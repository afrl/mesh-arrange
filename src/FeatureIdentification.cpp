#include "mesh_arrange/FeatureIdentification.h"

#include "detail/Containers.h"
#include "detail/DebugVTK.h"
#include "detail/ElementTables.h"
#include "detail/MathHelpers.h"
#include "igl/winding_number.h"
#include "mesh_arrange/GlobalOptions.h"
#include "mesh_arrange/Optimize.h"
#include "mesh_arrange/SplitEdges.h"

#include <Eigen/Dense>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <deque>
#include <iostream>
#include <limits>
#include <numeric>
#include <set>
#include <stdexcept>
#include <unordered_map>

namespace mesh_arrange
{

namespace
{

struct SortedTriFaceHash
{
    std::size_t operator()(const std::array<Index, 3>& k) const noexcept
    {
        std::size_t seed = 0;
        for (int i = 0; i < 3; ++i)
        {
            seed ^= std::hash<Index>{}(k[i]) + 0x9e3779b9ULL + (seed << 6) + (seed >> 2);
        }
        return seed;
    }
};

std::string join_output_path(const std::string& filename)
{
    const auto& dir = GlobalOptions::get().output_dir;
    if (dir.empty() || dir == ".")
    {
        return filename;
    }
    if (dir.back() == '/')
    {
        return dir + filename;
    }
    return dir + "/" + filename;
}

std::vector<detail::DynamicBitset> calculate_in_out_bitsets(
  const Triangulation& predicate_surface_mesh,
  const std::vector<std::vector<Index>>& object_to_surface_elements_map,
  const Eigen::MatrixXd& points)
{
    Eigen::MatrixX3d vertices(static_cast<Eigen::Index>(predicate_surface_mesh.vertices.size()), 3);
    for (Eigen::Index i = 0; i < vertices.rows(); ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            vertices(i, j) = predicate_surface_mesh.vertices[static_cast<std::size_t>(i)][j];
        }
    }

    std::vector<detail::DynamicBitset> in_out_bitsets(
      static_cast<std::size_t>(points.rows()),
      detail::DynamicBitset(object_to_surface_elements_map.size()));

    if (points.rows() == 0)
    {
        return in_out_bitsets;
    }

    Eigen::VectorXd winding_numbers(points.rows());
    for (std::size_t object_id = 0; object_id < object_to_surface_elements_map.size(); ++object_id)
    {
        const auto& elems = object_to_surface_elements_map[object_id];
        if (elems.empty())
        {
            continue;
        }

        Eigen::MatrixX3i triangles(static_cast<Eigen::Index>(elems.size()), 3);
        for (Eigen::Index i = 0; i < triangles.rows(); ++i)
        {
            const Index tri_id = elems[static_cast<std::size_t>(i)];
            for (int j = 0; j < 3; ++j)
            {
                triangles(i, j) = static_cast<int>(predicate_surface_mesh.triangles[tri_id][j]);
            }
        }

        igl::winding_number(vertices, triangles, points, winding_numbers);
        for (Eigen::Index i = 0; i < points.rows(); ++i)
        {
            if (std::abs(winding_numbers(i)) > 0.5)
            {
                in_out_bitsets[static_cast<std::size_t>(i)].set(object_id, true);
            }
        }
    }

    return in_out_bitsets;
}

void write_regions_vtk(const Tetrahedralization& tetrahedralization,
                       const RegionInfo& region_info,
                       const ObjectRegionInfo& object_region_info,
                       const std::string& prefix  = "",
                       const std::string& postfix = "")
{
    {
        std::vector<uint32_t> region_id(tetrahedralization.tetrahedra.size());
        std::vector<uint32_t> is_clipped(tetrahedralization.tetrahedra.size());
        for (std::size_t i = 0; i < tetrahedralization.tetrahedra.size(); ++i)
        {
            const RegionID rid = region_info.tet_to_region_map[i];
            region_id[i]       = rid;
            if (rid != NOT_SET_RegionID && rid < object_region_info.region_types.size())
            {
                is_clipped[i] =
                  (object_region_info.region_types[rid] == RegionType::CLIPPED) ? 1U : 0U;
            }
            else
            {
                is_clipped[i] = 1U;
            }
        }
        detail::write_vtk_tetrahedra(
          join_output_path(prefix + "regions" + postfix + ".vtk"),
          tetrahedralization.vertices,
          tetrahedralization.tetrahedra,
          { { "region_id", std::move(region_id) }, { "is_clipped", std::move(is_clipped) } });
    }

    {
        std::vector<std::array<Index, 4>> ov_tets;
        std::vector<uint32_t> ov_region_id;
        for (const RegionID overlap_region_id : object_region_info.regions_of_overlap)
        {
            for (const Index tid : region_info.region_to_tets_map[overlap_region_id])
            {
                ov_tets.push_back(tetrahedralization.tetrahedra[tid]);
                ov_region_id.push_back(region_info.tet_to_region_map[tid]);
            }
        }
        detail::write_vtk_tetrahedra(
          join_output_path(prefix + "overlap_regions" + postfix + ".vtk"),
          tetrahedralization.vertices,
          ov_tets,
          { { "region_id", std::move(ov_region_id) } });
    }
}

void write_region_boundaries_vtk(const Tetrahedralization& tetrahedralization,
                                 const RegionInfo& region_info,
                                 const std::vector<RegionBoundaryID>& triangle_to_boundary_id_map,
                                 const std::string& prefix  = "",
                                 const std::string& postfix = "")
{
    const std::size_t n_tris = tetrahedralization.triangles.size();
    std::vector<uint32_t> facet_id(n_tris);
    std::vector<uint32_t> boundary_id(n_tris);
    std::vector<uint32_t> bounded_region_1(n_tris);
    std::vector<uint32_t> bounded_region_2(n_tris);

    for (std::size_t i = 0; i < n_tris; ++i)
    {
        facet_id[i]         = region_info.tri_to_facet_map[i];
        boundary_id[i]      = triangle_to_boundary_id_map[i];
        bounded_region_1[i] = region_info.tri_to_bounded_regions_map[i].first;
        bounded_region_2[i] = region_info.tri_to_bounded_regions_map[i].second;
    }

    detail::write_vtk_triangles(join_output_path(prefix + "region_boundaries" + postfix + ".vtk"),
                                tetrahedralization.vertices,
                                tetrahedralization.triangles,
                                { { "facet_id", std::move(facet_id) },
                                  { "region_boundary_id", std::move(boundary_id) },
                                  { "bounded_region_1", std::move(bounded_region_1) },
                                  { "bounded_region_2", std::move(bounded_region_2) } });
}

void write_seams_vtk(const Tetrahedralization& tetrahedralization,
                     const SeamInfo& seam_info,
                     const std::string& prefix  = "",
                     const std::string& postfix = "")
{
    std::vector<std::array<Index, 2>> lines;
    std::vector<uint32_t> seam_ids;
    lines.reserve(seam_info.edges_on_seams.size());
    seam_ids.reserve(seam_info.edges_on_seams.size());

    for (const Index edge_id : seam_info.edges_on_seams)
    {
        lines.push_back(tetrahedralization.edges[edge_id]);
        seam_ids.push_back(seam_info.edge_id_to_seam_id_map[edge_id]);
    }

    detail::write_vtk_lines(join_output_path(prefix + "seams" + postfix + ".vtk"),
                            tetrahedralization.vertices,
                            lines,
                            { { "seam_id", std::move(seam_ids) } });
}

std::pair<SegmentRegionsResult, Eigen::MatrixXd> tetrahedralize_and_segment(
  const Triangulation& consistent_surface_mesh_to_segment,
  const PlcTetrahedralizerFn& plc_tetrahedralizer)
{
    if (GlobalOptions::get().verbose)
    {
        std::cout << GlobalOptions::get().curr_status_padding
                  << " - Tetrahedralizing surface for region segmentation" << std::endl;
    }

    const PlcTetrahedralizerFn& backend =
      plc_tetrahedralizer ? plc_tetrahedralizer : PlcTetrahedralizerFn(default_plc_tetrahedralizer);

    PlcTetOutput plc_out = backend(consistent_surface_mesh_to_segment);

    Tetrahedralization t;
    t.vertices = std::move(plc_out.vertices);

    t.tetrahedra.reserve(plc_out.tetrahedra.size());
    for (auto tet : plc_out.tetrahedra)
    {
        const double vol = detail::tetrahedron_signed_volume(t.vertices[tet[0]],
                                                             t.vertices[tet[1]],
                                                             t.vertices[tet[2]],
                                                             t.vertices[tet[3]]);
        if (vol < 0.0)
        {
            std::swap(tet[1], tet[2]);
        }
        t.tetrahedra.push_back(tet);
    }

    std::vector<FacetID> tri_to_facet_map;
    if (plc_out.surface_triangles.empty())
    {
        t.triangles      = consistent_surface_mesh_to_segment.triangles;
        tri_to_facet_map = consistent_surface_mesh_to_segment.tri_to_marker_map;
    }
    else
    {
        if (plc_out.surface_tri_parent_id.size() != plc_out.surface_triangles.size())
        {
            throw std::runtime_error(
              "PlcTetOutput::surface_tri_parent_id must match surface_triangles.size().");
        }
        t.triangles = std::move(plc_out.surface_triangles);
        tri_to_facet_map.reserve(t.triangles.size());
        for (const Index parent_id : plc_out.surface_tri_parent_id)
        {
            if (parent_id >= consistent_surface_mesh_to_segment.tri_to_marker_map.size())
            {
                throw std::runtime_error("PlcTetOutput::surface_tri_parent_id out of range.");
            }
            tri_to_facet_map.push_back(
              consistent_surface_mesh_to_segment.tri_to_marker_map[parent_id]);
        }
    }

    update_edges_and_maps(t);

    // Build tet face adjacency across all tetrahedra
    const std::size_t n_tets = t.tetrahedra.size();
    std::vector<std::array<Index, 4>> tet_adj(
      n_tets,
      { NOT_SET_Index, NOT_SET_Index, NOT_SET_Index, NOT_SET_Index });
    {
        std::unordered_map<std::array<Index, 3>, std::pair<Index, ConnIndex>, SortedTriFaceHash>
          first_tet_on_face;
        first_tet_on_face.reserve(n_tets * 3);

        for (std::size_t tid = 0; tid < n_tets; ++tid)
        {
            const auto& tc = t.tetrahedra[tid];
            for (ConnIndex f = 0; f < 4; ++f)
            {
                const auto& fc = detail::kTetFaces[f];
                std::array<Index, 3> key = { tc[fc[0]], tc[fc[1]], tc[fc[2]] };
                std::sort(key.begin(), key.end());
                auto [it, inserted] =
                  first_tet_on_face.emplace(key, std::make_pair(static_cast<Index>(tid), f));
                if (!inserted)
                {
                    const Index other_t        = it->second.first;
                    const ConnIndex other_f    = it->second.second;
                    tet_adj[tid][f]            = other_t;
                    tet_adj[other_t][other_f]  = static_cast<Index>(tid);
                }
            }
        }
    }

    // Segment connected tetrahedra without crossing constrained surface triangles
    std::vector<RegionID> tet_to_region_map(n_tets, NOT_SET_RegionID);
    std::vector<std::vector<Index>> region_to_tets_map;

    for (std::size_t seed = 0; seed < n_tets; ++seed)
    {
        if (tet_to_region_map[seed] != NOT_SET_RegionID)
        {
            continue;
        }
        const RegionID r_id = static_cast<RegionID>(region_to_tets_map.size());
        region_to_tets_map.emplace_back();

        std::deque<Index> bfs   = { static_cast<Index>(seed) };
        tet_to_region_map[seed] = r_id;
        while (!bfs.empty())
        {
            const Index cur = bfs.front();
            bfs.pop_front();
            region_to_tets_map[r_id].push_back(cur);

            for (ConnIndex f = 0; f < 4; ++f)
            {
                if (t.tet_to_tri_map[cur][f] == NOT_SET_Index)
                {
                    const Index nb = tet_adj[cur][f];
                    if (nb != NOT_SET_Index && tet_to_region_map[nb] == NOT_SET_RegionID)
                    {
                        tet_to_region_map[nb] = r_id;
                        bfs.push_back(nb);
                    }
                }
            }
        }
    }

    const RegionID n_regions = static_cast<RegionID>(region_to_tets_map.size());

    // Interior sample point per region: centroid of largest-volume tet in the region
    Eigen::MatrixXd regions_interior_point = Eigen::MatrixXd::Zero(n_regions, 3);
    std::vector<double> max_volume_tet_by_region(n_regions, -1.0);
    for (std::size_t i = 0; i < n_tets; ++i)
    {
        const RegionID r = tet_to_region_map[i];
        const auto& tc   = t.tetrahedra[i];
        const double vol = detail::tetrahedron_signed_volume(t.vertices[tc[0]],
                                                             t.vertices[tc[1]],
                                                             t.vertices[tc[2]],
                                                             t.vertices[tc[3]]);
        if (vol > max_volume_tet_by_region[r])
        {
            Eigen::Vector3d p = Eigen::Vector3d::Zero();
            for (const Index nid : tc)
            {
                p += detail::point_as_vector(t.vertices[nid]);
            }
            p /= 4.0;
            regions_interior_point.row(r) = p;
            max_volume_tet_by_region[r]   = vol;
        }
    }

    // Build region_to_bounding_triangles_map and tri_to_bounded_regions_map
    std::vector<std::vector<Index>> region_to_bounding_triangles_map(n_regions);
    std::vector<std::pair<RegionID, RegionID>> tri_to_bounded_regions_map(
      t.triangles.size(),
      { NOT_SET_RegionID, NOT_SET_RegionID });

    for (std::size_t tr = 0; tr < t.triangles.size(); ++tr)
    {
        const auto& tt = t.tri_to_tet_map[tr];
        if (tt[0] == NOT_SET_Index)
        {
            throw std::runtime_error(
              "PLC tetrahedralization is missing a constrained surface triangle (" +
              std::to_string(tr) + ").");
        }
        auto& br = tri_to_bounded_regions_map[tr];
        br.first = tet_to_region_map[tt[0]];
        region_to_bounding_triangles_map[br.first].push_back(static_cast<Index>(tr));

        if (tt[1] != NOT_SET_Index)
        {
            br.second = tet_to_region_map[tt[1]];
            if (br.second != br.first)
            {
                region_to_bounding_triangles_map[br.second].push_back(static_cast<Index>(tr));
            }
        }
        if (br.first > br.second)
        {
            std::swap(br.first, br.second);
        }
    }

    return { { std::move(t),
               { n_regions,
                 std::move(region_to_bounding_triangles_map),
                 std::move(tri_to_bounded_regions_map),
                 std::move(tet_to_region_map),
                 std::move(region_to_tets_map),
                 std::move(tri_to_facet_map) },
               ObjectRegionInfo{} },
             std::move(regions_interior_point) };
}

ObjectRegionInfo segment_regions_for_objects(
  const CombinedSurfaces& predicate_combined_surfaces,
  const Eigen::MatrixXd& regions_interior_point)
{
    const auto& predicate_surface_mesh = predicate_combined_surfaces.inconsistent_surface_mesh;
    const auto& object_to_predicate_surface_elements_map =
      predicate_combined_surfaces.object_to_surface_elements_map;
    const auto& object_types      = predicate_combined_surfaces.object_types;
    const auto& object_priorities = predicate_combined_surfaces.object_priorities;
    assert(object_types.size() == object_priorities.size());

    const std::size_t n_objects = object_to_predicate_surface_elements_map.size();
    const RegionID n_regions    = static_cast<RegionID>(regions_interior_point.rows());
    if (n_objects == 0 || n_regions == 0)
    {
        return { {},
                 std::vector<std::vector<RegionID>>(n_objects),
                 std::vector<std::vector<RegionID>>(n_regions),
                 std::vector<RegionType>(n_regions, RegionType::CLIPPED) };
    }

    const auto in_out_bitsets = calculate_in_out_bitsets(predicate_surface_mesh,
                                                         object_to_predicate_surface_elements_map,
                                                         regions_interior_point);

    std::vector<PriorityType> priorities_reversed = object_priorities;
    std::sort(priorities_reversed.begin(), priorities_reversed.end(), std::greater<PriorityType>());
    priorities_reversed.erase(std::unique(priorities_reversed.begin(), priorities_reversed.end()),
                              priorities_reversed.end());

    // Default any region outside all solid objects (e.g., cavities or exterior box air) to CLIPPED
    std::vector<std::vector<RegionID>> object_to_regions_map(n_objects);
    std::vector<RegionType> region_types(n_regions, RegionType::CLIPPED);
    std::vector<RegionID> regions_of_overlap;

    for (RegionID region_id = 0; region_id < n_regions; ++region_id)
    {
        for (const PriorityType p : priorities_reversed)
        {
            bool found_object = false;

            for (std::size_t object_id = 0; object_id < n_objects; ++object_id)
            {
                if (object_priorities[object_id] == p &&
                    ((object_types[object_id] == ObjectType::CLIP_INSIDE &&
                      in_out_bitsets[region_id].test(object_id)) ||
                     (object_types[object_id] == ObjectType::CLIP_OUTSIDE &&
                      !in_out_bitsets[region_id].test(object_id))))
                {
                    object_to_regions_map[object_id].push_back(region_id);
                    region_types[region_id] = RegionType::CLIPPED;
                    found_object            = true;
                }
            }

            if (found_object)
            {
                break;
            }

            int inside_solid_object_count = 0;
            for (std::size_t object_id = 0; object_id < n_objects; ++object_id)
            {
                if (object_priorities[object_id] == p &&
                    object_types[object_id] == ObjectType::SOLID &&
                    in_out_bitsets[region_id].test(object_id))
                {
                    object_to_regions_map[object_id].push_back(region_id);
                    region_types[region_id] = RegionType::KEPT;
                    found_object            = true;
                    ++inside_solid_object_count;
                }
            }

            if (found_object)
            {
                if (inside_solid_object_count > 1)
                {
                    regions_of_overlap.push_back(region_id);
                }
                break;
            }
        }
    }

    std::vector<std::vector<RegionID>> region_to_objects_map(n_regions);
    for (std::size_t object_id = 0; object_id < object_to_regions_map.size(); ++object_id)
    {
        for (const RegionID region_id : object_to_regions_map[object_id])
        {
            region_to_objects_map[region_id].push_back(static_cast<RegionID>(object_id));
        }
    }

    return { std::move(regions_of_overlap),
             std::move(object_to_regions_map),
             std::move(region_to_objects_map),
             std::move(region_types) };
}

void update_segment_regions_result(const RegionID& n_new_regions,
                                   const std::vector<RegionID>& old_to_new_region_id_map,
                                   const Index& n_new_tris,
                                   const std::vector<Index>& old_to_new_tri_map,
                                   Tetrahedralization& tetrahedralization,
                                   RegionInfo& region_info,
                                   ObjectRegionInfo& object_region_info)
{
    const std::size_t n_old_tris = region_info.tri_to_bounded_regions_map.size();

    region_info.n_regions = n_new_regions;
    std::vector<std::pair<RegionID, RegionID>> new_tri_to_bounded_regions_map(n_new_tris);
    std::vector<FacetID> new_tri_to_facet_map(n_new_tris);

    for (std::size_t old_tri_index = 0; old_tri_index < n_old_tris; ++old_tri_index)
    {
        if (old_to_new_tri_map[old_tri_index] != NOT_SET_Index)
        {
            const Index new_tri_index = old_to_new_tri_map[old_tri_index];
            new_tri_to_bounded_regions_map[new_tri_index] =
              region_info.tri_to_bounded_regions_map[old_tri_index];
            new_tri_to_facet_map[new_tri_index] = region_info.tri_to_facet_map[old_tri_index];

            auto& bounded_regions = new_tri_to_bounded_regions_map[new_tri_index];

            if (object_region_info.region_types[bounded_regions.first] == RegionType::CLIPPED)
            {
                if (bounded_regions.second != NOT_SET_RegionID &&
                    object_region_info.region_types[bounded_regions.second] != RegionType::CLIPPED)
                {
                    bounded_regions.first  = old_to_new_region_id_map[bounded_regions.second];
                    bounded_regions.second = NOT_SET_RegionID;
                }
            }
            else
            {
                bounded_regions.first = old_to_new_region_id_map[bounded_regions.first];
                if (bounded_regions.second != NOT_SET_RegionID &&
                    object_region_info.region_types[bounded_regions.second] == RegionType::CLIPPED)
                {
                    bounded_regions.second = NOT_SET_RegionID;
                }
                else if (bounded_regions.second != NOT_SET_RegionID)
                {
                    bounded_regions.second = old_to_new_region_id_map[bounded_regions.second];
                }
            }

            if (bounded_regions.first > bounded_regions.second)
            {
                std::swap(bounded_regions.first, bounded_regions.second);
            }
        }
    }

    std::swap(new_tri_to_bounded_regions_map, region_info.tri_to_bounded_regions_map);
    std::swap(new_tri_to_facet_map, region_info.tri_to_facet_map);

    region_info.region_to_bounding_triangles_map.assign(n_new_regions, {});
    for (std::size_t i = 0; i < region_info.tri_to_bounded_regions_map.size(); ++i)
    {
        assert(region_info.tri_to_bounded_regions_map[i].first != NOT_SET_RegionID);
        const RegionID region_id_1 = region_info.tri_to_bounded_regions_map[i].first;
        region_info.region_to_bounding_triangles_map[region_id_1].push_back(static_cast<Index>(i));
        if (region_info.tri_to_bounded_regions_map[i].second != NOT_SET_RegionID)
        {
            const RegionID region_id_2 = region_info.tri_to_bounded_regions_map[i].second;
            if (region_id_2 != region_id_1)
            {
                region_info.region_to_bounding_triangles_map[region_id_2].push_back(
                  static_cast<Index>(i));
            }
        }
    }

    Tetrahedralization new_t;
    new_t.vertices = tetrahedralization.vertices;

    new_t.triangles.resize(n_new_tris);
    for (std::size_t i = 0; i < tetrahedralization.triangles.size(); ++i)
    {
        if (old_to_new_tri_map[i] != NOT_SET_Index)
        {
            new_t.triangles[old_to_new_tri_map[i]] = tetrahedralization.triangles[i];
        }
    }

    Index curr_new_tet_id = 0;
    std::vector<Index> old_to_new_tet_id_map(tetrahedralization.tetrahedra.size(), NOT_SET_Index);
    std::vector<RegionID> new_tet_to_region_map;
    for (std::size_t old_tet_id = 0; old_tet_id < tetrahedralization.tetrahedra.size();
         ++old_tet_id)
    {
        const RegionID old_region_id = region_info.tet_to_region_map[old_tet_id];
        if (old_to_new_region_id_map[old_region_id] != NOT_SET_RegionID)
        {
            new_t.tetrahedra.push_back(tetrahedralization.tetrahedra[old_tet_id]);
            new_tet_to_region_map.push_back(old_to_new_region_id_map[old_region_id]);
            new_t.tet_to_tri_map.push_back(tetrahedralization.tet_to_tri_map[old_tet_id]);
            new_t.tet_to_edges_map.push_back(tetrahedralization.tet_to_edges_map[old_tet_id]);
            for (Index& t_idx : new_t.tet_to_tri_map.back())
            {
                if (t_idx != NOT_SET_Index)
                {
                    t_idx = old_to_new_tri_map[t_idx];
                }
            }
            old_to_new_tet_id_map[old_tet_id] = curr_new_tet_id++;
        }
    }
    std::swap(region_info.tet_to_region_map, new_tet_to_region_map);

    std::vector<std::vector<Index>> new_region_to_tets_map(n_new_regions);
    for (std::size_t old_region_id = 0; old_region_id < old_to_new_region_id_map.size();
         ++old_region_id)
    {
        const RegionID new_region_id = old_to_new_region_id_map[old_region_id];
        if (new_region_id != NOT_SET_RegionID)
        {
            for (const Index old_tet_id : region_info.region_to_tets_map[old_region_id])
            {
                const Index new_tet_id = old_to_new_tet_id_map[old_tet_id];
                if (new_tet_id != NOT_SET_Index)
                {
                    new_region_to_tets_map[new_region_id].push_back(new_tet_id);
                }
            }
        }
    }
    std::swap(region_info.region_to_tets_map, new_region_to_tets_map);

    for (std::size_t old_tri_id = 0; old_tri_id < n_old_tris; ++old_tri_id)
    {
        if (old_to_new_tri_map[old_tri_id] != NOT_SET_Index)
        {
            new_t.tri_to_tet_map.push_back(tetrahedralization.tri_to_tet_map[old_tri_id]);
            new_t.tri_to_edges_map.push_back(tetrahedralization.tri_to_edges_map[old_tri_id]);
            for (Index& t_idx : new_t.tri_to_tet_map.back())
            {
                if (t_idx != NOT_SET_Index)
                {
                    t_idx = old_to_new_tet_id_map[t_idx];
                }
            }
            std::sort(new_t.tri_to_tet_map.back().begin(), new_t.tri_to_tet_map.back().end());
        }
    }

    std::vector<bool> is_old_edge_kept(tetrahedralization.edges.size(), false);
    for (const auto& old_edges : new_t.tet_to_edges_map)
    {
        for (const Index old_e_id : old_edges)
        {
            is_old_edge_kept[old_e_id] = true;
        }
    }
    for (const auto& old_edges : new_t.tri_to_edges_map)
    {
        for (const Index old_e_id : old_edges)
        {
            is_old_edge_kept[old_e_id] = true;
        }
    }

    std::vector<Index> old_to_new_edges_map(tetrahedralization.edges.size(), NOT_SET_Index);
    Index curr_new_edge_id = 0;
    for (std::size_t old_e_id = 0; old_e_id < is_old_edge_kept.size(); ++old_e_id)
    {
        if (is_old_edge_kept[old_e_id])
        {
            new_t.edges.push_back(tetrahedralization.edges[old_e_id]);
            new_t.edge_to_tris_map.emplace_back();
            for (const auto& tp : tetrahedralization.edge_to_tris_map[old_e_id])
            {
                if (old_to_new_tri_map[tp.first] != NOT_SET_Index)
                {
                    new_t.edge_to_tris_map.back().emplace_back(old_to_new_tri_map[tp.first],
                                                               tp.second);
                }
            }
            new_t.edge_to_tets_map.emplace_back();
            for (const auto& tp : tetrahedralization.edge_to_tets_map[old_e_id])
            {
                if (old_to_new_tet_id_map[tp.first] != NOT_SET_Index)
                {
                    new_t.edge_to_tets_map.back().emplace_back(old_to_new_tet_id_map[tp.first],
                                                               tp.second);
                }
            }
            old_to_new_edges_map[old_e_id] = curr_new_edge_id++;
        }
    }

    for (auto& edges : new_t.tri_to_edges_map)
    {
        for (Index& e : edges)
        {
            e = old_to_new_edges_map[e];
        }
    }
    for (auto& edges : new_t.tet_to_edges_map)
    {
        for (Index& e : edges)
        {
            e = old_to_new_edges_map[e];
        }
    }

    tetrahedralization = std::move(new_t);

    ObjectRegionInfo new_object_region_info;
    new_object_region_info.object_to_regions_map.resize(
      object_region_info.object_to_regions_map.size());
    for (std::size_t object_id = 0; object_id < object_region_info.object_to_regions_map.size();
         ++object_id)
    {
        for (const RegionID old_region_id : object_region_info.object_to_regions_map[object_id])
        {
            if (old_to_new_region_id_map[old_region_id] != NOT_SET_RegionID)
            {
                new_object_region_info.object_to_regions_map[object_id].push_back(
                  old_to_new_region_id_map[old_region_id]);
            }
        }
    }
    for (auto& regions : new_object_region_info.object_to_regions_map)
    {
        std::sort(regions.begin(), regions.end());
        regions.erase(std::unique(regions.begin(), regions.end()), regions.end());
    }

    new_object_region_info.region_to_objects_map.resize(n_new_regions);
    new_object_region_info.region_types.resize(n_new_regions);
    for (std::size_t old_region_id = 0; old_region_id < old_to_new_region_id_map.size();
         ++old_region_id)
    {
        const RegionID new_region_id = old_to_new_region_id_map[old_region_id];
        if (new_region_id != NOT_SET_RegionID)
        {
            new_object_region_info.region_to_objects_map[new_region_id] =
              object_region_info.region_to_objects_map[old_region_id];
            new_object_region_info.region_types[new_region_id] =
              object_region_info.region_types[old_region_id];
        }
    }
    for (const RegionID old_region_id : object_region_info.regions_of_overlap)
    {
        if (old_to_new_region_id_map[old_region_id] != NOT_SET_RegionID)
        {
            new_object_region_info.regions_of_overlap.push_back(
              old_to_new_region_id_map[old_region_id]);
        }
    }
    std::sort(new_object_region_info.regions_of_overlap.begin(),
              new_object_region_info.regions_of_overlap.end());
    new_object_region_info.regions_of_overlap.erase(
      std::unique(new_object_region_info.regions_of_overlap.begin(),
                  new_object_region_info.regions_of_overlap.end()),
      new_object_region_info.regions_of_overlap.end());

    object_region_info = std::move(new_object_region_info);
}

void remove_clipped_regions(SegmentRegionsResult& srr)
{
    std::vector<RegionID> new_to_old_region_id_map;
    for (RegionID old_region_id = 0; old_region_id < srr.region_info.n_regions; ++old_region_id)
    {
        if (srr.object_region_info.region_types[old_region_id] != RegionType::CLIPPED)
        {
            new_to_old_region_id_map.push_back(old_region_id);
        }
    }

    std::vector<RegionID> old_to_new_region_id_map(srr.region_info.n_regions, NOT_SET_RegionID);
    for (std::size_t new_region_id = 0; new_region_id < new_to_old_region_id_map.size();
         ++new_region_id)
    {
        old_to_new_region_id_map[new_to_old_region_id_map[new_region_id]] =
          static_cast<RegionID>(new_region_id);
    }

    const std::size_t n_old_tris = srr.region_info.tri_to_bounded_regions_map.size();
    std::vector<Index> old_to_new_tri_map(n_old_tris, NOT_SET_Index);

    Index curr_new_tri_index = 0;
    for (std::size_t i = 0; i < n_old_tris; ++i)
    {
        const auto& bounded_regions = srr.region_info.tri_to_bounded_regions_map[i];
        if (srr.object_region_info.region_types[bounded_regions.first] == RegionType::CLIPPED)
        {
            if (bounded_regions.second != NOT_SET_RegionID &&
                srr.object_region_info.region_types[bounded_regions.second] != RegionType::CLIPPED)
            {
                old_to_new_tri_map[i] = curr_new_tri_index++;
            }
        }
        else
        {
            old_to_new_tri_map[i] = curr_new_tri_index++;
        }
    }

    update_segment_regions_result(static_cast<RegionID>(new_to_old_region_id_map.size()),
                                  old_to_new_region_id_map,
                                  curr_new_tri_index,
                                  old_to_new_tri_map,
                                  srr.tetrahedralization,
                                  srr.region_info,
                                  srr.object_region_info);
}

void remove_unneccessary_region_boundaries(
  const std::vector<PriorityType>& predicate_object_priorities,
  SegmentRegionsResult& srr,
  RegionBoundaryInfo& rbi,
  ObjectRegionBoundaryInfo& orbi)
{
    std::vector<RegionBoundaryID> region_boundaries_to_remove;
    for (std::size_t boundary_id = 0; boundary_id < orbi.boundary_id_to_objects_map.size();
         ++boundary_id)
    {
        const auto& bounded_regions = rbi.region_boundary_id_to_region_ids_map[boundary_id];
        if (bounded_regions.size() == 2)
        {
            PriorityType region_boundary_max_p = std::numeric_limits<PriorityType>::min();
            for (const RegionID object_id : orbi.boundary_id_to_objects_map[boundary_id])
            {
                if (predicate_object_priorities[object_id] > region_boundary_max_p)
                {
                    region_boundary_max_p = predicate_object_priorities[object_id];
                }
            }

            bool boundary_is_less_than_regions = true;
            for (const RegionID region_id : bounded_regions)
            {
                PriorityType region_max_p = std::numeric_limits<PriorityType>::min();
                for (const RegionID object_id :
                     srr.object_region_info.region_to_objects_map[region_id])
                {
                    assert(region_max_p == predicate_object_priorities[object_id] ||
                           region_max_p == std::numeric_limits<PriorityType>::min());
                    region_max_p = predicate_object_priorities[object_id];
                }
                if (region_boundary_max_p >= region_max_p)
                {
                    boundary_is_less_than_regions = false;
                }
            }
            if (boundary_is_less_than_regions)
            {
                region_boundaries_to_remove.push_back(static_cast<RegionBoundaryID>(boundary_id));
            }
        }
    }

    if (!region_boundaries_to_remove.empty())
    {
        remove_region_boundaries(region_boundaries_to_remove,
                                 srr.tetrahedralization,
                                 srr.region_info,
                                 srr.object_region_info,
                                 rbi,
                                 orbi);
    }
}

ObjectRegionBoundaryInfo find_object_and_region_boundary_info(
  const CombinedSurfaces& inconsistent_combined_surfaces,
  const ConsistentSurfaces& consistent_surfaces,
  const ObjectRegionInfo& object_region_info,
  const RegionBoundaryID& n_region_boundaries,
  const std::vector<std::vector<RegionBoundaryID>>& region_id_to_region_boundary_ids_map,
  const std::vector<FacetID>& region_boundary_id_to_consistent_facet_map)
{
    std::vector<std::vector<RegionBoundaryID>> object_to_external_boundary_ids_map;
    object_to_external_boundary_ids_map.reserve(object_region_info.object_to_regions_map.size());

    for (std::size_t object_id = 0; object_id < object_region_info.object_to_regions_map.size();
         ++object_id)
    {
        const auto& regions_in_object = object_region_info.object_to_regions_map[object_id];

        std::set<RegionBoundaryID> unmatched_boundary_ids;
        for (const RegionID region_id : regions_in_object)
        {
            for (const RegionBoundaryID boundary_id :
                 region_id_to_region_boundary_ids_map[region_id])
            {
                if (unmatched_boundary_ids.count(boundary_id) == 0)
                {
                    unmatched_boundary_ids.insert(boundary_id);
                }
                else
                {
                    unmatched_boundary_ids.erase(boundary_id);
                }
            }
        }

        object_to_external_boundary_ids_map.emplace_back();
        for (const RegionBoundaryID boundary_id : unmatched_boundary_ids)
        {
            const FacetID consistent_facet_id =
              region_boundary_id_to_consistent_facet_map[boundary_id];

            if (consistent_facet_id >=
                consistent_surfaces.constistent_facet_to_inconsistent_facets_map.size())
            {
                continue;
            }

            bool boundary_is_along_object = false;
            for (const FacetID predicate_facet_id :
                 consistent_surfaces
                   .constistent_facet_to_inconsistent_facets_map[consistent_facet_id])
            {
                for (const FacetID f : inconsistent_combined_surfaces.object_to_facets_map[object_id])
                {
                    if (f == predicate_facet_id)
                    {
                        boundary_is_along_object = true;
                        break;
                    }
                }
                if (boundary_is_along_object)
                {
                    break;
                }
            }
            if (boundary_is_along_object)
            {
                object_to_external_boundary_ids_map.back().push_back(boundary_id);
            }
        }
    }

    std::vector<std::vector<RegionID>> boundary_id_to_objects_map(n_region_boundaries);
    for (std::size_t object_id = 0; object_id < object_to_external_boundary_ids_map.size();
         ++object_id)
    {
        for (const RegionBoundaryID boundary_id : object_to_external_boundary_ids_map[object_id])
        {
            boundary_id_to_objects_map[boundary_id].push_back(static_cast<RegionID>(object_id));
        }
    }

    return { std::move(object_to_external_boundary_ids_map),
             std::move(boundary_id_to_objects_map) };
}

std::pair<bool, SegmentRegionsResult> attempt_optimize_then_segment_regions(
  const CombinedSurfaces& inconsistent_combined_surfaces,
  ConsistentSurfaces& consistent_surfaces,
  const PlcTetrahedralizerFn& plc_tetrahedralizer)
{
    try
    {
        if (GlobalOptions::get().verbose)
        {
            std::cout << " - Attempting to optimize surface before segmenting regions" << std::endl;
        }
        GlobalOptions::get().curr_status_padding = "   ";

        auto consistent_surface_copy = consistent_surfaces;
        auto seam_info               = segment_seams(consistent_surface_copy.consistent_surface_mesh);
        split_constrained_edges(consistent_surface_copy.consistent_surface_mesh, seam_info);
        optimize_tangential_smoothing(consistent_surface_copy.consistent_surface_mesh,
                                      seam_info,
                                      10,
                                      true);

        auto seg_region_result =
          segment_regions(inconsistent_combined_surfaces, consistent_surface_copy, plc_tetrahedralizer);

        std::swap(consistent_surfaces, consistent_surface_copy);
        GlobalOptions::get().curr_status_padding = "";

        if (GlobalOptions::get().verbose)
        {
            std::cout << " - Optimization strategy succeeded" << std::endl;
        }
        return { true, std::move(seg_region_result) };
    }
    catch (...)
    {
        GlobalOptions::get().curr_status_padding = "";
        if (GlobalOptions::get().verbose)
        {
            std::cout << " - Optimization strategy failed" << std::endl;
        }
        return { false, SegmentRegionsResult{} };
    }
}

} // namespace

void remove_region_boundaries(const std::vector<RegionBoundaryID>& region_boundaries_to_remove,
                              Tetrahedralization& tetrahedralization,
                              RegionInfo& region_info,
                              ObjectRegionInfo& object_region_info,
                              RegionBoundaryInfo& rbi,
                              ObjectRegionBoundaryInfo& orbi)
{
    if (region_boundaries_to_remove.empty())
    {
        return;
    }

    std::vector<bool> keep_region_boundary(rbi.n_region_boundaries, true);
    for (const RegionBoundaryID boundary_id : region_boundaries_to_remove)
    {
        keep_region_boundary[boundary_id] = false;
    }

    std::vector<RegionBoundaryID> old_to_new_region_boundary_map(rbi.n_region_boundaries,
                                                                 NOT_SET_RegionBoundaryID);
    RegionBoundaryID curr_new_boundary_id = 0;
    for (RegionBoundaryID boundary_id = 0; boundary_id < rbi.n_region_boundaries; ++boundary_id)
    {
        if (keep_region_boundary[boundary_id])
        {
            old_to_new_region_boundary_map[boundary_id] = curr_new_boundary_id++;
        }
    }

    std::vector<bool> triangle_is_kept(tetrahedralization.triangles.size(), true);
    for (const RegionBoundaryID old_boundary_id : region_boundaries_to_remove)
    {
        for (const Index t : rbi.region_boundary_to_triangles_map[old_boundary_id])
        {
            triangle_is_kept[t] = false;
        }
    }

    std::vector<Index> old_to_new_tri_map(tetrahedralization.triangles.size(), NOT_SET_Index);
    Index curr_new_tri_id = 0;
    for (std::size_t old_tri_id = 0; old_tri_id < triangle_is_kept.size(); ++old_tri_id)
    {
        if (triangle_is_kept[old_tri_id])
        {
            old_to_new_tri_map[old_tri_id] = curr_new_tri_id++;
        }
    }

    // Connected components of merged regions (using union-find so transitive merges unite properly)
    std::vector<RegionID> uf_parent(region_info.n_regions);
    std::iota(uf_parent.begin(), uf_parent.end(), 0);
    auto uf_find = [&](RegionID x)
    {
        while (uf_parent[x] != x)
        {
            uf_parent[x] = uf_parent[uf_parent[x]];
            x            = uf_parent[x];
        }
        return x;
    };
    std::vector<bool> is_region_being_merged(region_info.n_regions, false);
    for (const RegionBoundaryID boundary_id : region_boundaries_to_remove)
    {
        const auto& bounded_regions = rbi.region_boundary_id_to_region_ids_map[boundary_id];
        assert(bounded_regions.size() == 2);
        is_region_being_merged[bounded_regions[0]] = true;
        is_region_being_merged[bounded_regions[1]] = true;
        const RegionID ra                          = uf_find(bounded_regions[0]);
        const RegionID rb                          = uf_find(bounded_regions[1]);
        if (ra != rb)
        {
            uf_parent[std::max(ra, rb)] = std::min(ra, rb);
        }
    }

    std::vector<RegionID> old_to_new_region_map(region_info.n_regions, NOT_SET_RegionID);
    RegionID curr_new_region_id = 0;
    for (RegionID old_region_id = 0; old_region_id < region_info.n_regions; ++old_region_id)
    {
        if (!is_region_being_merged[old_region_id])
        {
            old_to_new_region_map[old_region_id] = curr_new_region_id++;
        }
    }
    std::unordered_map<RegionID, RegionID> root_to_new_id;
    for (RegionID old_region_id = 0; old_region_id < region_info.n_regions; ++old_region_id)
    {
        if (is_region_being_merged[old_region_id])
        {
            const RegionID root = uf_find(old_region_id);
            auto [it, ins]      = root_to_new_id.emplace(root, curr_new_region_id);
            if (ins)
            {
                ++curr_new_region_id;
            }
            old_to_new_region_map[old_region_id] = it->second;
        }
    }

    update_segment_regions_result(curr_new_region_id,
                                  old_to_new_region_map,
                                  curr_new_tri_id,
                                  old_to_new_tri_map,
                                  tetrahedralization,
                                  region_info,
                                  object_region_info);

    RegionBoundaryInfo new_rbi;
    new_rbi.n_region_boundaries = curr_new_boundary_id;

    new_rbi.region_boundary_id_to_region_ids_map.resize(curr_new_boundary_id);
    for (std::size_t old_boundary_id = 0; old_boundary_id < old_to_new_region_boundary_map.size();
         ++old_boundary_id)
    {
        const RegionBoundaryID new_boundary_id = old_to_new_region_boundary_map[old_boundary_id];
        if (new_boundary_id != NOT_SET_RegionBoundaryID)
        {
            auto& vec = new_rbi.region_boundary_id_to_region_ids_map[new_boundary_id];
            vec       = rbi.region_boundary_id_to_region_ids_map[old_boundary_id];
            for (RegionID& r : vec)
            {
                if (r != NOT_SET_RegionID)
                {
                    r = old_to_new_region_map[r];
                }
            }
            std::sort(vec.begin(), vec.end());
            vec.erase(std::unique(vec.begin(), vec.end()), vec.end());
        }
    }

    new_rbi.region_boundary_to_triangles_map.resize(curr_new_boundary_id);
    for (std::size_t old_boundary_id = 0; old_boundary_id < old_to_new_region_boundary_map.size();
         ++old_boundary_id)
    {
        const RegionBoundaryID new_boundary_id = old_to_new_region_boundary_map[old_boundary_id];
        if (new_boundary_id != NOT_SET_RegionBoundaryID)
        {
            new_rbi.region_boundary_to_triangles_map[new_boundary_id] =
              rbi.region_boundary_to_triangles_map[old_boundary_id];
            for (Index& t : new_rbi.region_boundary_to_triangles_map[new_boundary_id])
            {
                assert(old_to_new_tri_map[t] != NOT_SET_Index);
                t = old_to_new_tri_map[t];
            }
        }
    }

    new_rbi.triangle_to_boundary_id_map.resize(curr_new_tri_id, NOT_SET_RegionBoundaryID);
    for (std::size_t old_tri_id = 0; old_tri_id < old_to_new_tri_map.size(); ++old_tri_id)
    {
        const Index new_tri_id = old_to_new_tri_map[old_tri_id];
        if (new_tri_id != NOT_SET_Index)
        {
            const RegionBoundaryID old_boundary_id = rbi.triangle_to_boundary_id_map[old_tri_id];
            const RegionBoundaryID new_boundary_id =
              old_to_new_region_boundary_map[old_boundary_id];
            assert(new_boundary_id != NOT_SET_RegionBoundaryID);
            new_rbi.triangle_to_boundary_id_map[new_tri_id] = new_boundary_id;
        }
    }

    new_rbi.region_id_to_region_boundary_ids_map.resize(curr_new_region_id);
    for (std::size_t old_region_id = 0; old_region_id < old_to_new_region_map.size();
         ++old_region_id)
    {
        const RegionID new_region_id = old_to_new_region_map[old_region_id];
        assert(new_region_id != NOT_SET_RegionID);

        auto& new_boundaries = new_rbi.region_id_to_region_boundary_ids_map[new_region_id];
        for (const RegionBoundaryID old_boundary_id :
             rbi.region_id_to_region_boundary_ids_map[old_region_id])
        {
            const RegionBoundaryID new_boundary_id =
              old_to_new_region_boundary_map[old_boundary_id];
            if (new_boundary_id != NOT_SET_RegionBoundaryID)
            {
                new_boundaries.push_back(new_boundary_id);
            }
        }
    }
    for (auto& nb : new_rbi.region_id_to_region_boundary_ids_map)
    {
        std::sort(nb.begin(), nb.end());
        nb.erase(std::unique(nb.begin(), nb.end()), nb.end());
    }

    std::swap(rbi, new_rbi);

    std::vector<std::vector<RegionID>> new_boundary_id_to_objects_map(curr_new_boundary_id);
    for (std::size_t old_boundary_id = 0; old_boundary_id < old_to_new_region_boundary_map.size();
         ++old_boundary_id)
    {
        const RegionBoundaryID new_boundary_id = old_to_new_region_boundary_map[old_boundary_id];
        if (new_boundary_id != NOT_SET_RegionBoundaryID)
        {
            new_boundary_id_to_objects_map[new_boundary_id] =
              orbi.boundary_id_to_objects_map[old_boundary_id];
        }
    }
    std::swap(orbi.boundary_id_to_objects_map, new_boundary_id_to_objects_map);

    for (auto& boundary_ids : orbi.object_to_external_boundary_ids_map)
    {
        std::vector<RegionBoundaryID> kept_b;
        kept_b.reserve(boundary_ids.size());
        for (const RegionBoundaryID b : boundary_ids)
        {
            const RegionBoundaryID nb = old_to_new_region_boundary_map[b];
            if (nb != NOT_SET_RegionBoundaryID)
            {
                kept_b.push_back(nb);
            }
        }
        boundary_ids = std::move(kept_b);
    }
}

SegmentRegionsResult segment_regions(const CombinedSurfaces& inconsistent_combined_surfaces,
                                     const ConsistentSurfaces& consistent_surfaces,
                                     const PlcTetrahedralizerFn& plc_tetrahedralizer)
{
    auto seg_result = tetrahedralize_and_segment(consistent_surfaces.consistent_surface_mesh,
                                                 plc_tetrahedralizer);

    if (GlobalOptions::get().check_common_errors)
    {
        check_tri_edge_relations(seg_result.first.tetrahedralization);
        check_tet_edge_relations(seg_result.first.tetrahedralization);
        check_tet_tri_relations(seg_result.first.tetrahedralization);
    }

    seg_result.first.object_region_info =
      segment_regions_for_objects(inconsistent_combined_surfaces, seg_result.second);

    remove_clipped_regions(seg_result.first);

    if (GlobalOptions::get().write_intermediate_meshes)
    {
        write_regions_vtk(seg_result.first.tetrahedralization,
                          seg_result.first.region_info,
                          seg_result.first.object_region_info);
    }

    return std::move(seg_result.first);
}

SegmentRegionBoundariesResult segment_region_boundaries(
  const CombinedSurfaces& inconsistent_combined_surfaces,
  const ConsistentSurfaces& consistent_surfaces,
  const Tetrahedralization& tetrahedralization,
  const RegionInfo& region_info,
  const ObjectRegionInfo& object_region_info)
{
    std::unordered_map<std::size_t, RegionBoundaryID> boundary_hash_to_boundary_id_map;
    std::vector<std::vector<RegionBoundaryID>> region_id_to_boundary_ids_map(region_info.n_regions);
    std::vector<RegionBoundaryID> triangle_to_boundary_id_map(tetrahedralization.triangles.size());
    std::vector<FacetID> region_boundary_id_to_consistent_facet_map;

    for (std::size_t tri_id = 0; tri_id < tetrahedralization.triangles.size(); ++tri_id)
    {
        const FacetID facet_id      = region_info.tri_to_facet_map[tri_id];
        const auto& bounded_regions = region_info.tri_to_bounded_regions_map[tri_id];
        const std::size_t bounded_regions_hash =
          hash_linear_edge(bounded_regions.first, bounded_regions.second);
        const std::size_t boundary_hash =
          detail::hash_combine_pair(static_cast<std::size_t>(facet_id), bounded_regions_hash);

        const auto find_iter = boundary_hash_to_boundary_id_map.find(boundary_hash);
        if (find_iter == boundary_hash_to_boundary_id_map.end())
        {
            const RegionBoundaryID next_index =
              static_cast<RegionBoundaryID>(boundary_hash_to_boundary_id_map.size());

            if (bounded_regions.first != NOT_SET_RegionID)
            {
                region_id_to_boundary_ids_map[bounded_regions.first].push_back(next_index);
            }
            if (bounded_regions.second != NOT_SET_RegionID &&
                bounded_regions.second != bounded_regions.first)
            {
                region_id_to_boundary_ids_map[bounded_regions.second].push_back(next_index);
            }

            boundary_hash_to_boundary_id_map[boundary_hash] = next_index;
            triangle_to_boundary_id_map[tri_id]             = next_index;
            region_boundary_id_to_consistent_facet_map.push_back(facet_id);
        }
        else
        {
            triangle_to_boundary_id_map[tri_id] = find_iter->second;
        }
    }

    if (GlobalOptions::get().write_intermediate_meshes)
    {
        write_region_boundaries_vtk(tetrahedralization, region_info, triangle_to_boundary_id_map);
    }

    const RegionBoundaryID n_region_boundaries =
      static_cast<RegionBoundaryID>(boundary_hash_to_boundary_id_map.size());

    std::vector<std::vector<Index>> region_boundary_to_triangles_map(n_region_boundaries);
    for (std::size_t i = 0; i < triangle_to_boundary_id_map.size(); ++i)
    {
        region_boundary_to_triangles_map[triangle_to_boundary_id_map[i]].push_back(
          static_cast<Index>(i));
    }

    std::vector<std::vector<RegionID>> region_boundary_id_to_region_ids_map(n_region_boundaries);
    for (RegionID region_id = 0; region_id < region_info.n_regions; ++region_id)
    {
        for (const RegionBoundaryID region_boundary_id : region_id_to_boundary_ids_map[region_id])
        {
            region_boundary_id_to_region_ids_map[region_boundary_id].push_back(region_id);
        }
    }

    auto object_region_boundary_info =
      find_object_and_region_boundary_info(inconsistent_combined_surfaces,
                                           consistent_surfaces,
                                           object_region_info,
                                           n_region_boundaries,
                                           region_id_to_boundary_ids_map,
                                           region_boundary_id_to_consistent_facet_map);

    return { { n_region_boundaries,
               std::move(triangle_to_boundary_id_map),
               std::move(region_boundary_to_triangles_map),
               std::move(region_id_to_boundary_ids_map),
               std::move(region_boundary_id_to_region_ids_map) },
             std::move(object_region_boundary_info),
             std::move(region_boundary_id_to_consistent_facet_map) };
}

SeamInfo segment_seams(const Tetrahedralization& tetrahedralization,
                       const RegionBoundaryInfo& region_boundary_info)
{
    return segment_seams(tetrahedralization.edges,
                         tetrahedralization.tri_to_edges_map,
                         region_boundary_info.n_region_boundaries,
                         region_boundary_info.triangle_to_boundary_id_map);
}

SegmentedTetrahedralization identify_features(
  const CombinedSurfaces& inconsistent_combined_surfaces,
  ConsistentSurfaces& consistent_surfaces,
  const PlcTetrahedralizerFn& plc_tetrahedralizer)
{
    if (GlobalOptions::get().verbose)
    {
        std::cout << "* identify_features - Starting feature segmentation" << std::endl;
    }

    const auto stats = query_edge_length_stats(consistent_surfaces.consistent_surface_mesh);
    if (GlobalOptions::get().verbose)
    {
        std::cout << " - Edge length statistics:" << std::endl;
        std::cout << "     min: " << stats.min << std::endl;
        std::cout << "     max: " << stats.max << std::endl;
        std::cout << "     mean: " << stats.mean << std::endl;
        std::cout << "     std. dev.: " << stats.std_dev << std::endl;
    }

    SegmentRegionsResult seg_region_result;
    if (stats.min < 1e-12 && !consistent_surfaces.consistent_surface_mesh.triangles.empty())
    {
        auto result = attempt_optimize_then_segment_regions(inconsistent_combined_surfaces,
                                                            consistent_surfaces,
                                                            plc_tetrahedralizer);
        if (!result.first)
        {
            if (GlobalOptions::get().verbose)
            {
                std::cout << " - Segmenting consistent surface without modification" << std::endl;
            }
            seg_region_result =
              segment_regions(inconsistent_combined_surfaces, consistent_surfaces, plc_tetrahedralizer);
        }
        else
        {
            seg_region_result = std::move(result.second);
        }
    }
    else
    {
        if (GlobalOptions::get().verbose)
        {
            std::cout << " - Segmenting consistent surface without modification" << std::endl;
        }
        seg_region_result =
          segment_regions(inconsistent_combined_surfaces, consistent_surfaces, plc_tetrahedralizer);
    }

    if (GlobalOptions::get().verbose)
    {
        std::cout << " - Segmenting region boundaries" << std::endl;
    }
    auto region_boundary_result = segment_region_boundaries(inconsistent_combined_surfaces,
                                                            consistent_surfaces,
                                                            seg_region_result.tetrahedralization,
                                                            seg_region_result.region_info,
                                                            seg_region_result.object_region_info);

    if (GlobalOptions::get().write_intermediate_meshes)
    {
        write_region_boundaries_vtk(
          seg_region_result.tetrahedralization,
          seg_region_result.region_info,
          region_boundary_result.region_boundaries_info.triangle_to_boundary_id_map,
          "initial_");
        write_regions_vtk(seg_region_result.tetrahedralization,
                          seg_region_result.region_info,
                          seg_region_result.object_region_info,
                          "initial_");
    }

    if (GlobalOptions::get().verbose)
    {
        std::cout << " - Removing unnecessary region boundaries" << std::endl;
    }
    remove_unneccessary_region_boundaries(inconsistent_combined_surfaces.object_priorities,
                                          seg_region_result,
                                          region_boundary_result.region_boundaries_info,
                                          region_boundary_result.object_region_boundary_info);

    region_boundary_result = segment_region_boundaries(inconsistent_combined_surfaces,
                                                       consistent_surfaces,
                                                       seg_region_result.tetrahedralization,
                                                       seg_region_result.region_info,
                                                       seg_region_result.object_region_info);

    if (GlobalOptions::get().write_intermediate_meshes)
    {
        write_region_boundaries_vtk(
          seg_region_result.tetrahedralization,
          seg_region_result.region_info,
          region_boundary_result.region_boundaries_info.triangle_to_boundary_id_map,
          "final_");
        write_regions_vtk(seg_region_result.tetrahedralization,
                          seg_region_result.region_info,
                          seg_region_result.object_region_info,
                          "final_");
    }

    if (GlobalOptions::get().verbose)
    {
        std::cout << " - Segmenting seams" << std::endl;
    }
    auto seam_info = segment_seams(seg_region_result.tetrahedralization,
                                   region_boundary_result.region_boundaries_info);

    if (GlobalOptions::get().write_intermediate_meshes)
    {
        write_seams_vtk(seg_region_result.tetrahedralization, seam_info);
    }

    return { std::move(seg_region_result.tetrahedralization),
             std::move(seg_region_result.region_info),
             std::move(seg_region_result.object_region_info),
             std::move(region_boundary_result.region_boundaries_info),
             std::move(region_boundary_result.object_region_boundary_info),
             std::move(seam_info) };
}

void write_features_to_vtk(const SegmentedTetrahedralization& st,
                           const std::string prefix,
                           const std::string postfix)
{
    write_regions_vtk(st.tetrahedralization, st.region_info, st.object_region_info, prefix, postfix);
    write_region_boundaries_vtk(st.tetrahedralization,
                                st.region_info,
                                st.region_boundaries_info.triangle_to_boundary_id_map,
                                prefix,
                                postfix);
    write_seams_vtk(st.tetrahedralization, st.seams_info, prefix, postfix);
}

} // namespace mesh_arrange
