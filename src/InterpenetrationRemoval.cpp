#include "mesh_arrange/InterpenetrationRemoval.h"

#include "detail/DebugVTK.h"
#include "detail/ElementTables.h"
#include "detail/MathHelpers.h"
#include "igl/copyleft/cgal/remesh_self_intersections.h"
#include "mesh_arrange/GlobalOptions.h"
#include "mesh_arrange/QualityTriangulation.h"
#include "mesh_arrange/Refinement.h"

#ifdef MESHARRANGE_ENABLE_TETGEN
#include "tetgen.h"
#endif

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>
#include <Eigen/SparseLU>
#include <algorithm>
#include <atomic>
#include <cassert>
#include <cmath>
#include <deque>
#include <limits>
#include <map>
#include <numeric>
#include <set>
#include <stdexcept>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace mesh_arrange
{

namespace
{

struct OverlapRegionGroup
{
    RegionID group_id = 0;
    std::vector<RegionID> region_ids;
    std::set<RegionBoundaryID> internal_boundary_ids;
    std::set<RegionBoundaryID> external_boundary_ids;
};

std::vector<OverlapRegionGroup> find_overlap_region_groups(
  const RegionBoundaryInfo& region_boundary_info,
  const std::vector<RegionID>& all_overlap_regions)
{
    const auto& region_id_to_boundary_ids_map =
      region_boundary_info.region_id_to_region_boundary_ids_map;
    const auto& boundary_id_to_region_ids_map =
      region_boundary_info.region_boundary_id_to_region_ids_map;

    const std::set<RegionID> all_overlap_set(all_overlap_regions.begin(),
                                             all_overlap_regions.end());
    std::set<RegionID> remaining_overlap_regions = all_overlap_set;

    RegionID curr_group_id = 0;
    std::vector<OverlapRegionGroup> overlap_region_groups;
    while (!remaining_overlap_regions.empty())
    {
        std::vector<RegionID> overlap_regions_group;
        std::set<RegionBoundaryID> internal_boundary_ids;
        std::set<RegionBoundaryID> external_boundary_ids;

        std::set<RegionID> adjacency_to_search_queue = { *remaining_overlap_regions.begin() };
        while (!adjacency_to_search_queue.empty())
        {
            const RegionID curr_region_id = *adjacency_to_search_queue.begin();
            overlap_regions_group.push_back(curr_region_id);
            remaining_overlap_regions.erase(curr_region_id);
            adjacency_to_search_queue.erase(adjacency_to_search_queue.begin());

            for (const RegionBoundaryID boundary_id : region_id_to_boundary_ids_map[curr_region_id])
            {
                bool boundary_connects_to_overlap = false;
                for (const RegionID adjacent_region_id : boundary_id_to_region_ids_map[boundary_id])
                {
                    if (adjacent_region_id != curr_region_id &&
                        all_overlap_set.count(adjacent_region_id) > 0)
                    {
                        boundary_connects_to_overlap = true;
                        if (remaining_overlap_regions.count(adjacent_region_id) > 0)
                        {
                            adjacency_to_search_queue.insert(adjacent_region_id);
                        }
                    }
                }

                if (boundary_connects_to_overlap)
                {
                    internal_boundary_ids.insert(boundary_id);
                }
                else
                {
                    external_boundary_ids.insert(boundary_id);
                }
            }
        }

        overlap_region_groups.push_back({ curr_group_id++,
                                          std::move(overlap_regions_group),
                                          std::move(internal_boundary_ids),
                                          std::move(external_boundary_ids) });
    }
    return overlap_region_groups;
}

struct SortedTriKeyHash
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

double point_triangle_sq_dist(const Eigen::Vector3d& p,
                              const Eigen::Vector3d& a,
                              const Eigen::Vector3d& b,
                              const Eigen::Vector3d& c) noexcept
{
    const Eigen::Vector3d ab = b - a;
    const Eigen::Vector3d ac = c - a;
    const Eigen::Vector3d ap = p - a;

    const double d1 = ab.dot(ap);
    const double d2 = ac.dot(ap);
    if (d1 <= 0.0 && d2 <= 0.0)
    {
        return ap.squaredNorm();
    }

    const Eigen::Vector3d bp = p - b;
    const double d3          = ab.dot(bp);
    const double d4          = ac.dot(bp);
    if (d3 >= 0.0 && d4 <= d3)
    {
        return bp.squaredNorm();
    }

    const double vc = d1 * d4 - d3 * d2;
    if (vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0)
    {
        const double v = d1 / (d1 - d3);
        return (ap - v * ab).squaredNorm();
    }

    const Eigen::Vector3d cp = p - c;
    const double d5          = ab.dot(cp);
    const double d6          = ac.dot(cp);
    if (d6 >= 0.0 && d5 <= d6)
    {
        return cp.squaredNorm();
    }

    const double vb = d5 * d2 - d1 * d6;
    if (vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0)
    {
        const double w = d2 / (d2 - d6);
        return (ap - w * ac).squaredNorm();
    }

    const double va = d3 * d6 - d5 * d4;
    if (va <= 0.0 && (d4 - d3) >= 0.0 && (d5 - d6) >= 0.0)
    {
        const double w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
        return (bp - w * (c - b)).squaredNorm();
    }

    const double denom = 1.0 / (va + vb + vc);
    const double v     = vb * denom;
    const double w     = vc * denom;
    return (ap - ab * v - ac * w).squaredNorm();
}

std::array<Index, 4> make_positive_tet(const std::vector<Point3>& vertices,
                                       Index v0,
                                       Index v1,
                                       Index v2,
                                       Index v3)
{
    if (detail::tetrahedron_signed_volume(vertices[v0], vertices[v1], vertices[v2], vertices[v3]) <
        0.0)
    {
        return { v0, v2, v1, v3 };
    }
    return { v0, v1, v2, v3 };
}

void refine_overlap_interiors(SegmentedTetrahedralization& st)
{
    if (st.object_region_info.regions_of_overlap.empty())
    {
        return;
    }

    // Step 1: Pre-split only oversized edges on removable/internal overlap caps (and giant
    // clipper faces cutting an overlap), never touching surviving solid outer caps or seam curves.
    for (int pass = 0; pass < 3; ++pass)
    {
        const auto groups = find_overlap_region_groups(st.region_boundaries_info,
                                                       st.object_region_info.regions_of_overlap);
        std::vector<double> boundary_max_len(
          st.region_boundaries_info.n_region_boundaries,
          std::numeric_limits<double>::infinity());

        for (const auto& grp : groups)
        {
            std::set<RegionID> grp_objs;
            double grp_vol = 0.0;
            for (const RegionID rid : grp.region_ids)
            {
                for (const RegionID oid : st.object_region_info.region_to_objects_map[rid])
                {
                    grp_objs.insert(oid);
                }
                for (const Index tid : st.region_info.region_to_tets_map[rid])
                {
                    const auto& tc = st.tetrahedralization.tetrahedra[tid];
                    grp_vol += std::max(
                      0.0,
                      detail::tetrahedron_signed_volume(st.tetrahedralization.vertices[tc[0]],
                                                        st.tetrahedralization.vertices[tc[1]],
                                                        st.tetrahedralization.vertices[tc[2]],
                                                        st.tetrahedralization.vertices[tc[3]]));
                }
            }

            double sum_l2   = 0.0;
            double sum_l1   = 0.0;
            double grp_area = 0.0;
            for (const RegionBoundaryID bid : grp.external_boundary_ids)
            {
                for (const Index tri_id :
                     st.region_boundaries_info.region_boundary_to_triangles_map[bid])
                {
                    const auto& tr           = st.tetrahedralization.triangles[tri_id];
                    const Eigen::Vector3d p0 = detail::point_as_vector(
                      st.tetrahedralization.vertices[tr[0]]);
                    const Eigen::Vector3d p1 = detail::point_as_vector(
                      st.tetrahedralization.vertices[tr[1]]);
                    const Eigen::Vector3d p2 = detail::point_as_vector(
                      st.tetrahedralization.vertices[tr[2]]);
                    grp_area += 0.5 * (p1 - p0).cross(p2 - p0).norm();
                    const double l01 = (p1 - p0).norm();
                    const double l12 = (p2 - p1).norm();
                    const double l20 = (p0 - p2).norm();
                    sum_l2 += l01 * l01 + l12 * l12 + l20 * l20;
                    sum_l1 += l01 + l12 + l20;
                }
            }

            const double h_bnd    = (sum_l1 > 1e-14) ? (sum_l2 / sum_l1) : 1.0;
            const double h_lens   = (grp_area > 1e-14) ? (4.0 * grp_vol / grp_area) : h_bnd;
            const double h_floor  = std::cbrt(std::max(1e-24, 6.0 * grp_vol / 4500.0));
            const double h_target = std::max(std::min(h_bnd, 0.50 * h_lens), h_floor);

            std::map<RegionBoundaryID, int> ext_owner_count;
            for (const RegionBoundaryID bid : grp.external_boundary_ids)
            {
                ext_owner_count[bid] = 0;
            }
            for (const RegionID oid : grp_objs)
            {
                for (const RegionBoundaryID bid :
                     st.object_region_boundary_info.object_to_external_boundary_ids_map[oid])
                {
                    auto it = ext_owner_count.find(bid);
                    if (it != ext_owner_count.end())
                    {
                        it->second += 1;
                    }
                }
            }

            for (const RegionBoundaryID bid : grp.internal_boundary_ids)
            {
                boundary_max_len[bid] =
                  std::min(boundary_max_len[bid], 1.65 * h_target);
            }
            for (const auto& kv : ext_owner_count)
            {
                if (kv.second == 1)
                {
                    boundary_max_len[kv.first] =
                      std::min(boundary_max_len[kv.first], 1.65 * h_target);
                }
            }
        }

        std::vector<Index> edges_to_split;
        const auto& t = st.tetrahedralization;
        for (std::size_t eid = 0; eid < t.edges.size(); ++eid)
        {
            if (t.edge_to_tris_map[eid].empty())
            {
                continue;
            }
            if (eid < st.seams_info.edge_id_to_seam_id_map.size() &&
                st.seams_info.edge_id_to_seam_id_map[eid] != NOT_SET_SeamID)
            {
                continue;
            }

            double allowed = 0.0;
            bool eligible  = true;
            for (const auto& tp : t.edge_to_tris_map[eid])
            {
                const RegionBoundaryID bid =
                  st.region_boundaries_info.triangle_to_boundary_id_map[tp.first];
                if (bid >= boundary_max_len.size() || !std::isfinite(boundary_max_len[bid]))
                {
                    eligible = false;
                    break;
                }
                allowed = std::max(allowed, boundary_max_len[bid]);

                // Never midpoint-split an edge incident to a sub-tolerance CGAL rational sliver,
                // as double-precision rounding of 0.5*(u+w) can land collinearly on (v,w).
                const auto& tc           = t.triangles[tp.first];
                const Eigen::Vector3d p0 = detail::point_as_vector(t.vertices[tc[0]]);
                const Eigen::Vector3d p1 = detail::point_as_vector(t.vertices[tc[1]]);
                const Eigen::Vector3d p2 = detail::point_as_vector(t.vertices[tc[2]]);
                const double e01         = (p1 - p0).norm();
                const double e12         = (p2 - p1).norm();
                const double e20         = (p0 - p2).norm();
                const double emax        = std::max({ e01, e12, e20 });
                const double emin        = std::min({ e01, e12, e20 });
                const double dbl_a       = (p1 - p0).cross(p2 - p0).norm();
                if (emin <= 1e-5 * emax || dbl_a <= 1e-6 * emax * emax)
                {
                    eligible = false;
                    break;
                }
            }
            if (!eligible || allowed <= 0.0)
            {
                continue;
            }

            const double len =
              detail::distance(t.vertices[t.edges[eid][0]], t.vertices[t.edges[eid][1]]);
            if (len > allowed)
            {
                edges_to_split.push_back(static_cast<Index>(eid));
            }
        }

        if (edges_to_split.empty())
        {
            break;
        }
        split_edges(edges_to_split, st);
    }

    std::unordered_set<RegionID> tetgen_refined_regions;

#ifdef MESHARRANGE_ENABLE_TETGEN
    {
        struct ShellTriMeta
        {
            std::array<Index, 3> local_conn;
            int expected_sign = +1;
            Eigen::Vector3d a;
            Eigen::Vector3d b;
            Eigen::Vector3d c;
            Eigen::Vector3d bb_min;
            Eigen::Vector3d bb_max;
        };

        struct RegionReplacement
        {
            RegionID region_id = 0;
            std::vector<Point3> new_interior_pts;
            std::vector<std::array<Index, 4>> new_tets_encoded;
        };

        std::vector<RegionReplacement> replacements;
        replacements.reserve(st.object_region_info.regions_of_overlap.size());

        for (const RegionID rid : st.object_region_info.regions_of_overlap)
        {
            const auto& region_tets = st.region_info.region_to_tets_map[rid];
            if (region_tets.empty())
            {
                continue;
            }

            std::vector<Index> shell_tri_ids =
              st.region_info.region_to_bounding_triangles_map[rid];
            std::sort(shell_tri_ids.begin(), shell_tri_ids.end());
            shell_tri_ids.erase(std::unique(shell_tri_ids.begin(), shell_tri_ids.end()),
                                shell_tri_ids.end());
            if (shell_tri_ids.size() < 4)
            {
                continue;
            }

            std::vector<Index> local_to_global;
            std::unordered_map<Index, Index> global_to_local;
            for (const Index tr_id : shell_tri_ids)
            {
                for (const Index gv : st.tetrahedralization.triangles[tr_id])
                {
                    if (global_to_local.emplace(gv, static_cast<Index>(local_to_global.size()))
                          .second)
                    {
                        local_to_global.push_back(gv);
                    }
                }
            }

            std::vector<ShellTriMeta> shell_tris;
            shell_tris.reserve(shell_tri_ids.size());
            std::map<std::pair<Index, Index>, int> shell_edge_deg;
            bool shell_ok     = true;
            double shell_area = 0.0;
            double sum_l2     = 0.0;
            double sum_l1     = 0.0;

            for (const Index tr_id : shell_tri_ids)
            {
                const auto& gconn = st.tetrahedralization.triangles[tr_id];
                const auto& tt    = st.tetrahedralization.tri_to_tet_map[tr_id];
                Index adj_tet     = NOT_SET_Index;
                if (tt[0] != NOT_SET_Index && st.region_info.tet_to_region_map[tt[0]] == rid)
                {
                    adj_tet = tt[0];
                }
                else if (tt[1] != NOT_SET_Index && st.region_info.tet_to_region_map[tt[1]] == rid)
                {
                    adj_tet = tt[1];
                }
                if (adj_tet == NOT_SET_Index)
                {
                    shell_ok = false;
                    break;
                }

                Index v_opp = NOT_SET_Index;
                for (const Index cv : st.tetrahedralization.tetrahedra[adj_tet])
                {
                    if (cv != gconn[0] && cv != gconn[1] && cv != gconn[2])
                    {
                        v_opp = cv;
                        break;
                    }
                }
                if (v_opp == NOT_SET_Index)
                {
                    shell_ok = false;
                    break;
                }

                const double sv =
                  detail::tetrahedron_signed_volume(st.tetrahedralization.vertices[gconn[0]],
                                                    st.tetrahedralization.vertices[gconn[1]],
                                                    st.tetrahedralization.vertices[gconn[2]],
                                                    st.tetrahedralization.vertices[v_opp]);
                if (std::abs(sv) <= 1e-30)
                {
                    shell_ok = false;
                    break;
                }

                ShellTriMeta sm;
                sm.local_conn = { global_to_local.at(gconn[0]),
                                  global_to_local.at(gconn[1]),
                                  global_to_local.at(gconn[2]) };
                sm.expected_sign = (sv > 0.0) ? +1 : -1;
                sm.a = detail::point_as_vector(st.tetrahedralization.vertices[gconn[0]]);
                sm.b = detail::point_as_vector(st.tetrahedralization.vertices[gconn[1]]);
                sm.c = detail::point_as_vector(st.tetrahedralization.vertices[gconn[2]]);
                sm.bb_min = sm.a.cwiseMin(sm.b).cwiseMin(sm.c);
                sm.bb_max = sm.a.cwiseMax(sm.b).cwiseMax(sm.c);
                shell_tris.push_back(sm);

                shell_area += 0.5 * (sm.b - sm.a).cross(sm.c - sm.a).norm();
                const double l01 = (sm.b - sm.a).norm();
                const double l12 = (sm.c - sm.b).norm();
                const double l20 = (sm.a - sm.c).norm();
                sum_l2 += l01 * l01 + l12 * l12 + l20 * l20;
                sum_l1 += l01 + l12 + l20;

                for (int e = 0; e < 3; ++e)
                {
                    const Index u = sm.local_conn[e];
                    const Index v = sm.local_conn[(e + 1) % 3];
                    shell_edge_deg[{ std::min(u, v), std::max(u, v) }] += 1;
                }
            }

            if (!shell_ok)
            {
                continue;
            }
            for (const auto& kv : shell_edge_deg)
            {
                if (kv.second < 2 || (kv.second % 2) != 0)
                {
                    shell_ok = false;
                    break;
                }
            }
            if (!shell_ok)
            {
                continue;
            }

            double reg_vol = 0.0;
            for (const Index tid : region_tets)
            {
                const auto& tc = st.tetrahedralization.tetrahedra[tid];
                reg_vol += std::max(
                  0.0,
                  detail::tetrahedron_signed_volume(st.tetrahedralization.vertices[tc[0]],
                                                    st.tetrahedralization.vertices[tc[1]],
                                                    st.tetrahedralization.vertices[tc[2]],
                                                    st.tetrahedralization.vertices[tc[3]]));
            }
            if (reg_vol <= 1e-18 || shell_area <= 1e-16)
            {
                continue;
            }

            const double h_bnd   = (sum_l1 > 1e-14) ? (sum_l2 / sum_l1) : 1.0;
            const double h_lens  = 4.0 * reg_vol / shell_area;
            const double h_floor = std::cbrt(std::max(1e-24, 6.0 * reg_vol / 5500.0));
            const double h_int   = std::max(std::min(h_bnd, 0.42 * h_lens), h_floor);
            const double v_max =
              std::max((std::sqrt(2.0) / 12.0) * h_int * h_int * h_int, reg_vol / 7000.0);

            // Collect well-cleared interior seed points from cross-shell edges and tet centroids
            std::vector<Eigen::Vector3d> cand_seeds;
            std::unordered_set<Index> seen_int_edges;
            for (const Index tid : region_tets)
            {
                for (const Index eid : st.tetrahedralization.tet_to_edges_map[tid])
                {
                    if (st.tetrahedralization.edge_to_tris_map[eid].empty() &&
                        seen_int_edges.insert(eid).second)
                    {
                        const auto& ed = st.tetrahedralization.edges[eid];
                        cand_seeds.push_back(0.5 *
                                             (detail::point_as_vector(
                                                st.tetrahedralization.vertices[ed[0]]) +
                                              detail::point_as_vector(
                                                st.tetrahedralization.vertices[ed[1]])));
                    }
                }
                const auto& tc = st.tetrahedralization.tetrahedra[tid];
                cand_seeds.push_back(
                  0.25 * (detail::point_as_vector(st.tetrahedralization.vertices[tc[0]]) +
                          detail::point_as_vector(st.tetrahedralization.vertices[tc[1]]) +
                          detail::point_as_vector(st.tetrahedralization.vertices[tc[2]]) +
                          detail::point_as_vector(st.tetrahedralization.vertices[tc[3]])));
            }

            const double min_clear    = std::max(1e-9, 0.14 * std::min(h_int, h_lens));
            const double min_clear_sq = min_clear * min_clear;
            const double min_sep      = 0.35 * h_int;
            const double min_sep_sq   = min_sep * min_sep;

            std::vector<Eigen::Vector3d> accepted_seeds;
            accepted_seeds.reserve(std::min<std::size_t>(cand_seeds.size(), 4000));
            for (const auto& sp : cand_seeds)
            {
                if (accepted_seeds.size() >= 4000)
                {
                    break;
                }
                bool too_close = false;
                for (const auto& kp : accepted_seeds)
                {
                    if ((sp - kp).squaredNorm() < min_sep_sq)
                    {
                        too_close = true;
                        break;
                    }
                }
                if (too_close)
                {
                    continue;
                }

                for (const auto& sm : shell_tris)
                {
                    if (sp.x() < sm.bb_min.x() - min_clear || sp.x() > sm.bb_max.x() + min_clear ||
                        sp.y() < sm.bb_min.y() - min_clear || sp.y() > sm.bb_max.y() + min_clear ||
                        sp.z() < sm.bb_min.z() - min_clear || sp.z() > sm.bb_max.z() + min_clear)
                    {
                        continue;
                    }
                    if (point_triangle_sq_dist(sp, sm.a, sm.b, sm.c) < min_clear_sq)
                    {
                        too_close = true;
                        break;
                    }
                }
                if (!too_close)
                {
                    accepted_seeds.push_back(sp);
                }
            }

            try
            {
                using namespace mesh_arrange::tetgen_internal;

                const int n_shell_verts = static_cast<int>(local_to_global.size());
                const int n_seed_verts  = static_cast<int>(accepted_seeds.size());
                tetgenio in;
                in.firstnumber    = 0;
                in.numberofpoints = n_shell_verts + n_seed_verts;
                in.pointlist      = new REAL[in.numberofpoints * 3];
                for (int i = 0; i < n_shell_verts; ++i)
                {
                    const auto& pt =
                      st.tetrahedralization.vertices[local_to_global[static_cast<std::size_t>(i)]];
                    in.pointlist[3 * i + 0] = pt[0];
                    in.pointlist[3 * i + 1] = pt[1];
                    in.pointlist[3 * i + 2] = pt[2];
                }
                for (int i = 0; i < n_seed_verts; ++i)
                {
                    const int dst               = n_shell_verts + i;
                    in.pointlist[3 * dst + 0]   = accepted_seeds[static_cast<std::size_t>(i)].x();
                    in.pointlist[3 * dst + 1]   = accepted_seeds[static_cast<std::size_t>(i)].y();
                    in.pointlist[3 * dst + 2]   = accepted_seeds[static_cast<std::size_t>(i)].z();
                }

                in.numberoffacets  = static_cast<int>(shell_tris.size());
                in.facetlist       = new tetgenio::facet[in.numberoffacets];
                in.facetmarkerlist = new int[in.numberoffacets];
                for (int i = 0; i < in.numberoffacets; ++i)
                {
                    in.facetmarkerlist[i] = i + 1;
                    tetgenio::facet* f    = &in.facetlist[i];
                    f->numberofpolygons   = 1;
                    f->polygonlist        = new tetgenio::polygon[1];
                    f->numberofholes      = 0;
                    f->holelist           = nullptr;
                    tetgenio::polygon* p  = &f->polygonlist[0];
                    p->numberofvertices   = 3;
                    p->vertexlist         = new int[3];
                    const auto& lc        = shell_tris[static_cast<std::size_t>(i)].local_conn;
                    p->vertexlist[0]      = static_cast<int>(lc[0]);
                    p->vertexlist[1]      = static_cast<int>(lc[1]);
                    p->vertexlist[2]      = static_cast<int>(lc[2]);
                }

                tetgenbehavior flags;
                flags.plc            = 1;
                flags.object         = tetgenbehavior::POLY;
                flags.nobisect       = 1;
                flags.nomergefacet   = 1;
                flags.nomergevertex  = 1;
                flags.nojettison     = 1;
                flags.quality        = 1;
                flags.minratio       = 1.35;
                flags.fixedvolume    = 1;
                flags.maxvolume      = v_max;
                flags.smooth_maxiter = 0;
                flags.noimprove      = 0;
                flags.steinerleft    = 15000;
                flags.quiet          = 1;
                flags.verbose        = 0;
                flags.epsilon        = 1e-14;

                tetgenio out;
                tetrahedralize(&flags, &in, &out);

                if (out.numberofpoints < n_shell_verts || out.numberoftetrahedra <= 0)
                {
                    continue;
                }

                std::vector<Point3> out_pts(static_cast<std::size_t>(out.numberofpoints));
                bool coords_match = true;
                for (int i = 0; i < out.numberofpoints; ++i)
                {
                    out_pts[static_cast<std::size_t>(i)] = { out.pointlist[3 * i + 0],
                                                             out.pointlist[3 * i + 1],
                                                             out.pointlist[3 * i + 2] };
                    if (i < n_shell_verts)
                    {
                        const auto& orig_p =
                          st.tetrahedralization.vertices[local_to_global[static_cast<std::size_t>(i)]];
                        if (detail::distance(out_pts[static_cast<std::size_t>(i)], orig_p) > 1e-11)
                        {
                            coords_match = false;
                            break;
                        }
                    }
                }
                if (!coords_match)
                {
                    continue;
                }

                const std::size_t n_out_tets = static_cast<std::size_t>(out.numberoftetrahedra);
                std::vector<std::array<Index, 4>> out_tets(n_out_tets);
                for (std::size_t ti = 0; ti < n_out_tets; ++ti)
                {
                    out_tets[ti] = make_positive_tet(
                      out_pts,
                      static_cast<Index>(out.tetrahedronlist[4 * ti + 0]),
                      static_cast<Index>(out.tetrahedronlist[4 * ti + 1]),
                      static_cast<Index>(out.tetrahedronlist[4 * ti + 2]),
                      static_cast<Index>(out.tetrahedronlist[4 * ti + 3]));
                }

                std::unordered_map<std::array<Index, 3>, std::size_t, SortedTriKeyHash> shell_key_map;
                shell_key_map.reserve(shell_tris.size() * 2);
                for (std::size_t si = 0; si < shell_tris.size(); ++si)
                {
                    auto key = shell_tris[si].local_conn;
                    std::sort(key.begin(), key.end());
                    shell_key_map.emplace(key, si);
                }

                std::vector<std::array<Index, 4>> out_adj(
                  n_out_tets,
                  { NOT_SET_Index, NOT_SET_Index, NOT_SET_Index, NOT_SET_Index });
                std::unordered_map<std::array<Index, 3>,
                                   std::vector<std::pair<Index, ConnIndex>>,
                                   SortedTriKeyHash>
                  face_inc;
                face_inc.reserve(n_out_tets * 3);

                for (std::size_t ti = 0; ti < n_out_tets; ++ti)
                {
                    const auto& tc = out_tets[ti];
                    for (ConnIndex f = 0; f < 4; ++f)
                    {
                        const auto& fc           = detail::kTetFaces[f];
                        std::array<Index, 3> key = { tc[fc[0]], tc[fc[1]], tc[fc[2]] };
                        std::sort(key.begin(), key.end());
                        face_inc[key].emplace_back(static_cast<Index>(ti), f);
                    }
                }
                for (const auto& kv : face_inc)
                {
                    if (kv.second.size() == 2)
                    {
                        const auto [t0, f0] = kv.second[0];
                        const auto [t1, f1] = kv.second[1];
                        out_adj[t0][f0]     = t1;
                        out_adj[t1][f1]     = t0;
                    }
                }

                std::vector<int> tet_state(n_out_tets, 0); // 0=unvisited, 1=inside, -1=outside
                std::deque<Index> bfs;
                bool all_shell_matched = true;
                for (std::size_t si = 0; si < shell_tris.size(); ++si)
                {
                    const auto& sm = shell_tris[si];
                    auto key       = sm.local_conn;
                    std::sort(key.begin(), key.end());
                    const auto it = face_inc.find(key);
                    if (it == face_inc.end())
                    {
                        all_shell_matched = false;
                        break;
                    }
                    bool found_in = false;
                    for (const auto& [ti, lf] : it->second)
                    {
                        const auto& tc  = out_tets[ti];
                        Index v_opp_loc = NOT_SET_Index;
                        for (const Index cv : tc)
                        {
                            if (cv != sm.local_conn[0] && cv != sm.local_conn[1] &&
                                cv != sm.local_conn[2])
                            {
                                v_opp_loc = cv;
                                break;
                            }
                        }
                        if (v_opp_loc == NOT_SET_Index)
                        {
                            continue;
                        }
                        const double sv = detail::tetrahedron_signed_volume(
                          out_pts[sm.local_conn[0]],
                          out_pts[sm.local_conn[1]],
                          out_pts[sm.local_conn[2]],
                          out_pts[v_opp_loc]);
                        const int sgn = (sv > 0.0) ? +1 : -1;
                        if (sgn == sm.expected_sign)
                        {
                            found_in = true;
                            if (tet_state[ti] == -1)
                            {
                                all_shell_matched = false;
                                break;
                            }
                            if (tet_state[ti] == 0)
                            {
                                tet_state[ti] = 1;
                                bfs.push_back(ti);
                            }
                        }
                        else
                        {
                            if (tet_state[ti] == 1)
                            {
                                all_shell_matched = false;
                                break;
                            }
                            tet_state[ti] = -1;
                        }
                    }
                    if (!found_in || !all_shell_matched)
                    {
                        all_shell_matched = false;
                        break;
                    }
                }
                if (!all_shell_matched)
                {
                    continue;
                }

                while (!bfs.empty())
                {
                    const Index cur = bfs.front();
                    bfs.pop_front();
                    const auto& tc = out_tets[cur];
                    for (ConnIndex f = 0; f < 4; ++f)
                    {
                        const auto& fc           = detail::kTetFaces[f];
                        std::array<Index, 3> key = { tc[fc[0]], tc[fc[1]], tc[fc[2]] };
                        std::sort(key.begin(), key.end());
                        if (shell_key_map.count(key) > 0)
                        {
                            continue;
                        }
                        const Index nb = out_adj[cur][f];
                        if (nb == NOT_SET_Index || tet_state[nb] == -1)
                        {
                            all_shell_matched = false;
                            break;
                        }
                        if (tet_state[nb] == 0)
                        {
                            tet_state[nb] = 1;
                            bfs.push_back(nb);
                        }
                    }
                    if (!all_shell_matched)
                    {
                        break;
                    }
                }
                if (!all_shell_matched)
                {
                    continue;
                }

                std::vector<std::array<Index, 4>> kept_tets;
                kept_tets.reserve(n_out_tets);
                double kept_vol = 0.0;
                bool vols_pos   = true;
                for (std::size_t ti = 0; ti < n_out_tets; ++ti)
                {
                    if (tet_state[ti] != 1)
                    {
                        continue;
                    }
                    const auto& tc = out_tets[ti];
                    const double v = detail::tetrahedron_signed_volume(out_pts[tc[0]],
                                                                       out_pts[tc[1]],
                                                                       out_pts[tc[2]],
                                                                       out_pts[tc[3]]);
                    if (v <= 1e-22)
                    {
                        vols_pos = false;
                        break;
                    }
                    kept_vol += v;
                    kept_tets.push_back(tc);
                }
                if (!vols_pos || kept_tets.empty() ||
                    std::abs(kept_vol - reg_vol) > 1e-5 * std::max(1.0, reg_vol))
                {
                    continue;
                }

                // Smooth newly added strictly-interior vertices (`v >= n_shell_verts`) with
                // positive-volume guard to optimize tet aspect ratios inside the lens.
                const std::size_t n_pts_loc = out_pts.size();
                if (n_pts_loc > static_cast<std::size_t>(n_shell_verts))
                {
                    std::vector<std::vector<Index>> v_inc_tets(n_pts_loc);
                    std::vector<std::vector<Index>> v_nbrs(n_pts_loc);
                    for (std::size_t ti = 0; ti < kept_tets.size(); ++ti)
                    {
                        const auto& tc = kept_tets[ti];
                        for (int c = 0; c < 4; ++c)
                        {
                            if (tc[c] >= static_cast<Index>(n_shell_verts))
                            {
                                v_inc_tets[tc[c]].push_back(static_cast<Index>(ti));
                                for (int d = 0; d < 4; ++d)
                                {
                                    if (d != c)
                                    {
                                        v_nbrs[tc[c]].push_back(tc[d]);
                                    }
                                }
                            }
                        }
                    }
                    for (std::size_t vi = static_cast<std::size_t>(n_shell_verts); vi < n_pts_loc;
                         ++vi)
                    {
                        auto& nb = v_nbrs[vi];
                        std::sort(nb.begin(), nb.end());
                        nb.erase(std::unique(nb.begin(), nb.end()), nb.end());
                    }

                    for (int iter = 0; iter < 4; ++iter)
                    {
                        for (std::size_t vi = static_cast<std::size_t>(n_shell_verts);
                             vi < n_pts_loc;
                             ++vi)
                        {
                            if (v_nbrs[vi].empty() || v_inc_tets[vi].empty())
                            {
                                continue;
                            }
                            Point3 avg = { 0.0, 0.0, 0.0 };
                            for (const Index nb : v_nbrs[vi])
                            {
                                avg[0] += out_pts[nb][0];
                                avg[1] += out_pts[nb][1];
                                avg[2] += out_pts[nb][2];
                            }
                            const double inv_deg = 1.0 / static_cast<double>(v_nbrs[vi].size());
                            avg[0] *= inv_deg;
                            avg[1] *= inv_deg;
                            avg[2] *= inv_deg;

                            const Point3 old_p = out_pts[vi];
                            double min_before  = std::numeric_limits<double>::max();
                            for (const Index ti : v_inc_tets[vi])
                            {
                                const auto& tc = kept_tets[ti];
                                min_before     = std::min(min_before,
                                                      detail::tetrahedron_signed_volume(
                                                        out_pts[tc[0]],
                                                        out_pts[tc[1]],
                                                        out_pts[tc[2]],
                                                        out_pts[tc[3]]));
                            }

                            for (const double alpha : { 0.5, 0.25 })
                            {
                                out_pts[vi] = { (1.0 - alpha) * old_p[0] + alpha * avg[0],
                                                (1.0 - alpha) * old_p[1] + alpha * avg[1],
                                                (1.0 - alpha) * old_p[2] + alpha * avg[2] };
                                bool ok     = true;
                                for (const Index ti : v_inc_tets[vi])
                                {
                                    const auto& tc = kept_tets[ti];
                                    const double nv =
                                      detail::tetrahedron_signed_volume(out_pts[tc[0]],
                                                                        out_pts[tc[1]],
                                                                        out_pts[tc[2]],
                                                                        out_pts[tc[3]]);
                                    if (nv < 0.5 * min_before || nv <= 1e-22)
                                    {
                                        ok = false;
                                        break;
                                    }
                                }
                                if (ok)
                                {
                                    break;
                                }
                                out_pts[vi] = old_p;
                            }
                        }
                    }
                }

                // Compact used interior vertices and encode indices (global or offset-flagged)
                RegionReplacement rep;
                rep.region_id = rid;
                std::unordered_map<Index, Index> int_remap;
                rep.new_tets_encoded.reserve(kept_tets.size());
                constexpr Index kIntFlag = (Index{ 1 } << (sizeof(Index) * 8 - 1));
                for (const auto& tc : kept_tets)
                {
                    std::array<Index, 4> enc;
                    for (int c = 0; c < 4; ++c)
                    {
                        if (tc[c] < static_cast<Index>(n_shell_verts))
                        {
                            enc[c] = local_to_global[tc[c]];
                        }
                        else
                        {
                            auto [it_m, ins] = int_remap.emplace(
                              tc[c],
                              static_cast<Index>(rep.new_interior_pts.size()));
                            if (ins)
                            {
                                rep.new_interior_pts.push_back(out_pts[tc[c]]);
                            }
                            enc[c] = kIntFlag | it_m->second;
                        }
                    }
                    rep.new_tets_encoded.push_back(enc);
                }

                replacements.push_back(std::move(rep));
                tetgen_refined_regions.insert(rid);
            }
            catch (...)
            {
                // Fall through cleanly to edge-splitting fallback for this region
            }
        }

        if (!replacements.empty())
        {
            std::vector<std::array<Index, 4>> merged_tets;
            std::vector<RegionID> merged_tet_region;
            merged_tets.reserve(st.tetrahedralization.tetrahedra.size() * 2);
            merged_tet_region.reserve(st.tetrahedralization.tetrahedra.size() * 2);

            for (std::size_t ti = 0; ti < st.tetrahedralization.tetrahedra.size(); ++ti)
            {
                const RegionID r = st.region_info.tet_to_region_map[ti];
                if (tetgen_refined_regions.count(r) == 0)
                {
                    merged_tets.push_back(st.tetrahedralization.tetrahedra[ti]);
                    merged_tet_region.push_back(r);
                }
            }

            constexpr Index kIntFlag = (Index{ 1 } << (sizeof(Index) * 8 - 1));
            for (const auto& rep : replacements)
            {
                const Index base_v = static_cast<Index>(st.tetrahedralization.vertices.size());
                for (const auto& p : rep.new_interior_pts)
                {
                    st.tetrahedralization.vertices.push_back(p);
                }
                for (auto enc : rep.new_tets_encoded)
                {
                    for (int c = 0; c < 4; ++c)
                    {
                        if ((enc[c] & kIntFlag) != 0)
                        {
                            enc[c] = base_v + (enc[c] & ~kIntFlag);
                        }
                    }
                    merged_tets.push_back(enc);
                    merged_tet_region.push_back(rep.region_id);
                }
            }

            st.tetrahedralization.tetrahedra = std::move(merged_tets);
            update_edges_and_maps(st.tetrahedralization);

            // Rebuild seam edge IDs on updated edge indices
            st.seams_info = segment_seams(st.tetrahedralization, st.region_boundaries_info);

            st.region_info.tet_to_region_map = std::move(merged_tet_region);
            st.region_info.region_to_tets_map.assign(st.region_info.n_regions, {});
            for (std::size_t ti = 0; ti < st.tetrahedralization.tetrahedra.size(); ++ti)
            {
                st.region_info.region_to_tets_map[st.region_info.tet_to_region_map[ti]].push_back(
                  static_cast<Index>(ti));
            }
        }
    }
#endif

    // Step 3: Split any remaining strictly-interior cross-shell edges inside overlap regions
    // (touches ZERO surface triangles, ensuring thin seam wedges always have an interior DOF).
    {
        std::vector<bool> is_surface_vert(st.tetrahedralization.vertices.size(), false);
        for (const auto& tr : st.tetrahedralization.triangles)
        {
            is_surface_vert[tr[0]] = true;
            is_surface_vert[tr[1]] = true;
            is_surface_vert[tr[2]] = true;
        }

        std::set<Index> cross_shell_interior_edges;
        for (const RegionID rid : st.object_region_info.regions_of_overlap)
        {
            for (const Index tid : st.region_info.region_to_tets_map[rid])
            {
                for (const Index eid : st.tetrahedralization.tet_to_edges_map[tid])
                {
                    if (st.tetrahedralization.edge_to_tris_map[eid].empty())
                    {
                        const auto& ed = st.tetrahedralization.edges[eid];
                        if (is_surface_vert[ed[0]] && is_surface_vert[ed[1]])
                        {
                            cross_shell_interior_edges.insert(eid);
                        }
                    }
                }
            }
        }
        if (!cross_shell_interior_edges.empty())
        {
            const std::vector<Index> to_split(cross_shell_interior_edges.begin(),
                                              cross_shell_interior_edges.end());
            split_edges(to_split, st);
        }
    }

    // Step 4: Fallback for any overlap region that was not TetGen-refined and still lacks interior
    // vertices (e.g. MESHARRANGE_WITH_TETGEN=OFF on a minimal single-tet overlap).
    {
        std::vector<bool> is_surface_vert(st.tetrahedralization.vertices.size(), false);
        for (const auto& tr : st.tetrahedralization.triangles)
        {
            is_surface_vert[tr[0]] = true;
            is_surface_vert[tr[1]] = true;
            is_surface_vert[tr[2]] = true;
        }

        std::set<Index> fallback_edges;
        for (const RegionID rid : st.object_region_info.regions_of_overlap)
        {
            if (tetgen_refined_regions.count(rid) > 0)
            {
                continue;
            }
            bool has_interior_v = false;
            for (const Index tid : st.region_info.region_to_tets_map[rid])
            {
                for (const Index cv : st.tetrahedralization.tetrahedra[tid])
                {
                    if (!is_surface_vert[cv])
                    {
                        has_interior_v = true;
                        break;
                    }
                }
                if (has_interior_v)
                {
                    break;
                }
            }
            if (!has_interior_v)
            {
                for (const Index tid : st.region_info.region_to_tets_map[rid])
                {
                    for (const Index eid : st.tetrahedralization.tet_to_edges_map[tid])
                    {
                        fallback_edges.insert(eid);
                    }
                }
            }
        }
        if (!fallback_edges.empty())
        {
            const std::vector<Index> to_split(fallback_edges.begin(), fallback_edges.end());
            split_edges(to_split, st);
        }
    }
}

struct GroupHarmonicField
{
    RegionID group_id = 0;
    std::vector<Index> tets_in_overlap;
    std::vector<Index> overlap_to_global_vertex_map;
    std::unordered_map<Index, Index> global_to_overlap_vertex_map;
    std::vector<RegionID> local_id_to_object_id_map;
    std::unordered_map<RegionID, RegionID> object_id_to_local_id_map;
    std::vector<Eigen::VectorXd> weights;
    std::vector<std::vector<RegionID>> vertex_label_sets;
    std::set<RegionBoundaryID> removable_boundary_ids;
};

GroupHarmonicField compute_group_harmonic_field(SegmentedTetrahedralization& sf,
                                                const OverlapRegionGroup& overlap_group)
{
    static std::atomic<int> overlap_count{ 0 };
    const int overlap_index = overlap_count.fetch_add(1, std::memory_order_relaxed);

    GroupHarmonicField out;
    out.group_id = overlap_group.group_id;

    const auto& overlap_regions_group = overlap_group.region_ids;
    const auto& internal_boundary_ids = overlap_group.internal_boundary_ids;
    const auto& external_boundary_ids = overlap_group.external_boundary_ids;

    for (const RegionID region_id : overlap_regions_group)
    {
        for (const Index tet_id : sf.region_info.region_to_tets_map[region_id])
        {
            out.tets_in_overlap.push_back(tet_id);
        }
    }

    std::set<Index> vertices_in_overlap_set;
    for (const Index tet_index : out.tets_in_overlap)
    {
        for (const Index n_id : sf.tetrahedralization.tetrahedra[tet_index])
        {
            vertices_in_overlap_set.insert(n_id);
        }
    }

    out.overlap_to_global_vertex_map.reserve(vertices_in_overlap_set.size());
    out.global_to_overlap_vertex_map.reserve(vertices_in_overlap_set.size());
    Index local_v_index = 0;
    for (const Index v_id : vertices_in_overlap_set)
    {
        out.overlap_to_global_vertex_map.push_back(v_id);
        out.global_to_overlap_vertex_map.emplace(v_id, local_v_index++);
    }

    const int n_verts = static_cast<int>(out.overlap_to_global_vertex_map.size());
    const int n_tets  = static_cast<int>(out.tets_in_overlap.size());

    Eigen::MatrixXd V_ov(n_verts, 3);
    for (int i = 0; i < n_verts; ++i)
    {
        const auto& p = sf.tetrahedralization.vertices[out.overlap_to_global_vertex_map[i]];
        V_ov(i, 0)    = p[0];
        V_ov(i, 1)    = p[1];
        V_ov(i, 2)    = p[2];
    }

    Eigen::MatrixXi T_ov(n_tets, 4);
    for (int t = 0; t < n_tets; ++t)
    {
        const auto& conn = sf.tetrahedralization.tetrahedra[out.tets_in_overlap[t]];
        for (int c = 0; c < 4; ++c)
        {
            T_ov(t, c) = static_cast<int>(out.global_to_overlap_vertex_map.at(conn[c]));
        }
    }

    std::set<RegionID> overlapping_objects;
    for (const RegionID region_id : overlap_regions_group)
    {
        for (const RegionID object_id : sf.object_region_info.region_to_objects_map[region_id])
        {
            overlapping_objects.insert(object_id);
        }
    }
    out.local_id_to_object_id_map.assign(overlapping_objects.begin(), overlapping_objects.end());
    for (RegionID k = 0; k < static_cast<RegionID>(out.local_id_to_object_id_map.size()); ++k)
    {
        out.object_id_to_local_id_map.emplace(out.local_id_to_object_id_map[k], k);
    }
    const int K = static_cast<int>(out.local_id_to_object_id_map.size());

    const std::set<RegionID> group_region_set(overlap_regions_group.begin(),
                                              overlap_regions_group.end());
    std::unordered_map<RegionID, std::vector<bool>> region_has_local_obj;
    for (const RegionID region_id : overlap_regions_group)
    {
        std::vector<bool> mask(K, false);
        for (const RegionID obj_id : sf.object_region_info.region_to_objects_map[region_id])
        {
            mask[out.object_id_to_local_id_map.at(obj_id)] = true;
        }
        region_has_local_obj.emplace(region_id, std::move(mask));
    }

    std::vector<std::vector<bool>> vertex_present_objects(n_verts, std::vector<bool>(K, false));
    for (const RegionID region_id : overlap_regions_group)
    {
        const auto& r_mask = region_has_local_obj.at(region_id);
        for (const Index tet_id : sf.region_info.region_to_tets_map[region_id])
        {
            for (const Index g_vid : sf.tetrahedralization.tetrahedra[tet_id])
            {
                const int lv = static_cast<int>(out.global_to_overlap_vertex_map.at(g_vid));
                for (int k = 0; k < K; ++k)
                {
                    if (r_mask[k])
                    {
                        vertex_present_objects[lv][k] = true;
                    }
                }
            }
        }
    }

    std::map<RegionBoundaryID, std::vector<int>> boundary_to_local_owners;
    for (const RegionBoundaryID b_id : internal_boundary_ids)
    {
        boundary_to_local_owners[b_id] = {};
    }
    for (const RegionBoundaryID b_id : external_boundary_ids)
    {
        boundary_to_local_owners[b_id] = {};
    }
    for (int k = 0; k < K; ++k)
    {
        const RegionID obj_id = out.local_id_to_object_id_map[k];
        for (const RegionBoundaryID b_id :
             sf.object_region_boundary_info.object_to_external_boundary_ids_map[obj_id])
        {
            auto it = boundary_to_local_owners.find(b_id);
            if (it != boundary_to_local_owners.end())
            {
                it->second.push_back(k);
            }
        }
    }
    for (const auto& kv : boundary_to_local_owners)
    {
        if (kv.second.size() == 1)
        {
            out.removable_boundary_ids.insert(kv.first);
        }
    }

    std::vector<bool> has_active_cap(n_verts, false);
    std::vector<std::vector<bool>> terminates_at_v(n_verts, std::vector<bool>(K, false));
    std::vector<std::vector<bool>> continues_at_v(n_verts, std::vector<bool>(K, false));

    for (const auto& kv : boundary_to_local_owners)
    {
        const RegionBoundaryID b_id = kv.first;
        const auto& owners_on_b     = kv.second;
        if (owners_on_b.empty())
        {
            continue;
        }

        std::vector<bool> is_owner_on_b(K, false);
        for (int k : owners_on_b)
        {
            is_owner_on_b[k] = true;
        }

        for (const Index tri_id : sf.region_boundaries_info.region_boundary_to_triangles_map[b_id])
        {
            const auto& br = sf.region_info.tri_to_bounded_regions_map[tri_id];
            std::vector<bool> c_in(K, false);
            if (br.first != NOT_SET_RegionID && group_region_set.count(br.first) > 0)
            {
                const auto& m = region_has_local_obj.at(br.first);
                for (int k = 0; k < K; ++k)
                {
                    c_in[k] = c_in[k] || m[k];
                }
            }
            if (br.second != NOT_SET_RegionID && group_region_set.count(br.second) > 0)
            {
                const auto& m = region_has_local_obj.at(br.second);
                for (int k = 0; k < K; ++k)
                {
                    c_in[k] = c_in[k] || m[k];
                }
            }

            std::vector<int> term_list;
            std::vector<int> surv_list;
            for (int k = 0; k < K; ++k)
            {
                if (c_in[k])
                {
                    if (is_owner_on_b[k])
                    {
                        term_list.push_back(k);
                    }
                    else
                    {
                        surv_list.push_back(k);
                    }
                }
            }

            if (!term_list.empty() && !surv_list.empty())
            {
                for (const Index g_vid : sf.tetrahedralization.triangles[tri_id])
                {
                    const auto it_lv = out.global_to_overlap_vertex_map.find(g_vid);
                    if (it_lv == out.global_to_overlap_vertex_map.end())
                    {
                        continue;
                    }
                    const int lv       = static_cast<int>(it_lv->second);
                    has_active_cap[lv] = true;
                    for (int k : term_list)
                    {
                        terminates_at_v[lv][k] = true;
                    }
                    for (int k : surv_list)
                    {
                        continues_at_v[lv][k] = true;
                    }
                }
            }
        }
    }

    // Exact element-by-element P1 Galerkin stiffness assembly:
    //   K^{(t)}_{ab} = int_{tet} (grad phi_a . grad phi_b) dV = (M_a . M_b) / (36 * V_tet)
    // Since M_0 + M_1 + M_2 + M_3 = 0 identically on every tetrahedron, each element matrix
    // K^{(t)} is symmetric positive semidefinite and reproduces affine functions to machine
    // precision without clamping obtuse cotangents or adding valence-biased graph Laplacians.
    std::vector<Eigen::Triplet<double>> a_triplets;
    a_triplets.reserve(static_cast<std::size_t>(n_tets) * 16 + static_cast<std::size_t>(n_verts));
    std::vector<std::vector<int>> adj(n_verts);
    std::vector<double> diag_accum(n_verts, 0.0);

    for (int t = 0; t < n_tets; ++t)
    {
        int vid[4] = { T_ov(t, 0), T_ov(t, 1), T_ov(t, 2), T_ov(t, 3) };
        for (int e = 0; e < 6; ++e)
        {
            const int u = vid[detail::kTetEdges[e][0]];
            const int v = vid[detail::kTetEdges[e][1]];
            adj[u].push_back(v);
            adj[v].push_back(u);
        }

        Eigen::Vector3d q[4] = { V_ov.row(vid[0]).transpose(),
                                 V_ov.row(vid[1]).transpose(),
                                 V_ov.row(vid[2]).transpose(),
                                 V_ov.row(vid[3]).transpose() };

        double sum_edge_sq = 0.0;
        for (int e = 0; e < 6; ++e)
        {
            sum_edge_sq +=
              (q[detail::kTetEdges[e][1]] - q[detail::kTetEdges[e][0]]).squaredNorm();
        }
        const double l_rms_sq = sum_edge_sq / 6.0;
        if (l_rms_sq <= 1e-30)
        {
            continue;
        }
        const double l_rms_cubed = l_rms_sq * std::sqrt(l_rms_sq);

        double six_v = (q[1] - q[0]).dot((q[2] - q[0]).cross(q[3] - q[0]));
        if (six_v < 0.0)
        {
            std::swap(vid[1], vid[2]);
            std::swap(q[1], q[2]);
            six_v = -six_v;
        }
        const double six_v_eff = std::max(six_v, 1e-6 * l_rms_cubed);

        const Eigen::Vector3d M[4] = { (q[2] - q[1]).cross(q[3] - q[1]),
                                       (q[3] - q[0]).cross(q[2] - q[0]),
                                       (q[1] - q[0]).cross(q[3] - q[0]),
                                       (q[2] - q[0]).cross(q[1] - q[0]) };
        const double inv_36_vol    = 1.0 / (6.0 * six_v_eff);

        for (int a = 0; a < 4; ++a)
        {
            const double k_aa = M[a].squaredNorm() * inv_36_vol;
            if (!std::isfinite(k_aa))
            {
                continue;
            }
            a_triplets.emplace_back(vid[a], vid[a], k_aa);
            diag_accum[vid[a]] += k_aa;
            for (int b = a + 1; b < 4; ++b)
            {
                const double k_ab = M[a].dot(M[b]) * inv_36_vol;
                if (std::isfinite(k_ab))
                {
                    a_triplets.emplace_back(vid[a], vid[b], k_ab);
                    a_triplets.emplace_back(vid[b], vid[a], k_ab);
                }
            }
        }
    }

    for (int v = 0; v < n_verts; ++v)
    {
        auto& nb = adj[v];
        std::sort(nb.begin(), nb.end());
        nb.erase(std::unique(nb.begin(), nb.end()), nb.end());
    }

    double diag_scale = 1.0;
    {
        std::vector<double> pos_diags;
        pos_diags.reserve(n_verts);
        for (const double d : diag_accum)
        {
            if (d > 0.0 && std::isfinite(d))
            {
                pos_diags.push_back(d);
            }
        }
        if (!pos_diags.empty())
        {
            const std::size_t mid = pos_diags.size() / 2;
            std::nth_element(pos_diags.begin(),
                             pos_diags.begin() + static_cast<std::ptrdiff_t>(mid),
                             pos_diags.end());
            diag_scale = std::max(1e-14, pos_diags[mid]);
        }
    }

    Eigen::SparseMatrix<double> A(n_verts, n_verts);
    A.setFromTriplets(a_triplets.begin(), a_triplets.end());

    std::vector<std::vector<bool>> is_dirichlet(K, std::vector<bool>(n_verts, false));
    std::vector<Eigen::VectorXd> bc_values(K, Eigen::VectorXd::Zero(n_verts));
    std::vector<std::vector<int>> seam_tied_local_objs(n_verts);

    for (int v = 0; v < n_verts; ++v)
    {
        for (int k = 0; k < K; ++k)
        {
            if (!vertex_present_objects[v][k])
            {
                is_dirichlet[k][v] = true;
                bc_values[k](v)    = 0.0;
            }
        }

        if (!has_active_cap[v])
        {
            continue;
        }

        std::vector<int> strict_survivors;
        for (int k = 0; k < K; ++k)
        {
            if (vertex_present_objects[v][k] && continues_at_v[v][k] && !terminates_at_v[v][k])
            {
                strict_survivors.push_back(k);
            }
        }

        if (strict_survivors.size() == 1)
        {
            const int winner = strict_survivors[0];
            for (int k = 0; k < K; ++k)
            {
                is_dirichlet[k][v] = true;
                bc_values[k](v)    = (k == winner) ? 1.0 : 0.0;
            }
        }
        else if (strict_survivors.size() >= 2)
        {
            std::vector<bool> is_surv(K, false);
            for (int k : strict_survivors)
            {
                is_surv[k] = true;
            }
            for (int k = 0; k < K; ++k)
            {
                if (!is_surv[k])
                {
                    is_dirichlet[k][v] = true;
                    bc_values[k](v)    = 0.0;
                }
            }
        }
        else
        {
            std::vector<int> tied_seam;
            for (int k = 0; k < K; ++k)
            {
                if (vertex_present_objects[v][k] && (continues_at_v[v][k] || terminates_at_v[v][k]))
                {
                    tied_seam.push_back(k);
                }
            }
            if (!tied_seam.empty())
            {
                const double eq_val = 1.0 / static_cast<double>(tied_seam.size());
                std::vector<bool> in_seam(K, false);
                for (int k : tied_seam)
                {
                    in_seam[k]         = true;
                    is_dirichlet[k][v] = true;
                    bc_values[k](v)    = eq_val;
                }
                for (int k = 0; k < K; ++k)
                {
                    if (!in_seam[k])
                    {
                        is_dirichlet[k][v] = true;
                        bc_values[k](v)    = 0.0;
                    }
                }
                seam_tied_local_objs[v] = std::move(tied_seam);
            }
        }
    }

    std::vector<Eigen::VectorXd> raw_u(K, Eigen::VectorXd::Zero(n_verts));
    for (int k = 0; k < K; ++k)
    {
        std::vector<bool> visited(n_verts, false);
        for (int start = 0; start < n_verts; ++start)
        {
            if (visited[start] || is_dirichlet[k][start])
            {
                continue;
            }
            std::vector<int> comp;
            bool touches_dirichlet = false;
            std::deque<int> q      = { start };
            visited[start]         = true;
            while (!q.empty())
            {
                const int cur = q.front();
                q.pop_front();
                comp.push_back(cur);
                for (int nb : adj[cur])
                {
                    if (is_dirichlet[k][nb])
                    {
                        touches_dirichlet = true;
                    }
                    else if (!visited[nb])
                    {
                        visited[nb] = true;
                        q.push_back(nb);
                    }
                }
            }
            if (!touches_dirichlet)
            {
                for (int v : comp)
                {
                    int first_present = 0;
                    while (first_present < K && !vertex_present_objects[v][first_present])
                    {
                        ++first_present;
                    }
                    is_dirichlet[k][v] = true;
                    bc_values[k](v)    = (k == first_present) ? 1.0 : 0.0;
                }
            }
        }

        std::vector<int> free_to_vert;
        std::vector<int> vert_to_free(n_verts, -1);
        for (int v = 0; v < n_verts; ++v)
        {
            if (!is_dirichlet[k][v])
            {
                vert_to_free[v] = static_cast<int>(free_to_vert.size());
                free_to_vert.push_back(v);
            }
            else
            {
                raw_u[k](v) = bc_values[k](v);
            }
        }

        const int n_free = static_cast<int>(free_to_vert.size());
        if (n_free == 0)
        {
            continue;
        }

        std::vector<Eigen::Triplet<double>> sub_trips;
        sub_trips.reserve(static_cast<std::size_t>(A.nonZeros()));
        Eigen::VectorXd rhs = Eigen::VectorXd::Zero(n_free);

        for (int col = 0; col < A.outerSize(); ++col)
        {
            const int fc = vert_to_free[col];
            for (Eigen::SparseMatrix<double>::InnerIterator it(A, col); it; ++it)
            {
                const int row = static_cast<int>(it.row());
                const int fr  = vert_to_free[row];
                if (fr >= 0)
                {
                    if (fc >= 0)
                    {
                        sub_trips.emplace_back(fr, fc, it.value());
                    }
                    else
                    {
                        rhs(fr) -= it.value() * bc_values[k](col);
                    }
                }
            }
        }

        if (rhs.squaredNorm() <= 1e-30)
        {
            for (int f = 0; f < n_free; ++f)
            {
                raw_u[k](free_to_vert[f]) = 0.0;
            }
            continue;
        }

        for (int f = 0; f < n_free; ++f)
        {
            sub_trips.emplace_back(f, f, 1e-11 * diag_scale);
        }
        Eigen::SparseMatrix<double> A_in(n_free, n_free);
        A_in.setFromTriplets(sub_trips.begin(), sub_trips.end());

        Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> ldlt;
        ldlt.compute(A_in);
        Eigen::VectorXd u_in = Eigen::VectorXd::Zero(n_free);
        bool solved          = false;
        if (ldlt.info() == Eigen::Success)
        {
            u_in   = ldlt.solve(rhs);
            solved = (ldlt.info() == Eigen::Success && u_in.allFinite());
        }
        if (!solved)
        {
            Eigen::SparseLU<Eigen::SparseMatrix<double>> lu;
            lu.compute(A_in);
            if (lu.info() == Eigen::Success)
            {
                u_in   = lu.solve(rhs);
                solved = (lu.info() == Eigen::Success && u_in.allFinite());
            }
        }
        if (!solved)
        {
            Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower | Eigen::Upper> cg;
            cg.compute(A_in);
            if (cg.info() == Eigen::Success)
            {
                u_in = cg.solve(rhs);
            }
        }

        for (int f = 0; f < n_free; ++f)
        {
            raw_u[k](free_to_vert[f]) = std::clamp(u_in(f), 0.0, 1.0);
        }
    }

    out.weights.assign(K, Eigen::VectorXd::Zero(n_verts));
    for (int v = 0; v < n_verts; ++v)
    {
        double sum_u = 0.0;
        for (int k = 0; k < K; ++k)
        {
            if (vertex_present_objects[v][k])
            {
                sum_u += std::max(0.0, raw_u[k](v));
            }
        }
        if (sum_u > 1e-14)
        {
            for (int k = 0; k < K; ++k)
            {
                if (vertex_present_objects[v][k])
                {
                    out.weights[k](v) = std::max(0.0, raw_u[k](v)) / sum_u;
                }
            }
        }
        else
        {
            for (int k = 0; k < K; ++k)
            {
                if (vertex_present_objects[v][k])
                {
                    out.weights[k](v) = 1.0;
                    break;
                }
            }
        }
    }

    out.vertex_label_sets.resize(n_verts);
    auto pick_single_winner = [&](int v) -> RegionID
    {
        int best_k      = -1;
        double best_val = -1.0;
        for (int k = 0; k < K; ++k)
        {
            if (!vertex_present_objects[v][k])
            {
                continue;
            }
            const double val = out.weights[k](v);
            if (best_k < 0 || val > best_val + 1e-14)
            {
                best_val = val;
                best_k   = k;
            }
        }
        return out.local_id_to_object_id_map[std::max(0, best_k)];
    };

    // Identify strictly-interior overlap vertices (not on any tracked surface triangle)
    std::vector<bool> is_tracked_surf_local(n_verts, false);
    for (const auto& tri : sf.tetrahedralization.triangles)
    {
        for (const Index gv : tri)
        {
            const auto it = out.global_to_overlap_vertex_map.find(gv);
            if (it != out.global_to_overlap_vertex_map.end())
            {
                is_tracked_surf_local[static_cast<int>(it->second)] = true;
            }
        }
    }

    std::vector<bool> snapped_interior(n_verts, false);
    std::vector<double> snap_gap(n_verts, 1.0);
    for (int v = 0; v < n_verts; ++v)
    {
        if (!seam_tied_local_objs[v].empty())
        {
            for (int k : seam_tied_local_objs[v])
            {
                out.vertex_label_sets[v].push_back(out.local_id_to_object_id_map[k]);
            }
            std::sort(out.vertex_label_sets[v].begin(), out.vertex_label_sets[v].end());
        }
        else
        {
            bool did_snap = false;
            if (!is_tracked_surf_local[v] && !has_active_cap[v] && K >= 2)
            {
                int k1 = -1, k2 = -1;
                double w1 = -1.0, w2 = -1.0, w3 = -1.0;
                for (int k = 0; k < K; ++k)
                {
                    if (!vertex_present_objects[v][k])
                    {
                        continue;
                    }
                    const double wk = out.weights[k](v);
                    if (wk > w1)
                    {
                        w3 = w2;
                        w2 = w1;
                        k2 = k1;
                        w1 = wk;
                        k1 = k;
                    }
                    else if (wk > w2)
                    {
                        w3 = w2;
                        w2 = wk;
                        k2 = k;
                    }
                    else if (wk > w3)
                    {
                        w3 = wk;
                    }
                }
                if (k1 >= 0 && k2 >= 0 && w1 > 0.35 && (w1 - w2) < 0.014 && (w2 - w3) > 0.06)
                {
                    out.vertex_label_sets[v] = { out.local_id_to_object_id_map[k1],
                                                 out.local_id_to_object_id_map[k2] };
                    std::sort(out.vertex_label_sets[v].begin(), out.vertex_label_sets[v].end());
                    snapped_interior[v] = true;
                    snap_gap[v]         = w1 - w2;
                    did_snap            = true;
                }
            }
            if (!did_snap)
            {
                out.vertex_label_sets[v] = { pick_single_winner(v) };
            }
        }
    }

    // Prevent >= 3 non-seam snapped vertices on a single tet (avoids degenerate all-tied tets)
    for (int t = 0; t < n_tets; ++t)
    {
        std::vector<int> snap_vs;
        for (int c = 0; c < 4; ++c)
        {
            const int lv = T_ov(t, c);
            if (snapped_interior[lv])
            {
                snap_vs.push_back(lv);
            }
        }
        while (snap_vs.size() > 2)
        {
            auto worst_it = std::max_element(
              snap_vs.begin(),
              snap_vs.end(),
              [&](int a, int b) { return snap_gap[a] < snap_gap[b]; });
            const int drop_v             = *worst_it;
            snapped_interior[drop_v]     = false;
            out.vertex_label_sets[drop_v] = { pick_single_winner(drop_v) };
            snap_vs.erase(worst_it);
        }
    }

    static constexpr int kLocalTetFaces[4][3] = {
        { 0, 1, 2 }, { 0, 1, 3 }, { 0, 2, 3 }, { 1, 2, 3 }
    };
    auto sets_intersect = [](const std::vector<RegionID>& a, const std::vector<RegionID>& b)
    {
        std::size_t i = 0, j = 0;
        while (i < a.size() && j < b.size())
        {
            if (a[i] == b[j])
            {
                return true;
            }
            if (a[i] < b[j])
            {
                ++i;
            }
            else
            {
                ++j;
            }
        }
        return false;
    };
    auto sets_intersect3 = [](const std::vector<RegionID>& a,
                              const std::vector<RegionID>& b,
                              const std::vector<RegionID>& c)
    {
        for (const RegionID x : a)
        {
            if (std::binary_search(b.begin(), b.end(), x) &&
                std::binary_search(c.begin(), c.end(), x))
            {
                return true;
            }
        }
        return false;
    };

    for (int t = 0; t < n_tets; ++t)
    {
        for (int f = 0; f < 4; ++f)
        {
            const int u    = T_ov(t, kLocalTetFaces[f][0]);
            const int v    = T_ov(t, kLocalTetFaces[f][1]);
            const int w    = T_ov(t, kLocalTetFaces[f][2]);
            const auto& lu = out.vertex_label_sets[u];
            const auto& lv = out.vertex_label_sets[v];
            const auto& lw = out.vertex_label_sets[w];
            if (sets_intersect(lu, lv) && sets_intersect(lv, lw) && sets_intersect(lu, lw) &&
                !sets_intersect3(lu, lv, lw))
            {
                snapped_interior[u]      = false;
                snapped_interior[v]      = false;
                snapped_interior[w]      = false;
                out.vertex_label_sets[u] = { pick_single_winner(u) };
                out.vertex_label_sets[v] = { pick_single_winner(v) };
                out.vertex_label_sets[w] = { pick_single_winner(w) };
            }
        }
    }

    // Gently project snapped interior vertices onto the local P1 zero-isosurface (w_a - w_b = 0)
    // whenever all incident tetrahedra maintain >= 40% of their positive volume.
    {
        std::vector<std::vector<Index>> snapped_inc_tets(n_verts);
        for (int t = 0; t < n_tets; ++t)
        {
            const Index g_tid = out.tets_in_overlap[static_cast<std::size_t>(t)];
            for (int c = 0; c < 4; ++c)
            {
                const int lv = T_ov(t, c);
                if (snapped_interior[lv])
                {
                    snapped_inc_tets[lv].push_back(g_tid);
                }
            }
        }

        for (int lv = 0; lv < n_verts; ++lv)
        {
            if (!snapped_interior[lv] || out.vertex_label_sets[lv].size() != 2 ||
                snapped_inc_tets[lv].empty())
            {
                continue;
            }
            const int ka     = static_cast<int>(out.object_id_to_local_id_map.at(out.vertex_label_sets[lv][0]));
            const int kb     = static_cast<int>(out.object_id_to_local_id_map.at(out.vertex_label_sets[lv][1]));
            const double d_v = out.weights[ka](lv) - out.weights[kb](lv);
            if (std::abs(d_v) <= 1e-12)
            {
                continue;
            }

            Eigen::Vector3d grad_sum = Eigen::Vector3d::Zero();
            double vol_sum           = 0.0;
            double max_edge_len      = 0.0;
            double min_vol_before    = std::numeric_limits<double>::max();
            const Index gv           = out.overlap_to_global_vertex_map[static_cast<std::size_t>(lv)];

            for (const Index g_tid : snapped_inc_tets[lv])
            {
                const auto& tc = sf.tetrahedralization.tetrahedra[g_tid];
                const Eigen::Vector3d q0 =
                  detail::point_as_vector(sf.tetrahedralization.vertices[tc[0]]);
                const Eigen::Vector3d q1 =
                  detail::point_as_vector(sf.tetrahedralization.vertices[tc[1]]);
                const Eigen::Vector3d q2 =
                  detail::point_as_vector(sf.tetrahedralization.vertices[tc[2]]);
                const Eigen::Vector3d q3 =
                  detail::point_as_vector(sf.tetrahedralization.vertices[tc[3]]);
                const double v_tet = detail::tetrahedron_signed_volume(q0, q1, q2, q3);
                min_vol_before     = std::min(min_vol_before, v_tet);
                if (v_tet <= 1e-22)
                {
                    continue;
                }

                Eigen::Matrix3d em;
                em.col(0) = q1 - q0;
                em.col(1) = q2 - q0;
                em.col(2) = q3 - q0;
                max_edge_len =
                  std::max({ max_edge_len, em.col(0).norm(), em.col(1).norm(), em.col(2).norm() });

                double d_c[4];
                for (int c = 0; c < 4; ++c)
                {
                    const int clv = static_cast<int>(out.global_to_overlap_vertex_map.at(tc[c]));
                    d_c[c]        = out.weights[ka](clv) - out.weights[kb](clv);
                }
                const Eigen::Vector3d rhs_d(d_c[1] - d_c[0], d_c[2] - d_c[0], d_c[3] - d_c[0]);
                const Eigen::Vector3d g_tet = em.transpose().colPivHouseholderQr().solve(rhs_d);
                if (g_tet.allFinite())
                {
                    grad_sum += v_tet * g_tet;
                    vol_sum += v_tet;
                }
            }

            if (vol_sum <= 1e-22 || min_vol_before <= 1e-22)
            {
                continue;
            }
            const Eigen::Vector3d grad_d = grad_sum / vol_sum;
            const double g_norm2         = grad_d.squaredNorm();
            if (g_norm2 <= 1e-20)
            {
                continue;
            }

            Eigen::Vector3d step = -(d_v / g_norm2) * grad_d;
            const double step_n  = step.norm();
            if (step_n > 0.20 * max_edge_len)
            {
                step *= (0.20 * max_edge_len / step_n);
            }

            const Point3 old_p = sf.tetrahedralization.vertices[gv];
            sf.tetrahedralization.vertices[gv] = { old_p[0] + step.x(),
                                                   old_p[1] + step.y(),
                                                   old_p[2] + step.z() };
            bool ok = true;
            for (const Index g_tid : snapped_inc_tets[lv])
            {
                const auto& tc  = sf.tetrahedralization.tetrahedra[g_tid];
                const double nv = detail::tetrahedron_signed_volume(
                  sf.tetrahedralization.vertices[tc[0]],
                  sf.tetrahedralization.vertices[tc[1]],
                  sf.tetrahedralization.vertices[tc[2]],
                  sf.tetrahedralization.vertices[tc[3]]);
                if (nv < 0.45 * min_vol_before || nv <= 1e-22)
                {
                    ok = false;
                    break;
                }
            }
            if (!ok)
            {
                sf.tetrahedralization.vertices[gv] = old_p;
            }
            else
            {
                V_ov.row(lv) =
                  detail::point_as_vector(sf.tetrahedralization.vertices[gv]).transpose();
            }
        }
    }

    if (GlobalOptions::get().write_debug_meshes)
    {
        std::vector<Point3> v_pts(n_verts);
        for (int i = 0; i < n_verts; ++i)
        {
            v_pts[i] = { V_ov(i, 0), V_ov(i, 1), V_ov(i, 2) };
        }
        std::vector<std::array<Index, 4>> t_els(n_tets);
        for (int t = 0; t < n_tets; ++t)
        {
            t_els[t] = { static_cast<Index>(T_ov(t, 0)),
                         static_cast<Index>(T_ov(t, 1)),
                         static_cast<Index>(T_ov(t, 2)),
                         static_cast<Index>(T_ov(t, 3)) };
        }
        std::vector<detail::VtkScalarFieldF64> pt_fields;
        pt_fields.reserve(K);
        for (int k = 0; k < K; ++k)
        {
            std::vector<double> vals(out.weights[k].data(),
                                     out.weights[k].data() + out.weights[k].size());
            pt_fields.push_back({ "harmonic_weight_" + std::to_string(k), std::move(vals) });
        }
        detail::write_vtk_tetrahedra(GlobalOptions::get().output_dir + "/harmonic_weights_" +
                                       std::to_string(overlap_index) + ".vtk",
                                     v_pts,
                                     t_els,
                                     {},
                                     pt_fields);
    }

    return out;
}

void subdivide_face_canonically(
  const std::array<Index, 3>& orig_face,
  const std::unordered_map<std::size_t, Index>& edge_cut_map,
  std::unordered_map<std::array<Index, 3>, Index, SortedTriKeyHash>& face_steiner_map,
  std::vector<Point3>& vertices,
  std::vector<std::array<Index, 3>>& out_sub_tris)
{
    std::array<Index, 3> s = orig_face;
    bool flipped           = false;
    if (s[0] > s[1])
    {
        std::swap(s[0], s[1]);
        flipped = !flipped;
    }
    if (s[1] > s[2])
    {
        std::swap(s[1], s[2]);
        flipped = !flipped;
    }
    if (s[0] > s[1])
    {
        std::swap(s[0], s[1]);
        flipped = !flipped;
    }

    auto find_cut = [&](Index a, Index b) -> Index
    {
        const auto it = edge_cut_map.find(hash_linear_edge(a, b));
        return (it != edge_cut_map.end()) ? it->second : NOT_SET_Index;
    };

    const Index p01 = find_cut(s[0], s[1]);
    const Index p12 = find_cut(s[1], s[2]);
    const Index p02 = find_cut(s[0], s[2]);

    const int n_cuts = (p01 != NOT_SET_Index ? 1 : 0) + (p12 != NOT_SET_Index ? 1 : 0) +
                       (p02 != NOT_SET_Index ? 1 : 0);

    auto emit = [&](Index a, Index b, Index c)
    {
        if (flipped)
        {
            out_sub_tris.push_back({ a, c, b });
        }
        else
        {
            out_sub_tris.push_back({ a, b, c });
        }
    };

    if (n_cuts == 0)
    {
        out_sub_tris.push_back(orig_face);
        return;
    }

    if (n_cuts == 1)
    {
        if (p01 != NOT_SET_Index)
        {
            emit(s[0], p01, s[2]);
            emit(p01, s[1], s[2]);
        }
        else if (p12 != NOT_SET_Index)
        {
            emit(s[1], p12, s[0]);
            emit(p12, s[2], s[0]);
        }
        else
        {
            emit(s[0], s[1], p02);
            emit(p02, s[1], s[2]);
        }
        return;
    }

    if (n_cuts == 2)
    {
        Index c         = NOT_SET_Index;
        Index a         = NOT_SET_Index;
        Index b         = NOT_SET_Index;
        Index pca       = NOT_SET_Index;
        Index pcb       = NOT_SET_Index;
        bool local_flip = false;

        if (p12 == NOT_SET_Index)
        {
            c          = s[0];
            a          = s[1];
            b          = s[2];
            pca        = p01;
            pcb        = p02;
            local_flip = false;
        }
        else if (p02 == NOT_SET_Index)
        {
            c          = s[1];
            a          = s[0];
            b          = s[2];
            pca        = p01;
            pcb        = p12;
            local_flip = true;
        }
        else
        {
            c          = s[2];
            a          = s[0];
            b          = s[1];
            pca        = p02;
            pcb        = p12;
            local_flip = false;
        }

        auto emit_oriented = [&](Index x, Index y, Index z)
        {
            if (local_flip)
            {
                emit(x, z, y);
            }
            else
            {
                emit(x, y, z);
            }
        };

        emit_oriented(c, pca, pcb);
        emit_oriented(a, b, pcb);
        emit_oriented(a, pcb, pca);
        return;
    }

    Index pf        = NOT_SET_Index;
    const auto it_f = face_steiner_map.find(s);
    if (it_f != face_steiner_map.end())
    {
        pf = it_f->second;
    }
    else
    {
        pf             = static_cast<Index>(vertices.size());
        const auto& a0 = vertices[p01];
        const auto& a1 = vertices[p12];
        const auto& a2 = vertices[p02];
        vertices.push_back(
          { (a0[0] + a1[0] + a2[0]) / 3.0,
            (a0[1] + a1[1] + a2[1]) / 3.0,
            (a0[2] + a1[2] + a2[2]) / 3.0 });
        face_steiner_map.emplace(s, pf);
    }

    emit(s[0], p01, pf);
    emit(s[0], pf, p02);
    emit(s[1], p12, pf);
    emit(s[1], pf, p01);
    emit(s[2], p02, pf);
    emit(s[2], pf, p12);
}

RegionID find_max_priority_object_for_region(const std::vector<RegionID>& region_objects,
                                             const std::vector<PriorityType>& object_priorities)
{
    PriorityType best_p = std::numeric_limits<PriorityType>::min();
    RegionID best_obj   = NOT_SET_RegionID;
    for (const RegionID obj_id : region_objects)
    {
        const PriorityType p = object_priorities[obj_id];
        if (best_obj == NOT_SET_RegionID || p > best_p || (p == best_p && obj_id < best_obj))
        {
            best_p   = p;
            best_obj = obj_id;
        }
    }
    return best_obj;
}

void subdivide_and_partition_overlaps(const CombinedSurfaces& inconsistent_combined_surfaces,
                                      const ConsistentSurfaces& consistent_surfaces,
                                      const std::vector<OverlapRegionGroup>& /*overlap_groups*/,
                                      const std::vector<GroupHarmonicField>& group_fields,
                                      SegmentedTetrahedralization& st)
{
    const std::size_t n_orig_verts = st.tetrahedralization.vertices.size();
    const std::size_t n_orig_tets  = st.tetrahedralization.tetrahedra.size();
    const std::size_t n_orig_tris  = st.tetrahedralization.triangles.size();

    std::vector<int> tet_to_group_idx(n_orig_tets, -1);
    std::set<RegionBoundaryID> all_removable_boundaries;
    for (std::size_t g = 0; g < group_fields.size(); ++g)
    {
        for (const Index tet_id : group_fields[g].tets_in_overlap)
        {
            tet_to_group_idx[tet_id] = static_cast<int>(g);
        }
        all_removable_boundaries.insert(group_fields[g].removable_boundary_ids.begin(),
                                        group_fields[g].removable_boundary_ids.end());
    }

    struct EdgeCutAccum
    {
        Index v_min  = 0;
        Index v_max  = 0;
        double t_sum = 0.0;
        int count    = 0;
    };
    std::unordered_map<std::size_t, EdgeCutAccum> cut_accum_map;

    auto label_sets_disjoint = [](const std::vector<RegionID>& a, const std::vector<RegionID>& b)
    {
        std::size_t i = 0, j = 0;
        while (i < a.size() && j < b.size())
        {
            if (a[i] == b[j])
            {
                return false;
            }
            if (a[i] < b[j])
            {
                ++i;
            }
            else
            {
                ++j;
            }
        }
        return true;
    };

    for (const auto& gf : group_fields)
    {
        std::unordered_set<std::size_t> visited_group_edges;
        for (const Index tet_id : gf.tets_in_overlap)
        {
            const auto& conn = st.tetrahedralization.tetrahedra[tet_id];
            for (int e = 0; e < 6; ++e)
            {
                Index ga = conn[detail::kTetEdges[e][0]];
                Index gb = conn[detail::kTetEdges[e][1]];
                if (ga > gb)
                {
                    std::swap(ga, gb);
                }
                const std::size_t h = hash_linear_edge(ga, gb);
                if (!visited_group_edges.insert(h).second)
                {
                    continue;
                }

                const int la      = static_cast<int>(gf.global_to_overlap_vertex_map.at(ga));
                const int lb      = static_cast<int>(gf.global_to_overlap_vertex_map.at(gb));
                const auto& set_a = gf.vertex_label_sets[la];
                const auto& set_b = gf.vertex_label_sets[lb];
                if (!label_sets_disjoint(set_a, set_b))
                {
                    continue;
                }

                auto eval_max_group_w = [&](const std::vector<RegionID>& s, int lv) -> double
                {
                    double m = 0.0;
                    for (const RegionID obj : s)
                    {
                        const RegionID k = gf.object_id_to_local_id_map.at(obj);
                        m                = std::max(m, gf.weights[k](lv));
                    }
                    return m;
                };

                const double da    = eval_max_group_w(set_a, la) - eval_max_group_w(set_b, la);
                const double db    = eval_max_group_w(set_a, lb) - eval_max_group_w(set_b, lb);
                const double denom = da - db;
                const double t =
                  (std::abs(denom) > 1e-14) ? std::clamp(da / denom, 0.002, 0.998) : 0.5;

                auto& acc = cut_accum_map[h];
                acc.v_min = ga;
                acc.v_max = gb;
                acc.t_sum += t;
                acc.count += 1;
            }
        }
    }

    std::vector<EdgeCutAccum> sorted_cuts;
    sorted_cuts.reserve(cut_accum_map.size());
    for (const auto& kv : cut_accum_map)
    {
        sorted_cuts.push_back(kv.second);
    }
    std::sort(sorted_cuts.begin(),
              sorted_cuts.end(),
              [](const EdgeCutAccum& a, const EdgeCutAccum& b)
              { return (a.v_min != b.v_min) ? (a.v_min < b.v_min) : (a.v_max < b.v_max); });

    std::unordered_map<std::size_t, Index> edge_cut_map;
    edge_cut_map.reserve(sorted_cuts.size());
    struct EdgeSteinerOrigin
    {
        Index v0 = 0;
        Index v1 = 0;
        double t = 0.5;
    };
    std::unordered_map<Index, EdgeSteinerOrigin> edge_steiner_origin;
    edge_steiner_origin.reserve(sorted_cuts.size());

    for (const auto& c : sorted_cuts)
    {
        const double t   = c.t_sum / static_cast<double>(c.count);
        const auto& p0   = st.tetrahedralization.vertices[c.v_min];
        const auto& p1   = st.tetrahedralization.vertices[c.v_max];
        const Index vid  = static_cast<Index>(st.tetrahedralization.vertices.size());
        st.tetrahedralization.vertices.push_back({ (1.0 - t) * p0[0] + t * p1[0],
                                                   (1.0 - t) * p0[1] + t * p1[1],
                                                   (1.0 - t) * p0[2] + t * p1[2] });
        edge_cut_map.emplace(hash_linear_edge(c.v_min, c.v_max), vid);
        edge_steiner_origin.emplace(vid, EdgeSteinerOrigin{ c.v_min, c.v_max, t });
    }

    std::unordered_map<std::array<Index, 3>, Index, SortedTriKeyHash> face_steiner_map;

    std::vector<RegionID> orig_region_owner(st.region_info.n_regions, NOT_SET_RegionID);
    for (RegionID r = 0; r < st.region_info.n_regions; ++r)
    {
        orig_region_owner[r] = find_max_priority_object_for_region(
          st.object_region_info.region_to_objects_map[r],
          inconsistent_combined_surfaces.object_priorities);
    }

    auto eval_weight_at_vert = [&](const GroupHarmonicField& gf,
                                   int k,
                                   Index vid,
                                   const std::array<double, 4>& tet_corner_w,
                                   const Eigen::Vector3d& p0,
                                   const Eigen::Matrix3d& inv_edge_mat,
                                   bool has_inv_edge_mat) -> double
    {
        if (vid < n_orig_verts)
        {
            const auto it = gf.global_to_overlap_vertex_map.find(vid);
            return (it != gf.global_to_overlap_vertex_map.end()) ? gf.weights[k](it->second) : 0.0;
        }
        if (has_inv_edge_mat)
        {
            const Eigen::Vector3d pv   = detail::point_as_vector(st.tetrahedralization.vertices[vid]);
            const Eigen::Vector3d bary = inv_edge_mat * (pv - p0);
            const double l0            = 1.0 - bary.x() - bary.y() - bary.z();
            return l0 * tet_corner_w[0] + bary.x() * tet_corner_w[1] +
                   bary.y() * tet_corner_w[2] + bary.z() * tet_corner_w[3];
        }
        const auto it_e = edge_steiner_origin.find(vid);
        if (it_e != edge_steiner_origin.end())
        {
            const auto it0  = gf.global_to_overlap_vertex_map.find(it_e->second.v0);
            const auto it1  = gf.global_to_overlap_vertex_map.find(it_e->second.v1);
            const double w0 =
              (it0 != gf.global_to_overlap_vertex_map.end()) ? gf.weights[k](it0->second) : 0.0;
            const double w1 =
              (it1 != gf.global_to_overlap_vertex_map.end()) ? gf.weights[k](it1->second) : 0.0;
            return (1.0 - it_e->second.t) * w0 + it_e->second.t * w1;
        }
        return 0.25 * (tet_corner_w[0] + tet_corner_w[1] + tet_corner_w[2] + tet_corner_w[3]);
    };

    struct SubTetMeta
    {
        bool is_overlap      = false;
        RegionID group_id    = 0;
        RegionID owner       = NOT_SET_RegionID;
        RegionID orig_region = 0;
    };

    std::vector<std::array<Index, 4>> new_tetrahedra;
    std::vector<SubTetMeta> new_tet_meta;
    new_tetrahedra.reserve(n_orig_tets * 2);
    new_tet_meta.reserve(n_orig_tets * 2);

    for (std::size_t tid = 0; tid < n_orig_tets; ++tid)
    {
        const auto& conn      = st.tetrahedralization.tetrahedra[tid];
        const RegionID orig_r = st.region_info.tet_to_region_map[tid];
        const int g_idx       = tet_to_group_idx[tid];
        const bool is_ov      = (g_idx >= 0);

        std::vector<std::pair<int, Index>> tet_cuts;
        for (int e = 0; e < 6; ++e)
        {
            const Index va = conn[detail::kTetEdges[e][0]];
            const Index vb = conn[detail::kTetEdges[e][1]];
            const auto it  = edge_cut_map.find(hash_linear_edge(va, vb));
            if (it != edge_cut_map.end())
            {
                tet_cuts.emplace_back(e, it->second);
            }
        }

        const std::size_t sub_begin = new_tetrahedra.size();
        Index tet_steiner_vid       = NOT_SET_Index;

        if (tet_cuts.empty())
        {
            new_tetrahedra.push_back(conn);
        }
        else if (tet_cuts.size() == 1)
        {
            const int e    = tet_cuts[0].first;
            const Index pe = tet_cuts[0].second;
            const int ia   = detail::kTetEdges[e][0];
            const int ib   = detail::kTetEdges[e][1];
            int ic = -1, id = -1;
            for (int c = 0; c < 4; ++c)
            {
                if (c != ia && c != ib)
                {
                    if (ic < 0)
                    {
                        ic = c;
                    }
                    else
                    {
                        id = c;
                    }
                }
            }
            new_tetrahedra.push_back(
              make_positive_tet(st.tetrahedralization.vertices, conn[ia], pe, conn[ic], conn[id]));
            new_tetrahedra.push_back(
              make_positive_tet(st.tetrahedralization.vertices, pe, conn[ib], conn[ic], conn[id]));
        }
        else
        {
            int single_face_apex = -1;
            for (int apex_c = 0; apex_c < 4; ++apex_c)
            {
                bool any_incident_cut = false;
                for (const auto& ec : tet_cuts)
                {
                    if (detail::kTetEdges[ec.first][0] == apex_c ||
                        detail::kTetEdges[ec.first][1] == apex_c)
                    {
                        any_incident_cut = true;
                        break;
                    }
                }
                if (!any_incident_cut)
                {
                    single_face_apex = apex_c;
                    break;
                }
            }

            if (single_face_apex >= 0)
            {
                int i0 = -1, i1 = -1, i2 = -1;
                for (int c = 0; c < 4; ++c)
                {
                    if (c != single_face_apex)
                    {
                        if (i0 < 0)
                        {
                            i0 = c;
                        }
                        else if (i1 < 0)
                        {
                            i1 = c;
                        }
                        else
                        {
                            i2 = c;
                        }
                    }
                }
                std::vector<std::array<Index, 3>> base_sub_tris;
                subdivide_face_canonically({ conn[i0], conn[i1], conn[i2] },
                                           edge_cut_map,
                                           face_steiner_map,
                                           st.tetrahedralization.vertices,
                                           base_sub_tris);
                const Index v_apex = conn[single_face_apex];
                for (const auto& st_tri : base_sub_tris)
                {
                    new_tetrahedra.push_back(make_positive_tet(st.tetrahedralization.vertices,
                                                               st_tri[0],
                                                               st_tri[1],
                                                               st_tri[2],
                                                               v_apex));
                }
            }
            else
            {
                bool handled_without_tet_steiner = false;

                // Case 1-vs-3: 3 cut edges sharing a single corner `cut_apex_c` ->
                // exact 4-tet Steiner-free Dompelaar prism + tip decomposition (1 flat medial tri).
                if (tet_cuts.size() == 3)
                {
                    int cut_apex_c = -1;
                    for (int cand_c = 0; cand_c < 4; ++cand_c)
                    {
                        bool all_inc = true;
                        for (const auto& ec : tet_cuts)
                        {
                            if (detail::kTetEdges[ec.first][0] != cand_c &&
                                detail::kTetEdges[ec.first][1] != cand_c)
                            {
                                all_inc = false;
                                break;
                            }
                        }
                        if (all_inc)
                        {
                            cut_apex_c = cand_c;
                            break;
                        }
                    }

                    if (cut_apex_c >= 0)
                    {
                        const Index v_apex = conn[cut_apex_c];
                        std::array<Index, 3> b_verts;
                        int bi = 0;
                        for (int c = 0; c < 4; ++c)
                        {
                            if (c != cut_apex_c)
                            {
                                b_verts[bi++] = conn[c];
                            }
                        }
                        std::sort(b_verts.begin(), b_verts.end());

                        const Index q0 = edge_cut_map.at(hash_linear_edge(v_apex, b_verts[0]));
                        const Index q1 = edge_cut_map.at(hash_linear_edge(v_apex, b_verts[1]));
                        const Index q2 = edge_cut_map.at(hash_linear_edge(v_apex, b_verts[2]));

                        auto& V = st.tetrahedralization.vertices;
                        new_tetrahedra.push_back(make_positive_tet(V, v_apex, q0, q1, q2));
                        new_tetrahedra.push_back(
                          make_positive_tet(V, b_verts[0], b_verts[1], b_verts[2], q2));
                        new_tetrahedra.push_back(
                          make_positive_tet(V, b_verts[0], b_verts[1], q1, q2));
                        new_tetrahedra.push_back(make_positive_tet(V, b_verts[0], q0, q1, q2));
                        handled_without_tet_steiner = true;
                    }
                }
                // Case 2-vs-2: 4 cut edges leaving 2 opposite uncut edges ->
                // exact 6-tet Steiner-free Dompelaar bi-prism decomposition (2 flat medial tris).
                else if (tet_cuts.size() == 4)
                {
                    std::array<bool, 6> is_cut_edge = { false, false, false, false, false, false };
                    for (const auto& ec : tet_cuts)
                    {
                        is_cut_edge[ec.first] = true;
                    }
                    int ue0 = -1, ue1 = -1;
                    for (int e = 0; e < 6; ++e)
                    {
                        if (!is_cut_edge[e])
                        {
                            if (ue0 < 0)
                            {
                                ue0 = e;
                            }
                            else
                            {
                                ue1 = e;
                            }
                        }
                    }

                    if (ue0 >= 0 && ue1 >= 0)
                    {
                        const int u0a = detail::kTetEdges[ue0][0];
                        const int u0b = detail::kTetEdges[ue0][1];
                        const int u1a = detail::kTetEdges[ue1][0];
                        const int u1b = detail::kTetEdges[ue1][1];
                        if (u0a != u1a && u0a != u1b && u0b != u1a && u0b != u1b)
                        {
                            Index a0 = std::min(conn[u0a], conn[u0b]);
                            Index a1 = std::max(conn[u0a], conn[u0b]);
                            Index b0 = std::min(conn[u1a], conn[u1b]);
                            Index b1 = std::max(conn[u1a], conn[u1b]);
                            if (a0 > b0)
                            {
                                std::swap(a0, b0);
                                std::swap(a1, b1);
                            }

                            const Index q00 = edge_cut_map.at(hash_linear_edge(a0, b0));
                            const Index q01 = edge_cut_map.at(hash_linear_edge(a0, b1));
                            const Index q10 = edge_cut_map.at(hash_linear_edge(a1, b0));
                            const Index q11 = edge_cut_map.at(hash_linear_edge(a1, b1));

                            auto& V = st.tetrahedralization.vertices;
                            new_tetrahedra.push_back(make_positive_tet(V, a0, a1, q10, q11));
                            new_tetrahedra.push_back(make_positive_tet(V, a0, q00, q10, q11));
                            new_tetrahedra.push_back(make_positive_tet(V, a0, q00, q01, q11));
                            new_tetrahedra.push_back(make_positive_tet(V, b0, b1, q01, q11));
                            new_tetrahedra.push_back(make_positive_tet(V, b0, q00, q01, q11));
                            new_tetrahedra.push_back(make_positive_tet(V, b0, q00, q10, q11));
                            handled_without_tet_steiner = true;
                        }
                    }
                }

                if (!handled_without_tet_steiner)
                {
                    Point3 p_center = { 0.0, 0.0, 0.0 };
                    for (const auto& ec : tet_cuts)
                    {
                        const auto& pt = st.tetrahedralization.vertices[ec.second];
                        p_center[0] += pt[0];
                        p_center[1] += pt[1];
                        p_center[2] += pt[2];
                    }
                    const double inv_n = 1.0 / static_cast<double>(tet_cuts.size());
                    p_center[0] *= inv_n;
                    p_center[1] *= inv_n;
                    p_center[2] *= inv_n;

                    tet_steiner_vid = static_cast<Index>(st.tetrahedralization.vertices.size());
                    st.tetrahedralization.vertices.push_back(p_center);

                    std::vector<std::array<Index, 3>> face_sub_tris;
                    for (int f = 0; f < 4; ++f)
                    {
                        const auto& fc = detail::kTetFaces[f];
                        face_sub_tris.clear();
                        subdivide_face_canonically({ conn[fc[0]], conn[fc[1]], conn[fc[2]] },
                                                   edge_cut_map,
                                                   face_steiner_map,
                                                   st.tetrahedralization.vertices,
                                                   face_sub_tris);
                        for (const auto& st_tri : face_sub_tris)
                        {
                            new_tetrahedra.push_back(make_positive_tet(
                              st.tetrahedralization.vertices,
                              st_tri[0],
                              st_tri[1],
                              st_tri[2],
                              tet_steiner_vid));
                        }
                    }
                }
            }
        }

        const std::size_t sub_end = new_tetrahedra.size();
        if (!is_ov)
        {
            const RegionID owner = orig_region_owner[orig_r];
            for (std::size_t s = sub_begin; s < sub_end; ++s)
            {
                new_tet_meta.push_back(SubTetMeta{ false, 0, owner, orig_r });
            }
        }
        else
        {
            const auto& gf             = group_fields[g_idx];
            const auto& allowed_global = st.object_region_info.region_to_objects_map[orig_r];

            const Eigen::Vector3d tet_p0 =
              detail::point_as_vector(st.tetrahedralization.vertices[conn[0]]);
            Eigen::Matrix3d edge_mat;
            edge_mat.col(0) =
              detail::point_as_vector(st.tetrahedralization.vertices[conn[1]]) - tet_p0;
            edge_mat.col(1) =
              detail::point_as_vector(st.tetrahedralization.vertices[conn[2]]) - tet_p0;
            edge_mat.col(2) =
              detail::point_as_vector(st.tetrahedralization.vertices[conn[3]]) - tet_p0;
            Eigen::Matrix3d inv_edge_mat = Eigen::Matrix3d::Zero();
            bool has_inv_edge_mat        = false;
            const double edge_det        = edge_mat.determinant();
            if (std::abs(edge_det) > 1e-28)
            {
                inv_edge_mat     = edge_mat.inverse();
                has_inv_edge_mat = inv_edge_mat.allFinite();
            }

            for (std::size_t s = sub_begin; s < sub_end; ++s)
            {
                const auto& sc = new_tetrahedra[s];

                std::vector<RegionID> cand = allowed_global;
                std::sort(cand.begin(), cand.end());
                bool has_orig_corner = false;
                for (int c = 0; c < 4; ++c)
                {
                    if (sc[c] < n_orig_verts)
                    {
                        has_orig_corner  = true;
                        const int lv     = static_cast<int>(gf.global_to_overlap_vertex_map.at(sc[c]));
                        const auto& lset = gf.vertex_label_sets[lv];
                        std::vector<RegionID> next_cand;
                        std::set_intersection(cand.begin(),
                                              cand.end(),
                                              lset.begin(),
                                              lset.end(),
                                              std::back_inserter(next_cand));
                        cand = std::move(next_cand);
                    }
                }

                RegionID chosen_owner = NOT_SET_RegionID;
                if (has_orig_corner && cand.size() == 1)
                {
                    chosen_owner = cand[0];
                }
                else
                {
                    const auto& search_objs = (!cand.empty()) ? cand : allowed_global;
                    double best_score       = -1.0;
                    for (const RegionID obj_id : search_objs)
                    {
                        const int k = static_cast<int>(gf.object_id_to_local_id_map.at(obj_id));
                        std::array<double, 4> corner_w;
                        for (int c = 0; c < 4; ++c)
                        {
                            const int lv =
                              static_cast<int>(gf.global_to_overlap_vertex_map.at(conn[c]));
                            corner_w[c] = gf.weights[k](lv);
                        }
                        double score = 0.0;
                        for (int c = 0; c < 4; ++c)
                        {
                            score += eval_weight_at_vert(gf,
                                                         k,
                                                         sc[c],
                                                         corner_w,
                                                         tet_p0,
                                                         inv_edge_mat,
                                                         has_inv_edge_mat);
                        }
                        if (chosen_owner == NOT_SET_RegionID || score > best_score + 1e-14 ||
                            (std::abs(score - best_score) <= 1e-14 && obj_id < chosen_owner))
                        {
                            best_score   = score;
                            chosen_owner = obj_id;
                        }
                    }
                }

                new_tet_meta.push_back(SubTetMeta{ true, gf.group_id, chosen_owner, orig_r });
            }
        }
    }

    // 2-manifold edge-star and vertex-star repair
    {
        const std::size_t n_sub_tets = new_tetrahedra.size();
        std::vector<std::array<Index, 4>> sub_tet_adj(
          n_sub_tets,
          { NOT_SET_Index, NOT_SET_Index, NOT_SET_Index, NOT_SET_Index });
        std::vector<double> sub_tet_vol(n_sub_tets, 0.0);

        {
            std::unordered_map<std::array<Index, 3>, std::pair<Index, ConnIndex>, SortedTriKeyHash>
              first_on_face;
            first_on_face.reserve(n_sub_tets * 3);
            for (std::size_t tid = 0; tid < n_sub_tets; ++tid)
            {
                const auto& tc   = new_tetrahedra[tid];
                sub_tet_vol[tid] = std::max(
                  0.0,
                  detail::tetrahedron_signed_volume(st.tetrahedralization.vertices[tc[0]],
                                                    st.tetrahedralization.vertices[tc[1]],
                                                    st.tetrahedralization.vertices[tc[2]],
                                                    st.tetrahedralization.vertices[tc[3]]));

                for (ConnIndex f = 0; f < 4; ++f)
                {
                    const auto& fc           = detail::kTetFaces[f];
                    std::array<Index, 3> key = { tc[fc[0]], tc[fc[1]], tc[fc[2]] };
                    std::sort(key.begin(), key.end());
                    auto [it, ins] =
                      first_on_face.emplace(key, std::make_pair(static_cast<Index>(tid), f));
                    if (!ins)
                    {
                        sub_tet_adj[tid][f] = it->second.first;
                        sub_tet_adj[it->second.first][it->second.second] =
                          static_cast<Index>(tid);
                    }
                }
            }
        }

        std::unordered_map<std::size_t, std::vector<Index>> ov_edge_stars;
        std::unordered_map<Index, std::vector<Index>> ov_vert_stars;
        for (std::size_t tid = 0; tid < n_sub_tets; ++tid)
        {
            if (!new_tet_meta[tid].is_overlap)
            {
                continue;
            }
            const auto& tc = new_tetrahedra[tid];
            for (int e = 0; e < 6; ++e)
            {
                ov_edge_stars[hash_linear_edge(tc[detail::kTetEdges[e][0]],
                                               tc[detail::kTetEdges[e][1]])];
            }
            for (int c = 0; c < 4; ++c)
            {
                ov_vert_stars[tc[c]];
            }
        }
        for (std::size_t tid = 0; tid < n_sub_tets; ++tid)
        {
            const auto& tc = new_tetrahedra[tid];
            for (int e = 0; e < 6; ++e)
            {
                auto it = ov_edge_stars.find(
                  hash_linear_edge(tc[detail::kTetEdges[e][0]], tc[detail::kTetEdges[e][1]]));
                if (it != ov_edge_stars.end())
                {
                    it->second.push_back(static_cast<Index>(tid));
                }
            }
            for (int c = 0; c < 4; ++c)
            {
                auto it = ov_vert_stars.find(tc[c]);
                if (it != ov_vert_stars.end())
                {
                    it->second.push_back(static_cast<Index>(tid));
                }
            }
        }

        auto repair_star = [&](const std::vector<Index>& star_tets) -> bool
        {
            if (star_tets.size() < 2)
            {
                return false;
            }
            std::unordered_set<Index> star_set(star_tets.begin(), star_tets.end());
            std::set<RegionID> owners_present;
            for (const Index t : star_tets)
            {
                if (new_tet_meta[t].owner != NOT_SET_RegionID)
                {
                    owners_present.insert(new_tet_meta[t].owner);
                }
            }
            if (owners_present.size() <= 1)
            {
                return false;
            }

            bool modified = false;
            for (const RegionID owner_o : owners_present)
            {
                std::vector<std::vector<Index>> comps;
                std::unordered_set<Index> visited;
                for (const Index seed : star_tets)
                {
                    if (new_tet_meta[seed].owner != owner_o || visited.count(seed) > 0)
                    {
                        continue;
                    }
                    std::vector<Index> comp;
                    std::deque<Index> q = { seed };
                    visited.insert(seed);
                    while (!q.empty())
                    {
                        const Index cur = q.front();
                        q.pop_front();
                        comp.push_back(cur);
                        for (int f = 0; f < 4; ++f)
                        {
                            const Index nb = sub_tet_adj[cur][f];
                            if (nb != NOT_SET_Index && star_set.count(nb) > 0 &&
                                new_tet_meta[nb].owner == owner_o && visited.count(nb) == 0)
                            {
                                visited.insert(nb);
                                q.push_back(nb);
                            }
                        }
                    }
                    comps.push_back(std::move(comp));
                }

                if (comps.size() <= 1)
                {
                    continue;
                }

                int keep_idx    = -1;
                bool keep_fixed = false;
                double keep_vol = -1.0;
                for (std::size_t ci = 0; ci < comps.size(); ++ci)
                {
                    bool has_fixed = false;
                    double c_vol   = 0.0;
                    for (const Index t : comps[ci])
                    {
                        if (!new_tet_meta[t].is_overlap)
                        {
                            has_fixed = true;
                        }
                        c_vol += sub_tet_vol[t];
                    }
                    if (keep_idx < 0 || (has_fixed && !keep_fixed) ||
                        (has_fixed == keep_fixed && c_vol > keep_vol))
                    {
                        keep_idx   = static_cast<int>(ci);
                        keep_fixed = has_fixed;
                        keep_vol   = c_vol;
                    }
                }

                for (std::size_t ci = 0; ci < comps.size(); ++ci)
                {
                    if (static_cast<int>(ci) == keep_idx)
                    {
                        continue;
                    }
                    bool all_ov = true;
                    std::vector<RegionID> common_allowed;
                    for (std::size_t idx = 0; idx < comps[ci].size(); ++idx)
                    {
                        const Index t = comps[ci][idx];
                        if (!new_tet_meta[t].is_overlap)
                        {
                            all_ov = false;
                            break;
                        }
                        std::vector<RegionID> r_objs =
                          st.object_region_info.region_to_objects_map[new_tet_meta[t].orig_region];
                        std::sort(r_objs.begin(), r_objs.end());
                        if (idx == 0)
                        {
                            common_allowed = std::move(r_objs);
                        }
                        else
                        {
                            std::vector<RegionID> inter;
                            std::set_intersection(common_allowed.begin(),
                                                  common_allowed.end(),
                                                  r_objs.begin(),
                                                  r_objs.end(),
                                                  std::back_inserter(inter));
                            common_allowed = std::move(inter);
                        }
                    }
                    if (!all_ov || common_allowed.empty())
                    {
                        continue;
                    }

                    std::map<RegionID, int> adj_counts;
                    for (const Index t : comps[ci])
                    {
                        for (int f = 0; f < 4; ++f)
                        {
                            const Index nb = sub_tet_adj[t][f];
                            if (nb != NOT_SET_Index && star_set.count(nb) > 0)
                            {
                                const RegionID nb_o = new_tet_meta[nb].owner;
                                if (nb_o != owner_o && nb_o != NOT_SET_RegionID &&
                                    std::binary_search(common_allowed.begin(),
                                                       common_allowed.end(),
                                                       nb_o))
                                {
                                    adj_counts[nb_o]++;
                                }
                            }
                        }
                    }

                    RegionID best_alt = NOT_SET_RegionID;
                    int best_cnt      = -1;
                    for (const auto& kv : adj_counts)
                    {
                        if (kv.second > best_cnt)
                        {
                            best_cnt = kv.second;
                            best_alt = kv.first;
                        }
                    }
                    if (best_alt != NOT_SET_RegionID)
                    {
                        for (const Index t : comps[ci])
                        {
                            new_tet_meta[t].owner = best_alt;
                        }
                        modified = true;
                    }
                }
            }
            return modified;
        };

        for (int pass = 0; pass < 24; ++pass)
        {
            bool any_change = false;
            for (const auto& kv : ov_edge_stars)
            {
                if (repair_star(kv.second))
                {
                    any_change = true;
                }
            }
            for (const auto& kv : ov_vert_stars)
            {
                if (repair_star(kv.second))
                {
                    any_change = true;
                }
            }
            if (!any_change)
            {
                break;
            }
        }
    }

    struct OrigSubTriInfo
    {
        std::array<Index, 3> oriented_conn;
        FacetID facet_id             = 0;
        RegionBoundaryID boundary_id = 0;
    };
    std::unordered_map<std::array<Index, 3>, OrigSubTriInfo, SortedTriKeyHash> orig_sub_tri_lookup;
    orig_sub_tri_lookup.reserve(n_orig_tris * 2);

    std::vector<std::array<Index, 3>> sub_tris_buf;
    for (std::size_t tr = 0; tr < n_orig_tris; ++tr)
    {
        const auto& tri_conn       = st.tetrahedralization.triangles[tr];
        const FacetID fid          = st.region_info.tri_to_facet_map[tr];
        const RegionBoundaryID bid = st.region_boundaries_info.triangle_to_boundary_id_map[tr];

        sub_tris_buf.clear();
        subdivide_face_canonically(tri_conn,
                                   edge_cut_map,
                                   face_steiner_map,
                                   st.tetrahedralization.vertices,
                                   sub_tris_buf);
        for (const auto& st_conn : sub_tris_buf)
        {
            auto key = st_conn;
            std::sort(key.begin(), key.end());
            orig_sub_tri_lookup.emplace(key, OrigSubTriInfo{ st_conn, fid, bid });
        }
    }

    struct TetFaceIncidence
    {
        Index tet0                = NOT_SET_Index;
        Index tet1                = NOT_SET_Index;
        std::array<Index, 3> conn = { 0, 0, 0 };
    };
    std::unordered_map<std::array<Index, 3>, TetFaceIncidence, SortedTriKeyHash> all_tet_faces;
    all_tet_faces.reserve(new_tetrahedra.size() * 3);
    std::vector<std::array<Index, 3>> face_keys_in_order;
    face_keys_in_order.reserve(new_tetrahedra.size() * 3);

    for (std::size_t tid = 0; tid < new_tetrahedra.size(); ++tid)
    {
        const auto& tc = new_tetrahedra[tid];
        for (int f = 0; f < 4; ++f)
        {
            const auto& fc             = detail::kTetFaces[f];
            std::array<Index, 3> f_conn = { tc[fc[0]], tc[fc[1]], tc[fc[2]] };
            auto key                   = f_conn;
            std::sort(key.begin(), key.end());
            auto [it, inserted] = all_tet_faces.emplace(
              key,
              TetFaceIncidence{ static_cast<Index>(tid), NOT_SET_Index, f_conn });
            if (inserted)
            {
                face_keys_in_order.push_back(key);
            }
            else
            {
                it->second.tet1 = static_cast<Index>(tid);
            }
        }
    }

    const FacetID n_consistent_facets = static_cast<FacetID>(
      consistent_surfaces.consistent_facet_to_inconsistent_facets_map.size());
    FacetID max_existing_facet = 0;
    if (!st.region_info.tri_to_facet_map.empty())
    {
        max_existing_facet =
          *std::max_element(st.region_info.tri_to_facet_map.begin(),
                            st.region_info.tri_to_facet_map.end()) +
          1;
    }
    const FacetID medial_facet_base = std::max(n_consistent_facets, max_existing_facet);
    FacetID next_medial_facet       = medial_facet_base;
    std::map<std::tuple<RegionID, RegionID, RegionID>, FacetID> medial_pair_to_facet;

    std::vector<std::array<Index, 3>> final_triangles;
    std::vector<FacetID> final_tri_to_facet_map;
    final_triangles.reserve(orig_sub_tri_lookup.size() + 256);
    final_tri_to_facet_map.reserve(orig_sub_tri_lookup.size() + 256);

    for (const auto& key : face_keys_in_order)
    {
        const auto& inc     = all_tet_faces.at(key);
        const auto it_orig  = orig_sub_tri_lookup.find(key);
        const bool has_orig = (it_orig != orig_sub_tri_lookup.end());

        if (inc.tet1 == NOT_SET_Index)
        {
            if (has_orig)
            {
                final_triangles.push_back(it_orig->second.oriented_conn);
                final_tri_to_facet_map.push_back(it_orig->second.facet_id);
            }
            else
            {
                final_triangles.push_back(inc.conn);
                final_tri_to_facet_map.push_back(0);
            }
            continue;
        }

        const auto& m0 = new_tet_meta[inc.tet0];
        const auto& m1 = new_tet_meta[inc.tet1];

        if (!m0.is_overlap && !m1.is_overlap)
        {
            if (has_orig)
            {
                final_triangles.push_back(it_orig->second.oriented_conn);
                final_tri_to_facet_map.push_back(it_orig->second.facet_id);
            }
            continue;
        }

        if (m0.owner == m1.owner)
        {
            continue;
        }

        if (has_orig && all_removable_boundaries.count(it_orig->second.boundary_id) == 0 &&
            (!m0.is_overlap || !m1.is_overlap))
        {
            final_triangles.push_back(it_orig->second.oriented_conn);
            final_tri_to_facet_map.push_back(it_orig->second.facet_id);
        }
        else
        {
            const RegionID grp   = m0.is_overlap ? m0.group_id : m1.group_id;
            const RegionID obj_a = std::min(m0.owner, m1.owner);
            const RegionID obj_b = std::max(m0.owner, m1.owner);
            const auto mkey      = std::make_tuple(grp, obj_a, obj_b);
            auto [mit, mins]     = medial_pair_to_facet.emplace(mkey, next_medial_facet);
            if (mins)
            {
                ++next_medial_facet;
            }
            final_triangles.push_back(has_orig ? it_orig->second.oriented_conn : inc.conn);
            final_tri_to_facet_map.push_back(mit->second);
        }
    }

    st.tetrahedralization.tetrahedra = std::move(new_tetrahedra);
    st.tetrahedralization.triangles  = std::move(final_triangles);
    update_edges_and_maps(st.tetrahedralization);

    const std::size_t n_final_tets = st.tetrahedralization.tetrahedra.size();
    std::vector<std::array<Index, 4>> tet_adj(
      n_final_tets,
      { NOT_SET_Index, NOT_SET_Index, NOT_SET_Index, NOT_SET_Index });

    {
        std::unordered_map<std::array<Index, 3>, std::pair<Index, ConnIndex>, SortedTriKeyHash>
          first_tet_on_face;
        first_tet_on_face.reserve(n_final_tets * 3);
        for (std::size_t tid = 0; tid < n_final_tets; ++tid)
        {
            const auto& tc = st.tetrahedralization.tetrahedra[tid];
            for (ConnIndex f = 0; f < 4; ++f)
            {
                const auto& fc           = detail::kTetFaces[f];
                std::array<Index, 3> key = { tc[fc[0]], tc[fc[1]], tc[fc[2]] };
                std::sort(key.begin(), key.end());
                auto [it, inserted] =
                  first_tet_on_face.emplace(key, std::make_pair(static_cast<Index>(tid), f));
                if (!inserted)
                {
                    const Index other_t       = it->second.first;
                    const ConnIndex other_f   = it->second.second;
                    tet_adj[tid][f]           = other_t;
                    tet_adj[other_t][other_f] = static_cast<Index>(tid);
                }
            }
        }
    }

    std::vector<RegionID> new_tet_to_region(n_final_tets, NOT_SET_RegionID);
    std::vector<std::vector<Index>> new_region_to_tets;
    std::vector<RegionID> new_region_owner;

    for (std::size_t seed = 0; seed < n_final_tets; ++seed)
    {
        if (new_tet_to_region[seed] != NOT_SET_RegionID)
        {
            continue;
        }
        const RegionID r_id       = static_cast<RegionID>(new_region_to_tets.size());
        const RegionID seed_owner = new_tet_meta[seed].owner;
        new_region_to_tets.emplace_back();
        new_region_owner.push_back(seed_owner);

        std::deque<Index> bfs   = { static_cast<Index>(seed) };
        new_tet_to_region[seed] = r_id;
        while (!bfs.empty())
        {
            const Index cur = bfs.front();
            bfs.pop_front();
            new_region_to_tets[r_id].push_back(cur);

            for (ConnIndex f = 0; f < 4; ++f)
            {
                if (st.tetrahedralization.tet_to_tri_map[cur][f] == NOT_SET_Index)
                {
                    const Index nb = tet_adj[cur][f];
                    if (nb != NOT_SET_Index && new_tet_to_region[nb] == NOT_SET_RegionID)
                    {
                        assert(new_tet_meta[nb].owner == seed_owner);
                        new_tet_to_region[nb] = r_id;
                        bfs.push_back(nb);
                    }
                }
            }
        }
        std::sort(new_region_to_tets[r_id].begin(), new_region_to_tets[r_id].end());
    }

    const RegionID n_new_regions = static_cast<RegionID>(new_region_to_tets.size());

    st.region_info.n_regions          = n_new_regions;
    st.region_info.tet_to_region_map  = std::move(new_tet_to_region);
    st.region_info.region_to_tets_map = std::move(new_region_to_tets);
    st.region_info.tri_to_facet_map   = std::move(final_tri_to_facet_map);

    const std::size_t n_final_tris = st.tetrahedralization.triangles.size();
    st.region_info.tri_to_bounded_regions_map.assign(
      n_final_tris,
      { NOT_SET_RegionID, NOT_SET_RegionID });
    st.region_info.region_to_bounding_triangles_map.assign(n_new_regions, {});

    for (std::size_t tr = 0; tr < n_final_tris; ++tr)
    {
        const auto& tt = st.tetrahedralization.tri_to_tet_map[tr];
        auto& br       = st.region_info.tri_to_bounded_regions_map[tr];
        if (tt[0] != NOT_SET_Index)
        {
            br.first = st.region_info.tet_to_region_map[tt[0]];
            st.region_info.region_to_bounding_triangles_map[br.first].push_back(
              static_cast<Index>(tr));
        }
        if (tt[1] != NOT_SET_Index)
        {
            br.second = st.region_info.tet_to_region_map[tt[1]];
            if (br.second != br.first)
            {
                st.region_info.region_to_bounding_triangles_map[br.second].push_back(
                  static_cast<Index>(tr));
            }
        }
        if (br.first > br.second)
        {
            std::swap(br.first, br.second);
        }
    }

    const std::size_t n_objects = st.object_region_info.object_to_regions_map.size();
    st.object_region_info.regions_of_overlap.clear();
    st.object_region_info.region_types.assign(n_new_regions, RegionType::KEPT);
    st.object_region_info.region_to_objects_map.assign(n_new_regions, {});
    st.object_region_info.object_to_regions_map.assign(n_objects, {});

    for (RegionID r = 0; r < n_new_regions; ++r)
    {
        const RegionID owner                           = new_region_owner[r];
        st.object_region_info.region_to_objects_map[r] = { owner };
        st.object_region_info.object_to_regions_map[owner].push_back(r);
    }

    std::vector<std::set<FacetID>> object_pred_facet_sets(n_objects);
    for (std::size_t obj = 0; obj < n_objects; ++obj)
    {
        object_pred_facet_sets[obj].insert(
          inconsistent_combined_surfaces.object_to_facets_map[obj].begin(),
          inconsistent_combined_surfaces.object_to_facets_map[obj].end());
    }

    for (std::size_t tr = 0; tr < n_final_tris; ++tr)
    {
        FacetID& fid = st.region_info.tri_to_facet_map[tr];
        if (fid >= n_consistent_facets)
        {
            continue;
        }
        const auto& pred_facets =
          consistent_surfaces.consistent_facet_to_inconsistent_facets_map[fid];
        if (pred_facets.size() <= 1)
        {
            continue;
        }

        const auto& br = st.region_info.tri_to_bounded_regions_map[tr];
        std::set<FacetID> matched_preds;
        auto collect_for_region = [&](RegionID r)
        {
            if (r == NOT_SET_RegionID)
            {
                return;
            }
            const RegionID owner = new_region_owner[r];
            for (const FacetID pf : pred_facets)
            {
                if (object_pred_facet_sets[owner].count(pf) > 0)
                {
                    matched_preds.insert(pf);
                }
            }
        };
        collect_for_region(br.first);
        collect_for_region(br.second);

        if (matched_preds.size() == 1)
        {
            fid = *matched_preds.begin();
        }
    }

    if (!check_tri_edge_relations(st.tetrahedralization) ||
        !check_tet_edge_relations(st.tetrahedralization) ||
        !check_tet_tri_relations(st.tetrahedralization))
    {
        throw std::runtime_error(
          "Topological relation check failed after multi-material interpenetration removal.");
    }
}

double compute_tri_quality_2d(const Point3& p0, const Point3& p1, const Point3& p2) noexcept
{
    const Eigen::Vector3d a = detail::point_as_vector(p0);
    const Eigen::Vector3d b = detail::point_as_vector(p1);
    const Eigen::Vector3d c = detail::point_as_vector(p2);
    const double l2_sum =
      (b - a).squaredNorm() + (c - b).squaredNorm() + (a - c).squaredNorm();
    if (l2_sum <= 1e-30)
    {
        return 0.0;
    }
    const double dbl_area = (b - a).cross(c - a).norm();
    return (2.0 * std::sqrt(3.0) * dbl_area) / l2_sum;
}

Eigen::MatrixXi detect_surface_self_intersections(const Triangulation& t)
{
    Eigen::MatrixXd V(static_cast<Eigen::Index>(t.vertices.size()), 3);
    for (std::size_t i = 0; i < t.vertices.size(); ++i)
    {
        V(static_cast<Eigen::Index>(i), 0) = t.vertices[i][0];
        V(static_cast<Eigen::Index>(i), 1) = t.vertices[i][1];
        V(static_cast<Eigen::Index>(i), 2) = t.vertices[i][2];
    }
    Eigen::MatrixXi F(static_cast<Eigen::Index>(t.triangles.size()), 3);
    for (std::size_t i = 0; i < t.triangles.size(); ++i)
    {
        F(static_cast<Eigen::Index>(i), 0) = static_cast<int>(t.triangles[i][0]);
        F(static_cast<Eigen::Index>(i), 1) = static_cast<int>(t.triangles[i][1]);
        F(static_cast<Eigen::Index>(i), 2) = static_cast<int>(t.triangles[i][2]);
    }

    igl::copyleft::cgal::RemeshSelfIntersectionsParam params;
    params.detect_only = true;
    params.first_only  = false;
    params.stitch_all  = false;

    Eigen::MatrixXd VV;
    Eigen::MatrixXi FF;
    Eigen::MatrixXi IF;
    Eigen::VectorXi J;
    Eigen::VectorXi IM;
    igl::copyleft::cgal::remesh_self_intersections(V, F, params, VV, FF, IF, J, IM);
    return IF;
}

} // namespace

void smooth_medial_surfaces(const ConsistentSurfaces& consistent_surfaces,
                            SegmentedTetrahedralization& st)
{
    auto& tet_mesh             = st.tetrahedralization;
    const std::size_t n_verts  = tet_mesh.vertices.size();
    const std::size_t n_tris   = tet_mesh.triangles.size();
    const std::size_t n_tets   = tet_mesh.tetrahedra.size();
    const FacetID medial_start = static_cast<FacetID>(
      consistent_surfaces.consistent_facet_to_inconsistent_facets_map.size());

    if (n_verts == 0 || n_tris == 0 || n_tets == 0)
    {
        return;
    }

    std::vector<bool> is_medial_tri(n_tris, false);
    bool any_medial = false;
    for (std::size_t tr = 0; tr < n_tris; ++tr)
    {
        if (tr < st.region_info.tri_to_facet_map.size() &&
            st.region_info.tri_to_facet_map[tr] >= medial_start)
        {
            is_medial_tri[tr] = true;
            any_medial        = true;
        }
    }
    if (!any_medial)
    {
        return;
    }

    std::vector<bool> touches_non_medial(n_verts, false);
    std::vector<std::set<RegionBoundaryID>> vert_medial_bids(n_verts);
    std::vector<std::vector<Index>> medial_surf_nbrs(n_verts);

    for (std::size_t tr = 0; tr < n_tris; ++tr)
    {
        const auto& tc = tet_mesh.triangles[tr];
        if (!is_medial_tri[tr])
        {
            touches_non_medial[tc[0]] = true;
            touches_non_medial[tc[1]] = true;
            touches_non_medial[tc[2]] = true;
        }
        else
        {
            const RegionBoundaryID bid =
              (tr < st.region_boundaries_info.triangle_to_boundary_id_map.size())
                ? st.region_boundaries_info.triangle_to_boundary_id_map[tr]
                : 0;
            for (int c = 0; c < 3; ++c)
            {
                vert_medial_bids[tc[c]].insert(bid);
                medial_surf_nbrs[tc[c]].push_back(tc[(c + 1) % 3]);
                medial_surf_nbrs[tc[c]].push_back(tc[(c + 2) % 3]);
            }
        }
    }

    std::vector<std::vector<Index>> vert_seam_edges(n_verts);
    for (const Index eid : st.seams_info.edges_on_seams)
    {
        if (eid < tet_mesh.edges.size())
        {
            const auto& ed = tet_mesh.edges[eid];
            vert_seam_edges[ed[0]].push_back(eid);
            vert_seam_edges[ed[1]].push_back(eid);
        }
    }

    // Classify vertices:
    // 0 = locked, 1 = interior medial sheet (2D), 2 = interior multi-medial junction curve (1D),
    // 3 = strictly interior 3D volume vertex near medial sheet.
    std::vector<uint8_t> vert_mode(n_verts, 0);
    std::vector<std::array<Index, 2>> junc_curve_nbrs(n_verts, { NOT_SET_Index, NOT_SET_Index });
    std::vector<Index> medial_and_junc_verts;
    std::vector<bool> on_medial_closure(n_verts, false);

    for (std::size_t v = 0; v < n_verts; ++v)
    {
        if (!vert_medial_bids[v].empty())
        {
            on_medial_closure[v] = true;
        }
        if (touches_non_medial[v])
        {
            continue;
        }
        if (vert_seam_edges[v].empty() && vert_medial_bids[v].size() == 1)
        {
            auto& nb = medial_surf_nbrs[v];
            std::sort(nb.begin(), nb.end());
            nb.erase(std::unique(nb.begin(), nb.end()), nb.end());
            if (nb.size() >= 3)
            {
                vert_mode[v] = 1;
                medial_and_junc_verts.push_back(static_cast<Index>(v));
            }
        }
        else if (vert_seam_edges[v].size() == 2 && !vert_medial_bids[v].empty())
        {
            const Index e0  = vert_seam_edges[v][0];
            const Index e1  = vert_seam_edges[v][1];
            const SeamID s0 = st.seams_info.edge_id_to_seam_id_map[e0];
            const SeamID s1 = st.seams_info.edge_id_to_seam_id_map[e1];
            if (s0 != NOT_SET_SeamID && s0 == s1)
            {
                const Index n0 = (tet_mesh.edges[e0][0] == v) ? tet_mesh.edges[e0][1]
                                                              : tet_mesh.edges[e0][0];
                const Index n1 = (tet_mesh.edges[e1][0] == v) ? tet_mesh.edges[e1][1]
                                                              : tet_mesh.edges[e1][0];
                if (n0 != n1)
                {
                    vert_mode[v]       = 2;
                    junc_curve_nbrs[v] = { n0, n1 };
                    medial_and_junc_verts.push_back(static_cast<Index>(v));
                }
            }
        }
    }

    if (medial_and_junc_verts.empty())
    {
        return;
    }

    // Collect strictly-interior 3D volume vertices within 2 hops of medial sheets so surrounding
    // tetrahedra can relax smoothly out of the way during medial fairing.
    std::vector<bool> near_hop1(n_verts, false);
    for (std::size_t ti = 0; ti < n_tets; ++ti)
    {
        const auto& tc = tet_mesh.tetrahedra[ti];
        bool touches_m = false;
        for (int c = 0; c < 4; ++c)
        {
            if (on_medial_closure[tc[c]])
            {
                touches_m = true;
                break;
            }
        }
        if (touches_m)
        {
            for (int c = 0; c < 4; ++c)
            {
                near_hop1[tc[c]] = true;
            }
        }
    }
    std::vector<bool> near_hop2 = near_hop1;
    for (std::size_t ti = 0; ti < n_tets; ++ti)
    {
        const auto& tc  = tet_mesh.tetrahedra[ti];
        bool touches_h1 = false;
        for (int c = 0; c < 4; ++c)
        {
            if (near_hop1[tc[c]])
            {
                touches_h1 = true;
                break;
            }
        }
        if (touches_h1)
        {
            for (int c = 0; c < 4; ++c)
            {
                near_hop2[tc[c]] = true;
            }
        }
    }

    std::vector<Index> active_vol_verts;
    for (std::size_t v = 0; v < n_verts; ++v)
    {
        if (near_hop2[v] && !touches_non_medial[v] && vert_medial_bids[v].empty() &&
            vert_seam_edges[v].empty())
        {
            vert_mode[v] = 3;
            active_vol_verts.push_back(static_cast<Index>(v));
        }
    }

    std::vector<std::vector<Index>> vert_inc_tets(n_verts);
    std::vector<std::vector<Index>> vol_nbrs(n_verts);
    std::vector<double> orig_tet_vol(n_tets, 0.0);

    for (std::size_t ti = 0; ti < n_tets; ++ti)
    {
        const auto& tc   = tet_mesh.tetrahedra[ti];
        orig_tet_vol[ti] = detail::tetrahedron_signed_volume(tet_mesh.vertices[tc[0]],
                                                             tet_mesh.vertices[tc[1]],
                                                             tet_mesh.vertices[tc[2]],
                                                             tet_mesh.vertices[tc[3]]);
        for (int c = 0; c < 4; ++c)
        {
            const Index v = tc[c];
            if (vert_mode[v] != 0)
            {
                vert_inc_tets[v].push_back(static_cast<Index>(ti));
                if (vert_mode[v] == 3)
                {
                    for (int d = 0; d < 4; ++d)
                    {
                        if (d != c)
                        {
                            vol_nbrs[v].push_back(tc[d]);
                        }
                    }
                }
            }
        }
    }
    for (const Index v : active_vol_verts)
    {
        auto& nb = vol_nbrs[v];
        std::sort(nb.begin(), nb.end());
        nb.erase(std::unique(nb.begin(), nb.end()), nb.end());
    }

    auto try_move_vertex_3d = [&](Index v, const Point3& target, double step_weight)
    {
        const Point3 old_p = tet_mesh.vertices[v];
        const double dx    = target[0] - old_p[0];
        const double dy    = target[1] - old_p[1];
        const double dz    = target[2] - old_p[2];
        if (dx * dx + dy * dy + dz * dz <= 1e-30)
        {
            return;
        }

        double curr_min_vol = std::numeric_limits<double>::max();
        for (const Index ti : vert_inc_tets[v])
        {
            const auto& tc = tet_mesh.tetrahedra[ti];
            curr_min_vol   = std::min(curr_min_vol,
                                    detail::tetrahedron_signed_volume(tet_mesh.vertices[tc[0]],
                                                                      tet_mesh.vertices[tc[1]],
                                                                      tet_mesh.vertices[tc[2]],
                                                                      tet_mesh.vertices[tc[3]]));
        }
        if (curr_min_vol <= 1e-24)
        {
            return;
        }

        for (const double backoff : { 1.0, 0.5, 0.25, 0.125 })
        {
            const double w       = step_weight * backoff;
            tet_mesh.vertices[v] = { old_p[0] + w * dx, old_p[1] + w * dy, old_p[2] + w * dz };
            bool ok              = true;
            for (const Index ti : vert_inc_tets[v])
            {
                const auto& tc  = tet_mesh.tetrahedra[ti];
                const double nv = detail::tetrahedron_signed_volume(tet_mesh.vertices[tc[0]],
                                                                    tet_mesh.vertices[tc[1]],
                                                                    tet_mesh.vertices[tc[2]],
                                                                    tet_mesh.vertices[tc[3]]);
                const double floor_v =
                  std::min(0.20 * orig_tet_vol[ti], 0.80 * curr_min_vol);
                if (nv < floor_v || nv <= 1e-24)
                {
                    ok = false;
                    break;
                }
            }
            if (ok)
            {
                return;
            }
        }
        tet_mesh.vertices[v] = old_p;
    };

    constexpr double kLambda = +0.46;
    constexpr double kMu     = -0.482;
    std::vector<Point3> lap_targets(medial_and_junc_verts.size());

    for (int cycle = 0; cycle < 22; ++cycle)
    {
        if ((cycle % 2) == 0 && !active_vol_verts.empty())
        {
            for (const Index v : active_vol_verts)
            {
                const auto& nb = vol_nbrs[v];
                if (nb.empty())
                {
                    continue;
                }
                Point3 avg = { 0.0, 0.0, 0.0 };
                for (const Index u : nb)
                {
                    avg[0] += tet_mesh.vertices[u][0];
                    avg[1] += tet_mesh.vertices[u][1];
                    avg[2] += tet_mesh.vertices[u][2];
                }
                const double inv = 1.0 / static_cast<double>(nb.size());
                avg[0] *= inv;
                avg[1] *= inv;
                avg[2] *= inv;
                try_move_vertex_3d(v, avg, 0.45);
            }
        }

        for (const double taubin_step : { kLambda, kMu })
        {
            for (std::size_t idx = 0; idx < medial_and_junc_verts.size(); ++idx)
            {
                const Index v = medial_and_junc_verts[idx];
                if (vert_mode[v] == 1)
                {
                    const auto& nb = medial_surf_nbrs[v];
                    Point3 avg     = { 0.0, 0.0, 0.0 };
                    for (const Index u : nb)
                    {
                        avg[0] += tet_mesh.vertices[u][0];
                        avg[1] += tet_mesh.vertices[u][1];
                        avg[2] += tet_mesh.vertices[u][2];
                    }
                    const double inv = 1.0 / static_cast<double>(nb.size());
                    avg[0] *= inv;
                    avg[1] *= inv;
                    avg[2] *= inv;
                    lap_targets[idx] = avg;
                }
                else
                {
                    const auto& jn   = junc_curve_nbrs[v];
                    lap_targets[idx] = detail::midpoint(tet_mesh.vertices[jn[0]],
                                                        tet_mesh.vertices[jn[1]]);
                }
            }

            for (std::size_t idx = 0; idx < medial_and_junc_verts.size(); ++idx)
            {
                try_move_vertex_3d(medial_and_junc_verts[idx], lap_targets[idx], taubin_step);
            }
        }
    }
}

void polish_recovered_medial_surfaces(const ConsistentSurfaces& consistent_surfaces,
                                      const SegmentedTetrahedralization& st,
                                      RecoveredTriangulation& rt)
{
    if (rt.t.triangles.empty() || rt.t.n_markers == 0)
    {
        return;
    }

    const FacetID medial_start = static_cast<FacetID>(
      consistent_surfaces.consistent_facet_to_inconsistent_facets_map.size());
    const std::size_t n_bids = st.region_boundaries_info.n_region_boundaries;

    std::vector<bool> is_medial_marker(rt.t.n_markers, false);
    bool any_medial_marker = false;
    for (Marker m = 0; m < rt.t.n_markers && m < n_bids; ++m)
    {
        const auto& st_tris = st.region_boundaries_info.region_boundary_to_triangles_map[m];
        if (st_tris.empty())
        {
            continue;
        }
        bool all_med = true;
        for (const Index tr : st_tris)
        {
            if (tr >= st.region_info.tri_to_facet_map.size() ||
                st.region_info.tri_to_facet_map[tr] < medial_start)
            {
                all_med = false;
                break;
            }
        }
        if (all_med && m < st.object_region_boundary_info.boundary_id_to_objects_map.size() &&
            st.object_region_boundary_info.boundary_id_to_objects_map[m].empty())
        {
            is_medial_marker[m] = true;
            any_medial_marker   = true;
        }
    }
    if (!any_medial_marker)
    {
        return;
    }

    const Triangulation safe_backup = rt.t;

    auto run_medial_edge_flips = [&](Triangulation& surf, int max_passes)
    {
        for (int pass = 0; pass < max_passes; ++pass)
        {
            std::unordered_set<uint64_t> existing_edges;
            existing_edges.reserve(surf.edges.size() * 2);
            for (const auto& ed : surf.edges)
            {
                existing_edges.insert(hash_linear_edge(ed[0], ed[1]));
            }

            std::vector<bool> tri_busy(surf.triangles.size(), false);
            int flipped_in_pass = 0;

            for (std::size_t eid = 0; eid < surf.edges.size(); ++eid)
            {
                const auto& inc = surf.edge_to_tris_map[eid];
                if (inc.size() != 2)
                {
                    continue;
                }
                const Index t1  = inc[0].first;
                const Index t2  = inc[1].first;
                const Marker m1 = surf.tri_to_marker_map[t1];
                const Marker m2 = surf.tri_to_marker_map[t2];
                if (m1 != m2 || m1 >= is_medial_marker.size() || !is_medial_marker[m1] ||
                    tri_busy[t1] || tri_busy[t2])
                {
                    continue;
                }

                const auto& c1 = surf.triangles[t1];
                const auto& c2 = surf.triangles[t2];
                const Index ea = surf.edges[eid][0];
                const Index eb = surf.edges[eid][1];

                // Orient t1 as (a, b, c) along shared half-edge a -> b
                Index a = NOT_SET_Index, b = NOT_SET_Index, c = NOT_SET_Index, d = NOT_SET_Index;
                for (int i = 0; i < 3; ++i)
                {
                    const Index u = c1[i];
                    const Index v = c1[(i + 1) % 3];
                    if ((u == ea && v == eb) || (u == eb && v == ea))
                    {
                        a = u;
                        b = v;
                        c = c1[(i + 2) % 3];
                        break;
                    }
                }
                if (a == NOT_SET_Index)
                {
                    continue;
                }
                for (int i = 0; i < 3; ++i)
                {
                    const Index u = c2[i];
                    const Index v = c2[(i + 1) % 3];
                    if (u == b && v == a)
                    {
                        d = c2[(i + 2) % 3];
                        break;
                    }
                }
                if (d == NOT_SET_Index || c == d)
                {
                    continue;
                }
                if (existing_edges.count(hash_linear_edge(c, d)) > 0)
                {
                    continue;
                }

                const Eigen::Vector3d pa = detail::point_as_vector(surf.vertices[a]);
                const Eigen::Vector3d pb = detail::point_as_vector(surf.vertices[b]);
                const Eigen::Vector3d pc = detail::point_as_vector(surf.vertices[c]);
                const Eigen::Vector3d pd = detail::point_as_vector(surf.vertices[d]);

                const Eigen::Vector3d n1_raw = (pb - pa).cross(pc - pa);
                const Eigen::Vector3d n2_raw = (pa - pb).cross(pd - pb);
                const double a1              = n1_raw.norm();
                const double a2              = n2_raw.norm();
                if (a1 <= 1e-24 || a2 <= 1e-24)
                {
                    continue;
                }
                const Eigen::Vector3d n1 = n1_raw / a1;
                const Eigen::Vector3d n2 = n2_raw / a2;
                const double old_dot     = n1.dot(n2);
                if (old_dot < 0.50)
                {
                    continue;
                }

                // Candidate flipped triangles: t1' = (d, c, a), t2' = (c, d, b)
                const Eigen::Vector3d np1_raw = (pc - pd).cross(pa - pd);
                const Eigen::Vector3d np2_raw = (pd - pc).cross(pb - pc);
                const double ap1              = np1_raw.norm();
                const double ap2              = np2_raw.norm();
                if (ap1 <= 1e-24 || ap2 <= 1e-24)
                {
                    continue;
                }
                const Eigen::Vector3d np1   = np1_raw / ap1;
                const Eigen::Vector3d np2   = np2_raw / ap2;
                const Eigen::Vector3d n_avg = (n1_raw + n2_raw).normalized();

                if (np1.dot(n_avg) < 0.88 || np2.dot(n_avg) < 0.88)
                {
                    continue;
                }
                const double new_dot = np1.dot(np2);
                if (new_dot < std::min(0.95, old_dot - 0.01))
                {
                    continue;
                }

                const double q_old = std::min(compute_tri_quality_2d(surf.vertices[a],
                                                                     surf.vertices[b],
                                                                     surf.vertices[c]),
                                              compute_tri_quality_2d(surf.vertices[b],
                                                                     surf.vertices[a],
                                                                     surf.vertices[d]));
                const double q_new = std::min(compute_tri_quality_2d(surf.vertices[d],
                                                                     surf.vertices[c],
                                                                     surf.vertices[a]),
                                              compute_tri_quality_2d(surf.vertices[c],
                                                                     surf.vertices[d],
                                                                     surf.vertices[b]));

                if (q_new > q_old + 0.015 || (old_dot < 0.85 && new_dot > 0.96 && q_new > q_old - 0.05))
                {
                    surf.triangles[t1] = { d, c, a };
                    surf.triangles[t2] = { c, d, b };
                    tri_busy[t1]       = true;
                    tri_busy[t2]       = true;
                    existing_edges.erase(hash_linear_edge(a, b));
                    existing_edges.insert(hash_linear_edge(c, d));
                    ++flipped_in_pass;
                }
            }

            if (flipped_in_pass == 0)
            {
                break;
            }
            update_edges_and_maps(surf);
        }
    };

    auto run_medial_short_edge_collapses = [&](Triangulation& surf, int max_passes)
    {
        for (int pass = 0; pass < max_passes; ++pass)
        {
            const std::size_t nv = surf.vertices.size();
            const std::size_t nt = surf.triangles.size();

            std::vector<bool> is_locked_vert(nv, false);
            for (std::size_t eid = 0; eid < surf.edges.size(); ++eid)
            {
                const auto& inc     = surf.edge_to_tris_map[eid];
                bool interior_to_mm = false;
                if (inc.size() == 2)
                {
                    const Marker m0 = surf.tri_to_marker_map[inc[0].first];
                    const Marker m1 = surf.tri_to_marker_map[inc[1].first];
                    if (m0 == m1 && m0 < is_medial_marker.size() && is_medial_marker[m0])
                    {
                        interior_to_mm = true;
                    }
                }
                if (!interior_to_mm)
                {
                    is_locked_vert[surf.edges[eid][0]] = true;
                    is_locked_vert[surf.edges[eid][1]] = true;
                }
            }

            std::vector<std::vector<Index>> v_tris(nv);
            std::vector<std::vector<Index>> v_nbrs(nv);
            for (std::size_t ti = 0; ti < nt; ++ti)
            {
                const Marker m = surf.tri_to_marker_map[ti];
                if (m >= is_medial_marker.size() || !is_medial_marker[m])
                {
                    continue;
                }
                const auto& tc = surf.triangles[ti];
                for (int c = 0; c < 3; ++c)
                {
                    v_tris[tc[c]].push_back(static_cast<Index>(ti));
                    v_nbrs[tc[c]].push_back(tc[(c + 1) % 3]);
                    v_nbrs[tc[c]].push_back(tc[(c + 2) % 3]);
                }
            }
            for (std::size_t v = 0; v < nv; ++v)
            {
                auto& nb = v_nbrs[v];
                std::sort(nb.begin(), nb.end());
                nb.erase(std::unique(nb.begin(), nb.end()), nb.end());
            }

            struct CandEdge
            {
                Index eid  = 0;
                double len = 0.0;
            };
            std::vector<CandEdge> cands;
            for (std::size_t eid = 0; eid < surf.edges.size(); ++eid)
            {
                const auto& inc = surf.edge_to_tris_map[eid];
                if (inc.size() != 2)
                {
                    continue;
                }
                const Marker m0 = surf.tri_to_marker_map[inc[0].first];
                const Marker m1 = surf.tri_to_marker_map[inc[1].first];
                if (m0 != m1 || m0 >= is_medial_marker.size() || !is_medial_marker[m0])
                {
                    continue;
                }
                const Index a = surf.edges[eid][0];
                const Index b = surf.edges[eid][1];
                if (is_locked_vert[a] && is_locked_vert[b])
                {
                    continue;
                }
                cands.push_back(
                  { static_cast<Index>(eid), detail::distance(surf.vertices[a], surf.vertices[b]) });
            }
            std::sort(cands.begin(),
                      cands.end(),
                      [](const CandEdge& x, const CandEdge& y) { return x.len < y.len; });

            std::vector<bool> vert_busy(nv, false);
            std::vector<bool> tri_deleted(nt, false);
            std::vector<Index> collapse_target(nv);
            std::iota(collapse_target.begin(), collapse_target.end(), 0);
            int collapsed = 0;

            for (const auto& ce : cands)
            {
                const Index a = surf.edges[ce.eid][0];
                const Index b = surf.edges[ce.eid][1];
                if (vert_busy[a] || vert_busy[b])
                {
                    continue;
                }

                // Choose `u` (strictly interior medial vertex) -> `v`
                Index u = NOT_SET_Index;
                Index v = NOT_SET_Index;
                if (!is_locked_vert[a] && is_locked_vert[b])
                {
                    u = a;
                    v = b;
                }
                else if (!is_locked_vert[b] && is_locked_vert[a])
                {
                    u = b;
                    v = a;
                }
                else
                {
                    u = (v_nbrs[a].size() <= v_nbrs[b].size()) ? a : b;
                    v = (u == a) ? b : a;
                }

                bool ring_busy = false;
                double sum_ring_len = 0.0;
                for (const Index nb : v_nbrs[u])
                {
                    if (vert_busy[nb])
                    {
                        ring_busy = true;
                        break;
                    }
                    sum_ring_len += detail::distance(surf.vertices[u], surf.vertices[nb]);
                }
                if (ring_busy || v_nbrs[u].size() < 3)
                {
                    continue;
                }

                const double mean_ring_len = sum_ring_len / static_cast<double>(v_nbrs[u].size());
                if (ce.len > 0.28 * mean_ring_len)
                {
                    continue;
                }

                // Exact 2-manifold link condition: |N_m(u) intersect N_m(v)| == 2
                std::vector<Index> common_nb;
                std::set_intersection(v_nbrs[u].begin(),
                                      v_nbrs[u].end(),
                                      v_nbrs[v].begin(),
                                      v_nbrs[v].end(),
                                      std::back_inserter(common_nb));
                if (common_nb.size() != 2)
                {
                    continue;
                }

                // Prevent creating a new interior chord between two seam vertices
                if (is_locked_vert[v])
                {
                    bool creates_seam_chord = false;
                    for (const Index nb : v_nbrs[u])
                    {
                        if (nb != v && nb != common_nb[0] && nb != common_nb[1] &&
                            is_locked_vert[nb])
                        {
                            creates_seam_chord = true;
                            break;
                        }
                    }
                    if (creates_seam_chord)
                    {
                        continue;
                    }
                }

                const Index del_t0 = surf.edge_to_tris_map[ce.eid][0].first;
                const Index del_t1 = surf.edge_to_tris_map[ce.eid][1].first;
                if (tri_deleted[del_t0] || tri_deleted[del_t1])
                {
                    continue;
                }

                Eigen::Vector3d patch_n = Eigen::Vector3d::Zero();
                for (const Index ti : v_tris[u])
                {
                    const auto& tc = surf.triangles[ti];
                    patch_n += (detail::point_as_vector(surf.vertices[tc[1]]) -
                                detail::point_as_vector(surf.vertices[tc[0]]))
                                 .cross(detail::point_as_vector(surf.vertices[tc[2]]) -
                                        detail::point_as_vector(surf.vertices[tc[0]]));
                }
                if (patch_n.squaredNorm() <= 1e-30)
                {
                    continue;
                }
                patch_n.normalize();

                bool geom_ok      = true;
                int kept_fan_tris = 0;
                for (const Index ti : v_tris[u])
                {
                    if (ti == del_t0 || ti == del_t1)
                    {
                        continue;
                    }
                    ++kept_fan_tris;
                    auto nc = surf.triangles[ti];
                    const Eigen::Vector3d old_n =
                      (detail::point_as_vector(surf.vertices[nc[1]]) -
                       detail::point_as_vector(surf.vertices[nc[0]]))
                        .cross(detail::point_as_vector(surf.vertices[nc[2]]) -
                               detail::point_as_vector(surf.vertices[nc[0]]));
                    for (int c = 0; c < 3; ++c)
                    {
                        if (nc[c] == u)
                        {
                            nc[c] = v;
                        }
                    }
                    const Eigen::Vector3d new_n =
                      (detail::point_as_vector(surf.vertices[nc[1]]) -
                       detail::point_as_vector(surf.vertices[nc[0]]))
                        .cross(detail::point_as_vector(surf.vertices[nc[2]]) -
                               detail::point_as_vector(surf.vertices[nc[0]]));
                    const double new_a = new_n.norm();
                    if (new_a <= 1e-13 * mean_ring_len * mean_ring_len)
                    {
                        geom_ok = false;
                        break;
                    }
                    const Eigen::Vector3d new_nh = new_n / new_a;
                    if (new_nh.dot(patch_n) < 0.75)
                    {
                        geom_ok = false;
                        break;
                    }
                    const double old_a = old_n.norm();
                    if (old_a > 1e-24 && new_nh.dot(old_n / old_a) < 0.30)
                    {
                        geom_ok = false;
                        break;
                    }
                }
                if (!geom_ok || kept_fan_tris == 0)
                {
                    continue;
                }

                // Accept collapse u -> v
                collapse_target[u]  = v;
                tri_deleted[del_t0] = true;
                tri_deleted[del_t1] = true;
                vert_busy[u]        = true;
                vert_busy[v]        = true;
                for (const Index nb : v_nbrs[u])
                {
                    vert_busy[nb] = true;
                }
                ++collapsed;
            }

            if (collapsed == 0)
            {
                break;
            }

            // Rewrite surviving triangles and compact unused vertices
            std::vector<std::array<Index, 3>> next_tris;
            std::vector<Marker> next_markers;
            next_tris.reserve(nt);
            next_markers.reserve(nt);
            std::vector<Index> old_to_new_v(nv, NOT_SET_Index);
            std::vector<Point3> next_verts;
            next_verts.reserve(nv);

            for (std::size_t ti = 0; ti < nt; ++ti)
            {
                if (tri_deleted[ti])
                {
                    continue;
                }
                auto tc = surf.triangles[ti];
                for (int c = 0; c < 3; ++c)
                {
                    Index cur = tc[c];
                    while (collapse_target[cur] != cur)
                    {
                        cur = collapse_target[cur];
                    }
                    if (old_to_new_v[cur] == NOT_SET_Index)
                    {
                        old_to_new_v[cur] = static_cast<Index>(next_verts.size());
                        next_verts.push_back(surf.vertices[cur]);
                    }
                    tc[c] = old_to_new_v[cur];
                }
                next_tris.push_back(tc);
                next_markers.push_back(surf.tri_to_marker_map[ti]);
            }

            surf.vertices          = std::move(next_verts);
            surf.triangles         = std::move(next_tris);
            surf.tri_to_marker_map = std::move(next_markers);
            update_edges_and_maps(surf);
        }
    };

    auto run_medial_taubin_polish = [&](Triangulation& surf, int n_cycles)
    {
        const std::size_t nv = surf.vertices.size();
        const std::size_t nt = surf.triangles.size();

        std::vector<bool> is_locked(nv, false);
        for (std::size_t eid = 0; eid < surf.edges.size(); ++eid)
        {
            const auto& inc     = surf.edge_to_tris_map[eid];
            bool interior_to_mm = false;
            if (inc.size() == 2)
            {
                const Marker m0 = surf.tri_to_marker_map[inc[0].first];
                const Marker m1 = surf.tri_to_marker_map[inc[1].first];
                if (m0 == m1 && m0 < is_medial_marker.size() && is_medial_marker[m0])
                {
                    interior_to_mm = true;
                }
            }
            if (!interior_to_mm)
            {
                is_locked[surf.edges[eid][0]] = true;
                is_locked[surf.edges[eid][1]] = true;
            }
        }

        std::vector<std::vector<Index>> v_tris(nv);
        std::vector<std::vector<Index>> v_nbrs(nv);
        std::vector<Index> free_med_verts;
        for (std::size_t ti = 0; ti < nt; ++ti)
        {
            const Marker m = surf.tri_to_marker_map[ti];
            if (m >= is_medial_marker.size() || !is_medial_marker[m])
            {
                continue;
            }
            const auto& tc = surf.triangles[ti];
            for (int c = 0; c < 3; ++c)
            {
                if (!is_locked[tc[c]])
                {
                    v_tris[tc[c]].push_back(static_cast<Index>(ti));
                    v_nbrs[tc[c]].push_back(tc[(c + 1) % 3]);
                    v_nbrs[tc[c]].push_back(tc[(c + 2) % 3]);
                }
            }
        }
        for (std::size_t v = 0; v < nv; ++v)
        {
            if (!is_locked[v] && !v_nbrs[v].empty())
            {
                auto& nb = v_nbrs[v];
                std::sort(nb.begin(), nb.end());
                nb.erase(std::unique(nb.begin(), nb.end()), nb.end());
                if (nb.size() >= 3)
                {
                    free_med_verts.push_back(static_cast<Index>(v));
                }
            }
        }
        if (free_med_verts.empty())
        {
            return;
        }

        std::vector<Point3> targets(free_med_verts.size());
        constexpr double kLambda = +0.45;
        constexpr double kMu     = -0.472;

        for (int cyc = 0; cyc < n_cycles; ++cyc)
        {
            for (const double step_w : { kLambda, kMu })
            {
                for (std::size_t i = 0; i < free_med_verts.size(); ++i)
                {
                    const Index v  = free_med_verts[i];
                    const auto& nb = v_nbrs[v];
                    Point3 avg     = { 0.0, 0.0, 0.0 };
                    for (const Index u : nb)
                    {
                        avg[0] += surf.vertices[u][0];
                        avg[1] += surf.vertices[u][1];
                        avg[2] += surf.vertices[u][2];
                    }
                    const double inv = 1.0 / static_cast<double>(nb.size());
                    targets[i]       = { avg[0] * inv, avg[1] * inv, avg[2] * inv };
                }

                for (std::size_t i = 0; i < free_med_verts.size(); ++i)
                {
                    const Index v      = free_med_verts[i];
                    const Point3 old_p = surf.vertices[v];
                    const double dx    = targets[i][0] - old_p[0];
                    const double dy    = targets[i][1] - old_p[1];
                    const double dz    = targets[i][2] - old_p[2];
                    if (dx * dx + dy * dy + dz * dz <= 1e-30)
                    {
                        continue;
                    }

                    std::vector<Eigen::Vector3d> old_normals;
                    old_normals.reserve(v_tris[v].size());
                    double min_q_old = 1.0;
                    bool valid_old   = true;
                    for (const Index ti : v_tris[v])
                    {
                        const auto& tc = surf.triangles[ti];
                        const Eigen::Vector3d n_raw =
                          (detail::point_as_vector(surf.vertices[tc[1]]) -
                           detail::point_as_vector(surf.vertices[tc[0]]))
                            .cross(detail::point_as_vector(surf.vertices[tc[2]]) -
                                   detail::point_as_vector(surf.vertices[tc[0]]));
                        const double a = n_raw.norm();
                        if (a <= 1e-24)
                        {
                            valid_old = false;
                            break;
                        }
                        old_normals.push_back(n_raw / a);
                        min_q_old = std::min(min_q_old,
                                             compute_tri_quality_2d(surf.vertices[tc[0]],
                                                                    surf.vertices[tc[1]],
                                                                    surf.vertices[tc[2]]));
                    }
                    if (!valid_old)
                    {
                        continue;
                    }

                    bool accepted = false;
                    for (const double backoff : { 1.0, 0.5, 0.25 })
                    {
                        const double w   = step_w * backoff;
                        surf.vertices[v] = { old_p[0] + w * dx,
                                             old_p[1] + w * dy,
                                             old_p[2] + w * dz };
                        bool ok          = true;
                        for (std::size_t k = 0; k < v_tris[v].size(); ++k)
                        {
                            const auto& tc = surf.triangles[v_tris[v][k]];
                            const Eigen::Vector3d n_raw =
                              (detail::point_as_vector(surf.vertices[tc[1]]) -
                               detail::point_as_vector(surf.vertices[tc[0]]))
                                .cross(detail::point_as_vector(surf.vertices[tc[2]]) -
                                       detail::point_as_vector(surf.vertices[tc[0]]));
                            const double a = n_raw.norm();
                            if (a <= 1e-24 || (n_raw / a).dot(old_normals[k]) < 0.88)
                            {
                                ok = false;
                                break;
                            }
                            const double q_new = compute_tri_quality_2d(surf.vertices[tc[0]],
                                                                        surf.vertices[tc[1]],
                                                                        surf.vertices[tc[2]]);
                            if (q_new < std::min(0.18, 0.78 * min_q_old))
                            {
                                ok = false;
                                break;
                            }
                        }
                        if (ok)
                        {
                            accepted = true;
                            break;
                        }
                    }
                    if (!accepted)
                    {
                        surf.vertices[v] = old_p;
                    }
                }
            }
        }
    };

    run_medial_short_edge_collapses(rt.t, 3);
    run_medial_edge_flips(rt.t, 5);
    run_medial_taubin_polish(rt.t, 12);
    run_medial_edge_flips(rt.t, 3);

    if (detect_surface_self_intersections(rt.t).rows() > 0)
    {
        // Retry without collapses so vertex indices match `safe_backup` for per-patch fallback
        rt.t = safe_backup;
        run_medial_edge_flips(rt.t, 5);
        run_medial_taubin_polish(rt.t, 10);
        if (detect_surface_self_intersections(rt.t).rows() > 0)
        {
            rt.t = safe_backup;
        }
    }
}

void remove_interpenetrations(const CombinedSurfaces& inconsistent_combined_surfaces,
                              const ConsistentSurfaces& consistent_surfaces,
                              SegmentedTetrahedralization& st)
{
    if (st.object_region_info.regions_of_overlap.empty())
    {
        GlobalOptions::get().log_info("* No interpenetrations");
        return;
    }

    GlobalOptions::get().log_info("* remove_interpenetrations");
    GlobalOptions::get().log_info(" - Initial # of overlap regions: " +
                                  std::to_string(st.object_region_info.regions_of_overlap.size()));

    refine_overlap_interiors(st);

    const auto overlap_groups = find_overlap_region_groups(
      st.region_boundaries_info,
      st.object_region_info.regions_of_overlap);

    GlobalOptions::get().log_info(" - Discovered " + std::to_string(overlap_groups.size()) +
                                  " groups of contiguous overlap regions");

    std::vector<GroupHarmonicField> group_fields;
    group_fields.reserve(overlap_groups.size());
    for (const auto& overlap_group : overlap_groups)
    {
        group_fields.push_back(compute_group_harmonic_field(st, overlap_group));
    }

    subdivide_and_partition_overlaps(inconsistent_combined_surfaces,
                                     consistent_surfaces,
                                     overlap_groups,
                                     group_fields,
                                     st);

    auto region_boundary_result    = segment_region_boundaries(inconsistent_combined_surfaces,
                                                            consistent_surfaces,
                                                            st.tetrahedralization,
                                                            st.region_info,
                                                            st.object_region_info);
    st.region_boundaries_info      = std::move(region_boundary_result.region_boundaries_info);
    st.object_region_boundary_info = std::move(region_boundary_result.object_region_boundary_info);

    st.seams_info = segment_seams(st.tetrahedralization, st.region_boundaries_info);

    smooth_medial_surfaces(consistent_surfaces, st);

    if (GlobalOptions::get().write_intermediate_meshes)
    {
        write_features_to_vtk(st, "remove_interpenetrations_");
    }
}

} // namespace mesh_arrange
