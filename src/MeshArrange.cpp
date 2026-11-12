#include "mesh_arrange/MeshArrange.h"

#include "detail/MathHelpers.h"
#include "mesh_arrange/GlobalOptions.h"

#include <algorithm>
#include <cassert>
#include <set>
#include <unordered_map>
#include <utility>
#include <vector>

namespace mesh_arrange
{

namespace
{

struct ScopedGlobalOptions
{
    GlobalOptions saved;

    explicit ScopedGlobalOptions(const ArrangeOptions& opts)
      : saved(GlobalOptions::get())
    {
        auto& g                     = GlobalOptions::get();
        g.verbose                   = opts.verbose;
        g.log_callback              = opts.log_callback;
        g.check_common_errors       = opts.check_topology;
        g.write_intermediate_meshes = opts.write_debug_vtk;
        g.output_dir                = opts.debug_output_dir;
    }

    ~ScopedGlobalOptions()
    {
        GlobalOptions::get() = saved;
    }
};

/// Orient every tracked triangle in `st.tetrahedralization` so its right-hand normal points
/// strictly OUTWARD from `st.region_info.tri_to_bounded_regions_map[tr].first`.
void orient_segmented_surface_triangles_from_tets(SegmentedTetrahedralization& st)
{
    auto& tet_mesh = st.tetrahedralization;
    bool any_flip  = false;

    for (std::size_t tr = 0; tr < tet_mesh.triangles.size(); ++tr)
    {
        const auto& br = st.region_info.tri_to_bounded_regions_map[tr];
        if (br.first == NOT_SET_RegionID)
        {
            continue;
        }

        const auto& tt = tet_mesh.tri_to_tet_map[tr];
        Index tet_r0   = NOT_SET_Index;
        if (tt[0] != NOT_SET_Index && st.region_info.tet_to_region_map[tt[0]] == br.first)
        {
            tet_r0 = tt[0];
        }
        else if (tt[1] != NOT_SET_Index && st.region_info.tet_to_region_map[tt[1]] == br.first)
        {
            tet_r0 = tt[1];
        }
        if (tet_r0 == NOT_SET_Index)
        {
            continue;
        }

        auto& tri            = tet_mesh.triangles[tr];
        const auto& tet_conn = tet_mesh.tetrahedra[tet_r0];
        Index v_opp          = NOT_SET_Index;
        for (const Index tv : tet_conn)
        {
            if (tv != tri[0] && tv != tri[1] && tv != tri[2])
            {
                v_opp = tv;
                break;
            }
        }
        if (v_opp == NOT_SET_Index)
        {
            continue;
        }

        // Signed volume (p1 - p0) . ((p2 - p0) x (p_opp - p0)) / 6
        // equals ((p1 - p0) x (p2 - p0)) . (p_opp - p0) / 6.
        // Since (p_opp - p0) points INTO tet_r0, negative signed volume means
        // right-hand normal (p1 - p0) x (p2 - p0) points OUTWARD from tet_r0.
        const double sv = detail::tetrahedron_signed_volume(tet_mesh.vertices[tri[0]],
                                                            tet_mesh.vertices[tri[1]],
                                                            tet_mesh.vertices[tri[2]],
                                                            tet_mesh.vertices[v_opp]);
        if (sv > 0.0)
        {
            std::swap(tri[1], tri[2]);
            any_flip = true;
        }
    }

    if (any_flip)
    {
        update_edges_and_maps(tet_mesh);
    }
}

void extract_constituents_with_ids(const RecoveredTriangulation& rt,
                                   std::vector<Triangulation>& out_surfaces,
                                   std::vector<Index>& out_object_ids)
{
    out_surfaces = separate_objects(rt);
    out_object_ids.clear();
    out_object_ids.reserve(out_surfaces.size());

    std::vector<std::size_t> marker_tri_counts(rt.t.n_markers, 0);
    for (const Marker m : rt.t.tri_to_marker_map)
    {
        if (m < marker_tri_counts.size())
        {
            marker_tri_counts[m]++;
        }
    }

    for (std::size_t object_id = 0; object_id < rt.object_to_markers_map.size(); ++object_id)
    {
        std::size_t total_tris = 0;
        for (const Marker m : rt.object_to_markers_map[object_id])
        {
            if (m < marker_tri_counts.size())
            {
                total_tris += marker_tri_counts[m];
            }
        }
        if (total_tris > 0)
        {
            out_object_ids.push_back(static_cast<Index>(object_id));
        }
    }
    assert(out_object_ids.size() == out_surfaces.size());
}

void build_brep_topology(const SegmentedTetrahedralization& st, ArrangeResult& res)
{
    const auto& surf     = res.surface_mesh;
    const auto seam_info = segment_seams(surf);

    // -------------------------------------------------------------------------
    // 1. Decompose 1D seams into 0D BRepVertex corners and 1D BRepCurve arcs
    // -------------------------------------------------------------------------
    std::unordered_map<Index, std::vector<Index>> vert_to_seam_edges;
    for (const Index eid : seam_info.edges_on_seams)
    {
        const auto& e = surf.edges[eid];
        vert_to_seam_edges[e[0]].push_back(eid);
        vert_to_seam_edges[e[1]].push_back(eid);
    }

    std::vector<Index> sorted_seam_verts;
    sorted_seam_verts.reserve(vert_to_seam_edges.size());
    for (auto& kv : vert_to_seam_edges)
    {
        std::sort(kv.second.begin(), kv.second.end());
        sorted_seam_verts.push_back(kv.first);
    }
    std::sort(sorted_seam_verts.begin(), sorted_seam_verts.end());

    std::unordered_map<Index, Index> global_vid_to_brep_vid;
    for (const Index vid : sorted_seam_verts)
    {
        const auto& inc_edges = vert_to_seam_edges.at(vid);
        bool is_corner        = (inc_edges.size() != 2);
        if (!is_corner)
        {
            const SeamID s0 = seam_info.edge_id_to_seam_id_map[inc_edges[0]];
            const SeamID s1 = seam_info.edge_id_to_seam_id_map[inc_edges[1]];
            if (s0 != s1)
            {
                is_corner = true;
            }
        }
        if (is_corner)
        {
            const Index bvid            = static_cast<Index>(res.brep_vertices.size());
            global_vid_to_brep_vid[vid] = bvid;
            BRepVertex bv;
            bv.global_vertex_id = vid;
            bv.point            = surf.vertices[vid];
            res.brep_vertices.push_back(std::move(bv));
        }
    }

    // Map each seam edge `eid` -> {curve_id, +1 if curve traverses edge[0]->edge[1], -1 otherwise}
    std::unordered_map<Index, std::pair<Index, int>> edge_to_curve_dir;
    edge_to_curve_dir.reserve(seam_info.edges_on_seams.size());
    std::vector<bool> edge_visited(surf.edges.size(), false);

    auto collect_curve_faces = [&](const std::vector<Index>& curve_edges)
    {
        std::set<RegionBoundaryID> faces;
        for (const Index eid : curve_edges)
        {
            for (const auto& tp : surf.edge_to_tris_map[eid])
            {
                faces.insert(surf.tri_to_marker_map[tp.first]);
            }
        }
        return std::vector<RegionBoundaryID>(faces.begin(), faces.end());
    };

    auto walk_chain = [&](Index start_vid, Index first_eid)
    {
        const Index curve_id = static_cast<Index>(res.brep_curves.size());
        const SeamID sid     = seam_info.edge_id_to_seam_id_map[first_eid];

        std::vector<Index> poly = { start_vid };
        std::vector<Index> used_edges;

        Index curr_v = start_vid;
        Index curr_e = first_eid;
        while (true)
        {
            edge_visited[curr_e] = true;
            used_edges.push_back(curr_e);

            const auto& ed      = surf.edges[curr_e];
            const Index next_v  = (ed[0] == curr_v) ? ed[1] : ed[0];
            const int edge_sign = (ed[0] == curr_v) ? +1 : -1;
            edge_to_curve_dir.emplace(curr_e, std::make_pair(curve_id, edge_sign));
            poly.push_back(next_v);

            if (global_vid_to_brep_vid.count(next_v) > 0)
            {
                break;
            }

            const auto& next_inc = vert_to_seam_edges.at(next_v);
            Index next_e         = NOT_SET_Index;
            for (const Index cand_e : next_inc)
            {
                if (!edge_visited[cand_e] && seam_info.edge_id_to_seam_id_map[cand_e] == sid)
                {
                    next_e = cand_e;
                    break;
                }
            }
            if (next_e == NOT_SET_Index)
            {
                // Safety fallback: promote unexpected endpoint to BRepVertex
                const Index bvid               = static_cast<Index>(res.brep_vertices.size());
                global_vid_to_brep_vid[next_v] = bvid;
                BRepVertex bv;
                bv.global_vertex_id = next_v;
                bv.point            = surf.vertices[next_v];
                res.brep_vertices.push_back(std::move(bv));
                break;
            }
            curr_v = next_v;
            curr_e = next_e;
        }

        BRepCurve bc;
        bc.seam_id           = sid;
        bc.start_brep_vertex = global_vid_to_brep_vid.at(poly.front());
        bc.end_brep_vertex   = global_vid_to_brep_vid.at(poly.back());
        bc.is_closed_loop    = (poly.front() == poly.back());
        bc.incident_faces    = collect_curve_faces(used_edges);
        bc.polyline_vertices = std::move(poly);

        res.brep_vertices[bc.start_brep_vertex].incident_curves.push_back(curve_id);
        if (bc.end_brep_vertex != bc.start_brep_vertex)
        {
            res.brep_vertices[bc.end_brep_vertex].incident_curves.push_back(curve_id);
        }
        res.brep_curves.push_back(std::move(bc));
    };

    // Walk all open (or corner-anchored closed) curves starting at BRepVertex corners
    const std::size_t initial_n_corners = res.brep_vertices.size();
    for (std::size_t bi = 0; bi < initial_n_corners; ++bi)
    {
        const Index start_vid = res.brep_vertices[bi].global_vertex_id;
        const auto inc_copy   = vert_to_seam_edges.at(start_vid);
        for (const Index eid : inc_copy)
        {
            if (!edge_visited[eid])
            {
                walk_chain(start_vid, eid);
            }
        }
    }

    // Walk any remaining pure periodic closed loops (no degree!=2 corner) by anchoring at min vertex ID
    for (const Index vid : sorted_seam_verts)
    {
        const auto& inc = vert_to_seam_edges.at(vid);
        Index start_eid = NOT_SET_Index;
        for (const Index eid : inc)
        {
            if (!edge_visited[eid])
            {
                start_eid = eid;
                break;
            }
        }
        if (start_eid == NOT_SET_Index)
        {
            continue;
        }

        if (global_vid_to_brep_vid.count(vid) == 0)
        {
            const Index bvid            = static_cast<Index>(res.brep_vertices.size());
            global_vid_to_brep_vid[vid] = bvid;
            BRepVertex bv;
            bv.global_vertex_id = vid;
            bv.point            = surf.vertices[vid];
            res.brep_vertices.push_back(std::move(bv));
        }
        walk_chain(vid, start_eid);
    }

    // -------------------------------------------------------------------------
    // 2. Build 2D BRepFace patches (0 .. n_region_boundaries - 1)
    // -------------------------------------------------------------------------
    const RegionBoundaryID n_faces = st.region_boundaries_info.n_region_boundaries;
    res.brep_faces.resize(n_faces);

    for (RegionBoundaryID bid = 0; bid < n_faces; ++bid)
    {
        auto& bf       = res.brep_faces[bid];
        bf.boundary_id = bid;

        const auto& st_tris = st.region_boundaries_info.region_boundary_to_triangles_map[bid];
        if (!st_tris.empty())
        {
            bf.bounded_regions = st.region_info.tri_to_bounded_regions_map[st_tris[0]];
        }

        if (bid < st.object_region_boundary_info.boundary_id_to_objects_map.size())
        {
            const auto& objs = st.object_region_boundary_info.boundary_id_to_objects_map[bid];
            bf.origin_objects.assign(objs.begin(), objs.end());
        }
    }

    for (std::size_t tid = 0; tid < surf.triangles.size(); ++tid)
    {
        const RegionBoundaryID bid = surf.tri_to_marker_map[tid];
        if (bid < n_faces)
        {
            res.brep_faces[bid].triangles.push_back(surf.triangles[tid]);
            res.brep_faces[bid].surface_tri_ids.push_back(static_cast<Index>(tid));
        }
    }

    // Compute signed bounding curves and embedded curves/vertices for each BRepFace
    for (RegionBoundaryID bid = 0; bid < n_faces; ++bid)
    {
        auto& bf = res.brep_faces[bid];
        std::set<std::pair<Index, int>> signed_curves_set;
        std::set<Index> bound_curve_ids;
        std::set<Index> embedded_curves_set;
        std::set<Index> face_global_verts;

        for (std::size_t local_ti = 0; local_ti < bf.triangles.size(); ++local_ti)
        {
            const Index tid  = bf.surface_tri_ids[local_ti];
            const auto& conn = bf.triangles[local_ti];
            face_global_verts.insert(conn[0]);
            face_global_verts.insert(conn[1]);
            face_global_verts.insert(conn[2]);

            for (int le = 0; le < 3; ++le)
            {
                const Index eid = surf.tri_to_edges_map[tid][le];
                const auto it_c = edge_to_curve_dir.find(eid);
                if (it_c == edge_to_curve_dir.end())
                {
                    continue;
                }

                int same_face_deg = 0;
                for (const auto& tp : surf.edge_to_tris_map[eid])
                {
                    if (surf.tri_to_marker_map[tp.first] == bid)
                    {
                        ++same_face_deg;
                    }
                }
                if (same_face_deg == 1)
                {
                    const Index hu      = conn[le];
                    const Index hv      = conn[(le + 1) % 3];
                    const auto& ed      = surf.edges[eid];
                    const int he_sign   = (hu == ed[0] && hv == ed[1]) ? +1 : -1;
                    const int curve_dir = he_sign * it_c->second.second;
                    signed_curves_set.emplace(it_c->second.first, curve_dir);
                    bound_curve_ids.insert(it_c->second.first);
                }
                else if (same_face_deg >= 2)
                {
                    embedded_curves_set.insert(it_c->second.first);
                }
            }
        }

        for (const Index cid : bound_curve_ids)
        {
            embedded_curves_set.erase(cid);
        }

        bf.signed_bounding_curves.assign(signed_curves_set.begin(), signed_curves_set.end());
        bf.embedded_curves.assign(embedded_curves_set.begin(), embedded_curves_set.end());

        // Identify any isolated BRepVertex lying on this face with no curve attached to this face
        std::set<Index> curve_endpoint_bvs;
        for (const Index cid : bound_curve_ids)
        {
            const auto& bc = res.brep_curves[cid];
            if (bc.start_brep_vertex != NOT_SET_Index)
            {
                curve_endpoint_bvs.insert(bc.start_brep_vertex);
            }
            if (bc.end_brep_vertex != NOT_SET_Index)
            {
                curve_endpoint_bvs.insert(bc.end_brep_vertex);
            }
        }
        for (const Index cid : embedded_curves_set)
        {
            const auto& bc = res.brep_curves[cid];
            if (bc.start_brep_vertex != NOT_SET_Index)
            {
                curve_endpoint_bvs.insert(bc.start_brep_vertex);
            }
            if (bc.end_brep_vertex != NOT_SET_Index)
            {
                curve_endpoint_bvs.insert(bc.end_brep_vertex);
            }
        }
        for (std::size_t bvi = 0; bvi < res.brep_vertices.size(); ++bvi)
        {
            const Index idx = static_cast<Index>(bvi);
            if (face_global_verts.count(res.brep_vertices[bvi].global_vertex_id) > 0 &&
                curve_endpoint_bvs.count(idx) == 0)
            {
                bf.embedded_vertices.push_back(idx);
            }
        }
    }

    // -------------------------------------------------------------------------
    // 3. Build 3D BRepVolume regions (0 .. n_regions - 1)
    // -------------------------------------------------------------------------
    const RegionID n_regions = st.region_info.n_regions;
    res.brep_volumes.resize(n_regions);

    for (RegionID r = 0; r < n_regions; ++r)
    {
        auto& bvol     = res.brep_volumes[r];
        bvol.region_id = r;

        if (r < st.object_region_info.region_to_objects_map.size())
        {
            const auto& owners = st.object_region_info.region_to_objects_map[r];
            bvol.owner_objects.assign(owners.begin(), owners.end());
            std::sort(bvol.owner_objects.begin(), bvol.owner_objects.end());
        }
        bvol.is_overlap = (bvol.owner_objects.size() > 1);

        if (r < st.region_boundaries_info.region_id_to_region_boundary_ids_map.size())
        {
            for (const RegionBoundaryID bid :
                 st.region_boundaries_info.region_id_to_region_boundary_ids_map[r])
            {
                const auto& br = res.brep_faces[bid].bounded_regions;
                if (br.first == r && br.second != r)
                {
                    bvol.signed_bounding_faces.emplace_back(bid, +1);
                }
                else if (br.second == r && br.first != r)
                {
                    bvol.signed_bounding_faces.emplace_back(bid, -1);
                }
            }
        }
    }
}

} // namespace

ArrangeResult arrange_surfaces(const std::vector<InputSurface>& surfaces,
                               const ArrangeOptions& options)
{
    ScopedGlobalOptions scoped_opts(options);

    ArrangeResult result;
    if (surfaces.empty())
    {
        return result;
    }

    const CombinedSurfaces combined = combine_surfaces(surfaces);
    if (combined.inconsistent_surface_mesh.triangles.empty())
    {
        return result;
    }

    ConsistentSurfaces consistent = create_consistent_surface(combined);
    SegmentedTetrahedralization st =
      identify_features(combined, consistent, options.plc_tetrahedralizer);

    if (st.region_info.n_regions == 0 || st.tetrahedralization.triangles.empty())
    {
        result.segmented_tet = std::move(st);
        return result;
    }

    if (options.remove_overlaps && !st.object_region_info.regions_of_overlap.empty())
    {
        remove_interpenetrations(combined, consistent, st);
    }

    if (!options.optimize_surfaces)
    {
        if (options.split_constrained_seam_edges)
        {
            split_constrained_edges(st);
        }
        if (options.max_edge_length > 0.0)
        {
            split_large_edges(options.max_edge_length, st);
        }
    }

    // Orient all surface triangles outward from bounded_regions.first before recovery
    orient_segmented_surface_triangles_from_tets(st);

    RecoveredTriangulation rt = recover_triangulation(st);
    if (options.optimize_surfaces)
    {
        optimize_recovered_triangulation(rt);
        if (options.split_constrained_seam_edges)
        {
            auto si = segment_seams(rt.t);
            split_constrained_edges(rt.t, si);
        }
        if (options.max_edge_length > 0.0)
        {
            auto si = segment_seams(rt.t);
            split_large_edges(options.max_edge_length, rt.t, si);
        }
    }

    extract_constituents_with_ids(rt, result.constituent_surfaces, result.constituent_object_ids);
    result.exterior_boundary = find_global_boundary(rt);
    result.surface_mesh      = std::move(rt.t);

    build_brep_topology(st, result);
    result.segmented_tet = std::move(st);

    return result;
}

} // namespace mesh_arrange
