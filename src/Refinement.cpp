#include "mesh_arrange/Refinement.h"

#include "detail/MathHelpers.h"
#include "mesh_arrange/GlobalOptions.h"
#include "mesh_arrange/SplitEdges.h"

#include <iostream>
#include <set>
#include <vector>

namespace mesh_arrange
{

namespace
{

void update_features_for_edge_splits(const UpdateResult& update_result,
                                     SegmentedTetrahedralization& st)
{
    // Update tetrahedron region mappings
    st.region_info.tet_to_region_map.resize(st.tetrahedralization.tetrahedra.size(),
                                            NOT_SET_RegionID);
    for (const auto& update : update_result.updated_tetrahedra)
    {
        const RegionID orig_region_id = st.region_info.tet_to_region_map[update.orig_tet_id];
        for (const Index new_tet_id : update.new_tet_ids)
        {
            st.region_info.region_to_tets_map[orig_region_id].push_back(new_tet_id);
            st.region_info.tet_to_region_map[new_tet_id] = orig_region_id;
        }
    }

    // Update triangle bounded-region and facet mappings
    st.region_info.tri_to_bounded_regions_map.resize(st.tetrahedralization.triangles.size());
    st.region_info.tri_to_facet_map.resize(st.tetrahedralization.triangles.size());
    for (const auto& update : update_result.updated_triangles)
    {
        const auto orig_bounded_region_ids =
          st.region_info.tri_to_bounded_regions_map[update.orig_tri_id];
        const FacetID orig_facet_id = st.region_info.tri_to_facet_map[update.orig_tri_id];
        for (const Index new_tri_id : update.new_tri_ids)
        {
            st.region_info.tri_to_bounded_regions_map[new_tri_id] = orig_bounded_region_ids;
            if (orig_bounded_region_ids.first != NOT_SET_RegionID)
            {
                st.region_info.region_to_bounding_triangles_map[orig_bounded_region_ids.first]
                  .push_back(new_tri_id);
            }
            if (orig_bounded_region_ids.second != NOT_SET_RegionID)
            {
                st.region_info.region_to_bounding_triangles_map[orig_bounded_region_ids.second]
                  .push_back(new_tri_id);
            }
            st.region_info.tri_to_facet_map[new_tri_id] = orig_facet_id;
        }
    }

    // Update region boundary info for subdivided triangles
    st.region_boundaries_info.triangle_to_boundary_id_map.resize(
      st.tetrahedralization.triangles.size(),
      NOT_SET_RegionBoundaryID);
    for (const auto& t : update_result.updated_triangles)
    {
        const RegionBoundaryID orig_boundary_id =
          st.region_boundaries_info.triangle_to_boundary_id_map[t.orig_tri_id];
        for (const Index new_tri_id : t.new_tri_ids)
        {
            st.region_boundaries_info.triangle_to_boundary_id_map[new_tri_id] = orig_boundary_id;
            st.region_boundaries_info.region_boundary_to_triangles_map[orig_boundary_id].push_back(
              new_tri_id);
        }
    }

    // Update seam edge info
    st.seams_info.edge_id_to_seam_id_map.resize(st.tetrahedralization.edges.size(), NOT_SET_SeamID);
    for (const auto& update : update_result.updated_edges)
    {
        if (update.update_type == UpdatedEdgeInfo::SUBDIVIDED_EDGE &&
            st.seams_info.edge_id_to_seam_id_map[update.orig_edge_id] != NOT_SET_SeamID)
        {
            st.seams_info.edges_on_seams.push_back(update.new_edge_id);
            st.seams_info.edge_id_to_seam_id_map[update.new_edge_id] =
              st.seams_info.edge_id_to_seam_id_map[update.orig_edge_id];
        }
    }
}

} // namespace

void split_edges(const std::vector<Index>& edges_to_split, SegmentedTetrahedralization& st)
{
    if (edges_to_split.empty())
    {
        return;
    }
    const auto update_result = split_edges(edges_to_split, st.tetrahedralization);
    update_features_for_edge_splits(update_result, st);
}

void split_edges(const std::vector<Index>& edges_to_split,
                 const std::vector<Index>& existing_vertices_for_split_edges,
                 SegmentedTetrahedralization& st)
{
    if (edges_to_split.empty())
    {
        return;
    }
    const auto update_result =
      split_edges(edges_to_split, existing_vertices_for_split_edges, st.tetrahedralization);
    update_features_for_edge_splits(update_result, st);
}

void split_constrained_edges(SegmentedTetrahedralization& st)
{
    auto& t = st.tetrahedralization;

    std::vector<bool> is_node_on_seam(t.vertices.size(), false);
    for (const Index edge_id : st.seams_info.edges_on_seams)
    {
        is_node_on_seam[t.edges[edge_id][0]] = true;
        is_node_on_seam[t.edges[edge_id][1]] = true;
    }

    std::vector<Index> edges_to_split;
    for (std::size_t edge_id = 0; edge_id < t.edges.size(); ++edge_id)
    {
        const Index v_1                  = t.edges[edge_id][0];
        const Index v_2                  = t.edges[edge_id][1];
        const bool edge_pts_constrained  = is_node_on_seam[v_1] && is_node_on_seam[v_2];
        if (edge_pts_constrained &&
            st.seams_info.edge_id_to_seam_id_map[edge_id] == NOT_SET_SeamID)
        {
            edges_to_split.push_back(static_cast<Index>(edge_id));
        }
    }

    split_edges(edges_to_split, st);
}

void split_overlap_edges(SegmentedTetrahedralization& st)
{
    auto& t = st.tetrahedralization;

    std::set<Index> edges_to_split_set;
    for (const RegionID region_id : st.object_region_info.regions_of_overlap)
    {
        for (const Index tet_id : st.region_info.region_to_tets_map[region_id])
        {
            for (const Index edge_id : t.tet_to_edges_map[tet_id])
            {
                edges_to_split_set.insert(edge_id);
            }
        }
    }

    GlobalOptions::get().log_info(" - Splitting " + std::to_string(edges_to_split_set.size()) +
                                  " edges in overlap regions.");

    const std::vector<Index> edges_to_split(edges_to_split_set.begin(), edges_to_split_set.end());
    split_edges(edges_to_split, st);
}

void split_large_edges(const double max_edge_length, SegmentedTetrahedralization& st)
{
    if (max_edge_length <= 0.0)
    {
        return;
    }

    GlobalOptions::get().log_info("Refining via splitting edges...");

    while (true)
    {
        std::vector<Index> edges_to_split;
        for (std::size_t i = 0; i < st.tetrahedralization.edges.size(); ++i)
        {
            const auto& edge       = st.tetrahedralization.edges[i];
            const double edge_len  = detail::distance(st.tetrahedralization.vertices[edge[0]],
                                                     st.tetrahedralization.vertices[edge[1]]);
            if (edge_len > max_edge_length)
            {
                edges_to_split.push_back(static_cast<Index>(i));
            }
        }

        if (edges_to_split.empty())
        {
            break;
        }

        GlobalOptions::get().log_info(" - Splitting " + std::to_string(edges_to_split.size()) +
                                      " edges.");
        split_edges(edges_to_split, st);
    }
}

} // namespace mesh_arrange
