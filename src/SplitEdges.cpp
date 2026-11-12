#include "mesh_arrange/SplitEdges.h"

#include "detail/ElementTables.h"
#include "detail/MathHelpers.h"
#include "mesh_arrange/GlobalOptions.h"

#include <algorithm>
#include <cassert>
#include <iostream>
#include <map>
#include <stdexcept>

namespace mesh_arrange
{

struct UpdateState
{
    std::vector<Index> orig_edge_to_update_index;
    std::vector<Index> orig_tri_to_update_index;
    std::vector<Index> orig_tet_to_update_index;

    std::vector<Index> current_to_orig_edge_map;
    std::vector<Index> current_to_orig_tri_map;
    std::vector<Index> current_to_orig_tet_map;

    UpdateResult result;
};

namespace
{

struct TriangleNewEdgeInfo
{
    Index new_triangle_id;
    Index new_edge_id;
    bool is_edge_aligned;
};

template<typename T>
std::map<Index, TriangleNewEdgeInfo> split_triangles_on_edge(const Index& orig_edge_id,
                                                             const Index& new_split_edge_id,
                                                             T& t)
{
    auto& triangles        = t.triangles;
    auto& edges            = t.edges;
    auto& edge_to_tris_map = t.edge_to_tris_map;
    auto& tri_to_edges_map = t.tri_to_edges_map;

    const Index n_1_id         = edges[orig_edge_id][0];
    const Index new_node_index = edges[orig_edge_id][1];

    edge_to_tris_map.push_back({});

    std::map<Index, TriangleNewEdgeInfo> tri_updates;

    const auto parent_infos = edge_to_tris_map[orig_edge_id];
    for (const auto& parent_info : parent_infos)
    {
        const Index parent_triangle_index    = parent_info.first;
        const ConnIndex split_edge_pos_in_t  = parent_info.second;
        const auto& split_edge_conn          = detail::kTriEdges[split_edge_pos_in_t];

        const Index t_id     = parent_triangle_index;
        const Index new_t_id = static_cast<Index>(triangles.size());
        triangles.push_back(triangles[t_id]);
        auto& t_conn     = triangles[t_id];
        auto& new_t_conn = triangles.back();

        tri_to_edges_map.push_back(tri_to_edges_map[t_id]);

        ConnIndex external_conn_index = 0;
        for (ConnIndex i = 0; i < 3; ++i)
        {
            if (i != split_edge_conn[0] && i != split_edge_conn[1])
            {
                external_conn_index = i;
            }
        }
        const Index external_v = t_conn[external_conn_index];

        t_conn[split_edge_conn[1]]     = new_node_index;
        new_t_conn[split_edge_conn[0]] = new_node_index;

        const ConnIndex split_half_edge_position_in_t     = split_edge_pos_in_t;
        const ConnIndex internal_edge_position_in_t       = (split_edge_pos_in_t + 1) % 3;
        const ConnIndex split_half_edge_position_in_new_t = split_edge_pos_in_t;
        const ConnIndex external_edge_position_in_new_t   = (split_edge_pos_in_t + 1) % 3;
        const ConnIndex internal_edge_position_in_new_t   = (split_edge_pos_in_t + 2) % 3;

        const bool edge_aligned = (t_conn[split_edge_conn[0]] == n_1_id);
        if (edge_aligned)
        {
            tri_to_edges_map[new_t_id][split_half_edge_position_in_new_t] = new_split_edge_id;
            edge_to_tris_map[new_split_edge_id].push_back(
              { new_t_id, split_half_edge_position_in_new_t });

            for (auto& tri_for_edge : edge_to_tris_map[orig_edge_id])
            {
                if (tri_for_edge.first == t_id)
                {
                    tri_for_edge.second = split_half_edge_position_in_t;
                }
            }
        }
        else
        {
            tri_to_edges_map[t_id][split_half_edge_position_in_t] = new_split_edge_id;
            edge_to_tris_map[new_split_edge_id].push_back(
              { t_id, split_half_edge_position_in_t });

            for (auto& tri_for_edge : edge_to_tris_map[orig_edge_id])
            {
                if (tri_for_edge.first == t_id)
                {
                    tri_for_edge.first  = new_t_id;
                    tri_for_edge.second = split_half_edge_position_in_new_t;
                }
            }
        }

        const Index new_internal_edge_id = static_cast<Index>(edge_to_tris_map.size());
        edges.push_back({ new_node_index, external_v });
        edge_to_tris_map.push_back(
          { { t_id, internal_edge_position_in_t }, { new_t_id, internal_edge_position_in_new_t } });
        tri_to_edges_map[t_id][internal_edge_position_in_t]         = new_internal_edge_id;
        tri_to_edges_map[new_t_id][internal_edge_position_in_new_t] = new_internal_edge_id;

        const Index external_edge_id = tri_to_edges_map[new_t_id][external_edge_position_in_new_t];
        for (auto& p_info : edge_to_tris_map[external_edge_id])
        {
            if (p_info.first == t_id)
            {
                p_info.first  = new_t_id;
                p_info.second = external_edge_position_in_new_t;
            }
        }

        tri_updates[t_id] = { new_t_id, new_internal_edge_id, edge_aligned };
    }

    assert(tri_to_edges_map.size() == triangles.size());
    return tri_updates;
}

struct TetrahedraUpdateResult
{
    struct NewInternalEdgeInfo
    {
        Index new_edge_id;
        Index new_node_id;
    };

    std::map<Index, Index> original_to_new_tet_map;
    std::vector<NewInternalEdgeInfo> new_internal_edges;
};

TetrahedraUpdateResult update_tetrahedra_on_split_edge(
  const Index& orig_edge_id,
  const Index& new_split_edge_id,
  const std::map<Index, TriangleNewEdgeInfo>& tri_updates,
  Tetrahedralization& t)
{
    auto& triangles        = t.triangles;
    auto& tets             = t.tetrahedra;
    auto& edges            = t.edges;
    auto& edge_to_tris_map = t.edge_to_tris_map;
    auto& edge_to_tets_map = t.edge_to_tets_map;
    auto& tet_to_edges_map = t.tet_to_edges_map;
    auto& tet_to_tri_map   = t.tet_to_tri_map;
    auto& tri_to_tet_map   = t.tri_to_tet_map;

    TetrahedraUpdateResult r;

    static const ConnIndex split_face_table[6][2] = {
        { 0, 3 }, { 0, 1 }, { 1, 3 }, { 0, 2 }, { 2, 3 }, { 1, 2 }
    };
    static const ConnIndex recycled_tet_internal_edge_replace_table[6][2] = {
        { 3, 4 }, { 3, 5 }, { 5, 4 }, { 1, 5 }, { 5, 2 }, { 2, 4 }
    };
    static const ConnIndex new_tet_internal_edge_replace_table[6][2] = {
        { 1, 2 }, { 0, 2 }, { 1, 0 }, { 0, 4 }, { 3, 0 }, { 1, 3 }
    };
    static auto& recycled_tet_edges_to_keep_table = new_tet_internal_edge_replace_table;
    static auto& new_tet_edges_to_keep_table      = recycled_tet_internal_edge_replace_table;
    static const ConnIndex recycled_tet_face_to_keep[6] = { 1, 3, 0, 3, 0, 0 };
    static const ConnIndex new_tet_face_to_keep[6]      = { 2, 2, 2, 1, 1, 3 };
    static const ConnIndex common_orig_edge_table[6]    = { 5, 4, 3, 2, 1, 0 };

    const Index n_1_id   = edges[orig_edge_id][0];
    const Index n_new_id = edges[orig_edge_id][1];
    const Index n_2_id   = edges[new_split_edge_id][1];

    const auto orig_parent_tets = edge_to_tets_map[orig_edge_id];
    edge_to_tets_map[orig_edge_id].clear();
    edge_to_tets_map.push_back({});

    tri_to_tet_map.resize(triangles.size());
    for (const auto& u : tri_updates)
    {
        tri_to_tet_map[u.second.new_triangle_id] = tri_to_tet_map[u.first];
        edge_to_tets_map.push_back({});
    }

    assert(tri_to_tet_map.size() == triangles.size());
    assert(edge_to_tets_map.size() == edges.size());

    std::map<std::array<Index, 3>, Index> untracked_face_to_new_edge_map;

    for (const auto& parent_info : orig_parent_tets)
    {
        const Index orig_tet_index             = parent_info.first;
        const ConnIndex split_edge_pos_in_tet  = parent_info.second;

        const Index new_tet_index = static_cast<Index>(tets.size());
        tets.push_back(tets[orig_tet_index]);
        auto& recycled_tet_conn  = tets[orig_tet_index];
        auto& new_tet_conn       = tets.back();
        const auto orig_tet_conn = tets[orig_tet_index];
        const auto orig_tet_edges = tet_to_edges_map[orig_tet_index];

        r.original_to_new_tet_map[orig_tet_index] = new_tet_index;

        const auto& split_edge_c = detail::kTetEdges[split_edge_pos_in_tet];
        assert(orig_tet_conn[split_edge_c[0]] == n_1_id || orig_tet_conn[split_edge_c[0]] == n_2_id);
        assert(orig_tet_conn[split_edge_c[1]] == n_1_id || orig_tet_conn[split_edge_c[1]] == n_2_id);

        const bool edge_aligned = (orig_tet_conn[split_edge_c[0]] == n_1_id);

        recycled_tet_conn[split_edge_c[1]] = n_new_id;
        new_tet_conn[split_edge_c[0]]      = n_new_id;

        const auto orig_tet_to_tris    = tet_to_tri_map[orig_tet_index];
        tet_to_tri_map[orig_tet_index] = { NOT_SET_Index, NOT_SET_Index, NOT_SET_Index, NOT_SET_Index };
        tet_to_tri_map.push_back({ NOT_SET_Index, NOT_SET_Index, NOT_SET_Index, NOT_SET_Index });

        tet_to_edges_map[orig_tet_index] = {
            NOT_SET_Index, NOT_SET_Index, NOT_SET_Index, NOT_SET_Index, NOT_SET_Index, NOT_SET_Index
        };
        tet_to_edges_map.push_back({
            NOT_SET_Index, NOT_SET_Index, NOT_SET_Index, NOT_SET_Index, NOT_SET_Index, NOT_SET_Index
        });

        const ConnIndex common_edge_position                   = common_orig_edge_table[split_edge_pos_in_tet];
        const Index common_edge                                = orig_tet_edges[common_edge_position];
        tet_to_edges_map[orig_tet_index][common_edge_position] = common_edge;
        tet_to_edges_map[new_tet_index][common_edge_position]  = common_edge;
        edge_to_tets_map[common_edge].push_back({ new_tet_index, common_edge_position });

        if (edge_aligned)
        {
            tet_to_edges_map[orig_tet_index][split_edge_pos_in_tet] = orig_edge_id;
            tet_to_edges_map[new_tet_index][split_edge_pos_in_tet]  = new_split_edge_id;
            edge_to_tets_map[orig_edge_id].push_back({ orig_tet_index, split_edge_pos_in_tet });
            edge_to_tets_map[new_split_edge_id].push_back({ new_tet_index, split_edge_pos_in_tet });
        }
        else
        {
            tet_to_edges_map[new_tet_index][split_edge_pos_in_tet]  = orig_edge_id;
            tet_to_edges_map[orig_tet_index][split_edge_pos_in_tet] = new_split_edge_id;
            edge_to_tets_map[orig_edge_id].push_back({ new_tet_index, split_edge_pos_in_tet });
            edge_to_tets_map[new_split_edge_id].push_back({ orig_tet_index, split_edge_pos_in_tet });
        }

        tet_to_tri_map[orig_tet_index][recycled_tet_face_to_keep[split_edge_pos_in_tet]] =
          orig_tet_to_tris[recycled_tet_face_to_keep[split_edge_pos_in_tet]];

        const Index face_for_new_tet = orig_tet_to_tris[new_tet_face_to_keep[split_edge_pos_in_tet]];
        tet_to_tri_map[new_tet_index][new_tet_face_to_keep[split_edge_pos_in_tet]] = face_for_new_tet;
        if (face_for_new_tet != NOT_SET_Index)
        {
            for (auto& tet_ref : tri_to_tet_map[face_for_new_tet])
            {
                if (tet_ref == orig_tet_index)
                {
                    tet_ref = new_tet_index;
                }
            }
        }

        for (const ConnIndex edge_position : recycled_tet_edges_to_keep_table[split_edge_pos_in_tet])
        {
            tet_to_edges_map[orig_tet_index][edge_position] = orig_tet_edges[edge_position];
        }
        for (const ConnIndex edge_position : new_tet_edges_to_keep_table[split_edge_pos_in_tet])
        {
            const Index e_id = orig_tet_edges[edge_position];
            for (auto& edge_to_tet_info : edge_to_tets_map[e_id])
            {
                if (edge_to_tet_info.first == orig_tet_index)
                {
                    edge_to_tet_info.first  = new_tet_index;
                    edge_to_tet_info.second = edge_position;
                }
            }
            tet_to_edges_map[new_tet_index][edge_position] = e_id;
        }

        for (int i = 0; i < 2; ++i)
        {
            const ConnIndex split_face_position = split_face_table[split_edge_pos_in_tet][i];
            Index new_internal_edge_index       = NOT_SET_Index;

            const Index orig_tri_that_was_split = orig_tet_to_tris[split_face_position];
            if (orig_tri_that_was_split == NOT_SET_Index)
            {
                std::array<Index, 3> tri_face = {
                    orig_tet_conn[detail::kTetFaces[split_face_position][0]],
                    orig_tet_conn[detail::kTetFaces[split_face_position][1]],
                    orig_tet_conn[detail::kTetFaces[split_face_position][2]]
                };
                std::sort(tri_face.begin(), tri_face.end());

                auto it = untracked_face_to_new_edge_map.find(tri_face);
                if (it == untracked_face_to_new_edge_map.end())
                {
                    Index third_node = NOT_SET_Index;
                    for (Index v : tri_face)
                    {
                        if (v != n_1_id && v != n_2_id)
                        {
                            third_node = v;
                        }
                    }
                    new_internal_edge_index = static_cast<Index>(edges.size());
                    edges.push_back({ third_node, n_new_id });
                    edge_to_tets_map.push_back({});
                    edge_to_tris_map.push_back({});
                    r.new_internal_edges.push_back({ new_internal_edge_index, n_new_id });
                    untracked_face_to_new_edge_map[tri_face] = new_internal_edge_index;
                }
                else
                {
                    new_internal_edge_index = it->second;
                }
            }
            else
            {
                const auto it = tri_updates.find(orig_tri_that_was_split);
                if (it == tri_updates.end())
                {
                    throw std::runtime_error(
                      "Inconsistent tetrahedron-triangle update map while splitting edge.");
                }
                const Index new_split_tri_index = it->second.new_triangle_id;
                new_internal_edge_index         = it->second.new_edge_id;
                const bool is_tri_edge_aligned  = it->second.is_edge_aligned;

                if (edge_aligned == is_tri_edge_aligned)
                {
                    if (tri_to_tet_map[orig_tri_that_was_split][0] == orig_tet_index)
                    {
                        tri_to_tet_map[new_split_tri_index][0] = new_tet_index;
                    }
                    else if (tri_to_tet_map[orig_tri_that_was_split][1] == orig_tet_index)
                    {
                        tri_to_tet_map[new_split_tri_index][1] = new_tet_index;
                    }
                    tet_to_tri_map[new_tet_index][split_face_position]  = new_split_tri_index;
                    tet_to_tri_map[orig_tet_index][split_face_position] = orig_tri_that_was_split;
                }
                else
                {
                    if (tri_to_tet_map[orig_tri_that_was_split][0] == orig_tet_index)
                    {
                        tri_to_tet_map[orig_tri_that_was_split][0] = new_tet_index;
                    }
                    else if (tri_to_tet_map[orig_tri_that_was_split][1] == orig_tet_index)
                    {
                        tri_to_tet_map[orig_tri_that_was_split][1] = new_tet_index;
                    }
                    tet_to_tri_map[new_tet_index][split_face_position]  = orig_tri_that_was_split;
                    tet_to_tri_map[orig_tet_index][split_face_position] = new_split_tri_index;
                }
            }

            const ConnIndex edge_to_replace_in_recycled_tet =
              recycled_tet_internal_edge_replace_table[split_edge_pos_in_tet][i];
            tet_to_edges_map[orig_tet_index][edge_to_replace_in_recycled_tet] = new_internal_edge_index;
            edge_to_tets_map[new_internal_edge_index].push_back(
              { orig_tet_index, edge_to_replace_in_recycled_tet });

            const ConnIndex edge_to_replace_in_new_tet =
              new_tet_internal_edge_replace_table[split_edge_pos_in_tet][i];
            tet_to_edges_map[new_tet_index][edge_to_replace_in_new_tet] = new_internal_edge_index;
            edge_to_tets_map[new_internal_edge_index].push_back(
              { new_tet_index, edge_to_replace_in_new_tet });
        }
    }

    return r;
}

void apply_tri_updates(const Index& new_vertex_id,
                       const std::map<Index, TriangleNewEdgeInfo>& tri_updates,
                       UpdateState& state)
{
    state.current_to_orig_edge_map.resize(
      state.current_to_orig_edge_map.size() + tri_updates.size(),
      NOT_SET_Index);

    for (const auto& update : tri_updates)
    {
        state.result.updated_edges.push_back(
          { UpdatedEdgeInfo::NEW_EDGE, NOT_SET_Index, update.second.new_edge_id, new_vertex_id });
    }

    state.current_to_orig_tri_map.resize(
      state.current_to_orig_tri_map.size() + tri_updates.size(),
      NOT_SET_Index);

    for (const auto& update : tri_updates)
    {
        if (state.current_to_orig_tri_map[update.first] != NOT_SET_Index)
        {
            const Index orig_tri_id  = state.current_to_orig_tri_map[update.first];
            const Index new_tri_id   = update.second.new_triangle_id;
            const Index update_index = state.orig_tri_to_update_index[orig_tri_id];
            state.result.updated_triangles[update_index].new_tri_ids.push_back(new_tri_id);
            state.current_to_orig_tri_map[new_tri_id] = orig_tri_id;
        }
        else
        {
            const Index orig_tri_id = update.first;
            const Index new_tri_id  = update.second.new_triangle_id;
            if (state.orig_tri_to_update_index[orig_tri_id] == NOT_SET_Index)
            {
                const auto update_index                     = static_cast<Index>(state.result.updated_triangles.size());
                state.orig_tri_to_update_index[orig_tri_id] = update_index;
                state.current_to_orig_tri_map[new_tri_id]   = orig_tri_id;
                state.result.updated_triangles.push_back({ orig_tri_id, { new_tri_id } });
            }
            else
            {
                const Index update_index                  = state.orig_tri_to_update_index[orig_tri_id];
                state.current_to_orig_tri_map[new_tri_id] = orig_tri_id;
                state.result.updated_triangles[update_index].new_tri_ids.push_back(new_tri_id);
            }
        }
    }
}

void apply_tet_updates(const TetrahedraUpdateResult& tet_updates, UpdateState& state)
{
    state.current_to_orig_edge_map.resize(
      state.current_to_orig_edge_map.size() + tet_updates.new_internal_edges.size(),
      NOT_SET_Index);

    for (const auto& update : tet_updates.new_internal_edges)
    {
        state.result.updated_edges.push_back(
          { UpdatedEdgeInfo::NEW_EDGE, NOT_SET_Index, update.new_edge_id, update.new_node_id });
    }

    state.current_to_orig_tet_map.resize(
      state.current_to_orig_tet_map.size() + tet_updates.original_to_new_tet_map.size(),
      NOT_SET_Index);

    for (const auto& update : tet_updates.original_to_new_tet_map)
    {
        if (state.current_to_orig_tet_map[update.first] != NOT_SET_Index)
        {
            const Index orig_tet_id  = state.current_to_orig_tet_map[update.first];
            const Index new_tet_id   = update.second;
            const Index update_index = state.orig_tet_to_update_index[orig_tet_id];
            state.result.updated_tetrahedra[update_index].new_tet_ids.push_back(new_tet_id);
            state.current_to_orig_tet_map[new_tet_id] = orig_tet_id;
        }
        else
        {
            const Index orig_tet_id = update.first;
            const Index new_tet_id  = update.second;
            if (state.orig_tet_to_update_index[orig_tet_id] == NOT_SET_Index)
            {
                const auto update_index                     = static_cast<Index>(state.result.updated_tetrahedra.size());
                state.orig_tet_to_update_index[orig_tet_id] = update_index;
                state.current_to_orig_tet_map[new_tet_id]   = orig_tet_id;
                state.result.updated_tetrahedra.push_back({ orig_tet_id, { new_tet_id } });
            }
            else
            {
                const Index update_index                  = state.orig_tet_to_update_index[orig_tet_id];
                state.current_to_orig_tet_map[new_tet_id] = orig_tet_id;
                state.result.updated_tetrahedra[update_index].new_tet_ids.push_back(new_tet_id);
            }
        }
    }
}

void split_single_edge(const Index& edge_id,
                       const Index& existing_vertex_id_for_split,
                       Tetrahedralization& t,
                       UpdateState& state)
{
    if (edge_id < state.orig_edge_to_update_index.size())
    {
        assert(state.orig_edge_to_update_index[edge_id] == NOT_SET_Index);
    }

    auto& vertices = t.vertices;
    auto& edges    = t.edges;

    const Index v_1   = edges[edge_id][0];
    const Index v_2   = edges[edge_id][1];
    const Point3& p_1 = vertices[v_1];
    const Point3& p_2 = vertices[v_2];

    Index new_v = existing_vertex_id_for_split;
    if (new_v == NOT_SET_Index)
    {
        new_v = static_cast<Index>(vertices.size());
        vertices.push_back(detail::midpoint(p_1, p_2));
    }

    const Index new_n_n_2_edge_id = static_cast<Index>(edges.size());
    edges.push_back({ new_v, v_2 });
    edges[edge_id][1] = new_v;

    const Index orig_edge_id =
      (edge_id < state.current_to_orig_edge_map.size() &&
       state.current_to_orig_edge_map[edge_id] != NOT_SET_Index)
        ? state.current_to_orig_edge_map[edge_id]
        : edge_id;

    if (edge_id < state.orig_edge_to_update_index.size())
    {
        state.orig_edge_to_update_index[edge_id] =
          static_cast<Index>(state.result.updated_edges.size());
    }
    state.current_to_orig_edge_map.push_back(orig_edge_id);
    state.result.updated_edges.push_back(
      { UpdatedEdgeInfo::SUBDIVIDED_EDGE, orig_edge_id, new_n_n_2_edge_id, new_v });

    const auto tri_updates = split_triangles_on_edge(edge_id, new_n_n_2_edge_id, t);
    if (GlobalOptions::get().check_common_errors)
    {
        check_tri_edge_relations(t);
    }

    const auto tet_updates = update_tetrahedra_on_split_edge(edge_id, new_n_n_2_edge_id, tri_updates, t);
    if (GlobalOptions::get().check_common_errors)
    {
        check_tet_edge_relations(t);
        check_tri_edge_relations(t);
        check_tet_tri_relations(t);
    }

    apply_tri_updates(new_v, tri_updates, state);
    apply_tet_updates(tet_updates, state);
}

void split_single_edge(const Index& edge_id,
                       const Index& existing_vertex_id_for_split,
                       Triangulation& t,
                       UpdateState& state)
{
    if (edge_id < state.orig_edge_to_update_index.size())
    {
        assert(state.orig_edge_to_update_index[edge_id] == NOT_SET_Index);
    }

    auto& vertices = t.vertices;
    auto& edges    = t.edges;

    const Index v_1   = edges[edge_id][0];
    const Index v_2   = edges[edge_id][1];
    const Point3& p_1 = vertices[v_1];
    const Point3& p_2 = vertices[v_2];

    Index new_v = existing_vertex_id_for_split;
    if (new_v == NOT_SET_Index)
    {
        new_v = static_cast<Index>(vertices.size());
        vertices.push_back(detail::midpoint(p_1, p_2));
    }

    const Index new_n_n_2_edge_id = static_cast<Index>(edges.size());
    edges.push_back({ new_v, v_2 });
    edges[edge_id][1] = new_v;

    const Index orig_edge_id =
      (edge_id < state.current_to_orig_edge_map.size() &&
       state.current_to_orig_edge_map[edge_id] != NOT_SET_Index)
        ? state.current_to_orig_edge_map[edge_id]
        : edge_id;

    if (edge_id < state.orig_edge_to_update_index.size())
    {
        state.orig_edge_to_update_index[edge_id] =
          static_cast<Index>(state.result.updated_edges.size());
    }
    state.current_to_orig_edge_map.push_back(orig_edge_id);
    state.result.updated_edges.push_back(
      { UpdatedEdgeInfo::SUBDIVIDED_EDGE, orig_edge_id, new_n_n_2_edge_id, new_v });

    const auto tri_updates = split_triangles_on_edge(edge_id, new_n_n_2_edge_id, t);
    if (GlobalOptions::get().check_common_errors)
    {
        check_tri_edge_relations(t);
    }

    apply_tri_updates(new_v, tri_updates, state);
}

void update_triangulation_markers(const UpdateResult& update_result, Triangulation& t)
{
    t.tri_to_marker_map.resize(t.triangles.size(), NOT_SET_Marker);
    for (const auto& tri_update : update_result.updated_triangles)
    {
        const Marker orig_marker = t.tri_to_marker_map[tri_update.orig_tri_id];
        for (const Index new_tri_id : tri_update.new_tri_ids)
        {
            t.tri_to_marker_map[new_tri_id] = orig_marker;
        }
    }
}

} // namespace

UpdateResult split_edges(const std::vector<Index>& edges_to_split, Tetrahedralization& t)
{
    UpdateState state;
    state.orig_edge_to_update_index.resize(t.edges.size(), NOT_SET_Index);
    state.orig_tri_to_update_index.resize(t.triangles.size(), NOT_SET_Index);
    state.orig_tet_to_update_index.resize(t.tetrahedra.size(), NOT_SET_Index);

    state.current_to_orig_edge_map.resize(t.edges.size(), NOT_SET_Index);
    state.current_to_orig_tri_map.resize(t.triangles.size(), NOT_SET_Index);
    state.current_to_orig_tet_map.resize(t.tetrahedra.size(), NOT_SET_Index);

    for (const Index edge_id : edges_to_split)
    {
        split_single_edge(edge_id, NOT_SET_Index, t, state);
    }
    return state.result;
}

UpdateResult split_edges(const std::vector<Index>& edges_to_split,
                         const std::vector<Index>& existing_vertices_for_split_edges,
                         Tetrahedralization& t)
{
    UpdateState state;
    state.orig_edge_to_update_index.resize(t.edges.size(), NOT_SET_Index);
    state.orig_tri_to_update_index.resize(t.triangles.size(), NOT_SET_Index);
    state.orig_tet_to_update_index.resize(t.tetrahedra.size(), NOT_SET_Index);

    state.current_to_orig_edge_map.resize(t.edges.size(), NOT_SET_Index);
    state.current_to_orig_tri_map.resize(t.triangles.size(), NOT_SET_Index);
    state.current_to_orig_tet_map.resize(t.tetrahedra.size(), NOT_SET_Index);

    for (std::size_t i = 0; i < edges_to_split.size(); ++i)
    {
        split_single_edge(edges_to_split[i], existing_vertices_for_split_edges[i], t, state);
    }
    return state.result;
}

UpdateResult split_edges(const std::vector<Index>& edges_to_split, Triangulation& t)
{
    UpdateState state;
    state.orig_edge_to_update_index.resize(t.edges.size(), NOT_SET_Index);
    state.orig_tri_to_update_index.resize(t.triangles.size(), NOT_SET_Index);

    state.current_to_orig_edge_map.resize(t.edges.size(), NOT_SET_Index);
    state.current_to_orig_tri_map.resize(t.triangles.size(), NOT_SET_Index);

    for (const Index edge_id : edges_to_split)
    {
        split_single_edge(edge_id, NOT_SET_Index, t, state);
    }

    update_triangulation_markers(state.result, t);
    return state.result;
}

UpdateResult split_edges(const std::vector<Index>& edges_to_split,
                         const std::vector<Index>& existing_vertices_for_split_edges,
                         Triangulation& t)
{
    UpdateState state;
    state.orig_edge_to_update_index.resize(t.edges.size(), NOT_SET_Index);
    state.orig_tri_to_update_index.resize(t.triangles.size(), NOT_SET_Index);

    state.current_to_orig_edge_map.resize(t.edges.size(), NOT_SET_Index);
    state.current_to_orig_tri_map.resize(t.triangles.size(), NOT_SET_Index);

    for (std::size_t i = 0; i < edges_to_split.size(); ++i)
    {
        split_single_edge(edges_to_split[i], existing_vertices_for_split_edges[i], t, state);
    }

    update_triangulation_markers(state.result, t);
    return state.result;
}

UpdateResult split_constrained_edges(Triangulation& t, SeamInfo& seams)
{
    std::vector<bool> is_node_on_seam(t.vertices.size(), false);
    for (const Index edge_id : seams.edges_on_seams)
    {
        is_node_on_seam[t.edges[edge_id][0]] = true;
        is_node_on_seam[t.edges[edge_id][1]] = true;
    }

    std::vector<Index> edges_to_split;
    for (std::size_t edge_id = 0; edge_id < t.edges.size(); ++edge_id)
    {
        const Index v_1 = t.edges[edge_id][0];
        const Index v_2 = t.edges[edge_id][1];
        if (is_node_on_seam[v_1] && is_node_on_seam[v_2] &&
            seams.edge_id_to_seam_id_map[edge_id] == NOT_SET_SeamID)
        {
            edges_to_split.push_back(static_cast<Index>(edge_id));
        }
    }

    if (GlobalOptions::get().verbose)
    {
        std::cout << GlobalOptions::get().curr_status_padding
                  << "* split_constrained_edges - Splitting " << edges_to_split.size() << " edges."
                  << std::endl;
    }

    const auto update_result = split_edges(edges_to_split, t);

    seams.edge_id_to_seam_id_map.resize(t.edges.size(), NOT_SET_SeamID);
    for (const auto& update : update_result.updated_edges)
    {
        if (update.update_type == UpdatedEdgeInfo::SUBDIVIDED_EDGE &&
            seams.edge_id_to_seam_id_map[update.orig_edge_id] != NOT_SET_SeamID)
        {
            seams.edges_on_seams.push_back(update.new_edge_id);
            seams.edge_id_to_seam_id_map[update.new_edge_id] =
              seams.edge_id_to_seam_id_map[update.orig_edge_id];
        }
    }

    return update_result;
}

void split_large_edges(const double max_edge_length, Triangulation& t, SeamInfo& seams)
{
    if (GlobalOptions::get().verbose)
    {
        std::cout << GlobalOptions::get().curr_status_padding
                  << "* split_large_edges - Splitting edges larger than " << max_edge_length
                  << std::endl;
    }
    while (true)
    {
        std::vector<Index> edges_to_split;
        for (std::size_t edge_id = 0; edge_id < t.edges.size(); ++edge_id)
        {
            const auto& edge    = t.edges[edge_id];
            const double length = detail::distance(t.vertices[edge[0]], t.vertices[edge[1]]);
            if (length > max_edge_length)
            {
                edges_to_split.push_back(static_cast<Index>(edge_id));
            }
        }
        if (edges_to_split.empty())
        {
            break;
        }

        const auto update_result = split_edges(edges_to_split, t);

        seams.edge_id_to_seam_id_map.resize(t.edges.size(), NOT_SET_SeamID);
        for (const auto& update : update_result.updated_edges)
        {
            if (update.update_type == UpdatedEdgeInfo::SUBDIVIDED_EDGE &&
                seams.edge_id_to_seam_id_map[update.orig_edge_id] != NOT_SET_SeamID)
            {
                seams.edges_on_seams.push_back(update.new_edge_id);
                seams.edge_id_to_seam_id_map[update.new_edge_id] =
                  seams.edge_id_to_seam_id_map[update.orig_edge_id];
            }
        }
    }
}

} // namespace mesh_arrange
