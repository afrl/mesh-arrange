#include "mesh_arrange/Optimize.h"

#include "detail/Containers.h"
#include "detail/DebugVTK.h"
#include "detail/ElementTables.h"
#include "detail/IGLHelpers.h"
#include "detail/MathHelpers.h"
#include "mesh_arrange/GlobalOptions.h"
#include "mesh_arrange/Tetrahedralization.h"
#include "mesh_arrange/VertexMap.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>
#include <stdexcept>

namespace mesh_arrange
{

namespace
{

void write_seams_vtk(const Triangulation& t, const SeamInfo& seam_info, const std::string& prefix)
{
    std::vector<std::array<Index, 2>> lines;
    std::vector<uint32_t> seam_ids;
    for (const Index edge_id : seam_info.edges_on_seams)
    {
        lines.push_back(t.edges[edge_id]);
        seam_ids.push_back(seam_info.edge_id_to_seam_id_map[edge_id]);
    }
    detail::write_vtk_lines(prefix + "seams.vtk", t.vertices, lines, { { "seam_id", seam_ids } });
}

struct AdjacentNode
{
    Index node_id;
    Index edge_id;
};

int find_external_vertex_position(const Triangulation& t, const Index& tri_id, const Index& edge_id)
{
    const auto& edge = t.edges[edge_id];
    for (int i = 0; i < 3; ++i)
    {
        const Index v_id = t.triangles[tri_id][i];
        if (v_id != edge[0] && v_id != edge[1])
        {
            return i;
        }
    }
    throw std::runtime_error("Bad triangle and edge combination given. The triangle may be collapsed.");
}

int get_aligned_case(const Triangulation& t,
                     const bool is_t1,
                     const Index& tri_id,
                     const Index& edge_id,
                     const int edge_position)
{
    const Index first_v = t.triangles[tri_id][detail::kTriEdges[edge_position][0]];
    if (t.edges[edge_id][0] == first_v)
    {
        return is_t1 ? 0 : 1;
    }
    if (t.edges[edge_id][1] == first_v)
    {
        return is_t1 ? 1 : 0;
    }
    throw std::runtime_error("Bad triangle and edge combination given.");
}

bool can_safely_flip_edge(const Triangulation& t,
                          const Index edge_to_flip_id,
                          const std::vector<std::vector<Index>>& node_to_edge_map)
{
    const auto& connected_tris = t.edge_to_tris_map[edge_to_flip_id];
    if (connected_tris.size() != 2)
    {
        return false;
    }

    const Index t1_id  = connected_tris[0].first;
    const Index t2_id  = connected_tris[1].first;
    const int t1_case  = connected_tris[0].second;
    const int t2_case  = connected_tris[1].second;
    const int t1_align = get_aligned_case(t, true, t1_id, edge_to_flip_id, t1_case);
    const int t2_align = get_aligned_case(t, false, t2_id, edge_to_flip_id, t2_case);

    static constexpr int ext_v_pos[3] = { 2, 0, 1 };
    const Index v0_id                 = t.edges[edge_to_flip_id][0];
    const Index v2_id                 = t.edges[edge_to_flip_id][1];
    const Index v1_id                 = t.triangles[t2_id][ext_v_pos[t2_case]];
    const Index v3_id                 = t.triangles[t1_id][ext_v_pos[t1_case]];

    if (v1_id == v3_id || v1_id == v0_id || v1_id == v2_id || v3_id == v0_id || v3_id == v2_id)
    {
        return false;
    }

    for (const Index existing_e_id : node_to_edge_map[v1_id])
    {
        const auto& e = t.edges[existing_e_id];
        if ((e[0] == v1_id && e[1] == v3_id) || (e[0] == v3_id && e[1] == v1_id))
        {
            return false;
        }
    }

    auto tri_normal = [&](const std::array<Index, 3>& tri) -> Eigen::Vector3d
    {
        const Eigen::Vector3d p0 = detail::point_as_vector(t.vertices[tri[0]]);
        const Eigen::Vector3d p1 = detail::point_as_vector(t.vertices[tri[1]]);
        const Eigen::Vector3d p2 = detail::point_as_vector(t.vertices[tri[2]]);
        return (p1 - p0).cross(p2 - p0);
    };

    const Eigen::Vector3d n1_before = tri_normal(t.triangles[t1_id]);
    const Eigen::Vector3d n2_before = tri_normal(t.triangles[t2_id]);
    const double area1_b            = n1_before.norm();
    const double area2_b            = n2_before.norm();
    if (area1_b <= 1e-15 || area2_b <= 1e-15)
    {
        return false;
    }

    const Eigen::Vector3d u1_before = n1_before / area1_b;
    const Eigen::Vector3d u2_before = n2_before / area2_b;
    if (u1_before.dot(u2_before) < 0.5)
    {
        return false;
    }

    const std::array<Index, 3> new_t1 =
      (t1_align == 0) ? std::array<Index, 3>{ v1_id, v3_id, v0_id }
                      : std::array<Index, 3>{ v3_id, v1_id, v0_id };
    const std::array<Index, 3> new_t2 =
      (t2_align == 0) ? std::array<Index, 3>{ v3_id, v1_id, v2_id }
                      : std::array<Index, 3>{ v1_id, v3_id, v2_id };

    const Eigen::Vector3d n1_after = tri_normal(new_t1);
    const Eigen::Vector3d n2_after = tri_normal(new_t2);
    const double area1_a           = n1_after.norm();
    const double area2_a           = n2_after.norm();
    if (area1_a <= 1e-15 || area2_a <= 1e-15)
    {
        return false;
    }

    const Eigen::Vector3d u1_after   = n1_after / area1_a;
    const Eigen::Vector3d u2_after   = n2_after / area2_a;
    const Eigen::Vector3d avg_before = (u1_before + u2_before).normalized();

    if (u1_after.dot(avg_before) <= 0.5 || u2_after.dot(avg_before) <= 0.5 ||
        u1_after.dot(u2_after) < 0.5)
    {
        return false;
    }

    return true;
}

void flip_edge(Triangulation& t, const Index& edge_to_flip_id)
{
    auto& connected_tris = t.edge_to_tris_map[edge_to_flip_id];
    assert(connected_tris.size() == 2);
    assert(t.tri_to_marker_map[connected_tris[0].first] ==
           t.tri_to_marker_map[connected_tris[1].first]);

    static const int ext_v_pos[3]       = { 2, 0, 1 };
    static const int ext_e_pos[2][3]    = { { 1, 2, 0 }, { 2, 0, 1 } };
    static const int reused_e_pos[2][3] = { { 2, 0, 1 }, { 1, 2, 0 } };

    const Index t1_id = connected_tris[0].first;
    const Index t2_id = connected_tris[1].first;

    const int t1_case    = connected_tris[0].second;
    const int t1_align_i = get_aligned_case(t, true, t1_id, edge_to_flip_id, t1_case);
    const int t2_case    = connected_tris[1].second;
    const int t2_align_i = get_aligned_case(t, false, t2_id, edge_to_flip_id, t2_case);

    const Index v0_id = t.edges[edge_to_flip_id][0];
    const Index v1_id = t.triangles[t2_id][ext_v_pos[t2_case]];
    const Index v2_id = t.edges[edge_to_flip_id][1];
    const Index v3_id = t.triangles[t1_id][ext_v_pos[t1_case]];

    const Index e0_id = t.tri_to_edges_map[t2_id][ext_e_pos[t2_align_i][t2_case]];
    const Index e1_id = t.tri_to_edges_map[t2_id][reused_e_pos[t2_align_i][t2_case]];
    const Index e2_id = t.tri_to_edges_map[t1_id][ext_e_pos[t1_align_i][t1_case]];
    const Index e3_id = t.tri_to_edges_map[t1_id][reused_e_pos[t1_align_i][t1_case]];

    auto replace_edge = [&t](const Index edge_id,
                             const Index tri_to_find,
                             const Index tri_id,
                             const ConnIndex edge_position)
    {
        for (auto& tri : t.edge_to_tris_map[edge_id])
        {
            if (tri.first == tri_to_find)
            {
                tri.first  = tri_id;
                tri.second = edge_position;
                return;
            }
        }
        throw std::runtime_error("Triangle " + std::to_string(tri_to_find) +
                                 " did not appear in map for edge " + std::to_string(edge_id));
    };

    t.edges[edge_to_flip_id] = { v1_id, v3_id };

    if (t1_align_i == 0)
    {
        t.triangles[t1_id]        = { v1_id, v3_id, v0_id };
        t.tri_to_edges_map[t1_id] = { edge_to_flip_id, e3_id, e0_id };
        replace_edge(edge_to_flip_id, t1_id, t1_id, 0);
        replace_edge(e3_id, t1_id, t1_id, 1);
        replace_edge(e0_id, t2_id, t1_id, 2);
    }
    else
    {
        t.triangles[t1_id]        = { v3_id, v1_id, v0_id };
        t.tri_to_edges_map[t1_id] = { edge_to_flip_id, e0_id, e3_id };
        replace_edge(edge_to_flip_id, t1_id, t1_id, 0);
        replace_edge(e0_id, t2_id, t1_id, 1);
        replace_edge(e3_id, t1_id, t1_id, 2);
    }

    if (t2_align_i == 0)
    {
        t.triangles[t2_id]        = { v3_id, v1_id, v2_id };
        t.tri_to_edges_map[t2_id] = { edge_to_flip_id, e1_id, e2_id };
        replace_edge(edge_to_flip_id, t2_id, t2_id, 0);
        replace_edge(e1_id, t2_id, t2_id, 1);
        replace_edge(e2_id, t1_id, t2_id, 2);
    }
    else
    {
        t.triangles[t2_id]        = { v1_id, v3_id, v2_id };
        t.tri_to_edges_map[t2_id] = { edge_to_flip_id, e2_id, e1_id };
        replace_edge(edge_to_flip_id, t2_id, t2_id, 0);
        replace_edge(e2_id, t1_id, t2_id, 1);
        replace_edge(e1_id, t2_id, t2_id, 2);
    }
}

} // namespace

void optimize_tangential_smoothing(Triangulation& t,
                                   const SeamInfo& seams,
                                   const int n_iterations,
                                   const bool allow_smoothing_on_seams)
{
    if (t.vertices.empty() || t.triangles.empty() || t.tri_to_marker_map.empty() ||
        n_iterations <= 0)
    {
        return;
    }

    GlobalOptions::get().log_info(
      GlobalOptions::get().curr_status_padding +
      "* optimize_tangential_smoothing - Starting iterative optimization");

    if (GlobalOptions::get().write_debug_meshes)
    {
        detail::write_vtk_triangulation(GlobalOptions::get().output_dir + "/iteration_0.vtk", t);
        write_seams_vtk(t, seams, GlobalOptions::get().output_dir + "/iteration_0_");
    }

    const Marker n_facets =
      *std::max_element(t.tri_to_marker_map.begin(), t.tri_to_marker_map.end()) + 1;
    t.n_markers = std::max(t.n_markers, n_facets);

    const auto vert_to_tri_map    = build_vert_to_tri_map(t);
    const auto original_positions = t.vertices;

    std::vector<std::vector<Index>> facet_id_to_vert_id_map(n_facets);
    for (Index v_id = 0; v_id < static_cast<Index>(t.vertices.size()); ++v_id)
    {
        for (const Index tri_id : vert_to_tri_map.get_indices(v_id))
        {
            facet_id_to_vert_id_map[t.tri_to_marker_map[tri_id]].push_back(v_id);
        }
    }
    for (auto& verts : facet_id_to_vert_id_map)
    {
        std::sort(verts.begin(), verts.end());
        verts.erase(std::unique(verts.begin(), verts.end()), verts.end());
    }

    std::vector<bool> is_vert_on_seam(t.vertices.size(), false);
    std::vector<uint32_t> n_seam_elements_connected_to_vert(t.vertices.size(), 0);
    std::vector<std::vector<Index>> seam_id_to_vert_id_map(seams.n_seams);
    std::vector<std::vector<AdjacentNode>> vert_to_neighbors_on_seam(t.vertices.size());

    for (const Index feature_edge_id : seams.edges_on_seams)
    {
        const Marker seam_id   = seams.edge_id_to_seam_id_map[feature_edge_id];
        const Index v1_id      = t.edges[feature_edge_id][0];
        const Index v2_id      = t.edges[feature_edge_id][1];
        is_vert_on_seam[v1_id] = true;
        is_vert_on_seam[v2_id] = true;
        n_seam_elements_connected_to_vert[v1_id]++;
        n_seam_elements_connected_to_vert[v2_id]++;
        seam_id_to_vert_id_map[seam_id].push_back(v1_id);
        seam_id_to_vert_id_map[seam_id].push_back(v2_id);
        vert_to_neighbors_on_seam[v1_id].push_back({ v2_id, feature_edge_id });
        vert_to_neighbors_on_seam[v2_id].push_back({ v1_id, feature_edge_id });
    }
    for (auto& verts : seam_id_to_vert_id_map)
    {
        std::sort(verts.begin(), verts.end());
        verts.erase(std::unique(verts.begin(), verts.end()), verts.end());
    }

    std::vector<detail::VectorSet<Marker>> vert_id_to_input_facet_ids_map(t.vertices.size());
    for (std::size_t tri_id = 0; tri_id < t.triangles.size(); ++tri_id)
    {
        for (const Index v_id : t.triangles[tri_id])
        {
            vert_id_to_input_facet_ids_map[v_id].insert(t.tri_to_marker_map[tri_id]);
        }
    }

    std::vector<bool> is_vert_on_corner(t.vertices.size(), false);
    for (std::size_t n_id = 0; n_id < t.vertices.size(); ++n_id)
    {
        if (is_vert_on_seam[n_id])
        {
            if (n_seam_elements_connected_to_vert[n_id] != 2)
            {
                is_vert_on_corner[n_id] = true;
            }
            else
            {
                const Marker s0 =
                  seams.edge_id_to_seam_id_map[vert_to_neighbors_on_seam[n_id][0].edge_id];
                const Marker s1 =
                  seams.edge_id_to_seam_id_map[vert_to_neighbors_on_seam[n_id][1].edge_id];
                if (s0 != s1)
                {
                    is_vert_on_corner[n_id] = true;
                }
            }
        }
    }

    const auto igl_facets = detail::decompose_into_igl_facets(t);

    for (int iteration = 0; iteration < n_iterations; ++iteration)
    {
        std::vector<Point3> updated_coords(t.vertices.size(), { 0.0, 0.0, 0.0 });

        if (allow_smoothing_on_seams)
        {
            for (Marker seam_id = 0; seam_id < seams.n_seams; ++seam_id)
            {
                for (const Index v_id : seam_id_to_vert_id_map[seam_id])
                {
                    if (!is_vert_on_corner[v_id])
                    {
                        Eigen::Vector3d u       = Eigen::Vector3d::Zero();
                        Eigen::Vector3d tv      = Eigen::Vector3d::Zero();
                        double edge_length_sums = 0.0;
                        int n_connected_edges   = 0;
                        const Eigen::Vector3d p1 = detail::point_as_vector(t.vertices[v_id]);

                        for (const auto& neighbor : vert_to_neighbors_on_seam[v_id])
                        {
                            const Eigen::Vector3d p2 =
                              detail::point_as_vector(t.vertices[neighbor.node_id]);
                            const Eigen::Vector3d mid = 0.5 * (p1 + p2);
                            const double dist         = (p2 - p1).norm();
                            edge_length_sums += dist;
                            u += dist * mid;

                            if (n_connected_edges == 0)
                            {
                                if (dist > 0.0)
                                {
                                    tv += (p2 - p1) / dist;
                                }
                            }
                            else
                            {
                                if (dist > 0.0)
                                {
                                    tv -= (p2 - p1) / dist;
                                }
                            }
                            n_connected_edges++;
                        }

                        assert(n_connected_edges == 2);
                        const double tv_norm = tv.norm();
                        if (edge_length_sums > 0.0 && tv_norm > 0.0)
                        {
                            u /= edge_length_sums;
                            u -= p1;
                            tv /= tv_norm;
                            u                    = tv * u.dot(tv);
                            updated_coords[v_id] = detail::to_point(u);
                        }
                    }
                }
            }
        }

        for (Marker facet_id = 0; facet_id < n_facets; ++facet_id)
        {
            for (const Index v_id : facet_id_to_vert_id_map[facet_id])
            {
                if (!is_vert_on_corner[v_id] && !is_vert_on_seam[v_id])
                {
                    Eigen::Vector3d u            = Eigen::Vector3d::Zero();
                    double ww                    = 0.0;
                    const Eigen::Vector3d p_orig = detail::point_as_vector(t.vertices[v_id]);

                    for (const Index tri_id : vert_to_tri_map.get_indices(v_id))
                    {
                        if (t.tri_to_marker_map[tri_id] == facet_id)
                        {
                            const Eigen::Vector3d p1 =
                              detail::point_as_vector(t.vertices[t.triangles[tri_id][0]]);
                            const Eigen::Vector3d p2 =
                              detail::point_as_vector(t.vertices[t.triangles[tri_id][1]]);
                            const Eigen::Vector3d p3 =
                              detail::point_as_vector(t.vertices[t.triangles[tri_id][2]]);
                            const Eigen::Vector3d avg = (p1 + p2 + p3) / 3.0;
                            const double area         = (p2 - p1).cross(p3 - p1).norm();

                            ww += area;
                            u += area * avg;
                        }
                    }

                    if (ww > 0.0)
                    {
                        u /= ww;
                        u -= p_orig;
                        updated_coords[v_id] = detail::to_point(u);
                    }
                }
            }
        }

        for (std::size_t v_id = 0; v_id < t.vertices.size(); ++v_id)
        {
            if (!is_vert_on_corner[v_id])
            {
                t.vertices[v_id] = detail::add_points(t.vertices[v_id], updated_coords[v_id]);
            }
        }

        for (std::size_t v_id = 0; v_id < t.vertices.size(); ++v_id)
        {
            if (!is_vert_on_corner[v_id] && !is_vert_on_seam[v_id])
            {
                Eigen::Vector3d new_point = detail::point_as_vector(t.vertices[v_id]);
                for (const Marker facet_id : vert_id_to_input_facet_ids_map[v_id])
                {
                    if (facet_id < igl_facets.size() && igl_facets[facet_id].triangles.rows() > 0)
                    {
                        new_point = detail::igl_closest_point(igl_facets[facet_id], new_point);
                    }
                }
                t.vertices[v_id] = detail::to_point(new_point);
            }
        }

        if (allow_smoothing_on_seams)
        {
            const std::size_t max_seam_steps = seams.edges_on_seams.size() + 1;
            for (Marker seam_id = 0; seam_id < seams.n_seams; ++seam_id)
            {
                for (const Index v_id : seam_id_to_vert_id_map[seam_id])
                {
                    if (!is_vert_on_corner[v_id])
                    {
                        const Point3 p        = t.vertices[v_id];
                        const auto& neighbor1 = vert_to_neighbors_on_seam[v_id][0];
                        const auto result1    = detail::closest_point_on_line_segment(
                          p,
                          original_positions[t.edges[neighbor1.edge_id][0]],
                          original_positions[t.edges[neighbor1.edge_id][1]]);

                        const auto& neighbor2 = vert_to_neighbors_on_seam[v_id][1];
                        const auto result2    = detail::closest_point_on_line_segment(
                          p,
                          original_positions[t.edges[neighbor2.edge_id][0]],
                          original_positions[t.edges[neighbor2.edge_id][1]]);

                        Index curr_n_id    = v_id;
                        auto curr_neighbor = neighbor1;
                        auto curr_result   = result1;
                        if (result2.distance < result1.distance)
                        {
                            curr_neighbor = neighbor2;
                            curr_result   = result2;
                        }

                        std::size_t steps = 0;
                        while (steps++ < max_seam_steps)
                        {
                            if (is_vert_on_corner[curr_neighbor.node_id])
                            {
                                break;
                            }
                            const auto& next_neighbors =
                              vert_to_neighbors_on_seam[curr_neighbor.node_id];
                            if (next_neighbors.size() != 2)
                            {
                                break;
                            }

                            auto next_neighbor = next_neighbors[0];
                            if (next_neighbors[0].node_id == curr_n_id)
                            {
                                next_neighbor = next_neighbors[1];
                            }
                            const auto next_result = detail::closest_point_on_line_segment(
                              p,
                              original_positions[t.edges[next_neighbor.edge_id][0]],
                              original_positions[t.edges[next_neighbor.edge_id][1]]);

                            if (next_result.distance >= curr_result.distance)
                            {
                                break;
                            }

                            curr_n_id     = curr_neighbor.node_id;
                            curr_neighbor = next_neighbor;
                            curr_result   = next_result;
                        }

                        t.vertices[v_id] = curr_result.closest_point;
                    }
                }
            }
        }
    }

    if (GlobalOptions::get().write_debug_meshes)
    {
        detail::write_vtk_triangulation(
          GlobalOptions::get().output_dir + "/post_tangential_smoothing.vtk",
          t);
    }
}

void optimize_edge_flips(Triangulation& t, const int n_passes)
{
    check_tri_edge_relations(t);

    if (t.vertices.empty() || t.edges.empty() || n_passes <= 0)
    {
        return;
    }

    std::vector<int> vertex_valences(t.vertices.size(), 0);
    std::vector<std::vector<Index>> node_to_edge_map(t.vertices.size());
    for (std::size_t edge_id = 0; edge_id < t.edges.size(); ++edge_id)
    {
        const auto& edge = t.edges[edge_id];
        vertex_valences[edge[0]]++;
        vertex_valences[edge[1]]++;
        node_to_edge_map[edge[0]].push_back(static_cast<Index>(edge_id));
        node_to_edge_map[edge[1]].push_back(static_cast<Index>(edge_id));
    }

    auto remove_edge_from_node = [&](Index node_id, Index edge_id)
    {
        auto& vec = node_to_edge_map[node_id];
        auto it   = std::find(vec.begin(), vec.end(), edge_id);
        if (it != vec.end())
        {
            vec.erase(it);
        }
    };

    for (int pass = 0; pass < n_passes; ++pass)
    {
        auto sorted_permutation = detail::find_sorted_permutation(vertex_valences);
        std::reverse(sorted_permutation.begin(), sorted_permutation.end());

        int flips_in_pass = 0;
        for (std::size_t curr_i = 0; curr_i < sorted_permutation.size(); ++curr_i)
        {
            const std::size_t v_id = sorted_permutation[curr_i];
            int& v_valence         = vertex_valences[v_id];
            if (v_valence <= 7)
            {
                continue;
            }

            const auto orig_connected_edges = node_to_edge_map[v_id];
            const std::size_t n_conn_edges  = orig_connected_edges.size();
            std::vector<int> initial_max_external_valences(n_conn_edges, 1000);
            std::vector<std::array<Index, 2>> external_vertices(n_conn_edges, { 0, 0 });
            std::vector<bool> is_connected_edge_locked(n_conn_edges, false);

            for (std::size_t i = 0; i < n_conn_edges; ++i)
            {
                const Index edge_id = orig_connected_edges[i];
                if (t.edge_to_tris_map[edge_id].size() != 2 ||
                    t.tri_to_marker_map[t.edge_to_tris_map[edge_id][0].first] !=
                      t.tri_to_marker_map[t.edge_to_tris_map[edge_id][1].first])
                {
                    is_connected_edge_locked[i] = true;
                    continue;
                }

                const Index t1_id    = t.edge_to_tris_map[edge_id][0].first;
                const int i_1        = find_external_vertex_position(t, t1_id, edge_id);
                const Index t2_id    = t.edge_to_tris_map[edge_id][1].first;
                const int i_2        = find_external_vertex_position(t, t2_id, edge_id);
                external_vertices[i] = { t.triangles[t1_id][i_1], t.triangles[t2_id][i_2] };

                initial_max_external_valences[i] = std::max(
                  vertex_valences[external_vertices[i][0]],
                  vertex_valences[external_vertices[i][1]]);
            }

            const auto neighbor_sorted_p =
              detail::find_sorted_permutation(initial_max_external_valences);
            for (std::size_t i = 0; i < n_conn_edges; ++i)
            {
                const std::size_t sorted_i = neighbor_sorted_p[i];
                if (is_connected_edge_locked[sorted_i])
                {
                    continue;
                }
                const Index edge_id = orig_connected_edges[sorted_i];

                if (v_valence > 7 &&
                    v_valence > vertex_valences[external_vertices[sorted_i][0]] + 1 &&
                    v_valence > vertex_valences[external_vertices[sorted_i][1]] + 1 &&
                    can_safely_flip_edge(t, edge_id, node_to_edge_map))
                {
                    const auto orig_edge = t.edges[edge_id];
                    flip_edge(t, edge_id);
                    flips_in_pass++;

                    vertex_valences[orig_edge[0]]--;
                    vertex_valences[orig_edge[1]]--;
                    vertex_valences[t.edges[edge_id][0]]++;
                    vertex_valences[t.edges[edge_id][1]]++;

                    remove_edge_from_node(orig_edge[0], edge_id);
                    remove_edge_from_node(orig_edge[1], edge_id);
                    node_to_edge_map[t.edges[edge_id][0]].push_back(edge_id);
                    node_to_edge_map[t.edges[edge_id][1]].push_back(edge_id);
                }
            }
        }

        if (flips_in_pass == 0)
        {
            break;
        }
    }

    check_tri_edge_relations(t);
}

} // namespace mesh_arrange
