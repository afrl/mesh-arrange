#include "mesh_arrange/QualityTriangulation.h"

#include "detail/IGLHelpers.h"
#include "mesh_arrange/Optimize.h"
#include "mesh_arrange/SplitEdges.h"

#include <cassert>
#include <set>
#include <vector>

namespace mesh_arrange
{

RecoveredTriangulation recover_triangulation(SegmentedTetrahedralization& t)
{
    const std::size_t n_objects = t.object_region_info.object_to_regions_map.size();

    RecoveredTriangulation objects;
    objects.object_to_markers_map.resize(n_objects);

    std::vector<Index> t_to_objects_node_map(t.tetrahedralization.vertices.size(), NOT_SET_Index);
    std::vector<Index> objects_to_t_node_map;
    objects_to_t_node_map.reserve(t.tetrahedralization.vertices.size());

    Index curr_v_id = 0;
    for (std::size_t object_id = 0; object_id < n_objects; ++object_id)
    {
        std::set<Marker> object_markers;
        for (const RegionID region_id : t.object_region_info.object_to_regions_map[object_id])
        {
            const auto& tri_ids = t.region_info.region_to_bounding_triangles_map[region_id];
            for (const Index tri_id : tri_ids)
            {
                object_markers.insert(t.region_boundaries_info.triangle_to_boundary_id_map[tri_id]);
                const auto& conn = t.tetrahedralization.triangles[tri_id];
                for (int i = 0; i < 3; ++i)
                {
                    if (t_to_objects_node_map[conn[i]] == NOT_SET_Index)
                    {
                        t_to_objects_node_map[conn[i]] = curr_v_id++;
                        objects_to_t_node_map.push_back(conn[i]);
                    }
                }
            }
        }

        objects.object_to_markers_map[object_id].assign(object_markers.begin(),
                                                        object_markers.end());
    }

    // Also ensure any surface triangle vertex not reached via object_to_regions_map is mapped
    for (const auto& conn : t.tetrahedralization.triangles)
    {
        for (int i = 0; i < 3; ++i)
        {
            if (t_to_objects_node_map[conn[i]] == NOT_SET_Index)
            {
                t_to_objects_node_map[conn[i]] = curr_v_id++;
                objects_to_t_node_map.push_back(conn[i]);
            }
        }
    }

    objects.t.vertices.resize(curr_v_id);
    for (std::size_t i = 0; i < objects_to_t_node_map.size(); ++i)
    {
        objects.t.vertices[i] = t.tetrahedralization.vertices[objects_to_t_node_map[i]];
    }

    objects.t.triangles.resize(t.tetrahedralization.triangles.size());
    objects.t.tri_to_marker_map.resize(t.tetrahedralization.triangles.size());
    for (std::size_t tri_id = 0; tri_id < t.tetrahedralization.triangles.size(); ++tri_id)
    {
        const auto& conn = t.tetrahedralization.triangles[tri_id];
        for (int i = 0; i < 3; ++i)
        {
            objects.t.triangles[tri_id][i] = t_to_objects_node_map[conn[i]];
        }
        objects.t.tri_to_marker_map[tri_id] =
          t.region_boundaries_info.triangle_to_boundary_id_map[tri_id];
    }

    update_edges_and_maps(objects.t);

    return objects;
}

void optimize_recovered_triangulation(RecoveredTriangulation& rt)
{
    if (rt.t.triangles.empty())
    {
        return;
    }

    optimize_edge_flips(rt.t, 5);

    const auto stats              = query_edge_length_stats(rt.t);
    const double mean_edge_length = stats.mean;

    auto seam_info = segment_seams(rt.t);

    optimize_tangential_smoothing(rt.t, seam_info, 10, true);

    split_large_edges(1.25 * mean_edge_length, rt.t, seam_info);

    optimize_edge_flips(rt.t, 5);

    split_constrained_edges(rt.t, seam_info);

    optimize_tangential_smoothing(rt.t, seam_info, 10, true);
}

Triangulation find_global_boundary(const RecoveredTriangulation& rt)
{
    if (rt.t.triangles.empty())
    {
        return {};
    }

    std::vector<int> marker_object_count(rt.t.n_markers, 0);
    for (const auto& markers : rt.object_to_markers_map)
    {
        for (const Marker marker : markers)
        {
            if (marker < marker_object_count.size())
            {
                marker_object_count[marker]++;
            }
        }
    }

    std::vector<Index> t_to_object_node_map(rt.t.vertices.size(), NOT_SET_Index);
    std::vector<Index> object_to_t_node_map;
    object_to_t_node_map.reserve(rt.t.vertices.size());

    Index curr_v_id = 0;
    Index curr_t_id = 0;
    for (std::size_t tri_id = 0; tri_id < rt.t.triangles.size(); ++tri_id)
    {
        if (marker_object_count[rt.t.tri_to_marker_map[tri_id]] == 1)
        {
            ++curr_t_id;
            const auto& conn = rt.t.triangles[tri_id];
            for (int i = 0; i < 3; ++i)
            {
                if (t_to_object_node_map[conn[i]] == NOT_SET_Index)
                {
                    t_to_object_node_map[conn[i]] = curr_v_id++;
                    object_to_t_node_map.push_back(conn[i]);
                }
            }
        }
    }

    Triangulation boundary;
    if (curr_t_id == 0)
    {
        return boundary;
    }

    boundary.vertices.resize(curr_v_id);
    for (std::size_t i = 0; i < object_to_t_node_map.size(); ++i)
    {
        boundary.vertices[i] = rt.t.vertices[object_to_t_node_map[i]];
    }

    boundary.triangles.resize(curr_t_id);
    boundary.tri_to_marker_map.resize(curr_t_id);
    curr_t_id = 0;
    for (std::size_t tri_id = 0; tri_id < rt.t.triangles.size(); ++tri_id)
    {
        if (marker_object_count[rt.t.tri_to_marker_map[tri_id]] == 1)
        {
            const auto& conn = rt.t.triangles[tri_id];
            for (int i = 0; i < 3; ++i)
            {
                boundary.triangles[curr_t_id][i] = t_to_object_node_map[conn[i]];
            }
            boundary.tri_to_marker_map[curr_t_id] = rt.t.tri_to_marker_map[tri_id];
            ++curr_t_id;
        }
    }

    {
        Eigen::MatrixXd vertices(static_cast<Eigen::Index>(boundary.vertices.size()), 3);
        for (std::size_t i = 0; i < boundary.vertices.size(); ++i)
        {
            for (int j = 0; j < 3; ++j)
            {
                vertices(static_cast<Eigen::Index>(i), j) = boundary.vertices[i][j];
            }
        }
        Eigen::MatrixXi triangles(static_cast<Eigen::Index>(boundary.triangles.size()), 3);
        for (std::size_t i = 0; i < boundary.triangles.size(); ++i)
        {
            for (int j = 0; j < 3; ++j)
            {
                triangles(static_cast<Eigen::Index>(i), j) =
                  static_cast<int>(boundary.triangles[i][j]);
            }
        }
        detail::orient_volume(vertices, triangles);
        for (std::size_t i = 0; i < boundary.triangles.size(); ++i)
        {
            for (int j = 0; j < 3; ++j)
            {
                boundary.triangles[i][j] =
                  static_cast<Index>(triangles(static_cast<Eigen::Index>(i), j));
            }
        }
    }

    update_edges_and_maps(boundary);

    return boundary;
}

std::vector<Triangulation> separate_objects(const RecoveredTriangulation& rt)
{
    const std::size_t n_objects = rt.object_to_markers_map.size();

    std::vector<std::vector<Index>> marker_to_tris_map(rt.t.n_markers);
    for (std::size_t tri_id = 0; tri_id < rt.t.triangles.size(); ++tri_id)
    {
        marker_to_tris_map[rt.t.tri_to_marker_map[tri_id]].push_back(static_cast<Index>(tri_id));
    }

    std::vector<Triangulation> objects;
    objects.reserve(n_objects);
    for (std::size_t object_id = 0; object_id < n_objects; ++object_id)
    {
        if (rt.object_to_markers_map[object_id].empty())
        {
            continue;
        }

        std::vector<Index> object_to_tris_map;
        for (const Marker marker_id : rt.object_to_markers_map[object_id])
        {
            for (const Index tri_id : marker_to_tris_map[marker_id])
            {
                object_to_tris_map.push_back(tri_id);
            }
        }
        if (object_to_tris_map.empty())
        {
            continue;
        }

        std::vector<Index> t_to_object_node_map(rt.t.vertices.size(), NOT_SET_Index);
        std::vector<Index> object_to_t_node_map;
        object_to_t_node_map.reserve(rt.t.vertices.size());

        Index curr_v_id = 0;
        for (const Index tri_id : object_to_tris_map)
        {
            const auto& conn = rt.t.triangles[tri_id];
            for (int i = 0; i < 3; ++i)
            {
                if (t_to_object_node_map[conn[i]] == NOT_SET_Index)
                {
                    t_to_object_node_map[conn[i]] = curr_v_id++;
                    object_to_t_node_map.push_back(conn[i]);
                }
            }
        }

        Triangulation object_t;
        object_t.vertices.resize(curr_v_id);
        for (std::size_t i = 0; i < object_to_t_node_map.size(); ++i)
        {
            object_t.vertices[i] = rt.t.vertices[object_to_t_node_map[i]];
        }

        Index curr_t_id = 0;
        object_t.triangles.resize(object_to_tris_map.size());
        object_t.tri_to_marker_map.resize(object_to_tris_map.size());
        for (const Index tri_id : object_to_tris_map)
        {
            const auto& conn = rt.t.triangles[tri_id];
            for (int i = 0; i < 3; ++i)
            {
                object_t.triangles[curr_t_id][i] = t_to_object_node_map[conn[i]];
            }
            object_t.tri_to_marker_map[curr_t_id] = rt.t.tri_to_marker_map[tri_id];
            ++curr_t_id;
        }
        assert(curr_t_id == object_to_tris_map.size());

        {
            Eigen::MatrixXd vertices(static_cast<Eigen::Index>(object_t.vertices.size()), 3);
            for (std::size_t i = 0; i < object_t.vertices.size(); ++i)
            {
                for (int j = 0; j < 3; ++j)
                {
                    vertices(static_cast<Eigen::Index>(i), j) = object_t.vertices[i][j];
                }
            }
            Eigen::MatrixXi triangles(static_cast<Eigen::Index>(object_t.triangles.size()), 3);
            for (std::size_t i = 0; i < object_t.triangles.size(); ++i)
            {
                for (int j = 0; j < 3; ++j)
                {
                    triangles(static_cast<Eigen::Index>(i), j) =
                      static_cast<int>(object_t.triangles[i][j]);
                }
            }
            detail::orient_volume(vertices, triangles);
            for (std::size_t i = 0; i < object_t.triangles.size(); ++i)
            {
                for (int j = 0; j < 3; ++j)
                {
                    object_t.triangles[i][j] =
                      static_cast<Index>(triangles(static_cast<Eigen::Index>(i), j));
                }
            }
        }

        update_edges_and_maps(object_t);

        objects.push_back(std::move(object_t));
    }

    return objects;
}

} // namespace mesh_arrange
