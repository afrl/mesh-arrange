#include "mesh_arrange/ConsistentSurfaces.h"

#include "detail/DebugVTK.h"
#include "detail/IGLHelpers.h"
#include "detail/MathHelpers.h"
#include "igl/copyleft/cgal/remesh_self_intersections.h"
#include "mesh_arrange/GlobalOptions.h"

#include <algorithm>
#include <cassert>
#include <limits>
#include <unordered_map>
#include <vector>

namespace mesh_arrange
{

namespace
{

std::vector<Index> find_facet_to_object_map(
  FacetID n_facets,
  const std::vector<std::vector<FacetID>>& object_to_facets_map)
{
    std::vector<Index> facet_to_object_map(n_facets, 0);
    for (std::size_t object_id = 0; object_id < object_to_facets_map.size(); ++object_id)
    {
        for (const FacetID facet_id : object_to_facets_map[object_id])
        {
            if (facet_id >= facet_to_object_map.size())
            {
                facet_to_object_map.resize(static_cast<std::size_t>(facet_id) + 1, 0);
            }
            facet_to_object_map[facet_id] = static_cast<Index>(object_id);
        }
    }
    return facet_to_object_map;
}

} // namespace

ConsistentSurfaces create_consistent_surface(const CombinedSurfaces& combined_inconsistent_surfaces)
{
    GlobalOptions::get().log_info("* Creating consistent surface by remeshing intersections");

    const auto& inconsistent_surface = combined_inconsistent_surfaces.inconsistent_surface_mesh;
    if (inconsistent_surface.triangles.empty())
    {
        return {};
    }

    FacetID n_original_facets = inconsistent_surface.n_markers;
    for (const auto& fvec : combined_inconsistent_surfaces.object_to_facets_map)
    {
        for (const FacetID fid : fvec)
        {
            n_original_facets = std::max(n_original_facets, fid + 1);
        }
    }

    const std::vector<Index> original_facet_to_object_map =
      find_facet_to_object_map(n_original_facets,
                               combined_inconsistent_surfaces.object_to_facets_map);

    Eigen::MatrixXd new_vertices;
    Eigen::MatrixXi new_triangles;
    Eigen::MatrixXi intersecting_face_pairs;
    Eigen::VectorXi original_birth_triangles;
    Eigen::VectorXi unique_vertex_map;

    igl::copyleft::cgal::RemeshSelfIntersectionsParam params;
    params.detect_only = false;
    params.first_only  = false;
    params.stitch_all  = true;

    Eigen::MatrixXd vertices(static_cast<Eigen::Index>(inconsistent_surface.vertices.size()), 3);
    for (std::size_t i = 0; i < inconsistent_surface.vertices.size(); ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            vertices(static_cast<Eigen::Index>(i), j) = inconsistent_surface.vertices[i][j];
        }
    }

    Eigen::MatrixXi triangles(static_cast<Eigen::Index>(inconsistent_surface.triangles.size()), 3);
    for (std::size_t i = 0; i < inconsistent_surface.triangles.size(); ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            triangles(static_cast<Eigen::Index>(i), j) =
              static_cast<int>(inconsistent_surface.triangles[i][j]);
        }
    }

    GlobalOptions::get().log_info(" - Remeshing self intersections");

    igl::copyleft::cgal::remesh_self_intersections(vertices,
                                                   triangles,
                                                   params,
                                                   new_vertices,
                                                   new_triangles,
                                                   intersecting_face_pairs,
                                                   original_birth_triangles,
                                                   unique_vertex_map);

    // Apply duplicate vertex mapping unique_vertex_map to new_triangles
    std::for_each(new_triangles.data(),
                  new_triangles.data() + new_triangles.size(),
                  [&unique_vertex_map](auto& a) { a = unique_vertex_map(a); });

    // Compact referenced vertices in ascending index order
    std::vector<int> old_to_new(static_cast<std::size_t>(new_vertices.rows()), -1);
    for (Eigen::Index f = 0; f < new_triangles.rows(); ++f)
    {
        for (int c = 0; c < 3; ++c)
        {
            old_to_new[static_cast<std::size_t>(new_triangles(f, c))] = 0;
        }
    }
    int n_used = 0;
    for (Eigen::Index v = 0; v < new_vertices.rows(); ++v)
    {
        if (old_to_new[static_cast<std::size_t>(v)] == 0)
        {
            old_to_new[static_cast<std::size_t>(v)] = n_used++;
        }
    }
    vertices.resize(n_used, 3);
    for (Eigen::Index v = 0; v < new_vertices.rows(); ++v)
    {
        const int nv = old_to_new[static_cast<std::size_t>(v)];
        if (nv >= 0)
        {
            vertices.row(nv) = new_vertices.row(v);
        }
    }
    triangles.resize(new_triangles.rows(), 3);
    for (Eigen::Index f = 0; f < new_triangles.rows(); ++f)
    {
        for (int c = 0; c < 3; ++c)
        {
            triangles(f, c) = old_to_new[static_cast<std::size_t>(new_triangles(f, c))];
        }
    }

    GlobalOptions::get().log_info(" - Updating data structures");

    auto r                           = detail::igl_to_triangulation(vertices, triangles);
    auto& consistent_surface_mesh    = r.t;
    auto& element_to_origin_tris_map = r.element_to_origin_tris_map;

    std::vector<std::vector<FacetID>> new_facet_to_original_facets_map;
    new_facet_to_original_facets_map.reserve(n_original_facets);
    std::unordered_map<std::size_t, FacetID> facet_hash_to_new_facet_map;
    for (FacetID i = 0; i < n_original_facets; ++i)
    {
        new_facet_to_original_facets_map.push_back({ i });
        const std::vector<FacetID> facet_set    = { i };
        const std::size_t facet_hash            = detail::hash_container(facet_set);
        facet_hash_to_new_facet_map[facet_hash] = i;
    }

    consistent_surface_mesh.tri_to_marker_map.resize(consistent_surface_mesh.triangles.size());
    for (std::size_t i = 0; i < consistent_surface_mesh.triangles.size(); ++i)
    {
        const auto& origin_tris = element_to_origin_tris_map[i];
        assert(!origin_tris.empty());

        if (origin_tris.size() == 1)
        {
            const FacetID original_facet_id = inconsistent_surface.tri_to_marker_map
              [static_cast<std::size_t>(original_birth_triangles(origin_tris[0]))];
            consistent_surface_mesh.tri_to_marker_map[i] = original_facet_id;
        }
        else
        {
            std::vector<FacetID> original_facet_ids;
            original_facet_ids.reserve(origin_tris.size());
            for (const Index origin_tri : origin_tris)
            {
                const FacetID original_facet_id = inconsistent_surface.tri_to_marker_map
                  [static_cast<std::size_t>(original_birth_triangles(origin_tri))];
                original_facet_ids.push_back(original_facet_id);
            }

            PriorityType maximum_priority = std::numeric_limits<PriorityType>::min();
            for (const FacetID original_facet_id : original_facet_ids)
            {
                const Index object_id = original_facet_to_object_map[original_facet_id];
                const PriorityType priority =
                  combined_inconsistent_surfaces.object_priorities[object_id];
                if (priority > maximum_priority)
                {
                    maximum_priority = priority;
                }
            }

            std::vector<FacetID> max_priority_facet_ids;
            max_priority_facet_ids.reserve(original_facet_ids.size());
            for (const FacetID original_facet_id : original_facet_ids)
            {
                const Index object_id = original_facet_to_object_map[original_facet_id];
                const PriorityType priority =
                  combined_inconsistent_surfaces.object_priorities[object_id];
                if (priority == maximum_priority)
                {
                    max_priority_facet_ids.push_back(original_facet_id);
                }
            }

            const std::size_t facet_hash = detail::hash_sorted_container(max_priority_facet_ids);
            const auto find_iter         = facet_hash_to_new_facet_map.find(facet_hash);
            if (find_iter == facet_hash_to_new_facet_map.end())
            {
                const FacetID new_facet_id =
                  static_cast<FacetID>(new_facet_to_original_facets_map.size());
                consistent_surface_mesh.tri_to_marker_map[i] = new_facet_id;
                facet_hash_to_new_facet_map[facet_hash]      = new_facet_id;
                new_facet_to_original_facets_map.push_back(std::move(original_facet_ids));
            }
            else
            {
                consistent_surface_mesh.tri_to_marker_map[i] = find_iter->second;
            }
        }
    }

    update_edges_and_maps(consistent_surface_mesh);

    if (GlobalOptions::get().write_intermediate_meshes)
    {
        detail::write_vtk_triangulation(
          GlobalOptions::get().output_dir + "/consistent_surface.vtk",
          consistent_surface_mesh);
    }

    return { std::move(consistent_surface_mesh), std::move(new_facet_to_original_facets_map) };
}

} // namespace mesh_arrange
