#include "mesh_arrange/CombinedSurfaces.h"

#include "detail/DebugVTK.h"
#include "mesh_arrange/GlobalOptions.h"

#include <algorithm>
#include <numeric>
#include <set>
#include <stdexcept>

namespace mesh_arrange
{

CombinedSurfaces combine_surfaces(const std::vector<InputSurface>& surfaces)
{
    CombinedSurfaces combined_surfaces;
    combined_surfaces.object_to_surface_elements_map.reserve(surfaces.size());
    combined_surfaces.object_to_facets_map.reserve(surfaces.size());
    combined_surfaces.object_types.reserve(surfaces.size());
    combined_surfaces.object_priorities.reserve(surfaces.size());

    FacetID facet_count = 0;
    for (const auto& surface : surfaces)
    {
        if (!surface.triangle_facet_ids.empty() &&
            surface.triangle_facet_ids.size() != surface.triangles.size())
        {
            throw std::invalid_argument(
              "InputSurface::triangle_facet_ids must be empty or match triangles.size().");
        }

        const Index surf_start_index =
          static_cast<Index>(combined_surfaces.inconsistent_surface_mesh.triangles.size());
        const std::size_t surf_size = surface.triangles.size();

        combined_surfaces.object_to_surface_elements_map.emplace_back(surf_size, 0);
        std::iota(combined_surfaces.object_to_surface_elements_map.back().begin(),
                  combined_surfaces.object_to_surface_elements_map.back().end(),
                  surf_start_index);

        const Index v_offset =
          static_cast<Index>(combined_surfaces.inconsistent_surface_mesh.vertices.size());
        for (const auto& pt : surface.vertices)
        {
            combined_surfaces.inconsistent_surface_mesh.vertices.push_back(pt);
        }

        std::set<FacetID> facets_for_surface;
        FacetID max_local_facet = 0;
        for (std::size_t i = 0; i < surf_size; ++i)
        {
            const auto& conn = surface.triangles[i];
            combined_surfaces.inconsistent_surface_mesh.triangles.push_back(
              { conn[0] + v_offset, conn[1] + v_offset, conn[2] + v_offset });

            const FacetID local_facet =
              surface.triangle_facet_ids.empty() ? FacetID{ 0 } : surface.triangle_facet_ids[i];
            max_local_facet          = std::max(max_local_facet, local_facet);
            const FacetID global_fid = local_facet + facet_count;
            facets_for_surface.insert(global_fid);
            combined_surfaces.inconsistent_surface_mesh.tri_to_marker_map.push_back(global_fid);
        }

        combined_surfaces.object_to_facets_map.emplace_back(facets_for_surface.begin(),
                                                            facets_for_surface.end());
        if (surf_size > 0)
        {
            facet_count += (max_local_facet + 1);
        }

        combined_surfaces.object_types.push_back(surface.type);
        combined_surfaces.object_priorities.push_back(surface.priority);
    }

    update_edges_and_maps(combined_surfaces.inconsistent_surface_mesh);

    if (GlobalOptions::get().write_intermediate_meshes)
    {
        detail::write_vtk_triangulation(
          GlobalOptions::get().output_dir + "/combined_inconsistent_surface.vtk",
          combined_surfaces.inconsistent_surface_mesh);
    }

    return combined_surfaces;
}

} // namespace mesh_arrange
