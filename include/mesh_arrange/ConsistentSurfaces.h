#pragma once

#include "mesh_arrange/CombinedSurfaces.h"

#include <vector>

namespace mesh_arrange
{

/// @brief Conforming surface mesh after exact CGAL self-intersection remeshing and coplanar merging.
struct ConsistentSurfaces
{
    Triangulation consistent_surface_mesh;
    /// Maps each consistent facet ID to the 1-or-more original facet IDs that coincided there.
    std::vector<std::vector<FacetID>> consistent_facet_to_inconsistent_facets_map;
};

/// @brief Resolves all pairwise and self-intersections with exact predicates (`CGAL::Epeck`)
///        and merges duplicate/coplanar triangles with priority-aware facet provenance.
ConsistentSurfaces create_consistent_surface(const CombinedSurfaces& combined_inconsistent_surfaces);

} // namespace mesh_arrange
