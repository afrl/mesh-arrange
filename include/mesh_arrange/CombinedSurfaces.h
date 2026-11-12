#pragma once

#include "mesh_arrange/Common.h"
#include "mesh_arrange/Triangulation.h"

#include <array>
#include <vector>

namespace mesh_arrange
{

enum class ObjectType
{
    SOLID,        ///< Material volume; overlaps at equal priority are flagged or partitioned.
    CLIP_INSIDE,  ///< Boolean cutter: removes regions inside this closed surface at <= priority.
    CLIP_OUTSIDE  ///< Bounding domain cutter: removes regions outside this closed surface at <= priority.
};

using PriorityType = int;

/// @brief Single watertight input triangular surface passed to MeshArrange.
struct InputSurface
{
    /// Vertex coordinates of this closed triangular surface.
    std::vector<std::array<double, 3>> vertices;
    /// Outward-oriented (or consistently oriented) closed triangle connectivity.
    std::vector<std::array<Index, 3>> triangles;
    /// Optional per-triangle CAD/patch facet ID (0-based within this surface).
    /// If empty, all triangles on this surface are assigned local facet ID 0.
    std::vector<FacetID> triangle_facet_ids;
    /// Role of this surface in the arrangement (solid vs. interior/exterior clipper).
    ObjectType type = ObjectType::SOLID;
    /// Precedence level: higher priority cuts lower priority; clippers cut solids at equal priority.
    PriorityType priority = 0;
};

/// @brief Concatenated, un-intersected multi-object surface representation.
struct CombinedSurfaces
{
    Triangulation inconsistent_surface_mesh;
    std::vector<std::vector<Index>> object_to_surface_elements_map;
    std::vector<std::vector<FacetID>> object_to_facets_map;
    std::vector<ObjectType> object_types;
    std::vector<PriorityType> object_priorities;
};

/// @brief Combines a list of `InputSurface` objects into a unified `CombinedSurfaces` structure.
CombinedSurfaces combine_surfaces(const std::vector<InputSurface>& surfaces);

} // namespace mesh_arrange
