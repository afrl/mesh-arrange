#pragma once

#include "mesh_arrange/FeatureIdentification.h"

#include <vector>

namespace mesh_arrange
{

struct RecoveredTriangulation
{
    Triangulation t;
    std::vector<std::vector<Marker>> object_to_markers_map;
};

/// @brief Extracts the conforming multi-material surface triangulation from `st`.
RecoveredTriangulation recover_triangulation(SegmentedTetrahedralization& t);

/// @brief Applies edge flips, feature-preserving tangential smoothing, and edge refinement.
void optimize_recovered_triangulation(RecoveredTriangulation& rt);

/// @brief Extracts the outward-oriented global exterior boundary triangulation.
Triangulation find_global_boundary(const RecoveredTriangulation& rt);

/// @brief Extracts outward-oriented closed watertight triangulations for each surviving object.
std::vector<Triangulation> separate_objects(const RecoveredTriangulation& rt);

} // namespace mesh_arrange
