#pragma once

#include "mesh_arrange/Seams.h"
#include "mesh_arrange/Triangulation.h"

namespace mesh_arrange
{

/// @brief Feature-preserving tangential smoothing with surface/seam closest-point projection.
void optimize_tangential_smoothing(Triangulation& t,
                                   const SeamInfo& seams,
                                   const int n_iterations,
                                   const bool allow_smoothing_on_seams = true);

/// @brief Valence-reducingmanifold edge flipping preserving markers and sharp creases.
void optimize_edge_flips(Triangulation& t, const int n_passes);

} // namespace mesh_arrange
