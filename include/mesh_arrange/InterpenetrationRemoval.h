#pragma once

#include "mesh_arrange/CombinedSurfaces.h"
#include "mesh_arrange/ConsistentSurfaces.h"
#include "mesh_arrange/FeatureIdentification.h"

namespace mesh_arrange
{

/// @brief Resolves multi-material overlaps (`K >= 2`) via simplex harmonic fields and canonical
///        conforming tetrahedral subdivision, replacing overlapping regions with watertight
///        medial interface sheets.
void remove_interpenetrations(const CombinedSurfaces& inconsistent_combined_surfaces,
                              const ConsistentSurfaces& consistent_surfaces,
                              SegmentedTetrahedralization& st);

} // namespace mesh_arrange
