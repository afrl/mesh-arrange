#pragma once

#include "mesh_arrange/CombinedSurfaces.h"
#include "mesh_arrange/ConsistentSurfaces.h"
#include "mesh_arrange/FeatureIdentification.h"

namespace mesh_arrange
{

struct RecoveredTriangulation;

/// @brief Resolves multi-material overlaps (`K >= 2`) via simplex harmonic fields and canonical
///        conforming tetrahedral subdivision, replacing overlapping regions with watertight
///        medial interface sheets.
void remove_interpenetrations(const CombinedSurfaces& inconsistent_combined_surfaces,
                              const ConsistentSurfaces& consistent_surfaces,
                              SegmentedTetrahedralization& st);

/// @brief Coupled 3D volume-preserving Taubin fairing of synthesized medial interface sheets
///        and interior multi-sheet junction curves inside `st`.
void smooth_medial_surfaces(const ConsistentSurfaces& consistent_surfaces,
                            SegmentedTetrahedralization& st);

/// @brief Medial-only short-edge collapse, Delaunay/quality edge flips, and self-intersection-safe
///        Taubin polish on recovered surface triangulation `rt`.
void polish_recovered_medial_surfaces(const ConsistentSurfaces& consistent_surfaces,
                                      const SegmentedTetrahedralization& st,
                                      RecoveredTriangulation& rt);

} // namespace mesh_arrange
