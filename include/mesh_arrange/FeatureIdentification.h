#pragma once

#include "mesh_arrange/CombinedSurfaces.h"
#include "mesh_arrange/ConsistentSurfaces.h"
#include "mesh_arrange/PlcBackend.h"
#include "mesh_arrange/Seams.h"
#include "mesh_arrange/Tetrahedralization.h"
#include "mesh_arrange/Triangulation.h"

#include <limits>
#include <string>
#include <utility>
#include <vector>

namespace mesh_arrange
{

using RegionID         = Marker;
using RegionBoundaryID = Marker;

constexpr RegionID NOT_SET_RegionID                 = std::numeric_limits<RegionID>::max();
constexpr RegionBoundaryID NOT_SET_RegionBoundaryID = std::numeric_limits<RegionBoundaryID>::max();

enum class RegionType
{
    KEPT,
    CLIPPED
};

struct RegionInfo
{
    RegionID n_regions = 0;
    std::vector<std::vector<Index>> region_to_bounding_triangles_map;
    std::vector<std::pair<RegionID, RegionID>> tri_to_bounded_regions_map;
    std::vector<RegionID> tet_to_region_map;
    std::vector<std::vector<Index>> region_to_tets_map;
    std::vector<FacetID> tri_to_facet_map;
};

struct ObjectRegionInfo
{
    std::vector<RegionID> regions_of_overlap;
    std::vector<std::vector<RegionID>> object_to_regions_map;
    std::vector<std::vector<RegionID>> region_to_objects_map;
    std::vector<RegionType> region_types;
};

struct SegmentRegionsResult
{
    Tetrahedralization tetrahedralization;
    RegionInfo region_info;
    ObjectRegionInfo object_region_info;
};

struct RegionBoundaryInfo
{
    RegionBoundaryID n_region_boundaries = 0;
    std::vector<RegionBoundaryID> triangle_to_boundary_id_map;
    std::vector<std::vector<Index>> region_boundary_to_triangles_map;
    std::vector<std::vector<RegionBoundaryID>> region_id_to_region_boundary_ids_map;
    std::vector<std::vector<RegionID>> region_boundary_id_to_region_ids_map;
};

struct ObjectRegionBoundaryInfo
{
    std::vector<std::vector<RegionBoundaryID>> object_to_external_boundary_ids_map;
    std::vector<std::vector<RegionID>> boundary_id_to_objects_map;
};

struct SegmentRegionBoundariesResult
{
    RegionBoundaryInfo region_boundaries_info;
    ObjectRegionBoundaryInfo object_region_boundary_info;
    std::vector<FacetID> region_boundary_id_to_consistent_facet_map;
};

struct SegmentedTetrahedralization
{
    Tetrahedralization tetrahedralization;
    RegionInfo region_info;
    ObjectRegionInfo object_region_info;
    RegionBoundaryInfo region_boundaries_info;
    ObjectRegionBoundaryInfo object_region_boundary_info;
    SeamInfo seams_info;
};

void remove_region_boundaries(const std::vector<RegionBoundaryID>& region_boundaries_to_remove,
                              Tetrahedralization& tetrahedralization,
                              RegionInfo& region_info,
                              ObjectRegionInfo& object_region_info,
                              RegionBoundaryInfo& rbi,
                              ObjectRegionBoundaryInfo& orbi);

SegmentRegionsResult segment_regions(
  const CombinedSurfaces& inconsistent_combined_surfaces,
  const ConsistentSurfaces& consistent_surfaces,
  const PlcTetrahedralizerFn& plc_tetrahedralizer = {});

SegmentRegionBoundariesResult segment_region_boundaries(
  const CombinedSurfaces& inconsistent_combined_surfaces,
  const ConsistentSurfaces& consistent_surfaces,
  const Tetrahedralization& tetrahedralization,
  const RegionInfo& region_info,
  const ObjectRegionInfo& object_region_info);

SeamInfo segment_seams(const Tetrahedralization& tetrahedralization,
                       const RegionBoundaryInfo& region_boundary_info);

SegmentedTetrahedralization identify_features(
  const CombinedSurfaces& inconsistent_combined_surfaces,
  ConsistentSurfaces& consistent_surfaces,
  const PlcTetrahedralizerFn& plc_tetrahedralizer = {});

void write_features_to_vtk(const SegmentedTetrahedralization& st,
                           const std::string prefix  = "",
                           const std::string postfix = "");

} // namespace mesh_arrange
