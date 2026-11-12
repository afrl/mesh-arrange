#pragma once

#include "mesh_arrange/CombinedSurfaces.h"
#include "mesh_arrange/ConsistentSurfaces.h"
#include "mesh_arrange/FeatureIdentification.h"
#include "mesh_arrange/GlobalOptions.h"
#include "mesh_arrange/InterpenetrationRemoval.h"
#include "mesh_arrange/Optimize.h"
#include "mesh_arrange/PlcBackend.h"
#include "mesh_arrange/QualityTriangulation.h"
#include "mesh_arrange/Refinement.h"
#include "mesh_arrange/Seams.h"
#include "mesh_arrange/SplitEdges.h"
#include "mesh_arrange/Tetrahedralization.h"
#include "mesh_arrange/Triangulation.h"

#include <array>
#include <string>
#include <utility>
#include <vector>

namespace mesh_arrange
{

/// @brief Options controlling `arrange_surfaces(...)`.
struct ArrangeOptions
{
    /// When true, equal-priority overlapping solid volumes (`K >= 2`) are partitioned along
    /// smooth harmonic medial surfaces so all output regions become strictly disjoint.
    /// When false, overlap regions are preserved as distinct regions with multiple `owner_objects`.
    bool remove_overlaps = true;

    /// When true, applies feature-preserving edge flips and tangential smoothing to the
    /// recovered non-manifold surface triangulation.
    bool optimize_surfaces = false;

    /// When true, splits non-seam surface edges whose endpoints both lie on seams.
    bool split_constrained_seam_edges = false;

    /// If > 0.0, splits surface edges longer than `max_edge_length` while preserving seams.
    double max_edge_length = 0.0;

    /// Enables progress logging (writes to `log_callback` if set, otherwise `std::cout`).
    bool verbose = false;

    /// Optional message sink (e.g., forwarding to Gmsh's `Msg::Info` / `Msg::Debug`).
    /// When set, `mesh_arrange` routes all diagnostic output through this callback instead of `std::cout`.
    LogCallback log_callback = {};

    /// Runs bidirectional topological invariant checks on edges, faces, and tetrahedra.
    bool check_topology = false;

    /// Writes intermediate VTK diagnostic files into `debug_output_dir`.
    bool write_debug_vtk         = false;
    std::string debug_output_dir = ".";

    /// Optional custom coarse PLC tetrahedralizer (e.g., Gmsh's constrained 3D Delaunay/HXT).
    /// When empty (`nullptr`), uses `default_plc_tetrahedralizer` (built-in namespaced TetGen).
    PlcTetrahedralizerFn plc_tetrahedralizer = {};
};

/// @brief 0D topological corner / seam-junction vertex for B-Rep model construction (Gmsh `GVertex`).
struct BRepVertex
{
    /// 0-based vertex ID in `ArrangeResult::surface_mesh.vertices`.
    Index global_vertex_id = NOT_SET_Index;
    /// 3D coordinates.
    std::array<double, 3> point = { 0.0, 0.0, 0.0 };
    /// 0-based indices into `ArrangeResult::brep_curves` incident to this vertex.
    std::vector<Index> incident_curves;
};

/// @brief 1D manifold seam curve arc between two `BRepVertex` corners or forming a closed loop (Gmsh `GEdge`).
struct BRepCurve
{
    /// Underlying `SeamID` from `SeamInfo`.
    SeamID seam_id = NOT_SET_SeamID;
    /// Ordered polyline of vertex IDs in `ArrangeResult::surface_mesh.vertices`.
    /// For closed loops (`is_closed_loop == true`), `polyline_vertices.front() == polyline_vertices.back()`.
    std::vector<Index> polyline_vertices;
    /// Start `BRepVertex` index in `ArrangeResult::brep_vertices` (`NOT_SET_Index` only if loop has no corner,
    /// though closed loops also anchor a single `BRepVertex` so `start_brep_vertex == end_brep_vertex`).
    Index start_brep_vertex = NOT_SET_Index;
    /// End `BRepVertex` index in `ArrangeResult::brep_vertices`.
    Index end_brep_vertex = NOT_SET_Index;
    /// True if this curve forms a periodic loop (`start_brep_vertex == end_brep_vertex`).
    bool is_closed_loop = false;
    /// 0-based indices into `ArrangeResult::brep_faces` incident to this curve.
    std::vector<RegionBoundaryID> incident_faces;
};

/// @brief 2D manifold surface patch separating a specific pair of 3D regions/facets (Gmsh `GFace`).
struct BRepFace
{
    /// 0-based `RegionBoundaryID` matching `surface_mesh.tri_to_marker_map`.
    RegionBoundaryID boundary_id = NOT_SET_RegionBoundaryID;
    /// Consistently oriented triangles (indices into `ArrangeResult::surface_mesh.vertices`).
    /// When `bounded_regions.first != NOT_SET_RegionID`, triangle normals point OUTWARD from
    /// `bounded_regions.first` (and inward to `bounded_regions.second` when present).
    std::vector<std::array<Index, 3>> triangles;
    /// Corresponding triangle indices in `ArrangeResult::surface_mesh.triangles`.
    std::vector<Index> surface_tri_ids;
    /// Pair of bounded 3D region IDs `{r0, r1}` (`r1 == NOT_SET_RegionID` on exterior boundary).
    std::pair<RegionID, RegionID> bounded_regions = { NOT_SET_RegionID, NOT_SET_RegionID };
    /// Signed bounding curves `{curve_id, +1 or -1}` along the outer/inner boundary of this face patch.
    std::vector<std::pair<Index, int>> signed_bounding_curves;
    /// Internal / embedded seam curves (`curve_id` in `brep_curves`) whose two sides both belong to this face.
    std::vector<Index> embedded_curves;
    /// Isolated `BRepVertex` indices (`brep_vertices`) embedded in the interior of this face.
    std::vector<Index> embedded_vertices;
    /// Input `object_id`s whose original boundary coincided with this patch (empty for medial sheets).
    std::vector<Index> origin_objects;
};

/// @brief 3D volumetric region bounded by oriented `BRepFace` shells (Gmsh `GRegion`).
struct BRepVolume
{
    RegionID region_id = NOT_SET_RegionID;
    /// Input solid `object_id`s that contain this volume (size == 1 when `remove_overlaps == true`).
    std::vector<Index> owner_objects;
    /// True if `owner_objects.size() > 1` (only possible when `remove_overlaps == false`).
    bool is_overlap = false;
    /// Signed bounding faces `{boundary_id, +1 or -1}` where `+1` means `BRepFace::triangles`
    /// are oriented outward from this volume and `-1` means inward.
    std::vector<std::pair<RegionBoundaryID, int>> signed_bounding_faces;
};

/// @brief Complete arrangement, B-Rep topology, and separated constituent meshes returned by `arrange_surfaces`.
struct ArrangeResult
{
    /// Global conforming non-manifold surface triangulation (`tri_to_marker_map` = `RegionBoundaryID`).
    Triangulation surface_mesh;

    /// Explicit B-Rep hierarchy ready for Gmsh `discreteVertex` / `discreteEdge` / `discreteFace` / `discreteRegion`.
    std::vector<BRepVertex> brep_vertices;
    std::vector<BRepCurve> brep_curves;
    std::vector<BRepFace> brep_faces;
    std::vector<BRepVolume> brep_volumes;

    /// Outward-oriented closed watertight surface mesh per surviving input solid object.
    std::vector<Triangulation> constituent_surfaces;
    /// Input `object_id` corresponding to each entry in `constituent_surfaces`.
    std::vector<Index> constituent_object_ids;

    /// Outward-oriented global exterior domain boundary triangulation.
    Triangulation exterior_boundary;

    /// Underlying segmented tetrahedralization and region/boundary/seam tables.
    SegmentedTetrahedralization segmented_tet;
};

/// @brief Intersects a set of watertight triangular surfaces, applies boolean priority/clipping rules,
///        segments 3D regions, 2D region boundaries, and 1D/0D non-manifold seams, optionally removes
///        multi-material overlaps via harmonic medial surfaces, and returns a complete B-Rep arrangement.
ArrangeResult arrange_surfaces(const std::vector<InputSurface>& surfaces,
                               const ArrangeOptions& options = {});

} // namespace mesh_arrange
