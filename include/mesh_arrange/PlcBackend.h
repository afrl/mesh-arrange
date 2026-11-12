#pragma once

#include "mesh_arrange/Triangulation.h"

#include <array>
#include <functional>
#include <vector>

namespace mesh_arrange
{

/// @brief Output of a coarse constrained Piecewise-Linear Complex (PLC) tetrahedralizer.
///
/// A custom `PlcTetrahedralizerFn` (e.g., backed by Gmsh's 3D constrained Delaunay / HXT mesher)
/// takes a self-intersection-free `Triangulation` and returns a conforming tetrahedralization
/// preserving the input surface facets without splitting them unnecessarily.
struct PlcTetOutput
{
    /// Output vertex coordinates (includes input vertices plus any interior Steiner vertices).
    std::vector<std::array<double, 3>> vertices;
    /// Output tetrahedra (0-based indices into `vertices`; orientation is normalized automatically).
    std::vector<std::array<Index, 4>> tetrahedra;
    /// Optional constrained surface sub-triangles present as faces of `tetrahedra`.
    /// Leave empty if the tetrahedralizer preserved `plc_surface.triangles` 1-to-1 without
    /// splitting any input surface triangle.
    std::vector<std::array<Index, 3>> surface_triangles;
    /// For each entry in `surface_triangles`, its 0-based parent triangle index in
    /// `plc_surface.triangles` (required when `surface_triangles` is non-empty).
    std::vector<Index> surface_tri_parent_id;
};

using PlcTetrahedralizerFn = std::function<PlcTetOutput(const Triangulation& plc_surface)>;

/// @brief Returns true if MeshArrange was compiled with its built-in namespaced TetGen backend.
bool has_builtin_tetgen_backend() noexcept;

/// @brief Default coarse PLC tetrahedralizer (uses built-in TetGen when `MESHARRANGE_WITH_TETGEN=ON`,
///        otherwise throws `std::runtime_error` instructing caller to supply a `PlcTetrahedralizerFn`).
PlcTetOutput default_plc_tetrahedralizer(const Triangulation& plc_surface);

} // namespace mesh_arrange
