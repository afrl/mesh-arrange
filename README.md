# MeshArrange

**MeshArrange** is a standalone C++17 library for computing conforming 3D surface arrangements from collections of watertight triangular meshes, extracting complete non-manifold B-Rep topology (3D regions, 2D region-boundary patches, 1D seam curves, and 0D corner/junction vertices), applying priority-based boolean clipping, and optionally resolving multi-material overlaps (`K >= 2`) along smooth harmonic medial surfaces. It is designed for zero-leakage integration into mesh generators such as **Gmsh** as well as multi-material finite-element preprocessing pipelines.

---

## Capabilities

- **Exact Surface Arrangement (`CGAL::Epeck`):** Resolves all pairwise and self-intersections across arbitrary collections of closed triangular surfaces while merging coplanar facets with priority-aware facet provenance.
- **Priority & Boolean Clipping Rules:** Supports `ObjectType::SOLID`, `ObjectType::CLIP_INSIDE` (holes/voids), and `ObjectType::CLIP_OUTSIDE` (domain bounding cutters) across arbitrary integer priority levels.
- **Topological Feature & B-Rep Extraction:** Segments volumetric 3D regions (`BRepVolume`), consistently oriented 2D surface patches (`BRepFace`), 1D non-manifold seam arcs (`BRepCurve`), and 0D junction/loop-anchor corners (`BRepVertex`) ready for `discreteVertex` / `discreteEdge` / `discreteFace` / `discreteRegion` construction in Gmsh.
- **Multi-Material Interpenetration Removal (`K >= 2`):** Partitions equal-priority overlapping solid regions using a simplex-normalized clamped-cotangent harmonic field, canonical conforming tetrahedral subdivision, and 2-manifold star repair so every surviving object has a strictly disjoint, watertight manifold boundary.
- **Zero-Leakage Pure-STL Public API:** Public headers (`#include "mesh_arrange/MeshArrange.h"`) expose only standard C++17 STL types (`std::vector`, `std::array`, `std::function`, `uint32_t`)—no Eigen, CGAL, libigl, or TetGen headers leak into consumer translation units.
- **Pluggable Coarse PLC Tetrahedralizer:** Includes an optional built-in namespaced TetGen backend (`MESHARRANGE_WITH_TETGEN=ON`, isolated inside `mesh_arrange::tetgen_internal` to avoid symbol collisions with global `::tetgenio`) and accepts any custom `PlcTetrahedralizerFn` callback (such as Gmsh's constrained 3D Delaunay/HXT mesher with `MESHARRANGE_WITH_TETGEN=OFF`).

---

## Project Structure

- `include/mesh_arrange/` — Public pure-STL C++17 headers:
  - `MeshArrange.h` — High-level `arrange_surfaces(...)` API, `ArrangeOptions`, `ArrangeResult`, and B-Rep structs (`BRepVertex`, `BRepCurve`, `BRepFace`, `BRepVolume`).
  - `GmshModelBuilder.h` — Header-only `populate_gmsh_discrete_model(...)` adapter for Gmsh `GModel` discrete entities.
  - `PlcBackend.h` — Pluggable `PlcTetrahedralizerFn` interface for custom coarse 3D PLC tetrahedralizers.
  - Modular pipeline headers (`CombinedSurfaces.h`, `ConsistentSurfaces.h`, `FeatureIdentification.h`, `InterpenetrationRemoval.h`, `Refinement.h`, `QualityTriangulation.h`, `Triangulation.h`, `Tetrahedralization.h`, `Seams.h`, `SplitEdges.h`, `Optimize.h`).
- `src/` — Library implementation and private geometry helpers.
- `third_party/` — Vendored private header subset of `libigl` and optional namespaced TetGen 1.6.
- `tests/` — CTest unit and regression suites (`test_topology_and_seams`, `test_consistent_and_features`, `test_interpenetration_cases`, `test_custom_plc_and_gmsh_builder`).
- CMake package & target: `MeshArrange` / `MeshArrange::mesh_arrange`.

---

## Dependencies

All third-party mathematical dependencies are linked **`PRIVATE`** to `MeshArrange::mesh_arrange`:

- C++17 compiler (GCC 9+, Clang 10+, or MSVC 2019+) and CMake 3.16+
- [Eigen3](https://eigen.tuxfamily.org/) (>= 3.3)
- [CGAL](https://www.cgal.org/) (5.x+, header-only `CGAL::Epeck` exact predicates & constructions) with Boost headers, GMP (`libgmp`, `libgmpxx`), and MPFR (`libmpfr`)

---

## Quickstart Usage

```cpp
#include "mesh_arrange/MeshArrange.h"

std::vector<mesh_arrange::InputSurface> surfaces = { /* closed triangular meshes */ };

mesh_arrange::ArrangeOptions opts;
opts.remove_overlaps = true; // Partition equal-priority overlaps via harmonic medial sheets

mesh_arrange::ArrangeResult result = mesh_arrange::arrange_surfaces(surfaces, opts);
// result.brep_vertices, result.brep_curves, result.brep_faces, result.brep_volumes
// result.constituent_surfaces, result.exterior_boundary, result.surface_mesh
```

### Populating a Gmsh `GModel`

```cpp
#include "GModel.h"
#include "discreteVertex.h"
#include "discreteEdge.h"
#include "discreteFace.h"
#include "discreteRegion.h"
#include "MVertex.h"
#include "MLine.h"
#include "MTriangle.h"
#include "mesh_arrange/GmshModelBuilder.h"

auto result = mesh_arrange::arrange_surfaces(surfaces, opts);
mesh_arrange::populate_gmsh_discrete_model<
    GModel, discreteVertex, discreteEdge, discreteFace, discreteRegion,
    MVertex, MLine, MTriangle>(gmodel, result);
```

When embedding inside Gmsh without TetGen (`-DMESHARRANGE_WITH_TETGEN=OFF`), assign `opts.plc_tetrahedralizer` to a callback returning a coarse constrained tetrahedralization (`mesh_arrange::PlcTetOutput`) of the intersection-free `plc_surface`.

---

## Building & Running Tests

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DMESHARRANGE_BUILD_TESTS=ON
cmake --build build -j
ctest --test-dir build --output-on-failure
```

To build without the bundled TetGen backend (custom `PlcTetrahedralizerFn` only):

```bash
cmake -S . -B build-no-tetgen -DCMAKE_BUILD_TYPE=Release \
  -DMESHARRANGE_WITH_TETGEN=OFF -DMESHARRANGE_BUILD_TESTS=ON
cmake --build build-no-tetgen -j
ctest --test-dir build-no-tetgen --output-on-failure
```

---

## Public Release Information

Distribution Statement A. Approved for public release: distribution is unlimited. Case #: AFRL-2026-4513
