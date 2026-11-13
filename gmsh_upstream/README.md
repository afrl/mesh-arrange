# Upstream Gmsh Contribution Bundle (`Option C`)

This directory contains the drop-in files and integration snippets required to contribute `mesh-arrange` (**BSD-3-Clause**, AFRL Public Release `Case #: AFRL-2026-4513`) to upstream [Gmsh](https://gitlab.onelab.info/gmsh/gmsh) using the dual-entry-point architecture (**Option C**):

1. **Core `src/mesh/` Driver (`arrangeDiscreteSurfaces`):**
   - `src/mesh/meshArrange.h` $\rightarrow$ copy to Gmsh `src/mesh/meshArrange.h`
   - `src/mesh/meshArrange.cpp` $\rightarrow$ copy to Gmsh `src/mesh/meshArrange.cpp`
   - Uses Gmsh's native `meshGRegion` (`ALGO_3D_DELAUNAY` + `tetgenBR` with an embedded `discreteFace` inside a padded scratch `GRegion` and `renumber = 0`) as the `PlcTetrahedralizerFn` callback—**eliminating TetGen (AGPLv3) completely** (`MESHARRANGE_WITH_TETGEN=OFF`).
   - Populates `discreteVertex` (0D), `discreteEdge` (1D), `discreteFace` (2D, including `addEmbeddedEdge` / `addEmbeddedVertex`), and built-in `GEO` `Volume` (`GRegion`) 3D entities so `gmsh::model::mesh::generate(3)` tetrahedralizes all output regions out of the box.

2. **Multi-Language Public SDK API (`gmsh::model::mesh::arrangeSurfaces`):**
   - `api/gen_py_snippet.py` $\rightarrow$ insert into Gmsh `api/gen.py` inside `Module('mesh', ...)` and run `python3 gen.py` inside `api/` to regenerate `api/gmsh.h`, `api/gmshc.h`, `api/gmshc.cpp`, `api/gmsh.h_cwrap`, `api/gmsh.py`, `api/gmsh.jl`, `api/gmsh.f90`, `api/gmsh.mjs`, and `doc/texinfo/api.texi`.
   - `api/gmsh_cpp_snippet.cpp` $\rightarrow$ insert into Gmsh `src/common/gmsh.cpp` inside `gmsh::model::mesh`.
   - Mirrors `gmsh::model::occ::fragment(...)` by returning `(outDimTags, outDimTagsMap)` so scripts can map input constituents directly to Physical Volumes (`addPhysicalGroup(3, ...)`).

3. **Thin Interactive Plugin (`Plugin(ArrangeSurfaces)`):**
   - `src/plugin/ArrangeSurfaces.h` $\rightarrow$ copy to Gmsh `src/plugin/ArrangeSurfaces.h` (`GMSH_MeshPlugin`)
   - `src/plugin/ArrangeSurfaces.cpp` $\rightarrow$ copy to Gmsh `src/plugin/ArrangeSurfaces.cpp`
   - Register in `src/plugin/PluginManager.cpp`:
     ```cpp
     #include "ArrangeSurfaces.h"
     // Inside PluginManager::registerDefaultPlugins() under #if defined(HAVE_MESH):
     add(new GMSH_ArrangeSurfacesPlugin());
     ```
   - Exposes interactive FLTK GUI controls, `.geo` scripting (`Plugin(ArrangeSurfaces).Run;`), automatic 3D Physical Group propagation (`PropagatePhysicals = 1`), and CLI usage.

4. **CMake Build Wiring (`ENABLE_MESH_ARRANGE`):**
   - Copy `mesh-arrange` (`LICENSE`, `README.md`, `include/`, `src/`, `third_party/libigl/` — **without** `third_party/tetgen/`) to `contrib/MeshArrange/` alongside `CMakeLists_contrib_MeshArrange.cmake` as `contrib/MeshArrange/CMakeLists.txt`.
   - Add `#cmakedefine HAVE_MESH_ARRANGE` to `src/common/GmshConfig.h.in`, apply `CMakeLists_gmsh.cmake` in Gmsh's root `CMakeLists.txt`, add `meshArrange.cpp` to `src/mesh/CMakeLists.txt`, and add `ArrangeSurfaces.cpp` to `src/plugin/CMakeLists.txt`.

5. **Examples & Benchmarks:**
   - `demos/arrange_surfaces.py` $\rightarrow$ `examples/api/arrange_surfaces.py` (alongside `examples/api/arrange_surfaces.cpp` and `benchmarks/3d/arrange_surfaces_plugin.geo`).

---

## Ready-to-Post Upstream GitLab RFC (`gitlab.onelab.info/gmsh/gmsh/-/issues`)

**Title:** `RFC: Non-manifold discrete surface arrangement & multi-material overlap removal (gmsh::model::mesh::arrangeSurfaces + Plugin(ArrangeSurfaces))`

**Body:**
> Hello Christoph and the Gmsh team,
>
> We would like to contribute a standalone C++17 discrete surface arrangement and multi-material interpenetration-removal module (`mesh-arrange`, **BSD-3-Clause**, cleared for public release under AFRL Case `AFRL-2026-4513`) to Gmsh.
>
> ### Capabilities
> Given an arbitrary collection of closed 2D triangle shells or 3D volumes (STL/VTK micro-CT reconstructions, textile composites, granular packs, or independently meshed CAD assemblies):
> 1. **Exact Non-Manifold Arrangement:** Computes exact triangle-triangle intersections (via CGAL exact predicates/constructions + GMP) and segments 3D regions, 2D region-boundary patches, 1D non-manifold seam curves, and 0D junctions.
> 2. **Boolean Priorities & Clippers:** Supports integer constituent priorities plus `clipInside` (voids) and `clipOutside` (RVE windows).
> 3. **Gap-Free Multi-Material Overlap Removal ($K \ge 2$):** Partitions equal-priority overlapping solids along smooth harmonic medial surfaces ($\nabla^2 \phi_k = 0$) without leaving voids or slivers.
>
> ### Proposed Integration Architecture
> - **`contrib/mesh_arrange` (`ENABLE_MESH_ARRANGE`):** BSD-3-Clause library + minimal MPL-2.0 CGAL-intersection glue (49 header/inline files); reuses Gmsh's Eigen and GMP/MPFR (`MESHARRANGE_WITH_TETGEN=OFF`, **no TetGen dependency**).
> - **`src/mesh/meshArrange.{h,cpp}`:** Uses Gmsh's own `meshGRegion` (`ALGO_3D_DELAUNAY` + `tetgenBR`) for coarse PLC region recovery and populates `discreteVertex` / `discreteEdge` / `discreteFace` + `GEO` `Volume` entities ready for `gmsh::model::mesh::generate(3)`, with optional `createGeometryOfDiscreteEntities` surface parameterization.
> - **SDK API (`gmsh::model::mesh::arrangeSurfaces`) + `Plugin(ArrangeSurfaces)`:** Returns `(outDimTags, outDimTagsMap)` following `occ::fragment` conventions for Physical Group material assignment, alongside a thin `GMSH_PostPlugin` wrapper for `.geo` and GUI workflows.
>
> Please let us know if this layout and API signature align with your preferences before we open the Merge Request.
