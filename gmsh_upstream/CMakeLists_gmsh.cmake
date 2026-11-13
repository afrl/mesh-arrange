# Drop-in CMake integration snippet for Gmsh's root `CMakeLists.txt`
#
# Assumptions inside a Gmsh source checkout:
#   1. Copy `mesh-arrange/` (BSD-3-Clause) to `contrib/mesh_arrange/`
#      (excluding `contrib/mesh_arrange/third_party/tetgen/`, which is NEVER
#       compiled or needed inside Gmsh because `src/mesh/meshArrange.cpp`
#       supplies Gmsh's native `meshGRegion` + `tetgenBR` PLC callback).
#   2. Copy `gmsh_upstream/src/mesh/meshArrange.{h,cpp}` -> `src/mesh/`
#      and `gmsh_upstream/src/plugin/ArrangeSurfaces.{h,cpp}` -> `src/plugin/`
#      (registering `ArrangeSurfaces` in `src/plugin/PluginManager.cpp`).
#   3. Add `opt(MESH_ARRANGE "Enable discrete surface arrangement" ${DEFAULT})`
#      near `opt(CGAL ...)` / `opt(EIGEN ...)` in Gmsh's root `CMakeLists.txt`.

if(ENABLE_MESH_ARRANGE AND ENABLE_MESH)
  find_package(CGAL QUIET)
  find_package(Boost 1.66 QUIET)
  if(CGAL_FOUND AND Boost_FOUND AND (ENABLE_EIGEN OR EIGEN_INC))
    set(MESHARRANGE_WITH_TETGEN OFF CACHE BOOL "" FORCE)
    set(MESHARRANGE_BUILD_TESTS OFF CACHE BOOL "" FORCE)
    add_subdirectory(contrib/mesh_arrange EXCLUDE_FROM_ALL)

    # Reuse Gmsh's bundled or system Eigen and GMP/MPFR configuration
    if(EIGEN_INC)
      target_include_directories(mesh_arrange PUBLIC ${EIGEN_INC})
    endif()
    if(HAVE_GMP AND GMP_INC AND GMP_LIB)
      target_include_directories(mesh_arrange PRIVATE ${GMP_INC})
      target_link_libraries(mesh_arrange PUBLIC ${GMP_LIB})
    endif()

    list(APPEND EXTERNAL_INCLUDES
      ${CMAKE_CURRENT_SOURCE_DIR}/contrib/mesh_arrange/include)
    list(APPEND EXTERNAL_LIBRARIES mesh_arrange)
    set_config_option(HAVE_MESH_ARRANGE "MeshArrange")
  else()
    message(STATUS "CGAL, Boost, or Eigen not found: disabling ENABLE_MESH_ARRANGE")
  endif()
endif()
