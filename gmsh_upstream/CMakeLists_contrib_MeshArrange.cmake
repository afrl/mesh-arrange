# Gmsh - Copyright (C) 1997-2026 C. Geuzaine, J.-F. Remacle
#
# See the LICENSE.txt file in the Gmsh root directory for license information.
# Please report all issues on https://gitlab.onelab.info/gmsh/gmsh/issues.

set(SRC
  src/CombinedSurfaces.cpp
  src/ConsistentSurfaces.cpp
  src/FeatureIdentification.cpp
  src/InterpenetrationRemoval.cpp
  src/MeshArrange.cpp
  src/Optimize.cpp
  src/PlcBackend.cpp
  src/QualityTriangulation.cpp
  src/Refinement.cpp
  src/Seams.cpp
  src/SplitEdges.cpp
  src/Tetrahedralization.cpp
  src/Triangulation.cpp
  src/VertexMap.cpp
  src/detail/DebugVTK.cpp
  src/detail/IGLHelpers.cpp
)

file(GLOB_RECURSE HDR RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} include/*.h src/*.h)

add_library(mesh_arrange OBJECT ${SRC} ${HDR})
set_target_properties(mesh_arrange PROPERTIES
  POSITION_INDEPENDENT_CODE ON
  CXX_STANDARD 17
  CXX_STANDARD_REQUIRED ON
  CXX_VISIBILITY_PRESET hidden
)

if(HAVE_CONTRIB_EIGEN)
  set(_MA_EIGEN_INC ${PROJECT_SOURCE_DIR}/contrib/eigen)
else()
  set(_MA_EIGEN_INC ${EIGEN_INC})
endif()

target_include_directories(mesh_arrange PRIVATE
  ${CMAKE_CURRENT_SOURCE_DIR}/include
  ${CMAKE_CURRENT_SOURCE_DIR}/src
  ${CMAKE_CURRENT_SOURCE_DIR}/third_party/libigl/include
  ${_MA_EIGEN_INC}
  ${MESHARRANGE_EXTERNAL_INCLUDES}
)

target_compile_definitions(mesh_arrange PRIVATE CGAL_HEADER_ONLY=1)

if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
  target_compile_options(mesh_arrange PRIVATE -w)
endif()

set(GMSH_SRC ${GMSH_SRC};$<TARGET_OBJECTS:mesh_arrange> PARENT_SCOPE)
set(GMSH_DIRS ${GMSH_DIRS};contrib/MeshArrange PARENT_SCOPE)
