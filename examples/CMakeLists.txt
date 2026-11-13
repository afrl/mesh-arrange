cmake_minimum_required(VERSION 3.16)
project(MeshArrangeExamples LANGUAGES CXX)

if(NOT TARGET MeshArrange::mesh_arrange)
    # Allow standalone configuration: cmake -S examples -B build-examples
    add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/.." "${CMAKE_CURRENT_BINARY_DIR}/mesh_arrange_lib")
endif()

set(MESHARRANGE_DEMO_TARGETS
    demo_brep_arrangement
    demo_priority_clipping
    demo_harmonic_medial
    demo_multiway_junction
    demo_gmsh_rve_workflow
)

foreach(demo_target IN LISTS MESHARRANGE_DEMO_TARGETS)
    add_executable(${demo_target} "${demo_target}.cpp")
    target_compile_features(${demo_target} PRIVATE cxx_std_17)
    target_include_directories(${demo_target} PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}")
    target_link_libraries(${demo_target} PRIVATE MeshArrange::mesh_arrange)
endforeach()
