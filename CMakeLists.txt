cmake_minimum_required(VERSION 3.16)

project(MeshArrange VERSION 0.1.0 LANGUAGES CXX)

# ==============================================================================
# Options
# ==============================================================================
option(MESHARRANGE_WITH_TETGEN
       "Build with bundled namespaced TetGen backend for default coarse PLC tetrahedralization"
       ON)

set(_mesh_arrange_is_top_level OFF)
if(CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)
    set(_mesh_arrange_is_top_level ON)
endif()

option(MESHARRANGE_BUILD_TESTS
       "Build MeshArrange unit and regression tests"
       ${_mesh_arrange_is_top_level})
option(MESHARRANGE_BUILD_EXAMPLES
       "Build MeshArrange demonstration executables in examples/"
       OFF)
option(MESHARRANGE_INSTALL
       "Generate install and CMake package export rules for MeshArrange"
       ${_mesh_arrange_is_top_level})

set(MESHARRANGE_CGAL_INCLUDE_DIR "" CACHE PATH
    "Optional path to header-only CGAL include directory if find_package(CGAL) is unavailable")

# ==============================================================================
# Private Dependencies (Eigen3, CGAL + Boost + GMP/MPFR, Threads)
# ==============================================================================
find_package(Threads REQUIRED)

if(NOT TARGET Eigen3::Eigen)
    find_package(Eigen3 3.3 REQUIRED NO_MODULE)
endif()

if(NOT TARGET Boost::headers AND NOT TARGET Boost::boost)
    find_package(Boost REQUIRED)
endif()
if(TARGET Boost::headers)
    set(_mesh_arrange_boost_target Boost::headers)
elseif(TARGET Boost::boost)
    set(_mesh_arrange_boost_target Boost::boost)
endif()

set(_mesh_arrange_cgal_libs "")
set(_mesh_arrange_cgal_includes "")
set(_mesh_arrange_cgal_defs "")

if(TARGET CGAL::CGAL)
    list(APPEND _mesh_arrange_cgal_libs $<BUILD_INTERFACE:CGAL::CGAL>)
else()
    if(MESHARRANGE_CGAL_INCLUDE_DIR)
        list(APPEND _mesh_arrange_cgal_includes "${MESHARRANGE_CGAL_INCLUDE_DIR}")
        list(APPEND _mesh_arrange_cgal_defs CGAL_HEADER_ONLY=1)
    else()
        find_package(CGAL QUIET)
        if(TARGET CGAL::CGAL)
            list(APPEND _mesh_arrange_cgal_libs $<BUILD_INTERFACE:CGAL::CGAL>)
        elseif(CGAL_DIR AND EXISTS "${CGAL_DIR}/include/CGAL/version.h")
            list(APPEND _mesh_arrange_cgal_includes "${CGAL_DIR}/include")
            list(APPEND _mesh_arrange_cgal_defs CGAL_HEADER_ONLY=1)
        else()
            find_package(CGAL REQUIRED)
            list(APPEND _mesh_arrange_cgal_libs $<BUILD_INTERFACE:CGAL::CGAL>)
        endif()
    endif()
endif()

# Link GMP/MPFR when available; otherwise fall back cleanly to header-only Boost.Multiprecision
if(NOT WIN32)
    find_path(MESHARRANGE_GMP_INCLUDE_DIR NAMES gmp.h)
    find_path(MESHARRANGE_GMPXX_INCLUDE_DIR NAMES gmpxx.h)
    find_library(MESHARRANGE_GMP_LIBRARY NAMES gmp)
    find_library(MESHARRANGE_GMPXX_LIBRARY NAMES gmpxx)
    find_path(MESHARRANGE_MPFR_INCLUDE_DIR NAMES mpfr.h)
    find_library(MESHARRANGE_MPFR_LIBRARY NAMES mpfr)

    if(MESHARRANGE_GMP_INCLUDE_DIR AND MESHARRANGE_GMP_LIBRARY AND
       MESHARRANGE_MPFR_INCLUDE_DIR AND MESHARRANGE_MPFR_LIBRARY)
        list(APPEND _mesh_arrange_cgal_includes
             "${MESHARRANGE_GMP_INCLUDE_DIR}"
             "${MESHARRANGE_MPFR_INCLUDE_DIR}")
        list(APPEND _mesh_arrange_cgal_libs
             "${MESHARRANGE_MPFR_LIBRARY}"
             "${MESHARRANGE_GMP_LIBRARY}")
        if(MESHARRANGE_GMPXX_INCLUDE_DIR AND MESHARRANGE_GMPXX_LIBRARY)
            list(APPEND _mesh_arrange_cgal_includes "${MESHARRANGE_GMPXX_INCLUDE_DIR}")
            list(APPEND _mesh_arrange_cgal_libs "${MESHARRANGE_GMPXX_LIBRARY}")
            list(APPEND _mesh_arrange_cgal_defs CGAL_USE_GMPXX=1)
        endif()
    else()
        list(APPEND _mesh_arrange_cgal_defs CGAL_DISABLE_GMP=1 CGAL_USE_BOOST_MP=1)
    endif()
else()
    list(APPEND _mesh_arrange_cgal_defs CGAL_DISABLE_GMP=1 CGAL_USE_BOOST_MP=1)
endif()

# ==============================================================================
# Library Target: MeshArrange::mesh_arrange
# ==============================================================================
set(MESHARRANGE_PUBLIC_HEADERS
    include/mesh_arrange/Common.h
    include/mesh_arrange/GlobalOptions.h
    include/mesh_arrange/Triangulation.h
    include/mesh_arrange/Tetrahedralization.h
    include/mesh_arrange/VertexMap.h
    include/mesh_arrange/Seams.h
    include/mesh_arrange/SplitEdges.h
    include/mesh_arrange/Optimize.h
    include/mesh_arrange/PlcBackend.h
    include/mesh_arrange/CombinedSurfaces.h
    include/mesh_arrange/ConsistentSurfaces.h
    include/mesh_arrange/FeatureIdentification.h
    include/mesh_arrange/InterpenetrationRemoval.h
    include/mesh_arrange/Refinement.h
    include/mesh_arrange/QualityTriangulation.h
    include/mesh_arrange/MeshArrange.h
    include/mesh_arrange/GmshModelBuilder.h
)

set(MESHARRANGE_SOURCES
    src/detail/IGLHelpers.cpp
    src/detail/DebugVTK.cpp
    src/Triangulation.cpp
    src/Tetrahedralization.cpp
    src/VertexMap.cpp
    src/Seams.cpp
    src/SplitEdges.cpp
    src/Optimize.cpp
    src/PlcBackend.cpp
    src/CombinedSurfaces.cpp
    src/ConsistentSurfaces.cpp
    src/FeatureIdentification.cpp
    src/InterpenetrationRemoval.cpp
    src/Refinement.cpp
    src/QualityTriangulation.cpp
    src/MeshArrange.cpp
)

add_library(mesh_arrange STATIC ${MESHARRANGE_SOURCES} ${MESHARRANGE_PUBLIC_HEADERS})
add_library(MeshArrange::mesh_arrange ALIAS mesh_arrange)

target_compile_features(mesh_arrange PUBLIC cxx_std_17)
set_target_properties(mesh_arrange PROPERTIES POSITION_INDEPENDENT_CODE ON)

target_include_directories(mesh_arrange
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include>
    PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/src
        ${CMAKE_CURRENT_SOURCE_DIR}/third_party/libigl/include
        ${_mesh_arrange_cgal_includes}
)

target_compile_definitions(mesh_arrange
    PRIVATE
        ${_mesh_arrange_cgal_defs}
)

target_link_libraries(mesh_arrange
    PRIVATE
        $<BUILD_INTERFACE:Eigen3::Eigen>
        Threads::Threads
        $<BUILD_INTERFACE:${_mesh_arrange_boost_target}>
        ${_mesh_arrange_cgal_libs}
)

if(MESHARRANGE_WITH_TETGEN)
    target_sources(mesh_arrange
        PRIVATE
            third_party/tetgen/tetgen.cpp
            third_party/tetgen/predicates.cpp
    )
    target_include_directories(mesh_arrange PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/third_party/tetgen")
    target_compile_definitions(mesh_arrange PRIVATE TETLIBRARY MESHARRANGE_ENABLE_TETGEN=1)

    if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
        set_source_files_properties(
            third_party/tetgen/tetgen.cpp
            PROPERTIES COMPILE_OPTIONS "-w"
        )
        # TetGen predicates require strict IEEE-754 floating-point evaluation without contraction
        set_source_files_properties(
            third_party/tetgen/predicates.cpp
            PROPERTIES COMPILE_OPTIONS "-O0;-ffp-contract=off;-w"
        )
    elseif(MSVC)
        set_source_files_properties(
            third_party/tetgen/predicates.cpp
            PROPERTIES COMPILE_OPTIONS "/Od;/fp:strict"
        )
    endif()
endif()

# ==============================================================================
# Tests
# ==============================================================================
if(MESHARRANGE_BUILD_TESTS)
    include(CTest)
    enable_testing()

    function(mesh_arrange_add_test test_name test_src)
        add_executable(${test_name} ${test_src})
        target_compile_features(${test_name} PRIVATE cxx_std_17)
        target_link_libraries(${test_name} PRIVATE MeshArrange::mesh_arrange)
        add_test(NAME ${test_name} COMMAND ${test_name})
    endfunction()

    mesh_arrange_add_test(test_topology_and_seams tests/test_topology_and_seams.cpp)
    mesh_arrange_add_test(test_custom_plc_and_gmsh_builder tests/test_custom_plc_and_gmsh_builder.cpp)

    if(MESHARRANGE_WITH_TETGEN)
        mesh_arrange_add_test(test_consistent_and_features tests/test_consistent_and_features.cpp)
        mesh_arrange_add_test(test_interpenetration_cases tests/test_interpenetration_cases.cpp)
    endif()
endif()

# ==============================================================================
# Examples / Gallery Demos
# ==============================================================================
if(MESHARRANGE_BUILD_EXAMPLES AND MESHARRANGE_WITH_TETGEN)
    add_subdirectory(examples)
endif()

# ==============================================================================
# Installation & Package Config
# ==============================================================================
if(MESHARRANGE_INSTALL)
    include(GNUInstallDirs)
    include(CMakePackageConfigHelpers)

    install(TARGETS mesh_arrange
        EXPORT MeshArrangeTargets
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    )

    install(DIRECTORY include/mesh_arrange
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
    )

    install(EXPORT MeshArrangeTargets
        FILE MeshArrangeTargets.cmake
        NAMESPACE MeshArrange::
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/MeshArrange
    )

    export(EXPORT MeshArrangeTargets
        FILE "${CMAKE_CURRENT_BINARY_DIR}/MeshArrangeTargets.cmake"
        NAMESPACE MeshArrange::
    )

    configure_package_config_file(
        "${CMAKE_CURRENT_SOURCE_DIR}/cmake/MeshArrangeConfig.cmake.in"
        "${CMAKE_CURRENT_BINARY_DIR}/MeshArrangeConfig.cmake"
        INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/MeshArrange
    )

    write_basic_package_version_file(
        "${CMAKE_CURRENT_BINARY_DIR}/MeshArrangeConfigVersion.cmake"
        VERSION ${PROJECT_VERSION}
        COMPATIBILITY SameMajorVersion
    )

    install(FILES
        "${CMAKE_CURRENT_BINARY_DIR}/MeshArrangeConfig.cmake"
        "${CMAKE_CURRENT_BINARY_DIR}/MeshArrangeConfigVersion.cmake"
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/MeshArrange
    )
endif()
