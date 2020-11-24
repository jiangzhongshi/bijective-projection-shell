cmake_minimum_required(VERSION 3.1)
project(prism)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)


set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
set(PRISM_EXTERNAL ${CMAKE_CURRENT_SOURCE_DIR}/external/)

include(PrismDependencies)

if(UNIX)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fPIC")
  if(NOT ${U_CMAKE_BUILD_TYPE} MATCHES DEBUG)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fvisibility=hidden -flto")
  endif()
endif()

add_library(prism_library
    src/correspondence.cpp
    src/prism/PrismCage.cpp
    src/prism/extraction.cpp
    src/prism/feature_utils.cpp
    src/prism/cage_utils.cpp
    src/prism/bevel_utils.cpp
    src/prism/geogram/geogram_utils.cpp
    src/prism/cgal/polyhedron_self_intersect.cpp
    src/prism/cgal/triangle_triangle_intersection.cpp
    src/prism/cgal/tetrahedron_tetrahedron_intersection.cpp
    src/prism/cgal/QP.cpp
    src/prism/geogram/AABB.cpp
    src/prism/phong/projection.cpp
    src/prism/phong/tangent_orientation.cpp
    src/prism/phong/query_correspondence.cpp
    src/prism/predicates/inside_prism_tetra.cpp
    src/prism/predicates/inside_octahedron.cpp
    src/prism/predicates/positive_prism_volume_12.cpp
    src/prism/predicates/triangle_triangle_intersection.cpp
    src/prism/predicates/tetrahedron_overlap.cpp
    src/prism/energy/prism_quality.cpp
    src/prism/energy/map_distortion.cpp
    src/prism/energy/smoother_pillar.cpp
    src/prism/local_operations/remesh_pass.cpp
    src/prism/local_operations/smooth_pass.cpp
    src/prism/local_operations/remesh_pass_collapse.cpp
    src/prism/local_operations/section_remesh.cpp
    src/prism/local_operations/mesh_coloring.cpp
    src/prism/local_operations/local_mesh_edit.cpp
    src/prism/local_operations/validity_checks.cpp
    src/prism/local_operations/retain_triangle_adjacency.cpp
    src/prism/spatial-hash/AABB_hash.cpp
    src/prism/spatial-hash/self_intersection.cpp
    src/prism/osqp/osqp_normal.cpp
    src/prism/cage_check.cpp
    src/prism/intersections.cpp
  )

target_compile_features(prism_library PUBLIC cxx_std_17)
target_link_libraries(prism_library PUBLIC spdlog igl::core igl::cgal osqp highfive geogram mitsuba_autodiff)
target_include_directories(prism_library PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src/)

option(PYPRISM "Activate Python Bindings" OFF)
if(PYPRISM)
  find_package(PythonInterp QUIET)
  prism_download_pybind11()
  add_subdirectory(${PRISM_EXTERNAL}/pybind11/ pybind11)
  add_subdirectory(src/seism/python)
endif()

option(PRISM_TESTS "Unit Tests" ON)
if(PRISM_TESTS)
  include(CTest)
  add_subdirectory(tests)
endif()

add_executable(breve_bin src/main.cpp
    src/construct_shell.cpp
    src/remesh_section_in_shell.cpp
    src/getRSS.c
)
target_link_libraries(breve_bin prism_library CLI11::CLI11 json)