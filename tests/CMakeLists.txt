cmake_minimum_required(VERSION 3.1)
project(prism_tests)

list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/../cmake)


### Download unit test framework
if (NOT TARGET doctest)
  prism_download_doctest()
  add_library(doctest INTERFACE)
  target_include_directories(doctest SYSTEM INTERFACE ${PRISM_EXTERNAL}/doctest/doctest/)
endif()
# Add catch2
# Create test executable
add_executable(prism_tests test_main.cpp)
target_sources(prism_tests PRIVATE 
                cgal_AABB.cpp
                cgal_qp.cpp
                spatial_hash.cpp
                predicates.cpp
                phong.cpp
                numerical_self_intersection.cpp
                feature_parse.cpp
                )

target_link_libraries(prism_tests PUBLIC prism_library doctest)

target_compile_features(prism_tests PUBLIC cxx_std_17)
target_include_directories(prism_tests PUBLIC ${CMAKE_CURRENT_LIST_DIR} ${PROJECT_SOURCE_DIR}/src/prism)
# Register tests
# set(PARSE_CATCH_TESTS_ADD_TO_CONFIGURE_DEPENDS ON)
# include(Catch)
# catch_discover_tests(prism_tests
                    # WORKING_DIRECTORY ${PROJECT_SOURCE_DIR})