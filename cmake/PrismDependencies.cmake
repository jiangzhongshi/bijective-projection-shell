
include(PrismDownloadExternal)

if (NOT TARGET igl::core)
  prism_download_libigl()
  set(LIBIGL_EIGEN_VERSION 3.3.7)
  set(CGAL_DO_NOT_WARN_ABOUT_CMAKE_BUILD_TYPE TRUE)
  option(LIBIGL_USE_STATIC_LIBRARY "" OFF)
  option(LIBIGL_WITH_CGAL              "Use CGAL"           ON)
  option(LIBIGL_WITH_EMBREE              "Use embree"           ON)
  set(LIBIGL_INCLUDE_DIR ${PRISM_EXTERNAL}/libigl/include)
  find_package(LIBIGL REQUIRED)
endif()

if(NOT TARGET spdlog::spdlog)
  prism_download_spdlog()
  option(SPDLOG_BUILD_SHARED ON)
  add_subdirectory(${PRISM_EXTERNAL}/spdlog)
  set_property(TARGET spdlog PROPERTY POSITION_INDEPENDENT_CODE ON)
endif()

if (NOT TARGET json)
  prism_download_json()
  add_library(json INTERFACE)
  target_include_directories(json SYSTEM INTERFACE ${PRISM_EXTERNAL}/json/include)
endif()

if(NOT TARGET CLI11::CLI11)
    prism_download_cli11()
    add_subdirectory(${PRISM_EXTERNAL}/cli11)
endif()

if(NOT TARGET highfive)
  prism_download_HighFive()
  option(HIGHFIVE_USE_EIGEN ON)
  find_package(HDF5 REQUIRED)
  add_library(highfive INTERFACE)
  target_include_directories(highfive SYSTEM INTERFACE ${PRISM_EXTERNAL}/HighFive/include/ ${HDF5_INCLUDE_DIRS})
  target_link_libraries(highfive INTERFACE ${HDF5_LIBRARIES})
endif()

if(NOT TARGET geogram::geogram)
  set_property(GLOBAL PROPERTY ALLOW_DUPLICATE_CUSTOM_TARGETS 1) # geogram conflicts with embree
  prism_download_geogram()
  include(geogram)
endif()

if (NOT TARGET tbb::tbb)
  prism_download_tbb()
  add_subdirectory(${PRISM_EXTERNAL}/tbb)
endif()

if (NOT TARGET cvc3_rational)
  file(DOWNLOAD https://raw.githubusercontent.com/wildmeshing/fTetWild/master/src/external/Rational.h
  ${PRISM_EXTERNAL}/rational/Rational.h)
  add_library(cvc3_rational INTERFACE)
  target_include_directories(cvc3_rational INTERFACE ${PRISM_EXTERNAL}/rational/)
  target_link_directories(cvc3_rational INTERFACE ${GMP_LIBRARIES})
endif()

if (NOT TARGET mitsuba_autodiff)
  file(DOWNLOAD https://www.mitsuba-renderer.org/files/eigen/autodiff.h
  ${PRISM_EXTERNAL}/autodiff/autodiff_mitsuba.h)
  add_library(mitsuba_autodiff INTERFACE)
  target_include_directories(mitsuba_autodiff INTERFACE ${PRISM_EXTERNAL}/autodiff/)
endif()

if (NOT TARGET osqp)
  prism_download_project(osqp
  GIT_REPOSITORY https://github.com/oxfordcontrol/osqp.git
  GIT_TAG da403d4b41e86b7dc00237047ea4f00354d902ed
  # GIT_SHALLOW true
  # GIT_SUBMODULES qdldl
  )
  add_subdirectory(${PRISM_EXTERNAL}/osqp)
endif()