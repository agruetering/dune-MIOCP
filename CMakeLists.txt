# We require CMake version 3.1 to prevent issues
# with dune_enable_all_packages and older CMake versions.
cmake_minimum_required(VERSION 3.1)

option(CXX "enable C++ compilation" ON)
if(CXX)
    enable_language(CXX)
endif()

project(dune-MIOCP CXX)

if(NOT (dune-common_DIR OR dune-common_ROOT OR
      "${CMAKE_PREFIX_PATH}" MATCHES ".*dune-common.*"))
    string(REPLACE  ${CMAKE_PROJECT_NAME} dune-common dune-common_DIR
      ${PROJECT_BINARY_DIR})
endif()

#find dune-common and set the module path
find_package(dune-common REQUIRED)
list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake/modules"
  ${dune-common_MODULE_PATH})
  

find_package(GUROBI REQUIRED)

include_directories(${GUROBI_INCLUDE_DIRS})


option(DUNE_GRID_GRIDTYPE_SELECTOR "" ON)

#include the dune macros
include(DuneMacros)

# The directory where the example grid files used for the examples are located
set_property(DIRECTORY APPEND PROPERTY
  COMPILE_DEFINITIONS "PROJECT_SOURCE_DIR=\"${PROJECT_SOURCE_DIR}\"" "MIOCP_SOURCE_PATH=\"${PROJECT_SOURCE_DIR}/src/\"" "MIOCP_OUTPUT_PATH=\"${PROJECT_SOURCE_DIR}/output/\"")


# start a dune project with information from dune.module
dune_project()

dune_enable_all_packages()

add_subdirectory(src)
add_subdirectory(dune)
add_subdirectory(cmake/modules)



# finalize the dune project, e.g. generating config.h etc.
finalize_dune_project(GENERATE_CONFIG_H_CMAKE)


