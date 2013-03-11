# set up project
project("dune-geometry" C CXX)

# general stuff
cmake_minimum_required(VERSION 2.8.6)

#find dune-common and set the module path
find_package(dune-common)
list(APPEND CMAKE_MODULE_PATH ${dune-common_MODULE_PATH})
message("dune-common_MODULE_PATH=${dune-common_MODULE_PATH}")
#include the dune macros
include(DuneMacros)

# start a dune project with information from dune.module
dune_project()

add_subdirectory("m4")
add_subdirectory("dune")
add_subdirectory("doc")
add_subdirectory("lib")

# finalize the dune project, e.g. generating config.h etc.
finalize_dune_project(GENERATE_CONFIG_H_CMAKE)