dune_add_library(dunegeometry
  _DUNE_TARGET_OBJECTS:geometry_
  _DUNE_TARGET_OBJECTS:quadraturerules_
  ADD_LIBS dunecommon)

# create and install package config.hh header
dune_add_config_header(${PROJECT_SOURCE_DIR}/dune/geometry/config.hh.in ${PROJECT_BINARY_DIR}/dune/geometry/config.hh)
install(FILES ${PROJECT_BINARY_DIR}/dune/geometry/config.hh
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/geometry)

# force include of the generated config header
if(NOT DUNE_DISABLE_FORCED_INCLUDE)
  target_compile_options(dunegeometry PUBLIC
    $<BUILD_INTERFACE:-include${PROJECT_BINARY_DIR}/dune/geometry/config.hh>
    $<INSTALL_INTERFACE:-include$<INSTALL_PREFIX>/${CMAKE_INSTALL_INCLUDEDIR}/dune/geometry/config.hh>)
endif()
