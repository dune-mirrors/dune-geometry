add_subdirectory("quadraturerules")
add_subdirectory("refinement")
add_subdirectory("utility")
add_subdirectory("test")

install(FILES
  affinegeometry.hh
  axisalignedcubegeometry.hh
  dimension.hh
  generalvertexorder.hh
  multilineargeometry.hh
  quadraturerules.hh
  referenceelement.hh
  referenceelementimplementation.hh
  referenceelements.hh
  refinement.hh
  topologyfactory.hh
  type.hh
  typeindex.hh
  virtualrefinement.hh
  virtualrefinement.cc
  deprecated_topology.hh
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/geometry
)

# install the header as done for the auto-tools
install(FILES test/checkgeometry.hh
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/geometry/test)

dune_add_library(geometry OBJECT
  referenceelementimplementation.cc
  )
