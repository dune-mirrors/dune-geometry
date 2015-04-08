add_subdirectory("test" EXCLUDE_FROM_ALL)
install(FILES
  codimtable.hh
  conversion.hh
  geometrytraits.hh
  matrixhelper.hh
  referencedomain.hh
  referenceelements.hh
  subtopologies.hh
  topologytypes.hh
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/geometry/genericgeometry)

dune_library_add_sources(
  dunegeometry
  SOURCES
  referencedomain.cc
  subtopologies.cc
  )
