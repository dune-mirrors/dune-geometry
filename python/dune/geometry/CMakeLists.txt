add_python_targets(geometry
  __init__
  _referenceelements
  quadpy
)
dune_add_pybind11_module(NAME _geometry)
set_property(TARGET _geometry PROPERTY LINK_LIBRARIES dunecommon dunegeometry APPEND)
install(TARGETS _geometry LIBRARY DESTINATION python/dune/geometry)

set(PRECOMPILE
  _refelem0
  _refelem1
  _refelem2
  _refelem3
)

foreach(obj ${PRECOMPILE})
  dune_add_pybind11_module(NAME ${obj})
  set_property(TARGET ${obj} PROPERTY LINK_LIBRARIES dunecommon dunegeometry APPEND)
  install(TARGETS ${obj} LIBRARY DESTINATION python/dune/geometry)
endforeach()
