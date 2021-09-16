add_python_targets(geometry
  __init__
  _referenceelements
  quadpy
)
dune_add_pybind11_module(NAME _geometry)
target_link_libraries(_geometry PUBLIC Dune::geometry)
if(SKBUILD)
  install(TARGETS _geometry LIBRARY DESTINATION python/dune/geometry)
endif()
