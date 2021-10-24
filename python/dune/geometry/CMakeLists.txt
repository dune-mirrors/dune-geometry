# SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
# SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception

add_python_targets(geometry
  __init__
  _referenceelements
  quadpy
)

dune_add_pybind11_module(NAME _geometry)
set_property(TARGET _geometry PROPERTY LINK_LIBRARIES Dune::Geometry APPEND)
if(SKBUILD)
  install(TARGETS _geometry LIBRARY DESTINATION python/dune/geometry)
endif()

set( DIMS 0 1 2 3 )
foreach( DIM ${DIMS} )
  dune_add_pybind11_module(NAME _defaultreferenceelements_${DIM}
                           SOURCES _defaultreferenceelements.cc
                           COMPILE_DEFINITIONS PRIVATE "DIM=${DIM};MODNAME=_defaultreferenceelements_${DIM}" )
  set_property(TARGET _defaultreferenceelements_${DIM} PROPERTY LINK_LIBRARIES dunecommon dunegeometry APPEND)
  if(SKBUILD)
    install(TARGETS _defaultreferenceelements_${DIM} LIBRARY DESTINATION python/dune/geometry)
  endif()
endforeach()
