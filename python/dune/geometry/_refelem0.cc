// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <config.h>

#include <dune/python/geometry/referenceelements.hh>

#include <dune/python/common/typeregistry.hh>
#include <dune/python/geometry/type.hh>
#include <dune/python/pybind11/pybind11.h>

PYBIND11_MODULE( _refelem0, module )
{
  // register reference element (dim = 0)
  auto refElem0 = Dune::Python::insertClass< Dune::Geo::ReferenceElement<Dune::Geo::ReferenceElementImplementation<double,0>> >( module,
    "ReferenceElement0",
    Dune::Python::GenerateTypeName("Dune::Geo::ReferenceElement<Dune::Geo::ReferenceElementImplementation<double,0> >"),
    Dune::Python::IncludeFiles{"dune/python/geometry/referenceelements.hh"}
  ).first;
  Dune::Python::registerReferenceElements( module, refElem0 );
}
