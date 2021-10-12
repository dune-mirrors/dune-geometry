// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <config.h>

#include <dune/python/geometry/referenceelements.hh>

#include <dune/python/common/typeregistry.hh>
#include <dune/python/geometry/type.hh>
#include <dune/python/pybind11/pybind11.h>

PYBIND11_MODULE( _refelem1, module )
{
  // register reference element (dim = 1)
  auto refElem1 = Dune::Python::insertClass< Dune::Geo::ReferenceElement<Dune::Geo::ReferenceElementImplementation<double,1>> >( module,
    "ReferenceElement1",
    Dune::Python::GenerateTypeName("Dune::Geo::ReferenceElement<Dune::Geo::ReferenceElementImplementation<double,1> >"),
    Dune::Python::IncludeFiles{"dune/python/geometry/referenceelements.hh"}
  ).first;
  Dune::Python::registerReferenceElements( module, refElem1 );
}
