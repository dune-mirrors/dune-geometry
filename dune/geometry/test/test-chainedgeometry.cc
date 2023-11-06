// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#include <config.h>

#include <chrono>
#include <type_traits>
#include <vector>

#include <dune/geometry/chainedgeometry.hh>
#include <dune/geometry/referenceelements.hh>
#include <dune/geometry/type.hh>
#include <dune/geometry/test/checkgeometry.hh>

struct Timings
{
  double time1 = 0.0, time2 = 0.0, time3 = 0.0, time4 = 0.0;
  std::size_t counter = 0;
};

template <class ctype, int cdim, Dune::GeometryType::Id id>
static bool testChainedGeometry (Timings& timings)
{
  bool pass = true;

  constexpr auto gt = Dune::GeometryType{id};
  auto refElem = Dune::referenceElement<ctype,gt.dim()>(gt);

  // Chain two reference-element geometries
  auto start1 = std::chrono::high_resolution_clock::now();
  auto geo1 = refElem.template geometry<0>(0);
  auto geometry1 = Dune::ChainedGeometry{geo1, geo1};
  pass &= checkGeometry(geometry1);
  auto end1 = std::chrono::high_resolution_clock::now();
  auto elapsed_seconds1 = std::chrono::duration_cast<std::chrono::duration<double>>(end1 - start1);
  timings.time1 += elapsed_seconds1.count();

  // Chain two ReferenceElementGeometries
  auto start2 = std::chrono::high_resolution_clock::now();
  auto geometry2 = Dune::ChainedGeometry{geo1, geo1};
  pass &= checkGeometry(geometry2);
  auto end2 = std::chrono::high_resolution_clock::now();
  auto elapsed_seconds2 = std::chrono::duration_cast<std::chrono::duration<double>>(end2 - start2);
  timings.time2 += elapsed_seconds2.count();

  timings.counter++;

  return pass;
}


template <class ctype>
static bool testChainedGeometry (Timings& timings)
{
  bool pass = true;

  // pass &= testChainedGeometry<ctype, 0, 0, Dune::GeometryTypes::simplex(0)>();

  pass &= testChainedGeometry<ctype, 1, Dune::GeometryTypes::simplex(1)>(timings);
  pass &= testChainedGeometry<ctype, 2, Dune::GeometryTypes::simplex(1)>(timings);
  pass &= testChainedGeometry<ctype, 3, Dune::GeometryTypes::simplex(1)>(timings);
  pass &= testChainedGeometry<ctype, 4, Dune::GeometryTypes::simplex(1)>(timings);

  pass &= testChainedGeometry<ctype, 1, Dune::GeometryTypes::cube(1)>(timings);
  pass &= testChainedGeometry<ctype, 2, Dune::GeometryTypes::cube(1)>(timings);
  pass &= testChainedGeometry<ctype, 3, Dune::GeometryTypes::cube(1)>(timings);
  pass &= testChainedGeometry<ctype, 4, Dune::GeometryTypes::cube(1)>(timings);

  pass &= testChainedGeometry<ctype, 2, Dune::GeometryTypes::simplex(2)>(timings);
  pass &= testChainedGeometry<ctype, 3, Dune::GeometryTypes::simplex(2)>(timings);
  pass &= testChainedGeometry<ctype, 4, Dune::GeometryTypes::simplex(2)>(timings);

  pass &= testChainedGeometry<ctype, 2, Dune::GeometryTypes::cube(2)>(timings);
  pass &= testChainedGeometry<ctype, 3, Dune::GeometryTypes::cube(2)>(timings);
  pass &= testChainedGeometry<ctype, 4, Dune::GeometryTypes::cube(2)>(timings);

  pass &= testChainedGeometry<ctype, 3, Dune::GeometryTypes::simplex(3)>(timings);
  pass &= testChainedGeometry<ctype, 4, Dune::GeometryTypes::simplex(3)>(timings);

  pass &= testChainedGeometry<ctype, 3, Dune::GeometryTypes::cube(3)>(timings);
  pass &= testChainedGeometry<ctype, 4, Dune::GeometryTypes::cube(3)>(timings);

  pass &= testChainedGeometry<ctype, 4, Dune::GeometryTypes::simplex(4)>(timings);
  pass &= testChainedGeometry<ctype, 5, Dune::GeometryTypes::simplex(4)>(timings);

  pass &= testChainedGeometry<ctype, 4, Dune::GeometryTypes::cube(4)>(timings);
  pass &= testChainedGeometry<ctype, 5, Dune::GeometryTypes::cube(4)>(timings);

  pass &= testChainedGeometry<ctype, 3, Dune::GeometryTypes::pyramid>(timings);
  pass &= testChainedGeometry<ctype, 4, Dune::GeometryTypes::pyramid>(timings);

  pass &= testChainedGeometry<ctype, 3, Dune::GeometryTypes::prism>(timings);
  pass &= testChainedGeometry<ctype, 4, Dune::GeometryTypes::prism>(timings);

  return pass;
}

int main ( int argc, char **argv )
{
  bool pass = true;

  Timings timings{};

  std::cout << ">>> Checking ctype = double" << std::endl;
  for (std::size_t i = 0; i < 10; ++i)
    pass &= testChainedGeometry< double >(timings);

  // test single precision geometries
  // std::cout << ">>> Checking ctype = float" << std::endl;
  // pass &= testChainedGeometry< float >();

  std::cout << "timings:" << std::endl;
  std::cout << "time1: " << timings.time1/timings.counter << std::endl;
  std::cout << "time2: " << timings.time2/timings.counter << std::endl;

  return (pass ? 0 : 1);
}
