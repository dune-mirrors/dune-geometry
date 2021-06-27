// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <config.h>

#include <type_traits>
#include <vector>

#include <dune/geometry/parametrizedgeometry.hh>
#include <dune/geometry/referenceelements.hh>
#include <dune/geometry/type.hh>
#include <dune/geometry/test/checkgeometry.hh>
#include <dune/geometry/test/localfiniteelements.hh>


template <class ctype, int cdim, Dune::GeometryType::Id id>
static bool testParametrizedGeometry ()
{
  bool pass = true;

  constexpr auto gt = Dune::GeometryType{id};
  auto refElem = Dune::referenceElement<ctype,gt.dim()>(gt);

  using LFE = std::conditional_t<
    gt.isSimplex(), Dune::Impl::P1LocalFiniteElement<ctype,ctype,gt.dim()>, std::conditional_t<
    gt.isCube(),    Dune::Impl::Q1LocalFiniteElement<ctype,ctype,gt.dim()>, void>
    >;
  auto lfe = LFE{};

  // mapping to generate coordinates from reference-element corners
  auto f = [&](Dune::FieldVector<ctype,gt.dim()> const& x) {
    Dune::FieldVector<ctype,cdim> y;
    for (std::size_t i = 0; i < gt.dim(); ++i)
      y[i] = x[i] + i;
    return y;
  };

  auto corners = std::vector<Dune::FieldVector<ctype,cdim>>(refElem.size(gt.dim()));
  for (int i = 0; i < refElem.size(gt.dim()); ++i)
    corners[i] = f(refElem.position(i, gt.dim()));

  using Geometry = Dune::ParametrizedGeometry<LFE,cdim>;

  // construct a geometry using given parametrization coefficients
  auto geometry = Geometry{refElem, lfe, corners};
  pass &= checkGeometry(geometry);

  // construct a geometry using local interpolation
  auto geometry2 = Geometry{refElem, lfe, f};
  pass &= checkGeometry(geometry2);

  return pass;
}


template <class ctype>
static bool testParametrizedGeometry ()
{
  bool pass = true;

  // pass &= testParametrizedGeometry<ctype, 0, 0, Dune::GeometryTypes::simplex(0)>();

  pass &= testParametrizedGeometry<ctype, 1, Dune::GeometryTypes::simplex(1)>();
  pass &= testParametrizedGeometry<ctype, 2, Dune::GeometryTypes::simplex(1)>();
  pass &= testParametrizedGeometry<ctype, 3, Dune::GeometryTypes::simplex(1)>();
  pass &= testParametrizedGeometry<ctype, 4, Dune::GeometryTypes::simplex(1)>();

  pass &= testParametrizedGeometry<ctype, 1, Dune::GeometryTypes::cube(1)>();
  pass &= testParametrizedGeometry<ctype, 2, Dune::GeometryTypes::cube(1)>();
  pass &= testParametrizedGeometry<ctype, 3, Dune::GeometryTypes::cube(1)>();
  pass &= testParametrizedGeometry<ctype, 4, Dune::GeometryTypes::cube(1)>();

  pass &= testParametrizedGeometry<ctype, 2, Dune::GeometryTypes::simplex(2)>();
  pass &= testParametrizedGeometry<ctype, 3, Dune::GeometryTypes::simplex(2)>();
  pass &= testParametrizedGeometry<ctype, 4, Dune::GeometryTypes::simplex(2)>();

  pass &= testParametrizedGeometry<ctype, 2, Dune::GeometryTypes::cube(2)>();
  pass &= testParametrizedGeometry<ctype, 3, Dune::GeometryTypes::cube(2)>();
  pass &= testParametrizedGeometry<ctype, 4, Dune::GeometryTypes::cube(2)>();

  pass &= testParametrizedGeometry<ctype, 3, Dune::GeometryTypes::simplex(3)>();
  pass &= testParametrizedGeometry<ctype, 4, Dune::GeometryTypes::simplex(3)>();

  pass &= testParametrizedGeometry<ctype, 3, Dune::GeometryTypes::cube(3)>();
  pass &= testParametrizedGeometry<ctype, 4, Dune::GeometryTypes::cube(3)>();

  // pass &= testParametrizedGeometry<ctype, 3, 3, Dune::GeometryTypes::pyramid>();
  // pass &= testParametrizedGeometry<ctype, 3, 4, Dune::GeometryTypes::pyramid>();

  // pass &= testParametrizedGeometry<ctype, 3, 3, Dune::GeometryTypes::prism>();
  // pass &= testParametrizedGeometry<ctype, 3, 4, Dune::GeometryTypes::prism>();

  pass &= testParametrizedGeometry<ctype, 4, Dune::GeometryTypes::simplex(4)>();
  pass &= testParametrizedGeometry<ctype, 5, Dune::GeometryTypes::simplex(4)>();

  pass &= testParametrizedGeometry<ctype, 4, Dune::GeometryTypes::cube(4)>();
  pass &= testParametrizedGeometry<ctype, 5, Dune::GeometryTypes::cube(4)>();

  return pass;
}

int main ( int argc, char **argv )
{
  bool pass = true;

  std::cout << ">>> Checking ctype = double" << std::endl;
  pass &= testParametrizedGeometry< double >();
  //std::cout << ">>> Checking ctype = float" << std::endl;
  //pass &= testParametrizedGeometry< float >();

  return (pass ? 0 : 1);
}
