// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <config.h>

#include <type_traits>
#include <vector>

#include <dune/geometry/mappedgeometry.hh>
#include <dune/geometry/multilineargeometry.hh>
#include <dune/geometry/referenceelements.hh>
#include <dune/geometry/type.hh>
#include <dune/geometry/test/checkgeometry.hh>

// A mapping modelling a local element-function
template <class ctype, int mydim, int cdim>
struct AffineMapping
{
  // Minimal base element with a geometry
  struct LocalContext
  {
    const Dune::MultiLinearGeometry<ctype,mydim,cdim>& geometry_;
    const Dune::MultiLinearGeometry<ctype,mydim,cdim>& geometry () const { return geometry_; }
  };

  AffineMapping (const Dune::MultiLinearGeometry<ctype,mydim,cdim>& geometry)
    : localContext_{geometry}
  {
    for (int i = 0; i < cdim; ++i)
      b[i] = i + 1;
    for (int i = 0; i < mydim; ++i)
      A[i][i] = i + 1;
  }

  Dune::FieldVector<ctype,cdim> operator() (const Dune::FieldVector<ctype,mydim>& x) const
  {
    Dune::FieldVector<ctype,cdim> result;
    A.mv(x, result);
    result += b;
    return result;
  }

  // the element, this mapping is bound to
  const LocalContext& localContext () const { return localContext_; }

  // bind the mapping to a local element
  void bind (const LocalContext&) {}

  friend auto derivative (const AffineMapping& mapping)
  {
    return [&A=mapping.A](const Dune::FieldVector<ctype,mydim>& x) { return A; };
  }

private:
  LocalContext localContext_;
  Dune::FieldMatrix<ctype,cdim,mydim> A;
  Dune::FieldVector<ctype,cdim> b;
};


template <class ctype, int cdim, Dune::GeometryType::Id id>
static bool testMappedGeometry ()
{
  bool pass = true;

  constexpr auto gt = Dune::GeometryType{id};
  auto refElem = Dune::referenceElement<ctype,gt.dim()>(gt);

  std::vector<Dune::FieldVector<ctype,gt.dim()>> corners(refElem.size(gt.dim()));
  for (int i = 0; i < refElem.size(gt.dim()); ++i)
    corners[i] = refElem.position(i,gt.dim());

  auto localContext = Dune::MultiLinearGeometry<ctype,gt.dim(),gt.dim()>{refElem, corners};
  auto mapping = AffineMapping{localContext};

  // construct a geometry using given parametrization coefficients
  auto geometry = Dune::MappedGeometry{refElem, mapping};
  pass &= checkGeometry(geometry);

  return pass;
}


template <class ctype>
static bool testMappedGeometry ()
{
  bool pass = true;

  // pass &= testMappedGeometry<ctype, 0, 0, Dune::GeometryTypes::simplex(0)>();

  pass &= testMappedGeometry<ctype, 1, Dune::GeometryTypes::simplex(1)>();
  pass &= testMappedGeometry<ctype, 2, Dune::GeometryTypes::simplex(1)>();
  pass &= testMappedGeometry<ctype, 3, Dune::GeometryTypes::simplex(1)>();
  pass &= testMappedGeometry<ctype, 4, Dune::GeometryTypes::simplex(1)>();

  pass &= testMappedGeometry<ctype, 1, Dune::GeometryTypes::cube(1)>();
  pass &= testMappedGeometry<ctype, 2, Dune::GeometryTypes::cube(1)>();
  pass &= testMappedGeometry<ctype, 3, Dune::GeometryTypes::cube(1)>();
  pass &= testMappedGeometry<ctype, 4, Dune::GeometryTypes::cube(1)>();

  pass &= testMappedGeometry<ctype, 2, Dune::GeometryTypes::simplex(2)>();
  pass &= testMappedGeometry<ctype, 3, Dune::GeometryTypes::simplex(2)>();
  pass &= testMappedGeometry<ctype, 4, Dune::GeometryTypes::simplex(2)>();

  pass &= testMappedGeometry<ctype, 2, Dune::GeometryTypes::cube(2)>();
  pass &= testMappedGeometry<ctype, 3, Dune::GeometryTypes::cube(2)>();
  pass &= testMappedGeometry<ctype, 4, Dune::GeometryTypes::cube(2)>();

  pass &= testMappedGeometry<ctype, 3, Dune::GeometryTypes::simplex(3)>();
  pass &= testMappedGeometry<ctype, 4, Dune::GeometryTypes::simplex(3)>();

  pass &= testMappedGeometry<ctype, 3, Dune::GeometryTypes::cube(3)>();
  pass &= testMappedGeometry<ctype, 4, Dune::GeometryTypes::cube(3)>();

  pass &= testMappedGeometry<ctype, 4, Dune::GeometryTypes::simplex(4)>();
  pass &= testMappedGeometry<ctype, 5, Dune::GeometryTypes::simplex(4)>();

  pass &= testMappedGeometry<ctype, 4, Dune::GeometryTypes::cube(4)>();
  pass &= testMappedGeometry<ctype, 5, Dune::GeometryTypes::cube(4)>();

  pass &= testMappedGeometry<ctype, 3, Dune::GeometryTypes::pyramid>();
  pass &= testMappedGeometry<ctype, 4, Dune::GeometryTypes::pyramid>();

  pass &= testMappedGeometry<ctype, 3, Dune::GeometryTypes::prism>();
  pass &= testMappedGeometry<ctype, 4, Dune::GeometryTypes::prism>();

  return pass;
}

int main ( int argc, char **argv )
{
  bool pass = true;

  std::cout << ">>> Checking ctype = double" << std::endl;
  pass &= testMappedGeometry< double >();
  //std::cout << ">>> Checking ctype = float" << std::endl;
  //pass &= testMappedGeometry< float >();

  return (pass ? 0 : 1);
}
