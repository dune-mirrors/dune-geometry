// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <config.h>

#include <chrono>
#include <type_traits>
#include <vector>

#include <dune/geometry/chainedgeometry.hh>
#include <dune/geometry/mappedgeometry.hh>
#include <dune/geometry/multilineargeometry.hh>
#include <dune/geometry/referenceelements.hh>
#include <dune/geometry/type.hh>
#include <dune/geometry/test/checkgeometry.hh>

// A mapping modelling a local element-function
template <class ctype, int mydim, int cdim>
struct AffineMapping
{
  AffineMapping ()
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

  friend auto derivative (const AffineMapping& mapping)
  {
    return [&A=mapping.A](const Dune::FieldVector<ctype,mydim>& x) { return A; };
  }

private:
  Dune::FieldMatrix<ctype,cdim,mydim> A;
  Dune::FieldVector<ctype,cdim> b;
};

struct Timings
{
  double time1 = 0.0, time2 = 0.0, time3 = 0.0, time4 = 0.0, time5 = 0.0;
  std::size_t counter = 0;
};

template <class ctype, int cdim, Dune::GeometryType::Id id>
static bool testMappedGeometry (Timings& timings)
{
  bool pass = true;

  constexpr auto gt = Dune::GeometryType{id};
  auto refElem = Dune::referenceElement<ctype,gt.dim()>(gt);

  std::vector<Dune::FieldVector<ctype,gt.dim()>> corners(refElem.size(gt.dim()));
  for (int i = 0; i < refElem.size(gt.dim()); ++i)
    corners[i] = refElem.position(i,gt.dim());

  auto mapping = AffineMapping<ctype,gt.dim(),gt.dim()>{};

  // Transform using a MultiLinearGeometry
  auto start1 = std::chrono::high_resolution_clock::now();
  auto geo = Dune::MultiLinearGeometry<ctype,gt.dim(),gt.dim()>{refElem, corners};
  auto geometry1 = Dune::MappedGeometry{mapping, geo};
  pass &= checkGeometry(geometry1);
  auto end1 = std::chrono::high_resolution_clock::now();
  auto elapsed_seconds1 = std::chrono::duration_cast<std::chrono::duration<double>>(end1 - start1);
  timings.time1 += elapsed_seconds1.count();

  // Transform using the geometry of a reference element
  auto start2 = std::chrono::high_resolution_clock::now();
  auto refGeo = refElem.template geometry<0>(0);
  auto geometry2 = Dune::MappedGeometry{mapping, refGeo};
  pass &= checkGeometry(geometry2);
  auto end2 = std::chrono::high_resolution_clock::now();
  auto elapsed_seconds2 = std::chrono::duration_cast<std::chrono::duration<double>>(end2 - start2);
  timings.time2 += elapsed_seconds2.count();

  // Transform using the ReferenceElementGeometry wrapper
  auto start3 = std::chrono::high_resolution_clock::now();
  auto refElemWrapper = Dune::ReferenceElementGeometry{refElem};
  auto geometry3 = Dune::MappedGeometry{mapping, refElemWrapper};
  pass &= checkGeometry(geometry3);
  auto end3 = std::chrono::high_resolution_clock::now();
  auto elapsed_seconds3 = std::chrono::duration_cast<std::chrono::duration<double>>(end3 - start3);
  timings.time3 += elapsed_seconds3.count();

  // Transform using the LocalJacobianGeometry wrapper
  auto start4 = std::chrono::high_resolution_clock::now();
  auto jacobianTransform = Dune::LocalDerivativeGeometry{geometry1};
  auto geometry4 = Dune::MappedGeometry{mapping, jacobianTransform};
  pass &= checkGeometry(geometry4);
  auto end4 = std::chrono::high_resolution_clock::now();
  auto elapsed_seconds4 = std::chrono::duration_cast<std::chrono::duration<double>>(end4 - start4);
  timings.time4 += elapsed_seconds4.count();

  // Transform using the LocalJacobianGeometry wrapper
  auto start5 = std::chrono::high_resolution_clock::now();
  auto geometry5 = Dune::ChainedGeometry{geometry3, geo};
  pass &= checkGeometry(geometry4);
  auto end5 = std::chrono::high_resolution_clock::now();
  auto elapsed_seconds5 = std::chrono::duration_cast<std::chrono::duration<double>>(end5 - start5);
  timings.time5 += elapsed_seconds4.count();

  timings.counter++;

  return pass;
}


template <class ctype>
static bool testMappedGeometry (Timings& timings)
{
  bool pass = true;

  // pass &= testMappedGeometry<ctype, 0, 0, Dune::GeometryTypes::simplex(0)>();

  pass &= testMappedGeometry<ctype, 1, Dune::GeometryTypes::simplex(1)>(timings);
  pass &= testMappedGeometry<ctype, 2, Dune::GeometryTypes::simplex(1)>(timings);
  pass &= testMappedGeometry<ctype, 3, Dune::GeometryTypes::simplex(1)>(timings);
  pass &= testMappedGeometry<ctype, 4, Dune::GeometryTypes::simplex(1)>(timings);

  pass &= testMappedGeometry<ctype, 1, Dune::GeometryTypes::cube(1)>(timings);
  pass &= testMappedGeometry<ctype, 2, Dune::GeometryTypes::cube(1)>(timings);
  pass &= testMappedGeometry<ctype, 3, Dune::GeometryTypes::cube(1)>(timings);
  pass &= testMappedGeometry<ctype, 4, Dune::GeometryTypes::cube(1)>(timings);

  pass &= testMappedGeometry<ctype, 2, Dune::GeometryTypes::simplex(2)>(timings);
  pass &= testMappedGeometry<ctype, 3, Dune::GeometryTypes::simplex(2)>(timings);
  pass &= testMappedGeometry<ctype, 4, Dune::GeometryTypes::simplex(2)>(timings);

  pass &= testMappedGeometry<ctype, 2, Dune::GeometryTypes::cube(2)>(timings);
  pass &= testMappedGeometry<ctype, 3, Dune::GeometryTypes::cube(2)>(timings);
  pass &= testMappedGeometry<ctype, 4, Dune::GeometryTypes::cube(2)>(timings);

  pass &= testMappedGeometry<ctype, 3, Dune::GeometryTypes::simplex(3)>(timings);
  pass &= testMappedGeometry<ctype, 4, Dune::GeometryTypes::simplex(3)>(timings);

  pass &= testMappedGeometry<ctype, 3, Dune::GeometryTypes::cube(3)>(timings);
  pass &= testMappedGeometry<ctype, 4, Dune::GeometryTypes::cube(3)>(timings);

  pass &= testMappedGeometry<ctype, 4, Dune::GeometryTypes::simplex(4)>(timings);
  pass &= testMappedGeometry<ctype, 5, Dune::GeometryTypes::simplex(4)>(timings);

  pass &= testMappedGeometry<ctype, 4, Dune::GeometryTypes::cube(4)>(timings);
  pass &= testMappedGeometry<ctype, 5, Dune::GeometryTypes::cube(4)>(timings);

  pass &= testMappedGeometry<ctype, 3, Dune::GeometryTypes::pyramid>(timings);
  pass &= testMappedGeometry<ctype, 4, Dune::GeometryTypes::pyramid>(timings);

  pass &= testMappedGeometry<ctype, 3, Dune::GeometryTypes::prism>(timings);
  pass &= testMappedGeometry<ctype, 4, Dune::GeometryTypes::prism>(timings);

  return pass;
}

int main ( int argc, char **argv )
{
  bool pass = true;

  Timings timings{};

  std::cout << ">>> Checking ctype = double" << std::endl;
  for (std::size_t i = 0; i < 10; ++i)
    pass &= testMappedGeometry< double >(timings);
  //std::cout << ">>> Checking ctype = float" << std::endl;
  //pass &= testMappedGeometry< float >();

  std::cout << "timings:" << std::endl;
  std::cout << "time1: " << timings.time1/timings.counter << std::endl;
  std::cout << "time2: " << timings.time2/timings.counter << std::endl;
  std::cout << "time3: " << timings.time3/timings.counter << std::endl;
  std::cout << "time4: " << timings.time4/timings.counter << std::endl;
  std::cout << "time5: " << timings.time5/timings.counter << std::endl;

  return (pass ? 0 : 1);
}
