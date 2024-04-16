// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#include <config.h>

#include <dune/common/test/testsuite.hh>
#include <dune/geometry/type.hh>

#include <iostream>

template<Dune::GeometryType::Id gtid>
struct Foo
{
  static constexpr Dune::GeometryType gt = gtid;
  static unsigned int apply()
  {
    //return Foo<Dune::GeometryTypes::prismaticExtension(gt)>::gt.id(); //does not work for gcc < 10.2
    return Foo<Dune::GeometryTypes::prismaticExtension(gt).toId()>::gt.id();
  }
};

int main(int /* argc */, char** /* argv */)
{
  Dune::TestSuite suite;

  // make sure we can correctly roundtrip between GeometryType
  // and its Id in constexpr context
  constexpr Dune::GeometryType gt2a = Dune::GeometryTypes::triangle;

  Foo<gt2a> foo2;

  constexpr Dune::GeometryType gt2b = foo2.gt;

  static_assert(gt2a == gt2b, "The two geometry types have to compare equal");
  Foo<Dune::GeometryTypes::prismaticExtension(gt2b)> foo3;

  constexpr Dune::GeometryType gt3 = foo3.gt;
  static_assert(gt3 == Dune::GeometryTypes::prism, "The two geometry types have to compare equal");

  suite.check(foo2.apply() == foo3.gt.id());

  {
    using namespace Dune;
    GeometryType vertex = GeometryTypes::vertex;
    GeometryType line = GeometryTypes::line;
    GeometryType tri = GeometryTypes::simplex(2);
    GeometryType quad = GeometryTypes::cube(2);
    GeometryType tet = GeometryTypes::simplex(3);
    GeometryType hex = GeometryTypes::cube(3);
    GeometryType prism = GeometryTypes::prism;
    GeometryType pyramid = GeometryTypes::pyramid;

    suite.check(GeometryTypes::prismaticProduct(vertex,line) == line, "p(vertex,line)");
    suite.check(GeometryTypes::prismaticProduct(line,vertex) == line, "p(line,vertex)");
    suite.check(GeometryTypes::prismaticProduct(line,line) == quad, "p(line,line)");
    suite.check(GeometryTypes::prismaticProduct(tri,line) == prism, "p(tri,line)");
    // suite.check(GeometryTypes::prismaticProduct(line,tri) == prism, "p(line,tri)"); // ERROR: hex
    suite.check(GeometryTypes::prismaticProduct(quad,line) == hex, "p(quad,line)");
    suite.check(GeometryTypes::prismaticProduct(line,quad) == hex, "p(line,quad)");
    suite.check(GeometryTypes::prismaticProduct(hex,line) == GeometryTypes::cube(4));
    suite.check(GeometryTypes::prismaticProduct(quad,quad) == GeometryTypes::cube(4));

    suite.check(GeometryTypes::conicalProduct(vertex,line) == line);
    suite.check(GeometryTypes::conicalProduct(line,vertex) == line);
    suite.check(GeometryTypes::conicalProduct(line,line) == tri);
    suite.check(GeometryTypes::conicalProduct(tri,line) == tet);
    suite.check(GeometryTypes::conicalProduct(line,tri) == tet);
    suite.check(GeometryTypes::conicalProduct(quad,line) == pyramid);
    // suite.check(GeometryTypes::conicalProduct(line,quad) == pyramid); // ERROR: tet
    suite.check(GeometryTypes::conicalProduct(tet,line) == GeometryTypes::simplex(4));
    suite.check(GeometryTypes::conicalProduct(tri,tri) == GeometryTypes::simplex(4));
  }

  return suite.exit();

}
