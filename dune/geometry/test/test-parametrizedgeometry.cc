// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <config.h>

#include <vector>

#include <dune/geometry/parametrizedgeometry.hh>
#include <dune/geometry/referenceelements.hh>
#include <dune/geometry/test/checkgeometry.hh>

namespace Dune { namespace Impl
{
  template <class D, class R, unsigned int dim>
  struct DefaultLocalBasisTraits
  {
    using DomainFieldType = D;
    using RangeFieldType = R;
    using DomainType = FieldVector<D,dim>;
    using RangeType = FieldVector<R,1>;
    using JacobianType = FieldMatrix<R,1,dim>;

    enum {
      dimDomain = dim,
      dimRange = 1
    };
  };

  /// Lagrange shape functions of order 1 on the reference simplex
  template <class D, class R, unsigned int dim>
  class P1LocalBasis
  {
  public:
    using Traits = DefaultLocalBasisTraits<D,R,dim>;

    /// Number of shape functions
    static constexpr unsigned int size () { return dim+1; }

    /// Evaluate all shape functions
    void evaluateFunction (const typename Traits::DomainType& x,
                           std::vector<typename Traits::RangeType>& out) const
    {
      out.resize(size());
      out[0] = 1.0;
      for (unsigned int i=0; i<dim; i++)
      {
        out[0]  -= x[i];
        out[i+1] = x[i];
      }
    }

    /// Evaluate Jacobian of all shape functions
    void evaluateJacobian (const typename Traits::DomainType& x,
                           std::vector<typename Traits::JacobianType>& out) const
    {
      out.resize(size());
      std::fill(out[0][0].begin(), out[0][0].end(), -1);

      for (unsigned int i=0; i<dim; i++)
        for (unsigned int j=0; j<dim; j++)
          out[i+1][0][j] = (i==j);
    }

    /// Polynomial order of the shape functions
    static constexpr unsigned int order () {  return 1; }
  };

  /// Lagrange shape functions of order 1 on the reference cube
  template <class D, class R, unsigned int dim>
  class Q1LocalBasis
  {
  public:
    using Traits = DefaultLocalBasisTraits<D,R,dim>;

    /// Number of shape functions
    static constexpr unsigned int size () { return power(2, dim); }

    /// Evaluate all shape functions
    void evaluateFunction (const typename Traits::DomainType& x,
                           std::vector<typename Traits::RangeType>& out) const
    {
      out.resize(size());
      for (size_t i=0; i<size(); i++)
      {
        out[i] = 1;

        for (unsigned int j=0; j<dim; j++)
          // if j-th bit of i is set multiply with x[j], else with 1-x[j]
          out[i] *= (i & (1<<j)) ? x[j] :  1-x[j];
      }
    }

    /// Evaluate Jacobian of all shape functions
    void evaluateJacobian (const typename Traits::DomainType& x,
                           std::vector<typename Traits::JacobianType>& out) const
    {
      out.resize(size());

      // Loop over all shape functions
      for (unsigned int i=0; i<size(); i++)
      {
        // Loop over all coordinate directions
        for (unsigned int j=0; j<dim; j++)
        {
          // Initialize: the overall expression is a product
          // if j-th bit of i is set to 1, else -11
          out[i][0][j] = (i & (1<<j)) ? 1 : -1;

          for (unsigned int l=0; l<dim; l++)
          {
            if (j!=l)
              // if l-th bit of i is set multiply with x[l], else with 1-x[l]
              out[i][0][j] *= (i & (1<<l)) ? x[l] :  1-x[l];
          }
        }
      }
    }

    /// Polynomial order of the shape functions
    static constexpr unsigned int order () { return 1; }
  };


  template <class LB>
  class P1LocalInterpolation
  {
  public:
    /// Evaluate a given function at the Lagrange nodes
    template <class F, class C>
    void interpolate (F f, std::vector<C>& out) const
    {
      constexpr auto dim = LB::Traits::dimDomain;
      out.resize(LB::size());

      // vertex 0
      typename LB::Traits::DomainType x;
      std::fill(x.begin(), x.end(), 0);
      out[0] = f(x);

      // remaining vertices
      for (int i=0; i<dim; i++) {
        for (int j=0; j<dim; j++)
          x[j] = (i==j);

        out[i+1] = f(x);
      }
    }
  };

  template <class LB>
  class Q1LocalInterpolation
  {
  public:
    /// Evaluate a given function at the Lagrange nodes
    template <class F, class C>
    void interpolate (F f, std::vector<C>& out) const
    {
      constexpr auto dim = LB::Traits::dimDomain;
      out.resize(LB::size());

      typename LB::Traits::DomainType x;
      for (unsigned int i=0; i<LB::size(); i++)
      {
        // Generate coordinate of the i-th corner of the reference cube
        for (int j=0; j<dim; j++)
          x[j] = (i & (1<<j)) ? 1.0 : 0.0;

        out[i] = f(x);
      }
    }
  };


  /// Wrapper for local basis and local interpolation
  template <class LB, template <class> class LI>
  class LocalFiniteElement
  {
  public:
    struct Traits
    {
      using LocalBasisType = LB;
      using LocalInterpolationType = LI<LB>;
    };

    const auto& localBasis () const { return basis_; }
    const auto& localInterpolation () const { return interpolation_; }

    /// The number of shape functions
    static constexpr std::size_t size () { return LB::size(); }

  private:
    LB basis_;
    LI<LB> interpolation_;
  };

  template <class D, class R, int d>
  using P1LocalFiniteElement = LocalFiniteElement<P1LocalBasis<D,R,d>, P1LocalInterpolation>;

  template <class D, class R, int d>
  using Q1LocalFiniteElement = LocalFiniteElement<Q1LocalBasis<D,R,d>, Q1LocalInterpolation>;

}} // end namespace Dune::Impl


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
  auto f = [](Dune::FieldVector<ctype,gt.dim()> const& x) {
    Dune::FieldVector<ctype,cdim> y;
    for (int i = 0; i < mydim; ++i)
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
