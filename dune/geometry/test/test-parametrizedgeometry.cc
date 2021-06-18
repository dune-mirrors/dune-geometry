// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <config.h>

#include <vector>

#include <dune/geometry/parametrizedgeometry.hh>
#include <dune/geometry/referenceelements.hh>

#include <dune/geometry/test/checkgeometry.hh>

namespace Dune { namespace Impl
{
  /** \brief Lagrange shape functions of order 1 on the reference simplex
   *
   * \tparam D Type to represent the field in the domain
   * \tparam R Type to represent the field in the range
   * \tparam dim Dimension of the domain simplex
   */
  template <class D, class R, unsigned int dim>
  class P1LocalBasis
  {
  public:
    struct Traits
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

    /// Number of shape functions
    static constexpr unsigned int size () { return dim+1; }

    //! \brief Evaluate all shape functions
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

    /** \brief Evaluate Jacobian of all shape functions
     *
     * \param x Point in the reference simplex where to evaluation the Jacobians
     * \param[out] out The Jacobians of all shape functions at the point x
     */
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

  /** \brief Evaluate the degrees of freedom of a Lagrange basis
   *
   * \tparam LB  The corresponding set of shape functions
   */
  template <class LB>
  class P1LocalInterpolation
  {
  public:
    /** \brief Evaluate a given function at the Lagrange nodes
     *
     * \tparam F Type of function to evaluate
     * \tparam C Type used for the values of the function
     * \param[in] ff Function to evaluate
     * \param[out] out Array of function values
     */
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

  /** \brief Lagrange finite element for simplices with order 1
   *
   * \tparam D type used for domain coordinates
   * \tparam R type used for function values
   * \tparam d dimension of the reference element
   */
  template <class D, class R, int d>
  class P1LocalFiniteElement
  {
    using LB = P1LocalBasis<D,R,d>;
    using LI = P1LocalInterpolation<LB>;

  public:
    struct Traits
    {
      using LocalBasisType = LB;
      using LocalInterpolationType = LI;
    };

    const auto& localBasis () const { return basis_; }
    const auto& localInterpolation () const { return interpolation_; }

    /// The number of shape functions
    static constexpr std::size_t size () { return LB::size(); }

    /// The reference element that the local finite element is defined on
    static constexpr GeometryType type () { return GeometryTypes::simplex(d); }

  private:
    LB basis_;
    LI interpolation_;
  };

}} // end namespace Dune::Impl


template< class ctype, int mydim, int cdim >
static bool testParametrizedGeometry ( Dune::GeometryType gt )
{
  bool pass = true;

  auto refElem = Dune::referenceElement<double,mydim>(gt);

  using LFE = Dune::Impl::P1LocalFiniteElement<double,double,mydim>;
  auto lfe = LFE{};

  // mapping to generate coordinates from reference-element corners
  auto f = [](Dune::FieldVector<double,mydim> const& x) {
    Dune::FieldVector<double,cdim> y;
    for (int i = 0; i < mydim; ++i)
      y[i] = x[i] + i;
    return y;
  };

  auto corners = std::vector<Dune::FieldVector<double,cdim>>(refElem.size(mydim));
  for (int i = 0; i < refElem.size(mydim); ++i)
    corners[i] = f(refElem.position(i, mydim));

  using Geometry = Dune::ParametrizedGeometry<LFE,cdim>;

  // construct a geometry using given parametrization coefficients
  auto geometry = Geometry{refElem, lfe, corners};
  pass &= checkGeometry(geometry);

  // construct a geometry using local interpolation
  auto geometry2 = Geometry{refElem, lfe, f};
  pass &= checkGeometry(geometry2);

  return pass;
}


template< class ctype >
static bool testParametrizedGeometry ()
{
  bool pass = true;

  // pass &= testParametrizedGeometry< ctype, 0, 0 >( Dune::GeometryTypes::simplex(0) );

  pass &= testParametrizedGeometry< ctype, 1, 1 >( Dune::GeometryTypes::simplex(1) );
  pass &= testParametrizedGeometry< ctype, 1, 2 >( Dune::GeometryTypes::simplex(1) );
  pass &= testParametrizedGeometry< ctype, 1, 3 >( Dune::GeometryTypes::simplex(1) );
  pass &= testParametrizedGeometry< ctype, 1, 4 >( Dune::GeometryTypes::simplex(1) );

  // pass &= testParametrizedGeometry< ctype, 1, 3 >( Dune::GeometryTypes::cube(1) );
  // pass &= testParametrizedGeometry< ctype, 1, 1 >( Dune::GeometryTypes::cube(1) );
  // pass &= testParametrizedGeometry< ctype, 1, 2 >( Dune::GeometryTypes::cube(1) );
  // pass &= testParametrizedGeometry< ctype, 1, 4 >( Dune::GeometryTypes::cube(1) );

  pass &= testParametrizedGeometry< ctype, 2, 2 >( Dune::GeometryTypes::simplex(2) );
  pass &= testParametrizedGeometry< ctype, 2, 3 >( Dune::GeometryTypes::simplex(2) );
  pass &= testParametrizedGeometry< ctype, 2, 4 >( Dune::GeometryTypes::simplex(2) );

  // pass &= testParametrizedGeometry< ctype, 2, 2 >( Dune::GeometryTypes::cube(2) );
  // pass &= testParametrizedGeometry< ctype, 2, 3 >( Dune::GeometryTypes::cube(2) );
  // pass &= testParametrizedGeometry< ctype, 2, 4 >( Dune::GeometryTypes::cube(2) );

  pass &= testParametrizedGeometry< ctype, 3, 3 >( Dune::GeometryTypes::simplex(3) );
  pass &= testParametrizedGeometry< ctype, 3, 4 >( Dune::GeometryTypes::simplex(3) );

  /** \bug These tests currently fail. */
  //   pass &= testParametrizedGeometry< ctype, 3, 3 >( Dune::GeometryTypes::pyramid );
  //   pass &= testParametrizedGeometry< ctype, 3, 4 >( Dune::GeometryTypes::pyramid );

  // pass &= testParametrizedGeometry< ctype, 3, 3 >( Dune::GeometryTypes::prism );
  // pass &= testParametrizedGeometry< ctype, 3, 4 >( Dune::GeometryTypes::prism );

  /** \bug These tests currently fail. */
  //   pass &= testParametrizedGeometry< ctype, 3, 3 >( Dune::GeometryTypes::cube(3) );
  //   pass &= testParametrizedGeometry< ctype, 3, 4 >( Dune::GeometryTypes::cube(3) );

  pass &= testParametrizedGeometry< ctype, 4, 4 >( Dune::GeometryTypes::simplex(4) );
  pass &= testParametrizedGeometry< ctype, 4, 5 >( Dune::GeometryTypes::simplex(4) );

  /** \bug These tests currently fail. */
  //   pass &= testParametrizedGeometry< ctype, 4, 4 >( Dune::GeometryTypes::cube(4) );
  //   pass &= testParametrizedGeometry< ctype, 4, 5 >( Dune::GeometryTypes::cube(4) );

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
