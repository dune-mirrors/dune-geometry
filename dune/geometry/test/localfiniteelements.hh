// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_GEOMETRY_TEST_LOCALFINITEELEMENT_HH
#define DUNE_GEOMETRY_TEST_LOCALFINITEELEMENT_HH

#include <algorithm>
#include <array>
#include <numeric>
#include <vector>

#include <dune/common/exceptions.hh>
#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>
#include <dune/common/math.hh>

namespace Dune::Impl {

template <class D, class R, unsigned int dim>
struct ScalarLocalBasisTraits
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
  using Traits = ScalarLocalBasisTraits<D,R,dim>;

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

  /** \brief Evaluate partial derivatives of any order of all shape functions */
  void partial(const std::array<unsigned int,dim>& order,
                const typename Traits::DomainType& in,
                std::vector<typename Traits::RangeType>& out) const
  {
    auto totalOrder = std::accumulate(order.begin(), order.end(), 0);
    out.resize(size());

    if (totalOrder == 0)
      evaluateFunction(in, out);
    else if (totalOrder == 1)
    {
      size_t direction = std::distance(order.begin(), std::find(order.begin(), order.end(), 1u));
      for (size_t i=0; i<size(); i++)
      {
        out[i][0] = (i == direction);
      }
    }
    else if (totalOrder == 2)
    {
      for (size_t i=0; i<size(); i++)
      {
        out[i][0] = 0;
      }
    }
    else
      DUNE_THROW(Dune::NotImplemented, "Partial derivative of order " << totalOrder << " is not implemented!");
  }

  /// Polynomial order of the shape functions
  static constexpr unsigned int order () {  return 1; }
};

/// Lagrange shape functions of order 1 on the reference cube
template <class D, class R, unsigned int dim>
class Q1LocalBasis
{
public:
  using Traits = ScalarLocalBasisTraits<D,R,dim>;

  /// Number of shape functions
  static constexpr unsigned int size () { return Dune::power(2, dim); }

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

  /** \brief Evaluate partial derivatives of any order of all shape functions */
  void partial(const std::array<unsigned int,dim>& order,
                const typename Traits::DomainType& in,
                std::vector<typename Traits::RangeType>& out) const
  {
    auto totalOrder = std::accumulate(order.begin(), order.end(), 0);
    out.resize(size());

    if (totalOrder == 0)
      evaluateFunction(in, out);
    else if (totalOrder == 1)
    {
      auto direction = std::distance(order.begin(), std::find(order.begin(), order.end(), 1));
      if (direction >= dim)
        DUNE_THROW(Dune::RangeError, "Direction of partial derivative not found!");

      // Loop over all shape functions
      for (std::size_t i = 0; i < size(); ++i)
      {
        // Initialize: the overall expression is a product
        // if j-th bit of i is set to 1, otherwise to -1
        out[i] = (i & (1<<direction)) ? 1 : -1;

        for (unsigned int j = 0; j < dim; ++j)
        {
          if (direction != j)
            // if j-th bit of i is set multiply with in[j], else with 1-in[j]
            out[i] *= (i & (1<<j)) ? in[j] :  1-in[j];
        }
      }
    }
    else if (totalOrder == 2)
    {
      for (size_t i=0; i<size(); i++)
      {
        // convert index i to multiindex
        std::array<unsigned int,dim> alpha(multiindex(i));

        // Initialize: the overall expression is a product
        out[i][0] = 1.0;

        // rest of the product
        for (std::size_t l=0; l<dim; l++)
        {
          switch (order[l])
          {
            case 0:
              out[i][0] *= p(alpha[l],in[l]);
              break;
            case 1:
              out[i][0] *= dp(alpha[l],in[l]);
              break;
            case 2:
              out[i][0] *= 0;
              break;
            default:
              DUNE_THROW(Dune::NotImplemented, "Desired derivative order is not implemented");
          }
        }
      }
    }
    else
      DUNE_THROW(Dune::NotImplemented, "Partial derivative of order " << totalOrder << " is not implemented!");
  }

  /// Polynomial order of the shape functions
  static constexpr unsigned int order () { return 1; }

private:
  // i-th Lagrange polynomial of degree k in one dimension
  static R p(unsigned int i, D x)
  {
    R result(1.0);
    for (unsigned int j=0; j<=1; j++)
      if (j!=i) result *= (x-j)/((int)i-(int)j);
    return result;
  }

  // derivative of ith Lagrange polynomial of degree k in one dimension
  static R dp(unsigned int i, D x)
  {
    R result(0.0);

    for (unsigned int j=0; j<=1; j++)
    {
      if (j!=i)
      {
        R prod( (1.0)/((int)i-(int)j) );
        for (unsigned int l=0; l<=1; l++)
          if (l!=i && l!=j)
            prod *= (x-l)/((int)i-(int)l);
        result += prod;
      }
    }
    return result;
  }

  // Return i as a d-digit number in the (k+1)-nary system
  static std::array<unsigned int,dim> multiindex (unsigned int i)
  {
    std::array<unsigned int,dim> alpha;
    for (unsigned int j=0; j<dim; j++)
    {
      alpha[j] = i % 2;
      i = i/2;
    }
    return alpha;
  }
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



template <class LFE, int cdim,
          class R = typename LFE::Traits::LocalBasisType::Traits::RangeFieldType>
class LocalFiniteElementFunction
{
  using LocalFiniteElement = LFE;
  using LocalBasis = typename LocalFiniteElement::Traits::LocalBasisType;
  using LocalBasisRange = typename LocalBasis::Traits::RangeType;
  using LocalBasisJacobian = typename LocalBasis::Traits::JacobianType;
  using Domain = typename LocalBasis::Traits::DomainType;
  using Range = FieldVector<R,cdim>;
  using Jacobian = FieldMatrix<R,cdim,LocalBasis::Traits::dimDomain>;

  static_assert(LocalBasis::Traits::dimRange == 1);

public:
  LocalFiniteElementFunction () = default;
  LocalFiniteElementFunction (const LocalFiniteElement& lfe, std::vector<Range> coefficients)
    : lfe_(lfe)
    , coefficients_(std::move(coefficients))
  {}

  Range operator() (const Domain& local) const
  {
    thread_local std::vector<LocalBasisRange> shapeValues;
    lfe_.localBasis().evaluateFunction(local, shapeValues);
    assert(shapeValues.size() == coefficients_.size());
    Range range(0);
    for (std::size_t i = 0; i < shapeValues.size(); ++i)
      range.axpy(shapeValues[i], coefficients_[i]);
    return range;
  }

  friend auto derivative (const LocalFiniteElementFunction& f)
  {
    return [&lfe=f.lfe_,coefficients=f.coefficients_](const Domain& local) -> Jacobian
    {
      thread_local std::vector<LocalBasisJacobian> shapeJacobians;
      lfe.localBasis().evaluateJacobian(local, shapeJacobians);
      assert(shapeJacobians.size() == coefficients.size());
      Jacobian jacobian(0);
      for (std::size_t i = 0; i < shapeJacobians.size(); ++i) {
        for (int j = 0; j < Jacobian::rows; ++j) {
          shapeJacobians[i].umtv(coefficients[i][j], jacobian[j]);
        }
      }
      return jacobian;
    };
  }

private:
  LocalFiniteElement lfe_{};
  std::vector<Range> coefficients_{};
};

} // end namespace Dune::Impl

#endif // DUNE_GEOMETRY_TEST_LOCALFINITEELEMENT_HH
