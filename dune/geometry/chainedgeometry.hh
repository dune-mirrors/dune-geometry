// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_GEOMETRY_CHAINEDGEOMETRY_HH
#define DUNE_GEOMETRY_CHAINEDGEOMETRY_HH

#include <cassert>
#include <limits>
#include <optional>
#include <type_traits>
#include <vector>

#include <dune/common/exceptions.hh>
#include <dune/common/diagonalmatrix.hh>
#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>
#include <dune/geometry/affinegeometry.hh>
#include <dune/geometry/multilineargeometry.hh>
#include <dune/geometry/quadraturerules.hh>
#include <dune/geometry/referenceelements.hh>
#include <dune/geometry/type.hh>
#include <dune/geometry/utility/convergenceoptions.hh>

namespace Dune {

/// \brief A chained geometry `f o g`
/**
 * \tparam OuterGeo  The geometry mapping `f`
 * \tparam InnerGeo  The geometry mapping `g`
 *
 **/
template <class OuterGeo, class InnerGeo>
class ChainedGeometry
{
public:
  /// type of local coordinates
  using LocalCoordinate = typename InnerGeo::LocalCoordinate;

  /// type of global coordinates
  using GlobalCoordinate = typename OuterGeo::GlobalCoordinate;

  /// coordinate type
  using ctype = typename OuterGeo::ctype;

  /// geometry dimension
  static const int mydimension = InnerGeo::mydimension;

  /// coordinate dimension
  static const int coorddimension = OuterGeo::coorddimension;

  /// type of volume
  using Volume = typename OuterGeo::Volume;

  /// type of jacobian transposed
  using Jacobian = FieldMatrix<ctype, coorddimension, mydimension>;

  /// type of jacobian transposed
  using JacobianTransposed = FieldMatrix<ctype, mydimension, coorddimension>;

  /// type of jacobian inverse transposed
  using JacobianInverse = FieldMatrix<ctype, mydimension, coorddimension>;

  /// type of jacobian inverse transposed
  using JacobianInverseTransposed = FieldMatrix<ctype, coorddimension, mydimension>;

private:
  /// \brief helper structure containing some matrix routines. See affinegeometry.hh
  using MatrixHelper = Impl::FieldMatrixHelper<ctype>;

public:
  /// \brief Constructor a chained geometry mapping `f o g`
  /**
   *  \param[in]  outerGeo  Outer geometry `f`
   *  \param[in]  innerGeo  Inner geometry `g`
   **/
  ChainedGeometry (const OuterGeo& outerGeo, const InnerGeo& innerGeo)
    : outerGeo_(outerGeo)
    , innerGeo_(innerGeo)
  {}

  /// \brief Is this mapping affine? Only if inner and outer geometry are affine.
  bool affine () const
  {
    return innerGeo_.affine() && outerGeo_.affine();
  }

  /// \brief Obtain the element type from the reference element
  GeometryType type () const
  {
    return innerGeo_.type();
  }

  /// \brief Obtain number of corners of the corresponding reference element
  int corners () const
  {
    return innerGeo_.corners();
  }

  /// \brief Obtain coordinates of the i-th corner
  GlobalCoordinate corner (int i) const
  {
    return outerGeo_.global(innerGeo_.corner(i));
  }

  /// \brief Obtain the centroid of the mapping's image
  GlobalCoordinate center () const
  {
    return outerGeo_.global(innerGeo_.center());
  }

  /// \brief Evaluate the coordinate mapping
  /**
   *  Evaluate the local function in local coordinates
   *
   *  \param[in] local  local coordinates
   *  \returns          corresponding global coordinate
   **/
  GlobalCoordinate global (const LocalCoordinate& local) const
  {
    return outerGeo_.global(innerGeo_.global(local));
  }

  /// \brief Evaluate the inverse coordinate mapping
  /**
   *  \param[in] globalCoord  global coordinate to map
   *  \return                 corresponding local coordinate
   *
   *  \throws in case of an error indicating that the local coordinate can not be obtained,
   *          an exception is thrown. See \ref checkedLocal for a variant that returns
   *          an optional instead.
   *
   *  \note For given global coordinate `y` the returned local coordinate `x` that minimizes
   *  the following function over the local coordinate space spanned by the reference element.
   *  \code
   *  (global( x ) - y).two_norm()
   *  \endcode
   **/
  LocalCoordinate local (const GlobalCoordinate& global) const
  {
    return innerGeo_.local(outerGeo_.local(global));
  }

  ///  \brief Obtain the integration element
  /**
   *  If the Jacobian of the mapping is denoted by $J(x)$, the integration
   *  integration element \f$\mu(x)\f$ is given by
   *  \f[ \mu(x) = \sqrt{|\det (J^T(x) J(x))|}.\f]
   *
   *  \param[in]  local  local coordinate to evaluate the integration element in
   *
   *  \returns the integration element \f$\mu(x)\f$.
   **/
  ctype integrationElement (const LocalCoordinate& local) const
  {
    return MatrixHelper::sqrtDetAAT(jacobianTransposed(local));
  }

  /// \brief Obtain the volume of the mapping's image
  /**
   * Calculates the volume of the entity by numerical integration. Since the polynomial order of the
   * Volume element is not known, iteratively compute numerical integrals with increasing order
   * of the quadrature rules, until a tolerance is reached or the maximum number of iterations
   * is exceeded. The convergence properties are given in the convergence options argument.
   **/
  Volume volume (Impl::ConvergenceOptions<ctype> opts = {}) const
  {
    Volume vol0 = volume(QuadratureRules<ctype, mydimension>::rule(type(), 1));
    if (affine())
      return vol0;

    using std::abs;
    for (int p = 2; p < opts.maxIt; ++p) {
      Volume vol1 = volume(QuadratureRules<ctype, mydimension>::rule(type(), p));
      if (abs(vol1 - vol0) < opts.absTol)
        return vol1;

      vol0 = vol1;
    }
    return vol0;
  }

  /// \brief Obtain the volume of the mapping's image by given quadrature rule
  Volume volume (const QuadratureRule<ctype, mydimension>& quadRule) const
  {
    Volume vol(0);
    for (const auto& qp : quadRule)
      vol += integrationElement(qp.position()) * qp.weight();
    return vol;
  }

  /// \brief Obtain the Jacobian of the geometry
  /**
   *  \param[in]  local  local coordinate to evaluate Jacobian in
   *  \returns           the matrix corresponding to the Jacobian
   **/
  Jacobian jacobian (const LocalCoordinate& local) const
  {
    auto&& jInner = innerGeo_.jacobian(local);
    auto&& jOuter = outerGeo_.jacobian(innerGeo_.global(local));
    return jOuter * jInner;
  }

  /// \brief Obtain the transposed of the Jacobian
  /**
   *  \param[in]  local  local coordinate to evaluate Jacobian in
   *  \returns           the matrix corresponding to the transposed of the Jacobian
   **/
  JacobianTransposed jacobianTransposed (const LocalCoordinate& local) const
  {
    return jacobian(local).transposed();
  }

  /// \brief obtain the Jacobian's inverse
  /**
   *  The Jacobian's inverse is defined as a pseudo-inverse. If we denote
   *  the Jacobian by \f$J(x)\f$, the following condition holds:
   *  \f[ J^{-1}(x) J(x) = I. \f]
   **/
  JacobianInverse jacobianInverse (const LocalCoordinate& local) const
  {
    JacobianInverse out;
    MatrixHelper::leftInvA(jacobian(local), out);
    return out;
  }

  /// \brief obtain the transposed of the Jacobian's inverse
  /**
   *  The Jacobian's inverse is defined as a pseudo-inverse. If we denote
   *  the Jacobian by \f$J(x)\f$, the following condition holds:
   *  \f[ J^{-1}(x) J(x) = I. \f]
   **/
  JacobianInverseTransposed jacobianInverseTransposed (const LocalCoordinate& local) const
  {
    return jacobianInverse(local).transposed();
  }

  /// \brief obtain the reference-element related to this geometry
  friend auto referenceElement (const ChainedGeometry& geometry)
  {
    return referenceElement(geometry.innerGeo_);
  }

private:
  /// The outer geometry mapping
  OuterGeo outerGeo_;

  /// The inner geometry mapping
  InnerGeo innerGeo_;
};

} // end namespace Dune

#endif // DUNE_GEOMETRY_CHAINEDGEOMETRY_HH
