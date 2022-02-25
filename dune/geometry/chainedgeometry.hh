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
  using JacobianTransposed = FieldMatrix<ctype, mydimension, coorddimension>;

  /// type of jacobian inverse transposed
  using JacobianInverseTransposed = FieldMatrix<ctype, coorddimension, mydimension>;

private:
  /// \brief helper structure containing some matrix routines. See affinegeometry.hh
  using MatrixHelper = Impl::FieldMatrixHelper<ctype>;

  struct Traits
  {
    /// \brief tolerance to numerical algorithms
    static ctype tolerance () { return ctype(16) * std::numeric_limits<ctype>::epsilon(); }

    /// \brief maximal number of Newton iteration in `geometry.local(global)`
    static int maxIteration () { return 100; }
  };

public:
  /// \brief Constructor a chained geometry mapping
  /**
   *  \param[in]  outerGeo  Outer geometry `f`
   *  \param[in]  innerGeo  Inner geometry `g`
   **/
  ChainedGeometry (const OuterGeo& outerGeo, const InnerGeo& innerGeo)
    : outerGeo_(outerGeo)
    , innerGeo_(innerGeo)
  {}

  /// \brief Is this mapping affine? Not in general, since we don't know anything about the mapping.
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

  /// \brief Evaluate the inverse coordinate mapping
  /**
   *  \param[in] global  global coordinate to map
   *  \return            corresponding local coordinate or nothing, in case the local
   *                     coordinate could not be calculated. The return value is wrapped
   *                     in an optional.
   **/
  std::optional<LocalCoordinate> checkedLocal (const GlobalCoordinate& global) const
  {
    if (auto outerLocal = outerGeo_.checkedLocal(global); outerLocal.has_value())
      return innerGeo_.checkedLocal(*outerLocal);
    else
      return std::nullopt;
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
   * of the quadrature rules, until tolerance is reached.
   **/
  Volume volume () const
  {
    using std::abs;
    Volume vol0 = volume(QuadratureRules<ctype, mydimension>::rule(type(), 1));
    for (int p = 2; p < Traits::maxIteration(); ++p) {
      Volume vol1 = volume(QuadratureRules<ctype, mydimension>::rule(type(), p));
      if (abs(vol1 - vol0) < Traits::tolerance())
        return vol1;

      vol0 = vol1;
    }
    return vol0;
  }

  /// \brief Obtain the volume of the mapping's image by given quadrature rules
  Volume volume (const QuadratureRule<ctype, mydimension>& quadRule) const
  {
    Volume vol(0);
    for (const auto& qp : quadRule)
      vol += integrationElement(qp.position()) * qp.weight();
    return vol;
  }

  /// \brief Obtain the transposed of the Jacobian
  /**
   *  \param[in]  local  local coordinate to evaluate Jacobian in
   *  \returns           the matrix corresponding to the transposed of the Jacobian
   **/
  JacobianTransposed jacobianTransposed (const LocalCoordinate& local) const
  {
    auto&& jtInner = innerGeo_.jacobianTransposed(local);
    auto&& jtOuter = outerGeo_.jacobianTransposed(innerGeo_.global(local));
    return jtInner * jtOuter;
  }


  /// \brief obtain the transposed of the Jacobian's inverse
  /**
   *  The Jacobian's inverse is defined as a pseudo-inverse. If we denote
   *  the Jacobian by \f$J(x)\f$, the following condition holds:
   *  \f[ J^{-1}(x) J(x) = I. \f]
   **/
  JacobianInverseTransposed jacobianInverseTransposed (const LocalCoordinate& local) const
  {
    JacobianInverseTransposed out;
    MatrixHelper::rightInvA(jacobianTransposed(local), out);
    return out;
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
