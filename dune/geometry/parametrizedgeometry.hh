#ifndef DUNE_GEOMETRY_PARAMETRIZEDGEOMETRY_HH
#define DUNE_GEOMETRY_PARAMETRIZEDGEOMETRY_HH

#include <cassert>
#include <functional>
#include <iterator>
#include <limits>
#include <optional>
#include <type_traits>
#include <vector>

#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>
#include <dune/common/math.hh>
#include <dune/common/typetraits.hh>
#include <dune/common/std/type_traits.hh>

#include <dune/geometry/affinegeometry.hh>
#include <dune/geometry/multilineargeometry.hh>
#include <dune/geometry/quadraturerules.hh>
#include <dune/geometry/referenceelements.hh>
#include <dune/geometry/type.hh>

namespace Dune {

// ParametrizedGeometryTraits
// -------------------------

/// \brief default traits class for ParametrizedGeometry
/**
 *  The ParametrizedGeometry allow tweaking
 *  some implementation details through a traits class.
 *
 *  This structure provides the default values.
 *
 *  \tparam  ct        coordinate type
 */
template <class ct>
struct ParametrizedGeometryTraits
{
  using ctype = ct;

  /// \brief helper structure containing some matrix routines. See affinegeometry.hh
  using MatrixHelper = Impl::FieldMatrixHelper<ct>;

  /// \brief tolerance to numerical algorithms
  static ct tolerance () { return ct(16) * std::numeric_limits<ct>::epsilon(); }

  /// \brief maximal number of Newton iteration in `geometry.local(global)`
  static int maxIteration () { return 100; }

  /// \brief Geometry is associated to just one GeometryType
  template <int dim>
  struct hasSingleGeometryType
  {
    static const bool v = false;
    static const unsigned int topologyId = ~0u; //< optionally, the topologyId of the single GeometryType
  };
};



// ParametrizedGeometry
// -------------------

/// \brief Curved geometry implementation based on local-basis function parametrization
/**
 *  Parametrization of the geometry by any localfunction interpolated into a local finite-element space.
 *
 *  \tparam  LFE         Type of a local finite-element
 *  \tparam  cdim        coordinate dimension
 *  \tparam  TraitsType  Parameters of the geometry, see \ref ParametrizedGeometryTraits
 *
 *  The requirements on the traits are documented along with their default,
 *  ParametrizedGeometryTraits.
 */
template <class LFE, int cdim,
          class TraitsType = ParametrizedGeometryTraits<typename LFE::Traits::LocalBasisType::Traits::DomainFieldType>>
class ParametrizedGeometry
{
  using LocalFiniteElement = LFE;
  using LocalBasis = typename LFE::Traits::LocalBasisType;
  using LocalBasisTraits = typename LocalBasis::Traits;

public:
  /// coordinate type
  using ctype = typename LocalBasisTraits::DomainFieldType;

  /// geometry dimension
  static const int mydimension = LocalBasisTraits::dimDomain;

  /// coordinate dimension
  static const int coorddimension = cdim;

  /// type of local coordinates
  using LocalCoordinate = FieldVector<ctype, mydimension>;

  /// type of global coordinates
  using GlobalCoordinate = FieldVector<ctype, coorddimension>;

  /// type of volume
  using Volume = decltype(power(std::declval<ctype>(),mydimension));

  /// type of jacobian transposed
  using JacobianTransposed = FieldMatrix<ctype, mydimension, coorddimension>;

  /// type of jacobian inverse transposed
  using JacobianInverseTransposed = FieldMatrix<ctype, coorddimension, mydimension>;

  /// type of the extended Weingarten map
  using NormalGradient = FieldMatrix<ctype, coorddimension, coorddimension>;

public:
  /// type of reference element
  using ReferenceElements = Dune::ReferenceElements<ctype, mydimension>;
  using ReferenceElement = typename ReferenceElements::ReferenceElement;

  /// Parametrization of the geometry
  using Traits = TraitsType;

protected:
  using MatrixHelper = typename Traits::MatrixHelper;
  static const bool hasSingleGeometryType
    = Traits::template hasSingleGeometryType<mydimension>::v;
  static const bool isFlatAffine
    = hasSingleGeometryType && ((Traits::template hasSingleGeometryType<mydimension>::topologyId) >> 1 == 0);

public:
  /// \brief Constructor from a vector of coefficients of the LocalBasis parametrizing
  /// the geometry.
  /**
   *  \param[in]  refElement  reference element for the geometry
   *  \param[in]  localFE     Local finite-element to use for the parametrization
   *  \param[in]  vertices    Coefficients of the local interpolation into the basis
   *
   *  \note The vertices are stored internally, so if possible move an external vertex storage
   *        to this constructor
   **/
  ParametrizedGeometry (const ReferenceElement& refElement, const LocalFiniteElement& localFE,
                        std::vector<GlobalCoordinate> vertices)
    : refElement_(refElement)
    , localFE_(localFE)
    , vertices_(std::move(vertices))
  {
    assert(localFE_.size() == vertices_.size());
  }

  /// \brief Constructor from a local parametrization function, mapping local to (curved)
  /// global coordinates
  /**
   *  \param[in]  refElement  reference element for the geometry
   *  \param[in]  localFE     Local finite-element to use for the parametrization
   *  \param[in]  param       parametrization function with signature GlobalCoordinate(LocalCoordinate)`
   **/
  template <class Parametrization,
    std::enable_if_t<Std::is_callable<Parametrization(LocalCoordinate), GlobalCoordinate>::value, bool> = true>
  ParametrizedGeometry (const ReferenceElement& refElement, const LocalFiniteElement& localFE,
                        Parametrization&& param)
    : refElement_(refElement)
    , localFE_(localFE)
  {
    const auto& localInterpolation = localFE_.localInterpolation();
    localInterpolation.interpolate(param, vertices_);
  }

  /// \brief Constructor, forwarding to the other constructors that take a reference-element
  /**
   *  \param[in]  gt       geometry type
   *  \param[in]  args...  arguments passed to the other constructors
   **/
  template <class... Args>
  ParametrizedGeometry (GeometryType gt, Args&&... args)
    : ParametrizedGeometry(ReferenceElements::general(gt), std::forward<Args>(args)...)
  {}

  /// \brief Copy constructor
  ParametrizedGeometry (const ParametrizedGeometry& that)
    : ParametrizedGeometry(that.refElement_, that.localFE_, that.vertices_)
  {}

  /// \brief Move constructor
  ParametrizedGeometry (ParametrizedGeometry&& that)
    : ParametrizedGeometry(that.refElement_, std::move(that.localFE_), std::move(that.vertices_))
  {}

  /// \brief Copy assignment operator
  ParametrizedGeometry& operator= (const ParametrizedGeometry& that)
  {
    assert(refElement_ == that.refElement_);
    vertices_ = that.vertices_;
    return *this;
  }

  /// \brief Move assignment operator
  ParametrizedGeometry& operator= (ParametrizedGeometry&& that)
  {
    assert(refElement_ == that.refElement_);
    vertices_ = std::move(that.vertices_);
    return *this;
  }

  /// \brief Obtain the polynomial order of the parametrization
  int order () const
  {
    return localBasis().order();
  }

  /// \brief Is this mapping affine? This is only true for flat affine geometries.
  bool affine () const
  {
    if (!affine_)
      affine_ = (order() == 1 && (isFlatAffine || type().isSimplex() || flatGeometry().affine() ));
    return *affine_;
  }

  /// \brief Obtain the name of the reference element
  GeometryType type () const
  {
    return refElement_.type();
  }

  /// \brief Obtain number of corners of the corresponding reference element
  int corners () const
  {
    return refElement_.size(mydimension);
  }

  /// \brief Obtain coordinates of the i-th corner
  GlobalCoordinate corner (int i) const
  {
    assert( (i >= 0) && (i < corners()) );
    return global(refElement_.position(i, mydimension));
  }

  /// \brief Obtain the centroid of the mapping's image
  GlobalCoordinate center () const
  {
    return global(refElement_.position(0, 0));
  }

  /// \brief Evaluate the coordinate mapping
  /**
   *  Implements a linear combination of local basis functions scaled by
   *  the vertices as coefficients.
   *
   *  \f[ global = \sum_i v_i \psi_i(local) \f]
   *
   *  \param[in] local  local coordinate to map
   *  \returns          corresponding global coordinate
   **/
  GlobalCoordinate global (const LocalCoordinate& local) const
  {
    thread_local std::vector<typename LocalBasisTraits::RangeType> shapeValues;
    localBasis().evaluateFunction(local, shapeValues);
    assert(shapeValues.size() == vertices_.size());

    GlobalCoordinate out(0);
    for (std::size_t i = 0; i < shapeValues.size(); ++i)
      out.axpy(shapeValues[i], vertices_[i]);

    return out;
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
  LocalCoordinate local (const GlobalCoordinate& globalCoord) const
  {
    auto localCoord = checkedLocal(globalCoord);
    if (!localCoord)
      DUNE_THROW(Exception, "Local coordinate cannot be recovered from given global coordinate " << globalCoord);

    return *localCoord;
  }

  /// \brief Evaluate the inverse coordinate mapping
  /**
   *  \param[in] globalCoord  global coordinate to map
   *  \return                 optional wrapping the corresponding local coordinate
   *
   *  See \ref local() for some details.
   *
   *  The evaluation of local coordinates may fail if the jacobian is not invertible, or
   *  the Newton method to calculate the local coordinate fails to converge. Either the
   *  number of iteration or the tolerance in the \ref Traits class could be modified to
   *  control the convergence of the Newton method.
   **/
  std::optional<LocalCoordinate> checkedLocal (const GlobalCoordinate& globalCoord) const
  {
    const ctype tolerance = Traits::tolerance();
    LocalCoordinate x = flatGeometry().local(globalCoord);

    LocalCoordinate dx;
    const bool affineMapping = affine();

    for (int i = 0; i < Traits::maxIteration(); ++i)
    {
      // Newton's method: DF^n dx^n = F^n, x^{n+1} -= dx^n
      const GlobalCoordinate dglobal = global(x) - globalCoord;
      const bool invertible = MatrixHelper::xTRightInvA(jacobianTransposed(x), dglobal, dx);

      // break if jacobian is not invertible
      if (!invertible)
        return std::nullopt;

      // update x with correction
      x -= dx;

      // for affine mappings only one iteration is needed
      if (affineMapping)
        return x;

      // break if tolerance is reached.
      if (dx.two_norm2() < tolerance)
        return x;
    }

    if (dx.two_norm2() > tolerance)
      return std::nullopt;

    return x;
  }

  /// \brief Construct a normal vector of the curved element evaluated at
  /// a given local coordinate
  /**
   * \note Implemented for codim=1 entities only, i.e. edges in 2D and faces in 3D
   **/
  GlobalCoordinate normal (const LocalCoordinate& local) const
  {
    GlobalCoordinate n = normalDirection(local);
    return n / n.two_norm();
  }

  /// \brief Construct a normal direction (not normalized) of the curved element
  /// evaluated at a given local coordinate
  /**
   * \note Implemented for codim=1 entities only, i.e. edges in 2D and faces in 3D
   **/
  GlobalCoordinate normalDirection (const LocalCoordinate& local) const
  {
    assert(coorddimension == mydimension+1);
    return [&] {
      if constexpr ((mydimension == 1) && (coorddimension == 2))
        return normalDirection1D(local);
      else if constexpr ((mydimension == 2) && (coorddimension == 3))
        return normalDirection2D(local);
      else
        return GlobalCoordinate(0);
    }();
  }

  /// \brief Construct the surface gradient (extended Weingarten map) of the normal vector field
  /**
   * First, interpolate the normal vector field into a local Lagrange basis, then take the
   * derivative of this interpolated field, normalized it and project it into the tangential
   * plane.
   *
   * \param local   The local coordinate where to evaluate the normal-vector gradient
   **/
  NormalGradient normalGradient (const LocalCoordinate& local) const
  {
    return normalGradientImpl(local, jacobianInverseTransposed(local));
  }

  /// \brief Construct the surface gradient (extended Weingarten map) of the normal vector field
  /**
   * See \ref normalGradient() but with additional parameter.
   *
   * \param jiT   Evaluation of the JacobianInverseTransposed at the local coordinate `local`.
   *              This can be passed, if already computed elsewhere.
   **/
  NormalGradient normalGradientImpl (const LocalCoordinate& local,
                                     const JacobianInverseTransposed& jiT) const
  {
    if (nCoefficients_.empty()) {
      // create local discrete function of normal vectors by interpolation of the geometry normal
      localFE_.localInterpolation().interpolate(
        [&](const LocalCoordinate& l) { return this->normalDirection(l); }, nCoefficients_);
    }

    // Interpolated normal vector evaluated at local coordinate
    localFE_.localBasis().evaluateFunction(local, nShapeValues_);
    GlobalCoordinate nh(0);
    for (std::size_t j = 0; j < nShapeValues_.size(); ++j)
      nh.axpy(nShapeValues_[j], nCoefficients_[j]);
    auto nh_nrm = nh.two_norm();
    nh /= nh_nrm;

    // P = I - n x n
    NormalGradient Ph;
    for (int r = 0; r < coorddimension; ++r)
      for (int s = 0; s < coorddimension; ++s)
        Ph[r][s] = (r == s ? 1 : 0) - nh[r]*nh[s];

    // Compute the shape function gradients on the real element
    localFE_.localBasis().evaluateJacobian(local, nShapeGradients_);
    nGradients_.resize(nShapeGradients_.size());
    for (std::size_t j = 0; j < nGradients_.size(); ++j)
      jiT.mv(nShapeGradients_[j][0], nGradients_[j]);

    // Normal gradient evaluated at local coordinate
    NormalGradient H(0);
    for (std::size_t j = 0; j < nGradients_.size(); ++j)
      for (int r = 0; r < coorddimension; ++r)
        for (int s = 0; s < coorddimension; ++s)
          H[r][s] += nGradients_[j][s] * nCoefficients_[j][r];
    H /= nh_nrm;
    H.leftmultiply(Ph);
    H.rightmultiply(Ph);

    return H;
  }

  ///  \brief Obtain the integration element
  /**
   *  If the Jacobian of the mapping is denoted by \f$J(x)\f$, the integration
   *  element \f$\mu(x)\f$ is given by
   *
   *  \f[ \mu(x) = \sqrt{|\det (J^T(x) J(x))|}.\f]
   *
   *  \param[in]  local  local coordinate to evaluate the integration element in
   *  \returns           the integration element \f$\mu(x)\f$.
   **/
  ctype integrationElement (const LocalCoordinate& local) const
  {
    return MatrixHelper::sqrtDetAAT(jacobianTransposed(local));
  }

  /// \brief Obtain the volume of the mapping's image
  /**
   * Calculates the volume of the entity by numerical integration. Since the
   * polynomial order of the Volume element is not known, iteratively compute
   * numerical integrals with increasing order of the quadrature rules, until
   * tolerance is reached.
   **/
  Volume volume () const
  {
    using std::abs;
    Volume vol0 = volume(QuadratureRules<ctype, mydimension>::rule(type(), 1));
    for (int p = 2; p < 10; ++p) {
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
    thread_local std::vector<typename LocalBasisTraits::JacobianType> shapeJacobians;
    localBasis().evaluateJacobian(local, shapeJacobians);
    assert(shapeJacobians.size() == vertices_.size());

    JacobianTransposed out(0);
    for (std::size_t i = 0; i < shapeJacobians.size(); ++i)
      outerProductAccumulate(shapeJacobians[i], vertices_[i], out);

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
    JacobianInverseTransposed out;
    MatrixHelper::rightInvA(jacobianTransposed(local), out);
    return out;
  }

  /// \brief Obtain the reference-element related to this geometry
  friend ReferenceElement referenceElement (const ParametrizedGeometry& geometry)
  {
    return geometry.refElement();
  }

  /// \brief Obtain the coefficients of the parametrization
  const std::vector<GlobalCoordinate>& coefficients () const
  {
    return vertices_;
  }

protected:
  // the internal stored reference element
  const ReferenceElement& refElement () const
  {
    return refElement_;
  }

  // the local basis of the stored local finite-element
  const LocalBasis& localBasis () const
  {
    return localFE_.localBasis();
  }

  // normal vector to an edge line-element
  GlobalCoordinate normalDirection1D (const LocalCoordinate& local) const
  {
    auto J = jacobianTransposed(local);
    return GlobalCoordinate{
       J[0][1],
      -J[0][0]};
  }

  // normal vector to a triangle or quad face element
  GlobalCoordinate normalDirection2D (const LocalCoordinate& local) const
  {
    auto J = jacobianTransposed(local);
    return GlobalCoordinate{
      J[0][1] * J[1][2] - J[0][2] * J[1][1],
      J[0][2] * J[1][0] - J[0][0] * J[1][2],
      J[0][0] * J[1][1] - J[0][1] * J[1][0]};
  }

public:
  // type of a flat geometry build of the corner vertices
  using FlatGeometry = std::conditional_t<isFlatAffine,
    AffineGeometry<ctype, mydimension, coorddimension>,
    MultiLinearGeometry<ctype, mydimension, coorddimension>>;

  // construct a flat geometry from the corner vertices
  const FlatGeometry& flatGeometry () const
  {
    if (!flatGeometry_) {
      std::vector<GlobalCoordinate> corners;
      corners.reserve(refElement_.size(mydimension));
      for (int i = 0; i < refElement_.size(mydimension); ++i)
        corners.push_back(global(refElement_.position(i, mydimension)));

      flatGeometry_ = FlatGeometry{refElement_, corners};
    }

    return *flatGeometry_;
  }

private:
  // Let a and b be two column vectors then res += a^T*b
  template <class T, int n, int m>
  static void outerProductAccumulate (const FieldVector<T,n>& a, const FieldVector<T,m>& b,
                                      FieldMatrix<T,n,m>& res)
  {
    for (int i = 0; i < n; ++i)
      for (int j = 0; j < m; ++j)
        res[i][j] += a[i] * b[j];
  }

  // Let a and b be two column vectors then res += a^T*b
  template <class T, int n, int m>
  static void outerProductAccumulate (const FieldMatrix<T,n,1>& a, const FieldVector<T,m>& b,
                                      FieldMatrix<T,n,m>& res)
  {
    for (int i = 0; i < n; ++i)
      for (int j = 0; j < m; ++j)
        res[i][j] += a[i][0] * b[j];
  }

  // Let a be a row vector and b be a column vector then res += a*b
  template <class T, int n, int m,
    std::enable_if_t<(n > 1), int> = 0>
  static void outerProductAccumulate (const FieldMatrix<T,1,n>& a, const FieldVector<T,m>& b,
                                      FieldMatrix<T,n,m>& res)
  {
    outerProductAccumulate(a[0], b, res);
  }

private:
  /// Reference of the geometry
  ReferenceElement refElement_;

  /// A local finite-element
  LocalFiniteElement localFE_;

  /// The (lagrange) coefficients of the interpolating geometry
  std::vector<GlobalCoordinate> vertices_;

  // some data optionally provided
  mutable std::optional<bool> affine_;
  mutable std::optional<FlatGeometry> flatGeometry_;
  mutable std::vector<GlobalCoordinate> nCoefficients_;
  mutable std::vector<GlobalCoordinate> nGradients_;
  mutable std::vector<typename LocalBasisTraits::RangeType> nShapeValues_;
  mutable std::vector<typename LocalBasisTraits::JacobianType> nShapeGradients_;
};

namespace Impl {

// extract the LocalCoordinate type from a LocalFiniteElement
template <class LFE>
using LocalCoordinate_t
  = FieldVector<typename LFE::Traits::LocalBasisType::Traits::DomainFieldType,
                LFE::Traits::LocalBasisType::Traits::dimDomain>;

} // end namespace Impl


// deduction guides
template <class I, class LFE, class GlobalCoordinate>
ParametrizedGeometry (Geo::ReferenceElement<I>, const LFE&, std::vector<GlobalCoordinate>)
  -> ParametrizedGeometry<LFE, GlobalCoordinate::dimension>;

template <class I, class LFE, class F,
          class Range = std::result_of_t<F(Impl::LocalCoordinate_t<LFE>)>>
ParametrizedGeometry (Geo::ReferenceElement<I>, const LFE&, const F&)
  -> ParametrizedGeometry<LFE, Range::dimension>;

template <class LFE, class GlobalCoordinate>
ParametrizedGeometry (GeometryType, const LFE& localFE, std::vector<GlobalCoordinate>)
  -> ParametrizedGeometry<LFE, GlobalCoordinate::dimension>;

template <class LFE, class F,
          class Range = std::result_of_t<F(Impl::LocalCoordinate_t<LFE>)>>
ParametrizedGeometry (GeometryType, const LFE&, const F&)
  -> ParametrizedGeometry<LFE, Range::dimension>;

} // namespace Dune

#endif // DUNE_GEOMETRY_PARAMETRIZEDGEOMETRY_HH
