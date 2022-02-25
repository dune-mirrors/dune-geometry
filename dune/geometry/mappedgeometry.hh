#ifndef DUNE_GEOMETRY_MAPPEDGEOMETRY_HH
#define DUNE_GEOMETRY_MAPPEDGEOMETRY_HH

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
#include <dune/geometry/referenceelementgeometry.hh>
#include <dune/geometry/referenceelements.hh>
#include <dune/geometry/type.hh>
#include <dune/geometry/utility/identitymatrix.hh>

namespace Dune {

/// \brief Default traits class for MappedGeometryTraits
/**
 *  The MappedGeometry allows tweaking
 *  some implementation details through a traits class.
 *
 *  This structure provides the default values.
 *
 *  \tparam  ct  coordinate type
 **/
template <class ct>
struct MappedGeometryTraits
{
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


/// \brief Geometry parametrized by a LocalFunction and a LocalGeometry
/**
 * \tparam Map  Mapping of element coordinates to global coordinates, must be differentiable
 *              w.r.t. local coordinates
 * \tparam Geo  A geometry type defining the local coordinates of the geometry and its
 *              transformation into the domain of the mapping.
 * \tparam TraitsType  Parameters of the geometry, see \ref MappedGeometryTraits
 *
 * This class represents a geometry that is parametrized by the chained mapping of a geometry
 * and the given function.
 *
 * The function (or geometry mapping) is a differentiable function in terms of the dune-functions
 * concept. It requires the following interface:
 * - `operator()(LG::GlobalCoordinate)`: Evaluation in local-context coordinates.
 * - `derivative(Map)`: A free function returning a mapping that represents the first derivative
 *                      of the mapping w.r.t. local coordinates in the reference element
 *
 * The requirements on the traits are documented along with their default, `MappedGeometryTraits`.
 **/
template <class Map, class Geo,
          class TraitsType = MappedGeometryTraits<typename Geo::ctype>>
class MappedGeometry
{
public:
  /// type of local coordinates
  using LocalCoordinate = typename Geo::LocalCoordinate;

  /// type of global coordinates
  using GlobalCoordinate = std::result_of_t<Map(typename Geo::GlobalCoordinate)>;

  /// coordinate type
  using ctype = typename Geo::ctype;

  /// geometry dimension
  static const int mydimension = LocalCoordinate::size();

  /// coordinate dimension
  static const int coorddimension = GlobalCoordinate::size();

  /// type of volume
  using Volume = decltype(Dune::power(std::declval<ctype>(),mydimension));

  /// type of jacobian transposed
  using JacobianTransposed = FieldMatrix<ctype, mydimension, coorddimension>;

  /// type of jacobian inverse transposed
  using JacobianInverseTransposed = FieldMatrix<ctype, coorddimension, mydimension>;

private:
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

  // type of the geometry that is wrapped
  using Geometry = Geo;

  // type of the mapping representation the geometry parametrization
  using Mapping = Map;

  // type of a mapping representing the derivative of `Map`
  using DerivativeMapping = std::decay_t<decltype(derivative(std::declval<Map>()))>;

public:
  /// \brief Constructor from mapping to parametrize the geometry
  /**
   *  \param[in]  mapping   A differentiable function for the parametrization of the geometry
   *  \param[in]  geometry  The geometry that is mapped
   **/
  template <class Geo_, class Map_,
    std::enable_if_t<Dune::IsInteroperable<Map, Map_>::value, int> = 0,
    std::enable_if_t<Dune::IsInteroperable<Geo, Geo_>::value, int> = 0>
  MappedGeometry (Map_&& mapping, Geo_&& geometry)
    : mapping_(std::forward<Map_>(mapping))
    , geometry_(std::forward<Geo_>(geometry))
  {}

  /// \brief Is this mapping affine? Not in general, since we don't know anything about the mapping.
  bool affine () const
  {
    return false;
  }

  /// \brief Obtain the element type from the reference element
  GeometryType type () const
  {
    return geometry_.type();
  }

  /// \brief Obtain number of corners of the corresponding reference element
  int corners () const
  {
    return geometry_.corners();
  }

  /// \brief Obtain coordinates of the i-th corner
  GlobalCoordinate corner (int i) const
  {
    assert( (i >= 0) && (i < corners()) );
    return mapping()(geometry_.corner(i));
  }

  /// \brief Obtain the centroid of the mapping's image
  GlobalCoordinate center () const
  {
    return mapping()(geometry_.center());
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
    return mapping()(geometry_.global(local));
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
      DUNE_THROW(Dune::Exception,
        "local coordinate can not be recovered from global coordinate " << globalCoord);

    return *localCoord;
  }

  /// \brief Evaluate the inverse coordinate mapping
  /**
   *  \param[in] globalCoord  global coordinate to map
   *  \return                 corresponding local coordinate or nothing, in case the local
   *                          coordinate could not be calculated. The return value is wrapped
   *                          in an optional.
   **/
  std::optional<LocalCoordinate> checkedLocal (const GlobalCoordinate& globalCoord) const
  {
    LocalCoordinate xOld = flatGeometry().local(globalCoord), x, dx;
    GlobalCoordinate dglobal = global(xOld) - globalCoord;
    ctype alpha = 1.0, resOld = dglobal.two_norm(), res = 0.0;

    for (int i = 0; i < Traits::maxIteration(); ++i)
    {
      // Newton's method: DF^n dx^n = F^n, x^{n+1} -= dx^n
      const bool invertible = MatrixHelper::xTRightInvA(jacobianTransposed(x), dglobal, dx);

      // break if jacobian is not invertible
      if (!invertible)
        return std::nullopt;

      // update x with correction
      for (int j = 0; j < Traits::maxIteration(); ++j) {
        x = xOld - alpha * dx;
        dglobal = global(x) - globalCoord;
        res = dglobal.two_norm();
        alpha /= 2.0;

        if (res < resOld)
          break;
      }

      // break if tolerance is reached.
      if (dx.two_norm2() < Traits::tolerance())
        return x;

      alpha *= 2.0;
    }

    if (dx.two_norm2() > Traits::tolerance())
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
    return [&]() {
      if constexpr ((mydimension == 1) && (coorddimension == 2))
        return normalDirection1D(local);
      else if constexpr ((mydimension == 2) && (coorddimension == 3))
        return normalDirection2D(local);
      else
        return GlobalCoordinate(0);
    }();
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
    if (!dMapping_)
      dMapping_.emplace(derivative(mapping()));

    auto&& jtLocal = geometry_.jacobianTransposed(local);
    auto&& jMapping = (*dMapping_)(geometry_.global(local));
    return ABt(jtLocal, jMapping);
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
  friend ReferenceElement referenceElement (const MappedGeometry& geometry)
  {
    return geometry.refElement();
  }

protected:
  // the internal stored reference element
  ReferenceElement refElement () const
  {
    return referenceElement(geometry_);
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

private:
  // internal reference to the stored mapping
  const Mapping& mapping () const
  {
    return mapping_;
  }

  // internal reference to the stored localgeometry
  const Geometry& geometry () const
  {
    return geometry_;
  }

  // Construct a flat geometry from corner vertices
  using FlatGeometry = std::conditional_t<isFlatAffine,
    AffineGeometry<ctype, mydimension, coorddimension>,
    MultiLinearGeometry<ctype, mydimension, coorddimension>>;
  const FlatGeometry& flatGeometry () const
  {
    if (!flatGeometry_) {
      std::vector<GlobalCoordinate> corners;
      corners.reserve(this->corners());
      for (int i = 0; i < this->corners(); ++i)
        corners.push_back(this->corner(i));

      flatGeometry_.emplace(referenceElement(geometry_), std::move(corners));
    }

    return *flatGeometry_;
  }

private:
  // matrix-matrix multiplication A*B^t
  template <class K, int n, int m, int l>
  static FieldMatrix<K,n,m> ABt (const FieldMatrix<K,n,l>& A, const FieldMatrix<K,m,l>& B)
  {
    FieldMatrix<K,n,m> ABt;
    for (int i = 0; i < n; ++i)
      for (int j = 0; j < m; ++j) {
        ABt[i][j] = 0;
        for (int k = 0; k < l; ++k)
          ABt[i][j] += A[i][k] * B[j][k];
      }

    return ABt;
  }

  // matrix-matrix multiplication A*B^t where A is a diagonal matrix
  template <class K, int n, int m>
  static FieldMatrix<K,n,m> ABt (const DiagonalMatrix<K,n>& A, const FieldMatrix<K,m,n>& B)
  {
    FieldMatrix<K,n,m> ABt;
    for (int i = 0; i < n; ++i)
      for (int j = 0; j < m; ++j) {
        ABt[i][j] = A[i][i] * B[j][i];
      }

    return ABt;
  }

  // matrix-matrix multiplication A*B^t where A is an identity matrix
  template <class K, int n, int m>
  static FieldMatrix<K,n,m> ABt (const IdentityMatrix& A, const FieldMatrix<K,m,n>& B)
  {
    FieldMatrix<K,n,m> Bt;
    for (int i = 0; i < n; ++i)
      for (int j = 0; j < m; ++j) {
        Bt[i][j] = B[j][i];
      }

    return Bt;
  }

private:
  /// Parametrization of the element
  Mapping mapping_;

  /// The geometry that is wrapped
  Geometry geometry_;

  // some data optionally provided
  mutable std::optional<DerivativeMapping> dMapping_;
  mutable std::optional<FlatGeometry> flatGeometry_;
};

// deduction guides
template <class Map, class Geo>
MappedGeometry (const Map&, const Geo&)
  -> MappedGeometry<Map,Geo>;

} // end namespace Dune

#endif // DUNE_GEOMETRY_MAPPEDGEOMETRY_HH
