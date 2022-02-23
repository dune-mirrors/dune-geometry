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
#include <dune/geometry/referenceelements.hh>
#include <dune/geometry/type.hh>

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

// Placeholder type for a trivial identity matrix without any functionality
struct IdentityMatrix {};

/// \brief Efficient implementation of a trivial local geometry for elements
template <class ctype, int mydim, int cdim = mydim>
struct DefaultLocalGeometry
{
  using LocalCoordinate = FieldVector<ctype,mydim>;
  using GlobalCoordinate = FieldVector<ctype,cdim>;

  template <class LocalCoordinate>
  decltype(auto) global (LocalCoordinate&& local) const
  {
    return std::forward<LocalCoordinate>(local);
  }

  template <class LocalCoordinate>
  IdentityMatrix jacobianTransposed (LocalCoordinate&& local) const
  {
    return {};
  }
};


/// \brief Geometry parametrized by a LocalFunction and a LocalGeometry
/**
 * \tparam Map  Mapping of element local coordinates to world coordinates, must be differentiable
 * \tparam LG   A local geometry type defining the local coordinates of the geometry and its
 *              transformation into local coordinates of the element the mapping the bound to.
 *              See \ref DefaultLocalGeometry
 * \tparam TraitsType  Parameters of the geometry, see \ref MappedGeometryTraits
 *
 * This class represents a geometry that is parametrized by the chained mapping of a local geometry
 * mapping and the given function. Thereby, the local geometry mapping related coordinates in
 * subentities of an element to element coordinates and the function then maps element coordinates
 * to the world coordinates representing the actual geometry.
 *
 * Typical choices for the local geometry are
 * - An identity mapping, if the geometry represents an element geometry.
 * - An intersection geometry-in-inside or geometry-in-outside mapping.
 * - A reference-element sub-entity geometry mapping.
 *
 * The function (or geometry mapping) is a local function in terms of the dune-functions concept. It
 * requires the following members:
 * - `bind(LocalContext)`: Initialize the mapping on a given local context, i.e., an element or
 *                         entity or intersection.
 * - `localContext()`: Return the element (or entity or intersection) the mapping is bound to.
 * - `operator()(LG::GlobalCoordinate)`: Evaluation in local-context coordinates.
 * - `derivative(Map)`: A free function returning a local function that represents the first
 *                      derivative of the mapping.
 *
 * The requirements on the traits are documented along with their default, `MappedGeometryTraits`.
 **/
template <class Map, class LG,
          class TraitsType = MappedGeometryTraits<typename LG::LocalCoordinate::value_type>>
class MappedGeometry
{
public:
  /// type of local coordinates
  using LocalCoordinate = typename LG::LocalCoordinate;

  /// type of global coordinates
  using GlobalCoordinate = std::result_of_t<Map(typename LG::GlobalCoordinate)>;

  /// coordinate type
  using ctype = typename LocalCoordinate::value_type;

  /// geometry dimension
  static const int mydimension = LocalCoordinate::size();

  /// coordinate dimension
  static const int coorddimension = GlobalCoordinate::size();

  /// type of volume
  using Volume = ctype;

  /// type of jacobian transposed
  using JacobianTransposed = FieldMatrix<ctype, mydimension, coorddimension>;

  /// type of jacobian inverse transposed
  using JacobianInverseTransposed = FieldMatrix<ctype, coorddimension, mydimension>;

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

  // type of coordinate transformation for subEntities to codim=0 entities
  using LocalGeometry = LG;

  // type of the mapping representation the geometry parametrization
  using Mapping = Map;

  // type of a mapping representing the derivative of `Map`
  using DerivativeMapping = std::decay_t<decltype(derivative(std::declval<Map>()))>;

public:
  /// \brief Constructor from mapping to parametrize the geometry
  /**
   *  \param[in]  refElement     reference element for the geometry
   *  \param[in]  mapping        a differentiable local function for the parametrization of the
   *                             geometry (stored by value)
   *  \param[in]  localGeometry  local geometry for local coordinate transformation
   **/
  template <class Map_, class LG_>
  MappedGeometry (const ReferenceElement& refElement, Map_&& mapping, LG_&& localGeometry)
    : refElement_(refElement)
    , mapping_(std::forward<Map_>(mapping))
    , localGeometry_(std::forward<LG_>(localGeometry))
  {}

  /// \brief Constructor, forwarding to the other constructor that take a reference-element
  /**
   *  \param[in]  gt             geometry type
   *  \param[in]  mapping        a differentiable local function for the parametrization of the
   *                             geometry (stored by value)
   *  \param[in]  localGeometry  local geometry for local coordinate transformation
   **/
  template <class Map_, class LG_>
  MappedGeometry (GeometryType gt, Map_&& mapping, LG_&& localGeometry)
    : MappedGeometry(ReferenceElements::general(gt),
                     std::forward<Map_>(mapping),
                     std::forward<LG_>(localGeometry))
  {}

  /// \brief Constructor, forwarding to the other constructors, with LocalGeometry is
  /// DefaultLocalGeometry.
  /**
   *  \param[in]  refGeo   geometry type or reference element
   *  \param[in]  mapping  a differentiable local function for the parametrization of the
   *                       geometry (stored by value)
   **/
  template <class RefGeo, class Map_, class LG_ = LG,
    std::enable_if_t<std::is_same_v<LG_, DefaultLocalGeometry<ctype,mydimension>>, bool> = true>
  MappedGeometry (RefGeo&& refGeo, Map_&& mapping)
    : MappedGeometry(std::forward<RefGeo>(refGeo),
                     std::forward<Map_>(mapping),
                     LG_{})
  {}


  /// \brief Is this mapping affine? Not in general, since we don't know anything about the mapping.
  bool affine () const
  {
    return false;
  }

  /// \brief Obtain the element type from the reference element
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
   *  Evaluate the local function in local coordinates
   *
   *  \param[in] local  local coordinates
   *  \returns          corresponding global coordinate
   **/
  GlobalCoordinate global (const LocalCoordinate& local) const
  {
    return mapping()(localGeometry().global(local));
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
    LocalCoordinate x = flatGeometry().local(globalCoord);
    LocalCoordinate dx;

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

      // break if tolerance is reached.
      if (dx.two_norm2() < Traits::tolerance())
        return x;
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

    // coordinate in the localContext of the mapping
    auto&& elementLocal = localGeometry().global(local);

    auto&& jLocal = localGeometry().jacobianTransposed(local);
    auto&& jGlobal = hostGeometry().jacobianTransposed(elementLocal);
    auto&& jacobian = (*dMapping_)(elementLocal);
    return AB(jLocal, ABt(jGlobal,jacobian));
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
  const ReferenceElement& refElement () const
  {
    return refElement_;
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
  const LocalGeometry& localGeometry () const
  {
    return localGeometry_;
  }

  // Construct a flat geometry from corner vertices
  using FlatGeometry = std::conditional_t<isFlatAffine,
    AffineGeometry<ctype, mydimension, coorddimension>,
    MultiLinearGeometry<ctype, mydimension, coorddimension>>;
  const FlatGeometry& flatGeometry () const
  {
    if (!flatGeometry_) {
      std::vector<GlobalCoordinate> corners;
      corners.reserve(refElement_.size(mydimension));
      for (int i = 0; i < refElement_.size(mydimension); ++i)
        corners.push_back(global(refElement_.position(i, mydimension)));

      flatGeometry_.emplace(refElement_, std::move(corners));
    }

    return *flatGeometry_;
  }

  // Return the geometry of the element the mapping is bound to
  using HostGeometry = std::decay_t<decltype(std::declval<Map>().localContext().geometry())>;
  const HostGeometry& hostGeometry () const
  {
    if (!hostGeometry_)
      hostGeometry_.emplace(mapping().localContext().geometry());

    return *hostGeometry_;
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

  // matrix-matrix multiplication A*B
  template <class K, int n, int m, int l>
  static FieldMatrix<K,n,m> AB (const FieldMatrix<K,n,l>& A, const FieldMatrix<K,l,m>& B)
  {
    FieldMatrix<K,n,m> AB;
    for (int i = 0; i < n; ++i)
      for (int j = 0; j < m; ++j) {
        AB[i][j] = 0;
        for (int k = 0; k < l; ++k)
          AB[i][j] += A[i][k] * B[k][j];
      }

    return AB;
  }

  // matrix-matrix multiplication A*B where A is a diagonal matrix
  template <class K, int n, int m>
  static FieldMatrix<K,n,m> AB (const DiagonalMatrix<K,n>& A, const FieldMatrix<K,n,m>& B)
  {
    FieldMatrix<K,n,m> AB;
    for (int i = 0; i < n; ++i)
      for (int j = 0; j < m; ++j) {
        AB[i][j] = A[i][i] * B[i][j];
      }

    return AB;
  }

  // multiplication of an identity-matrix with any real matrix
  template <class Mat>
  static const Mat& AB (const IdentityMatrix& A, const Mat& B)
  {
    return B;
  }

private:
  /// Reference of the geometry
  ReferenceElement refElement_;

  /// local parametrization of the element
  Mapping mapping_;

  /// transformation of local coordinates to element-local coordinates
  LocalGeometry localGeometry_;

  // some data optionally provided
  mutable std::optional<DerivativeMapping> dMapping_;
  mutable std::optional<FlatGeometry> flatGeometry_;
  mutable std::optional<HostGeometry> hostGeometry_;
};

// deduction guides
template <class RefGeo, class Map, class LG>
MappedGeometry (RefGeo, const Map&, const LG&)
  -> MappedGeometry<Map,LG>;

template <class I, class Map>
MappedGeometry (Geo::ReferenceElement<I>, const Map&)
  -> MappedGeometry<Map,DefaultLocalGeometry<typename I::ctype, I::dimension>,
                        MappedGeometryTraits<typename I::ctype> >;

// typedef for geometries on grid elements with local geometry is identity
template <class Mapping, class ctype, int dim>
using ElementMappedGeometry = MappedGeometry<Mapping,
  DefaultLocalGeometry<ctype,dim>, MappedGeometryTraits<ctype> >;

} // end namespace Dune

#endif // DUNE_GEOMETRY_MAPPEDGEOMETRY_HH
