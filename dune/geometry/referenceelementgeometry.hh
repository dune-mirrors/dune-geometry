#ifndef DUNE_GEOMETRY_REFERENCEELEMENTGEOMETRY_HH
#define DUNE_GEOMETRY_REFERENCEELEMENTGEOMETRY_HH

#include <type_traits>

#include <dune/geometry/referenceelements.hh>
#include <dune/geometry/utility/identitymatrix.hh>

namespace Dune {

/// \brief Wrap a reference Element into a Geometry
template <class RefElem>
class ReferenceElementGeometry
    : public RefElem::template Codim<0>::Geometry
{
  using Super = typename RefElem::template Codim<0>::Geometry;

public:
  using LocalCoordinate = typename RefElem::Coordinate;
  using GlobalCoordinate = typename RefElem::Coordinate;
  using JacobianTransposed = IdentityMatrix;
  using JacobianInverseTransposed = IdentityMatrix;

public:
  ReferenceElementGeometry (const RefElem& refElem)
    : Super{refElem.template geometry<0>(0)}
  {}

  /// \brief Evaluate the inverse coordinate mapping, this is an identity mapping
  const LocalCoordinate& local (const GlobalCoordinate& global) const
  {
    return global;
  }

  /// \brief Evaluate the coordinate mapping, this is an identity mapping
  const GlobalCoordinate& global (const LocalCoordinate& local) const
  {
    return local;
  }

  /// \brief Obtain the transposed of the Jacobian, this is an identity matrix.
  JacobianTransposed jacobianTransposed (const LocalCoordinate& local) const
  {
    return {};
  }

  /// \brief obtain the transposed of the Jacobian's inverse, this is an identity matrix.
  JacobianInverseTransposed jacobianInverseTransposed (const LocalCoordinate& local) const
  {
    return {};
  }

private:
  RefElem refElem_;
};


/// \brief Wrapped Geometry that only transforms the derivatives
template <class Geometry>
class LocalDerivativeGeometry
    : public ReferenceElementGeometry<
        typename Dune::ReferenceElements<typename Geometry::ctype, Geometry::mydimension>::ReferenceElement >
{
  using ReferenceElements = Dune::ReferenceElements<typename Geometry::ctype, Geometry::mydimension>;
  using ReferenceElement = typename ReferenceElements::ReferenceElement;
  using Super = ReferenceElementGeometry<ReferenceElement>;

public:
  using LocalCoordinate = typename Geometry::LocalCoordinate;
  using JacobianTransposed = typename Geometry::JacobianTransposed;
  using JacobianInverseTransposed = typename Geometry::JacobianInverseTransposed;

public:
  LocalDerivativeGeometry (const Geometry& geometry)
    : Super{referenceElement(geometry)}
    , geometry_{geometry}
  {}

  /// \brief Obtain the transposed of the Jacobian
  JacobianInverseTransposed jacobianTransposed (const LocalCoordinate& local) const
  {
    return geometry_.jacobianTransposed(local);
  }

  /// \brief obtain the transposed of the Jacobian's inverse
  JacobianTransposed jacobianInverseTransposed (const LocalCoordinate& local) const
  {
    return geometry_.jacobianInverseTransposed(local);
  }

private:
  Geometry geometry_;
};

} // end namespace Dune

#endif // DUNE_GEOMETRY_REFERENCEELEMENTGEOMETRY_HH
