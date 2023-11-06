// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_GEOMETRY_UTILITY_DEFAULTGEOMETRYTRAITS_HH
#define DUNE_GEOMETRY_UTILITY_DEFAULTGEOMETRYTRAITS_HH

#include <limits>

#include <dune/geometry/utility/fmatrixhelper.hh>

namespace Dune
{

/** \file
 *  \brief An implementation of the Geometry interface for affine geometries
 *  \author Martin Nolte
 */

  /** \brief default traits class for geometries
   *
   *  Some geometries allow tweaking some implementation details through
   *  a traits class. This type provides a default implementation.
   *
   *  \tparam  ct  coordinate type
   */
  template< class ct >
  struct DefaultGeometryTraits
  {
    /** \brief helper structure containing some matrix routines
     *
     *  This helper allows exchanging the matrix inversion algorithms.
     *  It must provide the following static methods:
     *  \code
     *  template< int m, int n >
     *  static ctype sqrtDetAAT ( const FieldMatrix< ctype, m, n > &A );
     *
     *  template< int m, int n >
     *  static ctype rightInvA ( const FieldMatrix< ctype, m, n > &A,
     *                           FieldMatrix< ctype, n, m > &ret );
     *
     *  template< int m, int n >
     *  static void xTRightInvA ( const FieldMatrix< ctype, m, n > &A,
     *                            const FieldVector< ctype, n > &x,
     *                            FieldVector< ctype, m > &y );
     *  \endcode
     */
    typedef Impl::FieldMatrixHelper< ct > MatrixHelper;

    /** \brief tolerance to numerical algorithms */
    static ct tolerance () { return ct( 16 ) * std::numeric_limits< ct >::epsilon(); }

    /** \brief maximal number of Newton iteration in `geometry.local(global)` */
    static int maxIteration () { return 100; }

    /** \brief will there be only one geometry type for a dimension?
     *
     *  If there is only a single geometry type for a certain dimension,
     *  <em>hasSingleGeometryType::v</em> can be set to true.
     *  Supporting only one geometry type might yield a gain in performance.
     *
     *  If <em>hasSingleGeometryType::v</em> is set to true, an additional
     *  parameter <em>topologyId</em> is required.
     *  Here's an example:
     *  \code
     *  static const unsigned int topologyId = GeometryTypes::simplex(dim).id();
     *  \endcode
     */
    template< int dim >
    struct hasSingleGeometryType
    {
      static const bool v = false;
      static const unsigned int topologyId = ~0u;
    };

    /** \brief flag whether the geometry is affine */
    static bool affine () { return false; }
  };

} // namespace Dune

#endif // DUNE_GEOMETRY_UTILITY_DEFAULTGEOMETRYTRAITS_HH