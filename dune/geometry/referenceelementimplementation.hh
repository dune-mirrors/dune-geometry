// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_GEOMETRY_REFERENCEELEMENTIMPLEMENTATION_HH
#define DUNE_GEOMETRY_REFERENCEELEMENTIMPLEMENTATION_HH

#include <cassert>

#include <algorithm>
#include <limits>
#include <tuple>
#include <utility>
#include <vector>
#include <array>
#include <bitset>

#include <dune/common/std/span.hh>

#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>
#include <dune/common/hybridutilities.hh>
#include <dune/common/typetraits.hh>
#include <dune/common/iteratorrange.hh>
#include <dune/common/math.hh>
#include <dune/common/indices.hh>
#include <dune/common/tuplevector.hh>


#include <dune/geometry/referenceelement.hh>
#include <dune/geometry/affinegeometry.hh>
#include <dune/geometry/type.hh>

namespace Dune
{

  namespace Geo
  {

#ifndef DOXYGEN

    // Internal Forward Declarations
    // -----------------------------

    namespace Impl
    {
      template< class ctype, int dim >
      class ReferenceElementContainer;
    }

    template< class ctype, int dim >
    struct ReferenceElements;



    namespace Impl
    {

      using Dune::Impl::isPrism;
      using Dune::Impl::isPyramid;
      using Dune::Impl::baseTopologyId;
      using Dune::Impl::prismConstruction;
      using Dune::Impl::pyramidConstruction;
      using Dune::Impl::numTopologies;

      // /** \brief Compute the number of subentities of a given codimension */
      // unsigned int size ( unsigned int topologyId, int dim, int codim );



      // /** \brief Compute the topology id of a given subentity
      //  *
      //  * \param topologyId Topology id of entity
      //  * \param dim Dimension of entity
      //  * \param codim Codimension of the subentity that we are interested in
      //  * \param i Number of the subentity that we are interested in
      //  */
      // unsigned int subTopologyId ( unsigned int topologyId, int dim, int codim, unsigned int i );



      // size
      // ----

      constexpr unsigned int size ( unsigned int topologyId, int dim, int codim )
      {
        assert( (dim >= 0) && (topologyId < numTopologies( dim )) );
        assert( (0 <= codim) && (codim <= dim) );

        if( codim > 0 )
          {
            const unsigned int baseId = baseTopologyId( topologyId, dim );
            const unsigned int m = size( baseId, dim-1, codim-1 );

            if( isPrism( topologyId, dim ) )
              {
                const unsigned int n = (codim < dim ? size( baseId, dim-1, codim ) : 0);
                return n + 2*m;
              }
            else
              {
                assert( isPyramid( topologyId, dim ) );
                const unsigned int n = (codim < dim ? size( baseId, dim-1, codim ) : 1);
                return m+n;
              }
          }
        else
          return 1;
      }



      // subTopologyId
      // -------------

      constexpr unsigned int subTopologyId ( unsigned int topologyId, int dim, int codim, unsigned int i )
      {
        assert( i < size( topologyId, dim, codim ) );
        const int mydim = dim - codim;

        if( codim > 0 )
          {
            const unsigned int baseId = baseTopologyId( topologyId, dim );
            const unsigned int m = size( baseId, dim-1, codim-1 );

            if( isPrism( topologyId, dim ) )
              {
                const unsigned int n = (codim < dim ? size( baseId, dim-1, codim ) : 0);
                if( i < n )
                  return subTopologyId( baseId, dim-1, codim, i ) | ((unsigned int)prismConstruction << (mydim - 1));
                else
                  return subTopologyId( baseId, dim-1, codim-1, (i < n+m ? i-n : i-(n+m)) );
              }
            else
              {
                assert( isPyramid( topologyId, dim ) );
                if( i < m )
                  return subTopologyId( baseId, dim-1, codim-1, i );
                else if( codim < dim )
                  return subTopologyId( baseId, dim-1, codim, i-m ) | ((unsigned int)pyramidConstruction << (mydim - 1));
                else
                  return 0u;
              }
          }
        else
          return topologyId;
      }


      // // subTopologyNumbering
      // // --------------------

      // void subTopologyNumbering ( unsigned int topologyId, int dim, int codim, unsigned int i, int subcodim,
      //                             unsigned int *beginOut, unsigned int *endOut );



      // subTopologyNumbering
      // --------------------

      constexpr void subTopologyNumbering ( unsigned int topologyId, int dim, int codim, unsigned int i, int subcodim,
                                  unsigned int *beginOut, unsigned int *endOut )
      {
        assert( (codim >= 0) && (subcodim >= 0) && (codim + subcodim <= dim) );
        assert( i < size( topologyId, dim, codim ) );
        assert( (endOut - beginOut) == size( subTopologyId( topologyId, dim, codim, i ), dim-codim, subcodim ) );

        if( codim == 0 )
          {
            for( unsigned int j = 0; (beginOut + j) != endOut; ++j )
              *(beginOut + j) = j;
          }
        else if( subcodim == 0 )
          {
            assert( endOut == beginOut + 1 );
            *beginOut = i;
          }
        else
          {
            const unsigned int baseId = baseTopologyId( topologyId, dim );

            const unsigned int m = size( baseId, dim-1, codim-1 );

            const unsigned int mb = size( baseId, dim-1, codim+subcodim-1 );
            const unsigned int nb = (codim + subcodim < dim ? size( baseId, dim-1, codim+subcodim ) : 0);

            if( isPrism( topologyId, dim ) )
              {
                const unsigned int n = size( baseId, dim-1, codim );
                if( i < n )
                  {
                    const unsigned int subId = subTopologyId( baseId, dim-1, codim, i );

                    unsigned int *beginBase = beginOut;
                    if( codim + subcodim < dim )
                      {
                        beginBase = beginOut + size( subId, dim-codim-1, subcodim );
                        subTopologyNumbering( baseId, dim-1, codim, i, subcodim, beginOut, beginBase );
                      }

                    const unsigned int ms = size( subId, dim-codim-1, subcodim-1 );
                    subTopologyNumbering( baseId, dim-1, codim, i, subcodim-1, beginBase, beginBase+ms );
                    for( unsigned int j = 0; j < ms; ++j )
                      {
                        *(beginBase+j) += nb;
                        *(beginBase+j+ms) = *(beginBase+j) + mb;
                      }
                  }
                else
                  {
                    const unsigned int s = (i < n+m ? 0 : 1);
                    subTopologyNumbering( baseId, dim-1, codim-1, i-(n+s*m), subcodim, beginOut, endOut );
                    for( unsigned int *it = beginOut; it != endOut; ++it )
                      *it += nb + s*mb;
                  }
              }
            else
              {
                assert( isPyramid( topologyId, dim ) );

                if( i < m )
                  subTopologyNumbering( baseId, dim-1, codim-1, i, subcodim, beginOut, endOut );
                else
                  {
                    const unsigned int subId = subTopologyId( baseId, dim-1, codim, i-m );
                    const unsigned int ms = size( subId, dim-codim-1, subcodim-1 );

                    subTopologyNumbering( baseId, dim-1, codim, i-m, subcodim-1, beginOut, beginOut+ms );
                    if( codim+subcodim < dim )
                      {
                        subTopologyNumbering( baseId, dim-1, codim, i-m, subcodim, beginOut+ms, endOut );
                        for( unsigned int *it = beginOut + ms; it != endOut; ++it )
                          *it += mb;
                      }
                    else
                      *(beginOut + ms) = mb;
                  }
              }
          }
      }




      // checkInside
      // -----------

      template< class ct, int cdim >
      inline bool
      checkInside ( unsigned int topologyId, int dim, const FieldVector< ct, cdim > &x, ct tolerance, ct factor = ct( 1 ) )
      {
        assert( (dim >= 0) && (dim <= cdim) );
        assert( topologyId < numTopologies( dim ) );

        if( dim > 0 )
          {
            const ct baseFactor = (isPrism( topologyId, dim ) ? factor : factor - x[ dim-1 ]);
            if( (x[ dim-1 ] > -tolerance) && (factor - x[ dim-1 ] > -tolerance) )
              return checkInside< ct, cdim >( baseTopologyId( topologyId, dim ), dim-1, x, tolerance, baseFactor );
            else
              return false;
          }
        else
          return true;
      }



      // referenceCorners
      // ----------------

      template< class ct, int cdim >
      constexpr unsigned int
      referenceCorners ( unsigned int topologyId, int dim, FieldVector< ct, cdim > *corners )
      {
        assert( (dim >= 0) && (dim <= cdim) );
        assert( topologyId < numTopologies( dim ) );

        if( dim > 0 )
          {
            const unsigned int nBaseCorners
              = referenceCorners( baseTopologyId( topologyId, dim ), dim-1, corners );
            assert( nBaseCorners == size( baseTopologyId( topologyId, dim ), dim-1, dim-1 ) );
            if( isPrism( topologyId, dim ) )
              {
                // pre: first face of the prism has correct corners
                // copy corners of first face to the opposite side (second half of the corners)
                for(std::size_t i = 0; i != nBaseCorners; ++i)
                  corners[ i+nBaseCorners ] = corners[i];
                // move opposite face by one
                for( unsigned int i = 0; i < nBaseCorners; ++i )
                  corners[ i+nBaseCorners ][ dim-1 ] = ct( 1 );
                return 2*nBaseCorners;
              }
            else
              {
                corners[ nBaseCorners ] = FieldVector< ct, cdim >( ct( 0 ) );
                corners[ nBaseCorners ][ dim-1 ] = ct( 1 );
                return nBaseCorners+1;
              }
          }
        else
          {
            *corners = FieldVector< ct, cdim >( ct( 0 ) );
            return 1;
          }
      }



      // referenceVolume
      // ---------------

      // unsigned long referenceVolumeInverse ( unsigned int topologyId, int dim );

      // ReferenceVolumeInverse
      // ----------------------

      constexpr unsigned long referenceVolumeInverse ( unsigned int topologyId, int dim )
      {
        assert( (dim >= 0) && (topologyId < numTopologies( dim )) );

        if( dim > 0 )
          {
            unsigned int baseValue = referenceVolumeInverse( baseTopologyId( topologyId, dim ), dim-1 );
            return (isPrism( topologyId, dim ) ? baseValue : baseValue * (unsigned long)dim);
          }
        else
          return 1;
      }


      template< class ct >
      constexpr ct referenceVolume ( unsigned int topologyId, int dim )
      {
        return ct( 1 ) / ct( referenceVolumeInverse( topologyId, dim ) );
      }



      // referenceOrigins
      // ----------------

      template< class ct, int cdim >
      constexpr unsigned int
      referenceOrigins ( unsigned int topologyId, int dim, int codim, FieldVector< ct, cdim > *origins )
      {
        assert( (dim >= 0) && (dim <= cdim) );
        assert( topologyId < numTopologies( dim ) );
        assert( (codim >= 0) && (codim <= dim) );

        if( codim > 0 )
          {
            const unsigned int baseId = baseTopologyId( topologyId, dim );
            if( isPrism( topologyId, dim ) )
              {
                const unsigned int n = (codim < dim ? referenceOrigins( baseId, dim-1, codim, origins ) : 0);
                const unsigned int m = referenceOrigins( baseId, dim-1, codim-1, origins+n );
                for( unsigned int i = 0; i < m; ++i )
                  {
                    origins[ n+m+i ] = origins[ n+i ];
                    origins[ n+m+i ][ dim-1 ] = ct( 1 );
                  }
                return n+2*m;
              }
            else
              {
                const unsigned int m = referenceOrigins( baseId, dim-1, codim-1, origins );
                if( codim == dim )
                  {
                    origins[ m ] = FieldVector< ct, cdim >( ct( 0 ) );
                    origins[ m ][ dim-1 ] = ct( 1 );
                    return m+1;
                  }
                else
                  return m+referenceOrigins( baseId, dim-1, codim, origins+m );
              }
          }
        else
          {
            origins[ 0 ] = FieldVector< ct, cdim >( ct( 0 ) );
            return 1;
          }
      }



      // referenceEmbeddings
      // -------------------

      template< class ct, int cdim, int mydim >
      constexpr unsigned int
      referenceEmbeddings ( unsigned int topologyId, int dim, int codim,
                            FieldVector< ct, cdim > *origins,
                            FieldMatrix< ct, mydim, cdim > *jacobianTransposeds )
      {
        assert( (0 <= codim) && (codim <= dim) && (dim <= cdim) );
        assert( (dim - codim <= mydim) && (mydim <= cdim) );
        assert( topologyId < numTopologies( dim ) );

        if( (0 < codim) && (codim <= dim) )
          {
            const unsigned int baseId = baseTopologyId( topologyId, dim );
            if( isPrism( topologyId, dim ) )
              {
                const unsigned int n = (codim < dim ? referenceEmbeddings( baseId, dim-1, codim, origins, jacobianTransposeds ) : 0);
                for( unsigned int i = 0; i < n; ++i )
                  jacobianTransposeds[ i ][ dim-codim-1 ][ dim-1 ] = ct( 1 );

                const unsigned int m = referenceEmbeddings( baseId, dim-1, codim-1, origins+n, jacobianTransposeds+n );
                for(std::size_t i = 0; i != m; ++i){
                  origins[n+m+i] = origins[n+i];
                  jacobianTransposeds[n+m+i] = jacobianTransposeds[n+i];
                }
                for( unsigned int i = 0; i < m; ++i )
                  origins[ n+m+i ][ dim-1 ] = ct( 1 );

                return n+2*m;
              }
            else // !isPrism
              {
                const unsigned int m = referenceEmbeddings( baseId, dim-1, codim-1, origins, jacobianTransposeds );
                if( codim == dim )
                  {
                    origins[ m ] = FieldVector< ct, cdim >( ct( 0 ) );
                    origins[ m ][ dim-1 ] = ct( 1 );
                    jacobianTransposeds[ m ] = FieldMatrix< ct, mydim, cdim >( ct( 0 ) );
                    return m+1;
                  }
                else if( codim < dim )
                  {
                    const unsigned int n = referenceEmbeddings( baseId, dim-1, codim, origins+m, jacobianTransposeds+m );
                    for( unsigned int i = 0; i < n; ++i )
                      {
                        for( int k = 0; k < dim-1; ++k )
                          jacobianTransposeds[ m+i ][ dim-codim-1 ][ k ] = -origins[ m+i ][ k ];
                        jacobianTransposeds[ m+i ][ dim-codim-1 ][ dim-1 ] = ct( 1 );
                      }
                    return m+n;
                  }
              }
          }
        else if( codim == 0 )
          {
            origins[ 0 ] = FieldVector< ct, cdim >( ct( 0 ) );
            jacobianTransposeds[ 0 ] = FieldMatrix< ct, mydim, cdim >( ct( 0 ) );
            for( int k = 0; k < dim; ++k )
              jacobianTransposeds[ 0 ][ k ][ k ] = ct( 1 );
            return 1;
          }

        // this point should not be reached since all cases are handled before.
        std::abort();
        return 0;
      }



      // referenceIntegrationOuterNormals
      // --------------------------------

      template< class ct, int cdim >
      constexpr unsigned int
      referenceIntegrationOuterNormals ( unsigned int topologyId, int dim,
                                         const FieldVector< ct, cdim > *origins,
                                         FieldVector< ct, cdim > *normals )
      {
        assert( (dim > 0) && (dim <= cdim) );
        assert( topologyId < numTopologies( dim ) );

        if( dim > 1 )
          {
            const unsigned int baseId = baseTopologyId( topologyId, dim );
            if( isPrism( topologyId, dim ) )
              {
                const unsigned int numBaseFaces
                  = referenceIntegrationOuterNormals( baseId, dim-1, origins, normals );

                for( unsigned int i = 0; i < 2; ++i )
                  {
                    normals[ numBaseFaces+i ] = FieldVector< ct, cdim >( ct( 0 ) );
                    normals[ numBaseFaces+i ][ dim-1 ] = ct( 2*int( i )-1 );
                  }

                return numBaseFaces+2;
              }
            else
              {
                normals[ 0 ] = FieldVector< ct, cdim >( ct( 0 ) );
                normals[ 0 ][ dim-1 ] = ct( -1 );

                const unsigned int numBaseFaces
                  = referenceIntegrationOuterNormals( baseId, dim-1, origins+1, normals+1 );
                for( unsigned int i = 1; i <= numBaseFaces; ++i )
                  normals[ i ][ dim-1 ] = normals[ i ]*origins[ i ];

                return numBaseFaces+1;
              }
          }
        else
          {
            for( unsigned int i = 0; i < 2; ++i )
              {
                normals[ i ] = FieldVector< ct, cdim >( ct( 0 ) );
                normals[ i ][ 0 ] = ct( 2*int( i )-1 );
              }

            return 2;
          }
      }

      template< int dim, unsigned int topologyId, class ct, int cdim>
      constexpr unsigned int
      referenceIntegrationOuterNormals ( FieldVector< ct, cdim > *normals )
      {
        assert( (dim > 0) && (dim <= cdim) );

        constexpr std::size_t sz = size( topologyId, dim, 1 );
        std::array<FieldVector< ct, cdim >, sz> origins;
        referenceOrigins( topologyId, dim, 1, origins.data() );

        const unsigned int numFaces
          = referenceIntegrationOuterNormals( topologyId, dim, origins.data(), normals );
        assert( numFaces == size( topologyId, dim, 1 ) );

        return numFaces;
      }

    } // namespace Impl



    // ReferenceElement
    // ----------------

    /** \class ReferenceElementImplementation
     *  \ingroup GeometryReferenceElements
     *  \brief This class provides access to geometric and topological
     *  properties of a reference element.
     *
     *  This includes its type,
     *  the number of subentities, the volume, and a method for checking
     *  if a point is contained in the reference element.
     *  The embedding of each subentity into the reference element is also
     *  provided.
     *
     *  A singleton of this class for a given geometry type can be accessed
     *  through the ReferenceElements class.

     *  \tparam ctype  field type for coordinates
     *  \tparam dim    dimension of the reference element
     *
     */
    template< class ctype_, int dim >
    class ReferenceElementImplementation
    {

    public:

      //! The coordinate field type.
      using ctype = ctype_;

      //! The coordinate field type.
      using CoordinateField = ctype;

      //! The coordinate type.
      using Coordinate = Dune::FieldVector<ctype,dim>;

      //! The dimension of the reference element.
      static constexpr int dimension = dim;

      /** \brief Type used for volume */
      typedef ctype Volume;

    private:

      friend class Impl::ReferenceElementContainer< ctype, dim >;

      struct SubEntityInfo;

      template< int codim > struct CreateGeometries;

    public:
      /** \brief Collection of types depending on the codimension */
      template< int codim >
      struct Codim
      {
        //! type of geometry embedding a subentity into the reference element
        typedef AffineGeometry< ctype, dim-codim, dim > Geometry;
      };

      // ReferenceElement cannot be copied.
      constexpr ReferenceElementImplementation ( const ReferenceElementImplementation& ) = default;

      // ReferenceElementImplementation cannot be copied.
      constexpr ReferenceElementImplementation& operator= ( const ReferenceElementImplementation& ) = default;

      // ReferenceElementImplementation is default-constructible (required for storage in std::array)
      constexpr ReferenceElementImplementation () = default;

      /** \brief number of subentities of codimension c
       *
       *  \param[in]  c  codimension whose size is desired
       */
      constexpr int size ( int c ) const
      {
        assert( (c >= 0) && (c <= dim) );
        return info_[ c ].size();
      }

      /** \brief number of subentities of codimension cc of subentity (i,c)
       *
       *  Denote by E the i-th subentity of codimension c of the current
       *  reference element. This method returns the number of subentities
       *  of codimension cc of the current reference element, that are also
       *  a subentity of E. If cc<c this number is zero.
       *
       *  \param[in]  i   number of subentity E (0 <= i < size( c ))
       *  \param[in]  c   codimension of subentity E (0 <= c <= dim)
       *  \param[in]  cc  codimension whose size is desired (0 <= cc <= dim)
       */
      constexpr int size ( int i, int c, int cc ) const
      {
        assert( (i >= 0) && (i < size( c )) );
        return info_[ c ][ i ].size( cc );
      }

      /** \brief obtain number of ii-th subentity with codim cc of (i,c)
       *
       *  Denote by E the i-th subentity of codimension c of the current
       *  reference element. And denote by S the ii-th subentity of codimension
       *  (cc-c) of E. Then, S is a also a subentity of codimension cc of the current
       *  reference element. This method returns the number of S with respect
       *  to the current reference element.
       *
       *  \param[in]  i   number of subentity E (0 <= i < size( c ))
       *  \param[in]  c   codimension of subentity E
       *  \param[in]  ii  number of subentity S (with respect to E)
       *  \param[in]  cc  codimension of subentity S (c <= cc <= dim)
       */
      constexpr int subEntity ( int i, int c, int ii, int cc ) const
      {
        assert( (i >= 0) && (i < size( c )) );
        return info_[ c ][ i ].number( ii, cc );
      }

      /** \brief Obtain the range of numbers of subentities with codim cc of (i,c)
       *
       *  Denote by E the i-th subentity of codimension c of the current
       *  reference element. This method returns a range of numbers of
       *  all subentities of E with codimension cc. Notice that the sub-subentity
       *  codimension as well as the numbers in the returned range are
       *  given with respect to the reference element itself and not with
       *  respect to E. For 0<=cc<c this will return an empty range.
       *
       *  \param[in]  i   number of subentity E (0 <= i < size( c ))
       *  \param[in]  c   codimension of subentity E
       *  \param[in]  cc  codimension of subentity S (0 <= cc <= dim)
       *
       *  \returns An iterable range of numbers of the sub-subentities.
       */
      constexpr auto subEntities ( int i, int c, int cc ) const
      {
        assert( (i >= 0) && (i < size( c )) );
        return info_[ c ][ i ].numbers( cc );
      }

      /** \brief obtain the type of subentity (i,c)
       *
       *  Denote by E the i-th subentity of codimension c of the current
       *  reference element. This method returns the GeometryType of E.
       *
       *  \param[in]  i      number of subentity E (0 <= i < size( c ))
       *  \param[in]  c      codimension of subentity E
       */
      constexpr GeometryType type ( int i, int c ) const
      {
        assert( (i >= 0) && (i < size( c )) );
        return info_[ c ][ i ].type();
      }

      /** \brief obtain the type of this reference element */
      constexpr GeometryType type () const { return type( 0, 0 ); }

      /** \brief position of the barycenter of entity (i,c)
       *
       *  Denote by E the i-th subentity of codimension c of the current
       *  reference element. This method returns the coordinates of
       *  the center of gravity of E within the current reference element.
       *
       *  \param[in]  i   number of subentity E (0 <= i < size( c ))
       *  \param[in]  c   codimension of subentity E
       */
      constexpr const Coordinate &position( int i, int c ) const
      {
        assert( (c >= 0) && (c <= dim) );
        return baryCenters_[ c ][ i ];
      }

      /** \brief check if a coordinate is in the reference element
       *
       *  This method returns true if the given local coordinate is within this
       *  reference element.
       *
       *  \param[in]  local  coordinates of the point
       */
      constexpr bool checkInside ( const Coordinate &local ) const
      {
        const ctype tolerance = ctype( 64 ) * std::numeric_limits< ctype >::epsilon();
        return Impl::template checkInside< ctype, dim >( type().id(), dim, local, tolerance );
      }

      /** \brief obtain the embedding of subentity (i,codim) into the reference
       *         element
       *
       *  Denote by E the i-th subentity of codimension codim of the current
       *  reference element. This method returns a \ref Dune::AffineGeometry
       *  that maps the reference element of E into the current reference element.
       *
       *  \tparam     codim  codimension of subentity E
       *
       *  \param[in]  i      number of subentity E (0 <= i < size( codim ))
       */
      template< int codim >
      constexpr typename Codim< codim >::Geometry geometry ( int i ) const
      {
        return std::get< codim >( geometries_ )[ i ];
      }

      /** \brief obtain the volume of the reference element */
      constexpr Volume volume () const
      {
        return volume_;
      }

      /** \brief obtain the integration outer normal of the reference element
       *
       *  The integration outer normal is the outer normal whose length coincides
       *  with the face's integration element.
       *
       *  \param[in]  face  index of the face, whose normal is desired
       */
      const Coordinate &integrationOuterNormal ( int face ) const
      {
        assert( (face >= 0) && (face < int( integrationNormals_.size() )) );
        return integrationNormals_[ face ];
      }

    private:

      // creates an nested tuple for every possible topologyId
      template<class F>
      static constexpr auto makeContainer(F f){
        return unpackIntegerSequence([=](auto... topologyId){
          return makeTupleVector(f(topologyId)...);
        }, std::make_index_sequence<Impl::numTopologies(dim)>());
      }

      // creates a container for all (topologyId, codim, i)
      static constexpr auto infos = makeContainer([](auto topologyId){
        auto make_entity_info = [=](auto codim){
            constexpr unsigned int size = Impl::size( topologyId, dim, codim );
            std::array<SubEntityInfo, size> info{};
            for( unsigned int i = 0; i < size; ++i )
              info[ i ] = SubEntityInfo{ topologyId, codim, i };
            return info;
        };

        return unpackIntegerSequence([=](auto... codim){
          return makeTupleVector(make_entity_info(codim)...);
        }, std::make_index_sequence<dim+1>());
      });

      // compute integration outer normals
      static constexpr auto integrationNormals = makeContainer([](auto topologyId){
        if constexpr ( dim > 0 )
          {
            using namespace Dune::Indices::Literals;
            std::array<Coordinate, infos[topologyId][1_ic].size()> normals{};
            Impl::referenceIntegrationOuterNormals<dim, topologyId>( &(normals[ 0 ]) );
            return normals;
          } else {
            return std::array<Coordinate,0>{};
          }
      });

      // creates a container for all baricenters
      static constexpr auto baryCenters = makeContainer([](auto topologyId){

        // compute corners
        const std::size_t numVertices = infos[topologyId][index_constant<dim>()].size();
        std::array<Coordinate, numVertices> baryCenter0{};
        Impl::referenceCorners( topologyId, dim, baryCenter0.data() );


        auto make_codim_barycenter = [=](auto codim){
          constexpr std::size_t sz = infos[topologyId][codim].size();
          std::array<Coordinate, sz> baryCenter{};
            for( int i = 0; i < sz; ++i )
              {
                baryCenter[ i ] = Coordinate( ctype( 0 ) );
                auto& info = infos[topologyId][codim][i];
                const unsigned int numCorners = info.size( dim );
                for( unsigned int j = 0; j < numCorners; ++j )
                  baryCenter[ i ] += baryCenter0[ info.number( j, dim )];
                baryCenter[ i ] *= ctype( 1 ) / ctype( numCorners );
              }
            return baryCenter;
        };

        return unpackIntegerSequence([=](auto... codim){
          return makeTupleVector(make_codim_barycenter(codim)..., baryCenter0);
        }, std::make_index_sequence<dim>());
      });

      // creates a container for all (topologyId, codim, i)
      static constexpr auto geometries = makeContainer([](auto topologyId){

        auto make_codim_geometry = [=](auto codim){

          constexpr std::size_t sz = infos[topologyId][codim].size();
          using namespace Indices::Literals;
          std::array< FieldVector< ctype, dim > , sz> origins{};
          std::array< FieldMatrix< ctype, dim - codim, dim >, sz > jacobianTransposeds{};
          Impl::referenceEmbeddings( infos[topologyId][0_ic][ 0 ].type().id(), dim, codim, origins.data(), jacobianTransposeds.data() );

          using Geometry = AffineGeometry< ctype, dim-codim, dim >;
          std::array<Geometry, sz> geometry{};

          for( std::size_t i = 0; i != sz; ++i )
            geometry[i] = Geometry{ infos[topologyId][ codim ][ i ].type(), origins[ i ], jacobianTransposeds[ i ] };
          return geometry;
        };

        return unpackIntegerSequence([=](auto... codim){
          return makeTupleVector(make_codim_geometry(codim)...);
        }, std::make_index_sequence<dim+1>());
      });

      constexpr void initialize ( unsigned int topologyId )
      {
        assert( topologyId < Impl::numTopologies( dim ) );

        Hybrid::switchCases(range(index_constant<Impl::numTopologies(dim)>()), topologyId, [&](auto topologyId_){
          integrationNormals_ = integrationNormals[topologyId_];
          Hybrid::forEach(range(index_constant<dim+1>()), [&](auto codim){
            info_[codim] = infos[topologyId_][codim];
            baryCenters_[codim] = baryCenters[topologyId_][codim];
            geometries_[codim] = geometries[topologyId_][codim];
          });
        });

        // compute reference element volume
        volume_ = Impl::template referenceVolume< ctype >( topologyId, dim );
      }

      template< int... codim >
      static TupleVector< Std::span< const typename Codim< codim >::Geometry >... >
      makeGeometryTable ( std::integer_sequence< int, codim... > );

      /** \brief Type to store all subentities of all codimensions */
      typedef decltype( makeGeometryTable( std::make_integer_sequence< int, dim+1 >() ) ) GeometryTable;

      /** \brief The reference element volume */
      ctype volume_ = {};

      std::array<Std::span<const Coordinate>, dim+1> baryCenters_;
      Std::span< const Coordinate > integrationNormals_;

      /** \brief Stores all subentities of all codimensions */
      GeometryTable geometries_;

      std::array<Std::span<const SubEntityInfo>, dim+1> info_;
    };

    /** \brief topological information about the subentities of a reference element */
    template< class ctype, int dim >
    struct ReferenceElementImplementation< ctype, dim >::SubEntityInfo
    {
      // creates an nested tuple for every possible (topologyId, codim, i)
      template<class F>
      static constexpr auto makeContainer(F f){
        auto make_size = [=](auto topologyId, auto codim){
          constexpr std::size_t size = Impl::size( topologyId, dim, codim );
          return unpackIntegerSequence([=](auto... i){
            return makeTupleVector(f(topologyId, codim, i)...);
          }, std::make_index_sequence<size>());
        };

        auto make_codims = [=](auto topologyId) {
          return unpackIntegerSequence([=](auto... codim){
            return makeTupleVector(make_size(topologyId, codim)...);
          }, std::make_index_sequence<dim+1>());
        };

        return unpackIntegerSequence([=](auto... topologyId){
          return makeTupleVector(make_codims(topologyId)...);
        }, std::make_index_sequence<Impl::numTopologies(dim)>());
      }

      // construct all possible offsets for every (topologyId, codim, i) combination
      constexpr static auto offsets = makeContainer([](auto topologyId, auto codim, auto i) {
        const unsigned int subId = Impl::subTopologyId( topologyId, dim, codim, i );

        // compute offsets
        std::array<unsigned int, dim+2> offset{};
        for( int cc = 0; cc <= codim; ++cc )
          offset[ cc ] = 0;
        for( int cc = codim; cc <= dim; ++cc )
          offset[ cc+1 ] = offset[ cc ] + Impl::size( subId, dim-codim, cc-codim );
        return offset;
      });

      // construct all possible numberings for every (topologyId, codim, i) combination
      constexpr static auto numberings = makeContainer([](auto topologyId, auto codim, auto i) {
        constexpr std::array<unsigned int, dim+2> offset = offsets[topologyId][codim][i];
        std::array<unsigned int, offset[ dim+1 ]> numbering{};
        for( int cc = codim; cc <= dim; ++cc )
          Impl::subTopologyNumbering( topologyId, dim, codim, i, cc-codim, numbering.data()+offset[ cc ], numbering.data()+offset[ cc+1 ] );
        return numbering;
      });


      // Compute upper bound for the number of subsentities.
      // If someone knows an explicit formal feel free to
      // implement it here.
      static constexpr std::size_t maxSubEntityCount()
      {
        std::size_t maxCount=0;
        for(std::size_t codim=0; codim<=dim; ++codim)
          maxCount = std::max(maxCount, binomial(std::size_t(dim),codim)*(1 << codim));
        return maxCount;
      }

      using SubEntityFlags = std::bitset<maxSubEntityCount()>;

      class SubEntityRange
        : public Dune::IteratorRange<const unsigned int*>
      {
        using Base = typename Dune::IteratorRange<const unsigned int*>;

      public:

        using iterator = Base::iterator;
        using const_iterator = Base::const_iterator;

        SubEntityRange(const iterator& begin, const iterator& end, const SubEntityFlags& contains) :
          Base(begin, end),
          containsPtr_(&contains),
          size_(end-begin)
        {}

        SubEntityRange() :
          Base(),
          containsPtr_(nullptr),
          size_(0)
        {}

        std::size_t size() const
        {
          return size_;
        }

        bool contains(std::size_t i) const
        {
          return (*containsPtr_)[i];
        }

      private:
        const SubEntityFlags* containsPtr_;
        std::size_t size_;
        std::size_t offset_;
      };

      using NumberRange = typename Dune::IteratorRange<const unsigned int*>;

      constexpr SubEntityInfo ()
        : numbering_(),
        offset_{}
      {
      }

      constexpr SubEntityInfo (unsigned int topologyId, int codim, unsigned int i)
        : SubEntityInfo( )
      {
        const unsigned int subId = Impl::subTopologyId( topologyId, dim, codim, i );
        type_ = GeometryType( subId, dim-codim );

        Hybrid::switchCases(range(index_constant<Impl::numTopologies(dim)>()), topologyId, [&](auto topologyId_){
          Hybrid::switchCases(range(index_constant<dim+1>()), codim, [&](auto codim_){
            constexpr std::size_t size = Impl::size( topologyId_, dim, codim_ );
            Hybrid::switchCases(range(index_constant<size>()), i, [&](auto i_){
              offset_ = offsets[topologyId_][codim_][i_];
              numbering_ = numberings[topologyId_][codim_][i_];
            });
          });
        });

        // initialize containsSubentity lookup-table
        for(std::size_t cc=0; cc<= dim; ++cc)
        {
          unsigned long cc_bitset = 0;
          for(std::size_t idx=0; idx<std::size_t(size(cc)); ++idx)
            cc_bitset |= (1 << number(idx,cc));
          containsSubentity_[cc] = SubEntityFlags{cc_bitset};
        }
      }

      constexpr SubEntityInfo ( const SubEntityInfo &other )
        : numbering_(other.numbering_),
          offset_( other.offset_ ),
          type_( other.type_ ),
          containsSubentity_( other.containsSubentity_ ),
      {}

      ~SubEntityInfo () = default;

      constexpr const SubEntityInfo &operator= ( const SubEntityInfo &other )
      {
        numbering_ = other.numbering_;
        offset_ = other.offset_;
        type_ = other.type_;
        containsSubentity_ = other.containsSubentity_;

        return *this;
      }

      constexpr int size ( int cc ) const
      {
        assert( (cc >= 0) && (cc <= dim) );
        return (offset_[ cc+1 ] - offset_[ cc ]);
      }

      constexpr int number ( int ii, int cc ) const
      {
        assert( (ii >= 0) && (ii < size( cc )) );
        return numbering_[ offset_[ cc ] + ii ];
      }

      constexpr auto numbers ( int cc ) const
      {
        return SubEntityRange( numbering_.begin() + offset_[ cc ], numbering_.begin() + offset_[ cc+1 ], containsSubentity_[cc]);
      }

      constexpr GeometryType type () const { return type_; }

    protected:
      constexpr int codim () const { return dim - type().dim(); }
      constexpr unsigned int capacity () const { return offset_[ dim+1 ]; }

    private:
      Std::span<const unsigned int> numbering_;
      std::array< unsigned int, dim+2 > offset_;
      GeometryType type_;
      std::array< SubEntityFlags, dim+1> containsSubentity_;
    };

#endif // DOXYGEN

  } // namespace Geo

} // namespace Dune

#endif // #ifndef DUNE_GEOMETRY_REFERENCEELEMENTIMPLEMENTATION_HH
