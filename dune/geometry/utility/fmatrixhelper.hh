// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_GEOMETRY_UTILITY_FMATRIXHELPER_HH
#define DUNE_GEOMETRY_UTILITY_FMATRIXHELPER_HH

#include <cmath>

#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>

namespace Dune
{
  namespace Impl
  {

    // FieldMatrixHelper
    // -----------------

    template< class ct >
    struct FieldMatrixHelper
    {
      typedef ct ctype;

      template< int m, int n >
      static void Ax ( const FieldMatrix< ctype, m, n > &A, const FieldVector< ctype, n > &x, FieldVector< ctype, m > &ret )
      {
        for( int i = 0; i < m; ++i )
        {
          ret[ i ] = ctype( 0 );
          for( int j = 0; j < n; ++j )
            ret[ i ] += A[ i ][ j ] * x[ j ];
        }
      }

      template< int m, int n >
      static void ATx ( const FieldMatrix< ctype, m, n > &A, const FieldVector< ctype, m > &x, FieldVector< ctype, n > &ret )
      {
        for( int i = 0; i < n; ++i )
        {
          ret[ i ] = ctype( 0 );
          for( int j = 0; j < m; ++j )
            ret[ i ] += A[ j ][ i ] * x[ j ];
        }
      }

      template< int m, int n, int p >
      static void AB ( const FieldMatrix< ctype, m, n > &A, const FieldMatrix< ctype, n, p > &B, FieldMatrix< ctype, m, p > &ret )
      {
        for( int i = 0; i < m; ++i )
        {
          for( int j = 0; j < p; ++j )
          {
            ret[ i ][ j ] = ctype( 0 );
            for( int k = 0; k < n; ++k )
              ret[ i ][ j ] += A[ i ][ k ] * B[ k ][ j ];
          }
        }
      }

      template< int m, int n, int p >
      static void ATBT ( const FieldMatrix< ctype, m, n > &A, const FieldMatrix< ctype, p, m > &B, FieldMatrix< ctype, n, p > &ret )
      {
        for( int i = 0; i < n; ++i )
        {
          for( int j = 0; j < p; ++j )
          {
            ret[ i ][ j ] = ctype( 0 );
            for( int k = 0; k < m; ++k )
              ret[ i ][ j ] += A[ k ][ i ] * B[ j ][ k ];
          }
        }
      }

      template< int m, int n >
      static void ATA_L ( const FieldMatrix< ctype, m, n > &A, FieldMatrix< ctype, n, n > &ret )
      {
        for( int i = 0; i < n; ++i )
        {
          for( int j = 0; j <= i; ++j )
          {
            ret[ i ][ j ] = ctype( 0 );
            for( int k = 0; k < m; ++k )
              ret[ i ][ j ] += A[ k ][ i ] * A[ k ][ j ];
          }
        }
      }

      template< int m, int n >
      static void ATA ( const FieldMatrix< ctype, m, n > &A, FieldMatrix< ctype, n, n > &ret )
      {
        for( int i = 0; i < n; ++i )
        {
          for( int j = 0; j <= i; ++j )
          {
            ret[ i ][ j ] = ctype( 0 );
            for( int k = 0; k < m; ++k )
              ret[ i ][ j ] += A[ k ][ i ] * A[ k ][ j ];
            ret[ j ][ i ] = ret[ i ][ j ];
          }

          ret[ i ][ i ] = ctype( 0 );
          for( int k = 0; k < m; ++k )
            ret[ i ][ i ] += A[ k ][ i ] * A[ k ][ i ];
        }
      }

      template< int m, int n >
      static void AAT_L ( const FieldMatrix< ctype, m, n > &A, FieldMatrix< ctype, m, m > &ret )
      {
        /*
           if (m==2) {
           ret[0][0] = A[0]*A[0];
           ret[1][1] = A[1]*A[1];
           ret[1][0] = A[0]*A[1];
           }
           else
         */
        for( int i = 0; i < m; ++i )
        {
          for( int j = 0; j <= i; ++j )
          {
            ctype &retij = ret[ i ][ j ];
            retij = A[ i ][ 0 ] * A[ j ][ 0 ];
            for( int k = 1; k < n; ++k )
              retij += A[ i ][ k ] * A[ j ][ k ];
          }
        }
      }

      template< int m, int n >
      static void AAT ( const FieldMatrix< ctype, m, n > &A, FieldMatrix< ctype, m, m > &ret )
      {
        for( int i = 0; i < m; ++i )
        {
          for( int j = 0; j < i; ++j )
          {
            ret[ i ][ j ] = ctype( 0 );
            for( int k = 0; k < n; ++k )
              ret[ i ][ j ] += A[ i ][ k ] * A[ j ][ k ];
            ret[ j ][ i ] = ret[ i ][ j ];
          }
          ret[ i ][ i ] = ctype( 0 );
          for( int k = 0; k < n; ++k )
            ret[ i ][ i ] += A[ i ][ k ] * A[ i ][ k ];
        }
      }

      template< int n >
      static void Lx ( const FieldMatrix< ctype, n, n > &L, const FieldVector< ctype, n > &x, FieldVector< ctype, n > &ret )
      {
        for( int i = 0; i < n; ++i )
        {
          ret[ i ] = ctype( 0 );
          for( int j = 0; j <= i; ++j )
            ret[ i ] += L[ i ][ j ] * x[ j ];
        }
      }

      template< int n >
      static void LTx ( const FieldMatrix< ctype, n, n > &L, const FieldVector< ctype, n > &x, FieldVector< ctype, n > &ret )
      {
        for( int i = 0; i < n; ++i )
        {
          ret[ i ] = ctype( 0 );
          for( int j = i; j < n; ++j )
            ret[ i ] += L[ j ][ i ] * x[ j ];
        }
      }

      template< int n >
      static void LTL ( const FieldMatrix< ctype, n, n > &L, FieldMatrix< ctype, n, n > &ret )
      {
        for( int i = 0; i < n; ++i )
        {
          for( int j = 0; j < i; ++j )
          {
            ret[ i ][ j ] = ctype( 0 );
            for( int k = i; k < n; ++k )
              ret[ i ][ j ] += L[ k ][ i ] * L[ k ][ j ];
            ret[ j ][ i ] = ret[ i ][ j ];
          }
          ret[ i ][ i ] = ctype( 0 );
          for( int k = i; k < n; ++k )
            ret[ i ][ i ] += L[ k ][ i ] * L[ k ][ i ];
        }
      }

      template< int n >
      static void LLT ( const FieldMatrix< ctype, n, n > &L, FieldMatrix< ctype, n, n > &ret )
      {
        for( int i = 0; i < n; ++i )
        {
          for( int j = 0; j < i; ++j )
          {
            ret[ i ][ j ] = ctype( 0 );
            for( int k = 0; k <= j; ++k )
              ret[ i ][ j ] += L[ i ][ k ] * L[ j ][ k ];
            ret[ j ][ i ] = ret[ i ][ j ];
          }
          ret[ i ][ i ] = ctype( 0 );
          for( int k = 0; k <= i; ++k )
            ret[ i ][ i ] += L[ i ][ k ] * L[ i ][ k ];
        }
      }

      template< int n >
      static bool cholesky_L ( const FieldMatrix< ctype, n, n > &A, FieldMatrix< ctype, n, n > &ret, const bool checkSingular = false )
      {
        using std::sqrt;
        for( int i = 0; i < n; ++i )
        {
          ctype &rii = ret[ i ][ i ];

          ctype xDiag = A[ i ][ i ];
          for( int j = 0; j < i; ++j )
            xDiag -= ret[ i ][ j ] * ret[ i ][ j ];

          // in some cases A can be singular, e.g. when checking local for
          // outside points during checkInside
          if( checkSingular && ! ( xDiag > ctype( 0 )) )
            return false ;

          // otherwise this should be true always
          assert( xDiag > ctype( 0 ) );
          rii = sqrt( xDiag );

          ctype invrii = ctype( 1 ) / rii;
          for( int k = i+1; k < n; ++k )
          {
            ctype x = A[ k ][ i ];
            for( int j = 0; j < i; ++j )
              x -= ret[ i ][ j ] * ret[ k ][ j ];
            ret[ k ][ i ] = invrii * x;
          }
        }

        // return true for meaning A is non-singular
        return true;
      }

      template< int n >
      static ctype detL ( const FieldMatrix< ctype, n, n > &L )
      {
        ctype det( 1 );
        for( int i = 0; i < n; ++i )
          det *= L[ i ][ i ];
        return det;
      }

      template< int n >
      static ctype invL ( FieldMatrix< ctype, n, n > &L )
      {
        ctype det( 1 );
        for( int i = 0; i < n; ++i )
        {
          ctype &lii = L[ i ][ i ];
          det *= lii;
          lii = ctype( 1 ) / lii;
          for( int j = 0; j < i; ++j )
          {
            ctype &lij = L[ i ][ j ];
            ctype x = lij * L[ j ][ j ];
            for( int k = j+1; k < i; ++k )
              x += L[ i ][ k ] * L[ k ][ j ];
            lij = (-lii) * x;
          }
        }
        return det;
      }

      // calculates x := L^{-1} x
      template< int n >
      static void invLx ( FieldMatrix< ctype, n, n > &L, FieldVector< ctype, n > &x )
      {
        for( int i = 0; i < n; ++i )
        {
          for( int j = 0; j < i; ++j )
            x[ i ] -= L[ i ][ j ] * x[ j ];
          x[ i ] /= L[ i ][ i ];
        }
      }

      // calculates x := L^{-T} x
      template< int n >
      static void invLTx ( FieldMatrix< ctype, n, n > &L, FieldVector< ctype, n > &x )
      {
        for( int i = n; i > 0; --i )
        {
          for( int j = i; j < n; ++j )
            x[ i-1 ] -= L[ j ][ i-1 ] * x[ j ];
          x[ i-1 ] /= L[ i-1 ][ i-1 ];
        }
      }

      template< int n >
      static ctype spdDetA ( const FieldMatrix< ctype, n, n > &A )
      {
        // return A[0][0]*A[1][1]-A[1][0]*A[1][0];
        FieldMatrix< ctype, n, n > L;
        cholesky_L( A, L );
        return detL( L );
      }

      template< int n >
      static ctype spdInvA ( FieldMatrix< ctype, n, n > &A )
      {
        FieldMatrix< ctype, n, n > L;
        cholesky_L( A, L );
        const ctype det = invL( L );
        LTL( L, A );
        return det;
      }

      // calculate x := A^{-1} x
      template< int n >
      static bool spdInvAx ( FieldMatrix< ctype, n, n > &A, FieldVector< ctype, n > &x, const bool checkSingular = false )
      {
        FieldMatrix< ctype, n, n > L;
        const bool invertible = cholesky_L( A, L, checkSingular );
        if( ! invertible ) return invertible ;
        invLx( L, x );
        invLTx( L, x );
        return invertible;
      }

      template< int m, int n >
      static ctype detATA ( const FieldMatrix< ctype, m, n > &A )
      {
        if( m >= n )
        {
          FieldMatrix< ctype, n, n > ata;
          ATA_L( A, ata );
          return spdDetA( ata );
        }
        else
          return ctype( 0 );
      }

      /** \brief Compute the square root of the determinant of A times A transposed
       *
       *  This is the volume element for an embedded submanifold and needed to
       *  implement the method integrationElement().
       */
      template< int m, int n >
      static ctype sqrtDetAAT ( const FieldMatrix< ctype, m, n > &A )
      {
        using std::abs;
        using std::sqrt;
        // These special cases are here not only for speed reasons:
        // The general implementation aborts if the matrix is almost singular,
        // and the special implementation provide a stable way to handle that case.
        if( (n == 2) && (m == 2) )
        {
          // Special implementation for 2x2 matrices: faster and more stable
          return abs( A[ 0 ][ 0 ]*A[ 1 ][ 1 ] - A[ 1 ][ 0 ]*A[ 0 ][ 1 ] );
        }
        else if( (n == 3) && (m == 3) )
        {
          // Special implementation for 3x3 matrices
          const ctype v0 = A[ 0 ][ 1 ] * A[ 1 ][ 2 ] - A[ 1 ][ 1 ] * A[ 0 ][ 2 ];
          const ctype v1 = A[ 0 ][ 2 ] * A[ 1 ][ 0 ] - A[ 1 ][ 2 ] * A[ 0 ][ 0 ];
          const ctype v2 = A[ 0 ][ 0 ] * A[ 1 ][ 1 ] - A[ 1 ][ 0 ] * A[ 0 ][ 1 ];
          return abs( v0 * A[ 2 ][ 0 ] + v1 * A[ 2 ][ 1 ] + v2 * A[ 2 ][ 2 ] );
        }
        else if ( (n == 3) && (m == 2) )
        {
          // Special implementation for 2x3 matrices
          const ctype v0 = A[ 0 ][ 0 ] * A[ 1 ][ 1 ] - A[ 0 ][ 1 ] * A[ 1 ][ 0 ];
          const ctype v1 = A[ 0 ][ 0 ] * A[ 1 ][ 2 ] - A[ 1 ][ 0 ] * A[ 0 ][ 2 ];
          const ctype v2 = A[ 0 ][ 1 ] * A[ 1 ][ 2 ] - A[ 0 ][ 2 ] * A[ 1 ][ 1 ];
          return sqrt( v0*v0 + v1*v1 + v2*v2);
        }
        else if( n >= m )
        {
          // General case
          FieldMatrix< ctype, m, m > aat;
          AAT_L( A, aat );
          return spdDetA( aat );
        }
        else
          return ctype( 0 );
      }

      // A^{-1}_L = (A^T A)^{-1} A^T
      // => A^{-1}_L A = I
      template< int m, int n >
      static ctype leftInvA ( const FieldMatrix< ctype, m, n > &A, FieldMatrix< ctype, n, m > &ret )
      {
        static_assert((m >= n), "Matrix has no left inverse.");
        FieldMatrix< ctype, n, n > ata;
        ATA_L( A, ata );
        const ctype det = spdInvA( ata );
        ATBT( ata, A, ret );
        return det;
      }

      template< int m, int n >
      static void leftInvAx ( const FieldMatrix< ctype, m, n > &A, const FieldVector< ctype, m > &x, FieldVector< ctype, n > &y )
      {
        static_assert((m >= n), "Matrix has no left inverse.");
        FieldMatrix< ctype, n, n > ata;
        ATx( A, x, y );
        ATA_L( A, ata );
        spdInvAx( ata, y );
      }

      /** \brief Compute right pseudo-inverse of matrix A */
      template< int m, int n >
      static ctype rightInvA ( const FieldMatrix< ctype, m, n > &A, FieldMatrix< ctype, n, m > &ret )
      {
        static_assert((n >= m), "Matrix has no right inverse.");
        using std::abs;
        if( (n == 2) && (m == 2) )
        {
          const ctype det = (A[ 0 ][ 0 ]*A[ 1 ][ 1 ] - A[ 1 ][ 0 ]*A[ 0 ][ 1 ]);
          const ctype detInv = ctype( 1 ) / det;
          ret[ 0 ][ 0 ] = A[ 1 ][ 1 ] * detInv;
          ret[ 1 ][ 1 ] = A[ 0 ][ 0 ] * detInv;
          ret[ 1 ][ 0 ] = -A[ 1 ][ 0 ] * detInv;
          ret[ 0 ][ 1 ] = -A[ 0 ][ 1 ] * detInv;
          return abs( det );
        }
        else
        {
          FieldMatrix< ctype, m , m > aat;
          AAT_L( A, aat );
          const ctype det = spdInvA( aat );
          ATBT( A , aat , ret );
          return det;
        }
      }

      template< int m, int n >
      static bool xTRightInvA ( const FieldMatrix< ctype, m, n > &A, const FieldVector< ctype, n > &x, FieldVector< ctype, m > &y )
      {
        static_assert((n >= m), "Matrix has no right inverse.");
        FieldMatrix< ctype, m, m > aat;
        Ax( A, x, y );
        AAT_L( A, aat );
        // check whether aat is singular and return true if non-singular
        return spdInvAx( aat, y, true );
      }
    };

  } // namespace Impl
} // namespace Dune

#endif // DUNE_GEOMETRY_UTILITY_FMATRIXHELPER_HH
