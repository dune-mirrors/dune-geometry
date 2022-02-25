#ifndef DUNE_GEOMETRY_UTILITY_IDENTITYMATRIX_HH
#define DUNE_GEOMETRY_UTILITY_IDENTITYMATRIX_HH

#include <dune/common/fmatrix.hh>

namespace Dune {

// Placeholder type for a trivial identity matrix without any functionality
struct IdentityMatrix
{
  // multiply Id * A
  template <class A>
  friend const A& operator* (IdentityMatrix, const A& a) { return a; }

  // multiply A * Id
  template <class A>
  friend const A& operator* (const A& a, IdentityMatrix) { return a; }

  // multiply Id * Id
  friend IdentityMatrix operator* (IdentityMatrix, IdentityMatrix) { return {}; }

  // cast into FieldMatrix
  template <class K, int n>
  operator FieldMatrix<K,n,n> () const
  {
    FieldMatrix<K,n,n> I;
    for (int i = 0; i < n; ++i)
      I[i][i] = K(1);
    return I;
  }
};

} // end namespace Dune

#endif // DUNE_GEOMETRY_UTILITY_IDENTITYMATRIX_HH
