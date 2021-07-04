// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_GEOMETRY_QUADRATURERULES_NUMBERCONVERSION_HH
#define DUNE_GEOMETRY_QUADRATURERULES_NUMBERCONVERSION_HH

#include <type_traits>

//! expand the number to a value and to a string
#define DUNE_NUMBER(type,value) \
  Dune::Impl::number< type >(value, #value)

namespace Dune {
namespace Impl {

//! Construct the number type `ct` from `double` or from a character sequence
template<typename ct>
ct number([[maybe_unused]] const double value, [[maybe_unused]] const char* str)
{
  if constexpr(std::is_constructible<ct,const char*>::value)
    return ct{str};
  else
    return value;
}

}} // end namespace Dune::Impl

#endif // DUNE_GEOMETRY_QUADRATURERULES_NUMBERCONVERSION_HH
