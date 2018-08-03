/*This routine only used by Intel compler: */

#ifdef HAVE_CONFIG_H
   #include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#include "./isnan.h"

#if defined(_HAVE_ADOLC_) && !defined(_WRAPPERS_)
template <> int xIsNan<adouble> (const adouble& X){
  return std::isnan(X.getValue());
}
#endif

#if defined(_HAVE_ADOLC_) && !defined(_WRAPPERS_)
template <> int xIsInf<adouble> (const adouble& X){
  return std::isinf(X.getValue());
}
#endif
