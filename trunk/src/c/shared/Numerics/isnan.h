/*!\file:  isnan.h
 * \brief: only used for intel compiler.
 */ 

#ifndef _XISNAN_H_
#define _XISNAN_H_

#ifdef HAVE_CONFIG_H
   #include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

/*If include cmath instead of math, isnan on linux64 murdo does not work: */
#include <cmath>

template <class T> int xIsNan(const T& X) {
#ifdef _INTEL_WIN_
		return (X!=X)?1:0;
#else
		return std::isnan(X);
#endif
}

template <class T> int xIsInf(const T& X) {
	return std::isinf(X);
}

#if defined(_HAVE_ADOLC_) && !defined(_WRAPPERS_)
#include "./types.h"
template <> int xIsNan<adouble> (const adouble& X);
template <> int xIsInf<adouble> (const adouble& X);
#endif

#endif
