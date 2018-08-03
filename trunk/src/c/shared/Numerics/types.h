/*!\file: types.h
 * \brief prototypes for types.h
 */ 

#ifndef _TYPES_H_
#define  _TYPES_H_

#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#include <stdio.h>

/*here are our abstracted types: inspired on petsc */
#if ISSM_USE_64BIT_INDICES == 1
typedef long long IssmInt;
#else
typedef int IssmInt;
#endif  

#if defined(_HAVE_ADOLC_) &&  !defined(_WRAPPERS_)
#include "adolc/adolc.h"
// for active variables
typedef adouble IssmDouble;
// for passive variables
typedef double IssmPDouble;
#else 
// see above
typedef double IssmDouble; 
// see above
typedef IssmDouble IssmPDouble;
#endif

#endif //ifndef _TYPES_H_
