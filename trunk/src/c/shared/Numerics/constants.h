/*!\file: constants.h
 * \brief prototypes for constants.h
 */ 

#ifndef _ISSM_CONSTANTS_H_
#define _ISSM_CONSTANTS_H_

#define UNDEF -9999
#define ONETHIRD 0.333333333333333333333333333333333333333333333333333333333333
#define SQRT2 1.414213562373095048801688724209698078569671875376948073176679738
#define SQRT3 1.732050807568877293527446341505872366942805253810380628055806979
#define PI 3.141592653589793238462643383279502884197169399375105820974944592308

#define NDOF1 1
#define NDOF2 2
#define NDOF3 3
#define NDOF4 4

// /*Windows specific typefefs: */
// #ifdef _INTEL_WIN_
// 
// #ifndef NAN
// //For reference, for Intel compile on win64
// //#define NAN 0.0/0.0 
// #define NAN (INFINITY-INFINITY)
// #endif
// 
// #ifndef INFINITY
// //For reference, for Intel compile on win64
// //#define INFINITY 1.0/0.0
// #define INFINITY (DBL_MAX+DBL_MAX)
// #endif
// 
// #endif /*_INTEL_WIN_*/

#endif /*_ISSM_CONSTANTS_H_*/
