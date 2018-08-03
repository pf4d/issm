/*! \file OptionUtilities.h 
 *  \brief: header file for option object utilities
 */

#ifndef _OPTIONUTILITIES_H_
#define _OPTIONUTILITIES_H_

/*Headers:{{{*/
#include "../../shared/Exceptions/exceptions.h"
#include "../../shared/Enum/Enum.h"

#include "./Option.h"
/*}}}*/

int ColumnWiseDimsFromIndex(int* dims, int index, int* size, int ndims);
int IndexFromColumnWiseDims(int* dims, int* size, int ndims);
int IndexFromRowWiseDims(int* dims, int* size, int ndims);
int RowWiseDimsFromIndex(int* dims, int index, int* size, int ndims);
int StringFromDims(char* cstr, int* dims, int ndims);
int StringFromSize(char* cstr, int* size, int ndims);

#endif  /* _OPTIONUTILITIES_H */
