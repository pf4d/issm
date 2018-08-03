/*!\file:  TriMeshx.h
 * \brief header file for TriMeshx module
 */ 

#ifndef _TRIMESHX_H_
#define _TRIMESHX_H_

#include <string.h>
#include "../../classes/classes.h"

/* local prototypes: */
void TriMeshx(int** pindex,IssmPDouble** px,IssmPDouble** py,int** psegments,int** psegmentmarkerlist,int* pnels,int* pnods, int* pnseg,Contours* domain,Contours* rifts,double area);
#endif  /* _TRIMESHX_H */
