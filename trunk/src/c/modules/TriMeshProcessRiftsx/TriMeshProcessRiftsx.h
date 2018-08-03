/*!\file:  TriMeshProcessRiftsx.h
 * \brief header file for TriMeshProcessRifts module
 */ 

#ifndef _TRIMESHPROCESSRIFTX_H
#define _TRIMESHPROCESSRIFTX_H

class RiftStruct;

void TriMeshProcessRiftsx(int** pindex,int* pnel,double** px,double** py,int* pnods,int** psegments,int** psegmentmarkers,int *pnum_seg,RiftStruct **priftstruct);

#endif  /* _TRIMESHPROCESSRIFTX_H*/
