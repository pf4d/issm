/*\file metispatches.h
 * \brief: our own patches for metis. Mainly to work through new apis from 4.0 to 5.0 version.
 */

#ifndef _METIS_PATCHES_H_
#define _METIS_PATCHES_H_

#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

void METIS_PartMeshNodalPatch(int *, int *, int *, int *, int *, int *, int *, int *, int *); //Common interface we are using in ISSM.

extern "C" {

#if _METIS_VERSION_ == 4
void METIS_PartMeshNodal(int *, int *, idxtype *, int *, int *, int *, int *, idxtype *, idxtype *);
#endif
#if _METIS_VERSION_ == 5
int METIS_PartMeshNodal(idx_t*, idx_t*, idx_t*, idx_t*, idx_t*, idx_t*, idx_t*, real_t*, idx_t*, idx_t*, idx_t*, idx_t*);
int METIS_SetDefaultOptions(idx_t *options);
#endif

}

#endif
