/*!\file METIS_PartMeshNodalPatch
 * \brief common interface to Metis 4.0 and 5.0
 */

#include <config.h>
#include "../metisincludes.h"
#include "../../../shared/shared.h"

void METIS_PartMeshNodalPatch(int* pnumberofelements,int* pnumberofnodes, int* index, int* petype, int* pnumflag, int* pnum_procs, int* pedgecut, int* epart, int* npart){

	#if _METIS_VERSION_ == 4
	/*Our interface originates in the Metis 4.0 version, hence identical calls*/
	METIS_PartMeshNodal(pnumberofelements,pnumberofnodes, index, petype, pnumflag, pnum_procs, pedgecut, epart, npart); 
	#elif _METIS_VERSION_ == 5
	/*This interface is heavily changed. More options, different ways of meshing, etc ...: */ 
	int i;

	idx_t options[METIS_NOPTIONS];
	idx_t objval;
	idx_t* eptr=NULL;
	idx_t  k=0;
	real_t* tpwgts=NULL;

	/*setup options: */
	METIS_SetDefaultOptions(options);

	options[METIS_OPTION_PTYPE]   = 1;
	options[METIS_OPTION_OBJTYPE] = 0;
	options[METIS_OPTION_CTYPE]   = 1;
	options[METIS_OPTION_IPTYPE]  = 4;
	options[METIS_OPTION_RTYPE]   = 1;
	options[METIS_OPTION_DBGLVL]  = 0;
	options[METIS_OPTION_UFACTOR] = 30;
	options[METIS_OPTION_MINCONN] = 0;
	options[METIS_OPTION_CONTIG]  = 0;
	options[METIS_OPTION_SEED]    = -1;
	options[METIS_OPTION_NITER]   = 10;
	options[METIS_OPTION_NCUTS]   = 1;

	/*create eptr: */
	eptr=xNew<idx_t>((*pnumberofelements+1));
	eptr[0]=0;
	for(i=0;i<*pnumberofelements;i++){
		k+=3;
		eptr[i+1]=k;
	}

	/*create tpwgts: */
	tpwgts=xNew<real_t>(*pnum_procs);
	for(i=0;i<*pnum_procs;i++){
		tpwgts[i]=1.0/(*pnum_procs);
	}

	METIS_PartMeshNodal(pnumberofelements,pnumberofnodes, eptr, index,
			NULL, NULL, pnum_procs, tpwgts, options, &objval,epart, npart);

	/*clean-up*/
	xDelete<idx_t>(eptr);
	xDelete<real_t>(tpwgts);

	#else
	_error_("METIS version not supported yet");
	#endif
}
