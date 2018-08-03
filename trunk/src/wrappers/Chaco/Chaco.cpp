/*\file Chaco.c
 *\brief:  Chaco partitioner mex module
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#include "./Chaco.h"

void ChacoUsage(void){/*{{{*/
	_printf0_("\n");
	_printf0_("Usage: [assgn] = Chaco(A,vwgts,ewgts,x,y,z,options,nparts,goal);\n");
	_printf0_("\n");
}/*}}}*/
WRAPPER(Chaco_python){

	int i;
	int nterms;

	/*Inputs: */
	int     nvtxs;               /* number of vertices in graph           */
	int    *start;               /* start of edge list for each vertex    */
	int    *adjacency;           /* edge list data                        */
	int    *vwgts       = NULL;  /* weights for all vertices              */
	int     nedges;
	float  *ewgts       = NULL;  /* weights for all edges                 */
	float  *x           = NULL;
	float  *y           = NULL;
	float  *z           = NULL;  /* coordinates for inertial method       */
	double  options[10] = {1,1,0,0,1,1,50,0,.001,7654321}; /* architecture and partitioning options */
	double *in_options  = NULL;
	int    *nparts      = NULL;   /* number of parts options               */
	int     npart;
	double *goal        = NULL;   /* desired set sizes                     */

	/*intermediary pointers: */
	mwIndex *mwstart, *mwadjacency;
	double  *doublepointer;

	/*output: */
   short  *assignment       = NULL; /* set number of each vtx (length nvtxs+1)                */
   double *doubleassignment = NULL; /*holds assignment, in double format, to return to matlab */

	#ifndef _HAVE_CHACO_ //only works if dakota library has been compiled in.
	_error_("Chaco not available! Cannot carry out Chaco partitioning!");
	#endif

	/*Boot module: */
	MODULEBOOT();

	/*checks on arguments on the matlab side: */
	CheckNumMatlabArguments(nlhs,NLHS,nrhs,NRHS,__FUNCT__,&ChacoUsage);

	/*Fetch adjacency matrix: */
	nvtxs = mxGetN(A_IN);
	mwstart = mxGetJc(A_IN);
	start=xNew<int>((nvtxs+1));
	for (i=0; i<nvtxs+1;i++)start[i]=(int)mwstart[i];

	mwadjacency = mxGetIr(A_IN);
	adjacency = xNew<int>(mxGetNzmax(A_IN));
	for (i=0; i<mxGetNzmax(A_IN); i++) adjacency[i]= (int)mwadjacency[i];

	nedges = start[nvtxs];
	if(!mxIsEmpty(EWGTS_IN)){
		ewgts = xNewZeroInit<float>(nedges);
		doublepointer=mxGetPr(A_IN);
		for (i = 0; i < nedges; i++)ewgts[i] = (float)doublepointer[i];
	}
	else ewgts=NULL;

	/*Fetch rest of data: */
	FetchData(&vwgts,&nterms,VWGTS_IN); 
	FetchData(&x,&nterms,X_IN); 
	FetchData(&y,&nterms,Y_IN); 
	FetchData(&z,&nterms,Z_IN); 
	FetchData(&in_options,&nterms,OPTNS_IN); 
	for (i=0;i<(nterms<10?nterms:10);i++) options[i]=in_options[i]; //copy in_options into default options
	FetchData(&npart,NPARTS_IN); 
	nparts=xNew<int>(1); nparts[0]=npart; //weird Chacox interface ain't it?
	FetchData(&goal,&nterms,GOAL_IN); 

	/*Some debugging print: {{{*/
	#ifdef _DEBUG_
	_printf_("nvtxs: " << nvtxs << "\n");
	_printf_("options: [");
	for(i=0;i<10;i++)_printf_(options[i] << "|");
	_printf_("]\n");
	_printf_("start: \n");
	for (i=0; i<nvtxs+1;i++)_printf_(start[i] << " ");
	_printf_("\n");
	_printf_("adjacency: \n");
	for (i=0; i<mxGetNzmax(A_IN);i++)_printf_("" <<adjacency[i]<< " ");i++)
	_printf_("\n");
	_printf_("nedges: " << nedges << " " << ewgts << "\n");
	if(ewgts) for (i = 0; i < nedges; i++)_printf_(ewgts[i] << " ");
	_printf_("\n");
	_printf_("vwgts:\n");
	for (i = 0; i < nvtxs; i++)_printf_(vwgts[i] << " ");
	_printf_("\n");
	_printf_("nparts: " << nparts[0] << "\n");
	_printf_("goal: " << goal << "\n");
	#endif
	/*}}}*/

	/*Allocate output: */
	assignment = xNewZeroInit<short>(nvtxs);

    /*Call core: */
	Chacox(nvtxs, start, adjacency, vwgts, ewgts, x, y, z, assignment, options, nparts, goal);

    /*Output data: */
	doubleassignment=xNew<double>(nvtxs);
	for(i=0;i<nvtxs;i++) doubleassignment[i]=(double)assignment[i];
	WriteData(ASSGN_OUT,doubleassignment,nvtxs);

	/*Free ressources:*/
	xDelete<short>(assignment); 
	xDelete<double>(goal);
	xDelete<int>(nparts);
	xDelete<float>(z);
	xDelete<float>(y);
	xDelete<float>(x);
	xDelete<float>(ewgts);
	xDelete<int>(vwgts);
	xDelete<int>(adjacency);
	xDelete<int>(start);
	xDelete<double>(doubleassignment);

	/*end module: */
	MODULEEND();
}
