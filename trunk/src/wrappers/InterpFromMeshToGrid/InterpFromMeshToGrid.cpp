/*\file InterpFromMeshToGrid.c
 *\brief: compute diff between observed and modeled velocity
 */

#include "./InterpFromMeshToGrid.h"

void InterpFromMeshToGridUsage(void){/*{{{*/
	_printf0_("INTERPFROMMESHTOGRID - interpolation of a data defined on a mesh onto a grid\n");
	_printf0_("\n");
	_printf0_("   This function is a multi-threaded mex file that interpolates a field\n");
	_printf0_("   defined on a triangular mesh onto a regular grid\n");
	_printf0_("\n");
	_printf0_("   Usage:\n");
	_printf0_("      [x_m,y_m,griddata]=InterpFromMeshToGrid(index,x,y,data,xmin,ymax,xposting,yposting,nlines,ncols,default_value)\n");
	_printf0_("\n");
	_printf0_("      index,x,y: delaunay triangulation defining the mesh.\n");
	_printf0_("      meshdata: vertex values of data to be interpolated.\n");
	_printf0_("      xmin,ymax,posting,nlines,ncols: parameters that define the grid\n");
	_printf0_("      default_value: value of points located out of the mesh.\n");
	_printf0_("\n");
}/*}}}*/
WRAPPER(InterpFromMeshToGrid_python){

	/*input datasets: */
	double* index=NULL;
	int     nel;
	double* x=NULL;
	int     nods;
	double* y=NULL;
	double* meshdata=NULL;
	int     meshdata_length;
	double  xmin;
	double  ymax;
	double  xposting;
	double  yposting;
	int     nlines,ncols;
	double  default_value;

	/* output datasets: */
	double* griddata=NULL;
	double* x_m=NULL;
	double* y_m=NULL;

	/*Boot module: */
	MODULEBOOT();

	/*checks on arguments on the matlab side: */
	#ifdef _HAVE_MATLAB_MODULES_
	CheckNumMatlabArguments(nlhs,NLHS,nrhs,NRHS,__FUNCT__,&InterpFromMeshToGridUsage);
	#endif

	/*Input datasets: */
	FetchData(&index,&nel,NULL,INDEX);
	FetchData(&x,&nods,NULL,X);
	FetchData(&y,NULL,NULL,Y);
	FetchData(&meshdata,&meshdata_length,NULL,MESHDATA);
	FetchData(&xmin,XMIN);
	FetchData(&ymax,YMAX);
	FetchData(&xposting,XPOSTING);
	FetchData(&yposting,YPOSTING);
	FetchData(&nlines,NLINES);
	FetchData(&ncols,NCOLS);
	FetchData(&default_value,DEFAULTVALUE);

	/*Call core of computation: */
	InterpFromMeshToGridx(&x_m,&y_m,&griddata,index,x,y,nods,nel,meshdata,meshdata_length,xmin,ymax,xposting,yposting,nlines,ncols,default_value);

	/*Write results: */
	WriteData(XM,x_m,ncols);
	WriteData(YM,y_m,nlines);
	WriteData(GRIDDATA,griddata,nlines,ncols);

	/*Free ressources: */
	xDelete<double>(index);
	xDelete<double>(x);
	xDelete<double>(y);
	xDelete<double>(meshdata);
	xDelete<double>(griddata);
	xDelete<double>(x_m);
	xDelete<double>(y_m);

	/*end module: */
	MODULEEND();
}
