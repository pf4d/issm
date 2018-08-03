/*
	InterpFromMeshToGrid.h
*/

#ifndef _INTERPFROMMESHTOGRID_H
#define _INTERPFROMMESHTOGRID_H

#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
	#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

/*For python modules: needs to come before header files inclusion*/
#ifdef _HAVE_PYTHON_
#define PY_ARRAY_UNIQUE_SYMBOL PythonIOSymbol
#endif

#include "../bindings.h"
#include "../../c/main/globals.h"
#include "../../c/toolkits/toolkits.h"
#include "../../c/modules/modules.h"
#include "../../c/shared/shared.h"
#include "../../c/shared/io/io.h"

#undef __FUNCT__ 
#define __FUNCT__  "InterpFromMeshToGrid"

#ifdef _HAVE_MATLAB_MODULES_
/* serial input macros: */
#define INDEX        prhs[0]
#define X            prhs[1]
#define Y            prhs[2]
#define MESHDATA     prhs[3]
#define XMIN         prhs[4]
#define YMAX         prhs[5]
#define XPOSTING     prhs[6]
#define YPOSTING     prhs[7]
#define NLINES       prhs[8]
#define NCOLS        prhs[9]
#define DEFAULTVALUE prhs[10]
/* serial output macros: */
#define XM (mxArray**)&plhs[0]
#define YM (mxArray**)&plhs[1]
#define GRIDDATA (mxArray**)&plhs[2]
#endif

#ifdef _HAVE_PYTHON_MODULES_
/* serial input macros: */
#define INDEX        PyTuple_GetItem(args,0)
#define X            PyTuple_GetItem(args,1)
#define Y            PyTuple_GetItem(args,2)
#define MESHDATA     PyTuple_GetItem(args,3)
#define XMIN         PyTuple_GetItem(args,4)
#define YMAX         PyTuple_GetItem(args,5)
#define XPOSTING     PyTuple_GetItem(args,6)
#define YPOSTING     PyTuple_GetItem(args,7)
#define NLINES       PyTuple_GetItem(args,8)
#define NCOLS        PyTuple_GetItem(args,9)
#define DEFAULTVALUE PyTuple_GetItem(args,10)
/* serial output macros: */
#define XM output,0
#define YM output,1
#define GRIDDATA output,2
#endif

/* serial arg counts: */
#undef NLHS
#define NLHS  3
#undef NRHS
#define NRHS  11

#endif  /* _INTERPFROMMESHTOGRID_H*/
