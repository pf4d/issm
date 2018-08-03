/*!\file TriMeshx
 * \brief: x code for TriMesh mesher
 */

/*Header files*/
#include "./CoordinateSystemTransformx.h"
#include "../../shared/shared.h"
#include "../../toolkits/toolkits.h"
#include <proj_api.h>

void CoordinateSystemTransformx(IssmDouble** px_dest,IssmDouble** py_dest,IssmDouble* x_src,IssmDouble* y_src,int size,const char* str_src,const char* str_dst){

#if !defined(_HAVE_PROJ4_)
	_error_("proj.4 has not been installed");
#else

	/*Allocate output and initialize values as src*/
	_assert_(size>0);
	IssmDouble* x_dest = xNew<IssmDouble>(size);
	IssmDouble* y_dest = xNew<IssmDouble>(size);
	for(int i=0;i<size;i++){
		x_dest[i] = x_src[i];
		y_dest[i] = y_src[i];
	}

	/*Create proj.4 projection objects for src and dst*/
	projPJ pj_src  = pj_init_plus(str_src);
	projPJ pj_dst  = pj_init_plus(str_dst);
	if(!pj_src) _error_("Failed to initialize PROJ.4 with source projection \""    <<str_src<<"\"\n");
	if(!pj_dst) _error_("Failed to initialize PROJ.4 with destination projection\""<<str_dst<<"\"\n");

	/*Perform transformation*/
	int p = pj_transform(pj_src,pj_dst,size,1,x_dest,y_dest,NULL);
	if(p!=0) _error_("Reprojection failed, PROJ.4 error code: "<<p<<"\n");

	/*Output : */
	*px_dest=x_dest;
	*py_dest=y_dest;
#endif
}
