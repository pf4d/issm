/*!\file:  ExpToLevelSetxt.cpp
 * \brief  "thread" core code for figuring out level set value from a contour and a cloud of points.
 */ 

#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

/*Include files: {{{*/
#include "./ExpToLevelSetx.h"
double minimum_distance(double x1, double y1, double x2, double y2, double x0, double y0);
void ContourToLevelSet(double* distance,double* contourx, double* contoury, int contournods, double* x, double* y, int i0, int i10);
/*}}}*/

void* ExpToLevelSetxt(void* vpthread_handle){

	/*gate variables :*/
	ExpToLevelSetxThreadStruct *gate        = NULL;
	pthread_handle             *handle      = NULL;
	int  i,i1,i0;

	/*recover handle and gate: */
	handle          = (pthread_handle*)vpthread_handle;
	gate            = (ExpToLevelSetxThreadStruct*)handle->gate;
	int my_thread   = handle->id;
	int num_threads = handle->num;

	/*recover parameters :*/
	Contours* contours  = gate->contours;
	int       nods      = gate->nods;
	double   *distance    = gate->distance;
	double   *x         = gate->x;
	double   *y         = gate->y;

	/*distribute indices across threads :*/
	PartitionRange(&i0,&i1,nods,num_threads,my_thread);

	/*Loop through all contours: */
	for (i=0;i<contours->Size();i++){
		Contour<double>* contour=(Contour<double>*)contours->GetObjectByOffset(i);
		ContourToLevelSet(distance,contour->x,contour->y,contour->nods,x,y,i0,i1);
	}

	return NULL;
}

void ContourToLevelSet(double* dist,double* contourx, double* contoury, int contournods, double* x, double* y, int i0, int i1){/*{{{*/
	int i,j;
	double x0,y0;
	double x1,y1;
	double x2,y2;
	double mind;
	
	for(i=i0;i<i1;i++){
		x0=x[i]; y0=y[i];

		/*Figure out distance from (x0,y0) to contour: */
		mind=1e+50;
		for (j=0;j<contournods-1;j++){
			x1=contourx[j]; y1=contoury[j];
			x2=contourx[j+1]; y2=contoury[j+1];
			mind=min(mind,minimum_distance(x1,y1,x2,y2,x0,y0));
		}
		dist[i]=min(dist[i],mind);
	}

}
double ddistance(double x1,double y1,double x2,double y2){
	return sqrt(pow(x2-x1,2)+pow(y2-y1,2));
}
double ddot(double x1, double y1, double x2, double y2){
	return x1*x2+y1*y2;
}

bool isPointLeftOfRay(double x, double y, double raySx, double raySy, double rayEx, double rayEy) {
	  return (y-raySy)*(rayEx-raySx)
		    >      (x-raySx)*(rayEy-raySy); 
}

double minimum_distance(double x1, double y1, double x2, double y2, double x0, double y0){
	
	// Return minimum distance between line segment [(x1,y1) (x2,y2)] and point (x0,y0) (v=(x1,y1), w=(x2,y2) and p=(x0,y0)
	double projectionx; 
	double projectiony; 
	double l2;
	double t;

	l2 = pow(x2-x1,2)+pow(y2-y1,2); // i.e. |w-v|^2 -  avoid a sqrt

	if (l2 == 0.0) return ddistance(x0,y0, x1,y1); // v == w case
	// Consider the line extending the segment, parameterized as v + t (w - v).
	//         // We find projection of point p onto the line. 
	//           // It falls where t = [(p-v) . (w-v)] / |w-v|^2
	t = ddot(x0-x1,y0-y1, x2-x1, y2-y1) / l2;
	if (t < 0.0) return ddistance(x0,y0, x1, y1);       // Beyond the 'v' end of the segment
	else if (t > 1.0) return ddistance(x0,y0, x2,y2);  // Beyond the 'w' end of the segment
	
	projectionx= x1 + t* (x2-x1);  // Projection falls on the segment
	projectiony= y1 + t* (y2-y1);
	return ddistance(x0, y0, projectionx, projectiony);
}
/*}}}*/
