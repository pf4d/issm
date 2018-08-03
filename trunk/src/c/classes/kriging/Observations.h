#ifndef _CONTAINER_OBSERVATIONS_H_
#define  _CONTAINER_OBSERVATIONS_H_

class Quadtree;
class Covertree;
class Variogram;
class Options;
#include "../../datastructures/datastructures.h"
#include "../../shared/shared.h"

/*!\brief Declaration of Observations class.
 *
 * Declaration of Observations class.  Observations are vector lists (Containers) of Observation objects.
 */ 

class Observations: public DataSet{

	private:
		int        treetype;
		Quadtree*  quadtree;
		Covertree* covertree;

	public:

		/*constructors, destructors*/
		Observations();
		Observations(IssmDouble* observations_list,IssmDouble* x,IssmDouble* y,int n,Options* options);
		~Observations();

		/*Initialize data structures*/
		void InitCovertree(IssmDouble* observations_list,IssmDouble* x,IssmDouble* y,int n,Options* options);
		void InitQuadtree(IssmDouble* observations_list,IssmDouble* x,IssmDouble* y,int n,Options* options);

		/*Methods*/
		void ClosestObservation(IssmDouble *px,IssmDouble *py,IssmDouble *pobs,IssmDouble x_interp,IssmDouble y_interp,IssmDouble radius);
		void ClosestObservationCovertree(IssmDouble *px,IssmDouble *py,IssmDouble *pobs,IssmDouble x_interp,IssmDouble y_interp,IssmDouble radius);
		void ClosestObservationQuadtree(IssmDouble *px,IssmDouble *py,IssmDouble *pobs,IssmDouble x_interp,IssmDouble y_interp,IssmDouble radius);
		void Distances(IssmPDouble* distances,IssmPDouble *x,IssmPDouble *y,int n,IssmPDouble radius);
		void InterpolationIDW(IssmDouble *pprediction,IssmDouble x_interp,IssmDouble y_interp,IssmDouble radius,int mindata,int maxdata,IssmDouble power);
		void InterpolationV4(IssmDouble *pprediction,IssmDouble x_interp,IssmDouble y_interp,IssmDouble radius,int mindata,int maxdata);
		void InterpolationKriging(IssmDouble *pprediction,IssmDouble *perror,IssmDouble x_interp,IssmDouble y_interp,IssmDouble radius,int mindata,int maxdata,Variogram* variogram);
		void InterpolationNearestNeighbor(IssmDouble *pprediction,IssmDouble x_interp,IssmDouble y_interp,IssmDouble radius);
		void ObservationList(IssmDouble **px,IssmDouble **py,IssmDouble **pobs,int* pnobs);
		void ObservationList(IssmDouble **px,IssmDouble **py,IssmDouble **pobs,int* pnobs,IssmDouble x_interp,IssmDouble y_interp,IssmDouble radius,int maxdata);
		void ObservationListCovertree(IssmDouble **px,IssmDouble **py,IssmDouble **pobs,int* pnobs,IssmDouble x_interp,IssmDouble y_interp,IssmDouble radius,int maxdata);
		void ObservationListQuadtree(IssmDouble **px,IssmDouble **py,IssmDouble **pobs,int* pnobs,IssmDouble x_interp,IssmDouble y_interp,IssmDouble radius,int maxdata);
		void QuadtreeColoring(IssmDouble* A,IssmDouble *x,IssmDouble *y,int n);
		void Variomap(IssmDouble* gamma,IssmDouble *x,int n);

};
#endif //ifndef _OBSERVATIONS_H_
