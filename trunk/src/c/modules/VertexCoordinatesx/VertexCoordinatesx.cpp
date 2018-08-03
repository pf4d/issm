/*!\file VertexCoordinatesx
 * \brief: compute a vector x,y and z of vertex coordinates by 
 * marching through all our vertices. 
 */

#include "./VertexCoordinatesx.h"

#include "../../shared/shared.h"
#include "../../toolkits/toolkits.h"

void VertexCoordinatesx( IssmDouble** px, IssmDouble** py, IssmDouble** pz, Vertices* vertices,bool spherical) {

	/*output: */
	IssmDouble* x=NULL;
	IssmDouble* y=NULL;
	IssmDouble* z=NULL;

	Vector<IssmDouble>* vx=NULL;
	Vector<IssmDouble>* vy=NULL;
	Vector<IssmDouble>* vz=NULL;

	/*intermediary: */
	int  numberofvertices;
	int  i;

	/*figure out how many vertices we have: */
	numberofvertices=vertices->NumberOfVertices();

	vx=new Vector<IssmDouble>(numberofvertices);
	vy=new Vector<IssmDouble>(numberofvertices);
	vz=new Vector<IssmDouble>(numberofvertices);

	/*march through our vertices: */
	for(i=0;i<vertices->Size();i++){
		Vertex* vertex=(Vertex*)vertices->GetObjectByOffset(i);
		vertex->VertexCoordinates(vx,vy,vz,spherical);
	}

	/*Assemble*/
	vx->Assemble();
	vy->Assemble();
	vz->Assemble();

	/*serialize: */
	x=vx->ToMPISerial();
	y=vy->ToMPISerial();
	z=vz->ToMPISerial();

	/*Free ressources: */
	delete vx;
	delete vy;
	delete vz;

	/*output: */
	if (px)*px=x;
	else xDelete<IssmDouble>(x);
	if (py)*py=y;
	else xDelete<IssmDouble>(y);
	if (pz)*pz=z;
	else xDelete<IssmDouble>(z);
}
