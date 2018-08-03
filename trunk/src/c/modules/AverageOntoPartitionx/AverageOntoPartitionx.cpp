/*!\file AverageOntoPartitionx
 * \brief: average vertex vector values onto a sub-partition of the vertices
 * used by scaled responses in Qmu analysis. See DakotaResponses module.
 */

#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#include "./AverageOntoPartitionx.h"
#include "../../shared/shared.h"
#include "../../toolkits/toolkits.h"

void AverageOntoPartitionx(double** paverage, Elements* elements, Nodes* nodes, Vertices* vertices, Loads* loads, Materials* materials, Parameters* parameters,double* vertex_response){

	int     dummy;
	int     npart;
	double *qmu_part  = NULL;

	/*output: */
	Vector<IssmDouble>* partition_contributions=NULL;
	Vector<IssmDouble>* partition_areas=NULL;
	Vector<IssmDouble>* vec_average=NULL;
	double* average=NULL;

	/*First, recover qmu partition of vertices: */
	parameters->FindParam(&qmu_part,&dummy,QmuPartitionEnum);

	/*Some parameters: */
	parameters->FindParam(&npart,QmuNumberofpartitionsEnum);

	/*average onto the separate areas. The result will be a npart sized vector. */

	/*allocate: */
	partition_contributions=new Vector<IssmDouble>(npart);
	partition_areas=new Vector<IssmDouble>(npart);
	vec_average=new Vector<IssmDouble>(npart);

	/*loop on each element, and add contribution of the element to the partition (surface weighted average): */
	for(int i=0;i<elements->Size();i++){
		Element* element=xDynamicCast<Element*>(elements->GetObjectByOffset(i));
		element->AverageOntoPartition(partition_contributions,partition_areas,vertex_response,qmu_part);
	}

	/*Assemble: */
	partition_contributions->Assemble();
	partition_areas->Assemble();

	/*We have the partition_areas and the partition_contributions for each partition -> compute the surfae weighted average: */
	vec_average->PointwiseDivide(partition_contributions,partition_areas);

	/*serialize:*/
	average=vec_average->ToMPISerial();

	/*Free ressources:*/
	xDelete<double>(qmu_part);
	delete partition_contributions;
	delete partition_areas;
	delete vec_average;

	/*Assign output pointers:*/
	*paverage=average;
}
