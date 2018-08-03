/*!\file InputUpdateFromDakotax
 * \brief: update datasets using  parameter inputs
 */

#include "./InputUpdateFromDakotax.h"
#include "../../shared/shared.h"
#include "../../toolkits/toolkits.h"
#include "../InputUpdateFromMatrixDakotax/InputUpdateFromMatrixDakotax.h"
#include "../InputUpdateFromConstantx/InputUpdateFromConstantx.h"
#include "../InputUpdateFromVectorDakotax/InputUpdateFromVectorDakotax.h"

void InputUpdateFromDakotax(FemModel* femmodel,double* variables,char* *variables_descriptors,int numvariables){

	int     i,j,k,l;
	int     dummy;

	int     numberofvertices;
	int     nrows;
	int     ncols;
	int     npart;
	double *qmu_part  = NULL;

	double *distributed_values = NULL;
	double *parameter          = NULL;
	char   *descriptor         = NULL;
	char    root[50]; //root name of variable, ex: DragCoefficent, RhoIce, etc ...

	/*retrieve parameters: */
	femmodel->parameters->FindParam(&npart,QmuNumberofpartitionsEnum);
	femmodel->parameters->FindParam(&qmu_part,&dummy,QmuPartitionEnum);
	numberofvertices=femmodel->vertices->NumberOfVertices();

	/*Go through all dakota descriptors, ex: "rho_ice","thermal_conductivity","thickness1","thickness2", etc ..., and 
	 * for each descriptor, take the variable value and plug it into the inputs: */

	for(i=0;i<numvariables;i++){

		descriptor=variables_descriptors[i];

		/*From descriptor, figure out if the variable is scaled, indexed, nodal, or just a simple variable: */
		if (strncmp(descriptor,"scaled_",7)==0){

			/*Variable is scaled. Determine root name of variable (ex: scaled_DragCoefficient_1 -> DragCoefficient). Allocate distributed_values and fill the 
			 * distributed_values with the next npart variables: */

			//strcpy(root,strstr(descriptor,"_")+1); *strstr(root,"_")='\0';
			memcpy(root,strstr(descriptor,"_")+1,(strlen(strstr(descriptor,"_")+1)+1)*sizeof(char));
			*strstr(root,"_")='\0';

			distributed_values=xNew<double>(npart);
			for(j=0;j<npart;j++){
				distributed_values[j]=variables[i+j];
			}

			/*Now, pick up the parameter corresponding to root: */
			femmodel->parameters->FindParam(&parameter,&nrows,&ncols,StringToEnumx(root));

			/*We've got the parameter, we need to update it using qmu_part (a partitioning vector), 
			 * and the distributed_values. Two cases: we either have a nrows=numberofvertices, in 
			 * which case our parameter is a vector, or nrows=numberofvertices+1, in which case, 
			 * our parameter is a transient vector. Deal with both cases accordingly: */
			for(k=0;k<numberofvertices;k++){
				for(l=0;l<ncols;l++){
					*(parameter+ncols*k+l)=*(parameter+ncols*k+l)*distributed_values[(int)qmu_part[k]];
				}
			}

			#ifdef _DEBUG_
				PetscSynchronizedPrintf(IssmComm::GetComm(),"Parameter matrix:");
				PetscSynchronizedFlush(IssmComm::GetComm());
				for(l=0;l<ncols;l++){
					PetscSynchronizedPrintf(IssmComm::GetComm()," time %i\n",l);
					PetscSynchronizedFlush(IssmComm::GetComm());

					for(k=0;k<numberofvertices;k++){
						PetscSynchronizedPrintf(IssmComm::GetComm()," node %i value %g\n",k+1,*(parameter+k*ncols+l));
						PetscSynchronizedFlush(IssmComm::GetComm());
					}
				}
				PetscSynchronizedPrintf(IssmComm::GetComm()," descriptor: %s root %s enum: %i\n",descriptor,root,StringToEnumx(root));
				PetscSynchronizedFlush(IssmComm::GetComm());
			#endif

			/*Update inputs using the parameter matrix: */
			InputUpdateFromMatrixDakotax(femmodel, parameter, nrows,ncols,StringToEnumx(root), VertexEnum);

			/*increment i to skip the distributed values just collected: */
			i+=npart-1; //careful, the for loop will add 1.

			/*Free allocations: */
			xDelete<double>(parameter);
			xDelete<double>(distributed_values);
		}
		else if (strncmp(descriptor,"indexed_",8)==0){
			_error_("indexed variables not supported yet!");
		}
		else if (strncmp(descriptor,"nodal_",8)==0){
			_error_("nodal variables not supported yet!");
		}
		else{
			/*Ok, standard variable, just update inputs using the variable: */
			InputUpdateFromConstantx(femmodel,variables[i],StringToEnumx(descriptor));
		}
	}

	/*Free ressources:*/
	xDelete<double>(qmu_part);
}
