/*!\file SurfaceAverageVelMisfitx
 * \brief: compute misfit between observations and model
 */

#include "./SurfaceAverageVelMisfitx.h"

#include "../../shared/shared.h"
#include "../../toolkits/toolkits.h"
#include "../SurfaceAreax/SurfaceAreax.h"

void SurfaceAverageVelMisfitx(IssmDouble* pJ,FemModel* femmodel){

	/*output: */
	IssmDouble J=0.;
	IssmDouble J_sum;

	/*Compute Misfit: */
	for(int i=0;i<femmodel->elements->Size();i++){
		Element* element=xDynamicCast<Element*>(femmodel->elements->GetObjectByOffset(i));
		J+=SurfaceAverageVelMisfit(element);
	}

	/*Sum all J from all cpus of the cluster:*/
	ISSM_MPI_Reduce (&J,&J_sum,1,ISSM_MPI_DOUBLE,ISSM_MPI_SUM,0,IssmComm::GetComm() );
	ISSM_MPI_Bcast(&J_sum,1,ISSM_MPI_DOUBLE,0,IssmComm::GetComm());
	J=J_sum;

	/*Assign output pointers: */
	*pJ=J;
}

IssmDouble SurfaceAverageVelMisfit(Element* element){

	int        domaintype,numcomponents;
	IssmDouble Jelem=0.;
	IssmDouble misfit,S,Jdet;
	IssmDouble vx,vy,vxobs,vyobs,weight;
	IssmDouble* xyz_list = NULL;

	/*Get basal element*/
	if(!element->IsOnSurface()) return 0.;

	/*If on water, return 0: */
	if(!element->IsIceInElement()) return 0.;

	/*Get problem dimension*/
	element->FindParam(&domaintype,DomainTypeEnum);
	switch(domaintype){
		case Domain2DverticalEnum:
			numcomponents   = 1;
			break;
		case Domain3DEnum:
			numcomponents   = 2;
			break;
		case Domain2DhorizontalEnum:
			numcomponents   = 2;
			break;
		default: _error_("not supported yet");
	}

	/*Spawn surface element*/
	Element* topelement = element->SpawnTopElement();

	/* Get node coordinates*/
	topelement->GetVerticesCoordinates(&xyz_list);

	/*Retrieve all inputs we will be needing: */
	topelement->GetInputValue(&S,SurfaceAreaEnum);
	Input* weights_input=topelement->GetInput(InversionCostFunctionsCoefficientsEnum); _assert_(weights_input);
	Input* vx_input     =topelement->GetInput(VxEnum);                                 _assert_(vx_input);
	Input* vxobs_input  =topelement->GetInput(InversionVxObsEnum);                     _assert_(vxobs_input);
	Input* vy_input     = NULL;
	Input* vyobs_input  = NULL;
	if(numcomponents==2){
		vy_input    =topelement->GetInput(VyEnum);              _assert_(vy_input);
		vyobs_input =topelement->GetInput(InversionVyObsEnum);  _assert_(vyobs_input);
	}

	/* Start  looping on the number of gaussian points: */
	Gauss* gauss=topelement->NewGauss(3);
	for(int ig=gauss->begin();ig<gauss->end();ig++){

		gauss->GaussPoint(ig);

		/* Get Jacobian determinant: */
		topelement->JacobianDeterminant(&Jdet,xyz_list,gauss);

		/*Get all parameters at gaussian point*/
		weights_input->GetInputValue(&weight,gauss,SurfaceAverageVelMisfitEnum);
		vx_input->GetInputValue(&vx,gauss);
		vxobs_input->GetInputValue(&vxobs,gauss);
		if(numcomponents==2){
			vy_input->GetInputValue(&vy,gauss);
			vyobs_input->GetInputValue(&vyobs,gauss);
		}

		/*Compute SurfaceAverageVelMisfitEnum:
		 *
		 *      1                    2              2
		 * J = ---  sqrt(  (u - u   )  +  (v - v   )  )
		 *      S                obs            obs
		 */
		if(numcomponents==1){
			misfit=1/S*(vx-vxobs)*(vx-vxobs);
		}
		else{
			misfit=1/S*sqrt( pow(vx-vxobs,2) + pow(vy-vyobs,2));
		}

		/*Add to cost function*/
		Jelem+=misfit*weight*Jdet*gauss->weight;
	}

	/*clean up and Return: */
	if(domaintype!=Domain2DhorizontalEnum){topelement->DeleteMaterials(); delete topelement;};
	xDelete<IssmDouble>(xyz_list);
	delete gauss;
	return Jelem;
}
