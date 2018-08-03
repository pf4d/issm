/*!\file RheologyBbarAbsGradientx
 * \brief: compute misfit between observations and model
 */

#include "./RheologyBbarAbsGradientx.h"

#include "../../shared/shared.h"
#include "../../toolkits/toolkits.h"

void RheologyBbarAbsGradientx( IssmDouble* pJ, Elements* elements,Nodes* nodes, Vertices* vertices, Loads* loads, Materials* materials,Parameters* parameters){

	/*output: */
	IssmDouble J=0.;
	IssmDouble J_sum;

	/*Compute Misfit: */
	for(int i=0;i<elements->Size();i++){
		Element* element=xDynamicCast<Element*>(elements->GetObjectByOffset(i));
		J+=RheologyBbarAbsGradient(element);
	}

	/*Sum all J from all cpus of the cluster:*/
	ISSM_MPI_Reduce (&J,&J_sum,1,ISSM_MPI_DOUBLE,ISSM_MPI_SUM,0,IssmComm::GetComm() );
	ISSM_MPI_Bcast(&J_sum,1,ISSM_MPI_DOUBLE,0,IssmComm::GetComm());
	J=J_sum;

	/*Assign output pointers: */
	*pJ=J;
}

IssmDouble RheologyBbarAbsGradient(Element* element){

	int         domaintype,numcomponents;
	IssmDouble  Jelem=0.;
	IssmDouble  misfit,Jdet;
	IssmDouble  dp[2],weight;
	IssmDouble* xyz_list      = NULL;

	/*Get basal element*/
	if(!element->IsOnBase()) return 0.;

	/*If on water, return 0: */
	if(!element->IsIceInElement()) return 0.;

	/*Get problem dimension*/
	element->FindParam(&domaintype,DomainTypeEnum);
	switch(domaintype){
		case Domain2DverticalEnum:   numcomponents   = 1; break;
		case Domain3DEnum:           numcomponents   = 2; break;
		case Domain2DhorizontalEnum: numcomponents   = 2; break;
		default: _error_("not supported yet");
	}

	/*Spawn basal element*/
	Element* basalelement = element->SpawnBasalElement();

	/* Get node coordinates*/
	basalelement->GetVerticesCoordinates(&xyz_list);

	/*Retrieve all inputs we will be needing: */
	Input* weights_input=basalelement->GetInput(InversionCostFunctionsCoefficientsEnum); _assert_(weights_input);
	Input* rheologyb_input=basalelement->GetInput(MaterialsRheologyBbarEnum);                  _assert_(rheologyb_input);

	/* Start  looping on the number of gaussian points: */
	Gauss* gauss=basalelement->NewGauss(2);
	for(int ig=gauss->begin();ig<gauss->end();ig++){

		gauss->GaussPoint(ig);

		/* Get Jacobian determinant: */
		basalelement->JacobianDeterminant(&Jdet,xyz_list,gauss);

		/*Get all parameters at gaussian point*/
		weights_input->GetInputValue(&weight,gauss,RheologyBbarAbsGradientEnum);
		rheologyb_input->GetInputDerivativeValue(&dp[0],xyz_list,gauss);

		/*Tikhonov regularization: J = 1/2 ((dp/dx)^2 + (dp/dy)^2) */ 
		Jelem+=weight*1/2*(dp[0]*dp[0] + dp[1]*dp[1])*Jdet*gauss->weight;
	}

	/*clean up and Return: */
	if(domaintype!=Domain2DhorizontalEnum){basalelement->DeleteMaterials(); delete basalelement;};
	xDelete<IssmDouble>(xyz_list);
	delete gauss;
	return Jelem;
}
