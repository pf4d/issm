#include "./HydrologyShreveAnalysis.h"
#include "../toolkits/toolkits.h"
#include "../classes/classes.h"
#include "../shared/shared.h"
#include "../modules/modules.h"

/*Model processing*/
void HydrologyShreveAnalysis::CreateConstraints(Constraints* constraints,IoModel* iomodel){/*{{{*/

	/*retrieve some parameters: */
	int          hydrology_model;
	iomodel->FindConstant(&hydrology_model,"md.hydrology.model");

	if(hydrology_model!=HydrologyshreveEnum) return;

	IoModelToConstraintsx(constraints,iomodel,"md.hydrologyshreve.spcwatercolumn",HydrologyShreveAnalysisEnum,P1Enum);

}/*}}}*/
void HydrologyShreveAnalysis::CreateLoads(Loads* loads, IoModel* iomodel){/*{{{*/
	/*No loads*/
}/*}}}*/
void HydrologyShreveAnalysis::CreateNodes(Nodes* nodes,IoModel* iomodel){/*{{{*/

	/*Fetch parameters: */
	int  hydrology_model;
	iomodel->FindConstant(&hydrology_model,"md.hydrology.model");

	/*Now, do we really want Shreve?*/
	if(hydrology_model!=HydrologyshreveEnum) return;

	if(iomodel->domaintype==Domain3DEnum) iomodel->FetchData(2,"md.mesh.vertexonbase","md.mesh.vertexonsurface");
	::CreateNodes(nodes,iomodel,HydrologyShreveAnalysisEnum,P1Enum);
	iomodel->DeleteData(2,"md.mesh.vertexonbase","md.mesh.vertexonsurface");
}/*}}}*/
int  HydrologyShreveAnalysis::DofsPerNode(int** doflist,int domaintype,int approximation){/*{{{*/
	return 1;
}/*}}}*/
void HydrologyShreveAnalysis::UpdateElements(Elements* elements,IoModel* iomodel,int analysis_counter,int analysis_type){/*{{{*/

	/*Fetch data needed: */
	int    hydrology_model;
	iomodel->FindConstant(&hydrology_model,"md.hydrology.model");

	/*Now, do we really want Shreve?*/
	if(hydrology_model!=HydrologyshreveEnum) return;

	/*Update elements: */
	int counter=0;
	for(int i=0;i<iomodel->numberofelements;i++){
		if(iomodel->my_elements[i]){
			Element* element=(Element*)elements->GetObjectByOffset(counter);
			element->Update(i,iomodel,analysis_counter,analysis_type,P1Enum);
			counter++;
		}
	}

	iomodel->FetchDataToInput(elements,"md.geometry.thickness",ThicknessEnum);
	iomodel->FetchDataToInput(elements,"md.geometry.surface",SurfaceEnum);
	iomodel->FetchDataToInput(elements,"md.geometry.base",BaseEnum);
	iomodel->FetchDataToInput(elements,"md.slr.sealevel",SealevelEnum,0);
	if(iomodel->domaintype!=Domain2DhorizontalEnum){
		iomodel->FetchDataToInput(elements,"md.mesh.vertexonbase",MeshVertexonbaseEnum);
		iomodel->FetchDataToInput(elements,"md.mesh.vertexonsurface",MeshVertexonsurfaceEnum);
	}
	iomodel->FetchDataToInput(elements,"md.mask.ice_levelset",MaskIceLevelsetEnum);
	iomodel->FetchDataToInput(elements,"md.mask.groundedice_levelset",MaskGroundediceLevelsetEnum);
	iomodel->FetchDataToInput(elements,"md.basalforcings.groundedice_melting_rate",BasalforcingsGroundediceMeltingRateEnum);
	iomodel->FetchDataToInput(elements,"md.initialization.watercolumn",WatercolumnEnum);

	elements->InputDuplicate(WatercolumnEnum,WaterColumnOldEnum);
}/*}}}*/
void HydrologyShreveAnalysis::UpdateParameters(Parameters* parameters,IoModel* iomodel,int solution_enum,int analysis_enum){/*{{{*/

	/*retrieve some parameters: */
	int  hydrology_model;
	iomodel->FindConstant(&hydrology_model,"md.hydrology.model");

	/*Now, do we really want Shreve?*/
	if(hydrology_model!=HydrologyshreveEnum) return;

	parameters->AddObject(new IntParam(HydrologyModelEnum,hydrology_model));
	parameters->AddObject(iomodel->CopyConstantObject("md.hydrology.stabilization",HydrologyshreveStabilizationEnum));

}/*}}}*/

/*Finite Element Analysis*/
void           HydrologyShreveAnalysis::Core(FemModel* femmodel){/*{{{*/
	_error_("not implemented");
}/*}}}*/
ElementVector* HydrologyShreveAnalysis::CreateDVector(Element* element){/*{{{*/
	/*Default, return NULL*/
	return NULL;
}/*}}}*/
void           HydrologyShreveAnalysis::CreateHydrologyWaterVelocityInput(Element* element){/*{{{*/

	/*Intermediaries*/
	IssmDouble dsdx,dsdy,dbdx,dbdy,w;

	/*Retrieve all inputs and parameters*/
	IssmDouble  rho_ice   = element->GetMaterialParameter(MaterialsRhoIceEnum);
	IssmDouble  rho_water = element->GetMaterialParameter(MaterialsRhoSeawaterEnum);
	IssmDouble  g         = element->GetMaterialParameter(ConstantsGEnum);
	IssmDouble  mu_water  = element->GetMaterialParameter(MaterialsMuWaterEnum);
	Input* surfaceslopex_input = element->GetInput(SurfaceSlopeXEnum); _assert_(surfaceslopex_input);
	Input* surfaceslopey_input = element->GetInput(SurfaceSlopeYEnum); _assert_(surfaceslopey_input);
	Input* bedslopex_input     = element->GetInput(BedSlopeXEnum);     _assert_(bedslopex_input);
	Input* bedslopey_input     = element->GetInput(BedSlopeYEnum);     _assert_(bedslopey_input);
	Input* watercolumn_input   = element->GetInput(WatercolumnEnum);   _assert_(watercolumn_input);

	/*Fetch number of vertices and allocate output*/
	int numvertices = element->GetNumberOfVertices();
	IssmDouble* vx  = xNew<IssmDouble>(numvertices);
	IssmDouble* vy  = xNew<IssmDouble>(numvertices);

	Gauss* gauss=element->NewGauss();
	for(int iv=0;iv<numvertices;iv++){
		gauss->GaussVertex(iv);
		surfaceslopex_input->GetInputValue(&dsdx,gauss);
		surfaceslopey_input->GetInputValue(&dsdy,gauss);
		bedslopex_input->GetInputValue(&dbdx,gauss);
		bedslopey_input->GetInputValue(&dbdy,gauss);
		watercolumn_input->GetInputValue(&w,gauss);

		/* Water velocity x and y components */
		vx[iv]= - w*w/(12 * mu_water)*(rho_ice*g*dsdx+(rho_water-rho_ice)*g*dbdx);
		vy[iv]= - w*w/(12 * mu_water)*(rho_ice*g*dsdy+(rho_water-rho_ice)*g*dbdy);
	}

	/*clean-up*/
	delete gauss;

	/*Add to inputs*/
	element->AddInput(HydrologyWaterVxEnum,vx,P1Enum);
	element->AddInput(HydrologyWaterVyEnum,vy,P1Enum);
	xDelete<IssmDouble>(vx);
	xDelete<IssmDouble>(vy);
}/*}}}*/
ElementMatrix* HydrologyShreveAnalysis::CreateJacobianMatrix(Element* element){/*{{{*/
_error_("Not implemented");
}/*}}}*/
ElementMatrix* HydrologyShreveAnalysis::CreateKMatrix(Element* element){/*{{{*/

	/*Intermediaries */
	IssmDouble diffusivity;
	IssmDouble Jdet,D_scalar,dt,h;
	IssmDouble vx,vy,vel,dvxdx,dvydy;
	IssmDouble dvx[2],dvy[2];
	IssmDouble* xyz_list = NULL;

	/*Fetch number of nodes and dof for this finite element*/
	int numnodes = element->GetNumberOfNodes();

	/*Initialize Element vector and other vectors*/
	ElementMatrix* Ke     = element->NewElementMatrix();
	IssmDouble*    basis  = xNew<IssmDouble>(numnodes);
	IssmDouble*    B      = xNew<IssmDouble>(2*numnodes);
	IssmDouble*    Bprime = xNew<IssmDouble>(2*numnodes);
	IssmDouble     D[2][2]={0.};

	/*Create water velocity vx and vy from current inputs*/
	CreateHydrologyWaterVelocityInput(element);

	/*Retrieve all inputs and parameters*/
	element->GetVerticesCoordinates(&xyz_list);
	element->FindParam(&dt,TimesteppingTimeStepEnum);
	element->FindParam(&diffusivity,HydrologyshreveStabilizationEnum);
	Input* vx_input=element->GetInput(HydrologyWaterVxEnum); _assert_(vx_input);
	Input* vy_input=element->GetInput(HydrologyWaterVyEnum); _assert_(vy_input);
	h = element->CharacteristicLength();

	/* Start  looping on the number of gaussian points: */
	Gauss* gauss=element->NewGauss(2);
	for(int ig=gauss->begin();ig<gauss->end();ig++){
		gauss->GaussPoint(ig);

		element->JacobianDeterminant(&Jdet,xyz_list,gauss);
		element->NodalFunctions(basis,gauss);

		vx_input->GetInputValue(&vx,gauss);
		vy_input->GetInputValue(&vy,gauss);
		vx_input->GetInputDerivativeValue(&dvx[0],xyz_list,gauss);
		vy_input->GetInputDerivativeValue(&dvy[0],xyz_list,gauss);

		D_scalar=gauss->weight*Jdet;

		TripleMultiply(basis,1,numnodes,1,
					&D_scalar,1,1,0,
					basis,1,numnodes,0,
					Ke->values,1);

		GetB(B,element,xyz_list,gauss);
		GetBprime(Bprime,element,xyz_list,gauss);

		dvxdx=dvx[0];
		dvydy=dvy[1];
		D_scalar=dt*gauss->weight*Jdet;

		D[0][0]=D_scalar*dvxdx;
		D[1][1]=D_scalar*dvydy;
		TripleMultiply(B,2,numnodes,1,
					&D[0][0],2,2,0,
					B,2,numnodes,0,
					&Ke->values[0],1);

		D[0][0]=D_scalar*vx;
		D[1][1]=D_scalar*vy;
		TripleMultiply(B,2,numnodes,1,
					&D[0][0],2,2,0,
					Bprime,2,numnodes,0,
					&Ke->values[0],1);

		/*Artificial diffusivity*/
		vel=sqrt(vx*vx+vy*vy);
		D[0][0]=D_scalar*diffusivity*h/(2*vel)*vx*vx;
		D[1][0]=D_scalar*diffusivity*h/(2*vel)*vy*vx;
		D[0][1]=D_scalar*diffusivity*h/(2*vel)*vx*vy;
		D[1][1]=D_scalar*diffusivity*h/(2*vel)*vy*vy;
		TripleMultiply(Bprime,2,numnodes,1,
					&D[0][0],2,2,0,
					Bprime,2,numnodes,0,
					&Ke->values[0],1);
	}

	/*Clean up and return*/
	xDelete<IssmDouble>(xyz_list);
	xDelete<IssmDouble>(basis);
	xDelete<IssmDouble>(B);
	xDelete<IssmDouble>(Bprime);
	delete gauss;
	return Ke;
}/*}}}*/
ElementVector* HydrologyShreveAnalysis::CreatePVector(Element* element){/*{{{*/

	/*Skip if water or ice shelf element*/
	if(element->IsFloating()) return NULL;

	/*Intermediaries */
	IssmDouble  Jdet,dt;
	IssmDouble  mb,oldw;
	IssmDouble* xyz_list = NULL;

	/*Fetch number of nodes and dof for this finite element*/
	int numnodes = element->GetNumberOfNodes();

	/*Initialize Element vector and other vectors*/
	ElementVector* pe    = element->NewElementVector();
	IssmDouble*    basis = xNew<IssmDouble>(numnodes);

	/*Retrieve all inputs and parameters*/
	element->GetVerticesCoordinates(&xyz_list);
	element->FindParam(&dt,TimesteppingTimeStepEnum);
	Input* mb_input   = element->GetInput(BasalforcingsGroundediceMeltingRateEnum); _assert_(mb_input);
	Input* oldw_input = element->GetInput(WaterColumnOldEnum);                      _assert_(oldw_input);

	/*Initialize mb_correction to 0, do not forget!:*/
	/* Start  looping on the number of gaussian points: */
	Gauss* gauss=element->NewGauss(2);
	for(int ig=gauss->begin();ig<gauss->end();ig++){
		gauss->GaussPoint(ig);

		element->JacobianDeterminant(&Jdet,xyz_list,gauss);
		element->NodalFunctions(basis,gauss);

		mb_input->GetInputValue(&mb,gauss);
		oldw_input->GetInputValue(&oldw,gauss);

		if(dt!=0.){
			for(int i=0;i<numnodes;i++) pe->values[i]+=Jdet*gauss->weight*(oldw+dt*mb)*basis[i];
		}
		else{
			for(int i=0;i<numnodes;i++) pe->values[i]+=Jdet*gauss->weight*mb*basis[i];
		}
	}

	/*Clean up and return*/
	xDelete<IssmDouble>(xyz_list);
	xDelete<IssmDouble>(basis);
	delete gauss;
	return pe;
}/*}}}*/
void           HydrologyShreveAnalysis::GetB(IssmDouble* B,Element* element,IssmDouble* xyz_list,Gauss* gauss){/*{{{*/
	/*Compute B  matrix. B=[B1 B2 B3] where Bi is of size 3*NDOF2. 
	 * For node i, Bi can be expressed in the actual coordinate system
	 * by: 
	 *       Bi=[ N ]
	 *          [ N ]
	 * where N is the finiteelement function for node i.
	 *
	 * We assume B_prog has been allocated already, of size: 2x(NDOF1*numnodes)
	 */

	/*Fetch number of nodes for this finite element*/
	int numnodes = element->GetNumberOfNodes();

	/*Get nodal functions*/
	IssmDouble* basis=xNew<IssmDouble>(numnodes);
	element->NodalFunctions(basis,gauss);

	/*Build B: */
	for(int i=0;i<numnodes;i++){
		B[numnodes*0+i] = basis[i];
		B[numnodes*1+i] = basis[i];
	}

	/*Clean-up*/
	xDelete<IssmDouble>(basis);
}/*}}}*/
void           HydrologyShreveAnalysis::GetBprime(IssmDouble* Bprime,Element* element,IssmDouble* xyz_list,Gauss* gauss){/*{{{*/
	/*Compute B'  matrix. B'=[B1' B2' B3'] where Bi' is of size 3*NDOF2. 
	 * For node i, Bi' can be expressed in the actual coordinate system
	 * by: 
	 *       Bi_prime=[ dN/dx ]
	 *                [ dN/dy ]
	 * where N is the finiteelement function for node i.
	 *
	 * We assume B' has been allocated already, of size: 3x(NDOF2*numnodes)
	 */

	/*Fetch number of nodes for this finite element*/
	int numnodes = element->GetNumberOfNodes();

	/*Get nodal functions derivatives*/
	IssmDouble* dbasis=xNew<IssmDouble>(2*numnodes);
	element->NodalFunctionsDerivatives(dbasis,xyz_list,gauss);

	/*Build B': */
	for(int i=0;i<numnodes;i++){
		Bprime[numnodes*0+i] = dbasis[0*numnodes+i];
		Bprime[numnodes*1+i] = dbasis[1*numnodes+i];
	}

	/*Clean-up*/
	xDelete<IssmDouble>(dbasis);

}/*}}}*/
void           HydrologyShreveAnalysis::GetSolutionFromInputs(Vector<IssmDouble>* solution,Element* element){/*{{{*/
	element->GetSolutionFromInputsOneDof(solution,WatercolumnEnum);
}/*}}}*/
void           HydrologyShreveAnalysis::GradientJ(Vector<IssmDouble>* gradient,Element* element,int control_type,int control_index){/*{{{*/
	_error_("Not implemented yet");
}/*}}}*/
void           HydrologyShreveAnalysis::InputUpdateFromSolution(IssmDouble* solution,Element* element){/*{{{*/

	/*Intermediary*/
	int* doflist = NULL;

	/*Fetch number of nodes for this finite element*/
	int numnodes = element->GetNumberOfNodes();

	/*Fetch dof list and allocate solution vector*/
	element->GetDofList(&doflist,NoneApproximationEnum,GsetEnum);
	IssmDouble* values = xNew<IssmDouble>(numnodes);

	/*Use the dof list to index into the solution vector: */
	for(int i=0;i<numnodes;i++){
		values[i]=solution[doflist[i]];
		if(xIsNan<IssmDouble>(values[i])) _error_("NaN found in solution vector");
		if(xIsInf<IssmDouble>(values[i])) _error_("Inf found in solution vector");
		if (values[i]<10e-10) values[i]=10e-10; //correcting the water column to positive values
	}

	/*Add input to the element: */
	element->AddInput(WatercolumnEnum,values,element->GetElementType());

	/*Free ressources:*/
	xDelete<IssmDouble>(values);
	xDelete<int>(doflist);
}/*}}}*/
void           HydrologyShreveAnalysis::UpdateConstraints(FemModel* femmodel){/*{{{*/
	/*Default, do nothing*/
	return;
}/*}}}*/


/*Needed changes to switch to the Johnson formulation*//*{{{*/
/*All the changes are to be done in the velocity computation.
	The new velocity needs some new parameter that should be introduce in the hydrologyshreve class:
	'p' and 'q' which are the exponent of the Manning formula for laminar (p=2,q=1) or turbulent (p=2/3,q=1/2) flow
	'R' the hydraulic radius
	'n' the manning roughness coeficient

	With these, the velocity reads ;

	v= - (1/n)* pow(R,p)*pow((grad phi(rho_water*g)),q)

	you should also redefine the water pressure potential 'phi' with respect to the effective pressure deffinition given in Johson:
	phi=(rho_ice*g*( surface + ((rho_water/rho_ice)-1)*base - k_n*((thickness* grad(base))/omega) )

	where 
	'omega' is the fractional area of the base occupied by the water film
	'k_n' is a parameter
	This last equation derives from the effective pressure definition developped in Alley 1989
*/
/*}}}*/
