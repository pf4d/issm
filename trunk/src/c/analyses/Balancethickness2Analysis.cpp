#include "./Balancethickness2Analysis.h"
#include "../toolkits/toolkits.h"
#include "../classes/classes.h"
#include "../shared/shared.h"
#include "../modules/modules.h"

/*Model processing*/
void Balancethickness2Analysis::CreateConstraints(Constraints* constraints,IoModel* iomodel){/*{{{*/


	int finiteelement = P1Enum;
	IoModelToConstraintsx(constraints,iomodel,"md.balancethickness.spcthickness",Balancethickness2AnalysisEnum,finiteelement);

}/*}}}*/
void Balancethickness2Analysis::CreateLoads(Loads* loads, IoModel* iomodel){/*{{{*/

}/*}}}*/
void Balancethickness2Analysis::CreateNodes(Nodes* nodes,IoModel* iomodel){/*{{{*/

	int finiteelement = P1Enum;
	::CreateNodes(nodes,iomodel,Balancethickness2AnalysisEnum,finiteelement);
}/*}}}*/
int  Balancethickness2Analysis::DofsPerNode(int** doflist,int domaintype,int approximation){/*{{{*/
	return 1;
}/*}}}*/
void Balancethickness2Analysis::UpdateElements(Elements* elements,IoModel* iomodel,int analysis_counter,int analysis_type){/*{{{*/

	/*Finite element type*/
	int finiteelement = P1Enum;

	/*Load variables in element*/
	iomodel->FetchDataToInput(elements,"md.geometry.thickness",ThicknessEnum);
	iomodel->FetchDataToInput(elements,"md.geometry.surface",SurfaceEnum);
	iomodel->FetchDataToInput(elements,"md.geometry.base",BaseEnum);
	iomodel->FetchDataToInput(elements,"md.slr.sealevel",SealevelEnum,0);
	iomodel->FetchDataToInput(elements,"md.mask.ice_levelset",MaskIceLevelsetEnum);
	iomodel->FetchDataToInput(elements,"md.basalforcings.groundedice_melting_rate",BasalforcingsGroundediceMeltingRateEnum);
	iomodel->FetchDataToInput(elements,"md.smb.mass_balance",SmbMassBalanceEnum);
	iomodel->FetchDataToInput(elements,"md.balancethickness.thickening_rate",BalancethicknessThickeningRateEnum);
	iomodel->FetchDataToInput(elements,"md.balancethickness.omega",BalancethicknessOmegaEnum);

	/*Update elements: */
	int counter=0;
	for(int i=0;i<iomodel->numberofelements;i++){
		if(iomodel->my_elements[i]){
			Element* element=(Element*)elements->GetObjectByOffset(counter);
			element->Update(i,iomodel,analysis_counter,analysis_type,finiteelement);

			counter++;
		}
	}

}/*}}}*/
void Balancethickness2Analysis::UpdateParameters(Parameters* parameters,IoModel* iomodel,int solution_enum,int analysis_enum){/*{{{*/
}/*}}}*/

/*Finite Element Analysis*/
void           Balancethickness2Analysis::Core(FemModel* femmodel){/*{{{*/
	_error_("not implemented yet");
}/*}}}*/
ElementVector* Balancethickness2Analysis::CreateDVector(Element* element){/*{{{*/
	/*Default, return NULL*/
	return NULL;
}/*}}}*/
void           Balancethickness2Analysis::CreateD0(Element* element){/*{{{*/

	/*Intermediaries */
	IssmDouble       Gamma,h,mu0,ds[2],Cmu,B,k,s,b,normds;
	const int        n = 3;
	const IssmDouble Hstar = 500.;
	const IssmDouble Lstar = 500.e+3;
	IssmDouble *xyz_list  = NULL;

	/*Fetch number of vertices and allocate output*/
	int  numnodes = element->GetNumberOfNodes();
	IssmDouble* D0     = xNew<IssmDouble>(numnodes);

	/*retrieve what we need: */
	element->GetVerticesCoordinates(&xyz_list);
	Input* surfaceslopex_input = element->GetInput(SurfaceSlopeXEnum);
	Input* surfaceslopey_input = element->GetInput(SurfaceSlopeYEnum);
	Input* surface_input       = element->GetInput(SurfaceEnum);            _assert_(surface_input);
	Input* B_input             = element->GetInput(MaterialsRheologyBEnum); _assert_(B_input);
	IssmDouble rhog            = element->GetMaterialParameter(MaterialsRhoIceEnum)*element->GetMaterialParameter(ConstantsGEnum);

	/*Calculate damage evolution source term: */
	Gauss* gauss=element->NewGauss();
	for (int i=0;i<numnodes;i++){
		gauss->GaussNode(element->GetElementType(),i);
		
		B_input->GetInputValue(&B,gauss);
		if(surfaceslopex_input && surfaceslopey_input){
			surfaceslopex_input->GetInputValue(&ds[0],gauss);
			surfaceslopey_input->GetInputValue(&ds[1],gauss);
		}
		else{
			surface_input->GetInputDerivativeValue(&ds[0],xyz_list,gauss);
		}

		/*check slopes*/
		normds = sqrt(ds[0]*ds[0]+ds[1]*ds[1]);
		if (normds==0.){
			_error_("surface slope is zero");
		}
		if(normds<1.e-5){
			ds[0] = ds[0]/normds*1.e+5;
			ds[1] = ds[1]/normds*1.e+5;
			normds = 1.e-5;
		}

		mu0   = pow(2.,(1-3*n)/(2.*n)) * B;
		Gamma = pow(rhog,n) * pow(Hstar,2*(n+1)) * pow(Hstar/Lstar,2*(n+1)) * 1./pow(mu0,n);

		D0[i] = Gamma*pow(ds[0]*ds[0]+ds[1]*ds[1],(n-1)/2)/(n+2);
	}

	/*Add input*/
	element->AddInput(BalancethicknessD0Enum,D0,element->GetElementType());
	//if(surfaceslopex_input && surfaceslopey_input){
	//	element->DeleteInput(SurfaceSlopeXEnum);
	//	element->DeleteInput(SurfaceSlopeYEnum);
	//}
	
	/*Clean up and return*/
	xDelete<IssmDouble>(D0);
	xDelete<IssmDouble>(xyz_list);
	delete gauss;
}/*}}}*/
ElementMatrix* Balancethickness2Analysis::CreateJacobianMatrix(Element* element){/*{{{*/
_error_("Not implemented");
}/*}}}*/
ElementMatrix* Balancethickness2Analysis::CreateKMatrix(Element* element){/*{{{*/

	/*Intermediaries */
	IssmDouble  Jdet,D0,omega;
	IssmDouble* xyz_list = NULL;

	/*Fetch number of nodes and dof for this finite element*/
	int numnodes = element->GetNumberOfNodes();

	/*Initialize Element vector and other vectors*/
	ElementMatrix* Ke     = element->NewElementMatrix();
	IssmDouble*    dbasis = xNew<IssmDouble>(2*numnodes);

	/*Retrieve all inputs and parameters*/
	element->GetVerticesCoordinates(&xyz_list);
	Input* omega_input = element->GetInput(BalancethicknessOmegaEnum); _assert_(omega_input);
	Input* D0_input    = element->GetInput(BalancethicknessD0Enum);
	if(!D0_input){
		this->CreateD0(element);
		D0_input = element->GetInput(BalancethicknessD0Enum); _assert_(D0_input);
	}

	/* Start  looping on the number of gaussian points: */
	Gauss* gauss=element->NewGauss(2);
	for(int ig=gauss->begin();ig<gauss->end();ig++){
		gauss->GaussPoint(ig);
		element->NodalFunctionsDerivatives(dbasis,xyz_list,gauss);
		element->JacobianDeterminant(&Jdet,xyz_list,gauss);
		D0_input->GetInputValue(&D0,gauss);
		omega_input->GetInputValue(&omega,gauss);

		for(int i=0;i<numnodes;i++){
			for(int j=0;j<numnodes;j++){
				Ke->values[i*numnodes+j] += D0*exp(omega)*gauss->weight*Jdet*(dbasis[0*numnodes+i]*dbasis[0*numnodes+j] + dbasis[1*numnodes+i]*dbasis[1*numnodes+j]);
			}
		}
	}

	/*Clean up and return*/
	xDelete<IssmDouble>(xyz_list);
	xDelete<IssmDouble>(dbasis);
	delete gauss;
	return Ke;
}/*}}}*/
ElementVector* Balancethickness2Analysis::CreatePVector(Element* element){/*{{{*/

	/*Intermediaries */
	IssmDouble  dhdt,mb,ms,Jdet;
	IssmDouble* xyz_list = NULL;

	/*Fetch number of nodes and dof for this finite element*/
	int numnodes = element->GetNumberOfNodes();

	/*Initialize Element vector and other vectors*/
	ElementVector* pe    = element->NewElementVector();
	IssmDouble*    basis = xNew<IssmDouble>(numnodes);

	/*Retrieve all inputs and parameters*/
	element->GetVerticesCoordinates(&xyz_list);
	Input* ms_input   = element->GetInput(SmbMassBalanceEnum);                _assert_(ms_input);
	Input* mb_input   = element->GetInput(BasalforcingsGroundediceMeltingRateEnum);       _assert_(mb_input);
	Input* dhdt_input = element->GetInput(BalancethicknessThickeningRateEnum);            _assert_(dhdt_input);

	/* Start  looping on the number of gaussian points: */
	Gauss* gauss=element->NewGauss(2);
	for(int ig=gauss->begin();ig<gauss->end();ig++){
		gauss->GaussPoint(ig);

		element->JacobianDeterminant(&Jdet,xyz_list,gauss);
		element->NodalFunctions(basis,gauss);

		ms_input->GetInputValue(&ms,gauss);
		mb_input->GetInputValue(&mb,gauss);
		dhdt_input->GetInputValue(&dhdt,gauss);

		for(int i=0;i<numnodes;i++) pe->values[i]+=Jdet*gauss->weight*(
					(ms-mb-dhdt)*basis[i]
					);
	}

	/*Clean up and return*/
	xDelete<IssmDouble>(xyz_list);
	xDelete<IssmDouble>(basis);
	delete gauss;
	return pe;
}/*}}}*/
void           Balancethickness2Analysis::GetSolutionFromInputs(Vector<IssmDouble>* solution,Element* element){/*{{{*/
		element->GetSolutionFromInputsOneDof(solution,SurfaceEnum);
}/*}}}*/
void           Balancethickness2Analysis::GradientJ(Vector<IssmDouble>* gradient,Element* element,int control_type,int control_index){/*{{{*/
	_error_("Not implemented yet");
}/*}}}*/
void           Balancethickness2Analysis::InputUpdateFromSolution(IssmDouble* solution,Element* element){/*{{{*/

	/*Intermediaries*/
	IssmDouble  ds[2],s,b,D;
	IssmDouble* xyz_list = NULL;

	//element->InputUpdateFromSolutionOneDof(solution,ThicknessEnum);
	element->InputUpdateFromSolutionOneDof(solution,SurfaceEnum);

	/*Fetch number of vertices and allocate velocity vectors*/
	int numvertices = element->GetNumberOfVertices();
	IssmDouble* vel_list = xNew<IssmDouble>(numvertices);
	IssmDouble* vx_list  = xNew<IssmDouble>(numvertices);
	IssmDouble* vy_list  = xNew<IssmDouble>(numvertices);

	/*Retrieve all inputs and parameters*/
	element->GetVerticesCoordinates(&xyz_list);
	Input* D_input   = element->GetInput(BalancethicknessDiffusionCoefficientEnum);
	Input* H_input   = element->GetInput(ThicknessEnum);                            _assert_(H_input);
	Input* s_input   = element->GetInput(SurfaceEnum);                              _assert_(s_input);
	Input* b_input   = element->GetInput(BaseEnum);                                 _assert_(b_input);

	/*Calculate velocities*/
	Gauss* gauss=element->NewGauss();
	for(int iv=0;iv<numvertices;iv++){
		gauss->GaussVertex(iv);

		if(D_input){
			D_input->GetInputValue(&D,gauss);
		}
		else{
			D = 0.;
		}
		b_input->GetInputValue(&b,gauss);
		s_input->GetInputValue(&s,gauss);
		s_input->GetInputDerivativeValue(&ds[0],xyz_list,gauss);

		vx_list[iv] = -1./(s-b)*D*ds[0];
		vy_list[iv] = -1./(s-b)*D*ds[1];
		vel_list[iv] = sqrt(pow(vx_list[iv],2) + pow(vy_list[iv],2));
	}

	/*Add vx and vy as inputs to the tria element: */
	element->AddInput(VxEnum,vx_list,P1Enum);
	element->AddInput(VyEnum,vy_list,P1Enum);
	element->AddInput(VelEnum,vel_list,P1Enum);

	/*Free ressources:*/
	delete gauss;
	xDelete<IssmDouble>(vy_list);
	xDelete<IssmDouble>(vx_list);
	xDelete<IssmDouble>(vel_list);
	xDelete<IssmDouble>(xyz_list);
}/*}}}*/
void           Balancethickness2Analysis::UpdateConstraints(FemModel* femmodel){/*{{{*/
	/*Default, do nothing*/
	return;
}/*}}}*/
