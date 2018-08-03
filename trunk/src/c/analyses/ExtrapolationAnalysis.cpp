#include "./ExtrapolationAnalysis.h"
#include "../toolkits/toolkits.h"
#include "../classes/classes.h"
#include "../shared/shared.h"
#include "../modules/modules.h"
#include "../solutionsequences/solutionsequences.h"

void ExtrapolationAnalysis::CreateConstraints(Constraints* constraints,IoModel* iomodel){/*{{{*/
	// do nothing for now
	return;
}
/*}}}*/
void ExtrapolationAnalysis::CreateLoads(Loads* loads, IoModel* iomodel){/*{{{*/
	//	do nothing for now
	return;
}/*}}}*/
void ExtrapolationAnalysis::CreateNodes(Nodes* nodes,IoModel* iomodel){/*{{{*/
	int finiteelement=P1Enum;
	if(iomodel->domaintype!=Domain2DhorizontalEnum) iomodel->FetchData(2,"md.mesh.vertexonbase","md.mesh.vertexonsurface");
	::CreateNodes(nodes,iomodel,ExtrapolationAnalysisEnum,finiteelement);
	iomodel->DeleteData(2,"md.mesh.vertexonbase","md.mesh.vertexonsurface");
}
/*}}}*/
int  ExtrapolationAnalysis::DofsPerNode(int** doflist,int domaintype,int approximation){/*{{{*/
	return 1;
}
/*}}}*/
void ExtrapolationAnalysis::UpdateElements(Elements* elements,IoModel* iomodel,int analysis_counter,int analysis_type){/*{{{*/
	int    finiteelement;

	/*Finite element type*/
	finiteelement = P1Enum;

	/*Update elements: */
	int counter=0;
	for(int i=0;i<iomodel->numberofelements;i++){
		if(iomodel->my_elements[i]){
			Element* element=(Element*)elements->GetObjectByOffset(counter);
			element->Update(i,iomodel,analysis_counter,analysis_type,finiteelement);
			counter++;
		}
	}
	if(iomodel->domaintype!=Domain2DhorizontalEnum){
		iomodel->FetchDataToInput(elements,"md.mesh.vertexonbase",MeshVertexonbaseEnum);
		iomodel->FetchDataToInput(elements,"md.mesh.vertexonsurface",MeshVertexonsurfaceEnum);
	}
}
/*}}}*/
void ExtrapolationAnalysis::UpdateParameters(Parameters* parameters,IoModel* iomodel,int solution_enum,int analysis_enum){/*{{{*/
	//do nothing for now
	return;
}
/*}}}*/

/*Finite element Analysis*/
void           ExtrapolationAnalysis::Core(FemModel* femmodel){/*{{{*/

	/* Intermediaries */
	bool save_results;
	int extvar_enum; 
	femmodel->parameters->FindParam(&save_results,SaveResultsEnum);
	femmodel->parameters->FindParam(&extvar_enum, ExtrapolationVariableEnum);

	/*activate formulation: */
	femmodel->SetCurrentConfiguration(ExtrapolationAnalysisEnum);

	if(VerboseSolution()) _printf0_("extrapolation of " << EnumToStringx(extvar_enum) << ": call computational core:\n");
	solutionsequence_linear(femmodel);

	if(save_results){
		if(VerboseSolution()) _printf0_("   saving results\n");
		femmodel->RequestedOutputsx(&femmodel->results,&extvar_enum,1);
	}
}/*}}}*/
ElementVector* ExtrapolationAnalysis::CreateDVector(Element* element){/*{{{*/
	/*Default, return NULL*/
	return NULL;
}/*}}}*/
ElementMatrix* ExtrapolationAnalysis::CreateJacobianMatrix(Element* element){/*{{{*/
	/* Jacobian required for the Newton solver */
	_error_("not implemented yet");
}/*}}}*/
ElementMatrix* ExtrapolationAnalysis::CreateKMatrix(Element* element){/*{{{*/

	/*Intermediaries */
	bool	      extrapolatebydiffusion = true;
	int	      dim, domaintype, extrapolationcase;
	int	      i,row,col,stabilization;
	IssmDouble  Jdet,D_scalar,h;
	IssmDouble  norm_dlsf;
	IssmDouble  hx,hy,hz,kappa;
	IssmDouble *xyz_list    = NULL;
	Element    *workelement = NULL;

	/*Get problem case*/
	extrapolationcase=GetExtrapolationCase(element);
	switch(extrapolationcase){
		case 0:
			if(!element->IsOnBase()) return NULL; 
			workelement = element->SpawnBasalElement(); 
			break;
		case 1: case 2: case 3: workelement=element; break;
	}

	/* get extrapolation dimension */
	workelement->FindParam(&domaintype,DomainTypeEnum);
	switch(domaintype){
		case Domain2DverticalEnum:   dim=1; break;
		case Domain2DhorizontalEnum: dim=2; break;
		case Domain3DEnum:           dim=3; break;
      default: _error_("not supported yet");
	}

	/*Fetch number of nodes and dof for this finite element*/
	int numnodes = workelement->GetNumberOfNodes();

	/*Initialize Element vector and other vectors*/
	ElementMatrix* Ke     = workelement->NewElementMatrix();
	IssmDouble*    B      = xNew<IssmDouble>(dim*numnodes);
	IssmDouble*    Bprime = xNew<IssmDouble>(dim*numnodes);
	IssmDouble*		dlsf   = xNew<IssmDouble>(dim);
	IssmDouble*		normal = xNew<IssmDouble>(dim);
   IssmDouble*		D	    = xNewZeroInit<IssmDouble>(dim*dim);

	/*Retrieve all inputs and parameters*/
	Input* lsf_slopex_input=workelement->GetInput(LevelsetfunctionSlopeXEnum); _assert_(lsf_slopex_input);
	Input* lsf_slopey_input=workelement->GetInput(LevelsetfunctionSlopeYEnum); _assert_(lsf_slopey_input);
	workelement->GetVerticesCoordinates(&xyz_list);

	/* In 3d we are going to extrude using horizontal diffusion. Since the layers are 3d,
	 * we change the geometry of the element to make it flat. That way, the diffusion is
	 * handled along the layers
	 */
	if(element->ObjectEnum()==PentaEnum){
		for(i=0;i<3;i++) xyz_list[3*i+2] = 0.;
		for(i=3;i<6;i++) xyz_list[3*i+2] = 1.;
	}

	/* Start  looping on the number of gaussian points: */
	Gauss* gauss=workelement->NewGauss(2);
	for(int ig=gauss->begin();ig<gauss->end();ig++){
		gauss->GaussPoint(ig);

		workelement->JacobianDeterminant(&Jdet,xyz_list,gauss);
		GetB(B,workelement,xyz_list,gauss,dim);
		GetBprime(Bprime,workelement,xyz_list,gauss,dim);
		
		D_scalar=gauss->weight*Jdet;

		if(extrapolatebydiffusion){

			/* diffuse values outward only along the xy-plane*/
         for(int i=0;i<2;i++) D[i*dim+i] = D_scalar;

			TripleMultiply(Bprime,dim,numnodes,1,
					D,dim,dim,0,
					Bprime,dim,numnodes,0,
					&Ke->values[0],1);
		}
		else{
			/* extrapolate values along normal */
			/* Get normal on ice boundary */
			lsf_slopex_input->GetInputValue(&dlsf[0],gauss);
			if(dim>1)
				lsf_slopey_input->GetInputValue(&dlsf[1],gauss);
			if(dim>2)
				dlsf[2]=0.;
			norm_dlsf=0.;
			for(i=0;i<dim;i++)	norm_dlsf+=dlsf[i]*dlsf[i]; 
			norm_dlsf=sqrt(norm_dlsf); _assert_(norm_dlsf>0.);

			if(norm_dlsf>0.)
				for(i=0;i<dim;i++)	normal[i]=dlsf[i]/norm_dlsf;
			else
				for(i=0;i<dim;i++)	normal[i]=0.;
			
			for(row=0;row<dim;row++)
				for(col=0;col<dim;col++)
					if(row==col)
						D[row*dim+col]=D_scalar*normal[row];
					else
						D[row*dim+col]=0.;
			TripleMultiply(B,dim,numnodes,1,
						D,dim,dim,0,
						Bprime,dim,numnodes,0,
						&Ke->values[0],1);

			/* stabilization */
			/* do not use streamline upwinding for extrapolation: it yields oscillating results due to diffusion along normal direction, but none across */
			stabilization=1;
			if (stabilization==0){/* no stabilization, do nothing*/}
			else if(stabilization==1){
				/* Artificial Diffusion */
				workelement->ElementSizes(&hx,&hy,&hz);
				h=sqrt(pow(hx*normal[0],2) + pow(hy*normal[1],2));
				kappa=h/2.+1.e-14; 
				for(row=0;row<dim;row++)
					for(col=0;col<dim;col++)
						if(row==col)
							D[row*dim+col]=D_scalar*kappa;
						else
							D[row*dim+col]=0.;
				TripleMultiply(Bprime,dim,numnodes,1,
							D,dim,dim,0,
							Bprime,dim,numnodes,0,
							&Ke->values[0],1);
			}
		}
	}

	/*Clean up and return*/
	xDelete<IssmDouble>(xyz_list);
	xDelete<IssmDouble>(B);
	xDelete<IssmDouble>(Bprime);
	xDelete<IssmDouble>(D);
	xDelete<IssmDouble>(dlsf);
	xDelete<IssmDouble>(normal);
	delete gauss;
	if(extrapolationcase==0){workelement->DeleteMaterials(); delete workelement;};
	return Ke;

}/*}}}*/
ElementVector* ExtrapolationAnalysis::CreatePVector(Element* element){/*{{{*/
	return NULL;

}/*}}}*/
void           ExtrapolationAnalysis::GetB(IssmDouble* B,Element* element,IssmDouble* xyz_list,Gauss* gauss, int dim){/*{{{*/
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
	for(int i=0;i<numnodes;i++)
		for(int d=0;d<dim;d++)
			B[numnodes*d+i] = basis[i];

	/*Clean-up*/
	xDelete<IssmDouble>(basis);
}/*}}}*/
void           ExtrapolationAnalysis::GetBprime(IssmDouble* Bprime,Element* element,IssmDouble* xyz_list,Gauss* gauss, int dim){/*{{{*/
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
	IssmDouble* dbasis=xNew<IssmDouble>(dim*numnodes);
	element->NodalFunctionsDerivatives(dbasis,xyz_list,gauss);

	/*Build B': */
	for(int i=0;i<numnodes;i++)
		for(int d=0;d<dim;d++)
			Bprime[numnodes*d+i] = dbasis[d*numnodes+i];

	/*Clean-up*/
	xDelete<IssmDouble>(dbasis);

}/*}}}*/
void           ExtrapolationAnalysis::GetSolutionFromInputs(Vector<IssmDouble>* solution,Element* element){/*{{{*/
	_error_("not implemented yet");
}/*}}}*/
void           ExtrapolationAnalysis::GradientJ(Vector<IssmDouble>* gradient,Element* element,int control_type,int control_index){/*{{{*/
	_error_("Not implemented yet");
}/*}}}*/
void           ExtrapolationAnalysis::InputUpdateFromSolution(IssmDouble* solution,Element* element){/*{{{*/

	int extrapolationvariable, extrapolationcase;
	extrapolationcase=GetExtrapolationCase(element);
	element->FindParam(&extrapolationvariable, ExtrapolationVariableEnum);
	switch(extrapolationcase){
		case 0:
			element->InputUpdateFromSolutionOneDof(solution,extrapolationvariable);
			break;
		case 1:
			element->InputUpdateFromSolutionOneDof(solution,extrapolationvariable);
			break;
		case 2:
			element->InputUpdateFromSolutionOneDofCollapsed(solution,extrapolationvariable);
			break;
		case 3:
			element->InputUpdateFromSolutionOneDof(solution,extrapolationvariable);
			break;
	}
}/*}}}*/
int				ExtrapolationAnalysis::GetExtrapolationCase(Element* element){/*{{{*/

	/* Get case of extrapolation, depending on domain quality, and extrapolation variable */
	int domaintype, extrapolationvariable;
	int extrapolationcase;

	element->FindParam(&domaintype,DomainTypeEnum);
	switch(domaintype){
		case Domain2DverticalEnum:   extrapolationcase=0; break;
		case Domain2DhorizontalEnum: extrapolationcase=1; break;
		case Domain3DEnum: 
			element->FindParam(&extrapolationvariable, ExtrapolationVariableEnum);
			if(extrapolationvariable==ThicknessEnum){
            extrapolationcase=2; // scalar fields that are constant along z-axis
         }
			else{
            extrapolationcase=3; // scalar fields that vary along z-axis
         }
			break;
	}
	return extrapolationcase;
}/*}}}*/
void           ExtrapolationAnalysis::SetConstraintsOnIce(Element* element){/*{{{*/

	int numnodes=element->GetNumberOfNodes();	

	/* Intermediaries */
	int extvar_enum;
	IssmDouble phi,value;
	Node* node = NULL;

	/* Get parameters */
	element->FindParam(&extvar_enum, ExtrapolationVariableEnum);
	
	Input* levelset_input=element->GetInput(MaskIceLevelsetEnum); _assert_(levelset_input);
	Input* extvar_input=element->GetInput(extvar_enum); _assert_(extvar_input);

	Gauss* gauss=element->NewGauss();
	for(int in=0;in<numnodes;in++){
		gauss->GaussNode(element->GetElementType(),in);
		node=element->GetNode(in);
		levelset_input->GetInputValue(&phi,gauss);
		if(node->IsActive()){
			if(phi<=0.){
				/* if ice, set dirichlet BC */
				extvar_input->GetInputValue(&value,gauss);
				node->ApplyConstraint(0,value);
			}
			else {
				/* no ice, set no spc */
				node->DofInFSet(0); 
			}
		}
	}
	delete gauss;
}/*}}}*/
void           ExtrapolationAnalysis::UpdateConstraints(FemModel* femmodel){/*{{{*/

	for(int i=0;i<femmodel->elements->Size();i++){
		Element* element=xDynamicCast<Element*>(femmodel->elements->GetObjectByOffset(i));
		this->SetConstraintsOnIce(element);
	}
}/*}}}*/
