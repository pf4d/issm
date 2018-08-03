#include "./FreeSurfaceTopAnalysis.h"
#include "../toolkits/toolkits.h"
#include "../classes/classes.h"
#include "../shared/shared.h"
#include "../modules/modules.h"

/*Model processing*/
void FreeSurfaceTopAnalysis::CreateConstraints(Constraints* constraints,IoModel* iomodel){/*{{{*/
}/*}}}*/
void FreeSurfaceTopAnalysis::CreateLoads(Loads* loads, IoModel* iomodel){/*{{{*/

	/*Intermediaries*/
	int penpair_ids[2];
	int count=0;
	int numvertex_pairing;

	/*Create Penpair for vertex_pairing: */
	IssmDouble *vertex_pairing=NULL;
	IssmDouble *nodeonsurface=NULL;
	iomodel->FetchData(&vertex_pairing,&numvertex_pairing,NULL,"md.masstransport.vertex_pairing");
	if(iomodel->domaintype!=Domain2DhorizontalEnum) iomodel->FetchData(&nodeonsurface,NULL,NULL,"md.mesh.vertexonsurface");
	for(int i=0;i<numvertex_pairing;i++){

		if(iomodel->my_vertices[reCast<int>(vertex_pairing[2*i+0])-1]){

			/*In debugging mode, check that the second node is in the same cpu*/
			_assert_(iomodel->my_vertices[reCast<int>(vertex_pairing[2*i+1])-1]);

			/*Skip if one of the two is not on the bed*/
			if(iomodel->domaintype!=Domain2DhorizontalEnum){
				if(!(reCast<bool>(nodeonsurface[reCast<int>(vertex_pairing[2*i+0])-1])) || !(reCast<bool>(nodeonsurface[reCast<int>(vertex_pairing[2*i+1])-1]))) continue;
			}

			/*Get node ids*/
			penpair_ids[0]=iomodel->nodecounter+reCast<int>(vertex_pairing[2*i+0]);
			penpair_ids[1]=iomodel->nodecounter+reCast<int>(vertex_pairing[2*i+1]);

			/*Create Load*/
			loads->AddObject(new Penpair(
							iomodel->loadcounter+count+1,
							&penpair_ids[0],
							FreeSurfaceTopAnalysisEnum));
			count++;
		}
	}

	/*free ressources: */
	iomodel->DeleteData(vertex_pairing,"md.masstransport.vertex_pairing");
	iomodel->DeleteData(nodeonsurface,"md.mesh.vertexonsurface");
}/*}}}*/
void FreeSurfaceTopAnalysis::CreateNodes(Nodes* nodes,IoModel* iomodel){/*{{{*/

	if(iomodel->domaintype!=Domain2DhorizontalEnum) iomodel->FetchData(2,"md.mesh.vertexonbase","md.mesh.vertexonsurface");
	::CreateNodes(nodes,iomodel,FreeSurfaceTopAnalysisEnum,P1Enum);
	iomodel->DeleteData(2,"md.mesh.vertexonbase","md.mesh.vertexonsurface");
}/*}}}*/
int  FreeSurfaceTopAnalysis::DofsPerNode(int** doflist,int domaintype,int approximation){/*{{{*/
	return 1;
}/*}}}*/
void FreeSurfaceTopAnalysis::UpdateElements(Elements* elements,IoModel* iomodel,int analysis_counter,int analysis_type){/*{{{*/

	/*Now, is the model 3d? otherwise, do nothing: */
	if (iomodel->domaintype==Domain2DhorizontalEnum)return;

	int smb_model;
	int finiteelement = P1Enum;

	/*Fetch data needed: */
	iomodel->FindConstant(&smb_model,"md.smb.model");

	/*Update elements: */
	int counter=0;
	for(int i=0;i<iomodel->numberofelements;i++){
		if(iomodel->my_elements[i]){
			Element* element=(Element*)elements->GetObjectByOffset(counter);
			element->Update(i,iomodel,analysis_counter,analysis_type,finiteelement);
			counter++;
		}
	}

	iomodel->FetchDataToInput(elements,"md.geometry.surface",SurfaceEnum);
	iomodel->FetchDataToInput(elements,"md.slr.sealevel",SealevelEnum,0);
	iomodel->FetchDataToInput(elements,"md.mask.ice_levelset",MaskIceLevelsetEnum);
	iomodel->FetchDataToInput(elements,"md.initialization.vx",VxEnum);
	if(iomodel->domaintype!=Domain2DhorizontalEnum){
		iomodel->FetchDataToInput(elements,"md.mesh.vertexonsurface",MeshVertexonsurfaceEnum);
		iomodel->FetchDataToInput(elements,"md.mesh.vertexonbase",MeshVertexonbaseEnum);
	}
	if(iomodel->domaindim==3){
		iomodel->FetchDataToInput(elements,"md.initialization.vz",VzEnum);
	}
	switch(smb_model){
		case SMBforcingEnum:
			iomodel->FetchDataToInput(elements,"md.smb.mass_balance",SmbMassBalanceEnum,0.);
			break;
		default:
			/*Nothing for now*/
			;
	}
}/*}}}*/
void FreeSurfaceTopAnalysis::UpdateParameters(Parameters* parameters,IoModel* iomodel,int solution_enum,int analysis_enum){/*{{{*/
}/*}}}*/

/*Finite Element Analysis*/
void           FreeSurfaceTopAnalysis::Core(FemModel* femmodel){/*{{{*/
	_error_("not implemented");
}/*}}}*/
ElementVector* FreeSurfaceTopAnalysis::CreateDVector(Element* element){/*{{{*/
	/*Default, return NULL*/
	return NULL;
}/*}}}*/
ElementMatrix* FreeSurfaceTopAnalysis::CreateJacobianMatrix(Element* element){/*{{{*/
_error_("Not implemented");
}/*}}}*/
ElementMatrix* FreeSurfaceTopAnalysis::CreateKMatrix(Element* element){/*{{{*/

	/*Intermediaries*/
	int         domaintype,dim,stabilization;
	Element*    topelement = NULL;
	IssmDouble *xyz_list  = NULL;
	IssmDouble  Jdet,D_scalar,dt,h;
	IssmDouble  vel,vx,vy;

	/*Get top element*/
	element->FindParam(&domaintype,DomainTypeEnum);
	switch(domaintype){
		case Domain2DhorizontalEnum:
			topelement = element;
			dim = 2;
			break;
		case Domain2DverticalEnum:
			if(!element->IsOnSurface()) return NULL;
			topelement = element->SpawnTopElement();
			dim = 1;
			break;
		case Domain3DEnum:
			if(!element->IsOnSurface()) return NULL;
			topelement = element->SpawnTopElement();
			dim = 2;
			break;
		default: _error_("mesh "<<EnumToStringx(domaintype)<<" not supported yet");
	}

	/*Fetch number of nodes and dof for this finite element*/
	int numnodes = topelement->GetNumberOfNodes();

	/*Initialize Element vector*/
	ElementMatrix* Ke     = topelement->NewElementMatrix(NoneApproximationEnum);
	IssmDouble*    basis  = xNew<IssmDouble>(numnodes);
	IssmDouble*    B      = xNew<IssmDouble>(dim*numnodes);
	IssmDouble*    Bprime = xNew<IssmDouble>(dim*numnodes);
	IssmDouble*    D      = xNew<IssmDouble>(dim*dim);

	/*Retrieve all inputs and parameters*/
	topelement->GetVerticesCoordinates(&xyz_list);
	topelement->FindParam(&dt,TimesteppingTimeStepEnum);
	topelement->FindParam(&stabilization,MasstransportStabilizationEnum);
	Input* vx_input=topelement->GetInput(VxEnum); _assert_(vx_input);
	Input* vy_input=NULL;
	if(dim>1){vy_input = topelement->GetInput(VyEnum); _assert_(vy_input);}
	h = topelement->CharacteristicLength();

	/* Start  looping on the number of gaussian points: */
	Gauss* gauss=topelement->NewGauss(2);
	for(int ig=gauss->begin();ig<gauss->end();ig++){
		gauss->GaussPoint(ig);

		topelement->JacobianDeterminant(&Jdet,xyz_list,gauss);
		topelement->NodalFunctions(basis,gauss);

		vx_input->GetInputValue(&vx,gauss);
		if(dim==2) vy_input->GetInputValue(&vy,gauss);

		D_scalar=gauss->weight*Jdet;
		TripleMultiply(basis,1,numnodes,1,
					&D_scalar,1,1,0,
					basis,1,numnodes,0,
					&Ke->values[0],1);

		GetB(B,topelement,dim,xyz_list,gauss);
		GetBprime(Bprime,topelement,dim,xyz_list,gauss);

		D_scalar=dt*gauss->weight*Jdet;
		for(int i=0;i<dim*dim;i++) D[i]=0.;
		D[0] = D_scalar*vx;
		if(dim==2) D[1*dim+1]=D_scalar*vy;

		TripleMultiply(B,dim,numnodes,1,
					D,dim,dim,0,
					Bprime,dim,numnodes,0,
					&Ke->values[0],1);

		if(stabilization==2){
			/*Streamline upwinding*/
			if(dim==1){
				vel=fabs(vx)+1.e-8;
				D[0] = h/(2.*vel)*vx*vx;
			}
			else{
				vel=sqrt(vx*vx+vy*vy)+1.e-8;
				D[0*dim+0]=h/(2*vel)*vx*vx;
				D[1*dim+0]=h/(2*vel)*vy*vx;
				D[0*dim+1]=h/(2*vel)*vx*vy;
				D[1*dim+1]=h/(2*vel)*vy*vy;
			}
		}
		else if(stabilization==1){
			/*SSA*/
			if(dim==1){
				vx_input->GetInputAverage(&vx);
				D[0]=h/2.*fabs(vx);
			}
			else{
				vx_input->GetInputAverage(&vx);
				vy_input->GetInputAverage(&vy);
				D[0*dim+0]=h/2.0*fabs(vx);
				D[1*dim+1]=h/2.0*fabs(vy);
			}
		}
		if(stabilization==1 || stabilization==2){
			for(int i=0;i<dim*dim;i++) D[i]=D_scalar*D[i];
			TripleMultiply(Bprime,dim,numnodes,1,
						D,dim,dim,0,
						Bprime,dim,numnodes,0,
						&Ke->values[0],1);
		}
	}

	/*Clean up and return*/
	xDelete<IssmDouble>(xyz_list);
	xDelete<IssmDouble>(basis);
	xDelete<IssmDouble>(B);
	xDelete<IssmDouble>(Bprime);
	xDelete<IssmDouble>(D);
	delete gauss;
	if(domaintype!=Domain2DhorizontalEnum){topelement->DeleteMaterials(); delete topelement;};
	return Ke;
}/*}}}*/
ElementVector* FreeSurfaceTopAnalysis::CreatePVector(Element* element){/*{{{*/
	/*Intermediaries*/
	int         domaintype,dim;
	IssmDouble  Jdet,dt;
	IssmDouble  ms,surface,vz;
	Element*    topelement = NULL;
	IssmDouble *xyz_list  = NULL;

	/*Get top element*/
	element->FindParam(&domaintype,DomainTypeEnum);
	switch(domaintype){
		case Domain2DhorizontalEnum:
			topelement = element;
			dim = 2;
			break;
		case Domain2DverticalEnum:
			if(!element->IsOnSurface()) return NULL;
			topelement = element->SpawnTopElement();
			dim = 1;
			break;
		case Domain3DEnum:
			if(!element->IsOnSurface()) return NULL;
			topelement = element->SpawnTopElement();
			dim = 2;
			break;
		default: _error_("mesh "<<EnumToStringx(domaintype)<<" not supported yet");
	}

	/*Fetch number of nodes and dof for this finite element*/
	int numnodes = topelement->GetNumberOfNodes();

	/*Initialize Element vector and other vectors*/
	ElementVector* pe    = topelement->NewElementVector();
	IssmDouble*    basis = xNew<IssmDouble>(numnodes);

	/*Retrieve all inputs and parameters*/
	topelement->GetVerticesCoordinates(&xyz_list);
	topelement->FindParam(&dt,TimesteppingTimeStepEnum);
	Input* ms_input      = topelement->GetInput(SmbMassBalanceEnum);  _assert_(ms_input);
	Input* surface_input = topelement->GetInput(SurfaceEnum);                     _assert_(surface_input);
	Input* vz_input      = NULL;
	switch(dim){
		case 1: vz_input = topelement->GetInput(VyEnum); _assert_(vz_input); break;
		case 2: vz_input = topelement->GetInput(VzEnum); _assert_(vz_input); break;
		default: _error_("not implemented");
	}

	/*Initialize mb_correction to 0, do not forget!:*/
	/* Start  looping on the number of gaussian points: */
	Gauss* gauss=topelement->NewGauss(2);
	for(int ig=gauss->begin();ig<gauss->end();ig++){
		gauss->GaussPoint(ig);

		topelement->JacobianDeterminant(&Jdet,xyz_list,gauss);
		topelement->NodalFunctions(basis,gauss);

		ms_input->GetInputValue(&ms,gauss);
		vz_input->GetInputValue(&vz,gauss);
		surface_input->GetInputValue(&surface,gauss);

		for(int i=0;i<numnodes;i++) pe->values[i]+=Jdet*gauss->weight*(surface + dt*ms + dt*vz)*basis[i];
	}

	/*Clean up and return*/
	xDelete<IssmDouble>(xyz_list);
	xDelete<IssmDouble>(basis);
	delete gauss;
	if(domaintype!=Domain2DhorizontalEnum){topelement->DeleteMaterials(); delete topelement;};
	return pe;

}/*}}}*/
void           FreeSurfaceTopAnalysis::GetB(IssmDouble* B,Element* element,int dim,IssmDouble* xyz_list,Gauss* gauss){/*{{{*/
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
		for(int j=0;j<dim;j++){
			B[numnodes*j+i] = basis[i];
		}
	}

	/*Clean-up*/
	xDelete<IssmDouble>(basis);
}/*}}}*/
void           FreeSurfaceTopAnalysis::GetBprime(IssmDouble* Bprime,Element* element,int dim,IssmDouble* xyz_list,Gauss* gauss){/*{{{*/
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
	for(int i=0;i<numnodes;i++){
		for(int j=0;j<dim;j++){
			Bprime[numnodes*j+i] = dbasis[j*numnodes+i];
		}
	}

	/*Clean-up*/
	xDelete<IssmDouble>(dbasis);

}/*}}}*/
void           FreeSurfaceTopAnalysis::GetSolutionFromInputs(Vector<IssmDouble>* solution,Element* element){/*{{{*/
	   _error_("not implemented yet");
}/*}}}*/
void           FreeSurfaceTopAnalysis::GradientJ(Vector<IssmDouble>* gradient,Element* element,int control_type,int control_index){/*{{{*/
	_error_("Not implemented yet");
}/*}}}*/
void           FreeSurfaceTopAnalysis::InputUpdateFromSolution(IssmDouble* solution,Element* element){/*{{{*/

	element->InputUpdateFromSolutionOneDof(solution,SurfaceEnum);
}/*}}}*/
void           FreeSurfaceTopAnalysis::UpdateConstraints(FemModel* femmodel){/*{{{*/
	/*Default, do nothing*/
	return;
}/*}}}*/
