#include "./FreeSurfaceBaseAnalysis.h"
#include "../toolkits/toolkits.h"
#include "../classes/classes.h"
#include "../shared/shared.h"
#include "../modules/modules.h"

/*Model processing*/
void FreeSurfaceBaseAnalysis::CreateConstraints(Constraints* constraints,IoModel* iomodel){/*{{{*/
}/*}}}*/
void FreeSurfaceBaseAnalysis::CreateLoads(Loads* loads, IoModel* iomodel){/*{{{*/

	/*Intermediaries*/
	int penpair_ids[2];
	int count=0;
	int numvertex_pairing;

	/*Create Penpair for vertex_pairing: */
	IssmDouble *vertex_pairing=NULL;
	IssmDouble *nodeonbase=NULL;
	iomodel->FetchData(&vertex_pairing,&numvertex_pairing,NULL,"md.masstransport.vertex_pairing");
	if(iomodel->domaintype!=Domain2DhorizontalEnum) iomodel->FetchData(&nodeonbase,NULL,NULL,"md.mesh.vertexonbase");
	for(int i=0;i<numvertex_pairing;i++){

		if(iomodel->my_vertices[reCast<int>(vertex_pairing[2*i+0])-1]){

			/*In debugging mode, check that the second node is in the same cpu*/
			_assert_(iomodel->my_vertices[reCast<int>(vertex_pairing[2*i+1])-1]);

			/*Skip if one of the two is not on the bed*/
			if(iomodel->domaintype!=Domain2DhorizontalEnum){
				if(!(reCast<bool>(nodeonbase[reCast<int>(vertex_pairing[2*i+0])-1])) || !(reCast<bool>(nodeonbase[reCast<int>(vertex_pairing[2*i+1])-1]))) continue;
			}

			/*Get node ids*/
			penpair_ids[0]=iomodel->nodecounter+reCast<int>(vertex_pairing[2*i+0]);
			penpair_ids[1]=iomodel->nodecounter+reCast<int>(vertex_pairing[2*i+1]);

			/*Create Load*/
			loads->AddObject(new Penpair(
							iomodel->loadcounter+count+1,
							&penpair_ids[0],
							FreeSurfaceBaseAnalysisEnum));
			count++;
		}
	}

	/*free ressources: */
	iomodel->DeleteData(vertex_pairing,"md.masstransport.vertex_pairing");
	iomodel->DeleteData(nodeonbase,"md.mesh.vertexonbase");
}/*}}}*/
void FreeSurfaceBaseAnalysis::CreateNodes(Nodes* nodes,IoModel* iomodel){/*{{{*/

	if(iomodel->domaintype!=Domain2DhorizontalEnum) iomodel->FetchData(2,"md.mesh.vertexonbase","md.mesh.vertexonsurface");
	::CreateNodes(nodes,iomodel,FreeSurfaceBaseAnalysisEnum,P1Enum);
	iomodel->DeleteData(2,"md.mesh.vertexonbase","md.mesh.vertexonsurface");
}/*}}}*/
int  FreeSurfaceBaseAnalysis::DofsPerNode(int** doflist,int domaintype,int approximation){/*{{{*/
	return 1;
}/*}}}*/
void FreeSurfaceBaseAnalysis::UpdateElements(Elements* elements,IoModel* iomodel,int analysis_counter,int analysis_type){/*{{{*/

	/*Now, is the model 3d? otherwise, do nothing: */
	if (iomodel->domaintype==Domain2DhorizontalEnum)return;

	/*Finite element type*/
	int finiteelement = P1Enum;

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
	iomodel->FetchDataToInput(elements,"md.basalforcings.groundedice_melting_rate",BasalforcingsGroundediceMeltingRateEnum);
	iomodel->FetchDataToInput(elements,"md.basalforcings.floatingice_melting_rate",BasalforcingsFloatingiceMeltingRateEnum);
	iomodel->FetchDataToInput(elements,"md.initialization.vx",VxEnum);
	iomodel->FetchDataToInput(elements,"md.initialization.vy",VyEnum);
	if(iomodel->domaindim==3){
		iomodel->FetchDataToInput(elements,"md.initialization.vz",VzEnum);
	}
	if(iomodel->domaintype!=Domain2DhorizontalEnum){
		iomodel->FetchDataToInput(elements,"md.mesh.vertexonbase",MeshVertexonbaseEnum);
		iomodel->FetchDataToInput(elements,"md.mesh.vertexonsurface",MeshVertexonsurfaceEnum);
	}
}/*}}}*/
void FreeSurfaceBaseAnalysis::UpdateParameters(Parameters* parameters,IoModel* iomodel,int solution_enum,int analysis_enum){/*{{{*/
}/*}}}*/

/*Finite Element Analysis*/
void           FreeSurfaceBaseAnalysis::Core(FemModel* femmodel){/*{{{*/
	_error_("not implemented");
}/*}}}*/
ElementVector* FreeSurfaceBaseAnalysis::CreateDVector(Element* element){/*{{{*/
	/*Default, return NULL*/
	return NULL;
}/*}}}*/
ElementMatrix* FreeSurfaceBaseAnalysis::CreateJacobianMatrix(Element* element){/*{{{*/
_error_("Not implemented");
}/*}}}*/
ElementMatrix* FreeSurfaceBaseAnalysis::CreateKMatrix(Element* element){/*{{{*/

	/*Intermediaries*/
	int         domaintype,dim,stabilization;
	Element*    basalelement = NULL;
	IssmDouble *xyz_list  = NULL;
	IssmDouble  Jdet,D_scalar,dt,h;
	IssmDouble  vel,vx,vy;

	/*Get basal element*/
	element->FindParam(&domaintype,DomainTypeEnum);
	switch(domaintype){
		case Domain2DhorizontalEnum:
			basalelement = element;
			dim = 2;
			break;
		case Domain2DverticalEnum:
			if(!element->IsOnBase()) return NULL;
			basalelement = element->SpawnBasalElement();
			dim = 1;
			break;
		case Domain3DEnum:
			if(!element->IsOnBase()) return NULL;
			basalelement = element->SpawnBasalElement();
			dim = 2;
			break;
		default: _error_("mesh "<<EnumToStringx(domaintype)<<" not supported yet");
	}

	/*Fetch number of nodes and dof for this finite element*/
	int numnodes = basalelement->GetNumberOfNodes();

	/*Initialize Element vector*/
	ElementMatrix* Ke     = basalelement->NewElementMatrix(NoneApproximationEnum);
	IssmDouble*    basis  = xNew<IssmDouble>(numnodes);
	IssmDouble*    B      = xNew<IssmDouble>(dim*numnodes);
	IssmDouble*    Bprime = xNew<IssmDouble>(dim*numnodes);
	IssmDouble*    D      = xNew<IssmDouble>(dim*dim);

	/*Retrieve all inputs and parameters*/
	basalelement->GetVerticesCoordinates(&xyz_list);
	basalelement->FindParam(&dt,TimesteppingTimeStepEnum);
	basalelement->FindParam(&stabilization,MasstransportStabilizationEnum);
	Input* vx_input=basalelement->GetInput(VxEnum); _assert_(vx_input);
	Input* vy_input=NULL;
	if(dim>1){vy_input = basalelement->GetInput(VyEnum); _assert_(vy_input);}
	h = basalelement->CharacteristicLength();

	/* Start  looping on the number of gaussian points: */
	Gauss* gauss=basalelement->NewGauss(2);
	for(int ig=gauss->begin();ig<gauss->end();ig++){
		gauss->GaussPoint(ig);

		basalelement->JacobianDeterminant(&Jdet,xyz_list,gauss);
		basalelement->NodalFunctions(basis,gauss);

		vx_input->GetInputValue(&vx,gauss);
		if(dim==2) vy_input->GetInputValue(&vy,gauss);

		D_scalar=gauss->weight*Jdet;
		TripleMultiply(basis,1,numnodes,1,
					&D_scalar,1,1,0,
					basis,1,numnodes,0,
					&Ke->values[0],1);

		GetB(B,basalelement,dim,xyz_list,gauss);
		GetBprime(Bprime,basalelement,dim,xyz_list,gauss);

		D_scalar=dt*gauss->weight*Jdet;
		for(int i=0;i<dim*dim;i++) D[i]=0.;
		D[0] = D_scalar*vx;
		if(dim==2) D[1*dim+1] = D_scalar*vy;

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
	if(domaintype!=Domain2DhorizontalEnum){basalelement->DeleteMaterials(); delete basalelement;};
	return Ke;
}/*}}}*/
ElementVector* FreeSurfaceBaseAnalysis::CreatePVector(Element* element){/*{{{*/
	/*Intermediaries*/
	int         domaintype,dim;
	IssmDouble  Jdet,dt;
	IssmDouble  gmb,fmb,mb,bed,phi,vz;
	Element*    basalelement = NULL;
	IssmDouble *xyz_list  = NULL;

	/*Get basal element*/
	element->FindParam(&domaintype,DomainTypeEnum);
	switch(domaintype){
		case Domain2DhorizontalEnum:
			basalelement = element;
			dim = 2;
			break;
		case Domain2DverticalEnum:
			if(!element->IsOnBase()) return NULL;
			basalelement = element->SpawnBasalElement();
			dim = 1;
			break;
		case Domain3DEnum:
			if(!element->IsOnBase()) return NULL;
			basalelement = element->SpawnBasalElement();
			dim = 2;
			break;
		default: _error_("mesh "<<EnumToStringx(domaintype)<<" not supported yet");
	}

	/*Fetch number of nodes and dof for this finite element*/
	int numnodes = basalelement->GetNumberOfNodes();

	/*Initialize Element vector and other vectors*/
	ElementVector* pe    = basalelement->NewElementVector();
	IssmDouble*    basis = xNew<IssmDouble>(numnodes);

	/*Retrieve all inputs and parameters*/
	basalelement->GetVerticesCoordinates(&xyz_list);
	basalelement->FindParam(&dt,TimesteppingTimeStepEnum);
	Input* groundedice_input   = basalelement->GetInput(MaskGroundediceLevelsetEnum);              _assert_(groundedice_input);
	Input* gmb_input           = basalelement->GetInput(BasalforcingsGroundediceMeltingRateEnum);  _assert_(gmb_input);
	Input* fmb_input           = basalelement->GetInput(BasalforcingsFloatingiceMeltingRateEnum);  _assert_(fmb_input);
	Input* base_input          = basalelement->GetInput(BaseEnum);                                 _assert_(base_input);
	Input* vz_input      = NULL;
	switch(dim){
		case 1: vz_input = basalelement->GetInput(VyEnum); _assert_(vz_input); break;
		case 2: vz_input = basalelement->GetInput(VzEnum); _assert_(vz_input); break;
		default: _error_("not implemented");
	}

	/* Start  looping on the number of gaussian points: */
	Gauss* gauss=basalelement->NewGauss(2);
	for(int ig=gauss->begin();ig<gauss->end();ig++){
		gauss->GaussPoint(ig);

		basalelement->JacobianDeterminant(&Jdet,xyz_list,gauss);
		basalelement->NodalFunctions(basis,gauss);

		vz_input->GetInputValue(&vz,gauss);
		gmb_input->GetInputValue(&gmb,gauss);
		fmb_input->GetInputValue(&fmb,gauss);
		base_input->GetInputValue(&bed,gauss);
		groundedice_input->GetInputValue(&phi,gauss);
		if(phi>0) mb=gmb;
		else mb=fmb;

		for(int i=0;i<numnodes;i++) pe->values[i]+=Jdet*gauss->weight*(bed+dt*(mb) + dt*vz)*basis[i];
	}

	/*Clean up and return*/
	xDelete<IssmDouble>(xyz_list);
	xDelete<IssmDouble>(basis);
	delete gauss;
	if(domaintype!=Domain2DhorizontalEnum){basalelement->DeleteMaterials(); delete basalelement;};
	return pe;

}/*}}}*/
void           FreeSurfaceBaseAnalysis::GetB(IssmDouble* B,Element* element,int dim,IssmDouble* xyz_list,Gauss* gauss){/*{{{*/
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
void           FreeSurfaceBaseAnalysis::GetBprime(IssmDouble* Bprime,Element* element,int dim,IssmDouble* xyz_list,Gauss* gauss){/*{{{*/
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
void           FreeSurfaceBaseAnalysis::GetSolutionFromInputs(Vector<IssmDouble>* solution,Element* element){/*{{{*/
	   _error_("not implemented yet");
}/*}}}*/
void           FreeSurfaceBaseAnalysis::GradientJ(Vector<IssmDouble>* gradient,Element* element,int control_type,int control_index){/*{{{*/
	_error_("Not implemented yet");
}/*}}}*/
void           FreeSurfaceBaseAnalysis::InputUpdateFromSolution(IssmDouble* solution,Element* element){/*{{{*/
	element->InputUpdateFromSolutionOneDof(solution,BaseEnum);
}/*}}}*/
void           FreeSurfaceBaseAnalysis::UpdateConstraints(FemModel* femmodel){/*{{{*/

	/*Intermediary*/
	IssmDouble phi,isonbase,base;

	for(int i=0;i<femmodel->elements->Size();i++){

		Element* element=xDynamicCast<Element*>(femmodel->elements->GetObjectByOffset(i));
		if(!element->IsOnBase()) continue;

		int             numnodes = element->GetNumberOfNodes();
		Input* groundedice_input = element->GetInput(MaskGroundediceLevelsetEnum);  _assert_(groundedice_input);
		Input* onbase_input       = element->GetInput(MeshVertexonbaseEnum);          _assert_(onbase_input);
		Input* base_input        = element->GetInput(BaseEnum);                     _assert_(base_input);

		Gauss* gauss=element->NewGauss();
		for(int iv=0;iv<numnodes;iv++){
			gauss->GaussNode(element->GetElementType(),iv);
			onbase_input->GetInputValue(&isonbase,gauss);
			if(isonbase==1.){
				groundedice_input->GetInputValue(&phi,gauss);
				if(phi>=0.){
					base_input->GetInputValue(&base,gauss);
					element->nodes[iv]->ApplyConstraint(0,base);
				}
				else{
					element->nodes[iv]->DofInFSet(0);
				}
			}
		}
		delete gauss;
	}
}/*}}}*/
