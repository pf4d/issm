#include "./MasstransportAnalysis.h"
#include "../toolkits/toolkits.h"
#include "../classes/classes.h"
#include "../shared/shared.h"
#include "../modules/modules.h"

/*Model processing*/
void MasstransportAnalysis::CreateConstraints(Constraints* constraints,IoModel* iomodel){/*{{{*/

	/*Fetch parameters: */
	int stabilization;
	iomodel->FindConstant(&stabilization,"md.masstransport.stabilization");

	/*Do not add constraints in DG,  they are weakly imposed*/
	if(stabilization!=3){
		IoModelToConstraintsx(constraints,iomodel,"md.masstransport.spcthickness",MasstransportAnalysisEnum,P1Enum);
	}

	/*FCT, constraints are imposed using penalties*/
	if(stabilization==4){
		constraints->ActivatePenaltyMethod(MasstransportAnalysisEnum);
	}
}/*}}}*/
void MasstransportAnalysis::CreateLoads(Loads* loads, IoModel* iomodel){/*{{{*/

	/*Intermediaries*/
	int element;
	int penpair_ids[2];
	int count=0;
	int stabilization;
	int numvertex_pairing;

	/*Fetch parameters: */
	iomodel->FindConstant(&stabilization,"md.masstransport.stabilization");

	/*Loads only in DG*/
	if (stabilization==3){

		/*Get faces and elements*/
		CreateFaces(iomodel);
		iomodel->FetchData(1,"md.geometry.thickness");

		/*First load data:*/
		for(int i=0;i<iomodel->numberoffaces;i++){

			/*Get left and right elements*/
			element=iomodel->faces[4*i+2]-1; //faces are [node1 node2 elem1 elem2]

			/*Now, if this element is not in the partition, pass: */
			if(!iomodel->my_elements[element]) continue;

			/* Add load */
			loads->AddObject(new Numericalflux(iomodel->loadcounter+i+1,i,i,iomodel,MasstransportAnalysisEnum));
		}

		/*Free data: */
		iomodel->DeleteData(1,"md.geometry.thickness");
	}

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
							MasstransportAnalysisEnum));
			count++;
		}
	}

	/*free ressources: */
	iomodel->DeleteData(vertex_pairing,"md.masstransport.vertex_pairing");
	iomodel->DeleteData(nodeonbase,"md.mesh.vertexonbase");

}/*}}}*/
void MasstransportAnalysis::CreateNodes(Nodes* nodes,IoModel* iomodel){/*{{{*/

	/*Fetch parameters: */
	int  stabilization;
	iomodel->FindConstant(&stabilization,"md.masstransport.stabilization");

	/*Check in 3d*/
	if(stabilization==3 && iomodel->domaintype==Domain3DEnum) _error_("DG 3d not implemented yet");

	/*Create Nodes either DG or CG depending on stabilization*/
	if(iomodel->domaintype!=Domain2DhorizontalEnum) iomodel->FetchData(2,"md.mesh.vertexonbase","md.mesh.vertexonsurface");
	if(stabilization!=3){
		::CreateNodes(nodes,iomodel,MasstransportAnalysisEnum,P1Enum);
	}
	else{
		::CreateNodes(nodes,iomodel,MasstransportAnalysisEnum,P1DGEnum);
	}
	iomodel->DeleteData(2,"md.mesh.vertexonbase","md.mesh.vertexonsurface");
}/*}}}*/
int  MasstransportAnalysis::DofsPerNode(int** doflist,int domaintype,int approximation){/*{{{*/
	return 1;
}/*}}}*/
void MasstransportAnalysis::UpdateElements(Elements* elements,IoModel* iomodel,int analysis_counter,int analysis_type){/*{{{*/

	int    stabilization,finiteelement;
	bool   dakota_analysis;
	bool   isgroundingline;
	bool   ismovingfront;
	bool   issmb;

	/*Fetch data needed: */
	iomodel->FindConstant(&stabilization,"md.masstransport.stabilization");
	iomodel->FindConstant(&dakota_analysis,"md.qmu.isdakota");
	iomodel->FindConstant(&isgroundingline,"md.transient.isgroundingline");
	iomodel->FindConstant(&ismovingfront,"md.transient.ismovingfront");
	iomodel->FindConstant(&issmb,"md.transient.issmb");

	/*Finite element type*/
	finiteelement = P1Enum;
	if(stabilization==3){
		finiteelement = P1DGEnum;
	}

	/*Update elements: */
	int counter=0;
	for(int i=0;i<iomodel->numberofelements;i++){
		if(iomodel->my_elements[i]){
			Element* element=(Element*)elements->GetObjectByOffset(counter);
			element->Update(i,iomodel,analysis_counter,analysis_type,finiteelement);
			counter++;
		}
	}

	iomodel->FetchDataToInput(elements,"md.geometry.thickness",ThicknessEnum);
	iomodel->FetchDataToInput(elements,"md.geometry.surface",SurfaceEnum);
	iomodel->FetchDataToInput(elements,"md.geometry.base",BaseEnum);
	iomodel->FetchDataToInput(elements,"md.slr.sealevel",SealevelEnum,0);
	iomodel->FetchDataToInput(elements,"md.mask.ice_levelset",MaskIceLevelsetEnum);
	iomodel->FetchDataToInput(elements,"md.mask.groundedice_levelset",MaskGroundediceLevelsetEnum);
	iomodel->FetchDataToInput(elements,"md.basalforcings.groundedice_melting_rate",BasalforcingsGroundediceMeltingRateEnum);
	iomodel->FetchDataToInput(elements,"md.basalforcings.floatingice_melting_rate",BasalforcingsFloatingiceMeltingRateEnum);
	iomodel->FetchDataToInput(elements,"md.initialization.vx",VxEnum);
	iomodel->FetchDataToInput(elements,"md.initialization.vy",VyEnum);

	if(!issmb){
		iomodel->FetchDataToInput(elements,"md.smb.mass_balance",SmbMassBalanceEnum);
	}
	if(stabilization==3){
		iomodel->FetchDataToInput(elements,"md.masstransport.spcthickness",MasstransportSpcthicknessEnum); //for DG, we need the spc in the element
	}
	if(stabilization==4){
		iomodel->FetchDataToInput(elements,"md.masstransport.spcthickness",MasstransportSpcthicknessEnum); //for FCT, we need the spc in the element (penlaties)
	}

	if(iomodel->domaintype!=Domain2DhorizontalEnum){
		iomodel->FetchDataToInput(elements,"md.mesh.vertexonbase",MeshVertexonbaseEnum);
		iomodel->FetchDataToInput(elements,"md.mesh.vertexonsurface",MeshVertexonsurfaceEnum);
	}
}/*}}}*/
void MasstransportAnalysis::UpdateParameters(Parameters* parameters,IoModel* iomodel,int solution_enum,int analysis_enum){/*{{{*/

	int     numoutputs;
	char**  requestedoutputs = NULL;

	parameters->AddObject(iomodel->CopyConstantObject("md.flowequation.isFS",FlowequationIsFSEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.masstransport.isfreesurface",MasstransportIsfreesurfaceEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.masstransport.hydrostatic_adjustment",MasstransportHydrostaticAdjustmentEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.masstransport.stabilization",MasstransportStabilizationEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.masstransport.min_thickness",MasstransportMinThicknessEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.masstransport.penalty_factor",MasstransportPenaltyFactorEnum));

	iomodel->FindConstant(&requestedoutputs,&numoutputs,"md.masstransport.requested_outputs");
	parameters->AddObject(new IntParam(MasstransportNumRequestedOutputsEnum,numoutputs));
	if(numoutputs)parameters->AddObject(new StringArrayParam(MasstransportRequestedOutputsEnum,requestedoutputs,numoutputs));
	iomodel->DeleteData(&requestedoutputs,numoutputs,"md.masstransport.requested_outputs");
	
	

}/*}}}*/

/*Finite Element Analysis*/
void           MasstransportAnalysis::Core(FemModel* femmodel){/*{{{*/
	_error_("not implemented");
}/*}}}*/
ElementVector* MasstransportAnalysis::CreateDVector(Element* element){/*{{{*/
	/*Default, return NULL*/
	return NULL;
}/*}}}*/
ElementMatrix* MasstransportAnalysis::CreateJacobianMatrix(Element* element){/*{{{*/
_error_("Not implemented");
}/*}}}*/
ElementMatrix* MasstransportAnalysis::CreateKMatrix(Element* element){/*{{{*/

	/* Check if ice in element */
	if(!element->IsIceInElement()) return NULL;

	if(!element->IsOnBase()) return NULL;
	Element* basalelement = element->SpawnBasalElement();

	ElementMatrix* Ke = NULL;
	switch(element->FiniteElement()){
		case P1Enum: case P2Enum:
			Ke = CreateKMatrixCG(basalelement);
			break;
		case P1DGEnum:
			Ke = CreateKMatrixDG(basalelement);
			break;
		default:
			_error_("Element type " << EnumToStringx(element->FiniteElement()) << " not supported yet");
	}

	int domaintype;
	element->FindParam(&domaintype,DomainTypeEnum);
	if(domaintype!=Domain2DhorizontalEnum){basalelement->DeleteMaterials(); delete basalelement;};
	return Ke;
}/*}}}*/
ElementMatrix* MasstransportAnalysis::CreateKMatrixCG(Element* element){/*{{{*/

	/* Check if ice in element */
	if(!element->IsIceInElement()) return NULL;

	/*Intermediaries */
	int        stabilization;
	int        domaintype,dim;
	IssmDouble Jdet,D_scalar,dt,h;
	IssmDouble vel,vx,vy,dvxdx,dvydy;
	IssmDouble dvx[2],dvy[2];
	IssmDouble* xyz_list = NULL;

	/*Get problem dimension*/
	element->FindParam(&domaintype,DomainTypeEnum);
	switch(domaintype){
		case Domain2DverticalEnum:   dim = 1; break;
		case Domain2DhorizontalEnum: dim = 2; break;
		case Domain3DEnum:           dim = 2; break;
		default: _error_("mesh "<<EnumToStringx(domaintype)<<" not supported yet");
	}

	/*Fetch number of nodes and dof for this finite element*/
	int numnodes = element->GetNumberOfNodes();

	/*Initialize Element vector and other vectors*/
	ElementMatrix* Ke     = element->NewElementMatrix();
	IssmDouble*    basis  = xNew<IssmDouble>(numnodes);
	IssmDouble*    B      = xNew<IssmDouble>(dim*numnodes);
	IssmDouble*    Bprime = xNew<IssmDouble>(dim*numnodes);
	IssmDouble*    D      = xNewZeroInit<IssmDouble>(dim*dim);

	/*Retrieve all inputs and parameters*/
	element->GetVerticesCoordinates(&xyz_list);
	element->FindParam(&dt,TimesteppingTimeStepEnum);
	element->FindParam(&domaintype,DomainTypeEnum);
	element->FindParam(&stabilization,MasstransportStabilizationEnum);
	Input* vxaverage_input=element->GetInput(VxAverageEnum); _assert_(vxaverage_input);
	Input* vyaverage_input=NULL;
	if(dim==2) vyaverage_input=element->GetInput(VyAverageEnum); _assert_(vyaverage_input);

	h = element->CharacteristicLength();

	/* Start  looping on the number of gaussian points: */
	Gauss* gauss=element->NewGauss(2);
	for(int ig=gauss->begin();ig<gauss->end();ig++){
		gauss->GaussPoint(ig);

		element->JacobianDeterminant(&Jdet,xyz_list,gauss);
		element->NodalFunctions(basis,gauss);

		vxaverage_input->GetInputValue(&vx,gauss);
		vxaverage_input->GetInputDerivativeValue(&dvx[0],xyz_list,gauss);
		if(dim==2){
			vyaverage_input->GetInputValue(&vy,gauss);
			vyaverage_input->GetInputDerivativeValue(&dvy[0],xyz_list,gauss);
		}

		D_scalar=gauss->weight*Jdet;
		TripleMultiply(basis,1,numnodes,1,
					&D_scalar,1,1,0,
					basis,1,numnodes,0,
					&Ke->values[0],1);

		GetB(B,element,dim,xyz_list,gauss);
		GetBprime(Bprime,element,dim,xyz_list,gauss);

		dvxdx=dvx[0];
		if(dim==2) dvydy=dvy[1];
		D_scalar=dt*gauss->weight*Jdet;

		D[0*dim+0]=D_scalar*dvxdx;
		if(dim==2) D[1*dim+1]=D_scalar*dvydy;

		TripleMultiply(B,dim,numnodes,1,
					D,dim,dim,0,
					B,dim,numnodes,0,
					&Ke->values[0],1);

		D[0*dim+0]=D_scalar*vx;
		if(dim==2) D[1*dim+1]=D_scalar*vy;

		TripleMultiply(B,dim,numnodes,1,
					D,dim,dim,0,
					Bprime,dim,numnodes,0,
					&Ke->values[0],1);

		switch(stabilization){
			case 0:
				/*Nothing to be onde*/
				break;
			case 1:
				/*SSA*/
				vxaverage_input->GetInputAverage(&vx);
				if(dim==2) vyaverage_input->GetInputAverage(&vy);
				D[0*dim+0]=h/2.0*fabs(vx);
				if(dim==2) D[1*dim+1]=h/2.0*fabs(vy);
				break;
			case 2:
				if(dim==1){
					vel=fabs(vx)+1.e-8;
					D[0]=h/(2*vel)*vx*vx;
				}
				else{
					/*Streamline upwinding*/
					vel=sqrt(vx*vx+vy*vy)+1.e-8;
					D[0*dim+0]=h/(2*vel)*vx*vx;
					D[1*dim+0]=h/(2*vel)*vy*vx;
					D[0*dim+1]=h/(2*vel)*vx*vy;
					D[1*dim+1]=h/(2*vel)*vy*vy;
				}
				break;
			default:
				_error_("Stabilization "<<stabilization<<" not supported yet");
		}
		if(stabilization==1 || stabilization==2){
			if(dim==1) D[0]=D_scalar*D[0];
			else{
				D[0*dim+0]=D_scalar*D[0*dim+0];
				D[1*dim+0]=D_scalar*D[1*dim+0];
				D[0*dim+1]=D_scalar*D[0*dim+1];
				D[1*dim+1]=D_scalar*D[1*dim+1];
			}

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
	return Ke;
}/*}}}*/
ElementMatrix* MasstransportAnalysis::CreateKMatrixDG(Element* element){/*{{{*/

	/* Check if ice in element */
	if(!element->IsIceInElement()) return NULL;

	/*Intermediaries */
	int        domaintype;
	IssmDouble Jdet,D_scalar,dt,vx,vy;
	IssmDouble* xyz_list = NULL;

	/*Fetch number of nodes and dof for this finite element*/
	int numnodes = element->GetNumberOfNodes();

	/*Initialize Element vector and other vectors*/
	ElementMatrix* Ke     = element->NewElementMatrix();
	IssmDouble*    basis  = xNew<IssmDouble>(numnodes);
	IssmDouble*    B      = xNew<IssmDouble>(2*numnodes);
	IssmDouble*    Bprime = xNew<IssmDouble>(2*numnodes);
	IssmDouble     D[2][2];

	/*Retrieve all inputs and parameters*/
	element->GetVerticesCoordinates(&xyz_list);
	element->FindParam(&dt,TimesteppingTimeStepEnum);
	element->FindParam(&domaintype,DomainTypeEnum);
	Input* vxaverage_input=element->GetInput(VxAverageEnum); _assert_(vxaverage_input);
	Input* vyaverage_input=element->GetInput(VyAverageEnum); _assert_(vyaverage_input);

	/* Start  looping on the number of gaussian points: */
	Gauss* gauss=element->NewGauss(2);
	for(int ig=gauss->begin();ig<gauss->end();ig++){
		gauss->GaussPoint(ig);

		element->JacobianDeterminant(&Jdet,xyz_list,gauss);
		element->NodalFunctions(basis,gauss);

		vxaverage_input->GetInputValue(&vx,gauss);
		vyaverage_input->GetInputValue(&vy,gauss);

		D_scalar=gauss->weight*Jdet;
		TripleMultiply(basis,1,numnodes,1,
					&D_scalar,1,1,0,
					basis,1,numnodes,0,
					&Ke->values[0],1);

		/*WARNING: B and Bprime are inverted compared to CG*/
		GetB(Bprime,element,2,xyz_list,gauss);
		GetBprime(B,element,2,xyz_list,gauss);

		D_scalar = - dt*gauss->weight*Jdet;
		D[0][0]  = D_scalar*vx;
		D[0][1]  = 0.;
		D[1][0]  = 0.;
		D[1][1]  = D_scalar*vy;
		TripleMultiply(B,2,numnodes,1,
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
ElementVector* MasstransportAnalysis::CreatePVector(Element* element){/*{{{*/

	/* Check if ice in element */
	if(!element->IsIceInElement()) return NULL;

	if(!element->IsOnBase()) return NULL;
	Element* basalelement = element->SpawnBasalElement();

	ElementVector* pe = NULL;
	switch(element->FiniteElement()){
		case P1Enum: case P2Enum:
			pe = CreatePVectorCG(basalelement);
			break;
		case P1DGEnum:
			pe = CreatePVectorDG(basalelement);
			break;
		default:
			_error_("Element type " << EnumToStringx(element->FiniteElement()) << " not supported yet");
	}

	int domaintype;
	element->FindParam(&domaintype,DomainTypeEnum);
	if(domaintype!=Domain2DhorizontalEnum){basalelement->DeleteMaterials(); delete basalelement;};
	return pe;
}/*}}}*/
ElementVector* MasstransportAnalysis::CreatePVectorCG(Element* element){/*{{{*/

	/* Check if ice in element */
	if(!element->IsIceInElement()) return NULL;

	/*Intermediaries */
	IssmDouble  Jdet,dt;
	IssmDouble  ms,mb,gmb,fmb,thickness,phi;
	IssmDouble* xyz_list = NULL;

	/*Fetch number of nodes and dof for this finite element*/
	int numnodes = element->GetNumberOfNodes();

	/*Initialize Element vector and other vectors*/
	ElementVector* pe    = element->NewElementVector();
	IssmDouble*    basis = xNew<IssmDouble>(numnodes);

	/*Retrieve all inputs and parameters*/
	element->GetVerticesCoordinates(&xyz_list);
	element->FindParam(&dt,TimesteppingTimeStepEnum);
	Input* gmb_input           = element->GetInput(BasalforcingsGroundediceMeltingRateEnum);  _assert_(gmb_input);
	Input* fmb_input           = element->GetInput(BasalforcingsFloatingiceMeltingRateEnum);  _assert_(fmb_input);
	Input* groundedice_input   = element->GetInput(MaskGroundediceLevelsetEnum);              _assert_(groundedice_input);
	Input* ms_input            = element->GetInput(SmbMassBalanceEnum);                       _assert_(ms_input);
	Input* thickness_input     = element->GetInput(ThicknessEnum);                            _assert_(thickness_input);

	/* Start  looping on the number of gaussian points: */
	Gauss* gauss=element->NewGauss(2);
	for(int ig=gauss->begin();ig<gauss->end();ig++){
		gauss->GaussPoint(ig);

		element->JacobianDeterminant(&Jdet,xyz_list,gauss);
		element->NodalFunctions(basis,gauss);

		ms_input->GetInputValue(&ms,gauss);
		gmb_input->GetInputValue(&gmb,gauss);
		fmb_input->GetInputValue(&fmb,gauss);
		groundedice_input->GetInputValue(&phi,gauss);
		thickness_input->GetInputValue(&thickness,gauss);
		if(phi>0.) mb=gmb;
		else mb=fmb;

		for(int i=0;i<numnodes;i++) pe->values[i]+=Jdet*gauss->weight*(thickness+dt*(ms-mb))*basis[i];
	}

	/*Clean up and return*/
	xDelete<IssmDouble>(xyz_list);
	xDelete<IssmDouble>(basis);
	delete gauss;
	return pe;
}/*}}}*/
ElementVector* MasstransportAnalysis::CreatePVectorDG(Element* element){/*{{{*/

	/* Check if ice in element */
	if(!element->IsIceInElement()) return NULL;

	/*Intermediaries */
	IssmDouble  Jdet,dt;
	IssmDouble  ms,mb,gmb,fmb,thickness,phi;
	IssmDouble* xyz_list = NULL;

	/*Fetch number of nodes and dof for this finite element*/
	int numnodes = element->GetNumberOfNodes();

	/*Initialize Element vector and other vectors*/
	ElementVector* pe    = element->NewElementVector();
	IssmDouble*    basis = xNew<IssmDouble>(numnodes);

	/*Retrieve all inputs and parameters*/
	element->GetVerticesCoordinates(&xyz_list);
	element->FindParam(&dt,TimesteppingTimeStepEnum);
	Input* gmb_input           = element->GetInput(BasalforcingsGroundediceMeltingRateEnum); _assert_(gmb_input);
	Input* fmb_input           = element->GetInput(BasalforcingsFloatingiceMeltingRateEnum); _assert_(fmb_input);
	Input* ms_input            = element->GetInput(SmbMassBalanceEnum);          _assert_(ms_input);
	Input* groundedice_input   = element->GetInput(MaskGroundediceLevelsetEnum);             _assert_(groundedice_input);
	Input* thickness_input     = element->GetInput(ThicknessEnum);                           _assert_(thickness_input);

	/* Start  looping on the number of gaussian points: */
	Gauss* gauss=element->NewGauss(2);
	for(int ig=gauss->begin();ig<gauss->end();ig++){
		gauss->GaussPoint(ig);

		element->JacobianDeterminant(&Jdet,xyz_list,gauss);
		element->NodalFunctions(basis,gauss);

		ms_input->GetInputValue(&ms,gauss);
		gmb_input->GetInputValue(&gmb,gauss);
		fmb_input->GetInputValue(&fmb,gauss);
		groundedice_input->GetInputValue(&phi,gauss);
		thickness_input->GetInputValue(&thickness,gauss);
		if(phi>0) mb=gmb;
		else mb=fmb;

		for(int i=0;i<numnodes;i++) pe->values[i]+=Jdet*gauss->weight*(thickness+dt*(ms-mb))*basis[i];
	}

	/*Clean up and return*/
	xDelete<IssmDouble>(xyz_list);
	xDelete<IssmDouble>(basis);
	delete gauss;
	return pe;
}/*}}}*/
void           MasstransportAnalysis::GetB(IssmDouble* B,Element* element,int dim,IssmDouble* xyz_list,Gauss* gauss){/*{{{*/
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
void           MasstransportAnalysis::GetBprime(IssmDouble* Bprime,Element* element,int dim,IssmDouble* xyz_list,Gauss* gauss){/*{{{*/
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
void           MasstransportAnalysis::GetSolutionFromInputs(Vector<IssmDouble>* solution,Element* element){/*{{{*/
	element->GetSolutionFromInputsOneDof(solution,ThicknessEnum);
}/*}}}*/
void           MasstransportAnalysis::GradientJ(Vector<IssmDouble>* gradient,Element* element,int control_type,int control_index){/*{{{*/
	_error_("Not implemented yet");
}/*}}}*/
void           MasstransportAnalysis::InputUpdateFromSolution(IssmDouble* solution,Element* element){/*{{{*/

	int        i,hydroadjustment,domaintype;
	int*       doflist=NULL;
	IssmDouble rho_ice,rho_water,minthickness;
	Element*   basalelement=NULL;

	element->FindParam(&domaintype,DomainTypeEnum);
	if(domaintype!=Domain2DhorizontalEnum){
		if(!element->IsOnBase()) return;
		basalelement=element->SpawnBasalElement();
	}
	else{
		basalelement = element;
	}

	/*Fetch number of nodes for this finite element*/
	int numnodes = basalelement->GetNumberOfNodes();

	/*Fetch dof list and allocate solution vector*/
	basalelement->GetDofList(&doflist,NoneApproximationEnum,GsetEnum);
	IssmDouble* newthickness   = xNew<IssmDouble>(numnodes);
	IssmDouble* deltathickness = xNew<IssmDouble>(numnodes);
	IssmDouble* newbase        = xNew<IssmDouble>(numnodes);
	IssmDouble* newsurface     = xNew<IssmDouble>(numnodes);
	IssmDouble* oldthickness   = xNew<IssmDouble>(numnodes);
	IssmDouble* oldbase        = xNew<IssmDouble>(numnodes);
	IssmDouble* oldsurface     = xNew<IssmDouble>(numnodes);
	IssmDouble* phi            = xNew<IssmDouble>(numnodes);
	IssmDouble* sealevel       = xNew<IssmDouble>(numnodes);

	/*Use the dof list to index into the solution vector: */
	basalelement->FindParam(&minthickness,MasstransportMinThicknessEnum);
	for(i=0;i<numnodes;i++){
		newthickness[i]=solution[doflist[i]];
		if(xIsNan<IssmDouble>(newthickness[i])) _error_("NaN found in solution vector");
		if(xIsInf<IssmDouble>(newthickness[i])) _error_("Inf found in solution vector");
		/*Constrain thickness to be at least 1m*/
		if(newthickness[i]<minthickness) newthickness[i]=minthickness;
	}

	/*Get previous base, thickness, surfac and current sealevel:*/
	basalelement->GetInputListOnNodes(&oldbase[0],BaseEnum);
	basalelement->GetInputListOnNodes(&oldsurface[0],SurfaceEnum);
	basalelement->GetInputListOnNodes(&oldthickness[0],ThicknessEnum);
	basalelement->GetInputListOnNodes(&phi[0],MaskGroundediceLevelsetEnum);
	basalelement->GetInputListOnNodes(&sealevel[0],SealevelEnum);

	/*What is the delta thickness forcing the sea-level rise core: */
	for(i=0;i<numnodes;i++) deltathickness[i]=newthickness[i]-oldthickness[i];

	/*Find MasstransportHydrostaticAdjustment to figure out how to update the geometry:*/
	basalelement->FindParam(&hydroadjustment,MasstransportHydrostaticAdjustmentEnum);
	rho_ice   = basalelement->GetMaterialParameter(MaterialsRhoIceEnum);
	rho_water = basalelement->GetMaterialParameter(MaterialsRhoSeawaterEnum);

	for(i=0;i<numnodes;i++) {
		if (phi[i]>0.){ //this is an ice sheet: just add thickness to base.
			newsurface[i] = oldbase[i]+newthickness[i]; //surface = oldbase + newthickness
			newbase[i]     = oldbase[i];                 //same base: do nothing
		}
		else{ //this is an ice shelf: hydrostatic equilibrium*/
			if(hydroadjustment==AbsoluteEnum){
				newsurface[i] = newthickness[i]*(1-rho_ice/rho_water)+sealevel[i];
				newbase[i]     = newthickness[i]*(-rho_ice/rho_water)+sealevel[i];
			}
			else if(hydroadjustment==IncrementalEnum){
				newsurface[i] = oldsurface[i]+(1.0-rho_ice/rho_water)*(newthickness[i]-oldthickness[i])+sealevel[i]; //surface = oldsurface + (1-di) * dH
				newbase[i]     = oldbase[i]-rho_ice/rho_water*(newthickness[i]-oldthickness[i])+sealevel[i]; //base               = oldbed + di * dH
			}
			else _error_("Hydrostatic adjustment " << hydroadjustment << " (" << EnumToStringx(hydroadjustment) << ") not supported yet");
		}
	}

	/*Add input to the element: */
	element->AddBasalInput(ThicknessEnum,newthickness,P1Enum);
	element->AddBasalInput(SealevelriseDeltathicknessEnum,deltathickness,P1Enum);
	element->AddBasalInput(SurfaceEnum,newsurface,P1Enum);
	element->AddBasalInput(BaseEnum,newbase,P1Enum);

	/*Free ressources:*/
	xDelete<IssmDouble>(newthickness);
	xDelete<IssmDouble>(newbase);
	xDelete<IssmDouble>(newsurface);
	xDelete<IssmDouble>(oldthickness);
	xDelete<IssmDouble>(deltathickness);
	xDelete<IssmDouble>(oldbase);
	xDelete<IssmDouble>(oldsurface);
	xDelete<IssmDouble>(phi);
	xDelete<IssmDouble>(sealevel);
	xDelete<int>(doflist);
	if(domaintype!=Domain2DhorizontalEnum){basalelement->DeleteMaterials(); delete basalelement;};
}/*}}}*/
void           MasstransportAnalysis::UpdateConstraints(FemModel* femmodel){/*{{{*/
	SetActiveNodesLSMx(femmodel);
}/*}}}*/

/*Flux Correction Transport*/
ElementMatrix* MasstransportAnalysis::CreateFctKMatrix(Element* element){/*{{{*/

	/* Check if ice in element */
	if(!element->IsIceInElement()) return NULL;

	/*Intermediaries */
	IssmDouble Jdet;
	IssmDouble vx,vy;
	IssmDouble* xyz_list = NULL;

	/*Fetch number of nodes and dof for this finite element*/
	int numnodes = element->GetNumberOfNodes();
	int dim      = 2;

	/*Initialize Element vector and other vectors*/
	ElementMatrix* Ke     = element->NewElementMatrix();
	IssmDouble*    B      = xNew<IssmDouble>(dim*numnodes);
	IssmDouble*    Bprime = xNew<IssmDouble>(dim*numnodes);
	IssmDouble*    D      = xNewZeroInit<IssmDouble>(dim*dim);

	/*Retrieve all inputs and parameters*/
	element->GetVerticesCoordinates(&xyz_list);
	Input* vxaverage_input=element->GetInput(VxEnum); _assert_(vxaverage_input);
	Input* vyaverage_input=element->GetInput(VyEnum); _assert_(vyaverage_input);

	/* Start  looping on the number of gaussian points: */
	Gauss* gauss=element->NewGauss(2);
	for(int ig=gauss->begin();ig<gauss->end();ig++){
		gauss->GaussPoint(ig);

		element->JacobianDeterminant(&Jdet,xyz_list,gauss);
		GetB(B,element,dim,xyz_list,gauss);
		GetBprime(Bprime,element,dim,xyz_list,gauss);
		vxaverage_input->GetInputValue(&vx,gauss);
		vyaverage_input->GetInputValue(&vy,gauss);

		D[0*dim+0] = -gauss->weight*vx*Jdet;
		D[1*dim+1] = -gauss->weight*vy*Jdet;

		TripleMultiply(B,dim,numnodes,1,
					D,dim,dim,0,
					Bprime,dim,numnodes,0,
					&Ke->values[0],1);

	}

	/*Clean up and return*/
	xDelete<IssmDouble>(xyz_list);
	xDelete<IssmDouble>(B);
	xDelete<IssmDouble>(Bprime);
	xDelete<IssmDouble>(D);
	delete gauss;
	return Ke;
}/*}}}*/
ElementMatrix* MasstransportAnalysis::CreateMassMatrix(Element* element){/*{{{*/

	/* Check if ice in element */
	if(!element->IsIceInElement()) return NULL;

	/*Intermediaries*/
	IssmDouble  D,Jdet;
	IssmDouble* xyz_list = NULL;

	/*Fetch number of nodes and dof for this finite element*/
	int numnodes = element->GetNumberOfNodes();

	/*Initialize Element vector and other vectors*/
	ElementMatrix* Me     = element->NewElementMatrix();
	IssmDouble*    basis  = xNew<IssmDouble>(numnodes);

	/*Retrieve all inputs and parameters*/
	element->GetVerticesCoordinates(&xyz_list);

	/* Start  looping on the number of gaussian points: */
	Gauss* gauss=element->NewGauss(2);
	for(int ig=gauss->begin();ig<gauss->end();ig++){
		gauss->GaussPoint(ig);

		element->JacobianDeterminant(&Jdet,xyz_list,gauss);
		element->NodalFunctions(basis,gauss);

		D=gauss->weight*Jdet;
		TripleMultiply(basis,1,numnodes,1,
					&D,1,1,0,
					basis,1,numnodes,0,
					&Me->values[0],1);
	}

	/*Clean up and return*/
	xDelete<IssmDouble>(xyz_list);
	xDelete<IssmDouble>(basis);
	delete gauss;
	return Me;
}/*}}}*/
void           MasstransportAnalysis::FctKMatrix(Matrix<IssmDouble>** pKff,Matrix<IssmDouble>** pKfs,FemModel* femmodel){/*{{{*/

	/*Output*/
	Matrix<IssmDouble>* Kff = NULL;
	Matrix<IssmDouble>* Kfs = NULL;

	/*Initialize Jacobian Matrix*/
	AllocateSystemMatricesx(&Kff,&Kfs,NULL,NULL,femmodel);

	/*Create and assemble matrix*/
	for(int i=0;i<femmodel->elements->Size();i++){
		Element*       element = xDynamicCast<Element*>(femmodel->elements->GetObjectByOffset(i));
		ElementMatrix* Ke     = this->CreateFctKMatrix(element);
		if(Ke) Ke->AddToGlobal(Kff,Kfs);
		delete Ke;
	}
	Kff->Assemble();
	Kfs->Assemble();

	/*Assign output pointer*/
	*pKff=Kff;
	if(pKfs){
		*pKfs=Kfs;
	}
	else{
		delete Kfs;
	}
}/*}}}*/
void           MasstransportAnalysis::LumpedMassMatrix(Vector<IssmDouble>** pMlff,FemModel* femmodel){/*{{{*/

	/*Intermediaries*/
	int  configuration_type;

	/*Initialize Lumped mass matrix (actually we just save its diagonal)*/
	femmodel->parameters->FindParam(&configuration_type,ConfigurationTypeEnum);
	int fsize      = femmodel->nodes->NumberOfDofs(configuration_type,FsetEnum);
	int flocalsize = femmodel->nodes->NumberOfDofsLocal(configuration_type,FsetEnum);
	Vector<IssmDouble>* Mlff = new Vector<IssmDouble>(flocalsize,fsize);

	/*Create and assemble matrix*/
	for(int i=0;i<femmodel->elements->Size();i++){
		Element*       element = xDynamicCast<Element*>(femmodel->elements->GetObjectByOffset(i));
		ElementMatrix* MLe     = this->CreateMassMatrix(element);
		if(MLe){
			MLe->Lump();
			MLe->AddDiagonalToGlobal(Mlff);
		}
		delete MLe;
	}
	Mlff->Assemble();

	/*Assign output pointer*/
	*pMlff=Mlff;
}/*}}}*/
void           MasstransportAnalysis::MassMatrix(Matrix<IssmDouble>** pMff,FemModel* femmodel){/*{{{*/

	/*Initialize Mass matrix*/
	Matrix<IssmDouble> *Mff = NULL;
	AllocateSystemMatricesx(&Mff,NULL,NULL,NULL,femmodel);

	/*Create and assemble matrix*/
	for(int i=0;i<femmodel->elements->Size();i++){
		Element*       element = xDynamicCast<Element*>(femmodel->elements->GetObjectByOffset(i));
		ElementMatrix* MLe     = this->CreateMassMatrix(element);
		if(MLe){
			MLe->AddToGlobal(Mff);
		}
		delete MLe;
	}
	Mff->Assemble();

	/*Assign output pointer*/
	*pMff=Mff;
}/*}}}*/
