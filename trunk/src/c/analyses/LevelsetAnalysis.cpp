#ifdef HAVE_CONFIG_H
   #include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif
#include "./LevelsetAnalysis.h"
#include "../toolkits/toolkits.h"
#include "../classes/classes.h"
#include "../shared/shared.h"
#include "../modules/modules.h"
#include "../solutionsequences/solutionsequences.h"

void LevelsetAnalysis::CreateConstraints(Constraints* constraints,IoModel* iomodel){/*{{{*/
	IoModelToConstraintsx(constraints,iomodel,"md.levelset.spclevelset",LevelsetAnalysisEnum,P1Enum);
}
/*}}}*/
void LevelsetAnalysis::CreateLoads(Loads* loads, IoModel* iomodel){/*{{{*/
	return;
}/*}}}*/
void LevelsetAnalysis::CreateNodes(Nodes* nodes,IoModel* iomodel){/*{{{*/
	int finiteelement=P1Enum;
	if(iomodel->domaintype!=Domain2DhorizontalEnum) iomodel->FetchData(2,"md.mesh.vertexonbase","md.mesh.vertexonsurface");
	::CreateNodes(nodes,iomodel,LevelsetAnalysisEnum,finiteelement);
	iomodel->DeleteData(2,"md.mesh.vertexonbase","md.mesh.vertexonsurface");
}
/*}}}*/
int  LevelsetAnalysis::DofsPerNode(int** doflist,int domaintype,int approximation){/*{{{*/
	return 1;
}
/*}}}*/
void LevelsetAnalysis::UpdateElements(Elements* elements,IoModel* iomodel,int analysis_counter,int analysis_type){/*{{{*/

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
	
	iomodel->FetchDataToInput(elements,"md.mask.ice_levelset",MaskIceLevelsetEnum);
	iomodel->FetchDataToInput(elements,"md.initialization.vx",VxEnum);
	iomodel->FetchDataToInput(elements,"md.initialization.vy",VyEnum);

	/*Get moving front parameters*/
	int  calvinglaw;
	iomodel->FindConstant(&calvinglaw,"md.calving.law");
	switch(calvinglaw){
		case DefaultCalvingEnum:
			iomodel->FetchDataToInput(elements,"md.calving.calvingrate",CalvingCalvingrateEnum);
			iomodel->FetchDataToInput(elements,"md.calving.meltingrate",CalvingMeltingrateEnum);
			break;
		case CalvingLevermannEnum:
			iomodel->FetchDataToInput(elements,"md.calving.coeff",CalvinglevermannCoeffEnum);
			iomodel->FetchDataToInput(elements,"md.calving.meltingrate",CalvinglevermannMeltingrateEnum);
			break;
		case CalvingDevEnum:
			iomodel->FetchDataToInput(elements,"md.calving.meltingrate",CalvingMeltingrateEnum);
			break;
		case CalvingMinthicknessEnum:
			iomodel->FetchDataToInput(elements,"md.calving.meltingrate",CalvingMeltingrateEnum);
			break;
		default:
			_error_("Calving law "<<EnumToStringx(calvinglaw)<<" not supported yet");
	}
}
/*}}}*/
void LevelsetAnalysis::UpdateParameters(Parameters* parameters,IoModel* iomodel,int solution_enum,int analysis_enum){/*{{{*/
	parameters->AddObject(iomodel->CopyConstantObject("md.levelset.stabilization",LevelsetStabilizationEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.levelset.reinit_frequency",LevelsetReinitFrequencyEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.levelset.calving_max",CalvingMaxEnum));
	int  calvinglaw;
	iomodel->FindConstant(&calvinglaw,"md.calving.law");
	switch(calvinglaw){
		case DefaultCalvingEnum:
		case CalvingLevermannEnum:
			break;
		case CalvingDevEnum:
			parameters->AddObject(iomodel->CopyConstantObject("md.calving.stress_threshold_groundedice",CalvingStressThresholdGroundediceEnum));
			parameters->AddObject(iomodel->CopyConstantObject("md.calving.stress_threshold_floatingice",CalvingStressThresholdFloatingiceEnum));
			break;
		case CalvingMinthicknessEnum:
			parameters->AddObject(iomodel->CopyConstantObject("md.calving.min_thickness",CalvingMinthicknessEnum));
			break;
		default:
			_error_("Calving law "<<EnumToStringx(calvinglaw)<<" not supported yet");
	}
	return;
}
/*}}}*/

/*Finite element Analysis*/
void           LevelsetAnalysis::Core(FemModel* femmodel){/*{{{*/

	/*parameters: */
	int  stabilization;
	bool save_results;
	femmodel->parameters->FindParam(&save_results,SaveResultsEnum);
	femmodel->parameters->FindParam(&stabilization,LevelsetStabilizationEnum);

	/*activate formulation: */
	femmodel->SetCurrentConfiguration(LevelsetAnalysisEnum);

	if(VerboseSolution()) _printf0_("call computational core:\n");
	if(stabilization==4){
		solutionsequence_fct(femmodel);
	}
	else{
		solutionsequence_linear(femmodel);
	}

	if(save_results){
		if(VerboseSolution()) _printf0_("   saving results\n");
		int outputs[1] = {MaskIceLevelsetEnum};
		femmodel->RequestedOutputsx(&femmodel->results,&outputs[0],1);
	}
}/*}}}*/
ElementVector* LevelsetAnalysis::CreateDVector(Element* element){/*{{{*/
	/*Default, return NULL*/
	return NULL;
}/*}}}*/
ElementMatrix* LevelsetAnalysis::CreateJacobianMatrix(Element* element){/*{{{*/
	/* Jacobian required for the Newton solver */
	_error_("not implemented yet");
}/*}}}*/
ElementMatrix* LevelsetAnalysis::CreateKMatrix(Element* element){/*{{{*/

	if(!element->IsOnBase()) return NULL;
	Element* basalelement = element->SpawnBasalElement();

	/*Intermediaries */
	int  stabilization,dim, domaintype, calvinglaw;
	int i, row, col;
	IssmDouble kappa;
	IssmDouble Jdet, dt, D_scalar;
	IssmDouble h,hx,hy,hz;
	IssmDouble vel;
	IssmDouble norm_dlsf, norm_calving, calvingrate, meltingrate, groundedice;
	IssmDouble calvingmax;
	IssmDouble* xyz_list = NULL;

	/*Get problem dimension and whether there is moving front or not*/
	basalelement->FindParam(&domaintype,DomainTypeEnum);
	basalelement->FindParam(&calvinglaw,CalvingLawEnum);
	basalelement->FindParam(&stabilization,LevelsetStabilizationEnum);
	switch(domaintype){
		case Domain2DverticalEnum:   dim = 1; break;
		case Domain2DhorizontalEnum: dim = 2; break;
		case Domain3DEnum:           dim = 2; break;
		default: _error_("mesh "<<EnumToStringx(domaintype)<<" not supported yet");
	}

	/*Calving threshold*/
	
	/*Fetch number of nodes and dof for this finite element*/
	int numnodes    = basalelement->GetNumberOfNodes();

	/*Initialize Element vector and other vectors*/
	ElementMatrix* Ke       = basalelement->NewElementMatrix();
	IssmDouble*    basis    = xNew<IssmDouble>(numnodes);
	IssmDouble*    B        = xNew<IssmDouble>(dim*numnodes);
	IssmDouble*    Bprime   = xNew<IssmDouble>(dim*numnodes);
	IssmDouble*    D        = xNew<IssmDouble>(dim*dim);
	IssmDouble*    v        = xNew<IssmDouble>(dim);
	IssmDouble*    w        = xNew<IssmDouble>(dim);
	IssmDouble*    c        = xNewZeroInit<IssmDouble>(dim);
	IssmDouble*    m        = xNewZeroInit<IssmDouble>(dim);
	IssmDouble*    dlsf     = xNew<IssmDouble>(dim);

	/*Retrieve all inputs and parameters*/
	basalelement->GetVerticesCoordinates(&xyz_list);
	basalelement->FindParam(&dt,TimesteppingTimeStepEnum);
	basalelement->FindParam(&calvingmax,CalvingMaxEnum);
	Input* vx_input           = NULL;
	Input* vy_input           = NULL;
	Input* calvingratex_input = NULL;
	Input* calvingratey_input = NULL;
	Input* lsf_slopex_input   = NULL;
	Input* lsf_slopey_input   = NULL;
	Input* calvingrate_input  = NULL;
	Input* meltingrate_input  = NULL;
	Input* gr_input           = NULL;

	/*Load velocities*/
	switch(domaintype){
		case Domain2DverticalEnum:
			vx_input=basalelement->GetInput(VxEnum); _assert_(vx_input);
			break;
		case Domain2DhorizontalEnum:
			vx_input=basalelement->GetInput(VxEnum); _assert_(vx_input);
			vy_input=basalelement->GetInput(VyEnum); _assert_(vy_input);
			gr_input=basalelement->GetInput(MaskGroundediceLevelsetEnum); _assert_(gr_input);
			break;
		case Domain3DEnum:
			vx_input=basalelement->GetInput(VxAverageEnum); _assert_(vx_input);
			vy_input=basalelement->GetInput(VyAverageEnum); _assert_(vy_input);
			gr_input=basalelement->GetInput(MaskGroundediceLevelsetEnum); _assert_(gr_input);
			break;
		default: _error_("mesh "<<EnumToStringx(domaintype)<<" not supported yet");
	}

	/*Load calving inputs*/
	switch(calvinglaw){
		case DefaultCalvingEnum:
		case CalvingDevEnum:
			lsf_slopex_input  = basalelement->GetInput(LevelsetfunctionSlopeXEnum); _assert_(lsf_slopex_input);
			if(dim==2) lsf_slopey_input  = basalelement->GetInput(LevelsetfunctionSlopeYEnum); _assert_(lsf_slopey_input);
			calvingrate_input = basalelement->GetInput(CalvingCalvingrateEnum);     _assert_(calvingrate_input);
			meltingrate_input = basalelement->GetInput(CalvingMeltingrateEnum);     _assert_(meltingrate_input);
			break;
		case CalvingLevermannEnum:
			switch(domaintype){
				case Domain2DverticalEnum:
					calvingratex_input=basalelement->GetInput(CalvingratexEnum); _assert_(calvingratex_input);
					break;
				case Domain2DhorizontalEnum:
					calvingratex_input=basalelement->GetInput(CalvingratexEnum); _assert_(calvingratex_input);
					calvingratey_input=basalelement->GetInput(CalvingrateyEnum); _assert_(calvingratey_input);
					break;
				case Domain3DEnum:
					calvingratex_input=basalelement->GetInput(CalvingratexAverageEnum); _assert_(calvingratex_input);
					calvingratey_input=basalelement->GetInput(CalvingrateyAverageEnum); _assert_(calvingratey_input);
					break;
				default: _error_("mesh "<<EnumToStringx(domaintype)<<" not supported yet");
			}
			meltingrate_input = basalelement->GetInput(CalvinglevermannMeltingrateEnum);     _assert_(meltingrate_input);
			break;
		case CalvingMinthicknessEnum:
			lsf_slopex_input  = basalelement->GetInput(LevelsetfunctionSlopeXEnum); _assert_(lsf_slopex_input);
			if(dim==2) lsf_slopey_input  = basalelement->GetInput(LevelsetfunctionSlopeYEnum); _assert_(lsf_slopey_input);
			meltingrate_input = basalelement->GetInput(CalvingMeltingrateEnum);     _assert_(meltingrate_input);
			break;
		default:
			_error_("Calving law "<<EnumToStringx(calvinglaw)<<" not supported yet");
	}

	/* Start  looping on the number of gaussian points: */
	Gauss* gauss=basalelement->NewGauss(2);
	for(int ig=gauss->begin();ig<gauss->end();ig++){
		gauss->GaussPoint(ig);

		basalelement->JacobianDeterminant(&Jdet,xyz_list,gauss);
		D_scalar=gauss->weight*Jdet;

		/* Transient */
		if(dt!=0.){
			basalelement->NodalFunctions(basis,gauss);
			TripleMultiply(basis,numnodes,1,0,
						&D_scalar,1,1,0,
						basis,1,numnodes,0,
						&Ke->values[0],1);
			D_scalar*=dt;
		}

		/* Advection */
		GetB(B,basalelement,xyz_list,gauss); 
		GetBprime(Bprime,basalelement,xyz_list,gauss); 
		vx_input->GetInputValue(&v[0],gauss);
		vy_input->GetInputValue(&v[1],gauss); 
		gr_input->GetInputValue(&groundedice,gauss);

		/*Get calving speed*/
		switch(calvinglaw){
			case DefaultCalvingEnum:
			case CalvingDevEnum:
				lsf_slopex_input->GetInputValue(&dlsf[0],gauss);
				if(dim==2) lsf_slopey_input->GetInputValue(&dlsf[1],gauss);
				calvingrate_input->GetInputValue(&calvingrate,gauss);
				meltingrate_input->GetInputValue(&meltingrate,gauss);

				/*Limit calving rate to c <= v + 3 km/yr */
				vel=sqrt(v[0]*v[0] + v[1]*v[1]);
				if(calvingrate>calvingmax+vel) calvingrate = vel+calvingmax;
				if(groundedice<0) meltingrate = 0.;

				norm_dlsf=0.;
				for(i=0;i<dim;i++) norm_dlsf+=pow(dlsf[i],2);
				norm_dlsf=sqrt(norm_dlsf);

				if(norm_dlsf>1.e-10)
				 for(i=0;i<dim;i++){ 
					 c[i]=calvingrate*dlsf[i]/norm_dlsf; m[i]=meltingrate*dlsf[i]/norm_dlsf;
				 }
				else
				 for(i=0;i<dim;i++){
					 c[i]=0.; m[i]=0.;
				 }
				break;

			case CalvingLevermannEnum:
				calvingratex_input->GetInputValue(&c[0],gauss);
				if(dim==2) calvingratey_input->GetInputValue(&c[1],gauss);
				meltingrate_input->GetInputValue(&meltingrate,gauss);
				norm_calving=0.;
				for(i=0;i<dim;i++) norm_calving+=pow(c[i],2);
				norm_calving=sqrt(norm_calving)+1.e-14;
				for(i=0;i<dim;i++) m[i]=meltingrate*c[i]/norm_calving;
				break;

			case CalvingMinthicknessEnum:
				lsf_slopex_input->GetInputValue(&dlsf[0],gauss);
				if(dim==2) lsf_slopey_input->GetInputValue(&dlsf[1],gauss);
				meltingrate_input->GetInputValue(&meltingrate,gauss);

				norm_dlsf=0.;
				for(i=0;i<dim;i++) norm_dlsf+=pow(dlsf[i],2);
				norm_dlsf=sqrt(norm_dlsf);

				if(norm_dlsf>1.e-10)
				 for(i=0;i<dim;i++){ 
					 c[i]=0.;
					 m[i]=meltingrate*dlsf[i]/norm_dlsf;
				 }
				else
				 for(i=0;i<dim;i++){
					 c[i]=0.;
					 m[i]=0.;
				 }
				break;

			default:
				_error_("Calving law "<<EnumToStringx(calvinglaw)<<" not supported yet");
		}

		/*Levelset speed is ice velocity - calving rate*/
		for(i=0;i<dim;i++) w[i]=v[i]-c[i]-m[i];

		/*Compute D*/
		for(row=0;row<dim;row++){
			for(col=0;col<dim;col++){
				if(row==col)
				 D[row*dim+col]=D_scalar*w[row];
				else
				 D[row*dim+col]=0.;
			}
		}

		TripleMultiply(B,dim,numnodes,1,
					D,dim,dim,0,
					Bprime,dim,numnodes,0,
					&Ke->values[0],1);

		/* Stabilization */
		vel=0.;
		for(i=0;i<dim;i++) vel+=w[i]*w[i];
		vel=sqrt(vel)+1.e-14;
		switch(stabilization){
			case 0:
				// no stabilization, do nothing
				break;
			case 1:
				/* Artificial Diffusion */
				basalelement->ElementSizes(&hx,&hy,&hz);
				h=sqrt( pow(hx*w[0]/vel,2) + pow(hy*w[1]/vel,2) ); 
				kappa=h*vel/2.;
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
				break;	
			case 2:
				/* Streamline Upwinding */
				basalelement->ElementSizes(&hx,&hy,&hz);
				h=sqrt( pow(hx*w[0]/vel,2) + pow(hy*w[1]/vel,2) );
				for(row=0;row<dim;row++) 
					for(col=0;col<dim;col++) 
						D[row*dim+col] = D_scalar*h/(2.*vel)*w[row]*w[col];

				TripleMultiply(Bprime,dim,numnodes,1,
							D,dim,dim,0,
							Bprime,dim,numnodes,0,
							&Ke->values[0],1);
				break;
			default:
				_error_("unknown type of stabilization in LevelsetAnalysis.cpp");
		}
	}

	/*Clean up and return*/
	xDelete<IssmDouble>(xyz_list);
	xDelete<IssmDouble>(basis);
	xDelete<IssmDouble>(B);
	xDelete<IssmDouble>(D);
	xDelete<IssmDouble>(Bprime);
	xDelete<IssmDouble>(v);
	xDelete<IssmDouble>(w);
	xDelete<IssmDouble>(c);
	xDelete<IssmDouble>(m);
	xDelete<IssmDouble>(dlsf);
	delete gauss;
	if(domaintype!=Domain2DhorizontalEnum){basalelement->DeleteMaterials(); delete basalelement;};
	return Ke;
}/*}}}*/
ElementVector* LevelsetAnalysis::CreatePVector(Element* element){/*{{{*/
	
	if(!element->IsOnBase()) return NULL;
	Element* basalelement = element->SpawnBasalElement();

	/*Intermediaries */
	int i, ig, domaintype;
	IssmDouble  Jdet,dt;
	IssmDouble  lsf;
	IssmDouble* xyz_list = NULL;
	
	/*Fetch number of nodes and dof for this finite element*/
	int numnodes = basalelement->GetNumberOfNodes();

	/*Initialize Element vector*/
	ElementVector* pe = basalelement->NewElementVector();
	basalelement->FindParam(&dt,TimesteppingTimeStepEnum);
	
	if(dt!=0.){
		/*Initialize basis vector*/
		IssmDouble*    basis = xNew<IssmDouble>(numnodes);

		/*Retrieve all inputs and parameters*/
		basalelement->GetVerticesCoordinates(&xyz_list);
		Input* levelset_input     = basalelement->GetInput(MaskIceLevelsetEnum);                    _assert_(levelset_input);

		/* Start  looping on the number of gaussian points: */
		Gauss* gauss=basalelement->NewGauss(2);
		for(ig=gauss->begin();ig<gauss->end();ig++){
			gauss->GaussPoint(ig);

			basalelement->JacobianDeterminant(&Jdet,xyz_list,gauss);
			basalelement->NodalFunctions(basis,gauss);

			/* old function value */
			levelset_input->GetInputValue(&lsf,gauss);
			for(i=0;i<numnodes;i++) pe->values[i]+=Jdet*gauss->weight*lsf*basis[i];
		}

		/*Clean up and return*/
		xDelete<IssmDouble>(xyz_list);
		xDelete<IssmDouble>(basis);
		basalelement->FindParam(&domaintype,DomainTypeEnum);
		if(domaintype!=Domain2DhorizontalEnum){basalelement->DeleteMaterials(); delete basalelement;};
		delete gauss;
	}

	return pe;
}/*}}}*/
void           LevelsetAnalysis::GetB(IssmDouble* B,Element* element,IssmDouble* xyz_list,Gauss* gauss){/*{{{*/
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
void           LevelsetAnalysis::GetBprime(IssmDouble* Bprime,Element* element,IssmDouble* xyz_list,Gauss* gauss){/*{{{*/
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
IssmDouble     LevelsetAnalysis::GetDistanceToStraight(IssmDouble* q, IssmDouble* s0, IssmDouble* s1){/*{{{*/
	// returns distance d of point q to straight going through points s0, s1
	// d=|a x b|/|b|
	// with a=q-s0, b=s1-s0
	
	/* Intermediaries */
	const int dim=2;
	int i;
	IssmDouble a[dim], b[dim];
	IssmDouble norm_b;

	for(i=0;i<dim;i++){
		a[i]=q[i]-s0[i];
		b[i]=s1[i]-s0[i];
	}
	
	norm_b=0.;
	for(i=0;i<dim;i++)
		norm_b+=b[i]*b[i];
	norm_b=sqrt(norm_b);
	_assert_(norm_b>0.);

	return fabs(a[0]*b[1]-a[1]*b[0])/norm_b;
}/*}}}*/
void           LevelsetAnalysis::GetSolutionFromInputs(Vector<IssmDouble>* solution,Element* element){/*{{{*/
	_error_("not implemented yet");
}/*}}}*/
void           LevelsetAnalysis::GradientJ(Vector<IssmDouble>* gradient,Element* element,int control_type,int control_index){/*{{{*/
	_error_("Not implemented yet");
}/*}}}*/
void           LevelsetAnalysis::InputUpdateFromSolution(IssmDouble* solution,Element* element){/*{{{*/

	int domaintype;
	element->FindParam(&domaintype,DomainTypeEnum);
	switch(domaintype){
		case Domain2DhorizontalEnum:
			element->InputUpdateFromSolutionOneDof(solution,MaskIceLevelsetEnum);
			break;
		case Domain3DEnum:
			element->InputUpdateFromSolutionOneDofCollapsed(solution,MaskIceLevelsetEnum);
			break;
		default: _error_("mesh "<<EnumToStringx(domaintype)<<" not supported yet");
	}
}/*}}}*/
void           LevelsetAnalysis::SetDistanceOnIntersectedElements(FemModel* femmodel){/*{{{*/

	/* Intermediaries */
	int i,k;
	IssmDouble dmaxp=0.,dmaxm=0,val=0.;

	/*Initialize vector with number of vertices*/
	int numvertices=femmodel->vertices->NumberOfVertices();
	Element* element;

	Vector<IssmDouble>* vec_dist_zerolevelset = NULL;
	GetVectorFromInputsx(&vec_dist_zerolevelset, femmodel, MaskIceLevelsetEnum, VertexPIdEnum);
	
	/* set NaN on elements intersected by zero levelset */
	for(i=0;i<femmodel->elements->Size();i++){
		element=xDynamicCast<Element*>(femmodel->elements->GetObjectByOffset(i));
		if(element->IsZeroLevelset(MaskIceLevelsetEnum))
			for(k=0;k<element->GetNumberOfVertices();k++)
				vec_dist_zerolevelset->SetValue(element->vertices[k]->Sid(),NAN,INS_VAL); 
	}

	/* set distance on elements intersected by zero levelset */
	for(i=0;i<femmodel->elements->Size();i++){
		element=xDynamicCast<Element*>(femmodel->elements->GetObjectByOffset(i));
		if(element->IsZeroLevelset(MaskIceLevelsetEnum)){
			SetDistanceToZeroLevelsetElement(vec_dist_zerolevelset, element);

			/* Get maximum distance to interface along vertices */
			for(k=0;k<element->GetNumberOfVertices();k++){
					vec_dist_zerolevelset->GetValue(&val,element->vertices[k]->Sid()); 
					if((val>0.) && (val>dmaxp))
						 dmaxp=val;
					else if((val<0.) && (val<dmaxm))
						 dmaxm=val;
			}
		}
	}

	/* set all none intersected vertices to max/min distance */
	for(i=0;i<numvertices;i++){
		vec_dist_zerolevelset->GetValue(&val,i);
		if(val==1.) //FIXME: improve check
			vec_dist_zerolevelset->SetValue(i,3.*dmaxp,INS_VAL);
		else if(val==-1.)
			vec_dist_zerolevelset->SetValue(i,3.*dmaxm,INS_VAL);
	}

	/*Assemble vector and serialize */
	vec_dist_zerolevelset->Assemble();
	IssmDouble* dist_zerolevelset=vec_dist_zerolevelset->ToMPISerial();
	InputUpdateFromVectorx(femmodel,dist_zerolevelset,MaskIceLevelsetEnum,VertexSIdEnum);

	/*Clean up and return*/
	delete vec_dist_zerolevelset;
	delete dist_zerolevelset;
}/*}}}*/
void           LevelsetAnalysis::SetDistanceToZeroLevelsetElement(Vector<IssmDouble>* vec_signed_dist, Element* element){/*{{{*/

	if(!element->IsZeroLevelset(MaskIceLevelsetEnum))
		return;

	/* Intermediaries */
	const int dim=3;
	int i,d;
	int numvertices=element->GetNumberOfVertices();
	IssmDouble s0[dim], s1[dim], v[dim];
	IssmDouble dist,lsf_old;

	IssmDouble* lsf = xNew<IssmDouble>(numvertices);
	IssmDouble* sign_lsf = xNew<IssmDouble>(numvertices);
	IssmDouble* signed_dist = xNew<IssmDouble>(numvertices);
	IssmDouble* xyz_list = NULL;
	IssmDouble* xyz_list_zero = NULL;

	/* retrieve inputs and parameters */
	element->GetVerticesCoordinates(&xyz_list);
	element->GetInputListOnVertices(lsf,MaskIceLevelsetEnum);

	/* get sign of levelset function */
	for(i=0;i<numvertices;i++)
		sign_lsf[i]=(lsf[i]>=0.?1.:-1.);

	element->ZeroLevelsetCoordinates(&xyz_list_zero, xyz_list, MaskIceLevelsetEnum);
	for(d=0;d<dim;d++){
		s0[d]=xyz_list_zero[0+d];
		s1[d]=xyz_list_zero[3+d];
	}

	/* get signed_distance of vertices to zero levelset straight */
	for(i=0;i<numvertices;i++){
		for(d=0;d<dim;d++)
			v[d]=xyz_list[3*i+d];
		dist=GetDistanceToStraight(&v[0],&s0[0],&s1[0]);
		signed_dist[i]=sign_lsf[i]*dist;
	}
	
	/* insert signed_distance into vec_signed_dist, if computed distance is smaller */
	for(i=0;i<numvertices;i++){
		vec_signed_dist->GetValue(&lsf_old, element->vertices[i]->Sid());
		if(xIsNan<IssmDouble>(lsf_old) || fabs(signed_dist[i])<fabs(lsf_old))
			vec_signed_dist->SetValue(element->vertices[i]->Sid(),signed_dist[i],INS_VAL);
	}

	xDelete<IssmDouble>(lsf);
	xDelete<IssmDouble>(sign_lsf);
	xDelete<IssmDouble>(signed_dist);
	xDelete<IssmDouble>(xyz_list);
	xDelete<IssmDouble>(xyz_list_zero);
}/*}}}*/
void           LevelsetAnalysis::UpdateConstraints(FemModel* femmodel){/*{{{*/

	/*Intermediaries*/
	int         calvinglaw;
	IssmDouble  min_thickness,thickness; 
	femmodel->parameters->FindParam(&calvinglaw,CalvingLawEnum);

	if(calvinglaw==CalvingMinthicknessEnum){

		/*Get minimum thickness threshold*/
		femmodel->parameters->FindParam(&min_thickness,CalvingMinthicknessEnum);

		/*Loop over all elements of this partition*/
		for(int i=0;i<femmodel->elements->Size();i++){
			Element* element  = xDynamicCast<Element*>(femmodel->elements->GetObjectByOffset(i));

			int      numnodes = element->GetNumberOfNodes();	
			Gauss*   gauss    = element->NewGauss();
			Input*   H_input  = element->GetInput(ThicknessEnum); _assert_(H_input);

			/*Potentially constrain nodes of this element*/
			for(int in=0;in<numnodes;in++){
				gauss->GaussNode(element->GetElementType(),in);
				Node* node=element->GetNode(in);
				H_input->GetInputValue(&thickness,gauss);
				if(thickness<min_thickness){
					node->ApplyConstraint(0,+1.);
				}
				else {
					/* no ice, set no spc */
					node->DofInFSet(0); 
				}
			}
			delete gauss;
		}
	}
	/*Default, do nothing*/
	return;
}/*}}}*/
