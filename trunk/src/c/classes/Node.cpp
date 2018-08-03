/*!\file Node.c
 * \brief: implementation of the Node object
 */

/*Include files: {{{*/
#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#include "./classes.h"
#include "shared/shared.h"
#include "modules/ModelProcessorx/ModelProcessorx.h"
#include "../analyses/analyses.h"
/*}}}*/

/*Node constructors and destructors:*/
Node::Node(){/*{{{*/
	this->approximation=0;
}
/*}}}*/
Node::Node(int node_id,int node_sid,int node_lid,int io_index, IoModel* iomodel,int analysis_enum,int in_approximation){/*{{{*/

	/*Intermediary*/
	int k,l;
	int *doftypes = NULL;

	/*id: */
	this->id            = node_id;
	this->sid           = node_sid;
	this->lid           = node_lid;
	this->analysis_enum = analysis_enum;

	/*Initialize coord_system: Identity matrix by default*/
	for(k=0;k<3;k++) for(l=0;l<3;l++) this->coord_system[k][l]=0.0;
	for(k=0;k<3;k++) this->coord_system[k][k]=1.0;

	/*indexing:*/
	this->indexingupdate = true;

	Analysis* analysis = EnumToAnalysis(analysis_enum);
	int numdofs        = analysis->DofsPerNode(&doftypes,iomodel->domaintype,in_approximation);
	indexing.Init(numdofs,doftypes);
	xDelete<int>(doftypes);
	delete analysis;

	if(analysis_enum==StressbalanceAnalysisEnum)
	 this->approximation=in_approximation;
	else
	 this->approximation=0;

	/*Stressbalance Horiz*/
	if(analysis_enum==StressbalanceAnalysisEnum){

		/*Coordinate system provided, convert to coord_system matrix*/
		_assert_(iomodel->Data("md.stressbalance.referential")); 
		XZvectorsToCoordinateSystem(&this->coord_system[0][0],&iomodel->Data("md.stressbalance.referential")[io_index*6]);
		_assert_(sqrt( coord_system[0][0]*coord_system[0][0] + coord_system[1][0]*coord_system[1][0]) >1.e-4);

		if(iomodel->domaintype!=Domain2DhorizontalEnum && iomodel->domaintype!=Domain3DsurfaceEnum){
			/*We have a  3d mesh, we may have collapsed elements, hence dead nodes. Freeze them out: */
			_assert_(iomodel->Data("md.mesh.vertexonbase")); 
			_assert_(iomodel->Data("md.flowequation.vertex_equation"));
			if(in_approximation==SSAApproximationEnum && !reCast<int>(iomodel->Data("md.mesh.vertexonbase")[io_index])){
				this->HardDeactivate();
			}
			if(in_approximation==L1L2ApproximationEnum && !reCast<int>(iomodel->Data("md.mesh.vertexonbase")[io_index])){
				this->HardDeactivate();
			}
			if(in_approximation==SSAHOApproximationEnum && reCast<int>(iomodel->Data("md.flowequation.borderSSA")[io_index])){
				if(!reCast<int>(iomodel->Data("md.mesh.vertexonbase")[io_index])){
					this->HardDeactivate();
				}
			}
			if(in_approximation==SSAFSApproximationEnum && reCast<int>(iomodel->Data("md.flowequation.borderSSA")[io_index])){
				if(!reCast<int>(iomodel->Data("md.mesh.vertexonbase")[io_index])){
					for(k=0;k<=1;k++) this->FreezeDof(k);
				}
			}
		}
		/*spc all nodes on SIA*/
		if(in_approximation==SIAApproximationEnum){
			this->HardDeactivate();
		}
	}

	/*2d solutions in 3d, we need to constrain all the nodes that are not on base*/
	if(
				analysis_enum==FreeSurfaceBaseAnalysisEnum || 
				analysis_enum==MasstransportAnalysisEnum || 
				analysis_enum==MeltingAnalysisEnum || 
				analysis_enum==L2ProjectionBaseAnalysisEnum || 
				analysis_enum==BalancethicknessAnalysisEnum ||
				analysis_enum==HydrologyDCInefficientAnalysisEnum ||
				analysis_enum==HydrologyDCEfficientAnalysisEnum ||
				analysis_enum==LevelsetAnalysisEnum
				){
		if(iomodel->domaintype!=Domain2DhorizontalEnum & iomodel->domaintype!=Domain3DsurfaceEnum){
			/*On a 3d mesh, we may have collapsed elements, hence dead nodes. Freeze them out: */
			_assert_(iomodel->Data("md.mesh.vertexonbase"));
			if(!(reCast<bool>(iomodel->Data("md.mesh.vertexonbase")[io_index]))){
				this->HardDeactivate();
			}
		}
	}
	if(
				analysis_enum==FreeSurfaceTopAnalysisEnum
				){
		if(iomodel->domaintype!=Domain2DhorizontalEnum){
			/*On a 3d mesh, we may have collapsed elements, hence dead nodes. Freeze them out: */
			_assert_(iomodel->Data("md.mesh.vertexonsurface"));
			if(!(reCast<bool>(iomodel->Data("md.mesh.vertexonsurface")[io_index]))){
				this->HardDeactivate();
			}
		}
	}

}
/*}}}*/
Node::~Node(){/*{{{*/
	return;
}
/*}}}*/
Object* Node::copy(void){/*{{{*/

	int k,l;

	/*output: */
	Node* output=NULL;

	/*initalize output: */
	output=new Node();

	/*id: */
	output->id  = this->id;
	output->sid = this->sid;
	output->lid = this->lid;
	output->analysis_enum = this->analysis_enum;
	output->approximation = this->approximation;

	/*Initialize coord_system: */
	for(k=0;k<3;k++) for(l=0;l<3;l++) output->coord_system[k][l]=this->coord_system[k][l];

	/*indexing:*/
	output->indexingupdate = this->indexingupdate;
	output->indexing.copy(this->indexing);

	return (Object*)output; 
}
/*}}}*/
void Node::Marshall(char** pmarshalled_data,int* pmarshalled_data_size, int marshall_direction){ /*{{{*/

	MARSHALLING_ENUM(NodeEnum);
	MARSHALLING(id);
	MARSHALLING(sid);
	MARSHALLING(lid);
	MARSHALLING(indexingupdate);
	indexing.Marshall(pmarshalled_data,pmarshalled_data_size,marshall_direction);
	MARSHALLING(analysis_enum);
	MARSHALLING_ARRAY(coord_system,IssmDouble,9);

}
/*}}}*/

/*Object virtual functions definitions:*/
void Node::DeepEcho(void){/*{{{*/

	_printf_("Node:\n");
	_printf_("   id: " << id << "\n");
	_printf_("   sid: " << sid << "\n");
	_printf_("   analysis_enum: " << EnumToStringx(analysis_enum) << "\n");
	_printf_("   approximation: " << EnumToStringx(approximation) << "\n");
	_printf_("   indexingupdate: " << indexingupdate << "\n");
	indexing.DeepEcho();

}
/*}}}*/
void Node::Echo(void){/*{{{*/

	_printf_("Node:\n");
	_printf_("   id : " << id << "\n");
	_printf_("   sid: " << sid << "\n");
	_printf_("   analysis_enum: " << EnumToStringx(analysis_enum) << "\n");
	_printf_("   approximation: " << EnumToStringx(approximation) << "\n");
	_printf_("   indexingupdate: " << indexingupdate << "\n");
	indexing.Echo();

}
/*}}}*/
int  Node::Id(void){ return id; }/*{{{*/
/*}}}*/
int  Node::ObjectEnum(void){/*{{{*/

	return NodeEnum;

}
/*}}}*/

/*Node management:*/
void Node::GetCoordinateSystem(IssmDouble* coord_system_out){/*{{{*/

	/*Copy coord_system*/
	for(int k=0;k<3;k++) for(int l=0;l<3;l++) coord_system_out[3*k+l]=this->coord_system[k][l];

}
/*}}}*/
int  Node::GetDof(int dofindex,int setenum){/*{{{*/

	_assert_(!this->indexingupdate);
	if(setenum==GsetEnum){
		_assert_(dofindex>=0 && dofindex<indexing.gsize);
		return indexing.gdoflist[dofindex];
	}
	else if(setenum==FsetEnum){
		_assert_(dofindex>=0 && dofindex<indexing.fsize);
		return indexing.fdoflist[dofindex];
	}
	else if(setenum==SsetEnum){
		_assert_(dofindex>=0 && dofindex<indexing.ssize);
		return indexing.sdoflist[dofindex];
	}
	else _error_("set of enum type " << EnumToStringx(setenum) << " not supported yet!");

} /*}}}*/
void Node::GetDofList(int* outdoflist,int approximation_enum,int setenum){/*{{{*/
	int i;
	int count=0;
	int count2=0;

	_assert_(!this->indexingupdate);

	if(approximation_enum==NoneApproximationEnum){
		if(setenum==GsetEnum)for(i=0;i<this->indexing.gsize;i++) outdoflist[i]=indexing.gdoflist[i];
		if(setenum==FsetEnum)for(i=0;i<this->indexing.fsize;i++) outdoflist[i]=indexing.fdoflist[i];
		if(setenum==SsetEnum)for(i=0;i<this->indexing.ssize;i++) outdoflist[i]=indexing.sdoflist[i];
	}
	else{

		if(setenum==GsetEnum){
			if(indexing.doftype){
				count=0;
				for(i=0;i<this->indexing.gsize;i++){
					if(indexing.doftype[i]==approximation_enum){
						outdoflist[count]=indexing.gdoflist[i];
						count++;
					}
				}
				_assert_(count); //at least one dof should be the approximation requested
			}
			else for(i=0;i<this->indexing.gsize;i++) outdoflist[i]=indexing.gdoflist[i];
		}
		else if(setenum==FsetEnum){
			if(indexing.doftype){
				count=0;
				count2=0;
				for(i=0;i<this->indexing.gsize;i++){
					if(indexing.f_set[i]){
						if(indexing.doftype[i]==approximation_enum){
							outdoflist[count]=indexing.fdoflist[count2];
							count++;
						}
						count2++;
					}
				}
			}
			else for(i=0;i<this->indexing.fsize;i++) outdoflist[i]=indexing.fdoflist[i];
		}
		else if(setenum==SsetEnum){
			if(indexing.doftype){
				count=0;
				count2=0;
				for(i=0;i<this->indexing.gsize;i++){
					if(indexing.s_set[i]){
						if(indexing.doftype[i]==approximation_enum){
							outdoflist[count]=indexing.sdoflist[count2];
							count++;
						}
						count2++;
					}
				}
			}
			else for(i=0;i<this->indexing.ssize;i++) outdoflist[i]=indexing.sdoflist[i];
		}
		else _error_("set of enum type " << EnumToStringx(setenum) << " not supported yet!");
	}
}
/*}}}*/
void Node::GetLocalDofList(int* outdoflist,int approximation_enum,int setenum){/*{{{*/
	int i;
	int count=0;
	int count2=0;

	_assert_(!this->indexingupdate);

	if(approximation_enum==NoneApproximationEnum){
		if(setenum==GsetEnum)for(i=0;i<this->indexing.gsize;i++) outdoflist[i]=i;
		else if(setenum==FsetEnum){
			count=0;
			for(i=0;i<this->indexing.gsize;i++){
				if(indexing.f_set[i]){
					outdoflist[count]=i;
					count++;
				}
			}
		}
		else if(setenum==SsetEnum){
			count=0;
			for(i=0;i<this->indexing.gsize;i++){
				if(indexing.s_set[i]){
					outdoflist[count]=i;
					count++;
				}
			}
		}
		else _error_("set of enum type " << EnumToStringx(setenum) << " not supported yet!");
	}
	else{

		if(setenum==GsetEnum){
			if(indexing.doftype){
				count=0;
				for(i=0;i<this->indexing.gsize;i++){
					if(indexing.doftype[i]==approximation_enum){
						outdoflist[count]=count;
						count++;
					}
				}
				_assert_(count);
			}
			else for(i=0;i<this->indexing.gsize;i++) outdoflist[i]=i;
		}
		else if(setenum==FsetEnum){

			if(indexing.doftype){
				count=0;
				count2=0;
				for(i=0;i<this->indexing.gsize;i++){
					if(indexing.doftype[i]==approximation_enum){
						if(indexing.f_set[i]){
							outdoflist[count]=count2;
							count++;
						}
						count2++;
					}
				}
				_assert_(count2);
			}
			else{

				count=0;
				for(i=0;i<this->indexing.gsize;i++){
					if(indexing.f_set[i]){
						outdoflist[count]=i;
						count++;
					}
				}
			}
		}
		else if(setenum==SsetEnum){
			if(indexing.doftype){
				count=0;
				count2=0;
				for(i=0;i<this->indexing.gsize;i++){
					if(indexing.doftype[i]==approximation_enum){
						if(indexing.s_set[i]){
							outdoflist[count]=count2;
							count++;
						}
						count2++;
					}
				}
				_assert_(count2);
			}
			else{
				count=0;
				for(i=0;i<this->indexing.gsize;i++){
					if(indexing.s_set[i]){
						outdoflist[count]=i;
						count++;
					}
				}
			}
		}
		else _error_("set of enum type " << EnumToStringx(setenum) << " not supported yet!");
	}
}
/*}}}*/
bool Node::InAnalysis(int in_analysis_enum){/*{{{*/
	if (in_analysis_enum==this->analysis_enum) return true;
	else return false;
}
/*}}}*/
int  Node::Lid(void){/*{{{*/
	return lid; 
}
/*}}}*/
int  Node::Sid(void){/*{{{*/
	return sid; 
}
/*}}}*/

/*Node numerics:*/
void Node::Activate(void){/*{{{*/

	if(!IsActive() && !this->indexing.freeze){
		this->indexingupdate = true;
		indexing.Activate();
	}

}
/*}}}*/
void Node::ApplyConstraint(int dof,IssmDouble value){/*{{{*/

	/*Dof should be added in the s set, describing which 
	 * dofs are constrained to a certain value (dirichlet boundary condition*/
	DofInSSet(dof);
	this->indexing.svalues[dof]=value;
}
/*}}}*/
void Node::CreateNodalConstraints(Vector<IssmDouble>* ys){/*{{{*/

	int i;
	IssmDouble* values=NULL;
	int count;

	/*Recover values for s set and plug them in constraints vector: */
	if(this->indexing.ssize){
		values=xNew<IssmDouble>(this->indexing.ssize);
		count=0;
		for(i=0;i<this->indexing.gsize;i++){
			if(this->indexing.s_set[i]){
				values[count]=this->indexing.svalues[i];
				_assert_(!xIsNan<IssmDouble>(values[count]));
				count++;
			}
		}

		/*Add values into constraint vector: */
		ys->SetValues(this->indexing.ssize,this->indexing.sdoflist,values,INS_VAL);
	}

	/*Free ressources:*/
	xDelete<IssmDouble>(values);

}
/*}}}*/
void Node::Deactivate(void){/*{{{*/

	if(IsActive() && !this->indexing.freeze){
		this->indexingupdate = true;
		indexing.Deactivate();
	}

}
/*}}}*/
void Node::DofInFSet(int dof){/*{{{*/

	/*Put dof for this node into the f set (ie, this dof will NOT be constrained 
	 * to a fixed value during computations. Only do this for active nodes. */
	_assert_(dof<this->indexing.gsize);
	_assert_(this->indexing.active);

	if(this->indexing.f_set[dof] == 0){
		if(this->indexing.freeze) _error_("Cannot change dof of frozen node");
		this->indexingupdate = true;
		this->indexing.f_set[dof]=1; 
		this->indexing.s_set[dof]=0;
	}
}
/*}}}*/
void Node::DofInSSet(int dof){/*{{{*/

	/*Put dof for this node into the s set (ie, this dof will be constrained 
	 * to a fixed value during computations. */
	_assert_(dof<this->indexing.gsize);

	if(this->indexing.f_set[dof] == 1){
		//if(this->indexing.freeze) _error_("Cannot change dof of frozen node");
		this->indexingupdate = true;
		this->indexing.f_set[dof]=0; //n splits into f (for which we solve) and s (single point constraints)
		this->indexing.s_set[dof]=1;
	}
}
/*}}}*/
void Node::FreezeDof(int dof){/*{{{*/

	DofInSSet(dof); //with 0 displacement for this dof.
	//FIXME: for now we don't want this element to change so we use freeze
	this->indexing.freeze =true;

}
/*}}}*/
int  Node::GetApproximation(){/*{{{*/

	return approximation;
}
/*}}}*/
void  Node::SetApproximation(int in_approximation){/*{{{*/
	
	this->approximation = in_approximation;
}
/*}}}*/
int  Node::GetNumberOfDofs(int approximation_enum,int setenum){/*{{{*/

	/*Get number of degrees of freedom in a node, for a certain set (g,f or s-set)
	 *and for a certain approximation type: */

	int i;
	int numdofs=0;

	if(approximation_enum==NoneApproximationEnum){
		if      (setenum==GsetEnum) numdofs=this->indexing.gsize;
		else if (setenum==FsetEnum) numdofs=this->indexing.fsize;
		else if (setenum==SsetEnum) numdofs=this->indexing.ssize;
		else _error_("set of enum type " << EnumToStringx(setenum) << " not supported yet!");
	}
	else{
		if(setenum==GsetEnum){
			if(this->indexing.doftype){
				numdofs=0;
				for(i=0;i<this->indexing.gsize;i++){
					if(this->indexing.doftype[i]==approximation_enum) numdofs++;
				}
			}
			else numdofs=this->indexing.gsize;
		}
		else if (setenum==FsetEnum){
			if(this->indexing.doftype){
				numdofs=0;
				for(i=0;i<this->indexing.gsize;i++){
					if((this->indexing.doftype[i]==approximation_enum) && (this->indexing.f_set[i])) numdofs++;
				}
			}
			else numdofs=this->indexing.fsize;
		}
		else if (setenum==SsetEnum){
			if(this->indexing.doftype){
			numdofs=0;
				for(i=0;i<this->indexing.gsize;i++){
					if((this->indexing.doftype[i]==approximation_enum) && (this->indexing.s_set[i])) numdofs++;
				}
			}
			else numdofs=this->indexing.ssize;
		}
		else _error_("set of enum type " << EnumToStringx(setenum) << " not supported yet!");
	}
	return numdofs;
}
/*}}}*/
void Node::HardDeactivate(void){/*{{{*/

	this->indexing.Deactivate();
	this->indexing.freeze =true;

}
/*}}}*/
bool Node::IsActive(void){/*{{{*/

	return indexing.active;

}
/*}}}*/
int  Node::IsClone(){/*{{{*/

	return indexing.clone;

}
/*}}}*/
void Node::ReindexingDone(void){/*{{{*/

	this->indexingupdate = false;

}
/*}}}*/
void Node::RelaxConstraint(int dof){/*{{{*/

	/*Dof should be added to the f-set, and taken out of the s-set:*/
	DofInFSet(dof);
	this->indexing.svalues[dof]=0.;
}
/*}}}*/
bool Node::RequiresDofReindexing(void){/*{{{*/

	return this->indexingupdate;

}
/*}}}*/
void Node::VecMerge(Vector<IssmDouble>* ug, IssmDouble* vector_serial,int setenum){/*{{{*/

	IssmDouble *values  = NULL;
	int        *indices = NULL;
	int         count   = 0;
	int         i;

	if(setenum==FsetEnum){
		if(this->indexing.fsize){
			indices=xNew<int>(this->indexing.fsize);
 			values=xNew<IssmDouble>(this->indexing.fsize);

			for(i=0;i<this->indexing.gsize;i++){
				if(this->indexing.f_set[i]){
					_assert_(vector_serial);
					values[count]=vector_serial[this->indexing.fdoflist[count]];
					indices[count]=this->indexing.gdoflist[i];
					count++;
				}
			}

			/*Add values into ug: */
			ug->SetValues(this->indexing.fsize,indices,values,INS_VAL);
		}
	}
	else if(setenum==SsetEnum){
		if(this->indexing.ssize){
			indices=xNew<int>(this->indexing.ssize);
			values=xNew<IssmDouble>(this->indexing.ssize);

			for(i=0;i<this->indexing.gsize;i++){
				if(this->indexing.s_set[i]){
					_assert_(vector_serial);
					values[count]=vector_serial[this->indexing.sdoflist[count]];
					indices[count]=this->indexing.gdoflist[i];
					count++;
				}
			}

			/*Add values into ug: */
			ug->SetValues(this->indexing.ssize,indices,values,INS_VAL);
		}
	}
	else _error_("VecMerge can only merge from the s or f-set onto the g-set!");

	/*Free ressources:*/
	xDelete<IssmDouble>(values);
	xDelete<int>(indices);
}
/*}}}*/
void Node::VecReduce(Vector<IssmDouble>* vector, IssmDouble* ug_serial,int setenum){/*{{{*/

	IssmDouble* values=NULL;
	int     count=0;
	int     i;

	if(setenum==FsetEnum){
		if(this->indexing.fsize){
 			values=xNew<IssmDouble>(this->indexing.fsize);

			for(i=0;i<this->indexing.gsize;i++){
				if(this->indexing.f_set[i]){
					_assert_(ug_serial);
					values[count]=ug_serial[this->indexing.gdoflist[i]];
					count++;
				}
			}

			/*Add values into ug: */
			vector->SetValues(this->indexing.fsize,this->indexing.fdoflist,values,INS_VAL);
		}
	}
	else if(setenum==SsetEnum){
		if(this->indexing.ssize){
			values=xNew<IssmDouble>(this->indexing.ssize);

			for(i=0;i<this->indexing.gsize;i++){
				if(this->indexing.s_set[i]){
					_assert_(ug_serial);
					values[count]=ug_serial[this->indexing.gdoflist[i]];
					count++;
				}
			}

			/*Add values into ug: */
			vector->SetValues(this->indexing.ssize,this->indexing.sdoflist,values,INS_VAL);
		}
	}
	else _error_("VecReduce can only merge from the s or f-set onto the g-set!");

	/*Free ressources:*/
	xDelete<IssmDouble>(values);
}
/*}}}*/

/* indexing routines:*/
void Node::DistributeDofs(int* pdofcount,int setenum){/*{{{*/

	int i;
	int dofcount;

	dofcount=*pdofcount;

	/*Initialize: */
	if(setenum==FsetEnum) this->indexing.InitSet(setenum);
	if(setenum==SsetEnum) this->indexing.InitSet(setenum);

	/*For clone nodfs, don't distribute dofs, we will get them from another cpu in UpdateCloneDofs!*/
	if(indexing.clone){
		return;
	}

	/*This node should distribute dofs for setenum set (eg, f_set or s_set), go ahead: */
	if(setenum==GsetEnum){
		for(i=0;i<this->indexing.gsize;i++){
			indexing.gdoflist[i]=dofcount+i;
		}
		dofcount+=this->indexing.gsize;
	}
	else if(setenum==FsetEnum){
		for(i=0;i<this->indexing.fsize;i++){
			indexing.fdoflist[i]=dofcount+i;
		}
		dofcount+=this->indexing.fsize;
	}
	else if(setenum==SsetEnum){
		for(i=0;i<this->indexing.ssize;i++){
			indexing.sdoflist[i]=dofcount+i;
		}
		dofcount+=this->indexing.ssize;
	}
	else _error_("set of enum type " << EnumToStringx(setenum) << " not supported yet!");

	/*Assign output pointers: */
	*pdofcount=dofcount;
}
/*}}}*/
void Node::OffsetDofs(int dofcount,int setenum){/*{{{*/

	int i;

	if(indexing.clone){
		/*This node is a clone, don't off_set the dofs!: */
		return;
	}

	/*This node should off_set the dofs, go ahead: */
	if(setenum==GsetEnum){
		for(i=0;i<this->indexing.gsize;i++) indexing.gdoflist[i]+=dofcount;
	}
	else if(setenum==FsetEnum){
		for(i=0;i<this->indexing.fsize;i++) indexing.fdoflist[i]+=dofcount;
	}
	else if(setenum==SsetEnum){
		for(i=0;i<this->indexing.ssize;i++) indexing.sdoflist[i]+=dofcount;
	}
	else _error_("set of enum type " << EnumToStringx(setenum) << " not supported yet!");
}
/*}}}*/
void Node::SetClone(int* minranks){/*{{{*/

	int my_rank;

	/*recover my_rank:*/
	my_rank=IssmComm::GetRank();

	if (minranks[sid]==my_rank){
		indexing.clone=false;
	}
	else{
		/*!there is a cpu with lower rank that has the same node, 
		therefore, I am a clone*/
		indexing.clone=true;	
	}
}
/*}}}*/
void Node::ShowTrueDofs(int* truedofs, int ncols,int setenum){/*{{{*/

	int j;

	/*Are we a clone? : */
	if(indexing.clone) return;

	/*Ok, we are not a clone, just plug our dofs into truedofs: */
	switch(setenum){
		case GsetEnum:
			for(j=0;j<this->indexing.gsize;j++) truedofs[ncols*sid+j]=indexing.gdoflist[j];
			break;
		case FsetEnum:
			for(j=0;j<this->indexing.fsize;j++) truedofs[ncols*sid+j]=indexing.fdoflist[j];
			break;
		case SsetEnum:
			for(j=0;j<this->indexing.ssize;j++) truedofs[ncols*sid+j]=indexing.sdoflist[j];
			break;
		default:
			_error_("set of enum type " << EnumToStringx(setenum) << " not supported yet!");
	}

}
/*}}}*/
void Node::UpdateCloneDofs(int* alltruedofs,int ncols,int setenum){/*{{{*/

	int j;

	/*If we are not a clone, don't update, we already have dofs!: */
	if(!indexing.clone)return;

	/*Ok, we are a clone node, but we did not create the dofs for this node.
	 *Therefore, our doflist is garbage right now. Go pick it up in the alltruedofs: */
	switch(setenum){
		case GsetEnum:
			for(j=0;j<this->indexing.gsize;j++) indexing.gdoflist[j]=alltruedofs[ncols*sid+j];
			break;
		case FsetEnum:
			for(j=0;j<this->indexing.fsize;j++) indexing.fdoflist[j]=alltruedofs[ncols*sid+j];
			break;
		case SsetEnum:
			for(j=0;j<this->indexing.ssize;j++) indexing.sdoflist[j]=alltruedofs[ncols*sid+j];
			break;
		default:
			_error_("set of enum type " << EnumToStringx(setenum) << " not supported yet!");
	}
}
/*}}}*/

/*Methods inherent to Node: */
int* GetGlobalDofList(Node** nodes,int numnodes,int setenum,int approximation){/*{{{*/

	int  i,numdof,count;
	int* ndof_list=NULL;
	int *doflist = NULL;

	if(numnodes){

		/*Allocate:*/
		ndof_list=xNew<int>(numnodes);

		/*First, figure out size of doflist: */
		numdof=0;
		for(i=0;i<numnodes;i++){
			ndof_list[i]=nodes[i]->GetNumberOfDofs(approximation,setenum);
			numdof+=ndof_list[i];
		}

		if(numdof){
			/*Allocate: */
			doflist=xNew<int>(numdof);

			/*Populate: */
			count=0;
			for(i=0;i<numnodes;i++){
				nodes[i]->GetDofList(&doflist[count],approximation,setenum);
				count+=ndof_list[i];
			}
		}
		else doflist=NULL;
	}
	/*Free ressources:*/
	xDelete<int>(ndof_list);

	return doflist;
}
/*}}}*/
int* GetLocalDofList(Node** nodes,int numnodes,int setenum,int approximation){ /*{{{*/

	int  i,j,count,numdof,numgdof;
	int* ndof_list=NULL;
	int* ngdof_list_cumulative=NULL;
	int *doflist = NULL;

	if(numnodes){
		/*allocate: */
		ndof_list=xNew<int>(numnodes);
		ngdof_list_cumulative=xNew<int>(numnodes);

		/*Get number of dofs per node, and total for this given set*/
		numdof=0;
		numgdof=0;
		for(i=0;i<numnodes;i++){

			/*Cumulative list= number of dofs before node i*/
			ngdof_list_cumulative[i]=numgdof;

			/*Number of dofs for node i for given set and for the g set*/
			ndof_list[i]=nodes[i]->GetNumberOfDofs(approximation,setenum);
			numgdof    +=nodes[i]->GetNumberOfDofs(approximation,GsetEnum);
			numdof     +=ndof_list[i];
		}

		if(numdof){
			/*Allocate: */
			doflist=xNew<int>(numdof);

			/*Populate: */
			count=0;
			for(i=0;i<numnodes;i++){
				nodes[i]->GetLocalDofList(&doflist[count],approximation,setenum);
				count+=ndof_list[i];
			}

			/*We now have something like: [0 1 0 2 1 2]. Offset by gsize, to get something like: [0 1 2 4 6 7]:*/
			count=0;
			for(i=0;i<numnodes;i++){
				for(j=0;j<ndof_list[i];j++){
					doflist[count+j]+=ngdof_list_cumulative[i];
				}
				count+=ndof_list[i];
			}
		}
		else doflist=NULL;
	}

	/*Free ressources:*/
	xDelete<int>(ndof_list);
	xDelete<int>(ngdof_list_cumulative);

	/*CLean-up and return*/
	return doflist;
}
/*}}}*/
int GetNumberOfDofs(Node** nodes,int numnodes,int setenum,int approximation){/*{{{*/

	/*output: */
	int numberofdofs=0;

	for(int i=0;i<numnodes;i++){
		numberofdofs+=nodes[i]->GetNumberOfDofs(approximation,setenum);
	}

	return numberofdofs;
}
/*}}}*/
