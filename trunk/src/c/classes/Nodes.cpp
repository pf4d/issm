/*
 * \file Nodes.cpp
 * \brief: Implementation of Nodes class, derived from DataSet class.
 */

/*Headers: {{{*/
#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#include "../shared/io/Comm/IssmComm.h"
#include "./Nodes.h"
#include "./Node.h"

using namespace std;
/*}}}*/

/*Object constructors and destructor*/
Nodes::Nodes(){/*{{{*/
	enum_type=NodesEnum;
	return;
}
/*}}}*/
Nodes::~Nodes(){/*{{{*/
	return;
}
/*}}}*/

/*Numerics*/
void  Nodes::DistributeDofs(int analysis_type,int setenum){/*{{{*/

	int  i;
	int  dofcount=0;
	int  maxdofspernode=0;
	int* alldofcount=NULL;
	int* truedofs=NULL;
	int* alltruedofs=NULL;
	int  numnodes=0;

	/*recover my_rank:*/
	int my_rank   = IssmComm::GetRank();
	int num_procs = IssmComm::GetSize();

	/*some check: */
	_assert_(setenum==GsetEnum || setenum==FsetEnum || setenum==SsetEnum);

	/*Go through objects, and distribute dofs locally, from 0 to numberofdofsperobject*/
	for(i=0;i<this->Size();i++){
		Node* node=xDynamicCast<Node*>(this->GetObjectByOffset(i));

		/*Check that this node corresponds to our analysis currently being carried out: */
		if(node->InAnalysis(analysis_type)){
			node->DistributeDofs(&dofcount,setenum);
		}
	}

	/* Now every object has distributed dofs, but locally, and with a dof count starting from 
	 * 0. This means the dofs between all the cpus are not unique. We now offset the dofs of eache
	 * cpus by the total last dofs of the previus cpu, starting from 0.
	 * First: get number of dofs for each cpu*/
	alldofcount=xNew<int>(num_procs);
	ISSM_MPI_Gather(&dofcount,1,ISSM_MPI_INT,alldofcount,1,ISSM_MPI_INT,0,IssmComm::GetComm());
	ISSM_MPI_Bcast(alldofcount,num_procs,ISSM_MPI_INT,0,IssmComm::GetComm());

	/* Every cpu should start its own dof count at the end of the dofcount from cpu-1*/
	dofcount=0;
	for(i=0;i<my_rank;i++){
		dofcount+=alldofcount[i];
	}
	for(i=0;i<this->Size();i++){
		/*Check that this node corresponds to our analysis currently being carried out: */
		Node* node=xDynamicCast<Node*>(this->GetObjectByOffset(i));
		if (node->InAnalysis(analysis_type)){
			node->OffsetDofs(dofcount,setenum);
		}
	}

	/* Finally, remember that cpus may have skipped some objects, because they were clones. For every 
	 * object that is not a clone, tell them to show their dofs, so that later on, they can get picked 
	 * up by their clones: */
	maxdofspernode=this->MaxNumDofs(analysis_type,setenum);
	numnodes=this->NumberOfNodes(analysis_type);
	if(numnodes*maxdofspernode){
		truedofs=   xNewZeroInit<int>(numnodes*maxdofspernode); //initialize to 0, so that we can pick up the max
		alltruedofs=xNewZeroInit<int>(numnodes*maxdofspernode);
	}

	for(i=0;i<this->Size();i++){
		Node* node=xDynamicCast<Node*>(this->GetObjectByOffset(i));
		if (node->InAnalysis(analysis_type)){
			node->ShowTrueDofs(truedofs,maxdofspernode,setenum);//give maxdofspernode, column size, so that nodes can index into truedofs
		}
	}

	ISSM_MPI_Allreduce((void*)truedofs,(void*)alltruedofs,numnodes*maxdofspernode,ISSM_MPI_INT,ISSM_MPI_MAX,IssmComm::GetComm());

	/* Now every cpu knows the true dofs of everyone else that is not a clone*/
	for(i=0;i<this->Size();i++){
		Node* node=xDynamicCast<Node*>(this->GetObjectByOffset(i));
		if (node->InAnalysis(analysis_type)){
			node->UpdateCloneDofs(alltruedofs,maxdofspernode,setenum);
		}
	}

	/*Update indexingupdateflag*/
	for(i=0;i<this->Size();i++){
		Node* node=xDynamicCast<Node*>(this->GetObjectByOffset(i));
		if (node->InAnalysis(analysis_type)){
			node->ReindexingDone();
		}
	}

	/* Free ressources: */
	xDelete<int>(alldofcount);
	xDelete<int>(truedofs);
	xDelete<int>(alltruedofs);
}
/*}}}*/
void  Nodes::FlagClones(int analysis_type){/*{{{*/

	int i;
	int num_procs;
	int numnodes;

	/*recover num_procs: */
	num_procs=IssmComm::GetSize();

	/*Figure out number of nodes for this analysis: */
	numnodes=this->NumberOfNodes(analysis_type);

	/*Allocate ranks: */
	int* ranks    = xNew<int>(numnodes);
	int* minranks = xNew<int>(numnodes);
	for(i=0;i<numnodes;i++)ranks[i]=num_procs; //no cpu can have rank num_procs. This is the maximum limit.

	/*Now go through all our objects and ask them to report to who they belong (which rank): */
	Ranks(ranks,analysis_type);

	/*We need to take the minimum rank for each vertex, and every cpu needs to get that result. That way, 
	 * when we start building the dof list for all vertexs, a cpu can check whether its vertex already has been 
	 * dealt with by another cpu. We take the minimum because we are going to manage dof assignment in increasing 
	 * order of cpu rank. This is also why we initialized this array to num_procs.*/
	ISSM_MPI_Allreduce((void*)ranks,(void*)minranks,numnodes,ISSM_MPI_INT,ISSM_MPI_MIN,IssmComm::GetComm());

	/*Now go through all objects, and use minranks to flag which objects are cloned: */
	for(i=0;i<this->Size();i++){

		Node* node=xDynamicCast<Node*>(this->GetObjectByOffset(i));

		/*Check that this node corresponds to our analysis currently being carried out: */
		if (node->InAnalysis(analysis_type)){

			/*For this object, decide whether it is a clone: */
			node->SetClone(minranks);
		}
	}

	/*Free ressources: */
	xDelete<int>(ranks); 
	xDelete<int>(minranks);

}
/*}}}*/
int   Nodes::MaximumId(){/*{{{*/

	int max=-1;
	int id,allmax;

	/*Now go through all nodes, and get how many dofs they own, unless they are clone nodes: */
	if(!sorted){
		for(int i=0;i<this->Size();i++){
			Node* node=xDynamicCast<Node*>(this->GetObjectByOffset(i));
			id=node->Id();
			if(id>max)max=id;
		}
	}
	else{
		if(this->Size()==0){
			max = 0;
		}
		else{
			Node* node=xDynamicCast<Node*>(this->GetObjectByOffset(this->Size()-1));
			max = node->Id();
		}
	}

	/*Grab max of all cpus: */
	ISSM_MPI_Allreduce((void*)&max,(void*)&allmax,1,ISSM_MPI_INT,ISSM_MPI_MAX,IssmComm::GetComm());
	max=allmax;

	return max;
}
/*}}}*/
int   Nodes::MaxNumDofs(int analysis_type,int setenum){/*{{{*/

	int max=0;
	int allmax,numdofs;

	/*Now go through all nodes, and get how many dofs they own, unless they are clone nodes: */
	for(int i=0;i<this->Size();i++){

		Node* node=xDynamicCast<Node*>(this->GetObjectByOffset(i));

		/*Check that this node corresponds to our analysis currently being carried out: */
		if (node->InAnalysis(analysis_type)){

			numdofs=node->GetNumberOfDofs(NoneApproximationEnum,setenum);
			if(numdofs>max)max=numdofs;
		}
	}

	/*Grab max of all cpus: */
	ISSM_MPI_Allreduce((void*)&max,(void*)&allmax,1,ISSM_MPI_INT,ISSM_MPI_MAX,IssmComm::GetComm());
	max=allmax;

	return max;
}
/*}}}*/
int   Nodes::NumberOfDofs(int analysis_type,int setenum){/*{{{*/

	int   allnumdofs;

	/*Get number of dofs on current cpu (excluding clones)*/
	int numdofs=this->NumberOfDofsLocal(analysis_type,setenum);

	/*Gather from all cpus: */
	ISSM_MPI_Allreduce ( (void*)&numdofs,(void*)&allnumdofs,1,ISSM_MPI_INT,ISSM_MPI_SUM,IssmComm::GetComm());
	return allnumdofs;
}
/*}}}*/
int   Nodes::NumberOfDofsLocal(int analysis_type,int setenum){/*{{{*/

	int   numdofs=0;

	/*Now go through all nodes, and get how many dofs they own, unless they are clone nodes: */
	for(int i=0;i<this->Size();i++){

		Node* node=xDynamicCast<Node*>(this->GetObjectByOffset(i));

		/*Check that this node corresponds to our analysis currently being carried out: */
		if (node->InAnalysis(analysis_type)){

			/*Ok, this object is a node, ask it to plug values into partition: */
			if (!node->IsClone()){
				numdofs+=node->GetNumberOfDofs(NoneApproximationEnum,setenum);
			}
		}
	}

	return numdofs;
}
/*}}}*/
int   Nodes::NumberOfNodes(void){/*{{{*/

	/*Careful! only use once all clones have been setup for all nodes!: */

	int numnodes=0;
	int allnumnodes=0;

	/*Now go through all nodes, and get how many dofs they own, unless they are clone nodes: */
	for(int i=0;i<this->Size();i++){
		Node* node=xDynamicCast<Node*>(this->GetObjectByOffset(i));

		/*Ok, this object is a node, ask it to plug values into partition: */
		if (!node->IsClone()) numnodes++;
	}

	/*Gather from all cpus: */
	ISSM_MPI_Allreduce ( (void*)&numnodes,(void*)&allnumnodes,1,ISSM_MPI_INT,ISSM_MPI_SUM,IssmComm::GetComm());

	return allnumnodes;
}
/*}}}*/
int   Nodes::NumberOfNodes(int analysis_type){/*{{{*/

	int i;

	int max_sid=-1;
	int sid;
	int node_max_sid;

	for(i=0;i<this->Size();i++){

		Node* node=xDynamicCast<Node*>(this->GetObjectByOffset(i));

		/*Check that this node corresponds to our analysis currently being carried out: */
		if (node->InAnalysis(analysis_type)){

			sid=node->Sid();
			if (sid>max_sid)max_sid=sid;
		}
	}

	ISSM_MPI_Reduce (&max_sid,&node_max_sid,1,ISSM_MPI_INT,ISSM_MPI_MAX,0,IssmComm::GetComm() );
	ISSM_MPI_Bcast(&node_max_sid,1,ISSM_MPI_INT,0,IssmComm::GetComm());
	max_sid=node_max_sid;

	/*sid starts at 0*/
	max_sid++;

	/*return*/
	return max_sid;
}
/*}}}*/
void  Nodes::Ranks(int* ranks,int analysis_type){/*{{{*/

	int my_rank;
	int sid;

	/*recover my_rank:*/
	my_rank=IssmComm::GetRank();

	/*Go through nodes, and for each object, report it cpu: */
	for(int i=0;i<this->Size();i++){

		Node* node=xDynamicCast<Node*>(this->GetObjectByOffset(i));

		/*Check that this node corresponds to our analysis currently being carried out: */
		if (node->InAnalysis(analysis_type)){
			/*Plug rank into ranks, according to sid: */
			sid=node->Sid();
			ranks[sid]=my_rank; 
		}
	}
}
/*}}}*/
bool Nodes::RequiresDofReindexing(int analysis_type){/*{{{*/

	int flag = 0;
	int allflag;

	/*Now go through all nodes, and get how many dofs they own, unless they are clone nodes: */
	for(int i=0;i<this->Size();i++){

		Node* node=xDynamicCast<Node*>(this->GetObjectByOffset(i));

		/*Check that this node corresponds to our analysis currently being carried out: */
		if(node->InAnalysis(analysis_type)){
			if(node->RequiresDofReindexing()){
				flag = 1;
				break;
			}
		}
	}

	/*Grab max of all cpus: */
	ISSM_MPI_Allreduce((void*)&flag,(void*)&allflag,1,ISSM_MPI_INT,ISSM_MPI_MAX,IssmComm::GetComm());

	if(allflag){
		return true;
	}
	else{
		return false;
	}
}
/*}}}*/
