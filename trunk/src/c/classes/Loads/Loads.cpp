/*
 * \file Loads.cpp
 * \brief: Implementation of Loads class, derived from DataSet class.
 */

/*Headers: {{{*/
#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#include <vector>
#include <functional>
#include <algorithm>

#include "../../shared/io/Comm/IssmComm.h"
#include "./Loads.h"
#include "./Load.h"

using namespace std;
/*}}}*/

/*Object constructors and destructor*/
Loads::Loads(){/*{{{*/
	enum_type=LoadsEnum;
	return;
}
/*}}}*/
Loads::~Loads(){/*{{{*/
	return;
}
/*}}}*/

/*Numerics:*/
void Loads::Configure(Elements* elements,Loads* loads, Nodes* nodes, Vertices* vertices, Materials* materials,Parameters* parameters){/*{{{*/

	vector<Object*>::iterator object;
	Load* load=NULL;

	for ( object=objects.begin() ; object < objects.end(); object++ ){

		load=xDynamicCast<Load*>(*object);
		load->Configure(elements,loads,nodes,vertices,materials,parameters);

	}

}
/*}}}*/
bool Loads::IsPenalty(int analysis_type){/*{{{*/

	int ispenalty=0;
	int allispenalty=0;

	/*Now go through all loads, and get how many nodes they own, unless they are clone nodes: */
	for(int i=0;i<this->Size();i++){

		Load* load=xDynamicCast<Load*>(this->GetObjectByOffset(i));
		if (load->InAnalysis(analysis_type)){
			if(load->IsPenalty()) ispenalty++;
		}
	}

	/*Grab sum of all cpus: */
	ISSM_MPI_Allreduce((void*)&ispenalty,(void*)&allispenalty,1,ISSM_MPI_INT,ISSM_MPI_SUM,IssmComm::GetComm());
	ispenalty=allispenalty;

	if(ispenalty)
	 return true;
	else
	 return false;
}
/*}}}*/
int  Loads::MaxNumNodes(int analysis_type){/*{{{*/

	int max=0;
	int allmax;
	int numnodes=0;

	/*Now go through all loads, and get how many nodes they own, unless they are clone nodes: */
	for(int i=0;i<this->Size();i++){

		Load* load=xDynamicCast<Load*>(this->GetObjectByOffset(i));
		if (load->InAnalysis(analysis_type)){
			numnodes=load->GetNumberOfNodes();
			if(numnodes>max)max=numnodes;
		}
	}

	/*Grab max of all cpus: */
	ISSM_MPI_Allreduce((void*)&max,(void*)&allmax,1,ISSM_MPI_INT,ISSM_MPI_MAX,IssmComm::GetComm());
	max=allmax;

	return max;
}
/*}}}*/
int  Loads::NumberOfLoads(void){/*{{{*/

	int localloads;
	int numberofloads;

	/*Get number of local loads*/
	localloads=this->Size();

	/*figure out total number of loads combining all the cpus (no clones here)*/
	ISSM_MPI_Reduce(&localloads,&numberofloads,1,ISSM_MPI_INT,ISSM_MPI_SUM,0,IssmComm::GetComm() );
	ISSM_MPI_Bcast(&numberofloads,1,ISSM_MPI_INT,0,IssmComm::GetComm());

	return numberofloads;
}
/*}}}*/
int  Loads::NumberOfLoads(int analysis_type){/*{{{*/

	int localloads = 0;
	int numberofloads;

	/*Get number of local loads*/
	for(int i=0;i<this->Size();i++){

		Load* load=xDynamicCast<Load*>(this->GetObjectByOffset(i));

		/*Check that this load corresponds to our analysis currently being carried out: */
		if (load->InAnalysis(analysis_type)) localloads++;
	}

	/*figure out total number of loads combining all the cpus (no clones here)*/
	ISSM_MPI_Reduce(&localloads,&numberofloads,1,ISSM_MPI_INT,ISSM_MPI_SUM,0,IssmComm::GetComm() );
	ISSM_MPI_Bcast(&numberofloads,1,ISSM_MPI_INT,0,IssmComm::GetComm());

	return numberofloads;
}
/*}}}*/
void Loads::ResetHooks(){/*{{{*/

	vector<Object*>::iterator object;
	Load* load=NULL;

	for ( object=objects.begin() ; object < objects.end(); object++ ){

		load=xDynamicCast<Load*>((*object));
		load->ResetHooks();

	}

}
/*}}}*/
void Loads::SetCurrentConfiguration(Elements* elements,Loads* loads, Nodes* nodes, Vertices* vertices, Materials* materials,Parameters* parameters){/*{{{*/

	vector<Object*>::iterator object;
	Load* load=NULL;

	for ( object=objects.begin() ; object < objects.end(); object++ ){

		load=xDynamicCast<Load*>(*object);
		load->SetCurrentConfiguration(elements,loads,nodes,vertices,materials,parameters);

	}

}
/*}}}*/
int  Loads::Size(void){/*{{{*/

	return this->DataSet::Size();
}
/*}}}*/
int  Loads::Size(int analysis_type){/*{{{*/

	int localloads = 0;

	/*Get number of local loads*/
	for(int i=0;i<this->Size();i++){

		Load* load=xDynamicCast<Load*>(this->GetObjectByOffset(i));

		/*Check that this load corresponds to our analysis currently being carried out: */
		if (load->InAnalysis(analysis_type)) localloads++;
	}

	return localloads;
}
/*}}}*/
