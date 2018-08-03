/*
 * \file Vertices.cpp
 * \brief: Implementation of Vertices class, derived from DataSet class.
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
#include <iostream>

#include "./Vertices.h"
#include "../shared/shared.h"
#include "./Vertex.h"

using namespace std;
/*}}}*/

/*Object constructors and destructor*/
Vertices::Vertices(){/*{{{*/
	enum_type=VerticesEnum;
	return;
}
/*}}}*/
Vertices::~Vertices(){/*{{{*/
	return;
}
/*}}}*/

/*Numerics management*/
void  Vertices::DistributePids(int numberofobjects){/*{{{*/

	int num_procs;
	int my_rank;

	int  i;
	int  pidcount    = 0;
	int *allpidcount = NULL;
	int *truepids    = NULL;
	int *alltruepids = NULL;

	/*recover my_rank:*/
	my_rank=IssmComm::GetRank();
	num_procs=IssmComm::GetSize();

	/*Go through objects, and distribute pids locally, from 0 to numberofpidsperobject*/
	for (i=0;i<this->Size();i++){
		Vertex* vertex=xDynamicCast<Vertex*>(this->GetObjectByOffset(i));
		vertex->DistributePids(&pidcount);
	}

	/* Now every object has distributed pids, but locally, and with a pid count starting from 
	 * 0. This means the pids between all the cpus are not unique. We now offset the pids of each
	 * cpus by the total last pids of the previus cpu, starting from 0.
	 * First: get number of pids for each cpu*/
	allpidcount=xNew<int>(num_procs);
	ISSM_MPI_Gather(&pidcount,1,ISSM_MPI_INT,allpidcount,1,ISSM_MPI_INT,0,IssmComm::GetComm());
	ISSM_MPI_Bcast(allpidcount,num_procs,ISSM_MPI_INT,0,IssmComm::GetComm());

	/* Every cpu should start its own pid count at the end of the pidcount from cpu-1*/
	pidcount=0;
	if(my_rank!=0){
		for(i=0;i<my_rank;i++){
			pidcount+=allpidcount[i];
		}
	}
	for (i=0;i<this->Size();i++){
		Vertex* vertex=xDynamicCast<Vertex*>(this->GetObjectByOffset(i));
		vertex->OffsetPids(pidcount);
	}

	/* Finally, remember that cpus may have skipped some objects, because they were clones. For every 
	 * object that is not a clone, tell them to show their pids, so that later on, they can get picked 
	 * up by their clones: */
	truepids   =xNewZeroInit<int>(numberofobjects);
	alltruepids=xNewZeroInit<int>(numberofobjects);
	for (i=0;i<this->Size();i++){
		Vertex* vertex=xDynamicCast<Vertex*>(this->GetObjectByOffset(i));
		vertex->ShowTruePids(truepids);
	}
	ISSM_MPI_Allreduce((void*)truepids,(void*)alltruepids,numberofobjects,ISSM_MPI_INT,ISSM_MPI_MAX,IssmComm::GetComm());

	/* Now every cpu knows the true pids of everyone else that is not a clone*/
	for(i=0;i<this->Size();i++){
		Vertex* vertex=xDynamicCast<Vertex*>(this->GetObjectByOffset(i));
		vertex->UpdateClonePids(alltruepids);
	}

	/* Free ressources: */
	xDelete<int>(allpidcount);
	xDelete<int>(truepids);
	xDelete<int>(alltruepids);
}
/*}}}*/
void  Vertices::FlagClones(int numberofobjects){/*{{{*/

	int i;
	int num_procs;

	int* ranks=NULL;
	int* minranks=NULL;

	/*recover num_procs:*/
	num_procs=IssmComm::GetSize();

	/*Allocate ranks: */
	ranks=xNew<int>(numberofobjects);
	minranks=xNew<int>(numberofobjects);

	for(i=0;i<numberofobjects;i++)ranks[i]=num_procs; //no cpu can have rank num_procs. This is the maximum limit.

	/*Now go through all our objects and ask them to report to who they belong (which rank): */
	Ranks(ranks);

	/*We need to take the minimum rank for each vertex, and every cpu needs to get that result. That way, 
	 * when we start building the dof list for all vertexs, a cpu can check whether its vertex already has been 
	 * dealt with by another cpu. We take the minimum because we are going to manage dof assignment in increasing 
	 * order of cpu rank. This is also why we initialized this array to num_procs.*/
	ISSM_MPI_Allreduce ( (void*)ranks,(void*)minranks,numberofobjects,ISSM_MPI_INT,ISSM_MPI_MIN,IssmComm::GetComm());

	/*Now go through all objects, and use minranks to flag which objects are cloned: */
	for(i=0;i<this->Size();i++){
		/*For this object, decide whether it is a clone: */
		Vertex* vertex=xDynamicCast<Vertex*>(this->GetObjectByOffset(i));
		vertex->SetClone(minranks);
	}

	/*Free ressources: */
	xDelete<int>(ranks); 
	xDelete<int>(minranks);

}
/*}}}*/
int Vertices::NumberOfVertices(void){/*{{{*/

	int i,sid;
	int max_sid=0;
	int vertex_max_sid;

	for(i=0;i<this->Size();i++){
		Vertex* vertex=xDynamicCast<Vertex*>(this->GetObjectByOffset(i));
		sid=vertex->Sid();
		if (sid>max_sid)max_sid=sid;
	}

	ISSM_MPI_Reduce (&max_sid,&vertex_max_sid,1,ISSM_MPI_INT,ISSM_MPI_MAX,0,IssmComm::GetComm() );
	ISSM_MPI_Bcast(&vertex_max_sid,1,ISSM_MPI_INT,0,IssmComm::GetComm());
	max_sid=vertex_max_sid;

	/*sid starts at 0*/
	max_sid++;

	/*return:*/
	return max_sid;
}
/*}}}*/
void   Vertices::Ranks(int* ranks){/*{{{*/

	int my_rank;
	int        sid;

	/*recover my_rank:*/
	my_rank=IssmComm::GetRank();

	/*Go through a dataset, and for each object, report it cpu: */
	for(int i=0;i<this->Size();i++){
		/*Plug rank into ranks, according to id: */
		Vertex* vertex=xDynamicCast<Vertex*>(this->GetObjectByOffset(i));
		sid=vertex->Sid();
		ranks[sid]=my_rank; 
	}
}
/*}}}*/
IssmDouble* Vertices::ToXYZ(void){/*{{{*/

	/*intermediary: */
	int i;
	int my_rank;
	int num_vertices;

	/*output: */
	Matrix<IssmDouble>* xyz = NULL;
	IssmDouble* xyz_serial=NULL;

	/*recover my_rank:*/
	my_rank=IssmComm::GetRank();

	/*First, figure out number of vertices: */
	num_vertices=this->NumberOfVertices();

	/*Now, allocate matrix to hold all the vertices x,y and z values: */
	xyz= new Matrix<IssmDouble>(num_vertices,3);

	/*Go through vertices, and for each vertex, object, report it cpu: */
	for(i=0;i<this->Size();i++){

		/*let vertex fill matrix: */
		Vertex* vertex=xDynamicCast<Vertex*>(this->GetObjectByOffset(i));
		vertex->ToXYZ(xyz);
	}

	/*Assemble:*/
	xyz->Assemble();

	/*gather on cpu 0: */
	xyz_serial=xyz->ToSerial();

	/*free ressources: */
	delete xyz;
	if(my_rank!=0)delete xyz_serial;

	/*return matrix: */
	return xyz_serial;
}
/*}}}*/
