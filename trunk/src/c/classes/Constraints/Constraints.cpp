/*
 * \file Constraints.cpp
 * \brief: Implementation of Constraints class, derived from DataSet class.
 */

/*Headers: {{{*/
#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#include "./Constraints.h"
#include "./Constraint.h"
#include "../../shared/shared.h"
#include "../../toolkits/toolkits.h"

using namespace std;
/*}}}*/

/*Numerics: */
void Constraints::ActivatePenaltyMethod(int in_analysis){/*{{{*/

	for(int i=0;i<this->Size();i++){
		Constraint* constraint=(Constraint*)this->GetObjectByOffset(i);
		if(constraint->InAnalysis(in_analysis)){
			constraint->ActivatePenaltyMethod();
		}
	}

}
/*}}}*/
int  Constraints::NumberOfConstraints(void){/*{{{*/

	int localconstraints;
	int numberofconstraints;

	/*Get number of local constraints*/
	localconstraints=this->Size();

	/*figure out total number of constraints combining all the cpus (no clones here)*/
	ISSM_MPI_Reduce(&localconstraints,&numberofconstraints,1,ISSM_MPI_INT,ISSM_MPI_SUM,0,IssmComm::GetComm() );
	ISSM_MPI_Bcast(&numberofconstraints,1,ISSM_MPI_INT,0,IssmComm::GetComm());

	return numberofconstraints;
}
/*}}}*/
