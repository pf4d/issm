/*!\file ConstraintsStatex
 * \brief: set up penalty constraints on loads
 */

#include "./ConstraintsStatex.h"
#include "./ConstraintsStateLocal.h"
#include "../../shared/shared.h"
#include "../../toolkits/toolkits.h"

void ConstraintsStatex(int* pconverged, int* pnum_unstable_constraints,FemModel* femmodel){

	/*output: */
	int converged                     = 1;
	int num_unstable_constraints      = 0;
	int min_mechanical_constraints    = 0;
	int  unstable                     = 0;
	int  sum_num_unstable_constraints = 0;
	int analysis_type;

	/*Display message*/
	if(VerboseModule()) _printf0_("   Constraining penalties\n");

	/*recover parameters: */
	femmodel->parameters->FindParam(&min_mechanical_constraints,StressbalanceRiftPenaltyThresholdEnum);
	femmodel->parameters->FindParam(&analysis_type,AnalysisTypeEnum);

	/*Rift penalties first*/
	if(RiftIsPresent(femmodel->loads,analysis_type)){
		RiftConstraintsState(&converged,&num_unstable_constraints,femmodel->loads,min_mechanical_constraints,analysis_type);
	}

	/*Deal with pengrid*/
	for(int i=0;i<femmodel->loads->Size();i++){
		Load* load=(Load*)femmodel->loads->GetObjectByOffset(i);
		if(load->InAnalysis(analysis_type)){
			if(load->ObjectEnum()==PengridEnum){
				Pengrid* pengrid=(Pengrid*)load;
				pengrid->ConstraintActivate(&unstable);
				num_unstable_constraints += unstable;
			}
		}
	}
	ISSM_MPI_Reduce(&num_unstable_constraints,&sum_num_unstable_constraints,1,ISSM_MPI_INT,ISSM_MPI_SUM,0,IssmComm::GetComm() );
	ISSM_MPI_Bcast(&sum_num_unstable_constraints,1,ISSM_MPI_INT,0,IssmComm::GetComm());                
	num_unstable_constraints=sum_num_unstable_constraints;

	/*Have we converged? : */
	if(num_unstable_constraints) converged=0;

	/*Assign output pointers: */
	*pconverged                = converged;
	*pnum_unstable_constraints = num_unstable_constraints;
}
