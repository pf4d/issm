/*
 * cores.h: 
 */

#ifndef _CORES_H_
#define _CORES_H_

/*forward declarations: */
class FemModel;
class Parameters;
template <class doubletype> class Matrix;
template <class doubletype> class Vector;

#include "../shared/io/Comm/IssmComm.h"
#include "../shared/Numerics/types.h"

/*cores: */
void adjointstressbalance_core(FemModel* femmodel);
void adjointbalancethickness_core(FemModel* femmodel);
void adjointbalancethickness2_core(FemModel* femmodel);
void stressbalance_core(FemModel* femmodel);
void hydrology_core(FemModel* femmodel);
void thermal_core(FemModel* femmodel);
void surfaceslope_core(FemModel* femmodel);
void levelsetfunctionslope_core(FemModel* femmodel);
void movingfront_core(FemModel* femmodel);
void bedslope_core(FemModel* femmodel);
void meshdeformation_core(FemModel* femmodel);
void control_core(FemModel* femmodel);
void controltao_core(FemModel* femmodel);
void controlm1qn3_core(FemModel* femmodel);
void controlad_core(FemModel* femmodel);
void controlvalidation_core(FemModel* femmodel);
void masstransport_core(FemModel* femmodel);
void depthaverage_core(FemModel* femmodel);
void extrudefrombase_core(FemModel* femmodel);
void extrudefromtop_core(FemModel* femmodel);
void balancethickness_core(FemModel* femmodel);
void balancethickness2_core(FemModel* femmodel);
void balancevelocity_core(FemModel* femmodel);
void slopecompute_core(FemModel* femmodel);
void steadystate_core(FemModel* femmodel);
void transient_core(FemModel* femmodel);
void dakota_core(FemModel* femmodel);
void ad_core(FemModel* femmodel);
void adgradient_core(FemModel* femmodel);
void dummy_core(FemModel* femmodel);
void gia_core(FemModel* femmodel);
void esa_core(FemModel* femmodel);
void smb_core(FemModel* femmodel);
void damage_core(FemModel* femmodel);
void sealevelrise_core(FemModel* femmodel);
Vector<IssmDouble>* sealevelrise_core_eustatic(FemModel* femmodel);
Vector<IssmDouble>* sealevelrise_core_noneustatic(FemModel* femmodel,Vector<IssmDouble>* Sg_eustatic);
IssmDouble objectivefunction(IssmDouble search_scalar,FemModel* femmodel);

//optimization
int GradJSearch(IssmDouble* search_vector,FemModel* femmodel,int step);

//diverse
void ProcessArguments(int* solution,char** pbinname,char** poutbinname,char** ptoolkitsname,char** plockname,char** prestartname, char** prootpath,int argc,char **argv);
void WriteLockFile(char* filename);
void ResetBoundaryConditions(FemModel* femmodel, int analysis_type);
void PrintBanner(void);
void TransferForcing(FemModel* femmodel,int forcingenum);
void TransferSealevel(FemModel* femmodel,int forcingenum);

//solution configuration
void CorePointerFromSolutionEnum(void (**psolutioncore)(FemModel*),Parameters* parameters,int solutiontype);
void WrapperCorePointerFromSolutionEnum(void (**psolutioncore)(FemModel*),Parameters* parameters,int solutiontype,bool nodakotacore=false);
void AdjointCorePointerFromSolutionEnum(void (**padjointcore)(FemModel*),int solutiontype);

#endif
