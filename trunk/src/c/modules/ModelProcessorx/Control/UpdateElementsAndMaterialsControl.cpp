/*
 * UpdateElementsAndMaterialsControl:
 */

#include "../../../toolkits/toolkits.h"
#include "../../../classes/classes.h"
#include "../../../shared/shared.h"
#include "../../MeshPartitionx/MeshPartitionx.h"
#include "../ModelProcessorx.h"

void	UpdateElementsAndMaterialsControl(Elements* elements,Materials* materials, IoModel* iomodel){

	/*Intermediary*/
	bool       control_analysis;
	int        control,cost_function,domaintype;
	int        num_controls,num_cost_functions;
	Element   *element          = NULL;
	Material  *material         = NULL;
	int       *control_enums    = NULL;
	char     **controls         = NULL;
	char     **cost_functions   = NULL;

	/*Fetch parameters: */
	iomodel->FindConstant(&control_analysis,"md.inversion.iscontrol");
	if(control_analysis) iomodel->FindConstant(&num_controls,"md.inversion.num_control_parameters");

	/*Now, return if no control*/
	if(!control_analysis) return;

	/*Process controls and convert from string to enums*/
	iomodel->FindConstant(&controls,&num_controls,"md.inversion.control_parameters");
	if(num_controls<1) _error_("no controls found");
	control_enums=xNew<int>(num_controls);
	for(int i=0;i<num_controls;i++){
		control_enums[i]=StringToEnumx(controls[i]);
	}

	/*Process cost functions and convert from string to enums*/
	iomodel->FindConstant(&num_cost_functions,"md.inversion.num_cost_functions");
	iomodel->FindConstant(&cost_functions,&num_cost_functions,"md.inversion.cost_functions");
	if(num_cost_functions<1) _error_("No cost functions found");
	int* cost_function_enums=xNew<int>(num_cost_functions);
	for(int i=0;i<num_cost_functions;++i){
		cost_function_enums[i]=StringToEnumx(cost_functions[i]);
	}

	iomodel->FetchData(3,"md.inversion.cost_functions_coefficients","md.inversion.min_parameters","md.inversion.max_parameters");
	
	/*Fetch Observations */
	iomodel->FindConstant(&domaintype,"md.mesh.domain_type");
	for(int i=0;i<num_cost_functions;i++){
		cost_function=cost_function_enums[i];
		if(     cost_function==ThicknessAbsMisfitEnum) iomodel->FetchDataToInput(elements,"md.inversion.thickness_obs",InversionThicknessObsEnum);
		else if(cost_function==SurfaceAbsMisfitEnum)   iomodel->FetchDataToInput(elements,"md.inversion.surface_obs",InversionSurfaceObsEnum);
		else if(cost_function==SurfaceAbsVelMisfitEnum
			  || cost_function==SurfaceRelVelMisfitEnum
			  || cost_function==SurfaceLogVelMisfitEnum
			  || cost_function==SurfaceLogVxVyMisfitEnum
			  || cost_function==SurfaceAverageVelMisfitEnum){
			iomodel->FetchDataToInput(elements,"md.inversion.vx_obs",InversionVxObsEnum);
			if(domaintype!=Domain2DverticalEnum) iomodel->FetchDataToInput(elements,"md.inversion.vy_obs",InversionVyObsEnum); 
		}
	}

	for(int i=0;i<num_controls;i++){
		control = control_enums[i];
		switch(control){
			/*List of supported controls*/
			case BalancethicknessThickeningRateEnum:      iomodel->FetchData(1,"md.balancethickness.thickening_rate"); break;
			case VxEnum:                                  iomodel->FetchData(1,"md.initialization.vx"); break;
			case VyEnum:                                  iomodel->FetchData(1,"md.initialization.vy"); break;
			case ThicknessEnum:                           iomodel->FetchData(1,"md.geometry.thickness"); break;
			case FrictionCoefficientEnum:                 iomodel->FetchData(1,"md.friction.coefficient"); break;
			case FrictionAsEnum:                          iomodel->FetchData(1,"md.friction.As"); break;
			case BalancethicknessApparentMassbalanceEnum: iomodel->FetchData(1,"md.balancethickness.apparent_massbalance"); break;
			case BalancethicknessOmegaEnum:               iomodel->FetchData(1,"md.balancethickness.omega"); break;
			case MaterialsRheologyBEnum:                  iomodel->FetchData(1,"md.materials.rheology_B"); break;
			/*Special cases*/
			case MaterialsRheologyBbarEnum: iomodel->FetchData(1,"md.materials.rheology_B"); break;
			case DamageDbarEnum:            iomodel->FetchData(1,"md.damage.D");            break;
			default:
				_error_("Control " << EnumToStringx(control) << " not implemented yet");
		}
	}

	/*Update elements: */
	int counter=0;
	for(int i=0;i<iomodel->numberofelements;i++){
		if(iomodel->my_elements[i]){
			element=(Element*)elements->GetObjectByOffset(counter);
			element->InputUpdateFromIoModel(i,iomodel); //we need i to index into elements.
			counter++;
		}
	}

	/*Free data: */
	for(int i=0;i<num_controls;i++){
		switch(control_enums[i]){
			/*List of supported controls*/
			case BalancethicknessThickeningRateEnum:      iomodel->DeleteData(1,"md.balancethickness.thickening_rate"); break;
			case VxEnum:                                  iomodel->DeleteData(1,"md.initialization.vx"); break;
			case VyEnum:                                  iomodel->DeleteData(1,"md.initialization.vy"); break;
			case ThicknessEnum:                           iomodel->DeleteData(1,"md.geometry.thickness"); break;
			case FrictionCoefficientEnum:                 iomodel->DeleteData(1,"md.friction.coefficient"); break;
			case FrictionAsEnum:                          iomodel->DeleteData(1,"md.friction.As"); break;
			case BalancethicknessApparentMassbalanceEnum: iomodel->DeleteData(1,"md.balancethickness.apparent_massbalance"); break;
			case BalancethicknessOmegaEnum:               iomodel->DeleteData(1,"md.balancethickness.omega"); break;
			case MaterialsRheologyBEnum:                  iomodel->DeleteData(1,"md.materials.rheology_B"); break;
			/*Special cases*/
			case MaterialsRheologyBbarEnum: iomodel->DeleteData(1,"md.materials.rheology_B"); break;
			case DamageDbarEnum:            iomodel->DeleteData(1,"md.damage.D");            break;
			default:
				_error_("Control " << EnumToStringx(control_enums[i]) << " not implemented yet");
		}
	}

	iomodel->DeleteData(3,"md.inversion.cost_functions_coefficients","md.inversion.min_parameters","md.inversion.max_parameters");
	xDelete<int>(control_enums);
	xDelete<int>(cost_function_enums);
	for(int i=0;i<num_cost_functions;i++) xDelete<char>(cost_functions[i]);
	xDelete<char*>(cost_functions);
	for(int i=0;i<num_controls;i++) xDelete<char>(controls[i]);
	xDelete<char*>(controls);
}
