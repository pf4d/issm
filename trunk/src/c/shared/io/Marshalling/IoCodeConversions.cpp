#include <math.h>
#include <float.h>
#include <cstring>

#include "../../Enum/EnumDefinitions.h"
#include "../../MemOps/MemOps.h"
#include "../../Exceptions/exceptions.h"

void FieldAndEnumFromCode(int* out_enum,char** pfield,const char* string_in){/*{{{*/

	/*output*/
	char* fieldname = NULL;
	int   param_enum = -1;

	if(strcmp(string_in,"Thickness")==0){
		const char* field = "md.geometry.thickness";
		param_enum        = ThicknessEnum;
		fieldname=xNew<char>((strlen(field)+1)); xMemCpy<char>(fieldname,field,(strlen(field)+1));
	}
	else if(strcmp(string_in,"MaterialsRheologyB")==0){
		const char* field = "md.materials.rheology_B";
		param_enum        = MaterialsRheologyBEnum;
		fieldname=xNew<char>((strlen(field)+1)); xMemCpy<char>(fieldname,field,(strlen(field)+1));
	}
	else if(strcmp(string_in,"SmbMassBalance")==0){
		const char* field = "md.smb.mass_balance";
		param_enum        = SmbMassBalanceEnum;
		fieldname=xNew<char>((strlen(field)+1)); xMemCpy<char>(fieldname,field,(strlen(field)+1));
	}
	else if(strcmp(string_in,"SmbAccumulation")==0){
		const char* field = "md.smb.accumulation";
		param_enum        = SmbAccumulationEnum;
		fieldname=xNew<char>((strlen(field)+1)); xMemCpy<char>(fieldname,field,(strlen(field)+1));
	}
	else if(strcmp(string_in,"SmbMelt")==0){
		const char* field = "md.smb.melt";
		param_enum        = SmbMeltEnum;
		fieldname=xNew<char>((strlen(field)+1)); xMemCpy<char>(fieldname,field,(strlen(field)+1));
	}
	else if(strcmp(string_in,"SmbRefreeze")==0){
		const char* field = "md.smb.refreeze";
		param_enum        = SmbRefreezeEnum;
		fieldname=xNew<char>((strlen(field)+1)); xMemCpy<char>(fieldname,field,(strlen(field)+1));
	}
	else if(strcmp(string_in,"SmbRunoff")==0){
		const char* field = "md.smb.runoff";
		param_enum        = SmbRunoffEnum;
		fieldname=xNew<char>((strlen(field)+1)); xMemCpy<char>(fieldname,field,(strlen(field)+1));
	}
	else if(strcmp(string_in,"SmbEvaporation")==0){
		const char* field = "md.smb.evaporation";
		param_enum        = SmbEvaporationEnum;
		fieldname=xNew<char>((strlen(field)+1)); xMemCpy<char>(fieldname,field,(strlen(field)+1));
	}
	else if(strcmp(string_in,"BasalforcingsFloatingiceMeltingRate")==0){
		const char* field = "md.basalforcings.floatingice_melting_rate";
		param_enum        = BasalforcingsFloatingiceMeltingRateEnum;
		fieldname=xNew<char>((strlen(field)+1)); xMemCpy<char>(fieldname,field,(strlen(field)+1));
	}
	else if(strcmp(string_in,"FrictionCoefficient")==0){
		const char* field = "md.friction.coefficient";
		param_enum        = FrictionCoefficientEnum;
		fieldname=xNew<char>((strlen(field)+1)); xMemCpy<char>(fieldname,field,(strlen(field)+1));
	}
	else{
		_error_("Field \""<<string_in<<"\" not supported yet");
	}

	/*Assign output pointers*/
	*out_enum = param_enum;
	*pfield   = fieldname;
	return;
}/*}}}*/
int IoCodeToEnumSMB(int enum_in){/*{{{*/
	switch(enum_in){
		case 1: return SMBforcingEnum;
		case 2: return SMBcomponentsEnum;
		case 3: return SMBmeltcomponentsEnum;
		case 4: return SMBpddEnum;
		case 5: return SMBd18opddEnum;
		case 6: return SMBgradientsEnum;
		case 7: return SMBhenningEnum;
		case 8: return SMBgembEnum;
		case 9: return SMBgradientselaEnum;
		default: _error_("Marshalled SMB code \""<<enum_in<<"\" not supported yet");
	}
}/*}}}*/
int IoCodeToEnumBasal(int enum_in){/*{{{*/
	switch(enum_in){
		case 1: return FloatingMeltRateEnum;
		case 2: return LinearFloatingMeltRateEnum;
		case 3: return MismipFloatingMeltRateEnum;
		case 4: return MantlePlumeGeothermalFluxEnum;
		default: _error_("Marshalled Basal Forcings code \""<<enum_in<<"\" not supported yet");
	}
}/*}}}*/
int IoCodeToEnumCalving(int enum_in){/*{{{*/
	switch(enum_in){
		case 1: return DefaultCalvingEnum;
		case 2: return CalvingDevEnum;
		case 3: return CalvingLevermannEnum;
		case 4: return CalvingMinthicknessEnum;
		default: _error_("Marshalled Calving law code \""<<enum_in<<"\" not supported yet");
	}
}/*}}}*/
int IoCodeToEnumHydrology(int enum_in){/*{{{*/
	switch(enum_in){
		case 1: return HydrologydcEnum;
		case 2: return HydrologyshreveEnum;
		case 3: return HydrologysommersEnum;
		default: _error_("Marshalled hydrology code \""<<enum_in<<"\" not supported yet"); 
	}
}/*}}}*/
int IoCodeToEnumMaterials(int enum_in){/*{{{*/
	switch(enum_in){
		case 1: return MatdamageiceEnum;
		case 2: return MatestarEnum; 
		case 3: return MaticeEnum;
		case 4: return MatenhancediceEnum;
		default: _error_("Marshalled materials code \""<<enum_in<<"\" not supported yet"); 
	}
}/*}}}*/

int IoCodeToEnumVertexEquation(int enum_in){/*{{{*/
	switch(enum_in){
		case 0: return NoneApproximationEnum;
		case 1: return SIAApproximationEnum;
		case 2: return SSAApproximationEnum;
		case 3: return L1L2ApproximationEnum;
		case 4: return HOApproximationEnum;
		case 5: return FSApproximationEnum;
		case 6: return SSAHOApproximationEnum;
		case 7: return HOFSApproximationEnum;
		case 8: return SSAFSApproximationEnum;
		default: _error_("Marshalled vertex equation code \""<<enum_in<<"\" not supported yet.");
	}
}/*}}}*/
int IoCodeToEnumElementEquation(int enum_in){/*{{{*/
	switch(enum_in){
		case 0: return NoneApproximationEnum;
		case 1: return SIAApproximationEnum;
		case 2: return SSAApproximationEnum;
		case 3: return L1L2ApproximationEnum;
		case 4: return HOApproximationEnum;
		case 5: return FSApproximationEnum;
		case 6: return SSAHOApproximationEnum;
		case 7: return SSAFSApproximationEnum;
		case 8: return HOFSApproximationEnum;
		default: _error_("Marshalled element equation code \""<<enum_in<<"\" not supported yet.");
	}

}/*}}}*/

int IoRiftfillToEnum(int enum_in){/*{{{*/
	switch(enum_in){
		case 0: return AirEnum;
		case 1: return IceEnum;
		case 2: return MelangeEnum;
		case 3: return WaterEnum;
		default: _error_("Marshalled Riftfill enum \""<<enum_in<<"\" not supported yet.");
	}
}/*}}}*/
