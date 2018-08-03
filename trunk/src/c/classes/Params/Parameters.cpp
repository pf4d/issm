/*
 * \file Parameters.cpp
 * \brief: Implementation of the Parameters class, derived from DataSet class.
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

#include "./Parameters.h"
#include "./Param.h"

#include "./BoolParam.h"
#include "./DoubleMatParam.h"
#include "./DataSetParam.h"
#include "./DoubleParam.h"
#include "./DoubleVecParam.h"
#include "./IntParam.h"
#include "./IntVecParam.h"
#include "./IntMatParam.h"
#include "./FileParam.h"
#include "./MatrixParam.h"
#include "./VectorParam.h"
#include "./StringArrayParam.h"
#include "./StringParam.h"
#include "./DoubleMatArrayParam.h"
#include "./TransientParam.h"

#include "../../shared/shared.h"
#include "../../toolkits/toolkits.h"

using namespace std;
/*}}}*/

/*Object constructors and destructor*/
Parameters::Parameters(){/*{{{*/
	for(int i=0;i<NUMPARAMS;i++) this->params[i] = NULL;
	return;
}
/*}}}*/
Parameters::~Parameters(){/*{{{*/
	for(int i=0;i<NUMPARAMS;i++){
		if(this->params[i]) delete this->params[i];
	}
	return;
}
/*}}}*/

void Parameters::AddObject(Param* newparam){/*{{{*/

	/*Get Enum from Param*/
	_assert_(newparam);
	int param_enum = newparam->InstanceEnum();

	/*Get index in array*/
	#ifdef _ISSM_DEBUG_
	if(param_enum<=ParametersSTARTEnum) _error_("Enum "<<EnumToStringx(param_enum)<<" should appear after ParametersSTARTEnum");
	if(param_enum>=ParametersENDEnum)   _error_("Enum "<<EnumToStringx(param_enum)<<" should appear before ParametersENDEnum");
	#endif
	int index = param_enum - ParametersSTARTEnum -1;

	/*Delete param if it already exists*/
	if(this->params[index]){
		delete this->params[index];
		this->params[index] = NULL;
	}

	/*Add param to array*/
	this->params[index] = newparam;
}
/*}}}*/
Parameters* Parameters::Copy(void){/*{{{*/

	Parameters* output = new Parameters();

	for(int i=0;i<NUMPARAMS;i++){
		if(this->params[i]){
			output->params[i]=this->params[i]->copy();
		}
	}

	return output;
}
/*}}}*/
void Parameters::DeepEcho(void){/*{{{*/	
	for(int i=0;i<NUMPARAMS;i++) {
		if(this->params[i]) this->params[i]->DeepEcho();
	}
	return;
}
/*}}}*/
void Parameters::Echo(void){/*{{{*/	
	for(int i=0;i<NUMPARAMS;i++) {
		if(this->params[i]) this->params[i]->Echo();
	}
	return;
}
/*}}}*/
void Parameters::Marshall(char** pmarshalled_data, int* pmarshalled_data_size, int marshall_direction){/*{{{*/

	int obj_enum=-1;
	int num_params=0;

	MARSHALLING_ENUM(ParametersEnum);

	if(marshall_direction==MARSHALLING_FORWARD || marshall_direction==MARSHALLING_SIZE){

		/*Marshall num_params first*/
		for(int i=0;i<NUMPARAMS;i++){
			if(this->params[i]) num_params++;
		}
		MARSHALLING(num_params);

		/*Marshall Parameters one by one now*/
		for(int i=0;i<NUMPARAMS;i++){
			if(this->params[i]){
				obj_enum = this->params[i]->ObjectEnum();
				MARSHALLING(obj_enum);
				this->params[i]->Marshall(pmarshalled_data,pmarshalled_data_size,marshall_direction);
			}
		}
	}
	else{

		/*Get number of params marshalled*/
		MARSHALLING(num_params);

		/*Recover parameters one by one*/
		for(int i=0;i<num_params;i++){

			/*Recover enum of object first: */
			MARSHALLING(obj_enum); 

			if(obj_enum==DoubleParamEnum){
				DoubleParam* doubleparam=NULL;
				doubleparam=new DoubleParam();
				doubleparam->Marshall(pmarshalled_data,pmarshalled_data_size,marshall_direction);
				this->AddObject(doubleparam);
			}
			else if(obj_enum==IntParamEnum){
				IntParam* intparam=NULL;
				intparam=new IntParam();
				intparam->Marshall(pmarshalled_data,pmarshalled_data_size,marshall_direction);
				this->AddObject(intparam);
			}
			else if(obj_enum==IntMatParamEnum){
				IntMatParam* intmparam=NULL;
				intmparam=new IntMatParam();
				intmparam->Marshall(pmarshalled_data,pmarshalled_data_size,marshall_direction);
				this->AddObject(intmparam);
			}
			else if(obj_enum==IntVecParamEnum){
				IntVecParam* intvparam=NULL;
				intvparam=new IntVecParam();
				intvparam->Marshall(pmarshalled_data,pmarshalled_data_size,marshall_direction);
				this->AddObject(intvparam);
			}
			else if(obj_enum==BoolParamEnum){
				BoolParam* boolparam=NULL;
				boolparam=new BoolParam();
				boolparam->Marshall(pmarshalled_data,pmarshalled_data_size,marshall_direction);
				this->AddObject(boolparam);
			}
			else if(obj_enum==DataSetParamEnum){
				DataSetParam* dsparam=NULL;
				dsparam=new DataSetParam();
				dsparam->Marshall(pmarshalled_data,pmarshalled_data_size,marshall_direction);
				this->AddObject(dsparam);
			}
			else if(obj_enum==DoubleMatArrayParamEnum){
				DoubleMatArrayParam* dmaparam=NULL;
				dmaparam=new DoubleMatArrayParam();
				dmaparam->Marshall(pmarshalled_data,pmarshalled_data_size,marshall_direction);
				this->AddObject(dmaparam);
			}
			else if(obj_enum==DoubleMatParamEnum){
				DoubleMatParam* dmparam=NULL;
				dmparam=new DoubleMatParam();
				dmparam->Marshall(pmarshalled_data,pmarshalled_data_size,marshall_direction);
				this->AddObject(dmparam);
			}
			else if(obj_enum==DoubleVecParamEnum){
				DoubleVecParam* dvparam=NULL;
				dvparam=new DoubleVecParam();
				dvparam->Marshall(pmarshalled_data,pmarshalled_data_size,marshall_direction);
				this->AddObject(dvparam);
			}
			else if(obj_enum==FileParamEnum){
				FileParam* fileparam=NULL;
				fileparam=new FileParam();
				fileparam->Marshall(pmarshalled_data,pmarshalled_data_size,marshall_direction);
				delete fileparam;
				/* No need to add this object, the pointer is not valid             
					The FemModel should reset all FileParams in the restart function */
			}
			else if(obj_enum==StringParamEnum){
				StringParam* sparam=NULL;
				sparam=new StringParam();
				sparam->Marshall(pmarshalled_data,pmarshalled_data_size,marshall_direction);
				this->AddObject(sparam);
			}
			else if(obj_enum==StringArrayParamEnum){
				StringArrayParam* saparam=NULL;
				saparam=new StringArrayParam();
				saparam->Marshall(pmarshalled_data,pmarshalled_data_size,marshall_direction);
				this->AddObject(saparam);
			}
			else if(obj_enum==TransientParamEnum){
				TransientParam* transparam=NULL;
				transparam=new TransientParam();
				transparam->Marshall(pmarshalled_data,pmarshalled_data_size,marshall_direction);
				this->AddObject(transparam);
			}
			else if(obj_enum==GenericParamEnum){
				/*Skip for now (we don't want to Marhsall Comms*/
			}
		}
	}
}
/*}}}*/

/*Object management*/
void Parameters::Delete(int param_enum){/*{{{*/

	_assert_(param_enum>ParametersSTARTEnum);
	_assert_(param_enum<ParametersENDEnum);

	int index = param_enum - ParametersSTARTEnum -1;
	if(this->params[index]){
		delete this->params[index];
		this->params[index] = NULL;
	}

	return;
}
/*}}}*/
bool Parameters::Exist(int param_enum){/*{{{*/

	_assert_(param_enum>ParametersSTARTEnum);
	_assert_(param_enum<ParametersENDEnum);

	int index = param_enum - ParametersSTARTEnum -1;
	if(this->params[index]) return true;

	return false;
}
/*}}}*/
void Parameters::FindParam(bool* pbool,int param_enum){ _assert_(this);/*{{{*/

	_assert_(param_enum>ParametersSTARTEnum);
	_assert_(param_enum<ParametersENDEnum);

	int index = param_enum - ParametersSTARTEnum -1;
	if(!this->params[index]) _error_("Parameter " << EnumToStringx(param_enum) <<" not set");
	this->params[index]->GetParameterValue(pbool);
}
/*}}}*/
void Parameters::FindParam(int* pinteger,int param_enum){ _assert_(this);/*{{{*/

	_assert_(param_enum>ParametersSTARTEnum);
	_assert_(param_enum<ParametersENDEnum);

	int index = param_enum - ParametersSTARTEnum -1;
	if(!this->params[index]) _error_("Parameter " << EnumToStringx(param_enum) <<" not set");
	this->params[index]->GetParameterValue(pinteger);
}
/*}}}*/
void Parameters::FindParam(IssmDouble* pscalar,int param_enum){ _assert_(this);/*{{{*/

	_assert_(param_enum>ParametersSTARTEnum);
	_assert_(param_enum<ParametersENDEnum);

	int index = param_enum - ParametersSTARTEnum -1;
	if(!this->params[index]) _error_("Parameter " << EnumToStringx(param_enum) <<" not set");
	this->params[index]->GetParameterValue(pscalar);
}
/*}}}*/
void Parameters::FindParam(IssmDouble* pscalar, int param_enum,IssmDouble time){ _assert_(this);/*{{{*/

	_assert_(param_enum>ParametersSTARTEnum);
	_assert_(param_enum<ParametersENDEnum);

	int index = param_enum - ParametersSTARTEnum -1;
	if(!this->params[index]) _error_("Parameter " << EnumToStringx(param_enum) <<" not set");
	this->params[index]->GetParameterValue(pscalar,time);
}
/*}}}*/
void Parameters::FindParam(char** pstring,int param_enum){ _assert_(this);/*{{{*/

	_assert_(param_enum>ParametersSTARTEnum);
	_assert_(param_enum<ParametersENDEnum);

	int index = param_enum - ParametersSTARTEnum -1;
	if(!this->params[index]) _error_("Parameter " << EnumToStringx(param_enum) <<" not set");
	this->params[index]->GetParameterValue(pstring);

}
/*}}}*/
void Parameters::FindParam(char*** pstringarray,int* pM,int param_enum){ _assert_(this);/*{{{*/

	_assert_(param_enum>ParametersSTARTEnum);
	_assert_(param_enum<ParametersENDEnum);

	int index = param_enum - ParametersSTARTEnum -1;
	if(!this->params[index]) _error_("Parameter " << EnumToStringx(param_enum) <<" not set");
	this->params[index]->GetParameterValue(pstringarray,pM);
}
/*}}}*/
void Parameters::FindParam(int** pintarray,int* pM, int param_enum){ _assert_(this);/*{{{*/

	_assert_(param_enum>ParametersSTARTEnum);
	_assert_(param_enum<ParametersENDEnum);

	int index = param_enum - ParametersSTARTEnum -1;
	if(!this->params[index]) _error_("Parameter " << EnumToStringx(param_enum) <<" not set");
	this->params[index]->GetParameterValue(pintarray,pM);

}
/*}}}*/
void Parameters::FindParam(int** pintarray,int* pM,int *pN,int param_enum){ _assert_(this);/*{{{*/

	_assert_(param_enum>ParametersSTARTEnum);
	_assert_(param_enum<ParametersENDEnum);

	int index = param_enum - ParametersSTARTEnum -1;
	if(!this->params[index]) _error_("Parameter " << EnumToStringx(param_enum) <<" not set");
	this->params[index]->GetParameterValue(pintarray,pM,pN);

}
/*}}}*/
void Parameters::FindParam(IssmDouble** pIssmDoublearray,int* pM, int param_enum){ _assert_(this);/*{{{*/

	_assert_(param_enum>ParametersSTARTEnum);
	_assert_(param_enum<ParametersENDEnum);

	int index = param_enum - ParametersSTARTEnum -1;
	if(!this->params[index]) _error_("Parameter " << EnumToStringx(param_enum) <<" not set");
	this->params[index]->GetParameterValue(pIssmDoublearray,pM);
}
/*}}}*/
void Parameters::FindParam(IssmDouble** pIssmDoublearray,int* pM, int* pN,int param_enum){ _assert_(this);/*{{{*/

	_assert_(param_enum>ParametersSTARTEnum);
	_assert_(param_enum<ParametersENDEnum);

	int index = param_enum - ParametersSTARTEnum -1;
	if(!this->params[index]) _error_("Parameter " << EnumToStringx(param_enum) <<" not set");
	this->params[index]->GetParameterValue(pIssmDoublearray,pM,pN);
}
/*}}}*/
void Parameters::FindParam(IssmDouble*** parray,int* pM,int** pmdims_array,int** pndims_array,int param_enum){ _assert_(this);/*{{{*/

	_assert_(param_enum>ParametersSTARTEnum);
	_assert_(param_enum<ParametersENDEnum);

	int index = param_enum - ParametersSTARTEnum -1;
	if(!this->params[index]) _error_("Parameter " << EnumToStringx(param_enum) <<" not set");
	this->params[index]->GetParameterValue(parray,pM,pmdims_array,pndims_array);
}
/*}}}*/
void Parameters::FindParam(Vector<IssmDouble>** pvec,int param_enum){ _assert_(this);/*{{{*/

	_assert_(param_enum>ParametersSTARTEnum);
	_assert_(param_enum<ParametersENDEnum);

	int index = param_enum - ParametersSTARTEnum -1;
	if(!this->params[index]) _error_("Parameter " << EnumToStringx(param_enum) <<" not set");
	this->params[index]->GetParameterValue(pvec);
}
/*}}}*/
void Parameters::FindParam(Matrix<IssmDouble>** pmat,int param_enum){ _assert_(this);/*{{{*/

	_assert_(param_enum>ParametersSTARTEnum);
	_assert_(param_enum<ParametersENDEnum);

	int index = param_enum - ParametersSTARTEnum -1;
	if(!this->params[index]) _error_("Parameter " << EnumToStringx(param_enum) <<" not set");
	this->params[index]->GetParameterValue(pmat);
}
/*}}}*/
void Parameters::FindParam(FILE** pfid,int param_enum){ _assert_(this);/*{{{*/

	_assert_(param_enum>ParametersSTARTEnum);
	_assert_(param_enum<ParametersENDEnum);

	int index = param_enum - ParametersSTARTEnum -1;
	if(!this->params[index]) _error_("Parameter " << EnumToStringx(param_enum) <<" not set");
	this->params[index]->GetParameterValue(pfid);
}
/*}}}*/
void Parameters::FindParam(DataSet** pdataset,int param_enum){ /*{{{*/
	_assert_(this);

	_assert_(param_enum>ParametersSTARTEnum);
	_assert_(param_enum<ParametersENDEnum);

	int index = param_enum - ParametersSTARTEnum -1;
	if(!this->params[index]) _error_("Parameter " << EnumToStringx(param_enum) <<" not set");
	this->params[index]->GetParameterValue(pdataset);
}
/*}}}*/

void   Parameters::SetParam(bool boolean,int enum_type){/*{{{*/

	Param* param=NULL;

	/*first, figure out if the param has already been created: */
	param=xDynamicCast<Param*>(this->FindParamObject(enum_type));

	if(param) param->SetValue(boolean); //already exists, just set it.
	else this->AddObject(new BoolParam(enum_type,boolean)); //just add the new parameter.
}
/*}}}*/
void   Parameters::SetParam(int integer,int enum_type){/*{{{*/

	Param* param=NULL;

	/*first, figure out if the param has already been created: */
	param=xDynamicCast<Param*>(this->FindParamObject(enum_type));

	if(param) param->SetValue(integer); //already exists, just set it.
	else this->AddObject(new IntParam(enum_type,integer)); //just add the new parameter.
}
/*}}}*/
void   Parameters::SetParam(IssmDouble scalar,int enum_type){/*{{{*/

	Param* param=NULL;

	/*first, figure out if the param has already been created: */
	param=xDynamicCast<Param*>(this->FindParamObject(enum_type));

	if(param) param->SetValue(scalar); //already exists, just set it.
	else this->AddObject(new DoubleParam(enum_type,scalar)); //just add the new parameter.
}
/*}}}*/
void   Parameters::SetParam(char* string,int enum_type){/*{{{*/

	Param* param=NULL;

	/*first, figure out if the param has already been created: */
	param=xDynamicCast<Param*>(this->FindParamObject(enum_type));

	if(param) param->SetValue(string); //already exists, just set it.
	else this->AddObject(new StringParam(enum_type,string)); //just add the new parameter.
}
/*}}}*/
void   Parameters::SetParam(char** stringarray,int M, int enum_type){/*{{{*/

	Param* param=NULL;

	/*first, figure out if the param has already been created: */
	param=xDynamicCast<Param*>(this->FindParamObject(enum_type));

	if(param) param->SetValue(stringarray,M); //already exists, just set it.
	else this->AddObject(new StringArrayParam(enum_type,stringarray,M)); //just add the new parameter.
}
/*}}}*/
void   Parameters::SetParam(IssmDouble* IssmDoublearray,int M, int enum_type){/*{{{*/

	Param* param=NULL;

	/*first, figure out if the param has already been created: */
	param=xDynamicCast<Param*>(this->FindParamObject(enum_type));

	if(param) param->SetValue(IssmDoublearray,M); //already exists, just set it.
	else this->AddObject(new DoubleVecParam(enum_type,IssmDoublearray,M)); //just add the new parameter.
}
/*}}}*/
void   Parameters::SetParam(IssmDouble* IssmDoublearray,int M, int N, int enum_type){/*{{{*/

	Param* param=NULL;

	/*first, figure out if the param has already been created: */
	param=xDynamicCast<Param*>(this->FindParamObject(enum_type));

	if(param) param->SetValue(IssmDoublearray,M,N); //already exists, just set it.
	else this->AddObject(new DoubleMatParam(enum_type,IssmDoublearray,M,N)); //just add the new parameter.
}
/*}}}*/
void   Parameters::SetParam(int* intarray,int M, int enum_type){/*{{{*/

	Param* param=NULL;

	/*first, figure out if the param has already been created: */
	param=xDynamicCast<Param*>(this->FindParamObject(enum_type));

	if(param) param->SetValue(intarray,M); //already exists, just set it.
	else this->AddObject(new IntVecParam(enum_type,intarray,M)); //just add the new parameter.
}
/*}}}*/
void   Parameters::SetParam(int* intarray,int M, int N, int enum_type){/*{{{*/

	Param* param=NULL;

	/*first, figure out if the param has already been created: */
	param=xDynamicCast<Param*>(this->FindParamObject(enum_type));

	if(param) param->SetValue(intarray,M,N); //already exists, just set it.
	else this->AddObject(new IntMatParam(enum_type,intarray,M,N)); //just add the new parameter.
}
/*}}}*/
void   Parameters::SetParam(Vector<IssmDouble>* vector,int enum_type){/*{{{*/

	Param* param=NULL;

	/*first, figure out if the param has already been created: */
	param=xDynamicCast<Param*>(this->FindParamObject(enum_type));

	if(param) param->SetValue(vector); //already exists, just set it.
	else this->AddObject(new VectorParam(enum_type,vector)); //just add the new parameter.
}
/*}}}*/
void   Parameters::SetParam(Matrix<IssmDouble>* matrix,int enum_type){/*{{{*/

	Param* param=NULL;

	/*first, figure out if the param has already been created: */
	param=xDynamicCast<Param*>(this->FindParamObject(enum_type));

	if(param) param->SetValue(matrix); //already exists, just set it.
	else this->AddObject(new MatrixParam(enum_type,matrix)); //just add the new parameter.
}
/*}}}*/
void   Parameters::SetParam(FILE* fid,int enum_type){/*{{{*/

	Param* param=NULL;

	/*first, figure out if the param has already been created: */
	param=xDynamicCast<Param*>(this->FindParamObject(enum_type));

	if(param) param->SetValue(fid); //already exists, just set it.
	else this->AddObject(new FileParam(enum_type,fid)); //just add the new parameter.
}
/*}}}*/
void   Parameters::SetParam(DataSet* dataset,int enum_type){/*{{{*/

	Param* param=NULL;

	/*first, figure out if the param has already been created: */
	param=xDynamicCast<Param*>(this->FindParamObject(enum_type));

	if(param) param->SetValue(dataset); //already exists, just set it.
	else this->AddObject(new DataSetParam(enum_type,dataset)); //just add the new parameter.
}
/*}}}*/

Param* Parameters::FindParamObject(int param_enum){/*{{{*/

	#ifdef _ISSM_DEBUG_
	if(param_enum<=ParametersSTARTEnum) _error_("Enum "<<EnumToStringx(param_enum)<<" should appear after ParametersSTARTEnum");
	if(param_enum>=ParametersENDEnum)   _error_("Enum "<<EnumToStringx(param_enum)<<" should appear before ParametersENDEnum");
	#endif

	int index = param_enum - ParametersSTARTEnum -1;
	return this->params[index];
}
/*}}}*/

/*Methods relating to parameters: */
char* OptionsFromAnalysis(Parameters* parameters,int analysis_type){ /*{{{*/

	/* figure out ISSM options for current analysis, return a string. */ 

	/*output: */
	char*   outstring=NULL;

	/*intermediary: */
	int          dummy;
	IssmDouble  *analyses    = NULL;
	char       **strings     = NULL;
	char        *string      = NULL;
	int          numanalyses;
	int          found       = -1;
	int          i;

	numanalyses=0;
	parameters->FindParam(&strings,&numanalyses,ToolkitsOptionsStringsEnum);

	parameters->FindParam(&analyses,&dummy,ToolkitsOptionsAnalysesEnum);

	if(numanalyses==0)return NULL; //we did not find petsc options, don't bother.

	/*ok, go through analyses and figure out if it corresponds to our analysis_type: */
	for(i=0;i<numanalyses;i++){
		if(analyses[i]==analysis_type){
			found=i;
			break;
		}
	}
	if(found==-1){
		/*still haven't found a list of petsc options, go find the default one, for analysis type DefaultAnalysisEnum: */
		for(i=0;i<numanalyses;i++){
			if(analyses[i]==DefaultAnalysisEnum){
				found=i;
				break;
			}
		}
	}
	if (found==-1){
		_error_("could find neither a default analysis nor analysis " << EnumToStringx(analysis_type));
	}

	/*ok, grab the option string: */
	outstring=xNew<char>(strlen(strings[found])+1);
	strcpy(outstring,strings[found]);

	/*Free ressources*/
	xDelete<IssmDouble>(analyses);
	for(i=0;i<numanalyses;i++){
		string=strings[i];
		xDelete<char>(string);
	}
	xDelete<char*>(strings);
	return outstring;
} 
/*}}}*/
void ToolkitsOptionsFromAnalysis(Parameters* parameters,int analysis_type){ /*{{{*/

	/*!\file:  ToolkitsOptionsFromAnalysis.cpp
	 * \brief: for each analysis, setup the issmoptions string. 
	 * This is mainly for the case where we run our toolkits using petsc. In this case, we need to 
	 * plug our toolkits options directly into the petsc options database. This is the case for each analysis type 
	 * and parameters
	 */ 

	char* options=NULL;

	/*Recover first the options string for this analysis: */
	options=OptionsFromAnalysis(parameters,analysis_type);

	/*Initialize our Toolkit Options: */
	ToolkitOptions::Init(options);

	#ifdef _HAVE_PETSC_
		/*In case we are using PETSC, we do not rely on issmoptions. Instead, we dump issmoptions into the Petsc 
		 * options database: */

		#if _PETSC_MAJOR_ == 2 
		PetscOptionsDestroy();
		PetscOptionsCreate();
		//PetscOptionsCheckInitial_Private();
		//PetscOptionsCheckInitial_Components();
		PetscOptionsSetFromOptions();
		PetscOptionsInsertMultipleString(options); //our patch
		#else
		#if (_PETSC_MINOR_>=7)
		PetscOptionsSetFromOptions(NULL);
		PetscOptionsClear(NULL);
		#else
		PetscOptionsSetFromOptions();
		PetscOptionsClear();
		#endif
		//PetscOptionsSetFromOptions();
		PetscOptionsInsertMultipleString(options); //our patch
		#endif

	#endif

	xDelete<char>(options);
}
/*}}}*/
