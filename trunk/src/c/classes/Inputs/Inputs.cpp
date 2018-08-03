/*
 * \file Inputs.c
 * \brief: implementation of the Inputs class, derived from DataSet class
 */

/*Headers: {{{*/
#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#include "./Input.h"
#include "./Inputs.h"
#include "../../shared/shared.h"

using namespace std;
/*}}}*/

/*Object constructors and destructor*/
Inputs::Inputs(){/*{{{*/
	return;
}
/*}}}*/
Inputs::~Inputs(){/*{{{*/
	return;
}
/*}}}*/

/*Object management*/
int  Inputs::AddInput(Input* in_input){/*{{{*/

	/*First, go through dataset of inputs and check whether any input 
	 * with the same name is already in. If so, erase the corresponding 
	 * object before adding this new one: */
	vector<Object*>::iterator object;
	Input* input=NULL;

	/*In debugging mode, check that the input is not a NULL pointer*/
	_assert_(in_input);

	for ( object=objects.begin() ; object < objects.end(); object++ ){

		input=xDynamicCast<Input*>(*object);

		if (input->InstanceEnum()==in_input->InstanceEnum()){
			this->DeleteObject(input);
			break;
		}
	}
	this->AddObject(in_input);

	return 1;
}
/*}}}*/
void  Inputs::AXPY(int inputy_enum, IssmDouble scalar, int inputx_enum){/*{{{*/

	/*Find x and y inputs: */
	Input* xinput=xDynamicCast<Input*>(this->GetInput(inputx_enum));
	Input* yinput=xDynamicCast<Input*>(this->GetInput(inputy_enum));

	/*some checks: */
	if(!xinput) _error_("input " << EnumToStringx(inputx_enum) << " could not be found!");
	if(!yinput) _error_("input " << EnumToStringx(inputy_enum) << " could not be found!");

	/*Apply AXPY: */
	yinput->AXPY(xinput,scalar);
}
/*}}}*/
void  Inputs::ChangeEnum(int oldenumtype,int newenumtype){/*{{{*/

	/*Go through dataset of inputs and look for input with 
	 * same enum as input enum, once found, just change its name */
	vector<Object*>::iterator object;
	Input* input=NULL;

	/*Delete existing input of newenumtype if it exists*/
	for ( object=objects.begin() ; object < objects.end(); object++ ){
		input=xDynamicCast<Input*>(*object);

		if (input->InstanceEnum()==newenumtype){
			this->DeleteObject(input);
			break;
		}
	}

	/*Change enum_type of input of oldenumtype*/
	for ( object=objects.begin() ; object < objects.end(); object++ ){

		input=xDynamicCast<Input*>(*object);

		if (input->InstanceEnum()==oldenumtype){
			input->ChangeEnum(newenumtype);
			break;
		}
	}
}
/*}}}*/
void Inputs::Configure(Parameters* parameters){/*{{{*/

	vector<Object*>::iterator object;
	Input* input=NULL;

	for ( object=objects.begin() ; object < objects.end(); object++ ){

		input=xDynamicCast<Input*>(*object);
		input->Configure(parameters);

	}

}
/*}}}*/
void  Inputs::ConstrainMin(int constrain_enum, IssmDouble minimum){/*{{{*/

	/*Find x and y inputs: */
	Input* constrain_input=xDynamicCast<Input*>(this->GetInput(constrain_enum));

	/*some checks: */
	if(!constrain_input) _error_("input " << EnumToStringx(constrain_enum) << " could not be found!");

	/*Apply ContrainMin: */
	constrain_input->ConstrainMin(minimum);
}
/*}}}*/
int  Inputs::DeleteInput(int enum_type){/*{{{*/

	vector<Object*>::iterator object;
	Input* input=NULL;

	for ( object=objects.begin() ; object < objects.end(); object++ ){

		input=xDynamicCast<Input*>(*object);

		if (input->InstanceEnum()==enum_type){
			this->DeleteObject(input);
			break;
		}
	}

	return 1;

}
/*}}}*/
void  Inputs::DuplicateInput(int original_enum,int new_enum){/*{{{*/

	/*Make a copy of the original input: */
	Input* original=xDynamicCast<Input*>(this->GetInput(original_enum));
	if(!original)_error_("could not find input with enum: " << EnumToStringx(original_enum)); 
	Input* copy=xDynamicCast<Input*>(original->copy());

	/*Change copy enum to reinitialized_enum: */
	copy->ChangeEnum(new_enum);

	/*Add copy into inputs, it will wipe off the one already there: */
	this->AddInput(xDynamicCast<Input*>(copy));
}
/*}}}*/
Input* Inputs::GetInput(int enum_name){/*{{{*/

	vector<Object*>::iterator object;
	Input* input=NULL;

	for ( object=objects.begin() ; object < objects.end(); object++ ){

		input=xDynamicCast<Input*>(*object);

		if (input->InstanceEnum()==enum_name){
			return input;
		}
	}
	return NULL;
}
/*}}}*/
void Inputs::GetInputAverage(IssmDouble* pvalue,int enum_type){/*{{{*/

	vector<Object*>::iterator object;
	Input* input=NULL;
	bool   found=false;

	/*Go through inputs and check whether any input with the same name is already in: */
	for ( object=objects.begin() ; object < objects.end(); object++ ){

		input=xDynamicCast<Input*>(*object);
		if (input->InstanceEnum()==enum_type){
			found=true;
			break;
		}
	}

	if (!found){
		/*we could not find an input with the correct enum type. No defaults values were provided, 
		 * error out: */
		_error_("could not find input with enum type " << enum_type << " (" << EnumToStringx(enum_type) << ")");
	}

	/*Ok, we have an input if we made it here, request the input to return the value: */
	input->GetInputAverage(pvalue);

}
/*}}}*/
void Inputs::GetInputValue(bool* pvalue,int enum_type){/*{{{*/

	vector<Object*>::iterator object;
	Input* input=NULL;
	bool   found=false;

	/*Go through inputs and check whether any input with the same name is already in: */
	for ( object=objects.begin() ; object < objects.end(); object++ ){

		input=xDynamicCast<Input*>(*object);
		if (input->InstanceEnum()==enum_type){
			found=true;
			break;
		}
	}

	if (!found){
		/*we could not find an input with the correct enum type. No defaults values were provided, 
		 * error out: */
		_error_("could not find input with enum type " << enum_type << " (" << EnumToStringx(enum_type) << ")");
	}

	/*Ok, we have an input if we made it here, request the input to return the value: */
	input->GetInputValue(pvalue);

}
/*}}}*/
void Inputs::GetInputValue(int* pvalue,int enum_type){/*{{{*/

	vector<Object*>::iterator object;
	Input* input=NULL;
	bool   found=false;

	/*Go through inputs and check whether any input with the same name is already in: */
	for ( object=objects.begin() ; object < objects.end(); object++ ){

		input=xDynamicCast<Input*>(*object);
		if (input->InstanceEnum()==enum_type){
			found=true;
			break;
		}
	}

	if (!found){
		/*we could not find an input with the correct enum type. No defaults values were provided, 
		 * error out: */
		_error_("could not find input with enum type " << enum_type << " (" << EnumToStringx(enum_type) << ")");
	}

	/*Ok, we have an input if we made it here, request the input to return the value: */
	input->GetInputValue(pvalue);

}
/*}}}*/
void Inputs::GetInputValue(IssmDouble* pvalue,int enum_type){/*{{{*/

	vector<Object*>::iterator object;
	Input* input=NULL;
	bool   found=false;

	/*Go through inputs and check whether any input with the same name is already in: */
	for ( object=objects.begin() ; object < objects.end(); object++ ){

		input=xDynamicCast<Input*>(*object); 
		if (input->InstanceEnum()==enum_type){
			found=true;
			break;
		}
	}

	if (!found){
		/*we could not find an input with the correct enum type. No defaults values were provided, 
		 * error out: */
		_error_("could not find input with enum type " << enum_type << " (" << EnumToStringx(enum_type) << ")");
	}

	/*Ok, we have an input if we made it here, request the input to return the value: */
	input->GetInputValue(pvalue);

}
/*}}}*/
IssmDouble Inputs::InfinityNorm(int enumtype){/*{{{*/

	/*Output*/
	IssmDouble norm;

	/*Get input*/
	Input* input=xDynamicCast<Input*>(this->GetInput(enumtype));

	/*Apply ContrainMin: */
	if (input){
		norm=input->InfinityNorm();
	}
	else{
		norm=0;
	}

	/*Return output*/
	return norm;
}
/*}}}*/
IssmDouble Inputs::Max(int enumtype){/*{{{*/

	/*Output*/
	IssmDouble max;

	/*Get input*/
	Input* input=xDynamicCast<Input*>(this->GetInput(enumtype));

	/*Apply ContrainMin: */
	if (input){
		max=input->Max();
	}
	else{
		_error_("Input " << EnumToStringx(enumtype) << " not found");
	}

	/*Return output*/
	return max;
}
/*}}}*/
IssmDouble Inputs::MaxAbs(int enumtype){/*{{{*/

	/*Output*/
	IssmDouble max;

	/*Get input*/
	Input* input=xDynamicCast<Input*>(this->GetInput(enumtype));

	/*Apply ContrainMin: */
	if (input){
		max=input->MaxAbs();
	}
	else{
		_error_("Input " << EnumToStringx(enumtype) << " not found");
	}

	/*Return output*/
	return max;
}
/*}}}*/
IssmDouble Inputs::Min(int enumtype){/*{{{*/

	/*Output*/
	IssmDouble min;

	/*Get input*/
	Input* input=xDynamicCast<Input*>(this->GetInput(enumtype));

	/*Apply ContrainMin: */
	if (input){
		min=input->Min();
	}
	else{
		_error_("Input " << EnumToStringx(enumtype) << " not found");
	}

	/*Return output*/
	return min;
}
/*}}}*/
IssmDouble Inputs::MinAbs(int enumtype){/*{{{*/

	/*Output*/
	IssmDouble min;

	/*Get input*/
	Input* input=xDynamicCast<Input*>(this->GetInput(enumtype));

	/*Apply ContrainMin: */
	if (input){
		min=input->MinAbs();
	}
	else{
		_error_("Input " << EnumToStringx(enumtype) << " not found");
	}

	/*Return output*/
	return min;
}
/*}}}*/
Inputs* Inputs::SpawnSegInputs(int index1,int index2){/*{{{*/

	/*Intermediary*/
	vector<Object*>::iterator object;
	Input* inputin=NULL;
	Input* inputout=NULL;

	/*Output*/
	Inputs* newinputs=new Inputs();

	/*Go through inputs and call Spawn function*/
	for ( object=objects.begin() ; object < objects.end(); object++ ){

		/*Create new input*/
		inputin=xDynamicCast<Input*>(*object);
		inputout=inputin->SpawnSegInput(index1,index2);

		/*Add input to new inputs*/
		newinputs->AddObject(inputout);
	}

	/*Assign output pointer*/
	return newinputs;
}
/*}}}*/
Inputs* Inputs::SpawnTriaInputs(int index1,int index2,int index3){/*{{{*/

	/*Intermediary*/
	vector<Object*>::iterator object;
	Input* inputin=NULL;
	Input* inputout=NULL;

	/*Output*/
	Inputs* newinputs=new Inputs();

	/*Go through inputs and call Spawn function*/
	for ( object=objects.begin() ; object < objects.end(); object++ ){

		/*Create new input*/
		inputin=xDynamicCast<Input*>(*object);
		inputout=inputin->SpawnTriaInput(index1,index2,index3);

		/*Add input to new inputs*/
		newinputs->AddObject(inputout);
	}

	/*Assign output pointer*/
	return newinputs;
}
/*}}}*/
