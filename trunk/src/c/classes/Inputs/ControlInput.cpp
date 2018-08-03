/*!\file ControlInput.c
 * \brief: implementation of the ControlInput object
 */

#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#include "../classes.h"
#include "../../shared/shared.h"

/*ControlInput constructors and destructor*/
ControlInput::ControlInput(){/*{{{*/
	control_id  = 0;
	values      = NULL;
	savedvalues = NULL;
	minvalues   = NULL;
	maxvalues   = NULL;
	gradient    = NULL;
}
/*}}}*/
ControlInput::ControlInput(int in_enum_type,int enum_input,IssmDouble* pvalues,IssmDouble* pmin,IssmDouble* pmax,int id){/*{{{*/

	control_id=id;
	enum_type=in_enum_type;

	switch(enum_input){
		case TriaInputEnum:
			values     =new TriaInput(enum_type,pvalues,P1Enum);
			savedvalues=new TriaInput(enum_type,pvalues,P1Enum);
			minvalues  =new TriaInput(enum_type,pmin,P1Enum);
			maxvalues  =new TriaInput(enum_type,pmax,P1Enum);
			break;
		case PentaInputEnum:
			values     =new PentaInput(enum_type,pvalues,P1Enum);
			savedvalues=new PentaInput(enum_type,pvalues,P1Enum);
			minvalues  =new PentaInput(enum_type,pmin,P1Enum);
			maxvalues  =new PentaInput(enum_type,pmax,P1Enum);
			break;
		default:
			_error_("Input of Enum " << EnumToStringx(enum_input) << " not supported yet by ControlInput");
	}
	gradient   =NULL;
}
/*}}}*/
ControlInput::~ControlInput(){/*{{{*/
	delete values;
	delete savedvalues;
	delete minvalues;
	delete maxvalues;
	delete gradient;
}
/*}}}*/

/*Object virtual functions definitions:*/
Object* ControlInput::copy() {/*{{{*/

	ControlInput* output=NULL;

	output = new ControlInput();
	output->enum_type=this->enum_type;
	output->control_id=this->control_id;

	if(values)      output->values      = xDynamicCast<Input*>(this->values->copy());
	if(savedvalues) output->savedvalues = xDynamicCast<Input*>(this->savedvalues->copy());
	if(minvalues)   output->minvalues   = xDynamicCast<Input*>(this->minvalues->copy());
	if(maxvalues)   output->maxvalues   = xDynamicCast<Input*>(this->maxvalues->copy());
	if(gradient)    output->gradient    = xDynamicCast<Input*>(this->gradient->copy());

	return output;
}
/*}}}*/
void ControlInput::DeepEcho(void){/*{{{*/

	_printf_("ControlInput:\n");
	_printf_(setw(15)<<"   ControlInput "<<setw(25)<<left<<EnumToStringx(this->enum_type)<<"\n");
	_printf_("---values: \n");     if (values)      values->Echo();
	_printf_("---savedvalues: \n");if (savedvalues) savedvalues->Echo();
	_printf_("---minvalues: \n");  if (minvalues)   minvalues->Echo();
	_printf_("---maxvalues: \n");  if (maxvalues)   maxvalues->Echo();
	_printf_("---gradient: \n");   if (gradient){    gradient->Echo();} else{_printf_("     Not set yet\n");}
}
/*}}}*/
void ControlInput::Echo(void){/*{{{*/
	this->DeepEcho();
}
/*}}}*/
int  ControlInput::Id(void){ return -1; }/*{{{*/
/*}}}*/
void ControlInput::Marshall(char** pmarshalled_data,int* pmarshalled_data_size, int marshall_direction){ /*{{{*/

	MARSHALLING_ENUM(ControlInputEnum);

	MARSHALLING(enum_type);
	MARSHALLING(control_id);

	if (marshall_direction == MARSHALLING_BACKWARD){
		switch(enum_type){
			case TriaInputEnum:
				values     =new TriaInput();
				savedvalues=new TriaInput();
				minvalues  =new TriaInput();
				maxvalues  =new TriaInput();
				gradient   =new TriaInput();
				break;
			case PentaInputEnum:
				values     =new PentaInput();
				savedvalues=new PentaInput();
				minvalues  =new PentaInput();
				maxvalues  =new PentaInput();
				gradient   =new PentaInput();
				break;
			default:
				_error_("Input of Enum " << EnumToStringx(enum_type) << " not supported yet by ControlInput");
		}
	}
	if(values) this->values->Marshall(pmarshalled_data,pmarshalled_data_size,marshall_direction);
	if(savedvalues) this->savedvalues->Marshall(pmarshalled_data,pmarshalled_data_size,marshall_direction);
	if(minvalues) this->minvalues->Marshall(pmarshalled_data,pmarshalled_data_size,marshall_direction);
	if(maxvalues) this->maxvalues->Marshall(pmarshalled_data,pmarshalled_data_size,marshall_direction);
	if(gradient) this->gradient->Marshall(pmarshalled_data,pmarshalled_data_size,marshall_direction);
}
/*}}}*/
int  ControlInput::ObjectEnum(void){/*{{{*/

	return ControlInputEnum;

}
/*}}}*/

/*ControlInput management*/
int ControlInput::InstanceEnum(void){/*{{{*/

	return this->enum_type;

}
/*}}}*/

/*Object functions*/
void ControlInput::AXPY(Input* xinput,IssmDouble scalar){/*{{{*/
	values->AXPY(xinput,scalar);
}/*}}}*/
void ControlInput::Configure(Parameters* parameters){/*{{{*/
	/*do nothing: */
}
/*}}}*/
void ControlInput::Constrain(void){/*{{{*/

	Input* newvalues=NULL;

	newvalues=this->values->PointwiseMin(maxvalues);
	delete values; this->values=newvalues;
	newvalues=this->values->PointwiseMax(minvalues);
	delete values; this->values=newvalues;
}/*}}}*/
void ControlInput::Constrain(IssmDouble min, IssmDouble max){/*{{{*/
	   values->Constrain(min,max);
}/*}}}*/
void ControlInput::Extrude(int start){/*{{{*/
	values->Extrude(start);
	savedvalues->Extrude(start);
	//gradient->Extrude();
}/*}}}*/
void ControlInput::GetGradient(Vector<IssmDouble>* gradient_vec,int* doflist){/*{{{*/
	if(gradient) gradient->GetVectorFromInputs(gradient_vec,doflist);
}/*}}}*/
void ControlInput::GetGradientValue(IssmDouble* pvalue,Gauss* gauss){/*{{{*/
	gradient->GetInputValue(pvalue,gauss);
}/*}}}*/
void ControlInput::GetInputAverage(IssmDouble* pvalue){/*{{{*/
	values->GetInputAverage(pvalue);
}/*}}}*/
void ControlInput::GetInputDerivativeValue(IssmDouble* derivativevalues, IssmDouble* xyz_list, Gauss* gauss){/*{{{*/
	values->GetInputDerivativeValue(derivativevalues,xyz_list,gauss);
}/*}}}*/
void ControlInput::GetInputValue(bool* pvalue){/*{{{*/
	values->GetInputValue(pvalue);
}/*}}}*/
void ControlInput::GetInputValue(int* pvalue){/*{{{*/
	values->GetInputValue(pvalue);
}/*}}}*/
void ControlInput::GetInputValue(IssmDouble* pvalue){/*{{{*/
	values->GetInputValue(pvalue);
}/*}}}*/
void ControlInput::GetInputValue(IssmDouble* pvalue,Gauss* gauss){/*{{{*/
	values->GetInputValue(pvalue,gauss);
}/*}}}*/
int  ControlInput::GetResultInterpolation(void){/*{{{*/

	return values->GetResultInterpolation();

}
/*}}}*/
int  ControlInput::GetResultNumberOfNodes(void){/*{{{*/

	return values->GetResultNumberOfNodes();

}
/*}}}*/
void ControlInput::GetVectorFromInputs(Vector<IssmDouble>* vector,int* doflist){/*{{{*/
	values->GetVectorFromInputs(vector,doflist);
}/*}}}*/
void ControlInput::GetVectorFromInputs(Vector<IssmDouble>* vector,int* doflist,const char* data){/*{{{*/
	 if(strcmp(data,"value")==0){
		 _assert_(values);
		 values->GetVectorFromInputs(vector,doflist);
	 }
	 else if (strcmp(data,"lowerbound")==0){
		 _assert_(minvalues);
		 minvalues->GetVectorFromInputs(vector,doflist);
	 }
	 else if (strcmp(data,"upperbound")==0){
		 _assert_(maxvalues);
		 maxvalues->GetVectorFromInputs(vector,doflist);
	 }
	 else if (strcmp(data,"gradient")==0){
		 _assert_(gradient);
		 gradient->GetVectorFromInputs(vector,doflist);
	 }
	 else{
		 _error_("Data " << data << " not supported yet");
	 }
}/*}}}*/
IssmDouble ControlInput::Min(void){/*{{{*/

	return values->Min();

}
/*}}}*/
void ControlInput::SaveValue(void){/*{{{*/
	if(!values) _error_("Values of " << EnumToStringx(this->enum_type) << " not found");

	if(savedvalues) delete this->savedvalues;
	this->savedvalues=xDynamicCast<Input*>(this->values->copy());
}/*}}}*/
void ControlInput::ScaleGradient(IssmDouble scaling_factor){/*{{{*/
	if(!gradient) _error_("Gradient of ControlInput " << EnumToStringx(enum_type) << " not found");
	gradient->Scale(scaling_factor);
}/*}}}*/
void ControlInput::SetGradient(Input* gradient_in){/*{{{*/

	/*Get enum for current gradient*/
	switch(this->control_id){
		case 1:
			gradient_in->ChangeEnum(Gradient1Enum);
			break;
		case 2:
			gradient_in->ChangeEnum(Gradient2Enum);
			break;
		case 3:
			gradient_in->ChangeEnum(Gradient3Enum);
			break;
		default:
			_error_("more than 3 controls not implemented yet (Gradient " << this->control_id << " was requested). EnumDefinitions.h needs to be updated.");
	}

	/*Delete old gradient and assign new gradient*/
	if(gradient) delete gradient;
	gradient=gradient_in;

}/*}}}*/
void ControlInput::SetInput(Input* in_input){/*{{{*/

	delete values; this->values=in_input;
	this->SaveValue(); //because this is what SpawnResult saves FIXME

}/*}}}*/
Input* ControlInput::SpawnSegInput(int index1,int index2){/*{{{*/
	return values->SpawnSegInput(index1,index2);
}/*}}}*/
Input* ControlInput::SpawnTriaInput(int index1,int index2,int index3){/*{{{*/
	return values->SpawnTriaInput(index1,index2,index3);
}/*}}}*/
void ControlInput::UpdateValue(IssmDouble scalar){/*{{{*/
	if(!gradient)    _error_("Gradient of " << EnumToStringx(this->enum_type) << " not found");
	if(!savedvalues) _error_("Values of " << EnumToStringx(this->enum_type) << " not found");

	if(values) delete this->values;
	this->values=xDynamicCast<Input*>(this->savedvalues->copy());
	this->values->AXPY(gradient,scalar);
}/*}}}*/
void ControlInput::VerticallyIntegrate(Input* thickness_input){/*{{{*/
	values->VerticallyIntegrate(thickness_input);
}/*}}}*/
