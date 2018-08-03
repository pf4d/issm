/*!\file DoubleInput.c
 * \brief: implementation of the DoubleInput object
 */

#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#include "../classes.h"
#include "../../shared/shared.h"

/*DoubleInput constructors and destructor*/
DoubleInput::DoubleInput(){/*{{{*/
	return;
}
/*}}}*/
DoubleInput::DoubleInput(int in_enum_type,IssmDouble in_value){/*{{{*/

	enum_type=in_enum_type;
	value=in_value;
}
/*}}}*/
DoubleInput::~DoubleInput(){/*{{{*/
	return;
}
/*}}}*/

/*Object virtual functions definitions:*/
Object* DoubleInput::copy() {/*{{{*/

	return new DoubleInput(this->enum_type,this->value);

}
/*}}}*/
void DoubleInput::DeepEcho(void){/*{{{*/

	_printf_(setw(15)<<"   DoubleInput "<<setw(25)<<left<<EnumToStringx(this->enum_type)<<" "<<this->value<<"\n");
}
/*}}}*/
void DoubleInput::Echo(void){/*{{{*/
	this->DeepEcho();
}
/*}}}*/
int    DoubleInput::Id(void){ return -1; }/*{{{*/
/*}}}*/
void DoubleInput::Marshall(char** pmarshalled_data,int* pmarshalled_data_size, int marshall_direction){ /*{{{*/

	MARSHALLING_ENUM(DoubleInputEnum);

	MARSHALLING(enum_type);
	MARSHALLING(value);

}
/*}}}*/
int DoubleInput::ObjectEnum(void){/*{{{*/

	return DoubleInputEnum;

}
/*}}}*/

/*DoubleInput management*/
int DoubleInput::InstanceEnum(void){/*{{{*/

	return this->enum_type;

}
/*}}}*/
Input* DoubleInput::SpawnSegInput(int index1,int index2){/*{{{*/

	/*output*/
	DoubleInput* outinput=new DoubleInput();

	/*only copy current value*/
	outinput->enum_type=this->enum_type;
	outinput->value=this->value;

	/*Assign output*/
	return outinput;

}
/*}}}*/
Input* DoubleInput::SpawnTriaInput(int index1,int index2,int index3){/*{{{*/

	/*output*/
	DoubleInput* outinput=new DoubleInput();

	/*only copy current value*/
	outinput->enum_type=this->enum_type;
	outinput->value=this->value;

	/*Assign output*/
	return outinput;

}
/*}}}*/

/*Object functions*/
void DoubleInput::AXPY(Input* xinput,IssmDouble scalar){/*{{{*/

	DoubleInput*  xIssmDoubleinput=NULL;

	/*xinput is of the same type, so cast it: */
	xIssmDoubleinput=(DoubleInput*)xinput;

	/*Carry out the AXPY operation depending on type:*/
	switch(xinput->ObjectEnum()){

		case DoubleInputEnum:
			this->value=this->value+scalar*xIssmDoubleinput->value;
			return;

		default:
			_error_("not implemented yet");
	}

}
/*}}}*/
void DoubleInput::ChangeEnum(int newenumtype){/*{{{*/
	this->enum_type=newenumtype;
}
/*}}}*/
void DoubleInput::Configure(Parameters* parameters){/*{{{*/
	/*do nothing: */
}
/*}}}*/
void DoubleInput::Constrain(IssmDouble cm_min, IssmDouble cm_max){/*{{{*/

	if(!xIsNan<IssmDouble>(cm_min)) if (this->value<cm_min)this->value=cm_min;
	if(!xIsNan<IssmDouble>(cm_max)) if (this->value>cm_max)this->value=cm_max;

}
/*}}}*/
void DoubleInput::ConstrainMin(IssmDouble minimum){/*{{{*/
	if (value<minimum) value=minimum;
}
/*}}}*/
void DoubleInput::GetInputAverage(IssmDouble* pvalue){/*{{{*/
	*pvalue=value;
}
/*}}}*/
void DoubleInput::GetInputValue(bool* pvalue){/*{{{*/
	_error_("Double input of enum " << EnumToStringx(enum_type) << " cannot return a boolean");

}
/*}}}*/
void DoubleInput::GetInputValue(int* pvalue){/*{{{*/
	_error_("Double input of enum " << enum_type << " (" << EnumToStringx(enum_type) << ") cannot return an integer");

}
/*}}}*/
void DoubleInput::GetInputValue(IssmDouble* pvalue){/*{{{*/

	/*return value*/
	*pvalue=value;
}
/*}}}*/
void DoubleInput::GetInputValue(IssmDouble* pvalue,Gauss* gauss){*pvalue=this->value;}/*{{{*/
/*}}}*/
void DoubleInput::GetVectorFromInputs(Vector<IssmDouble>* vector,int* doflist){/*{{{*/

	_error_("not supporte yet!");

}
/*}}}*/
IssmDouble DoubleInput::Max(void){/*{{{*/
	return this->value;
}
/*}}}*/
IssmDouble DoubleInput::MaxAbs(void){/*{{{*/
	return fabs(this->value);
}
/*}}}*/
IssmDouble DoubleInput::Min(void){/*{{{*/
	return this->value;
}
/*}}}*/
IssmDouble DoubleInput::MinAbs(void){/*{{{*/
	return fabs(this->value);
}
/*}}}*/
Input* DoubleInput::PointwiseDivide(Input* inputB){/*{{{*/

	/*Ouput*/
	DoubleInput* outinput=NULL;

	/*Intermediaries*/
	IssmDouble       Bvalue;

	/*Check that inputB is of the same type*/
	inputB->GetInputAverage(&Bvalue);

	/*Create new DoubleInput*/
	outinput=new DoubleInput(this->enum_type,this->value/Bvalue);

	/*Return output pointer*/
	return outinput;

}
/*}}}*/
Input* DoubleInput::PointwiseMax(Input* input){/*{{{*/

	/*Ouput*/
	DoubleInput* outinput=NULL;

	/*Intermediaries*/
	IssmDouble       max;

	/*Check that inputB is of the same type*/
	if (input->Max() > this->Max()) max=input->Max();
	else max=this->Max();

	/*Create new DoubleInput*/
	outinput=new DoubleInput(this->enum_type,max);

	/*Return output pointer*/
	return outinput;

}
/*}}}*/
Input* DoubleInput::PointwiseMin(Input* input){/*{{{*/

	/*Ouput*/
	DoubleInput* outinput=NULL;

	/*Intermediaries*/
	IssmDouble       min;

	/*Check that inputB is of the same type*/
	if (input->Min() < this->Min()) min=input->Min();
	else min=this->Min();

	/*Create new DoubleInput*/
	outinput=new DoubleInput(this->enum_type,min);

	/*Return output pointer*/
	return outinput;

}
/*}}}*/
void DoubleInput::Scale(IssmDouble scale_factor){/*{{{*/
	value=value*scale_factor;
}
/*}}}*/
void DoubleInput::SquareMin(IssmDouble* psquaremin,Parameters* parameters){/*{{{*/

	/*square min of a IssmDouble is the square of the IssmDouble itself: */
	*psquaremin=pow(value,2);
}
/*}}}*/
void DoubleInput::VerticallyIntegrate(Input* thickness_input){/*{{{*/

	/*Intermediaries*/
	IssmDouble thickness_value;

	/*Check that input provided is a thickness*/
	if (thickness_input->InstanceEnum()!=ThicknessEnum) _error_("Input provided is not a Thickness (enum_type is " << EnumToStringx(thickness_input->InstanceEnum()) << ")");

	/*vertically integrate depending on type:*/
	switch(thickness_input->ObjectEnum()){

		case PentaInputEnum:
			thickness_input->GetInputAverage(&thickness_value);
			this->value=this->value*thickness_value;
			return;

		default:
			_error_("not implemented yet");
	}
}
/*}}}*/
