/*!\file IntInput.c
 * \brief: implementation of the IntInput object
 */

#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#include "../classes.h"
#include "../../shared/shared.h"

/*IntInput constructors and destructor*/
IntInput::IntInput(){/*{{{*/
	return;
}
/*}}}*/
IntInput::IntInput(int in_enum_type,IssmInt in_value){/*{{{*/

	enum_type=in_enum_type;
	value=in_value;
}
/*}}}*/
IntInput::~IntInput(){/*{{{*/
	return;
}
/*}}}*/

/*Object virtual functions definitions:*/
Object* IntInput::copy() {/*{{{*/

	return new IntInput(this->enum_type,this->value);

}
/*}}}*/
void IntInput::DeepEcho(void){/*{{{*/

	_printf_(setw(15)<<"   IntInput "<<setw(25)<<left<<EnumToStringx(this->enum_type)<<" "<<this->value<<"\n");
}
/*}}}*/
int  IntInput::Id(void){ return -1; }/*{{{*/
/*}}}*/
void IntInput::Marshall(char** pmarshalled_data,int* pmarshalled_data_size, int marshall_direction){ /*{{{*/

	MARSHALLING_ENUM(IntInputEnum);

	MARSHALLING(enum_type);
	MARSHALLING(value);

}
/*}}}*/
int  IntInput::ObjectEnum(void){/*{{{*/

	return IntInputEnum;

}
/*}}}*/

/*IntInput management*/
void IntInput::Echo(void){/*{{{*/
	this->DeepEcho();
}
/*}}}*/
int IntInput::InstanceEnum(void){/*{{{*/

	return this->enum_type;

}
/*}}}*/
Input* IntInput::SpawnSegInput(int index1,int index2){/*{{{*/

	/*output*/
	IntInput* outinput=new IntInput();

	/*only copy current value*/
	outinput->enum_type=this->enum_type;
	outinput->value=this->value;

	/*Assign output*/
	return outinput;
}
/*}}}*/
Input* IntInput::SpawnTriaInput(int index1,int index2,int index3){/*{{{*/

	/*output*/
	IntInput* outinput=new IntInput();

	/*only copy current value*/
	outinput->enum_type=this->enum_type;
	outinput->value=this->value;

	/*Assign output*/
	return outinput;
}
/*}}}*/

/*Object functions*/
void IntInput::AXPY(Input* xinput,IssmDouble scalar){/*{{{*/

	IssmDouble dvalue;
	IntInput*  xintinput=NULL;

	/*xinput is of the same type, so cast it: */
	xintinput=(IntInput*)xinput;

	/*Carry out the AXPY operation depending on type:*/
	switch(xinput->ObjectEnum()){

		case IntInputEnum:
			dvalue=(IssmDouble)this->value+scalar*(IssmDouble)xintinput->value;
			this->value=reCast<int>(dvalue);
			return;

		default:
			_error_("not implemented yet");
	}

}
/*}}}*/
void IntInput::ChangeEnum(int newenumtype){/*{{{*/
	this->enum_type=newenumtype;
}
/*}}}*/
void IntInput::Configure(Parameters* parameters){/*{{{*/
	/*do nothing: */
}
/*}}}*/
void IntInput::Constrain(IssmDouble cm_min, IssmDouble cm_max){/*{{{*/

	if(!xIsNan<IssmDouble>(cm_min)) if (this->value<cm_min)this->value=reCast<int>(cm_min);
	if(!xIsNan<IssmDouble>(cm_max)) if (this->value>cm_max)this->value=reCast<int>(cm_max);

}
/*}}}*/
void IntInput::GetInputAverage(IssmDouble* pvalue){/*{{{*/
	*pvalue=reCast<IssmDouble>(value);
}
/*}}}*/
void IntInput::GetInputValue(bool* pvalue){_error_("not supported yet!");}/*{{{*/
/*}}}*/
void IntInput::GetInputValue(int* pvalue){/*{{{*/
	*pvalue=value;
}
/*}}}*/
void IntInput::GetInputValue(IssmDouble* pvalue){/*{{{*/
	_error_("IntInput cannot return a IssmDouble in parallel");
}
/*}}}*/
void IntInput::GetInputValue(IssmDouble* pvalue,Gauss* gauss){_error_("not supported yet!");}/*{{{*/
/*}}}*/
void IntInput::GetVectorFromInputs(Vector<IssmDouble>* vector,int* doflist){/*{{{*/

	_error_("not supporte yet!");

}
/*}}}*/
void IntInput::Scale(IssmDouble scale_factor){/*{{{*/
	IssmDouble dvalue=(IssmDouble)value*scale_factor;
	value=reCast<int>(dvalue);
}
/*}}}*/
void IntInput::SquareMin(IssmDouble* psquaremin,Parameters* parameters){/*{{{*/

	/*square min of an integer is the square of the integer itself: */
	*psquaremin=pow((IssmDouble)value,2);
}
/*}}}*/
