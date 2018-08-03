/*!\file BoolInput.c
 * \brief: implementation of the BoolInput object
 */

#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#include "../classes.h"
#include "../../shared/shared.h"

/*BoolInput constructors and destructor*/
BoolInput::BoolInput(){/*{{{*/
	return;
}
/*}}}*/
BoolInput::BoolInput(int in_enum_type,bool in_value){/*{{{*/

	enum_type=in_enum_type;
	value=in_value;
}
/*}}}*/
BoolInput::~BoolInput(){/*{{{*/
	return;
}
/*}}}*/

/*Object virtual functions definitions:*/
Object* BoolInput::copy() {/*{{{*/

	return new BoolInput(this->enum_type,this->value);

}
/*}}}*/
void BoolInput::DeepEcho(void){/*{{{*/

	_printf_(setw(15)<<"   BoolInput "<<setw(25)<<left<<EnumToStringx(this->enum_type)<<" "<<(value?"true":"false") << "\n");
}
/*}}}*/
void BoolInput::Echo(void){/*{{{*/
	this->DeepEcho();
}
/*}}}*/
int  BoolInput::Id(void){ return -1; }/*{{{*/
/*}}}*/
void BoolInput::Marshall(char** pmarshalled_data,int* pmarshalled_data_size, int marshall_direction){ /*{{{*/

	MARSHALLING_ENUM(BoolInputEnum);

	MARSHALLING(enum_type);
	MARSHALLING(value);

}
/*}}}*/
int  BoolInput::ObjectEnum(void){/*{{{*/

	return BoolInputEnum;

}
/*}}}*/

/*BoolInput management*/
int BoolInput::InstanceEnum(void){/*{{{*/

	return this->enum_type;

}
/*}}}*/
Input* BoolInput::SpawnSegInput(int index1,int index2){/*{{{*/

		/*output*/
		BoolInput* outinput=new BoolInput();

		/*only copy current value*/
		outinput->enum_type=this->enum_type;
		outinput->value=this->value;

		/*Assign output*/
		return outinput;

}
/*}}}*/
Input* BoolInput::SpawnTriaInput(int index1,int index2,int index3){/*{{{*/

		/*output*/
		BoolInput* outinput=new BoolInput();

		/*only copy current value*/
		outinput->enum_type=this->enum_type;
		outinput->value=this->value;

		/*Assign output*/
		return outinput;

}
/*}}}*/

/*Object functions*/
void BoolInput::AXPY(Input* xinput,IssmDouble scalar){/*{{{*/

	BoolInput*  xboolinput=NULL;

	/*xinput is of the same type, so cast it: */
	xboolinput=(BoolInput*)xinput;

	/*Carry out the AXPY operation depending on type:*/
	switch(xinput->ObjectEnum()){

		case BoolInputEnum:
			this->value=reCast<bool,IssmDouble>(this->value+scalar*xboolinput->value);
			return;

		default:
			_error_("not implemented yet");
	}

}
/*}}}*/
void BoolInput::ChangeEnum(int newenumtype){/*{{{*/
	this->enum_type=newenumtype;
}
/*}}}*/
void BoolInput::Configure(Parameters* parameters){/*{{{*/
	/*do nothing: */
}
/*}}}*/
void BoolInput::Extrude(int start){/*{{{*/

	/*do nothing*/
	return;
}
/*}}}*/
void BoolInput::GetInputAverage(IssmDouble* pvalue){/*{{{*/
	*pvalue=reCast<IssmDouble>(value);
}
/*}}}*/
void BoolInput::GetInputValue(bool* pvalue){/*{{{*/
	*pvalue=value;
}
/*}}}*/
void BoolInput::GetInputValue(int* pvalue){_error_("not supported yet!");}/*{{{*/
/*}}}*/
void BoolInput::GetInputValue(IssmDouble* pvalue){_error_("not supported yet!");}/*{{{*/
/*}}}*/
void BoolInput::GetInputValue(IssmDouble* pvalue,Gauss* gauss){_error_("not supported yet!");}/*{{{*/
/*}}}*/
void BoolInput::GetVectorFromInputs(Vector<IssmDouble>* vector,int* doflist){/*{{{*/

	_error_("not supporte yet!");

}
/*}}}*/
void BoolInput::Scale(IssmDouble scale_factor){/*{{{*/
	/*a bool cannot be scaled: */
}
/*}}}*/
void BoolInput::SquareMin(IssmDouble* psquaremin,Parameters* parameters){/*{{{*/
	/*square of a bool is the bool itself: */
	*psquaremin=value;
}
/*}}}*/
