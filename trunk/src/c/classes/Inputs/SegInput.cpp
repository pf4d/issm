/*!\file SegInput.c
 * \brief: implementation of the SegInput object
 */

#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#include "../classes.h"
#include "../../shared/shared.h"

/*SegInput constructors and destructor*/
SegInput::SegInput(){/*{{{*/
	values = NULL;
}
/*}}}*/
SegInput::SegInput(int in_enum_type,IssmDouble* in_values,int interpolation_type_in){/*{{{*/

	/*Set Enum*/
	enum_type=in_enum_type;
	this->interpolation_type=interpolation_type_in;

	/*Set values*/
	this->values=xNew<IssmDouble>(this->NumberofNodes(this->interpolation_type));
	for(int i=0;i<this->NumberofNodes(this->interpolation_type);i++) values[i]=in_values[i];
}
/*}}}*/
SegInput::~SegInput(){/*{{{*/
	xDelete<IssmDouble>(this->values);
}
/*}}}*/

/*Object virtual functions definitions:*/
Object* SegInput::copy() {/*{{{*/

	return new SegInput(this->enum_type,this->values,this->interpolation_type);

}
/*}}}*/
void SegInput::DeepEcho(void){/*{{{*/

	_printf_(setw(15)<<"   SegInput "<<setw(25)<<left<<EnumToStringx(this->enum_type)<<" [");
	for(int i=0;i<this->NumberofNodes(this->interpolation_type);i++) _printf_(" "<<this->values[i]);
	_printf_("]\n");
}
/*}}}*/
void SegInput::Echo(void){/*{{{*/
	this->DeepEcho();
}
/*}}}*/
int  SegInput::Id(void){ return -1; }/*{{{*/
/*}}}*/
void SegInput::Marshall(char** pmarshalled_data,int* pmarshalled_data_size, int marshall_direction){ /*{{{*/

	MARSHALLING_ENUM(SegInputEnum);

	MARSHALLING(enum_type);
	MARSHALLING(interpolation_type);

	int numnodes = this->NumberofNodes(this->interpolation_type);
	if(numnodes > 0){
		MARSHALLING_DYNAMIC(this->values,IssmDouble,numnodes)
	}
	else this->values = NULL;
}
/*}}}*/
int  SegInput::ObjectEnum(void){/*{{{*/

	return SegInputEnum;

}
/*}}}*/

/*SegInput management*/
int SegInput::InstanceEnum(void){/*{{{*/

	return this->enum_type;

}
/*}}}*/

/*Object functions*/
void SegInput::Configure(Parameters* parameters){/*{{{*/
	/*do nothing: */
}
/*}}}*/
void SegInput::GetInputAverage(IssmDouble* pvalue){/*{{{*/

	int        numnodes  = this->NumberofNodes(this->interpolation_type);
	IssmDouble numnodesd = reCast<int,IssmDouble>(numnodes);
	IssmDouble value     = 0.;

	for(int i=0;i<numnodes;i++) value+=values[i];
	value = value/numnodesd;

	*pvalue=value;
}
/*}}}*/
void SegInput::GetInputDerivativeValue(IssmDouble* p, IssmDouble* xyz_list,Gauss* gauss){/*{{{*/

	/*Call SegRef function*/
	_assert_(gauss->Enum()==GaussSegEnum);
	SegRef::GetInputDerivativeValue(p,&values[0],xyz_list,(GaussSeg*)gauss,this->interpolation_type);
}
/*}}}*/
void SegInput::GetInputValue(IssmDouble* pvalue,Gauss* gauss){/*{{{*/

	/*Call SegRef function*/
	_assert_(gauss->Enum()==GaussSegEnum);
	SegRef::GetInputValue(pvalue,&values[0],(GaussSeg*)gauss,this->interpolation_type);

}
/*}}}*/
IssmDouble SegInput::Min(void){/*{{{*/

	const int  numnodes=this->NumberofNodes(this->interpolation_type);
	IssmDouble min=values[0];

	for(int i=1;i<numnodes;i++){
		if(values[i]<min) min=values[i];
	}
	return min;
}
/*}}}*/
