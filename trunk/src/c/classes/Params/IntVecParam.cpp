/*!\file IntVecParam.c
 * \brief: implementation of the IntVecParam object
 */

/*header files: */
/*{{{*/
#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#include "../classes.h"
#include "../../shared/shared.h"
/*}}}*/

/*IntVecParam constructors and destructor*/
IntVecParam::IntVecParam(){/*{{{*/
	return;
}
/*}}}*/
IntVecParam::IntVecParam(int in_enum_type,int* in_values, int in_M){/*{{{*/

	enum_type=in_enum_type;
	M=in_M;

	if(M){
		values=xNew<int>(M);
		xMemCpy<int>(values,in_values,M);
	}
	else values=NULL;
}
/*}}}*/
IntVecParam::IntVecParam(int in_enum_type,IssmDouble* in_values, int in_M){/*{{{*/

	enum_type=in_enum_type;
	M=in_M;

	if(M){
		values=xNew<int>(M);
		for(int i=0;i<in_M;i++) values[i]=reCast<int>(in_values[i]);
	}
	else values=NULL;
}
/*}}}*/
IntVecParam::~IntVecParam(){/*{{{*/
	xDelete<int>(values);
	return;
}
/*}}}*/

/*Object virtual functions definitions:*/
Param* IntVecParam::copy() {/*{{{*/

	return new IntVecParam(this->enum_type,this->values,this->M);

}
/*}}}*/
void IntVecParam::DeepEcho(void){/*{{{*/

	int i;

	_printf_("IntVecParam:\n");
	_printf_("   enum: " << this->enum_type << " (" << EnumToStringx(this->enum_type) << ")\n");
	_printf_("   vector size: " << this->M << "\n");
	for(i=0;i<this->M;i++){
		_printf_(i << " " << this->values[i] << "\n");
	}
}
/*}}}*/
void IntVecParam::Echo(void){/*{{{*/

	_printf_("IntVecParam:\n");
	_printf_("   enum: " << this->enum_type << " (" << EnumToStringx(this->enum_type) << ")\n");
	_printf_("   vector size: " << this->M << "\n");

}
/*}}}*/
int  IntVecParam::Id(void){ return -1; }/*{{{*/
/*}}}*/
void IntVecParam::Marshall(char** pmarshalled_data,int* pmarshalled_data_size, int marshall_direction){ /*{{{*/

	MARSHALLING_ENUM(IntVecParamEnum);

	MARSHALLING(enum_type);
	MARSHALLING(M);
	if(M) { 
		MARSHALLING_DYNAMIC(values,int,M);
	}
	else values=NULL;

}
/*}}}*/
int  IntVecParam::ObjectEnum(void){/*{{{*/

	return IntVecParamEnum;

}
/*}}}*/

/*IntVecParam virtual functions definitions: */
void  IntVecParam::GetParameterValue(int** pintarray,int* pM){/*{{{*/
	int* output=NULL;

	if(M){
		output=xNew<int>(M);
		xMemCpy<int>(output,values,M);
	}

	/*Assign output pointers:*/
	if(pM) *pM=M;
	*pintarray=output;
}
/*}}}*/
void  IntVecParam::SetValue(int* intarray,int in_M){/*{{{*/

	/*avoid leak: */
	xDelete<int>(this->values);

	if(in_M){
		this->values=xNew<int>(in_M);
		xMemCpy<int>(this->values,intarray,in_M);
	}
	else this->values=NULL;

	this->M=in_M;
}
/*}}}*/
