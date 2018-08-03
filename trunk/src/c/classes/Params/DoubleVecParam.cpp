/*!\file DoubleVecParam.c
 * \brief: implementation of the DoubleVecParam object
 */

/*header files: */
/*{{{*/
#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#include "../classes.h"
#include "shared/shared.h"
/*}}}*/

/*DoubleVecParam constructors and destructor*/
DoubleVecParam::DoubleVecParam(){/*{{{*/
	return;
}
/*}}}*/
DoubleVecParam::DoubleVecParam(int in_enum_type,IssmDouble* in_values, int in_M){/*{{{*/

	enum_type=in_enum_type;
	M=in_M;

	if(M){
		values=xNew<IssmDouble>(M);
		xMemCpy<IssmDouble>(values,in_values,M);
	}
	else values=NULL;
}
/*}}}*/
DoubleVecParam::~DoubleVecParam(){/*{{{*/
	xDelete<IssmDouble>(values);
	return;
}
/*}}}*/

/*Object virtual functions definitions:*/
Param* DoubleVecParam::copy() {/*{{{*/

	return new DoubleVecParam(this->enum_type,this->values,this->M);

}
/*}}}*/
void DoubleVecParam::DeepEcho(void){/*{{{*/
	_printf_(setw(22)<<"   DoubleVecParam "<<setw(35)<<left<<EnumToStringx(this->enum_type)<<" "<<"[");
	for(int i=0;i<this->M;i++) _printf_(" "<< this->values[i]);
	_printf_("\n");
}
/*}}}*/
void DoubleVecParam::Echo(void){/*{{{*/

	_printf_(setw(22)<<"   DoubleVecParam "<<setw(35)<<left<<EnumToStringx(this->enum_type)<<" size: "<<this->M<<"\n");

}
/*}}}*/
int    DoubleVecParam::Id(void){ return -1; }/*{{{*/
/*}}}*/
void DoubleVecParam::Marshall(char** pmarshalled_data,int* pmarshalled_data_size, int marshall_direction){ /*{{{*/

	MARSHALLING_ENUM(DoubleVecParamEnum);

	MARSHALLING(enum_type);
	MARSHALLING(M);
	MARSHALLING_DYNAMIC(values,IssmDouble,M);

}
/*}}}*/
int DoubleVecParam::ObjectEnum(void){/*{{{*/

	return DoubleVecParamEnum;

}
/*}}}*/

/*DoubleVecParam virtual functions definitions: */
void  DoubleVecParam::GetParameterValue(IssmDouble** pIssmDoublearray,int* pM){/*{{{*/
	IssmDouble* output=NULL;
	int M;

	M=this->M;
	output=xNew<IssmDouble>(M);
	xMemCpy<IssmDouble>(output,values,M);

	/*Assign output pointers:*/
	if(pM) *pM=M;
	*pIssmDoublearray=output;
}
/*}}}*/
void  DoubleVecParam::GetParameterValue(IssmDouble** pIssmDoublearray,int* pM,int* pN){/*{{{*/
	IssmDouble* output=NULL;
	int M;
	int N;

	N=1;
	M=this->M;
	output=xNew<IssmDouble>(M);
	xMemCpy<IssmDouble>(output,values,M);

	/*Assign output pointers:*/
	if(pM) *pM=M;
	if(pN) *pN=N;
	*pIssmDoublearray=output;
}
/*}}}*/
void  DoubleVecParam::GetParameterValue(int** pintarray,int* pM){/*{{{*/
	_error_("DoubleVec param of enum " << enum_type << " (" << EnumToStringx(enum_type) << ") cannot return an array of int");
}
/*}}}*/
void  DoubleVecParam::SetValue(IssmDouble* IssmDoublearray,int in_M){/*{{{*/

	/*avoid leak: */
	xDelete<IssmDouble>(this->values);

	this->values=xNew<IssmDouble>(in_M);
	xMemCpy<IssmDouble>(this->values,IssmDoublearray,in_M);

	this->M=in_M;
}
/*}}}*/

/*DoubleVecParam specific routines:*/
void  DoubleVecParam::GetParameterValueByPointer(IssmDouble** pIssmDoublearray,int* pM){/*{{{*/
	
	/*Assign output pointers:*/
	if(pM) *pM=M;
	*pIssmDoublearray=values;
}
/*}}}*/
