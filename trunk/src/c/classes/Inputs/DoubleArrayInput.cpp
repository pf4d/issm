/*!\file DoubleArrayInput.c
 * \brief: implementation of the DoubleArrayInput object
 */

#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#include "../classes.h"
#include "../../shared/shared.h"

/*DoubleArrayInput constructors and destructor*/
DoubleArrayInput::DoubleArrayInput(){/*{{{*/
	return;
}
/*}}}*/
DoubleArrayInput::DoubleArrayInput(int in_enum_type,IssmDouble* in_values,  int in_m){/*{{{*/

	enum_type=in_enum_type;
	m=in_m;
	values=xNew<IssmDouble>(m);
	xMemCpy<IssmDouble>(values,in_values,m);

}
/*}}}*/
DoubleArrayInput::~DoubleArrayInput(){/*{{{*/

	if(values)xDelete<IssmDouble>(values);

	return;
}
/*}}}*/

/*Object virtual functions definitions:*/
Object* DoubleArrayInput::copy() {/*{{{*/

	return new DoubleArrayInput(this->enum_type,this->values,this->m);

}
/*}}}*/
void DoubleArrayInput::DeepEcho(void){/*{{{*/

	_printf_(setw(15)<<"   DoubleArrayInput "<<setw(25)<<left<<EnumToStringx(this->enum_type)<<" Size: " << m << "\n");
	for (int i=0;i<m;i++) _printf_(setw(20) << this->values[i]<<"\n");

}
/*}}}*/
void DoubleArrayInput::Echo(void){/*{{{*/
	this->DeepEcho();
}
/*}}}*/
int DoubleArrayInput::Id(void){ return -1; }/*{{{*/
/*}}}*/
void DoubleArrayInput::Marshall(char** pmarshalled_data,int* pmarshalled_data_size, int marshall_direction){ /*{{{*/

	MARSHALLING_ENUM(DoubleArrayInputEnum);

	MARSHALLING(enum_type);
	MARSHALLING(m);
	MARSHALLING_DYNAMIC(this->values,IssmDouble,m);
}
/*}}}*/
int DoubleArrayInput::ObjectEnum(void){/*{{{*/

	return DoubleArrayInputEnum;

}
/*}}}*/

/*DoubleArrayInput management*/
void DoubleArrayInput::GetValues(IssmDouble** pvalues, int *pm){ /*{{{*/

	/*output: */
	IssmDouble*  outvalues= NULL;

	outvalues=xNew<IssmDouble>(m);

	xMemCpy<IssmDouble>(outvalues,values,m);

	/*assign output pointers: */
	*pm=m;
	*pvalues=outvalues;
}
/*}}}*/
int DoubleArrayInput::InstanceEnum(void){/*{{{*/

	return this->enum_type;

}
/*}}}*/
void DoubleArrayInput::ResultToMatrix(IssmDouble* values,int ncols,int sid){/*{{{*/

	int ncols_local = this->GetResultArraySize();

	/*Some checks*/
	_assert_(values);
	_assert_(ncols_local<=ncols);

	/*Fill in arrays*/
	for(int i=0;i<ncols_local;i++) values[sid*ncols + i] = this->values[i];
}
/*}}}*/

/*Object functions*/
void DoubleArrayInput::ChangeEnum(int newenumtype){/*{{{*/
	this->enum_type=newenumtype;
}
/*}}}*/
void DoubleArrayInput::Configure(Parameters* parameters){/*{{{*/
	/*do nothing: */
}
/*}}}*/
