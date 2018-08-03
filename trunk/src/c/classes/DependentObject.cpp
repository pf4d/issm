/*!\file DependentObject.c
 * \brief: implementation of the DependentObject object
 */

#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#include "./classes.h"
#include "shared/shared.h"
#include "../modules/modules.h"

/*DependentObject constructors and destructor*/
DependentObject::DependentObject(){/*{{{*/
	this->name=NULL;
	this->type=0;
	this->index=-1;
}
/*}}}*/
DependentObject::DependentObject(char* in_name, int in_type,int in_index){/*{{{*/

	this->name=xNew<char>(strlen(in_name)+1); xMemCpy<char>(this->name,in_name,strlen(in_name)+1);
	this->type=in_type;
	this->index=in_index;
	if(in_type!=0 && in_type!=1)_error_("cannot create an DependentObject of type " << in_type);
	if(in_type==1)_error_("not implemented yet!");

}
/*}}}*/
DependentObject::~DependentObject(){ //destructor/*{{{*/
	xDelete<char>(this->name);
}
/*}}}*/

/*Object virtual functions definitions:*/
Object* DependentObject::copy(void) { /*{{{*/
	return new DependentObject(name,type,index);
} /*}}}*/
void DependentObject::DeepEcho(void){/*{{{*/
	this->Echo();
}
/*}}}*/
void DependentObject::Echo(void){/*{{{*/

	_printf_("DependentObject:\n");
	_printf_("   name: " << this->name << "\n");
	if(this->type==0)
		_printf_("   type: scalar\n");
	else if(this->type==1)
		_printf_("   type: vertex\n");
	else
		_error_(" unknown type: " << this->type);
	if(this->index>=0) _printf_("   index: " << this->index << "\n");
}
/*}}}*/
int    DependentObject::Id(void){ return -1; }/*{{{*/
/*}}}*/
int DependentObject::ObjectEnum(void){/*{{{*/

	return DependentObjectEnum;

}
/*}}}*/

/*DependentObject methods: */
int  DependentObject::NumDependents(void){/*{{{*/

	/*Branch according to the type of variable: */
	if(type==0){ /*scalar:*/
		return 1;
	}
	else if(type==1){ /* vector:*/
		_error_("not implemented yet!");
	}
	else _error_("should not have a type of " << type);
}
/*}}}*/
void  DependentObject::Responsex(IssmDouble* poutput_value,FemModel* femmodel){/*{{{*/

	/*Is this some special type of response for which we need to go in the output definitions? :*/
	if (StringToEnumx(this->name,false)==-1){
		*poutput_value=OutputDefinitionsResponsex(femmodel,this->name);
	}
	else femmodel->Responsex(poutput_value,this->name);

}
/*}}}*/
