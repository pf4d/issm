/*!\file SpcDynamic.c
 * \brief: implementation of the SpcDynamic object
 */

#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#include "../classes.h"
#include "./Constraint.h"
#include "../../shared/shared.h"

/*SpcDynamic constructors and destructor*/
SpcDynamic::SpcDynamic(){/*{{{*/
	return;
}
/*}}}*/
SpcDynamic::SpcDynamic(int spc_sid,int spc_nodeid, int spc_dof,int spc_analysis_type){/*{{{*/

	sid           = spc_sid;
	nodeid        = spc_nodeid;
	dof           = spc_dof;
	value         = 0;
	analysis_type = spc_analysis_type;
	isset         = false;
	penalty       = false;

	return;
}
/*}}}*/
SpcDynamic::~SpcDynamic(){/*{{{*/
	return;
}
/*}}}*/

/*Object virtual functions definitions:*/
Object* SpcDynamic::copy() {/*{{{*/

	SpcDynamic* spcdyn = new SpcDynamic(*this); 

	spcdyn->sid=this->sid;
	spcdyn->nodeid=this->nodeid;
	spcdyn->dof=this->dof;
	spcdyn->value=this->value;
	spcdyn->analysis_type=this->analysis_type;
	spcdyn->isset=this->isset;

	return (Object*) spcdyn;
}
/*}}}*/
void    SpcDynamic::DeepEcho(void){/*{{{*/

	this->Echo();
	return;
}		
/*}}}*/
void    SpcDynamic::Echo(void){/*{{{*/

	_printf_("SpcDynamic:\n");
	_printf_("   sid: " << sid << "\n");
	_printf_("   nodeid: " << nodeid << "\n");
	_printf_("   dof: " << dof << "\n");
	_printf_("   value: " << value << "\n");
	_printf_("   isset: " <<(isset?"true":"false") << "\n");
	_printf_("   analysis_type: " << EnumToStringx(analysis_type) << "\n");
	return;
}
/*}}}*/
int     SpcDynamic::Id(void){ return sid; }/*{{{*/
/*}}}*/
void    SpcDynamic::Marshall(char** pmarshalled_data,int* pmarshalled_data_size, int marshall_direction){ /*{{{*/

	MARSHALLING_ENUM(SpcDynamicEnum);

	MARSHALLING(sid);
	MARSHALLING(nodeid);
	MARSHALLING(dof);
	MARSHALLING(value);
	MARSHALLING(analysis_type);
	MARSHALLING(isset);
	MARSHALLING(penalty);

}
/*}}}*/
int     SpcDynamic::ObjectEnum(void){/*{{{*/

	return SpcDynamicEnum;

}
/*}}}*/

/*Constraint virtual functions definitions: */
void SpcDynamic::ActivatePenaltyMethod(void){/*{{{*/
	this->penalty = true;
}
/*}}}*/
void SpcDynamic::ConstrainNode(Nodes* nodes,Parameters* parameters){/*{{{*/

	Node* node=NULL;

	/*Chase through nodes and find the node to which this SpcDynamic applys: */
	node=(Node*)nodes->GetObjectById(NULL,nodeid);

	/*Apply constraint: */
	if(node){ //in case the spc is dealing with a node on another cpu

		/*We should first check that the value has been set... (test306)*/
		node->ApplyConstraint(dof,value);
	}
}
/*}}}*/
bool SpcDynamic::InAnalysis(int in_analysis_type){/*{{{*/
	if (in_analysis_type==this->analysis_type) return true;
	else return false;
}
/*}}}*/

/*SpcDynamic functions*/
int        SpcDynamic::GetDof(){/*{{{*/
	return dof;
}
/*}}}*/
int        SpcDynamic::GetNodeId(){/*{{{*/

	return nodeid;
}
/*}}}*/
IssmDouble SpcDynamic::GetValue(){/*{{{*/
	_assert_(this->isset);
	_assert_(!xIsNan<IssmDouble>(value));
	return value;
}
/*}}}*/
void       SpcDynamic::SetDynamicConstraint(Nodes* nodes,IssmDouble* yg_serial){/*{{{*/

	int pos;

	Node* node=(Node*)nodes->GetObjectById(NULL,nodeid);
	pos=node->GetDof(dof,GsetEnum);

	this->value=yg_serial[pos];
	this->isset=true;
}
/*}}}*/
