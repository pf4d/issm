/*!\file SpcStatic.c
 * \brief: implementation of the SpcStatic object
 */

#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#include "../classes.h"
#include "./Constraint.h"
#include "../../shared/shared.h"

/*SpcStatic constructors and destructor*/
SpcStatic::SpcStatic(){/*{{{*/
	return;
}
/*}}}*/
SpcStatic::SpcStatic(int spc_sid,int spc_nodeid, int spc_dof,IssmDouble spc_value,int spc_analysis_type){/*{{{*/

	sid           = spc_sid;
	nodeid        = spc_nodeid;
	dof           = spc_dof;
	value         = spc_value;
	analysis_type = spc_analysis_type;
	penalty       = false;

	return;
}
/*}}}*/
SpcStatic::~SpcStatic(){/*{{{*/
	return;
}
/*}}}*/

/*Object virtual functions definitions:*/
Object* SpcStatic::copy() {/*{{{*/
	
	SpcStatic* spcstat = new SpcStatic(*this); 

	spcstat->sid=this->sid;
	spcstat->nodeid=this->nodeid;
	spcstat->dof=this->dof;
	spcstat->value=this->value;
	spcstat->analysis_type=this->analysis_type;

	return (Object*) spcstat;
}
/*}}}*/
void    SpcStatic::DeepEcho(void){/*{{{*/

	_printf_("SpcStatic:\n");
	_printf_("   sid: " << sid << "\n");
	_printf_("   nodeid: " << nodeid << "\n");
	_printf_("   dof: " << dof << "\n");
	_printf_("   value: " << value << "\n");
	_printf_("   analysis_type: " << EnumToStringx(analysis_type) << "\n");
	return;
}		
/*}}}*/
void    SpcStatic::Echo(void){/*{{{*/

	_printf_("SpcStatic:\n");
	_printf_("   sid: " << sid << "\n");
	_printf_("   nodeid: " << nodeid << "\n");
	_printf_("   dof: " << dof << "\n");
	_printf_("   value: " << value << "\n");
	_printf_("   analysis_type: " << EnumToStringx(analysis_type) << "\n");
	return;
}
/*}}}*/
int     SpcStatic::Id(void){ return sid; }/*{{{*/
/*}}}*/
void    SpcStatic::Marshall(char** pmarshalled_data,int* pmarshalled_data_size, int marshall_direction){ /*{{{*/

	MARSHALLING_ENUM(SpcStaticEnum);

	MARSHALLING(sid);
	MARSHALLING(nodeid);
	MARSHALLING(dof);
	MARSHALLING(value);
	MARSHALLING(analysis_type);
	MARSHALLING(penalty);

}
/*}}}*/
int     SpcStatic::ObjectEnum(void){/*{{{*/

	return SpcStaticEnum;

}
/*}}}*/

/*Constraint virtual functions definitions: */
void SpcStatic::ActivatePenaltyMethod(void){/*{{{*/
	   this->penalty = true;
}
/*}}}*/
void SpcStatic::ConstrainNode(Nodes* nodes,Parameters* parameters){/*{{{*/

	Node* node=NULL;

	/*Chase through nodes and find the node to which this SpcStatic applys: */
	node=(Node*)nodes->GetObjectById(NULL,nodeid);

	/*Apply constraint: */
	if(node){ //in case the spc is dealing with a node on another cpu
		node->ApplyConstraint(dof,value);
	}
}
/*}}}*/
bool SpcStatic::InAnalysis(int in_analysis_type){/*{{{*/
	if (in_analysis_type==this->analysis_type) return true;
	else return false;
}
/*}}}*/

/*SpcStatic functions*/
int        SpcStatic::GetDof(){/*{{{*/
	return dof;
}
/*}}}*/
int        SpcStatic::GetNodeId(){/*{{{*/

	return nodeid;
}
/*}}}*/
IssmDouble SpcStatic::GetValue(){/*{{{*/
	_assert_(!xIsNan<IssmDouble>(value));
	return value;
}
/*}}}*/
