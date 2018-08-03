/*!\file Nodalvalue.h
 * \brief: header file for Nodalvalue object
 */

#ifndef _NODALVALUE_H_
#define _NODALVALUE_H_

/*Headers:*/
/*{{{*/
#include "./Definition.h"
#include "../datastructures/datastructures.h"
#include "./Elements/Element.h"
#include "./Elements/Elements.h"
#include "./FemModel.h"
#include "../modules/SurfaceAreax/SurfaceAreax.h"
#include "../classes/Params/Parameters.h"
#include "../classes/Inputs/Input.h"
#include "../classes/gauss/Gauss.h"
/*}}}*/

void NodalValuex( IssmDouble* pnodalvalue, int natureofdataenum,Elements* elements,Nodes* nodes, Vertices* vertices, Loads* loads, Materials* materials, Parameters* parameters);
IssmDouble OutputDefinitionsResponsex(FemModel* femmodel,int output_enum);

class Nodalvalue: public Object, public Definition{

	public: 

		int         definitionenum;
		int         model_enum;
		char*       name;
		int         node;
		
		/*Nodalvalue constructors, destructors :*/
		Nodalvalue(){/*{{{*/

			this->definitionenum = -1;
			this->name = NULL;
			this->model_enum = UNDEF;
			this->node = -1;

		}
		/*}}}*/
		Nodalvalue(char* in_name, int in_definitionenum, int in_model_enum, int in_node){/*{{{*/

			this->definitionenum=in_definitionenum;
			this->name   = xNew<char>(strlen(in_name)+1);
			xMemCpy<char>(this->name,in_name,strlen(in_name)+1);

			this->model_enum=in_model_enum;
			this->node=in_node;
		}
		/*}}}*/
		~Nodalvalue(){/*{{{*/
			if(this->name)xDelete(this->name);
		}
		/*}}}*/
		/*Object virtual function resolutoin: */
		Object* copy() {/*{{{*/
			Nodalvalue* mf = new Nodalvalue(this->name,this->definitionenum, this->model_enum,this->node);
			return (Object*) mf;
		}
		/*}}}*/
		void DeepEcho(void){/*{{{*/
			this->Echo();
		}
		/*}}}*/
		void Echo(void){/*{{{*/
			_printf_(" Nodalvalue: " << name << " " << this->definitionenum << "\n");
			_printf_("    model_enum: " << model_enum << " " << EnumToStringx(model_enum) << "\n");
			_printf_("    node: " << node << "\n");
		}
		/*}}}*/
		int Id(void){/*{{{*/
			return -1;
		}
		/*}}}*/
		void Marshall(char** pmarshalled_data,int* pmarshalled_data_size, int marshall_direction){/*{{{*/
			_error_("not implemented yet!"); 
		} 
		/*}}}*/
		int ObjectEnum(void){/*{{{*/
			return NodalvalueEnum;
		}
		/*}}}*/
		/*Definition virtual function resolutoin: */
		int DefinitionEnum(){/*{{{*/

			return this->definitionenum;
		}
		/*}}}*/
		char* Name(){/*{{{*/

			char* name2=xNew<char>(strlen(this->name)+1);
			xMemCpy(name2,this->name,strlen(this->name)+1);

			return name2;
		}
		/*}}}*/
		 IssmDouble Response(FemModel* femmodel){/*{{{*/
			
			 /*output:*/
			 IssmDouble value;

			 /*set index, which will be used by the NodalValue module: */
			 femmodel->parameters->SetParam(node,IndexEnum);

			 /*call Nodalvalue:*/
			 NodalValuex(&value, model_enum, femmodel->elements, femmodel->nodes, femmodel->vertices, femmodel->loads, 
					 femmodel->materials, femmodel->parameters);

			 /*done:*/
			 return value;
		 }
		 /*}}}*/
};

#endif  /* _NODALVALUE_H_ */
