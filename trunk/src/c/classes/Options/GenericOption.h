/*! \file GenericOption.h 
 *  \brief: header file for generic option object
 */

#ifndef _GENERIC_OPTION_
#define _GENERIC_OPTION_

/*Headers:{{{*/
#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#include <cstring>
#include "../../shared/shared.h"
#include "../../datastructures/datastructures.h"
#include "./OptionUtilities.h"
/*}}}*/

template <class OptionType> 
class GenericOption: public Option {

	public:

		char       *name;
		OptionType  value;

		int         numel;   //in case OptionType is an array
		int         ndims;   //in case OptionType is a multi-dimensional array: */
		int        *size;

		/*GenericOption constructors, destructors*/
		GenericOption(){ /*{{{*/

			name   = NULL;
			numel  = 0;
			ndims  = 0;
			size   = NULL;

		} /*}}}*/
		~GenericOption(){ /*{{{*/

			if(name)   xDelete<char>(name);
			if(size)   xDelete<int>(size);

		} /*}}}*/

		/*Object virtual functions definitions:*/
		Object* copy(){/*{{{*/
			_error_("Not implemented yet");
		};/*}}}*/
		void DeepEcho(){ /*{{{*/

			char  indent[81]="";
			this->DeepEcho(indent);

		} /*}}}*/
		void DeepEcho(char* indent){ /*{{{*/

			char  cstr[81];
			bool  flag=true;

			if(flag) _printf0_(indent << "         name: \"" << name << "\"\n");
			if(flag) _printf0_(indent << "         numel: " << numel << "\n");
			if(flag) _printf0_(indent << "         ndims: " << ndims << "\n");
			if(size){
				StringFromSize(cstr,size,ndims);
				if(flag) _printf0_(indent << "          size: " << cstr << "\n");
			}
			else if(flag) _printf0_(indent << "          size: [empty]\n");
			_printf_(indent << "         value: " << value << "\n");;
		} /*}}}*/
		void Echo(){ /*{{{*/

			this->DeepEcho();

		} /*}}}*/
		int  Id(){/*{{{*/
			_error_("Not implemented yet");
		};/*}}}*/
		int  ObjectEnum(){/*{{{*/
			return GenericOptionEnum;
		};/*}}}*/

		/*GenericOption functions: */
		void  Get(OptionType* pvalue){/*{{{*/
			*pvalue=value; 
		};/*}}}*/
		char* Name(){/*{{{*/
			return name;
		};/*}}}*/
		int   NDims(){/*{{{*/
			return ndims;
		};/*}}}*/
		int   NumEl(){/*{{{*/
			return numel;
		};/*}}}*/
		int*  Size(){/*{{{*/
			return size;
		};/*}}}*/
};

#if defined(_HAVE_ADOLC_) && !defined(_WRAPPERS_)  //We hook off this specific specialization when not running ADOLC, otherwise we get a redeclaration with the next specialization. 
template <> inline void GenericOption<IssmPDouble*>::Get(IssmPDouble** pvalue){ /*{{{*/

	/*Copy vector*/
	IssmPDouble* outvalue=xNew<IssmPDouble>(this->NumEl());
	for(int i=0;i<this->NumEl();i++) outvalue[i]=this->value[i];

	/*Assign output pointer*/
	*pvalue=outvalue;
} /*}}}*/
#endif
template <> inline void GenericOption<IssmDouble*>::Get(IssmDouble** pvalue){ /*{{{*/

	/*Copy vector*/
	IssmDouble* outvalue=xNew<IssmDouble>(this->NumEl());
	for(int i=0;i<this->NumEl();i++) outvalue[i]=this->value[i];

	/*Assign output pointer*/
	*pvalue=outvalue;
} /*}}}*/
template <> inline void GenericOption<char*>::Get(char** pvalue){ /*{{{*/

	int   stringsize=strlen(this->value)+1;
	char* outstring=xNew<char>(stringsize);
	xMemCpy<char>(outstring,this->value,stringsize);

	*pvalue=outstring;
} 
/*}}}*/

/*Special destructors when there is a pointer*/
template <> inline GenericOption<char*>::~GenericOption(){ /*{{{*/

	if(name)   xDelete<char>(name);
	if(size)   xDelete<int>(size);
	if(value)  xDelete<char>(value);
} 
/*}}}*/

#endif  /* _OPTIONOBJECT_H */
