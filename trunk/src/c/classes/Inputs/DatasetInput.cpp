/*!\file DatasetInput.c
 * \brief: implementation of the datasetinput object
 */
/*Headers{{{*/
#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#include "../classes.h"
#include "../../shared/shared.h"
/*}}}*/

/*DatasetInput constructors and destructor*/
DatasetInput::DatasetInput(){/*{{{*/
	enum_type = UNDEF;
	inputs    = NULL;
	numids    = 0;
	ids       = NULL;
}
/*}}}*/
DatasetInput::DatasetInput(int in_enum_type){/*{{{*/

	enum_type = in_enum_type;
	numids    = 0;
	inputs    = new Inputs();
	ids       = NULL;
}
/*}}}*/
DatasetInput::~DatasetInput(){/*{{{*/
	xDelete<int>(this->ids);
	delete inputs;
}
/*}}}*/

/*Object virtual functions definitions:*/
Object* DatasetInput::copy() {/*{{{*/

	DatasetInput* output=NULL;

	output = new DatasetInput();
	output->enum_type=this->enum_type;
	output->numids=this->numids;
	output->ids=xNew<int>(output->numids);
	xMemCpy(output->ids,this->ids,output->numids);
	output->inputs=static_cast<Inputs*>(this->inputs->Copy());

	return (Object*)output;
}
/*}}}*/
void DatasetInput::DeepEcho(void){/*{{{*/

	_printf_("DatasetInput:\n");
	_printf_("   enum: " << this->enum_type << " (" << EnumToStringx(this->enum_type) << ")\n");
	_printf_("   numids:"<< this->numids<< "\n");
	_printf_("      ids: ");
	for(int i=0;i<this->numids;i++) _printf_(this->ids[i]<<" ");
	_printf_("\n");
	_printf_("   inputs: \n"); inputs->Echo();
}
/*}}}*/
void DatasetInput::Echo(void){/*{{{*/
	this->DeepEcho();
}
/*}}}*/
int    DatasetInput::Id(void){ return -1; }/*{{{*/
/*}}}*/
void DatasetInput::Marshall(char** pmarshalled_data,int* pmarshalled_data_size, int marshall_direction){ /*{{{*/

	MARSHALLING_ENUM(DatasetInputEnum);

	MARSHALLING(enum_type);
	MARSHALLING(numids);
	MARSHALLING_DYNAMIC(ids,int,numids)
	if (marshall_direction == MARSHALLING_BACKWARD) inputs = new Inputs();
	inputs->Marshall(pmarshalled_data,pmarshalled_data_size,marshall_direction);

}
/*}}}*/
int DatasetInput::ObjectEnum(void){/*{{{*/

	return DatasetInputEnum;

}
/*}}}*/
Input* DatasetInput::SpawnSegInput(int index1,int index2){/*{{{*/

	/*output*/
	DatasetInput* outinput=NULL;

	/*Create new Datasetinput (copy of current input)*/
	outinput=new DatasetInput();
	outinput->enum_type=this->enum_type;
	outinput->inputs=xDynamicCast<Inputs*>(this->inputs->SpawnSegInputs(index1,index2));
	outinput->numids=this->numids;
	outinput->ids=xNew<int>(this->numids);
	xMemCpy(outinput->ids,this->ids,this->numids);

	/*Assign output*/
	return outinput;
}
/*}}}*/
Input* DatasetInput::SpawnTriaInput(int index1,int index2,int index3){/*{{{*/

	/*output*/
	DatasetInput* outinput=NULL;

	/*Create new Datasetinput (copy of current input)*/
	outinput=new DatasetInput();
	outinput->enum_type=this->enum_type;
	outinput->inputs=xDynamicCast<Inputs*>(this->inputs->SpawnTriaInputs(index1,index2,index3));
	outinput->numids=this->numids;
	outinput->ids=xNew<int>(this->numids);
	xMemCpy(outinput->ids,this->ids,this->numids);

	/*Assign output*/
	return outinput;
}
/*}}}*/

/*DatasetInput management*/
void DatasetInput::AddInput(Input* input,int id){/*{{{*/

	_assert_(this->numids == this->inputs->Size());

	int *old_ids = NULL;

	if(this->numids>0){
		old_ids=xNew<int>(this->numids);
		xMemCpy(old_ids,this->ids,this->numids);
		xDelete<int>(this->ids);
	}

	this->numids=this->numids+1;
	this->ids=xNew<int>(this->numids);

	if(this->numids>1){
		xMemCpy(this->ids,old_ids,this->numids-1);
		xDelete<int>(old_ids);
	}

	/*go ahead and plug: */
	this->ids[this->numids-1]=id;
	inputs->AddObject(input);

	_assert_(this->numids == this->inputs->Size());
}
/*}}}*/
int DatasetInput::InstanceEnum(void){/*{{{*/

	return this->enum_type;

}
/*}}}*/

/*Object functions*/
void DatasetInput::Configure(Parameters* parameters){/*{{{*/
	/*do nothing: */
}
/*}}}*/
void DatasetInput::GetInputValue(IssmDouble* pvalue,Gauss* gauss,int id){/*{{{*/

	int  offset = -1;
	_assert_(this->numids == this->inputs->Size());

	/*Get requested input within dataset*/
	for(int i=0;i<this->numids;i++) if(this->ids[i]==id) offset=i;
	if(offset<0) _error_("Could not find input of id "<<id );

	Input* input=xDynamicCast<Input*>(this->inputs->GetObjectByOffset(offset));
	input->GetInputValue(pvalue,gauss);
}
/*}}}*/
