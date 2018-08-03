/*!\file TransientInput.c
 * \brief: implementation of the TransientInput object
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

/*TransientInput constructors and destructor*/
TransientInput::TransientInput(){/*{{{*/

	enum_type=UNDEF;
	inputs=NULL;
	this->numtimesteps=0;
	this->parameters=NULL;
	this->timesteps=NULL;

}
/*}}}*/
TransientInput::TransientInput(int in_enum_type)/*{{{*/
{
	/*Set Enum*/
	enum_type=in_enum_type;

	/*Allocate values and timesteps, and copy: */
	this->numtimesteps=0;
	this->timesteps=NULL;
	inputs = new Inputs();
	this->parameters=NULL;

}
/*}}}*/
TransientInput::TransientInput(int in_enum_type,IssmDouble* timesin,int N){/*{{{*/

	/*Set Enum*/
	enum_type=in_enum_type;

	/*Allocate values and timesteps, and copy: */
	this->numtimesteps=N;
	this->timesteps=xNew<IssmDouble>(N);
	xMemCpy(this->timesteps,timesin,N);

	inputs = new Inputs();
	this->parameters=NULL;
}
/*}}}*/
TransientInput::~TransientInput(){/*{{{*/
	xDelete(this->timesteps);
	this->timesteps=NULL;
	this->numtimesteps=0;
	parameters=NULL;
	delete this->inputs;
	return;
}
/*}}}*/

/*Object virtual functions definitions:*/
Object* TransientInput::copy() {/*{{{*/

	TransientInput* output=NULL;

	output = new TransientInput();
	output->enum_type=this->enum_type;
	output->numtimesteps=this->numtimesteps;
	output->timesteps=xNew<IssmDouble>(this->numtimesteps);
	xMemCpy(output->timesteps,this->timesteps,this->numtimesteps);
	output->inputs=static_cast<Inputs*>(this->inputs->Copy());
	output->parameters=this->parameters;

	return (Object*)output;

}
/*}}}*/
void TransientInput::DeepEcho(void){/*{{{*/

	int i;

	_printf_("TransientInput:\n");
	_printf_("   enum: " << this->enum_type << " (" << EnumToStringx(this->enum_type) << ")\n");
	_printf_("   numtimesteps: " << this->numtimesteps << "\n");
	_printf_("---inputs: \n"); 
	for(i=0;i<this->numtimesteps;i++){
		_printf_("   time: " << this->timesteps[i]<<"  ");
		((Input*)this->inputs->GetObjectByOffset(i))->Echo();
	}
}
/*}}}*/
void TransientInput::Echo(void){/*{{{*/
	this->DeepEcho();
}
/*}}}*/
int  TransientInput::Id(void){ return -1; }/*{{{*/
/*}}}*/
void TransientInput::Marshall(char** pmarshalled_data,int* pmarshalled_data_size, int marshall_direction){ /*{{{*/

	if (marshall_direction == MARSHALLING_BACKWARD) inputs = new Inputs();

	MARSHALLING_ENUM(TransientInputEnum);

	MARSHALLING(enum_type);
	MARSHALLING(numtimesteps);
	MARSHALLING_DYNAMIC(this->timesteps,IssmDouble,numtimesteps);
	inputs->Marshall(pmarshalled_data,pmarshalled_data_size,marshall_direction);
}
/*}}}*/
int  TransientInput::ObjectEnum(void){/*{{{*/

	return TransientInputEnum;

}
/*}}}*/

/*TransientInput management*/
int TransientInput::InstanceEnum(void){/*{{{*/

	return this->enum_type;

}
/*}}}*/
Input* TransientInput::SpawnSegInput(int index1,int index2){/*{{{*/

	/*output*/
	TransientInput* outinput=NULL;

	/*Create new Transientinput (copy of current input)*/
	outinput=new TransientInput();
	outinput->enum_type=this->enum_type;
	outinput->numtimesteps=this->numtimesteps;
	outinput->timesteps=xNew<IssmDouble>(this->numtimesteps);
	xMemCpy(outinput->timesteps,this->timesteps,this->numtimesteps);
	outinput->inputs=(Inputs*)this->inputs->SpawnSegInputs(index1,index2);
	outinput->parameters=this->parameters;

	/*Assign output*/
	return outinput;

}
/*}}}*/
Input* TransientInput::SpawnTriaInput(int index1,int index2,int index3){/*{{{*/

	/*output*/
	TransientInput* outinput=NULL;

	/*Create new Transientinput (copy of current input)*/
	outinput=new TransientInput();
	outinput->enum_type=this->enum_type;
	outinput->numtimesteps=this->numtimesteps;
	outinput->timesteps=xNew<IssmDouble>(this->numtimesteps);
	xMemCpy(outinput->timesteps,this->timesteps,this->numtimesteps);
	outinput->inputs=(Inputs*)this->inputs->SpawnTriaInputs(index1,index2,index3);
	outinput->parameters=this->parameters;

	/*Assign output*/
	return outinput;

}
/*}}}*/

/*Object functions*/
void TransientInput::ChangeEnum(int newenumtype){/*{{{*/
	this->enum_type=newenumtype;
}
/*}}}*/
void TransientInput::GetInputAllTimeAverages(IssmDouble** pvalues,IssmDouble** ptimes, int* pnumtimes){/*{{{*/

	int i;
	IssmDouble* times=NULL;
	IssmDouble* values=NULL;

	/*allocate: */
	times=xNew<IssmDouble>(this->numtimesteps);
	values=xNew<IssmDouble>(this->numtimesteps);

	for(i=0;i<numtimesteps;i++){
		Input* input=(Input*)this->inputs->GetObjectByOffset(i);
		input->GetInputAverage(values+i);
		times[i]=this->timesteps[i];
	}

	*pvalues=values;
	*ptimes=times;
	*pnumtimes=numtimesteps;
}
/*}}}*/
void TransientInput::GetInputAverage(IssmDouble* pvalue){/*{{{*/

	IssmDouble time;

	/*First, recover current time from parameters: */
	parameters->FindParam(&time,TimeEnum);

	/*Retrieve interpolated values for this time step: */
	Input* input=GetTimeInput(time);

	/*Call input function*/
	input->GetInputAverage(pvalue);

	delete input;

}
/*}}}*/
void TransientInput::GetInputDerivativeValue(IssmDouble* p, IssmDouble* xyz_list, Gauss* gauss){/*{{{*/

	IssmDouble time;

	/*First, recover current time from parameters: */
	parameters->FindParam(&time,TimeEnum);

	/*Retrieve interpolated values for this time step: */
	Input* input=GetTimeInput(time);

	/*Call input function*/
	input->GetInputDerivativeValue(p,xyz_list,gauss);

	delete input;
}
/*}}}*/
void TransientInput::GetInputValue(IssmDouble* pvalue,Gauss* gauss){/*{{{*/
	IssmDouble time;

	/*First, recover current time from parameters: */
	this->parameters->FindParam(&time,TimeEnum);

	/*Retrieve interpolated values for this time step: */
	Input* input=GetTimeInput(time);

	/*Call input function*/
	input->GetInputValue(pvalue,gauss);

	delete input;
}
/*}}}*/
void TransientInput::GetInputValue(IssmDouble* pvalue,Gauss* gauss,IssmDouble time){/*{{{*/

	/*Retrieve interpolated values for this time step: */
	Input* input=GetTimeInput(time);

	/*Call input function*/
	input->GetInputValue(pvalue,gauss);

	delete input;
}
/*}}}*/
void TransientInput::GetInputUpToCurrentTimeAverages(IssmDouble** pvalues, IssmDouble** ptimes, int* pnumtimes, IssmDouble currenttime){/*{{{*/

	int         i;
	IssmDouble *times  = NULL;
	IssmDouble *values = NULL;
	int         numsteps;
	bool        iscurrenttime_included = false;

	/*Figure out how many time steps we are going to return: */
	numsteps=0;
	for(i=0;i<numtimesteps;i++){
		if(this->timesteps[i]==currenttime)iscurrenttime_included=true;
		if (this->timesteps[i]>currenttime)break;
		else numsteps++;
	}
	if(iscurrenttime_included==false)numsteps++;

	/*allocate: */
	times=xNew<IssmDouble>(numsteps);
	values=xNew<IssmDouble>(numsteps);

	for(i=0;i<numsteps;i++){

		if((iscurrenttime_included==false) && (i==(numsteps-1))){

			/*Retrieve interpolated values for current time step: */
			Input* input=GetTimeInput(currenttime);
			input->GetInputAverage(values+i);
			times[i]=currenttime;
		}
		else{
			Input* input=(Input*)this->inputs->GetObjectByOffset(i);
			input->GetInputAverage(values+i);
			times[i]=this->timesteps[i];
		}
	}

	*pvalues=values;
	*ptimes=times;
	*pnumtimes=numtimesteps;
}
/*}}}*/

/*Intermediary*/
void TransientInput::AddTimeInput(Input* input,IssmDouble time){/*{{{*/

	/*insert values at time step: */
	if (this->numtimesteps>0 && time<=this->timesteps[this->numtimesteps-1]) _error_("timestep values must increase sequentially");

	//copy timesteps, add the new time, delete previous timesteps, and add the new input: inputs->AddObject(input);
	IssmDouble* old_timesteps=NULL;

	if (this->numtimesteps > 0){
		old_timesteps=xNew<IssmDouble>(this->numtimesteps);
		xMemCpy(old_timesteps,this->timesteps,this->numtimesteps);
		xDelete(this->timesteps);
	}

	this->numtimesteps=this->numtimesteps+1;
	this->timesteps=xNew<IssmDouble>(this->numtimesteps);

	if (this->numtimesteps > 1){
		xMemCpy(this->timesteps,old_timesteps,this->numtimesteps-1);
		xDelete(old_timesteps);
	}

	/*go ahead and plug: */
	this->timesteps[this->numtimesteps-1]=time;
	inputs->AddObject(input);

}
/*}}}*/
void TransientInput::AddTimeInput(Input* input){/*{{{*/

	_assert_(this->inputs->Size()<this->numtimesteps);
	inputs->AddObject(input);

}
/*}}}*/
void TransientInput::Configure(Parameters* parameters){/*{{{*/
	this->parameters=parameters;
}
/*}}}*/
void TransientInput::Extrude(int start){/*{{{*/

	for(int i=0;i<this->numtimesteps;i++){
		((Input*)this->inputs->GetObjectByOffset(i))->Extrude(start);
	}
}
/*}}}*/
int  TransientInput::GetResultArraySize(void){/*{{{*/

	return 1;
}
/*}}}*/
int  TransientInput::GetResultInterpolation(void){/*{{{*/

	IssmDouble time;
	int        output;

	parameters->FindParam(&time,TimeEnum);
	Input* input=GetTimeInput(time);
	output = input->GetResultInterpolation();

	/*Clean up and return*/
	delete input;
	return output;

}
/*}}}*/
int  TransientInput::GetResultNumberOfNodes(void){/*{{{*/

	IssmDouble time;
	int        output;

	parameters->FindParam(&time,TimeEnum);
	Input* input=GetTimeInput(time);
	output = input->GetResultNumberOfNodes();

	/*Clean up and return*/
	delete input;
	return output;

}
/*}}}*/
Input* TransientInput::GetTimeInput(IssmDouble intime){/*{{{*/

	IssmDouble deltat;
	IssmDouble alpha1,alpha2;
	int        found;
	int        offset;
	bool       interp;

	/*First, recover interp bool: */
	this->parameters->FindParam(&interp,TimesteppingInterpForcingsEnum);

	Input *input  = NULL;
	Input *input1 = NULL;
	Input *input2 = NULL;

	/*go through the timesteps, and figure out which interval we 
	 *fall within. Then interpolate the values on this interval: */
	found=binary_search(&offset,intime,this->timesteps,this->numtimesteps);
	if(!found) _error_("Input not found (is TransientInput sorted ?)");

	if (offset==-1){
		/*get values for the first time: */
		_assert_(intime<this->timesteps[0]);
		input=(Input*)((Input*)this->inputs->GetObjectByOffset(0))->copy();
	}
	else if(offset==(this->numtimesteps-1) || !interp){
		/*get values for the last time: */
		_assert_(intime>=this->timesteps[offset]);
		input=(Input*)((Input*)this->inputs->GetObjectByOffset(offset))->copy();
	}
	else {
		/*get values between two times [offset:offset+1[, Interpolate linearly*/
		_assert_(intime>=this->timesteps[offset] && intime<this->timesteps[offset+1]);
		deltat=this->timesteps[offset+1]-this->timesteps[offset];
		alpha2=(intime-this->timesteps[offset])/deltat;
		alpha1=(1.0-alpha2);

		input1=(Input*)this->inputs->GetObjectByOffset(offset); 
		input2=(Input*)this->inputs->GetObjectByOffset(offset+1);

		input=(Input*)input1->copy();
		input->Scale(alpha1);
		input->AXPY(input2,alpha2);
	}

	/*Assign output pointer*/
	return input;
}
/*}}}*/
void TransientInput::GetVectorFromInputs(Vector<IssmDouble>* vector,int* doflist){/*{{{*/

	IssmDouble time;

	/*First, recover current time from parameters: */
	parameters->FindParam(&time,TimeEnum);

	/*Retrieve interpolated values for this time step: */
	Input* input=GetTimeInput(time);

	/*Call input function*/
	input->GetVectorFromInputs(vector,doflist);

	delete input;

} /*}}}*/
IssmDouble TransientInput::InfinityNorm(void){/*{{{*/

	IssmDouble time;
	IssmDouble infnorm;

	/*First, recover current time from parameters: */
	parameters->FindParam(&time,TimeEnum);

   /*Retrieve interpolated values for this time step: */
	Input* input=GetTimeInput(time);

	/*Call input function*/
	infnorm=input->InfinityNorm();

	/*Clean-up and return*/
	delete input;
	return infnorm;
}
/*}}}*/
IssmDouble TransientInput::Max(void){/*{{{*/

	IssmDouble time;
	IssmDouble max;

	/*First, recover current time from parameters: */
	parameters->FindParam(&time,TimeEnum);

   /*Retrieve interpolated values for this time step: */
	Input* input=GetTimeInput(time);

	/*Call input function*/
	max=input->Max();

	delete input;

	return max;
}
/*}}}*/
IssmDouble TransientInput::MaxAbs(void){/*{{{*/

	IssmDouble time;
	IssmDouble maxabs;

	/*First, recover current time from parameters: */
	parameters->FindParam(&time,TimeEnum);

	/*Retrieve interpolated values for this time step: */
	Input* input=GetTimeInput(time);

	/*Call input function*/
	maxabs=input->MaxAbs();

	/*Clean-up and return*/
	delete input;
	return maxabs;

}
/*}}}*/
IssmDouble TransientInput::Min(void){/*{{{*/

	IssmDouble time;
	IssmDouble min;

	/*First, recover current time from parameters: */
	parameters->FindParam(&time,TimeEnum);

   /*Retrieve interpolated values for this time step: */
	Input* input=GetTimeInput(time);

	/*Call input function*/
	min=input->Min();

	/*Clean-up and return*/
	delete input;
	return min;

}
/*}}}*/
IssmDouble TransientInput::MinAbs(void){/*{{{*/

	IssmDouble time;
	IssmDouble minabs;

	/*First, recover current time from parameters: */
	parameters->FindParam(&time,TimeEnum);

	/*Retrieve interpolated values for this time step: */
	Input* input=GetTimeInput(time);

	/*Call input function*/
	minabs=input->MinAbs();

	/*Clean-up and return*/
	delete input;
	return minabs;
}
/*}}}*/
void TransientInput::SquareMin(IssmDouble* psquaremin,Parameters* parameters){/*{{{*/

	IssmDouble time;

	/*First, recover current time from parameters: */
	parameters->FindParam(&time,TimeEnum);

   /*Retrieve interpolated values for this time step: */
	Input* input=GetTimeInput(time);

	/*Call input function*/
	input->SquareMin(psquaremin,parameters);

	delete input;

}
/*}}}*/
