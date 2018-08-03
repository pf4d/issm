/*!\file PentaInput.c
 * \brief: implementation of the PentaInput object
 */

#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#include "../classes.h"
#include "../../shared/shared.h"

/*PentaInput constructors and destructor*/
PentaInput::PentaInput(){/*{{{*/
	values = NULL;
}
/*}}}*/
PentaInput::PentaInput(int in_enum_type,IssmDouble* in_values,int interpolation_type_in){/*{{{*/

	/*Set Enum*/
	this->enum_type=in_enum_type;
	this->interpolation_type=interpolation_type_in;

	int numnodes = this->NumberofNodes(this->interpolation_type);

	/*Set values*/
	if (numnodes > 0){
		this->values=xNew<IssmDouble>((unsigned int)numnodes);
		for(int i=0;i<numnodes;i++) values[i]=in_values[i];
	}
	else{
		this->values = NULL;
	}

}
/*}}}*/
PentaInput::~PentaInput(){/*{{{*/
	xDelete<IssmDouble>(this->values);
}
/*}}}*/

/*Object virtual functions definitions:*/
Object* PentaInput::copy() {/*{{{*/

	return new PentaInput(this->enum_type,this->values,this->interpolation_type);

}
/*}}}*/
void PentaInput::DeepEcho(void){/*{{{*/

	_printf_(setw(15)<<"   PentaInput "<<setw(25)<<left<<EnumToStringx(this->enum_type)<<" [");
	for(int i=0;i<this->NumberofNodes(this->interpolation_type);i++) _printf_(" "<<this->values[i]);
	_printf_("]\n");
}
/*}}}*/
void PentaInput::Echo(void){/*{{{*/
	this->DeepEcho();
}
/*}}}*/
int  PentaInput::Id(void){ return -1; }/*{{{*/
/*}}}*/
void PentaInput::Marshall(char** pmarshalled_data,int* pmarshalled_data_size, int marshall_direction){ /*{{{*/

	MARSHALLING_ENUM(PentaInputEnum);

	MARSHALLING(enum_type);
	MARSHALLING(interpolation_type);

	int numnodes = this->NumberofNodes(this->interpolation_type);
	if(numnodes > 0){
		MARSHALLING_DYNAMIC(this->values,IssmDouble,numnodes)
	}
	else this->values = NULL;
}
/*}}}*/
int  PentaInput::ObjectEnum(void){/*{{{*/

	return PentaInputEnum;

}
/*}}}*/

/*PentaInput management*/
int  PentaInput::GetResultInterpolation(void){/*{{{*/

	if(this->interpolation_type==P0Enum){
		return P0Enum;
	}
	return P1Enum;

}
/*}}}*/
int  PentaInput::GetResultNumberOfNodes(void){/*{{{*/

	return this->NumberofNodes(this->interpolation_type);;

}
/*}}}*/
int PentaInput::InstanceEnum(void){/*{{{*/

	return this->enum_type;

}
/*}}}*/
void PentaInput::ResultToPatch(IssmDouble* values,int nodesperelement,int sid){/*{{{*/

	int numnodes = this->NumberofNodes(this->interpolation_type);

	/*Some checks*/
	_assert_(values);
	_assert_(numnodes==nodesperelement);

	/*Fill in arrays*/
	for(int i=0;i<numnodes;i++) values[sid*numnodes + i] = this->values[i];
}
/*}}}*/
Input* PentaInput::SpawnSegInput(int index1,int index2){/*{{{*/

	_error_("not supported");
}
/*}}}*/
Input* PentaInput::SpawnTriaInput(int index1,int index2,int index3){/*{{{*/

	/*output*/
	TriaInput* outinput=NULL;

	if(this->interpolation_type==P0Enum){ 
		outinput=new TriaInput(this->enum_type,&this->values[0],P0Enum);
	}
	else{
		/*Assume P1 interpolation only for now*/
		IssmDouble newvalues[3]; 

		/*Create array of indices depending on location (0=base 1=surface)*/
		int indices[3];
		indices[0] = index1;
		indices[1] = index2;
		indices[2] = index3;

		/*Create new input*/
		for(int i=0;i<3;i++){
			_assert_(indices[i]>=0 && indices[i]<6);
			newvalues[i]=this->values[indices[i]];
		}
		outinput=new TriaInput(this->enum_type,&newvalues[0],P1Enum);
	}

	/*Assign output*/
	return outinput;
}
/*}}}*/

/*Object functions*/
void PentaInput::ChangeEnum(int newenumtype){/*{{{*/
	this->enum_type=newenumtype;
}
/*}}}*/
void PentaInput::GetInputAverage(IssmDouble* pvalue){/*{{{*/

	int        numnodes  = this->NumberofNodes(this->interpolation_type);
	IssmDouble numnodesd = reCast<int,IssmDouble>(numnodes);
	IssmDouble value     = 0.;

	for(int i=0;i<numnodes;i++) value+=values[i];
	value = value/numnodesd;

	*pvalue=value;
}
/*}}}*/
void PentaInput::GetInputDerivativeValue(IssmDouble* p, IssmDouble* xyz_list,Gauss* gauss){/*{{{*/

	/*Call PentaRef function*/
	_assert_(gauss->Enum()==GaussPentaEnum);
	PentaRef::GetInputDerivativeValue(p,&values[0],xyz_list,(GaussPenta*)gauss,this->interpolation_type);
}
/*}}}*/
void PentaInput::GetInputValue(IssmDouble* pvalue){/*{{{*/

	if(this->interpolation_type==P0Enum){
		pvalue=&values[0];
	}
	else _error_("not implemented yet");
}
/*}}}*/
void PentaInput::GetInputValue(IssmDouble* pvalue,Gauss* gauss){/*{{{*/

	/*Call PentaRef function*/
	_assert_(gauss->Enum()==GaussPentaEnum);
	PentaRef::GetInputValue(pvalue,&values[0],(GaussPenta*)gauss,this->interpolation_type);

}
/*}}}*/

/*Intermediary*/
void PentaInput::AXPY(Input* xinput,IssmDouble scalar){/*{{{*/

	const int numnodes=this->NumberofNodes(this->interpolation_type);
	PentaInput* xpentainput=NULL;

	/*If xinput is a ControlInput, take its values directly*/
	if(xinput->ObjectEnum()==ControlInputEnum){
		xinput=((ControlInput*)xinput)->values;
	}

	/*xinput is of the same type, so cast it: */
	if(xinput->ObjectEnum()!=PentaInputEnum)
	  _error_("Operation not permitted because xinput is of type " << EnumToStringx(xinput->ObjectEnum()));
	xpentainput=(PentaInput*)xinput;
	if(xpentainput->interpolation_type!=this->interpolation_type) _error_("Operation not permitted because xinput is of type " << EnumToStringx(xpentainput->interpolation_type));

	/*Carry out the AXPY operation depending on type:*/
	for(int i=0;i<numnodes;i++)this->values[i]=this->values[i]+scalar*xpentainput->values[i];

}
/*}}}*/
void PentaInput::Configure(Parameters* parameters){/*{{{*/
	/*do nothing: */
}
/*}}}*/
void PentaInput::Constrain(IssmDouble cm_min, IssmDouble cm_max){/*{{{*/

	int i;
	const int numnodes=this->NumberofNodes(this->interpolation_type);

	if(!xIsNan<IssmDouble>(cm_min)) for(i=0;i<numnodes;i++)if (this->values[i]<cm_min)this->values[i]=cm_min;
	if(!xIsNan<IssmDouble>(cm_max)) for(i=0;i<numnodes;i++)if (this->values[i]>cm_max)this->values[i]=cm_max;

}
/*}}}*/
void PentaInput::ConstrainMin(IssmDouble minimum){/*{{{*/

	int numnodes = this->NumberofNodes(this->interpolation_type);
	for(int i=0;i<numnodes;i++) if (values[i]<minimum) values[i]=minimum;
}
/*}}}*/
void PentaInput::Extrude(int start){/*{{{*/

	switch(this->interpolation_type){
		case P1Enum:
			if(start==-1){
				for(int i=0;i<3;i++) this->values[3+i]=this->values[i];
			}
			else{
				for(int i=0;i<3;i++) this->values[i]  =this->values[3+i];
			}
			break;
		default:
			_error_("not supported yet for type "<<EnumToStringx(this->interpolation_type));
	}
}
/*}}}*/
void PentaInput::GetVectorFromInputs(Vector<IssmDouble>* vector,int* doflist){/*{{{*/
	const int numvertices=6;
	vector->SetValues(numvertices,doflist,this->values,INS_VAL);
} /*}}}*/
IssmDouble PentaInput::InfinityNorm(void){/*{{{*/

	/*Output*/
	IssmDouble norm=0.;
	int numnodes=this->NumberofNodes(this->interpolation_type);

	for(int i=0;i<numnodes;i++) if(fabs(values[i])>norm) norm=fabs(values[i]);
	return norm;
}
/*}}}*/
IssmDouble PentaInput::Max(void){/*{{{*/

	int  numnodes=this->NumberofNodes(this->interpolation_type);
	IssmDouble max=values[0];

	for(int i=1;i<numnodes;i++){
		if(values[i]>max) max=values[i];
	}
	return max;
}
/*}}}*/
IssmDouble PentaInput::MaxAbs(void){/*{{{*/

	int  numnodes=this->NumberofNodes(this->interpolation_type);
	IssmDouble max=fabs(values[0]);

	for(int i=1;i<numnodes;i++){
		if(fabs(values[i])>max) max=fabs(values[i]);
	}
	return max;
}
/*}}}*/
IssmDouble PentaInput::Min(void){/*{{{*/

	const int  numnodes=this->NumberofNodes(this->interpolation_type);
	IssmDouble min=values[0];

	for(int i=1;i<numnodes;i++){
		if(values[i]<min) min=values[i];
	}
	return min;
}
/*}}}*/
IssmDouble PentaInput::MinAbs(void){/*{{{*/

	const int  numnodes=this->NumberofNodes(this->interpolation_type);
	IssmDouble min=fabs(values[0]);

	for(int i=1;i<numnodes;i++){
		if(fabs(values[i])<min) min=fabs(values[i]);
	}
	return min;
}
/*}}}*/
Input* PentaInput::PointwiseDivide(Input* inputB){/*{{{*/

	/*Ouput*/
	PentaInput* outinput=NULL;

	/*Intermediaries*/
	PentaInput *xinputB  = NULL;
	const int   numnodes = this->NumberofNodes(this->interpolation_type);

	/*Check that inputB is of the same type*/
	if(inputB->ObjectEnum()!=PentaInputEnum)     _error_("Operation not permitted because inputB is of type " << EnumToStringx(inputB->ObjectEnum()));
	xinputB=(PentaInput*)inputB;
	if(xinputB->interpolation_type!=this->interpolation_type) _error_("Operation not permitted because inputB is of type " << EnumToStringx(xinputB->interpolation_type));

	/*Allocate intermediary*/
	IssmDouble* AdotBvalues=xNew<IssmDouble>(numnodes);

	/*Create point wise sum*/
	for(int i=0;i<numnodes;i++){
		_assert_(xinputB->values[i]!=0);
		AdotBvalues[i]=this->values[i]/xinputB->values[i];
	}

	/*Create new Penta vertex input (copy of current input)*/
	outinput=new PentaInput(this->enum_type,AdotBvalues,this->interpolation_type);

	/*Return output pointer*/
	xDelete<IssmDouble>(AdotBvalues);
	return outinput;

}
/*}}}*/
Input* PentaInput::PointwiseMax(Input* inputB){/*{{{*/

	/*Ouput*/
	PentaInput* outinput=NULL;

	/*Intermediaries*/
	int         i;
	PentaInput  *xinputB   = NULL;
	const int   numnodes  = this->NumberofNodes(this->interpolation_type);
	IssmDouble *maxvalues = xNew<IssmDouble>(numnodes);

	/*Check that inputB is of the same type*/
	if(inputB->ObjectEnum()!=PentaInputEnum) _error_("Operation not permitted because inputB is of type " << EnumToStringx(inputB->ObjectEnum()));
	xinputB=(PentaInput*)inputB;
	if(xinputB->interpolation_type!=this->interpolation_type) _error_("Operation not permitted because inputB is of type " << EnumToStringx(xinputB->interpolation_type));

	/*Create point wise max*/
	for(i=0;i<numnodes;i++){
		if(this->values[i] < xinputB->values[i]) maxvalues[i]=xinputB->values[i];
		else maxvalues[i]=this->values[i];
	}

	/*Create new Penta vertex input (copy of current input)*/
	outinput=new PentaInput(this->enum_type,&maxvalues[0],this->interpolation_type);

	/*Return output pointer*/
	xDelete<IssmDouble>(maxvalues);
	return outinput;
}
/*}}}*/
Input* PentaInput::PointwiseMin(Input* inputB){/*{{{*/

	/*Ouput*/
	PentaInput* outinput=NULL;

	/*Intermediaries*/
	int         i;
	PentaInput  *xinputB   = NULL;
	const int   numnodes  = this->NumberofNodes(this->interpolation_type);
	IssmDouble *minvalues = xNew<IssmDouble>(numnodes);

	/*Check that inputB is of the same type*/
	if(inputB->ObjectEnum()!=PentaInputEnum)       _error_("Operation not permitted because inputB is of type " << EnumToStringx(inputB->ObjectEnum()));
	xinputB=(PentaInput*)inputB;
	if(xinputB->interpolation_type!=this->interpolation_type) _error_("Operation not permitted because inputB is of type " << EnumToStringx(xinputB->interpolation_type));

	/*Create point wise min*/
	for(i=0;i<numnodes;i++){
		if(this->values[i] > xinputB->values[i]) minvalues[i]=xinputB->values[i];
		else minvalues[i]=this->values[i];
	}

	/*Create new Penta vertex input (copy of current input)*/
	outinput=new PentaInput(this->enum_type,&minvalues[0],this->interpolation_type);

	/*Return output pointer*/
	xDelete<IssmDouble>(minvalues);
	return outinput;
}
/*}}}*/
void PentaInput::Scale(IssmDouble scale_factor){/*{{{*/

	const int numnodes=this->NumberofNodes(this->interpolation_type);
	for(int i=0;i<numnodes;i++)values[i]=values[i]*scale_factor;
}
/*}}}*/
void PentaInput::SquareMin(IssmDouble* psquaremin,Parameters* parameters){/*{{{*/

	int        numnodes=this->NumberofNodes(this->interpolation_type);
	IssmDouble squaremin;

	/*Now, figure out minimum of valuescopy: */
	squaremin=pow(this->values[0],2);
	for(int i=1;i<numnodes;i++){
		if(pow(this->values[i],2)<squaremin)squaremin=pow(this->values[i],2);
	}
	/*Assign output pointers:*/
	*psquaremin=squaremin;
}
/*}}}*/
void PentaInput::VerticallyIntegrate(Input* thickness_input){/*{{{*/

	IssmDouble thickness;
	IssmDouble value=0.;

	/*Check that input provided is a thickness*/
	if (thickness_input->InstanceEnum()!=ThicknessEnum) _error_("Input provided is not a Thickness (enum_type is " << EnumToStringx(thickness_input->InstanceEnum()) << ")");

	/*vertically integrate depending on type (and use P1 interpolation from now on)*/
	switch(this->interpolation_type){
		case P1Enum:
		case P1bubbleEnum:
		case P1xP2Enum:
		case P1xP3Enum:
		case P2Enum:
			  {
				this->interpolation_type=P1Enum;
				GaussPenta *gauss=new GaussPenta();
				for(int iv=0;iv<3;iv++){
					gauss->GaussVertex(iv);
					thickness_input->GetInputValue(&thickness,gauss);
					this->values[iv]=0.5*(this->values[iv]+this->values[iv+3]) * thickness;
					this->values[iv+3]=this->values[iv];
				}
				delete gauss;
				return; 
			  }
		default:
			_error_("not supported yet for type "<<EnumToStringx(this->interpolation_type));
	}
}
/*}}}*/
