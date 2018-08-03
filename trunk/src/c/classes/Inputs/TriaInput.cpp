/*!\file TriaInput.c
 * \brief: implementation of the TriaInput object
 */

#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#include "../classes.h"
#include "../../shared/shared.h"

/*TriaInput constructors and destructor*/
TriaInput::TriaInput(){/*{{{*/
	values = NULL;
}
/*}}}*/
TriaInput::TriaInput(int in_enum_type,IssmDouble* in_values,int interpolation_type_in){/*{{{*/

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
TriaInput::~TriaInput(){/*{{{*/
	xDelete<IssmDouble>(this->values);
}
/*}}}*/

/*Object virtual functions definitions:*/
Object* TriaInput::copy() {/*{{{*/

	return new TriaInput(this->enum_type,this->values,this->interpolation_type);

}
/*}}}*/
void TriaInput::DeepEcho(void){/*{{{*/

	_printf_(setw(15)<<"   TriaInput "<<setw(25)<<left<<EnumToStringx(this->enum_type)<<" [");
	for(int i=0;i<this->NumberofNodes(this->interpolation_type);i++) _printf_(" "<<this->values[i]);
	_printf_("] ("<<EnumToStringx(this->interpolation_type)<<")\n");
}
/*}}}*/
void TriaInput::Echo(void){/*{{{*/
	this->DeepEcho();
}
/*}}}*/
int  TriaInput::Id(void){ return -1; }/*{{{*/
/*}}}*/
void TriaInput::Marshall(char** pmarshalled_data,int* pmarshalled_data_size, int marshall_direction){ /*{{{*/

	MARSHALLING_ENUM(TriaInputEnum);

	MARSHALLING(enum_type);
	MARSHALLING(interpolation_type);

	int numnodes = this->NumberofNodes(this->interpolation_type);
	if(numnodes > 0){
		MARSHALLING_DYNAMIC(this->values,IssmDouble,numnodes)
	}
	else this->values = NULL;
}
/*}}}*/
int  TriaInput::ObjectEnum(void){/*{{{*/

	return TriaInputEnum;

}
/*}}}*/

/*TriaInput management*/
int TriaInput::InstanceEnum(void){/*{{{*/

	return this->enum_type;

}
/*}}}*/
int  TriaInput::GetResultArraySize(void){/*{{{*/

	return 1;

}
/*}}}*/
int  TriaInput::GetResultInterpolation(void){/*{{{*/

	if(this->interpolation_type==P0Enum){
		return P0Enum;
	}
	return P1Enum;

}
/*}}}*/
int  TriaInput::GetResultNumberOfNodes(void){/*{{{*/

	return this->NumberofNodes(this->interpolation_type);

}
/*}}}*/
void TriaInput::ResultToPatch(IssmDouble* values,int nodesperelement,int sid){/*{{{*/

	int numnodes = this->NumberofNodes(this->interpolation_type);

	/*Some checks*/
	_assert_(values);
	_assert_(numnodes==nodesperelement);

	/*Fill in arrays*/
	for(int i=0;i<numnodes;i++) values[sid*numnodes + i] = this->values[i];
}
/*}}}*/
Input* TriaInput::SpawnSegInput(int index1,int index2){/*{{{*/

	/*output*/
	SegInput* outinput=NULL;

	if(this->interpolation_type==P0Enum){ 
		outinput=new SegInput(this->enum_type,&this->values[0],P0Enum);
	}
	else{
		/*Assume P1 interpolation only for now*/
		IssmDouble newvalues[2];

		/*Create array of indices depending on location (0=base 1=surface)*/
		newvalues[0]=this->values[index1];
		newvalues[1]=this->values[index2];

		/*Create new Seg input*/
		outinput=new SegInput(this->enum_type,&newvalues[0],P1Enum);
	}

	/*Assign output*/
	return outinput;

}
/*}}}*/
Input* TriaInput::SpawnTriaInput(int index1,int index2,int index3){/*{{{*/

	/*output*/
	TriaInput* outinput=NULL;

	/*Create new Tria input (copy of current input)*/
	outinput=new TriaInput(this->enum_type,&this->values[0],this->interpolation_type);

	/*Assign output*/
	return outinput;

}
/*}}}*/

/*Object functions*/
void TriaInput::ChangeEnum(int newenumtype){/*{{{*/
	this->enum_type=newenumtype;
}
/*}}}*/
void TriaInput::GetInputAllTimeAverages(IssmDouble** pvalues,IssmDouble** ptimes, int* pnumtimes){/*{{{*/

	IssmDouble* outvalues=NULL;
	IssmDouble* times=NULL;
	int         numtimes;

	/*this is not a transient forcing, so we only have 1 value, steady state: */
	numtimes=1;
	outvalues=xNew<IssmDouble>(1);
	times=xNew<IssmDouble>(1);

	this->GetInputAverage(&outvalues[0]);
	times[0]=0.; /*we don't have a time*/

	*pvalues=outvalues;
	*ptimes=times;
	*pnumtimes=numtimes;
}
/*}}}*/
void TriaInput::GetInputAverage(IssmDouble* pvalue){/*{{{*/

	int        numnodes  = this->NumberofNodes(this->interpolation_type);
	IssmDouble numnodesd = reCast<int,IssmDouble>(numnodes);
	IssmDouble value     = 0.;

	for(int i=0;i<numnodes;i++) value+=values[i];
	value = value/numnodesd;

	*pvalue=value;
}
/*}}}*/
void TriaInput::GetInputDerivativeAverageValue(IssmDouble* derivativevalues, IssmDouble* xyz_list){/*{{{*/

	int        numnodes  = this->NumberofNodes(this->interpolation_type);
	IssmDouble numnodesd = reCast<int,IssmDouble>(numnodes);
	IssmDouble dvalue[3];

	derivativevalues[0] = 0.;
	derivativevalues[1] = 0.;

	GaussTria* gauss=new GaussTria();
	for(int iv=0;iv<numnodes;iv++){
		gauss->GaussNode(this->interpolation_type,iv);
		this->GetInputDerivativeValue(&dvalue[0],xyz_list,gauss);

		derivativevalues[0] += dvalue[0]/numnodesd;
		derivativevalues[1] += dvalue[1]/numnodesd;
	}
	delete gauss;
}
/*}}}*/
void TriaInput::GetInputDerivativeValue(IssmDouble* p, IssmDouble* xyz_list, Gauss* gauss){/*{{{*/

	/*Call TriaRef function*/
	_assert_(gauss->Enum()==GaussTriaEnum);
	TriaRef::GetInputDerivativeValue(p,&values[0],xyz_list,(GaussTria*)gauss,this->interpolation_type);
}
/*}}}*/
void TriaInput::GetInputUpToCurrentTimeAverages(IssmDouble** pvalues, IssmDouble** ptimes, int* pnumtimes, IssmDouble currenttime){/*{{{*/

	IssmDouble* outvalues=NULL;
	IssmDouble* times=NULL;
	int         numtimes;

	/*this is not a transient forcing, so we only have 1 value, steady state: */
	numtimes=1;
	outvalues=xNew<IssmDouble>(1);
	times=xNew<IssmDouble>(1);

	this->GetInputAverage(&outvalues[0]);
	times[0]=currenttime; /*we don't have a time*/

	*pvalues=outvalues;
	*ptimes=times;
	*pnumtimes=numtimes;
}
/*}}}*/
void TriaInput::GetInputValue(IssmDouble* pvalue,Gauss* gauss){/*{{{*/

	/*Call TriaRef function*/
	_assert_(gauss->Enum()==GaussTriaEnum);
	TriaRef::GetInputValue(pvalue,&values[0],(GaussTria*)gauss,this->interpolation_type);

}
/*}}}*/

/*Intermediary*/
void TriaInput::AXPY(Input* xinput,IssmDouble scalar){/*{{{*/

	const int numnodes=this->NumberofNodes(this->interpolation_type);
	TriaInput*  xtriainput=NULL;

	/*xinput is of the same type, so cast it: */
	if(xinput->ObjectEnum()!=TriaInputEnum) _error_("Operation not permitted because xinput is of type " << EnumToStringx(xinput->ObjectEnum()));
	xtriainput=(TriaInput*)xinput;
	if(xtriainput->interpolation_type!=this->interpolation_type) _error_("Operation not permitted because xinput is of type " << EnumToStringx(xinput->ObjectEnum()));

	/*Carry out the AXPY operation depending on type:*/
	for(int i=0;i<numnodes;i++)this->values[i]=this->values[i]+scalar*xtriainput->values[i];

}
/*}}}*/
void TriaInput::Configure(Parameters* parameters){/*{{{*/
	/*do nothing: */
}
/*}}}*/
void TriaInput::Constrain(IssmDouble cm_min, IssmDouble cm_max){/*{{{*/

	int i;
	const int numnodes=this->NumberofNodes(this->interpolation_type);

	if(!xIsNan<IssmDouble>(cm_min)) for(i=0;i<numnodes;i++)if (this->values[i]<cm_min)this->values[i]=cm_min;
	if(!xIsNan<IssmDouble>(cm_max)) for(i=0;i<numnodes;i++)if (this->values[i]>cm_max)this->values[i]=cm_max;

}
/*}}}*/
void TriaInput::ConstrainMin(IssmDouble minimum){/*{{{*/

	int numnodes = this->NumberofNodes(this->interpolation_type);
	for(int i=0;i<numnodes;i++) if (values[i]<minimum) values[i]=minimum;
}
/*}}}*/
void TriaInput::GetVectorFromInputs(Vector<IssmDouble>* vector,int* doflist){/*{{{*/
	const int numvertices=3;
	vector->SetValues(numvertices,doflist,this->values,INS_VAL);
} /*}}}*/
IssmDouble TriaInput::InfinityNorm(void){/*{{{*/

	/*Output*/
	IssmDouble norm=0.;
	int numnodes=this->NumberofNodes(this->interpolation_type);

	for(int i=0;i<numnodes;i++) if(fabs(values[i])>norm) norm=fabs(values[i]);
	return norm;
}
/*}}}*/
IssmDouble TriaInput::Max(void){/*{{{*/

	int  numnodes=this->NumberofNodes(this->interpolation_type);
	IssmDouble max=values[0];

	for(int i=1;i<numnodes;i++){
		if(values[i]>max) max=values[i];
	}
	return max;
}
/*}}}*/
IssmDouble TriaInput::MaxAbs(void){/*{{{*/

	int  numnodes=this->NumberofNodes(this->interpolation_type);
	IssmDouble max=fabs(values[0]);

	for(int i=1;i<numnodes;i++){
		if(fabs(values[i])>max) max=fabs(values[i]);
	}
	return max;
}
/*}}}*/
IssmDouble TriaInput::Min(void){/*{{{*/

	const int  numnodes=this->NumberofNodes(this->interpolation_type);
	IssmDouble min=values[0];

	for(int i=1;i<numnodes;i++){
		if(values[i]<min) min=values[i];
	}
	return min;
}
/*}}}*/
IssmDouble TriaInput::MinAbs(void){/*{{{*/

	const int  numnodes=this->NumberofNodes(this->interpolation_type);
	IssmDouble min=fabs(values[0]);

	for(int i=1;i<numnodes;i++){
		if(fabs(values[i])<min) min=fabs(values[i]);
	}
	return min;
}
/*}}}*/
Input* TriaInput::PointwiseMax(Input* inputB){/*{{{*/

	/*Ouput*/
	TriaInput* outinput=NULL;

	/*Intermediaries*/
	int         i;
	TriaInput  *xinputB   = NULL;
	const int   numnodes  = this->NumberofNodes(this->interpolation_type);
	IssmDouble *maxvalues = xNew<IssmDouble>(numnodes);

	/*Check that inputB is of the same type*/
	if(inputB->ObjectEnum()!=TriaInputEnum) _error_("Operation not permitted because inputB is of type " << EnumToStringx(inputB->ObjectEnum()));
	xinputB=(TriaInput*)inputB;
	if(xinputB->interpolation_type!=this->interpolation_type) _error_("Operation not permitted because inputB is of type " << EnumToStringx(xinputB->interpolation_type));

	/*Create point wise max*/
	for(i=0;i<numnodes;i++){
		if(this->values[i] < xinputB->values[i]) maxvalues[i]=xinputB->values[i];
		else maxvalues[i]=this->values[i];
	}

	/*Create new Tria vertex input (copy of current input)*/
	outinput=new TriaInput(this->enum_type,&maxvalues[0],this->interpolation_type);

	/*Return output pointer*/
	xDelete<IssmDouble>(maxvalues);
	return outinput;

}
/*}}}*/
Input* TriaInput::PointwiseMin(Input* inputB){/*{{{*/

	/*Ouput*/
	TriaInput* outinput=NULL;

	/*Intermediaries*/
	int         i;
	TriaInput  *xinputB   = NULL;
	const int   numnodes  = this->NumberofNodes(this->interpolation_type);
	IssmDouble *minvalues = xNew<IssmDouble>(numnodes);

	/*Check that inputB is of the same type*/
	if(inputB->ObjectEnum()!=TriaInputEnum)       _error_("Operation not permitted because inputB is of type " << EnumToStringx(inputB->ObjectEnum()));
	xinputB=(TriaInput*)inputB;
	if(xinputB->interpolation_type!=this->interpolation_type) _error_("Operation not permitted because inputB is of type " << EnumToStringx(xinputB->interpolation_type));

	/*Create point wise min*/
	for(i=0;i<numnodes;i++){
		if(this->values[i] > xinputB->values[i]) minvalues[i]=xinputB->values[i];
		else minvalues[i]=this->values[i];
	}

	/*Create new Tria vertex input (copy of current input)*/
	outinput=new TriaInput(this->enum_type,&minvalues[0],this->interpolation_type);

	/*Return output pointer*/
	xDelete<IssmDouble>(minvalues);
	return outinput;

}
/*}}}*/
Input* TriaInput::PointwiseDivide(Input* inputB){/*{{{*/

	/*Ouput*/
	TriaInput* outinput=NULL;

	/*Intermediaries*/
	TriaInput *xinputB  = NULL;
	const int   numnodes = this->NumberofNodes(this->interpolation_type);

	/*Check that inputB is of the same type*/
	if(inputB->ObjectEnum()!=TriaInputEnum)     _error_("Operation not permitted because inputB is of type " << EnumToStringx(inputB->ObjectEnum()));
	xinputB=(TriaInput*)inputB;
	if(xinputB->interpolation_type!=this->interpolation_type) _error_("Operation not permitted because inputB is of type " << EnumToStringx(xinputB->interpolation_type));

	/*Allocate intermediary*/
	IssmDouble* AdotBvalues=xNew<IssmDouble>(numnodes);

	/*Create point wise division*/
	for(int i=0;i<numnodes;i++){
		_assert_(xinputB->values[i]!=0);
		AdotBvalues[i]=this->values[i]/xinputB->values[i];
	}

	/*Create new Tria vertex input (copy of current input)*/
	outinput=new TriaInput(this->enum_type,AdotBvalues,this->interpolation_type);

	/*Return output pointer*/
	xDelete<IssmDouble>(AdotBvalues);
	return outinput;

}
/*}}}*/
void TriaInput::Scale(IssmDouble scale_factor){/*{{{*/

	const int numnodes=this->NumberofNodes(this->interpolation_type);
	for(int i=0;i<numnodes;i++)values[i]=values[i]*scale_factor;
}
/*}}}*/
void TriaInput::Set(IssmDouble setvalue){/*{{{*/

	const int numnodes=this->NumberofNodes(this->interpolation_type);
	for(int i=0;i<numnodes;i++)values[i]=setvalue;
}
/*}}}*/
void TriaInput::SquareMin(IssmDouble* psquaremin,Parameters* parameters){/*{{{*/

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
