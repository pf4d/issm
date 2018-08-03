/*!\file TetraInput.c
 * \brief: implementation of the TetraInput object
 */

#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#include "../classes.h"
#include "../../shared/shared.h"

/*TetraInput constructors and destructor*/
TetraInput::TetraInput(){/*{{{*/
	values = NULL;
}
/*}}}*/
TetraInput::TetraInput(int in_enum_type,IssmDouble* in_values,int interpolation_type_in){/*{{{*/

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
TetraInput::~TetraInput(){/*{{{*/
	xDelete<IssmDouble>(this->values);
}
/*}}}*/

/*Object virtual functions definitions:*/
Object* TetraInput::copy() {/*{{{*/

	return new TetraInput(this->enum_type,this->values,this->interpolation_type);

}
/*}}}*/
void TetraInput::DeepEcho(void){/*{{{*/

	_printf_(setw(15)<<"   TetraInput "<<setw(25)<<left<<EnumToStringx(this->enum_type)<<" [");
	for(int i=0;i<this->NumberofNodes(this->interpolation_type);i++) _printf_(" "<<this->values[i]);
	_printf_("] ("<<EnumToStringx(this->interpolation_type)<<")\n");
}
/*}}}*/
void TetraInput::Echo(void){/*{{{*/
	this->DeepEcho();
}
/*}}}*/
int  TetraInput::Id(void){ return -1; }/*{{{*/
/*}}}*/
void TetraInput::Marshall(char** pmarshalled_data,int* pmarshalled_data_size, int marshall_direction){ /*{{{*/

	MARSHALLING_ENUM(TetraInputEnum);

	MARSHALLING(enum_type);
	MARSHALLING(interpolation_type);

	int numnodes = this->NumberofNodes(this->interpolation_type);
	if(numnodes > 0){
		MARSHALLING_DYNAMIC(this->values,IssmDouble,numnodes)
	}
	else this->values = NULL;
}
/*}}}*/
int  TetraInput::ObjectEnum(void){/*{{{*/

	return TetraInputEnum;

}
/*}}}*/

/*TetraInput management*/
int TetraInput::InstanceEnum(void){/*{{{*/

	return this->enum_type;

}
/*}}}*/
int  TetraInput::GetResultInterpolation(void){/*{{{*/

	if(this->interpolation_type==P0Enum){
		return P0Enum;
	}
	return P1Enum;

}
/*}}}*/
int  TetraInput::GetResultNumberOfNodes(void){/*{{{*/

	return this->NumberofNodes(this->interpolation_type);

}
/*}}}*/
void TetraInput::ResultToPatch(IssmDouble* values,int nodesperelement,int sid){/*{{{*/

	int numnodes = this->NumberofNodes(this->interpolation_type);

	/*Some checks*/
	_assert_(values);
	_assert_(numnodes==nodesperelement);

	/*Fill in arrays*/
	for(int i=0;i<numnodes;i++) values[sid*numnodes + i] = this->values[i];
}
/*}}}*/

/*Object functions*/
void TetraInput::ChangeEnum(int newenumtype){/*{{{*/
	this->enum_type=newenumtype;
}
/*}}}*/
void TetraInput::GetInputAllTimeAverages(IssmDouble** pvalues,IssmDouble** ptimes, int* pnumtimes){/*{{{*/

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
void TetraInput::GetInputAverage(IssmDouble* pvalue){/*{{{*/

	int        numnodes  = this->NumberofNodes(this->interpolation_type);
	IssmDouble numnodesd = reCast<int,IssmDouble>(numnodes);
	IssmDouble value     = 0.;

	for(int i=0;i<numnodes;i++) value+=values[i];
	value = value/numnodesd;

	*pvalue=value;
}
/*}}}*/
void TetraInput::GetInputDerivativeValue(IssmDouble* p, IssmDouble* xyz_list, Gauss* gauss){/*{{{*/

	/*Call TetraRef function*/
	_assert_(gauss->Enum()==GaussTetraEnum);
	TetraRef::GetInputDerivativeValue(p,&values[0],xyz_list,(GaussTetra*)gauss,this->interpolation_type);
}
/*}}}*/
void TetraInput::GetInputUpToCurrentTimeAverages(IssmDouble** pvalues, IssmDouble** ptimes, int* pnumtimes, IssmDouble currenttime){/*{{{*/

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
void TetraInput::GetInputValue(IssmDouble* pvalue,Gauss* gauss){/*{{{*/

	/*Call TetraRef function*/
	_assert_(gauss->Enum()==GaussTetraEnum);
	TetraRef::GetInputValue(pvalue,&values[0],(GaussTetra*)gauss,this->interpolation_type);

}
/*}}}*/
Input* TetraInput::SpawnTriaInput(int index1,int index2,int index3){/*{{{*/

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
			_assert_(indices[i]>=0 && indices[i]<4);
			newvalues[i]=this->values[indices[i]];
		}
		outinput=new TriaInput(this->enum_type,&newvalues[0],P1Enum);
	}

	/*Assign output*/
	return outinput;
}
/*}}}*/

/*Intermediary*/
void TetraInput::AXPY(Input* xinput,IssmDouble scalar){/*{{{*/

	const int numnodes=this->NumberofNodes(this->interpolation_type);
	TetraInput*  xtriainput=NULL;

	/*xinput is of the same type, so cast it: */
	if(xinput->ObjectEnum()!=TetraInputEnum) _error_("Operation not permitted because xinput is of type " << EnumToStringx(xinput->ObjectEnum()));
	xtriainput=(TetraInput*)xinput;
	if(xtriainput->interpolation_type!=this->interpolation_type) _error_("Operation not permitted because xinput is of type " << EnumToStringx(xinput->ObjectEnum()));

	/*Carry out the AXPY operation depending on type:*/
	for(int i=0;i<numnodes;i++)this->values[i]=this->values[i]+scalar*xtriainput->values[i];

}
/*}}}*/
void TetraInput::Configure(Parameters* parameters){/*{{{*/
	/*do nothing: */
}
/*}}}*/
void TetraInput::Constrain(IssmDouble cm_min, IssmDouble cm_max){/*{{{*/

	int i;
	const int numnodes=this->NumberofNodes(this->interpolation_type);

	if(!xIsNan<IssmDouble>(cm_min)) for(i=0;i<numnodes;i++)if (this->values[i]<cm_min)this->values[i]=cm_min;
	if(!xIsNan<IssmDouble>(cm_max)) for(i=0;i<numnodes;i++)if (this->values[i]>cm_max)this->values[i]=cm_max;

}
/*}}}*/
void TetraInput::ConstrainMin(IssmDouble minimum){/*{{{*/

	int numnodes = this->NumberofNodes(this->interpolation_type);
	for(int i=0;i<numnodes;i++) if (values[i]<minimum) values[i]=minimum;
}
/*}}}*/
void TetraInput::GetVectorFromInputs(Vector<IssmDouble>* vector,int* doflist){/*{{{*/
	const int numvertices=4;
	vector->SetValues(numvertices,doflist,this->values,INS_VAL);
} /*}}}*/
IssmDouble TetraInput::InfinityNorm(void){/*{{{*/

	/*Output*/
	IssmDouble norm=0.;
	int numnodes=this->NumberofNodes(this->interpolation_type);

	for(int i=0;i<numnodes;i++) if(fabs(values[i])>norm) norm=fabs(values[i]);
	return norm;
}
/*}}}*/
IssmDouble TetraInput::MinAbs(void){/*{{{*/

	const int  numnodes=this->NumberofNodes(this->interpolation_type);
	IssmDouble min=fabs(values[0]);

	for(int i=1;i<numnodes;i++){
		if(fabs(values[i])<min) min=fabs(values[i]);
	}
	return min;
}
/*}}}*/
IssmDouble TetraInput::Max(void){/*{{{*/

	int  numnodes=this->NumberofNodes(this->interpolation_type);
	IssmDouble max=values[0];

	for(int i=1;i<numnodes;i++){
		if(values[i]>max) max=values[i];
	}
	return max;
}
/*}}}*/
IssmDouble TetraInput::MaxAbs(void){/*{{{*/

	int  numnodes=this->NumberofNodes(this->interpolation_type);
	IssmDouble max=fabs(values[0]);

	for(int i=1;i<numnodes;i++){
		if(fabs(values[i])>max) max=fabs(values[i]);
	}
	return max;
}
/*}}}*/
IssmDouble TetraInput::Min(void){/*{{{*/

	const int  numnodes=this->NumberofNodes(this->interpolation_type);
	IssmDouble min=values[0];

	for(int i=1;i<numnodes;i++){
		if(values[i]<min) min=values[i];
	}
	return min;
}
/*}}}*/
Input* TetraInput::PointwiseDivide(Input* inputB){/*{{{*/

	/*Ouput*/
	TetraInput* outinput=NULL;

	/*Intermediaries*/
	TetraInput *xinputB  = NULL;
	const int   numnodes = this->NumberofNodes(this->interpolation_type);

	/*Check that inputB is of the same type*/
	if(inputB->ObjectEnum()!=TetraInputEnum)     _error_("Operation not permitted because inputB is of type " << EnumToStringx(inputB->ObjectEnum()));
	xinputB=(TetraInput*)inputB;
	if(xinputB->interpolation_type!=this->interpolation_type) _error_("Operation not permitted because inputB is of type " << EnumToStringx(xinputB->interpolation_type));

	/*Allocate intermediary*/
	IssmDouble* AdotBvalues=xNew<IssmDouble>(numnodes);

	/*Create point wise division*/
	for(int i=0;i<numnodes;i++){
		_assert_(xinputB->values[i]!=0);
		AdotBvalues[i]=this->values[i]/xinputB->values[i];
	}

	/*Create new Tetra vertex input (copy of current input)*/
	outinput=new TetraInput(this->enum_type,AdotBvalues,this->interpolation_type);

	/*Return output pointer*/
	xDelete<IssmDouble>(AdotBvalues);
	return outinput;

}
/*}}}*/
Input* TetraInput::PointwiseMax(Input* inputB){/*{{{*/

	/*Ouput*/
	TetraInput* outinput=NULL;

	/*Intermediaries*/
	int         i;
	TetraInput  *xinputB   = NULL;
	const int   numnodes  = this->NumberofNodes(this->interpolation_type);
	IssmDouble *maxvalues = xNew<IssmDouble>(numnodes);

	/*Check that inputB is of the same type*/
	if(inputB->ObjectEnum()!=TetraInputEnum) _error_("Operation not permitted because inputB is of type " << EnumToStringx(inputB->ObjectEnum()));
	xinputB=(TetraInput*)inputB;
	if(xinputB->interpolation_type!=this->interpolation_type) _error_("Operation not permitted because inputB is of type " << EnumToStringx(xinputB->interpolation_type));

	/*Create point wise max*/
	for(i=0;i<numnodes;i++){
		if(this->values[i] < xinputB->values[i]) maxvalues[i]=xinputB->values[i];
		else maxvalues[i]=this->values[i];
	}

	/*Create new Tetra vertex input (copy of current input)*/
	outinput=new TetraInput(this->enum_type,&maxvalues[0],this->interpolation_type);

	/*Return output pointer*/
	xDelete<IssmDouble>(maxvalues);
	return outinput;

}
/*}}}*/
Input* TetraInput::PointwiseMin(Input* inputB){/*{{{*/

	/*Ouput*/
	TetraInput* outinput=NULL;

	/*Intermediaries*/
	int         i;
	TetraInput  *xinputB   = NULL;
	const int   numnodes  = this->NumberofNodes(this->interpolation_type);
	IssmDouble *minvalues = xNew<IssmDouble>(numnodes);

	/*Check that inputB is of the same type*/
	if(inputB->ObjectEnum()!=TetraInputEnum)       _error_("Operation not permitted because inputB is of type " << EnumToStringx(inputB->ObjectEnum()));
	xinputB=(TetraInput*)inputB;
	if(xinputB->interpolation_type!=this->interpolation_type) _error_("Operation not permitted because inputB is of type " << EnumToStringx(xinputB->interpolation_type));

	/*Create point wise min*/
	for(i=0;i<numnodes;i++){
		if(this->values[i] > xinputB->values[i]) minvalues[i]=xinputB->values[i];
		else minvalues[i]=this->values[i];
	}

	/*Create new Tetra vertex input (copy of current input)*/
	outinput=new TetraInput(this->enum_type,&minvalues[0],this->interpolation_type);

	/*Return output pointer*/
	xDelete<IssmDouble>(minvalues);
	return outinput;

}
/*}}}*/
void TetraInput::Scale(IssmDouble scale_factor){/*{{{*/

	const int numnodes=this->NumberofNodes(this->interpolation_type);
	for(int i=0;i<numnodes;i++)values[i]=values[i]*scale_factor;
}
/*}}}*/
void TetraInput::Set(IssmDouble setvalue){/*{{{*/

	const int numnodes=this->NumberofNodes(this->interpolation_type);
	for(int i=0;i<numnodes;i++)values[i]=setvalue;
}
/*}}}*/
void TetraInput::SquareMin(IssmDouble* psquaremin,Parameters* parameters){/*{{{*/

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
