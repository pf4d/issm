/*!\file PetscVec.cpp
 * \brief: implementation of the PetscVec object
 */

/*Headers:*/
/*{{{*/
#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#include <stdio.h>
#include <string.h>
#include "../petscincludes.h"
#include "../../../shared/shared.h"

/*}}}*/

/*PetscVec constructors and destructor*/
PetscVec::PetscVec(){/*{{{*/
	this->vector=NULL;
	#ifdef _HAVE_ADOLC_
	this->avector=NULL;
	#endif
}
/*}}}*/
PetscVec::PetscVec(int M,bool fromlocalsize){/*{{{*/

	this->vector=NewVec(M,IssmComm::GetComm(),fromlocalsize);

}
/*}}}*/
PetscVec::PetscVec(int m,int M){/*{{{*/

	VecCreate(IssmComm::GetComm(),&this->vector);
	VecSetSizes(this->vector,m,M);
	VecSetFromOptions(this->vector);
	VecSet(this->vector,0.0);

}
/*}}}*/
PetscVec::PetscVec(Vec petsc_vec){/*{{{*/

	if(petsc_vec==NULL){
		this->vector=NewVec(0,IssmComm::GetComm());
	}
	else{
		/*copy vector*/
		VecDuplicate(petsc_vec,&this->vector);
		VecCopy(petsc_vec,this->vector);
	}

}
/*}}}*/
PetscVec::PetscVec(IssmDouble* serial_vec,int M){/*{{{*/

	int* idxm=NULL;
	if(M)idxm=xNew<int>(M);
	for(int i=0;i<M;i++) idxm[i]=i;

	this->vector=NewVec(M,IssmComm::GetComm());
	VecSetValues(this->vector,M,idxm,serial_vec,INSERT_VALUES);
	VecAssemblyBegin(this->vector);
	VecAssemblyEnd(this->vector);

	xDelete<int>(idxm);
}
/*}}}*/
PetscVec::~PetscVec(){/*{{{*/
    VecFree(&this->vector);
}
/*}}}*/

/*PetscVec specific routines: */
void PetscVec::Echo(void){/*{{{*/

	_assert_(this->vector);
	VecView(this->vector,PETSC_VIEWER_STDOUT_WORLD);
}
/*}}}*/
void PetscVec::Assemble(void){/*{{{*/

	_assert_(this->vector);
	VecAssemblyBegin(this->vector); 
	VecAssemblyEnd(this->vector);

}
/*}}}*/
void PetscVec::SetValues(int ssize, int* list, IssmDouble* values, InsMode mode){/*{{{*/

	_assert_(this->vector);
	VecSetValues(this->vector,ssize,list,values,ISSMToPetscInsertMode(mode));

}
/*}}}*/
void PetscVec::SetValue(int dof, IssmDouble value, InsMode mode){/*{{{*/

	_assert_(this->vector);
	VecSetValues(this->vector,1,&dof,&value,ISSMToPetscInsertMode(mode));

}
/*}}}*/
void PetscVec::GetValue(IssmDouble* pvalue,int dof){/*{{{*/

	_assert_(this->vector);
	VecGetValues(this->vector,1,&dof,pvalue);

}
/*}}}*/
void PetscVec::GetSize(int* pM){/*{{{*/

	_assert_(this->vector);
	VecGetSize(this->vector,pM);

}
/*}}}*/
void PetscVec::GetLocalSize(int* pM){/*{{{*/

	_assert_(this->vector);
	VecGetLocalSize(this->vector,pM);

}
/*}}}*/
PetscVec* PetscVec::Duplicate(void){/*{{{*/

	_assert_(this->vector);
	Vec vec_output=NULL;
	VecDuplicate(this->vector,&vec_output);
	PetscVec* output=new PetscVec(vec_output);
	VecFree(&vec_output);

	return output;
}
/*}}}*/
void PetscVec::Set(IssmDouble value){/*{{{*/

	_assert_(this->vector);
	VecSet(this->vector,value);

}
/*}}}*/
void PetscVec::AXPY(PetscVec* X, IssmDouble a){/*{{{*/

	_assert_(this->vector);
	VecAXPY(this->vector,a,X->vector);

}
/*}}}*/
void PetscVec::AYPX(PetscVec* X, IssmDouble a){/*{{{*/

	_assert_(this->vector);
	VecAYPX(this->vector,a,X->vector);

}
/*}}}*/
IssmDouble* PetscVec::ToMPISerial(void){/*{{{*/

	IssmDouble* vec_serial=NULL;
	VecToMPISerial(&vec_serial, this->vector,IssmComm::GetComm());
	return vec_serial;

}
/*}}}*/
void PetscVec::Shift(IssmDouble shift){/*{{{*/

	if(this->vector) VecShift(this->vector,shift);

}
/*}}}*/
void PetscVec::Copy(PetscVec* to){/*{{{*/

	if(this->vector) VecCopy(this->vector,to->vector);

}
/*}}}*/
IssmDouble PetscVec::Max(void){/*{{{*/

	_assert_(this->vector);

	IssmDouble max;
	VecMax(this->vector,NULL,&max);
	return max;

}
/*}}}*/
IssmDouble PetscVec::Norm(NormMode mode){/*{{{*/

	IssmDouble norm=0;
	_assert_(this->vector);
	VecNorm(this->vector,ISSMToPetscNormMode(mode),&norm);
	return norm;

}
/*}}}*/
void PetscVec::Scale(IssmDouble scale_factor){/*{{{*/

	_assert_(this->vector);
	VecScale(this->vector,scale_factor); 

}
/*}}}*/
IssmDouble PetscVec::Dot(PetscVec* input){/*{{{*/

	IssmDouble dot;
	_assert_(this->vector);
	VecDot(this->vector,input->vector,&dot);
	return dot;

}
/*}}}*/
void PetscVec::PointwiseDivide(PetscVec* x,PetscVec* y){/*{{{*/

	_assert_(this->vector);
	VecPointwiseDivide(this->vector,x->vector,y->vector);

}
/*}}}*/
