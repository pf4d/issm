/*!\file ElementMatrix.cpp
 * \brief: implementation of the ElementMatrix object, used to plug values from element into global stiffness matrix
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
#include "../classes.h"
#include "../../shared/shared.h"
/*}}}*/

/*ElementMatrix constructors and destructor*/
ElementMatrix::ElementMatrix(){/*{{{*/

	this->nrows=0;
	this->ncols=0;
	this->values=NULL;
	this->dofsymmetrical=false;

	this->row_fsize=0;
	this->row_flocaldoflist=NULL;
	this->row_fglobaldoflist=NULL;
	this->row_ssize=0;
	this->row_slocaldoflist=NULL;
	this->row_sglobaldoflist=NULL;

	this->col_fsize=0;
	this->col_flocaldoflist=NULL;
	this->col_fglobaldoflist=NULL;
	this->col_ssize=0;
	this->col_slocaldoflist=NULL;
	this->col_sglobaldoflist=NULL;

}
/*}}}*/
ElementMatrix::ElementMatrix(ElementMatrix* Ke){/*{{{*/

	if(!Ke) _error_("Input Element Matrix is a NULL pointer");
	this->Init(Ke);
	return;
}
/*}}}*/
ElementMatrix::ElementMatrix(ElementMatrix* Ke1, ElementMatrix* Ke2){/*{{{*/

	/*intermediaries*/
	int i,j,counter;
	int gsize,fsize,ssize;
	int* P=NULL;
	bool found;

	/*If one of the two matrix is NULL, we copy the other one*/
	if(!Ke1 && !Ke2){
		_error_("Two input element matrices are NULL");
	}
	else if(!Ke1){
		this->Init(Ke2);
		return;
	}
	else if(!Ke2){
		this->Init(Ke1);
		return;
	}

	/*General Case: Ke1 and Ke2 are not empty*/
	if(!Ke1->dofsymmetrical || !Ke2->dofsymmetrical) _error_("merging 2 non dofsymmetrical matrices not implemented yet");

	/*Initialize itransformation matrix Ke[P[i]] = Ke2[i]*/
	P=xNew<int>(Ke2->nrows);

	/*1: Get the new numbering of Ke2 and get size of the new matrix*/
	gsize=Ke1->nrows;
	for(i=0;i<Ke2->nrows;i++){
		found=false;
		for(j=0;j<Ke1->nrows;j++){
			if(Ke2->gglobaldoflist[i]==Ke1->gglobaldoflist[j]){
				found=true; P[i]=j; break;
			}
		}
		if(!found){
			P[i]=gsize; gsize++;
		}
	}

	/*2: Initialize static fields*/
	this->nrows=gsize;
	this->ncols=gsize;
	this->dofsymmetrical=true;

	/*Gset and values*/
	this->gglobaldoflist=xNew<int>(this->nrows);
	this->values=xNewZeroInit<IssmDouble>(this->nrows*this->ncols);
	for(i=0;i<Ke1->nrows;i++){
		for(j=0;j<Ke1->ncols;j++){
			this->values[i*this->ncols+j] += Ke1->values[i*Ke1->ncols+j];
		}
		this->gglobaldoflist[i]=Ke1->gglobaldoflist[i];
	}
	for(i=0;i<Ke2->nrows;i++){
		for(j=0;j<Ke2->ncols;j++){
			this->values[P[i]*this->ncols+P[j]] += Ke2->values[i*Ke2->ncols+j];
		}
		this->gglobaldoflist[P[i]]=Ke2->gglobaldoflist[i];
	}

	/*Fset*/
	fsize=Ke1->row_fsize;
	for(i=0;i<Ke2->row_fsize;i++){
		if(P[Ke2->row_flocaldoflist[i]] >= Ke1->nrows) fsize++;
	}
	this->row_fsize=fsize;
	if(fsize){
		this->row_flocaldoflist =xNew<int>(fsize);
		this->row_fglobaldoflist=xNew<int>(fsize);
		for(i=0;i<Ke1->row_fsize;i++){
			this->row_flocaldoflist[i] =Ke1->row_flocaldoflist[i];
			this->row_fglobaldoflist[i]=Ke1->row_fglobaldoflist[i];
		}
		counter=Ke1->row_fsize;
		for(i=0;i<Ke2->row_fsize;i++){
			if(P[Ke2->row_flocaldoflist[i]] >= Ke1->nrows){
				this->row_flocaldoflist[counter] =P[Ke2->row_flocaldoflist[i]];
				this->row_fglobaldoflist[counter]=Ke2->row_fglobaldoflist[i];
				counter++;
			}
		}
	}
	else{
		this->row_flocaldoflist=NULL;
		this->row_fglobaldoflist=NULL;
	}

	/*Sset*/
	ssize=Ke1->row_ssize;
	for(i=0;i<Ke2->row_ssize;i++){
		if(P[Ke2->row_slocaldoflist[i]] >= Ke1->nrows) ssize++;
	}
	this->row_ssize=ssize;
	if(ssize){
		this->row_slocaldoflist =xNew<int>(ssize);
		this->row_sglobaldoflist=xNew<int>(ssize);
		for(i=0;i<Ke1->row_ssize;i++){
			this->row_slocaldoflist[i] =Ke1->row_slocaldoflist[i];
			this->row_sglobaldoflist[i]=Ke1->row_sglobaldoflist[i];
		}
		counter=Ke1->row_ssize;
		for(i=0;i<Ke2->row_ssize;i++){
			if(P[Ke2->row_slocaldoflist[i]] >= Ke1->nrows){
				this->row_slocaldoflist[counter] =P[Ke2->row_slocaldoflist[i]];
				this->row_sglobaldoflist[counter]=Ke2->row_sglobaldoflist[i];
				counter++;
			}
		}
	}
	else{
		this->row_slocaldoflist=NULL;
		this->row_sglobaldoflist=NULL;
	}

	/*don't do cols, we can pick them up from the rows: */
	this->col_fsize=0;
	this->col_flocaldoflist=NULL;
	this->col_fglobaldoflist=NULL;
	this->col_ssize=0;
	this->col_slocaldoflist=NULL;
	this->col_sglobaldoflist=NULL;

	/*clean-up*/
	xDelete<int>(P);
}
/*}}}*/
ElementMatrix::ElementMatrix(ElementMatrix* Ke1, ElementMatrix* Ke2,ElementMatrix* Ke3){/*{{{*/

	/*Concatenate all matrices*/
	ElementMatrix* Ke12 =new ElementMatrix(Ke1,Ke2);
	ElementMatrix* Ke123=new ElementMatrix(Ke12,Ke3);

	/*Initialize current object with this matrix*/
	this->Init(Ke123);

	/*clean-up*/
	delete Ke12;
	delete Ke123;
}
/*}}}*/
ElementMatrix::ElementMatrix(Node** nodes,int numnodes,Parameters* parameters,int approximation){/*{{{*/

	/*get Matrix size and properties*/
	this->dofsymmetrical=true;
	this->nrows=GetNumberOfDofs(nodes,numnodes,GsetEnum,approximation);
	this->ncols=this->nrows;

	/*fill values with 0: */
	this->values=xNewZeroInit<IssmDouble>(this->nrows*this->ncols);

	/*g list*/
	this->gglobaldoflist=GetGlobalDofList(nodes,numnodes,GsetEnum,approximation);

	/*get dof lists for f and s set: */
	this->row_fsize=GetNumberOfDofs(nodes,numnodes,FsetEnum,approximation);
	this->row_flocaldoflist =GetLocalDofList( nodes,numnodes,FsetEnum,approximation);
	this->row_fglobaldoflist=GetGlobalDofList(nodes,numnodes,FsetEnum,approximation);
	this->row_ssize=GetNumberOfDofs(nodes,numnodes,SsetEnum,approximation);
	this->row_slocaldoflist =GetLocalDofList( nodes,numnodes,SsetEnum,approximation);
	this->row_sglobaldoflist=GetGlobalDofList(nodes,numnodes,SsetEnum,approximation);

	/*Because this matrix is "dofsymmetrical" don't do cols, we can pick them up from the rows: */
	this->col_fsize=0;
	this->col_flocaldoflist=NULL;
	this->col_fglobaldoflist=NULL;
	this->col_ssize=0;
	this->col_slocaldoflist=NULL;
	this->col_sglobaldoflist=NULL;
}
/*}}}*/
ElementMatrix::~ElementMatrix(){/*{{{*/

	xDelete<IssmDouble>(this->values);
	xDelete<int>(this->gglobaldoflist);
	xDelete<int>(this->row_flocaldoflist);
	xDelete<int>(this->row_fglobaldoflist);
	xDelete<int>(this->row_slocaldoflist);
	xDelete<int>(this->row_sglobaldoflist);
	xDelete<int>(this->col_flocaldoflist);
	xDelete<int>(this->col_fglobaldoflist);
	xDelete<int>(this->col_slocaldoflist);
	xDelete<int>(this->col_sglobaldoflist);
}
/*}}}*/

/*ElementMatrix specific routines: */
void ElementMatrix::AddDiagonalToGlobal(Vector<IssmDouble>* pf){/*{{{*/

	IssmDouble* localvalues=NULL;

	/*Check that pf is not NULL*/
	_assert_(pf); 

	/*In debugging mode, check consistency (no NaN, and values not too big)*/
	this->CheckConsistency();

	if(this->dofsymmetrical){
		/*only use row dofs to add values into global matrices: */

		if(this->row_fsize){
			/*first, retrieve values that are in the f-set from the g-set values matrix: */
			localvalues=xNew<IssmDouble>(this->row_fsize);
			for(int i=0;i<this->row_fsize;i++){
				localvalues[i] = this->values[this->ncols*this->row_flocaldoflist[i]+ this->row_flocaldoflist[i]];
			}

			/*add local values into global  matrix, using the fglobaldoflist: */
			pf->SetValues(this->row_fsize,this->row_fglobaldoflist,localvalues,ADD_VAL);

			/*Free ressources:*/
			xDelete<IssmDouble>(localvalues);
		}
	}
	else{
		_error_("non dofsymmetrical matrix AddToGlobal routine not support yet!");
	}

}
/*}}}*/
void ElementMatrix::AddToGlobal(Matrix<IssmDouble>* Kff, Matrix<IssmDouble>* Kfs){/*{{{*/

	int i,j;
	IssmDouble* localvalues=NULL;

	/*Check that Kff has been alocated in debugging mode*/
	_assert_(Kff);

	/*If Kfs is not provided, call the other function*/
	if(!Kfs){
		this->AddToGlobal(Kff);
		return;
	}

	/*In debugging mode, check consistency (no NaN, and values not too big)*/
	this->CheckConsistency();

	if(this->dofsymmetrical){
		/*only use row dofs to add values into global matrices: */

		if(this->row_fsize){
			/*first, retrieve values that are in the f-set from the g-set values matrix: */
			localvalues=xNew<IssmDouble>(this->row_fsize*this->row_fsize);
			for(i=0;i<this->row_fsize;i++){
				for(j=0;j<this->row_fsize;j++){
					localvalues[this->row_fsize*i+j]=this->values[this->ncols*this->row_flocaldoflist[i]+this->row_flocaldoflist[j]];
				}
			}

			/*add local values into global  matrix, using the fglobaldoflist: */
			Kff->SetValues(this->row_fsize,this->row_fglobaldoflist,this->row_fsize,this->row_fglobaldoflist,localvalues,ADD_VAL);

			/*Free ressources:*/
			xDelete<IssmDouble>(localvalues);
		}

		if((this->row_ssize!=0) && (this->row_fsize!=0)){
			/*first, retrieve values that are in the f and s-set from the g-set values matrix: */
			localvalues=xNew<IssmDouble>(this->row_fsize*this->row_ssize);
			for(i=0;i<this->row_fsize;i++){
				for(j=0;j<this->row_ssize;j++){
					localvalues[this->row_ssize*i+j]=this->values[this->ncols*this->row_flocaldoflist[i]+this->row_slocaldoflist[j]];
				}
			}
			/*add local values into global  matrix, using the fglobaldoflist: */
			Kfs->SetValues(this->row_fsize,this->row_fglobaldoflist,this->row_ssize,this->row_sglobaldoflist,localvalues,ADD_VAL);

			/*Free ressources:*/
			xDelete<IssmDouble>(localvalues);
		}
	}
	else{
		_error_("non dofsymmetrical matrix AddToGlobal routine not support yet!");
	}

}
/*}}}*/
void ElementMatrix::AddToGlobal(Matrix<IssmDouble>* Jff){/*{{{*/

	int i,j;
	IssmDouble* localvalues=NULL;

	/*Check that Jff is not NULL*/
	_assert_(Jff); 

	/*In debugging mode, check consistency (no NaN, and values not too big)*/
	this->CheckConsistency();

	if(this->dofsymmetrical){
		/*only use row dofs to add values into global matrices: */

		if(this->row_fsize){
			/*first, retrieve values that are in the f-set from the g-set values matrix: */
			localvalues=xNew<IssmDouble>(this->row_fsize*this->row_fsize);
			for(i=0;i<this->row_fsize;i++){
				for(j=0;j<this->row_fsize;j++){
					*(localvalues+this->row_fsize*i+j)=*(this->values+this->ncols*this->row_flocaldoflist[i]+this->row_flocaldoflist[j]);
				}
			}
			/*add local values into global  matrix, using the fglobaldoflist: */
			Jff->SetValues(this->row_fsize,this->row_fglobaldoflist,this->row_fsize,this->row_fglobaldoflist,localvalues,ADD_VAL);

			/*Free ressources:*/
			xDelete<IssmDouble>(localvalues);
		}

	}
	else{
		_error_("non dofsymmetrical matrix AddToGlobal routine not support yet!");
	}

}
/*}}}*/
void ElementMatrix::CheckConsistency(void){/*{{{*/
	/*Check element matrix values, only in debugging mode*/
	#ifdef _ISSM_DEBUG_ 
	for (int i=0;i<this->nrows;i++){
		for(int j=0;j<this->ncols;j++){
			if (xIsNan<IssmDouble>(this->values[i*this->ncols+j])) _error_("NaN found in Element Matrix");
			if (xIsInf<IssmDouble>(this->values[i*this->ncols+j])) _error_("Inf found in Element Matrix");
			if (fabs(this->values[i*this->ncols+j])>1.e+50) _error_("Element Matrix values exceeds 1.e+50");
		}
	}
	#endif
}
/*}}}*/
void ElementMatrix::Echo(void){/*{{{*/

	int i,j;
	_printf_("Element Matrix echo:\n");
	_printf_("   nrows: " << this->nrows << "\n");
	_printf_("   ncols: " << this->ncols << "\n");
	_printf_("   dofsymmetrical: " << (dofsymmetrical?"true":"false") << "\n");

	_printf_("   values:\n");
	for(i=0;i<nrows;i++){
		_printf_(setw(4) << right << i << ": ");
		for(j=0;j<ncols;j++) _printf_( " " << setw(11) << setprecision (5) << right << values[i*ncols+j]);
		_printf_("\n");
	}

	_printf_("   gglobaldoflist (" << gglobaldoflist << "): ");
	if(gglobaldoflist) for(i=0;i<nrows;i++) _printf_(" " << gglobaldoflist[i]); _printf_("\n");

	_printf_("   row_fsize: " << row_fsize << "\n");
	_printf_("   row_flocaldoflist  (" << row_flocaldoflist << "): ");
	if(row_flocaldoflist) for(i=0;i<row_fsize;i++) _printf_(" " << row_flocaldoflist[i] ); _printf_(" \n");
	_printf_("   row_fglobaldoflist  (" << row_fglobaldoflist << "): ");
	if(row_fglobaldoflist) for(i=0;i<row_fsize;i++)_printf_(" " << row_fglobaldoflist[i]); _printf_(" \n");

	_printf_("   row_ssize: " << row_ssize << "\n");
	_printf_("   row_slocaldoflist  (" << row_slocaldoflist << "): ");
	if(row_slocaldoflist) for(i=0;i<row_ssize;i++) _printf_(" " << row_slocaldoflist[i] ); _printf_(" \n");
	_printf_("   row_sglobaldoflist  (" << row_sglobaldoflist << "): ");
	if(row_sglobaldoflist) for(i=0;i<row_ssize;i++)_printf_(" " << row_sglobaldoflist[i]); _printf_(" \n");

	if(!dofsymmetrical){
		_printf_("   col_fsize: " << col_fsize << "\n");
		_printf_("   col_flocaldoflist  (" << col_flocaldoflist << "): ");
		if(col_flocaldoflist) for(i=0;i<col_fsize;i++) _printf_(" " << col_flocaldoflist[i] ); _printf_(" \n");
		_printf_("   col_fglobaldoflist  (" << col_fglobaldoflist << "): ");
		if(col_fglobaldoflist) for(i=0;i<col_fsize;i++)_printf_(" " << col_fglobaldoflist[i]); _printf_(" \n");

		_printf_("   col_ssize: " << col_ssize << "\n");
		_printf_("   col_slocaldoflist  (" << col_slocaldoflist << "): ");
		if(col_slocaldoflist) for(i=0;i<col_ssize;i++) _printf_(" " << col_slocaldoflist[i] ); _printf_(" \n");
		_printf_("   col_sglobaldoflist  (" << col_sglobaldoflist << "): ");
		if(col_sglobaldoflist) for(i=0;i<col_ssize;i++)_printf_(" " << col_sglobaldoflist[i]); _printf_(" \n");
	}
}
/*}}}*/
void ElementMatrix::Init(ElementMatrix* Ke){/*{{{*/

	_assert_(Ke);
	_assert_(this);

	this->nrows =Ke->nrows;
	this->ncols =Ke->ncols;
	this->dofsymmetrical=Ke->dofsymmetrical;

	this->values=xNew<IssmDouble>(this->nrows*this->ncols);
	xMemCpy<IssmDouble>(this->values,Ke->values,this->nrows*this->ncols);

	this->gglobaldoflist=xNew<int>(this->nrows);
	xMemCpy<int>(this->gglobaldoflist,Ke->gglobaldoflist,this->nrows);

	this->row_fsize=Ke->row_fsize;
	if(this->row_fsize){
		this->row_flocaldoflist=xNew<int>(this->row_fsize);
		xMemCpy<int>(this->row_flocaldoflist,Ke->row_flocaldoflist,this->row_fsize);
		this->row_fglobaldoflist=xNew<int>(this->row_fsize);
		xMemCpy<int>(this->row_fglobaldoflist,Ke->row_fglobaldoflist,this->row_fsize);
	}
	else{
		this->row_flocaldoflist=NULL;
		this->row_fglobaldoflist=NULL;
	}

	this->row_ssize=Ke->row_ssize;
	if(this->row_ssize){
		this->row_slocaldoflist=xNew<int>(this->row_ssize);
		xMemCpy<int>(this->row_slocaldoflist,Ke->row_slocaldoflist,this->row_ssize);
		this->row_sglobaldoflist=xNew<int>(this->row_ssize);
		xMemCpy<int>(this->row_sglobaldoflist,Ke->row_sglobaldoflist,this->row_ssize);
	}
	else{
		this->row_slocaldoflist=NULL;
		this->row_sglobaldoflist=NULL;
	}

	this->col_fsize=Ke->col_fsize;
	if(this->col_fsize){
		this->col_flocaldoflist=xNew<int>(this->col_fsize);
		xMemCpy<int>(this->col_flocaldoflist,Ke->col_flocaldoflist,this->col_fsize);
		this->col_fglobaldoflist=xNew<int>(this->col_fsize);
		xMemCpy<int>(this->col_fglobaldoflist,Ke->col_fglobaldoflist,this->col_fsize);
	}
	else{
		this->col_flocaldoflist=NULL;
		this->col_fglobaldoflist=NULL;
	}

	this->col_ssize=Ke->col_ssize;
	if(this->col_ssize){
		this->col_slocaldoflist=xNew<int>(this->col_ssize);
		xMemCpy<int>(this->col_slocaldoflist,Ke->col_slocaldoflist,this->col_ssize);
		this->col_sglobaldoflist=xNew<int>(this->col_ssize);
		xMemCpy<int>(this->col_sglobaldoflist,Ke->col_sglobaldoflist,this->col_ssize);
	}
	else{
		this->col_slocaldoflist=NULL;
		this->col_sglobaldoflist=NULL;
	}
}
/*}}}*/
void ElementMatrix::Lump(void){/*{{{*/

	if(!dofsymmetrical) _error_("not supported yet");

	for(int i=0;i<this->nrows;i++){
		for(int j=0;j<this->ncols;j++){
			if(i!=j){
				this->values[i*this->ncols+i] += this->values[i*this->ncols+j];
				this->values[i*this->ncols+j]  = 0.;
			}
		}
	}

	return;
}
/*}}}*/
void ElementMatrix::Transpose(void){/*{{{*/

	/*Intermediaries*/
	ElementMatrix* Ke_copy=new ElementMatrix(this);

	/*Update sizes*/
	this->nrows=Ke_copy->ncols;
	this->ncols=Ke_copy->nrows;

	/*Transpose values*/
	for (int i=0;i<this->nrows;i++) for(int j=0;j<this->ncols;j++) this->values[i*this->ncols+j]=Ke_copy->values[j*Ke_copy->ncols+i];

	/*Transpose indices*/
	if(!dofsymmetrical){
		_error_("not supported yet");
	}

	/*Clean up and return*/
	delete Ke_copy;
	return;
}
/*}}}*/
void ElementMatrix::StaticCondensation(int bsize,int* bindices){/*{{{*/
	/* 
	 * | Kii  Kib | | Ui |    |Fi|
	 * | Kbi  Kbb | | Ub |  = |Fb|
	 *
	 * Kii Ui + Kib Ub = Fi
	 * Kbi Ui + Kbb Ub = Fb
	 *
	 * We want to remove Ub from the equation:
	 *
	 * Kii Ui + Kib inv(Kbb) (Fb - Kbi Ui) = Fi
	 *
	 * which gives:
	 *
	 * (Kii - Kib inv(Kbb) Kbi) Ui = Fi - Kib inv(Kbb) Fb
	 */

	/*Checks in debugging mode*/
	_assert_(this->nrows==this->ncols && bsize>0 && bsize<this->ncols && this->values); 

	/*Intermediaries*/
	int         counter,i,j,isize;
	IssmDouble *Kii         = NULL;
	IssmDouble *Kib         = NULL;
	IssmDouble *Kbi         = NULL;
	IssmDouble *Kbb         = NULL;
	IssmDouble *Kbbinv      = NULL;
	IssmDouble *Ktemp       = NULL;
	int        *iindices    = NULL;
	bool        flag;

	/*Get new sizes and indices*/
	isize    = this->nrows - bsize;
	iindices = xNew<int>(isize);
	counter  = 0;
	for(i=0;i<this->nrows;i++){
		flag = true;
		for(j=0;j<bsize;j++){
			if(i==bindices[j]){
				flag = false;
				break;
			}
		}
		if(flag){
			_assert_(counter<isize);
			iindices[counter++] = i;
		}
	}
	_assert_(counter == isize);

	/*Get submatrices*/
	Kii = xNew<IssmDouble>(isize*isize);
	Kib = xNew<IssmDouble>(isize*bsize);
	Kbi = xNew<IssmDouble>(bsize*isize);
	Kbb = xNew<IssmDouble>(bsize*bsize);
	for(i=0;i<isize;i++) for(j=0;j<isize;j++) Kii[i*isize+j] = this->values[iindices[i]*this->ncols + iindices[j]];
	for(i=0;i<isize;i++) for(j=0;j<bsize;j++) Kib[i*bsize+j] = this->values[iindices[i]*this->ncols + bindices[j]];
	for(i=0;i<bsize;i++) for(j=0;j<isize;j++) Kbi[i*isize+j] = this->values[bindices[i]*this->ncols + iindices[j]];
	for(i=0;i<bsize;i++) for(j=0;j<bsize;j++) Kbb[i*bsize+j] = this->values[bindices[i]*this->ncols + bindices[j]];

	/*Invert Kbb*/
	Kbbinv = xNew<IssmDouble>(bsize*bsize);
	switch(bsize){
		case 1:
			Kbbinv[0] = 1./Kbb[0];
			break;
		case 2:
			Matrix2x2Invert(Kbbinv,Kbb);
			break;
		case 3:
			Matrix3x3Invert(Kbbinv,Kbb);
			break;
		default:
			MatrixInverse(Kbbinv,bsize,bsize,NULL,0,NULL);
			break;
	}

	/*Calculate  Kib inv(Kbb) Kbi*/
	Ktemp = xNew<IssmDouble>(isize*isize);
	TripleMultiply(Kib,isize,bsize,0, Kbbinv,bsize,bsize,0, Kbi,bsize,isize,0, Ktemp,0);

	/*New Ke*/
	for(i=0;i<isize*isize;i++) Ktemp[i] = Kii[i] - Ktemp[i];

	/*Update matrix values*/
	for(i=0;i<this->nrows*this->ncols;i++) this->values[i]=0.;
	for(i=0;i<isize;i++){
		for(j=0;j<isize;j++){
			this->values[iindices[i]*this->ncols + iindices[j]] = Ktemp[i*isize+j];
		}
	}

	/*Clean up and return*/
	xDelete<IssmDouble>(Kii);
	xDelete<IssmDouble>(Kib);
	xDelete<IssmDouble>(Kbi);
	xDelete<IssmDouble>(Kbb);
	xDelete<IssmDouble>(Kbbinv);
	xDelete<IssmDouble>(Ktemp);
	xDelete<int>(iindices);
	return;
}
/*}}}*/
