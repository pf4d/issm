/*!\file DofIndexing.c
 * \brief: implementation of the DofIndexing object
 */

#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#include <string.h>

#include "./DofIndexing.h"
#include "../shared/Numerics/types.h"
#include "../shared/Numerics/constants.h"
#include "../shared/io/Print/Print.h"
#include "../shared/io/Marshalling/Marshalling.h"
#include "../shared/Exceptions/exceptions.h"
#include "../shared/MemOps/MemOps.h"
#include "../shared/Enum/Enum.h"

/*DofIndexing constructors and destructor*/
DofIndexing::DofIndexing(){/*{{{*/

	this->gsize    = UNDEF;
	this->fsize    = UNDEF;
	this->ssize    = UNDEF;
	this->clone    = false;
	this->active   = true;
	this->freeze   = false;
	this->f_set    = NULL;
	this->s_set    = NULL;
	this->svalues  = NULL;
	this->doftype  = NULL;
	this->gdoflist = NULL;
	this->fdoflist = NULL;
	this->sdoflist = NULL;

}
/*}}}*/
DofIndexing::DofIndexing(int in_gsize){/*{{{*/

	this->Init(in_gsize,NULL);

}
/*}}}*/
DofIndexing::DofIndexing(DofIndexing* in){ //copy constructor/*{{{*/

	this->gsize  = in->gsize;
	this->fsize  = in->fsize;
	this->ssize  = in->ssize;
	this->clone  = in->clone;
	this->active = in->active;
	this->freeze = in->freeze;

	if(this->gsize>0){
		this->f_set=xNew<bool>(this->gsize);
		this->s_set=xNew<bool>(this->gsize);
		this->svalues=xNew<IssmDouble>(this->gsize);
		if(in->doftype){
			this->doftype=xNew<int>(this->gsize); 
		}
		else{
			this->doftype=NULL;
		}
		this->gdoflist=xNew<int>(this->gsize); 
	}
	else{
		this->f_set    = NULL;
		this->s_set    = NULL;
		this->svalues  = NULL;
		this->doftype  = NULL;
		this->gdoflist = NULL;
	}
	if(this->fsize>0)this->fdoflist=xNew<int>(this->fsize); else this->fdoflist=NULL;
	if(this->ssize>0)this->sdoflist=xNew<int>(this->ssize); else this->sdoflist=NULL;

	if(this->gsize>0){
		memcpy(this->f_set,in->f_set,this->gsize*sizeof(bool));
		memcpy(this->s_set,in->s_set,this->gsize*sizeof(bool));
		xMemCpy<IssmDouble>(this->svalues,in->svalues,this->gsize);
		if(this->doftype)memcpy(this->doftype,in->doftype,this->gsize*sizeof(int));
		memcpy(this->gdoflist,in->gdoflist,this->gsize*sizeof(int));
	}
	if(this->fsize>0)memcpy(this->fdoflist,in->fdoflist,this->fsize*sizeof(int));
	if(this->ssize>0)memcpy(this->sdoflist,in->sdoflist,this->ssize*sizeof(int));

}
/*}}}*/
DofIndexing::~DofIndexing(){ //destructor/*{{{*/

	if(this->f_set) xDelete<bool>(f_set); 
	if(this->s_set) xDelete<bool>(s_set); 
	if(this->svalues) xDelete<IssmDouble>(svalues);
	if(this->doftype) xDelete<int>(doftype); 
	if(this->gdoflist) xDelete<int>(gdoflist);
	if(this->fdoflist) xDelete<int>(fdoflist);
	if(this->sdoflist) xDelete<int>(sdoflist);

}
/*}}}*/
DofIndexing DofIndexing::operator=( const DofIndexing& in ){/*{{{*/

	this->copy(in);

	return this;
}
/*}}}*/
void DofIndexing::copy(const DofIndexing& in ){/*{{{*/

	this->gsize  = in.gsize;
	this->fsize  = in.fsize;
	this->ssize  = in.ssize;
	this->clone  = in.clone;
	this->active = in.active;
	this->freeze = in.freeze;

	if(this->gsize>0){
		this->f_set=xNew<bool>(this->gsize);
		this->s_set=xNew<bool>(this->gsize);
		this->svalues=xNew<IssmDouble>(this->gsize);
		if(in.doftype){
			this->doftype=xNew<int>(this->gsize);
		}
		else{
			this->doftype=NULL;
		}
		this->gdoflist=xNew<int>(this->gsize);
	}
	else{
		this->f_set    = NULL;
		this->s_set    = NULL;
		this->svalues  = NULL;
		this->doftype  = NULL;
		this->gdoflist = NULL;
	}
	if(this->fsize>0)this->fdoflist=xNew<int>(this->fsize); else this->fdoflist=NULL;
	if(this->ssize>0)this->sdoflist=xNew<int>(this->ssize); else this->sdoflist=NULL;

	if(this->gsize>0){
		memcpy(this->f_set,in.f_set,this->gsize*sizeof(bool));
		memcpy(this->s_set,in.s_set,this->gsize*sizeof(bool));
		xMemCpy<IssmDouble>(this->svalues,in.svalues,this->gsize);
		if(this->doftype)memcpy(this->doftype,in.doftype,this->gsize*sizeof(int));
		memcpy(this->gdoflist,in.gdoflist,this->gsize*sizeof(int));
	}
	if(this->fsize>0)memcpy(this->fdoflist,in.fdoflist,this->fsize*sizeof(int));
	if(this->ssize>0)memcpy(this->sdoflist,in.sdoflist,this->ssize*sizeof(int));

	return;
}
/*}}}*/
void DofIndexing::Init(int in_gsize,int* in_doftype){/*{{{*/

	this->gsize = in_gsize;

	/*At this point, assume this is not a clone (will be dealt with later)*/
	this->clone = false;

	/*memory allocation */
	if(this->gsize>0){
		this->f_set    = xNew<bool>((unsigned int)in_gsize);
		this->s_set    = xNew<bool>((unsigned int)in_gsize);
		this->svalues  = xNew<IssmDouble>((unsigned int)in_gsize);
		this->gdoflist = xNew<int>((unsigned int)in_gsize);

		if(in_doftype) this->doftype = xNew<int>((unsigned int)in_gsize);
		else this->doftype = NULL;
	}

	/*Assign values assuming no Dirichlet at this point*/
	for(int i=0;i<this->gsize;i++){
		this->f_set[i]    = true;
		this->s_set[i]    = false;
		this->svalues[i]  = 0.;      //0 constraint is the default value
		this->gdoflist[i] = UNDEF;

		if(this->doftype) this->doftype[i]=in_doftype[i];
	}
}
/*}}}*/
void DofIndexing::InitSet(int setenum){/*{{{*/

	int i;
	int size=0;

	/*go through sets, and figure out how many dofs belong to this set, except for g-set, 
	 * which has already been initialized: */
	if(setenum==FsetEnum){
		size=0;
		for(i=0;i<this->gsize;i++) if(f_set[i])size++;
		this->fsize=size;
		xDelete<int>(this->fdoflist);

		if(this->fsize)
		 this->fdoflist=xNew<int>(size);
		else
		 this->fdoflist=NULL;
	}
	else if(setenum==SsetEnum){
		size=0;
		for(i=0;i<this->gsize;i++) if(s_set[i])size++;
		this->ssize=size;
		xDelete<int>(this->sdoflist);

		if(this->ssize)
		 this->sdoflist=xNew<int>(size);
		else
		 this->sdoflist=NULL;
	}
	else _error_("set of enum type " << EnumToStringx(setenum) << " not supported yet!");
}
/*}}}*/

/*Some of the Object functionality: */
void DofIndexing::Activate(void){/*{{{*/

	this->active = true;

	/*Constrain to 0. at this point*/
	for(int i=0;i<this->gsize;i++){
		this->f_set[i]    = true;
		this->s_set[i]    = false;
		this->svalues[i]  = 0.; 
	}
	return;
}
/*}}}*/
void DofIndexing::Deactivate(void){/*{{{*/
	this->active = false;

	/*Constrain to 0. at this point*/
	for(int i=0;i<this->gsize;i++){
		this->f_set[i]    = false;
		this->s_set[i]    = true;
		this->svalues[i]  = 0.; 
	}
	return;
}
/*}}}*/
void DofIndexing::DeepEcho(void){/*{{{*/

	int i;

	_printf_("DofIndexing:\n");
	_printf_("   gsize:  " << gsize << "\n");
	_printf_("   fsize:  " << fsize << "\n");
	_printf_("   ssize:  " << ssize << "\n");
	_printf_("   clone:  " << clone << "\n");
	_printf_("   active: " << active << "\n");
	_printf_("   freeze: " << freeze << "\n");

	_printf_("   f_set = [ ");
	for(i=0;i<gsize;i++) _printf_((f_set[i]?1:0)<< " ");
	_printf_("]\n");
	_printf_("   s_set = [ ");
	for(i=0;i<gsize;i++) _printf_((s_set[i]?1:0)<< " ");
	_printf_("]\n");

	_printf_("   svalues (" << this->ssize << "): |");
	for(i=0;i<this->gsize;i++){
		if(this->s_set[i])_printf_(" " << svalues[i] << " |");
	}
	_printf_("\n");

	if(doftype){
		_printf_("   doftype: |");
		for(i=0;i<gsize;i++){
			_printf_(" " << doftype[i] << " |");
		}
		_printf_("\n");
	}
	else _printf_("   doftype: NULL\n");

	_printf_("   g_doflist (" << this->gsize << "): |");
	for(i=0;i<this->gsize;i++){
		_printf_(" " << gdoflist[i] << " |");
	}
	_printf_("\n");

	_printf_("   f_doflist (" << this->fsize << "): |");
	for(i=0;i<this->fsize;i++){
		_printf_(" " << fdoflist[i] << " |");
	}
	_printf_("\n");

	_printf_("   s_doflist (" << this->ssize << "): |");
	for(i=0;i<this->ssize;i++){
		_printf_(" " << sdoflist[i] << " |");
	}
	_printf_("\n");
}		
/*}}}*/
void DofIndexing::Echo(void){/*{{{*/

	_printf_("DofIndexing:\n");
	_printf_("   gsize:  " << gsize << "\n");
	_printf_("   fsize:  " << fsize << "\n");
	_printf_("   ssize:  " << ssize << "\n");
	_printf_("   clone:  " << clone << "\n");
	_printf_("   active: " << active << "\n");
	_printf_("   freeze: " << freeze << "\n");
}
/*}}}*/
void DofIndexing::Marshall(char** pmarshalled_data,int* pmarshalled_data_size, int marshall_direction){ /*{{{*/

	MARSHALLING(gsize);
	MARSHALLING(fsize);
	MARSHALLING(ssize);
	MARSHALLING(clone);
	MARSHALLING(active);
	MARSHALLING(freeze);
	MARSHALLING_DYNAMIC(f_set,bool,gsize);
	MARSHALLING_DYNAMIC(s_set,bool,gsize);
	MARSHALLING_DYNAMIC(svalues,IssmDouble,gsize);
	MARSHALLING_DYNAMIC(doftype,int,gsize);
	MARSHALLING_DYNAMIC(gdoflist,int,gsize);
	MARSHALLING_DYNAMIC(fdoflist,int,fsize);
	MARSHALLING_DYNAMIC(sdoflist,int,ssize);

}
/*}}}*/
