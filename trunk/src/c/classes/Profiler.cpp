/*!\file Profiler.c
 * \brief: implementation of the Profiler object
 */

/*Include files: {{{*/
#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#include "./Profiler.h"
#include "./Params/DoubleParam.h"
#include "../toolkits/toolkits.h"
/*}}}*/

/*Profiler constructors and destructors:*/
Profiler::Profiler(){/*{{{*/
	for(int i=0;i<MAXIMUMSIZE;i++){
		this->time[i]  =NAN;
		this->flops[i] =NAN;
		this->memory[i]=NAN;
	}
} /*}}}*/
Profiler::~Profiler(){/*{{{*/
	/*Nothing to delete, everything is statically allocated*/
} /*}}}*/
Object* Profiler::copy(){/*{{{*/
	/*First do simple copy: */
	Profiler* output=new Profiler();

	for(int i=0;i<MAXIMUMSIZE;i++){
		output->time[i]  =this->time[i];
		output->flops[i] =this->flops[i];
		output->memory[i]=this->memory[i];
	}

	return (Object*)output;
}
/*}}}*/

/*Object virtual functions definitions:*/
void Profiler::DeepEcho(void){/*{{{*/

	this->Echo();

}
/*}}}*/
void Profiler::Echo(void){/*{{{*/

	_printf_("Profiler:\n");
	for(int i=0;i<MAXIMUMSIZE;i++){
		_printf_("    Tag "<<i<<":\n");
		_printf_("       flops:  "<<this->flops[i]<<"\n");
		_printf_("       memory: "<<this->memory[i]<<"\n");
		_printf_("       time:   "<<this->time[i]<<"\n");
	}

}
/*}}}*/
int  Profiler::Id(void){ return -1; }/*{{{*/
/*}}}*/
void Profiler::Marshall(char** pmarshalled_data,int* pmarshalled_data_size, int marshall_direction){ /*{{{*/

	IssmDouble* pointer = NULL;

	MARSHALLING_ENUM(ProfilerEnum);
	pointer = &this->time[0];
	MARSHALLING_DYNAMIC(pointer,IssmDouble,MAXIMUMSIZE);
	pointer = &this->flops[0];
	MARSHALLING_DYNAMIC(pointer,IssmDouble,MAXIMUMSIZE);
	pointer = &this->memory[0];
	MARSHALLING_DYNAMIC(pointer,IssmDouble,MAXIMUMSIZE);

} /*}}}*/
int  Profiler::ObjectEnum(void){/*{{{*/

	return ProfilerEnum;

}
/*}}}*/

/*Profiler routines:*/
IssmDouble  Profiler::DeltaFlops(int inittag, int finaltag){/*{{{*/

	/*Get initial flops*/
	_assert_(inittag>=0); 
	_assert_(inittag<MAXIMUMSIZE); 
	if(xIsNan<IssmDouble>(this->flops[inittag])) _error_("Tag not set");
	IssmDouble init = this->flops[inittag];

	/*Get final flops*/
	_assert_(finaltag>=0); 
	_assert_(finaltag<MAXIMUMSIZE); 
	if(xIsNan<IssmDouble>(this->flops[finaltag])) _error_("Tag not set");
	IssmDouble final = this->flops[finaltag];

	return final-init;
}
/*}}}*/
IssmDouble  Profiler::DeltaTime(int inittag, int finaltag){/*{{{*/

	/*Get initial time*/
	_assert_(inittag>=0); 
	_assert_(inittag<MAXIMUMSIZE); 
	if(xIsNan<IssmDouble>(this->time[inittag])) _error_("Tag "<<inittag<<" not set");
	IssmDouble init = this->time[inittag];

	/*Get final time*/
	_assert_(finaltag>=0); 
	_assert_(finaltag<MAXIMUMSIZE); 
	if(xIsNan<IssmDouble>(this->time[finaltag])) _error_("Tag "<<finaltag<<" not set");
	IssmDouble final = this->time[finaltag];

	#ifdef _HAVE_MPI_
	return final-init;
	#else
	return (final-init)/CLOCKS_PER_SEC;
	#endif
}
/*}}}*/
int Profiler::DeltaTimeModHour(int inittag, int finishtag){/*{{{*/

	IssmDouble delta = this->DeltaTime(inittag,finishtag);
	return int((reCast<int,IssmDouble>(delta))/3600);

}
/*}}}*/
int Profiler::DeltaTimeModMin(int inittag, int finishtag){/*{{{*/

	IssmDouble delta = this->DeltaTime(inittag,finishtag);
	return int(int(reCast<int,IssmDouble>(delta))%3600/60);
}
/*}}}*/
int Profiler::DeltaTimeModSec(int inittag, int finishtag){/*{{{*/

	IssmDouble delta = this->DeltaTime(inittag,finishtag);
	return int(reCast<int,IssmDouble>(delta)%60);
}
/*}}}*/
IssmDouble  Profiler::Memory(int tag){/*{{{*/

	/*Get initial flops*/
	_assert_(tag>=0); 
	_assert_(tag<MAXIMUMSIZE); 
	if(xIsNan<IssmDouble>(this->flops[tag])) _error_("Tag not set");
	return this->memory[tag];
}
/*}}}*/
void  Profiler::Tag(int tagenum,bool dontmpisync){/*{{{*/

	IssmDouble t;
	IssmDouble f;
	IssmDouble m;

	/*If mpisync requested, make sure all the cpus are at the same point 
	 *in the execution: */
	if(!dontmpisync){
		ISSM_MPI_Barrier(IssmComm::GetComm()); 
	}

	/*Capture time: */
	#ifdef _HAVE_MPI_
	t=ISSM_MPI_Wtime();
	#else
	t=(IssmPDouble)clock();
	#endif

	/*Capture flops: */
	#ifdef _HAVE_PETSC_
		PetscGetFlops(&f);
		PetscMemoryGetCurrentUsage(&m);
	#else
		/*do nothing for now:*/
	#endif

	/*Plug into this->time: */
	_assert_(tagenum>=0); 
	_assert_(tagenum<MAXIMUMSIZE); 
	if(!xIsNan<IssmDouble>(this->time[tagenum])) _error_("Tag already exists");
	this->time[tagenum]  = t;
	if(!xIsNan<IssmDouble>(this->flops[tagenum])) _error_("Tag already exists");
	this->flops[tagenum] = f;
	if(!xIsNan<IssmDouble>(this->memory[tagenum])) _error_("Tag already exists");
	this->memory[tagenum]= m;

}
/*}}}*/
