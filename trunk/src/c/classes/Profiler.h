/*!\file Profiler.h
 * \brief: header file for node object
 */

#ifndef _PROFILER_H_
#define _PROFILER_H_

/*Headers:*/
/*{{{*/
#include "../datastructures/datastructures.h"
#include "../shared/shared.h"
/*}}}*/

class DoubleParam;

#define START 0
#define STARTINIT 1
#define FINISHINIT 2
#define STARTCORE 3
#define FINISHCORE 4
#define STARTADCORE 5
#define FINISHADCORE 6
#define FINISH 7
#define MAXIMUMSIZE 8 

class Profiler: public Object{

	public: 
		IssmDouble flops[MAXIMUMSIZE];
		IssmDouble memory[MAXIMUMSIZE];
		IssmDouble time[MAXIMUMSIZE];

		/*Profiler constructors, destructors {{{*/
		Profiler();
		~Profiler();
		/*}}}*/
		/*Object virtual functions definitions:{{{ */
		Object *copy();
		void    DeepEcho();
		void    Echo();
		int     Id();
		void Marshall(char** pmarshalled_data,int* pmarshalled_data_size, int marshall_direction);
		int     ObjectEnum();
		/*}}}*/
		/*Profiler routines {{{*/
		IssmDouble  DeltaFlops(int inittag, int finaltag);
		IssmDouble  DeltaTime(int inittag, int finaltag);
		int     DeltaTimeModHour(int inittag, int finaltag);
		int     DeltaTimeModMin(int inittag, int finaltag);
		int     DeltaTimeModSec(int inittag, int finaltag);
		IssmDouble  Memory(int tag);
		void    Tag(int tagenum,bool dontmpisync=false);
		/*}}}*/
};

#endif  /* _PROFILER_H_ */
