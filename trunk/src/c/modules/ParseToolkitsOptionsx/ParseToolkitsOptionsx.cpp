/*!\file ParseToolkitsOptionsx
 * * \brief: parse options present in a petsc file, and create petsc options 
 * objects accordingly. This will be used to drive the behaviour of Toolkits for 
 * each analysis type.
 */

#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#include <cstring>

#include "./ParseToolkitsOptionsx.h"
#include "../../shared/shared.h"
#include "../../toolkits/toolkits.h"

void ParseToolkitsOptionsx(Parameters* parameters,FILE* fid){

	char line [1000];
	int my_rank;
	int i;

	/*intermediary: */
	IssmDouble* analyses=NULL;
	char** strings=NULL;
	int numanalyses;
	char* string=NULL;
	char* newstring=NULL;
	char* catstring=NULL;
	int   stringlength;

	/*Get my_rank:*/
	my_rank=IssmComm::GetRank();

	if(my_rank==0){

		/*Now, go through lines and figure out how many analyses we have: */
		numanalyses=0;
		while ( fgets(line, sizeof line, fid) ){
			/*skip comments and empty lines: */
			if ((line[0]=='%') || (line[0]=='\n') || (line[0]==' ') || (line[0]=='\t') || (line[0]=='\r'))continue;
			/*ok, first time, we should get an analysis enum, starting with a +: */
			if (line[0]=='+')numanalyses++;
			else continue;
		}

		/*Now, allocate analyses and strings: */
		analyses=xNew<IssmDouble>(numanalyses);
		strings=xNew<char*>(numanalyses);
		for(i=0;i<numanalyses;i++)strings[i]=NULL; 

		/*Go back to beginning of file:*/
		fseek(fid,0,SEEK_SET);
		numanalyses=0;
		while ( fgets(line, sizeof line, fid) ){

			/*skip comments and empty lines: */
			if ((line[0]=='%') || (line[0]=='\n') || (line[0]==' ') || (line[0]=='\t') || (line[0]=='\r'))continue;

			/*Get rid of end of line: */
			line[strlen(line)-1]='\0';

			if (line[0]=='+'){ /*this is the analysis line: */
				analyses[numanalyses]=StringToEnumx(&line[1]);  //skip the '+'
				numanalyses++;
			}
			else{ /*this is an option corresponding to analysis numanalyses-1. Add it 
			to the already existing options*/
				if(strings[numanalyses-1]==NULL){
					string=xNew<char>((strlen(line)+1)); 
					xMemCpy<char>(string,line,(strlen(line)+1));

					strings[numanalyses-1]=string;
				}
				else{
					string=strings[numanalyses-1];
					newstring=xNew<char>((strlen(line)+1));
					xMemCpy<char>(newstring,line,(strlen(line)+1));

					/*concatenate:*/
					catstring=xNew<char>(strlen(string)+1+strlen(newstring)+1+1); //fit in a space " "
					xMemCpy<char>(catstring,string,(strlen(string)+1));

					strcat(catstring," ");
					strcat(catstring,newstring);
					strings[numanalyses-1]=catstring;
					xDelete<char>(newstring);
					xDelete<char>(string);
				}
			}
		}
	}

	/*Ok, broadcast to other cpus: */
	ISSM_MPI_Bcast(&numanalyses,1,ISSM_MPI_INT,0,IssmComm::GetComm());
	if(my_rank!=0){
		analyses=xNew<IssmDouble>(numanalyses);
		strings=xNew<char*>(numanalyses);
	}
	ISSM_MPI_Bcast(analyses,numanalyses,ISSM_MPI_DOUBLE,0,IssmComm::GetComm());
	for(i=0;i<numanalyses;i++){
		char* string=strings[i];
		if(my_rank==0){
			if(string==NULL) _error_("PETSc options for analysis " << EnumToStringx(reCast<int>(analyses[i])) << " have been declared but were not found");
		}
		if(my_rank==0)stringlength=(strlen(string)+1)*sizeof(char);
		ISSM_MPI_Bcast(&stringlength,1,ISSM_MPI_INT,0,IssmComm::GetComm());
		if(my_rank!=0)string=xNew<char>(stringlength);
		ISSM_MPI_Bcast(string,stringlength,ISSM_MPI_CHAR,0,IssmComm::GetComm());
		if(my_rank!=0)strings[i]=string;
	}

	/*Ok, out of strings and analyses and numanalyses, create parameters, and plug them into parameters container: */
	parameters->AddObject(new StringArrayParam(ToolkitsOptionsStringsEnum,strings,numanalyses));
	parameters->AddObject(new DoubleVecParam(ToolkitsOptionsAnalysesEnum,analyses,numanalyses));

	/*Clean up and return*/
	for(i=0;i<numanalyses;i++) xDelete<char>(strings[i]);
	xDelete<char*>(strings);
	xDelete<IssmDouble>(analyses);
	return;
}
