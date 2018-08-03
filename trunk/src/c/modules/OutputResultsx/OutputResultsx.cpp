/*!\file:  OutputResultsx.cpp
 * \brief: go through our finite elements, and see what results they have stored within. 
 * Then output them into serialized patch arrays, and dump to disk.
 */ 

#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#include <stdio.h>
#include "./OutputResultsx.h"
#include "../../shared/io/io.h"
#include "../../classes/classes.h"

void OutputResultsx(FemModel* femmodel){

	int         my_rank;
	FILE       *fid                     = NULL;
	char       *outputfilename          = NULL;
	char        cpu_outputfilename[100];        //easier to convert an integer with sprintf
	bool        io_gather;
	int         solutiontype;
	char*       solutiontypestring      = NULL;
	bool        dakota_analysis         = false;

	/*retrieve parameters: */
	femmodel->parameters->FindParam(&dakota_analysis,QmuIsdakotaEnum);

	/*recover my_rank:*/
	my_rank=IssmComm::GetRank();

	if(dakota_analysis){
		//no need to output anything, Dakota analysis has different outputs
		return; 
	}

	/*Results do not include the type of solution being run	. In parallel, we output results to a filename, 
	 *therefore, we need to include the solutiontype into the filename: */
	if(my_rank==0){
		femmodel->parameters->FindParam(&solutiontype,SolutionTypeEnum);
		EnumToStringx(&solutiontypestring,solutiontype);
		femmodel->results->AddResult(new GenericExternalResult<char*>(femmodel->results->Size()+1,SolutionTypeEnum,solutiontypestring));
		xDelete<char>(solutiontypestring);
	}

#ifdef _HAVE_JAVASCRIPT_
	femmodel->parameters->FindParam(&fid,OutputFilePointerEnum);
#else

	/*Now, open file for writing*/
	_assert_(!femmodel->parameters->Exist(OutputFilePointerEnum));
	femmodel->parameters->FindParam(&outputfilename,OutputFileNameEnum);
	femmodel->parameters->FindParam(&io_gather,SettingsIoGatherEnum);

	if(io_gather){
		/*Just open the file for output on cpu 0. We are gathering the data on cpu 0 from all other cpus: */
		if(my_rank==0) fid=pfopen0(outputfilename ,"ab+");
	}
	else{
		/*We are opening different  files for output on all cpus. Append the  rank to the filename, and open: */
		sprintf(cpu_outputfilename,"%s.%i",outputfilename,my_rank);
		fid=pfopen(cpu_outputfilename ,"ab+");
	}

	/*Add file pointer in parameters for further calls to OutputResultsx: */
	femmodel->parameters->SetParam(fid,OutputFilePointerEnum);
#endif

	/*Write results to disk: */
	femmodel->results->Write(femmodel->parameters);

	femmodel->parameters->Delete(OutputFilePointerEnum);

	/*Delete and reinitialize results, in parallel: */
	femmodel->results->clear();

#ifndef _HAVE_JAVASCRIPT_
	/*Close output file? :*/
	if(io_gather){
		if(my_rank==0) pfclose(fid,outputfilename);
	}
	else pfclose(fid,cpu_outputfilename);
#endif

	/*Clean up and return*/
	xDelete<char>(outputfilename);
}
