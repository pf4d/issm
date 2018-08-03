/*!\file:  ProcessArguments.cpp
 * \brief: process arguments
 */ 

#include <stdio.h>
#include <cstring>

#include "../shared/shared.h"

void ProcessArguments(int* solution_type,char** pbinfilename,char** poutbinfilename,char** ptoolkitsfilename,char** plockfilename,char** prestartfilename, char** prootpath, int argc,char **argv){

	/*Check input arguments*/
	if(argc<2) _error_("Usage error: no solution requested");
	if(argc<3) _error_("Usage error: missing execution directory");
	if(argc<4) _error_("Usage error: missing model name");

	/*Get some arguments*/
	*solution_type = StringToEnumx(argv[1]);
	char* rootpatharg = argv[2];
	char* modelname   = argv[3];

	/*Recover myrank and length of string "my_rank" */
	int my_rank     = IssmComm::GetRank();
	int rank_length = (my_rank == 0 ? 1 : (int)(log10(static_cast<double>(my_rank))+1)); 

	/*Create rootpath from argument*/
	char* rootpath = xNew<char>(strlen(rootpatharg)+2); sprintf(rootpath,"%s/",rootpatharg);

	/*Create all file paths*/
	char* binfilename      = xNew<char>(strlen(rootpath)+strlen(modelname)+strlen(".bin")      +1); sprintf(binfilename,   "%s%s%s",rootpath,modelname,".bin");
	char* outbinfilename   = xNew<char>(strlen(rootpath)+strlen(modelname)+strlen(".outbin")   +1); sprintf(outbinfilename,"%s%s%s",rootpath,modelname,".outbin");
	char* toolkitsfilename = xNew<char>(strlen(rootpath)+strlen(modelname)+strlen(".toolkits") +1); sprintf(toolkitsfilename,"%s%s%s",rootpath,modelname,".toolkits");
	char* lockfilename     = xNew<char>(strlen(rootpath)+strlen(modelname)+strlen(".lock")     +1); sprintf(lockfilename,  "%s%s%s",rootpath,modelname,".lock");
	char* restartfilename  = xNew<char>(strlen(rootpath)+strlen(modelname)+strlen(".rst.")+rank_length+1);
	sprintf(restartfilename,"%s%s%s%i",rootpath,modelname,".rst.",my_rank);

	/*Clean up and assign output pointer*/
	*pbinfilename=binfilename;
	*poutbinfilename=outbinfilename;
	*ptoolkitsfilename=toolkitsfilename;
	*plockfilename=lockfilename;
	*prestartfilename=restartfilename;
	*prootpath=rootpath;

}
