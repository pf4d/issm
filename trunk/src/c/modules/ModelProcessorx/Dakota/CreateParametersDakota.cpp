/*!\file: CreateParametersDakota.cpp
 * \brief general driver for creating parameters dataset
 */ 

#include "../../../toolkits/toolkits.h"
#include "../../../classes/classes.h"
#include "../../../shared/shared.h"
#include "../../MeshPartitionx/MeshPartitionx.h"
#include "../ModelProcessorx.h"

void CreateParametersDakota(Parameters* parameters,IoModel* iomodel,char* rootpath){

	/*variable declarations*/
	int          i;
	int         *part                   = NULL;
	double      *dpart                  = NULL;
	char       **responsedescriptors    = NULL;
	int          numresponsedescriptors;
	char       **variabledescriptors    = NULL;
	int          numvariabledescriptors;
	char        *descriptor             = NULL;
	double      *dakota_parameter       = NULL;

	//qmu files
	char *qmuinname  = NULL;
	char *qmuerrname = NULL;
	char *qmuoutname = NULL;

	//descriptors:
	char tag[50];

	bool  dakota_analysis   = false;
	char *name              = NULL;
	int   numberofresponses;
	int   nrows,ncols;

	/*recover parameters: */
	iomodel->FindConstant(&dakota_analysis,"md.qmu.isdakota");

	if(dakota_analysis){

		iomodel->FindConstant(&name,"md.miscellaneous.name");
		iomodel->FindConstant(&numberofresponses,"md.qmu.numberofresponses");

		/*name of qmu input, error and output files*/
		qmuinname=xNew<char>((strlen(rootpath)+strlen(name)+strlen(".qmu.in")+1));
		sprintf(qmuinname,"%s%s%s",rootpath,name,".qmu.in");
		parameters->AddObject(new   StringParam(QmuInNameEnum,qmuinname));

		qmuoutname=xNew<char>((strlen(rootpath)+strlen(name)+strlen(".qmu.out")+1));
		sprintf(qmuoutname,"%s%s%s",rootpath,name,".qmu.out");
		parameters->AddObject(new   StringParam(QmuOutNameEnum,qmuoutname));

		qmuerrname=xNew<char>((strlen(rootpath)+strlen(name)+strlen(".qmu.err")+1));
		sprintf(qmuerrname,"%s%s%s",rootpath,name,".qmu.err");
		parameters->AddObject(new   StringParam(QmuErrNameEnum,qmuerrname));

		/*Fetch variable descriptors*/
		iomodel->FindConstant(&variabledescriptors,&numvariabledescriptors,"md.qmu.variabledescriptors");

		/*Ok, we have all the variable descriptors. Build a parameter with it: */
		parameters->AddObject(new StringArrayParam(QmuVariabledescriptorsEnum,variabledescriptors,numvariabledescriptors));

		/*Fetch response descriptors*/
		iomodel->FindConstant(&responsedescriptors,&numresponsedescriptors,"md.qmu.responsedescriptors");

		/*Ok, we have all the response descriptors. Build a parameter with it: */
		parameters->AddObject(new StringArrayParam(QmuResponsedescriptorsEnum,responsedescriptors,numresponsedescriptors));
		parameters->AddObject(new    IntParam(QmuNumberofresponsesEnum,numberofresponses));

		/*Deal with partitioning*/
		/*partition vertices in iomodel->qmu_npart parts, unless a partition is already present: */
		parameters->AddObject(iomodel->CopyConstantObject("md.qmu.numberofpartitions",QmuNumberofpartitionsEnum));
		iomodel->FetchData(&dpart,NULL,NULL,"md.qmu.partition");
		if(!dpart){
			/*Partition elements and vertices and nodes: */
			ElementsAndVerticesPartitioning(&iomodel->my_elements,&iomodel->my_vertices,iomodel);

			dpart=xNew<double>(iomodel->numberofvertices);
			for(i=0;i<iomodel->numberofvertices;i++)dpart[i]=iomodel->my_vertices[i];
		}
		parameters->AddObject(new DoubleVecParam(QmuPartitionEnum,dpart,iomodel->numberofvertices));

		/*Deal with data needed because of qmu variables*/
		for(i=0;i<numvariabledescriptors;i++){
			if (strncmp(variabledescriptors[i],"scaled_",7)==0){
				/*Ok, we are dealing with a variable that is distributed over nodes. Recover the name of the variable (ex: scaled_Thickness): */
				sscanf(variabledescriptors[i],"scaled_%s",tag);

				/*Get field name and input enum from tag*/
				char* fieldname  = NULL;
				int   param_enum = -1;
				FieldAndEnumFromCode(&param_enum,&fieldname,tag);

				/*Recover data: */
				iomodel->FetchData(&dakota_parameter,&nrows,&ncols,fieldname);
				if(nrows==iomodel->numberofvertices){
					parameters->AddObject(new DoubleMatParam(param_enum,dakota_parameter,nrows,ncols));
				}
				else{
					parameters->AddObject(new DoubleTransientMatParam(param_enum,dakota_parameter,nrows,ncols));
				}
				xDelete<double>(dakota_parameter);
				xDelete<char>(fieldname);
			}
		}

		/*clean-up*/
		for(i=0;i<numresponsedescriptors;i++){
			descriptor=responsedescriptors[i];
			xDelete<char>(descriptor);
		}
		xDelete<char*>(responsedescriptors);
		for(i=0;i<numvariabledescriptors;i++){
			descriptor=variabledescriptors[i];
			xDelete<char>(descriptor);
		}
		xDelete<char*>(variabledescriptors);
		xDelete<int>(part);
		xDelete<double>(dpart);
		xDelete<char>(qmuinname);
		xDelete<char>(qmuerrname);
		xDelete<char>(qmuoutname);
	}

	/*Free data*/
	xDelete<char>(name);
}
