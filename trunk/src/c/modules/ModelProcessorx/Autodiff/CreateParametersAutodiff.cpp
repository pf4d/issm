/*!\file: CreateParametersAutodiff.cpp
 * \brief driver for creating parameters dataset, for autodiff analysis.
 */ 

#include "../../../toolkits/toolkits.h"
#include "../../../classes/classes.h"
#include "../../../shared/shared.h"
#include "../ModelProcessorx.h"

void CreateParametersAutodiff(Parameters* parameters,IoModel* iomodel){

	int         i;
	bool        isautodiff;
	int         num_dependent_objects;
	int         num_dep=0;
	char**      names=NULL;
	int*        types=NULL;
	int         dummy;
	char*       autodiff_driver=NULL;
	int*        indices=NULL;
	int         num_indices;
	char* options=NULL;

	IssmDouble* xp=NULL;
	IssmDouble* xp_backup=NULL;
	int         num_ind,local_num_ind;
	DataSet*    dependent_objects=NULL;

	/*retrieve some parameters: */
	iomodel->FindConstant(&isautodiff,"md.autodiff.isautodiff");

	#ifdef _HAVE_ADOLC_
		/*initialize a placeholder to store solver pointers: {{{*/
		GenericParam<Adolc_edf> *theAdolcEDF_p=new GenericParam<Adolc_edf>(AdolcParamEnum);

		/*Solver pointers depend on what type of solver we are implementing: */
		options=OptionsFromAnalysis(parameters,DefaultAnalysisEnum); //options database is not filled in yet, use default.
		ToolkitOptions::Init(options);

		switch(IssmSolverTypeFromToolkitOptions()){
			case MumpsEnum:{
				#ifdef _HAVE_MUMPS_
				theAdolcEDF_p->GetParameterValue().myEDF_for_solverx_p=reg_ext_fct(mumpsSolveEDF);
				#else
				_error_("requesting mumps solver without MUMPS being compiled in!");
				#endif
				break;
							}
			case GslEnum: {
				#ifdef _HAVE_GSL_
				theAdolcEDF_p->GetParameterValue().myEDF_for_solverx_p=reg_ext_fct(EDF_for_solverx);
				#else
				_error_("requesting GSL solver without GSL being compiled in!");
				#endif
			    break;
						  }
			default:
				_error_("solver type not supported yet!");
		}

		// to save some space:
		// we know we won't use adolc inside of  the solver:
		theAdolcEDF_p->GetParameterValue().myEDF_for_solverx_p->nestedAdolc=false;
		// the solution vector is just allocated and doesn't have a meaningful prior value
		theAdolcEDF_p->GetParameterValue().myEDF_for_solverx_p->dp_y_priorRequired=false;
		// the solver wrapper makes sure the matrix and the right hand side don't change
		theAdolcEDF_p->GetParameterValue().myEDF_for_solverx_p->dp_x_changes=false;
		parameters->AddObject(theAdolcEDF_p);

		/*Free ressources: */
		xDelete<char>(options);
		/*}}}*/

	if(isautodiff){
		/*Copy some parameters from IoModel to parameters dataset: {{{*/
		parameters->AddObject(iomodel->CopyConstantObject("md.autodiff.keep",AutodiffKeepEnum));
		parameters->AddObject(iomodel->CopyConstantObject("md.autodiff.obufsize",AutodiffObufsizeEnum));
		parameters->AddObject(iomodel->CopyConstantObject("md.autodiff.cbufsize",AutodiffCbufsizeEnum));
		parameters->AddObject(iomodel->CopyConstantObject("md.autodiff.lbufsize",AutodiffLbufsizeEnum));
		parameters->AddObject(iomodel->CopyConstantObject("md.autodiff.tbufsize",AutodiffTbufsizeEnum));
		parameters->AddObject(iomodel->CopyConstantObject("md.autodiff.gcTriggerRatio",AutodiffGcTriggerRatioEnum));
		parameters->AddObject(iomodel->CopyConstantObject("md.autodiff.gcTriggerMaxSize",AutodiffGcTriggerMaxSizeEnum));
		/*}}}*/
		/*retrieve driver: {{{*/
		iomodel->FindConstant(&autodiff_driver,"md.autodiff.driver");
		parameters->AddObject(iomodel->CopyConstantObject("md.autodiff.driver",AutodiffDriverEnum));

		if(strcmp(autodiff_driver,"fos_forward")==0){
			parameters->AddObject(iomodel->CopyConstantObject("md.autodiff.fos_forward_index",AutodiffFosForwardIndexEnum));
		}
		else if(strcmp(autodiff_driver,"fos_reverse")==0){
			parameters->AddObject(iomodel->CopyConstantObject("md.autodiff.fos_reverse_index",AutodiffFosReverseIndexEnum));
		}
		else if(strcmp(autodiff_driver,"fov_forward")==0){
			/*Retrieve list of indices: */
			iomodel->FetchData(&indices,&num_indices,&dummy,"md.autodiff.fov_forward_indices");
			parameters->AddObject(new IntMatParam(AutodiffFovForwardIndicesEnum,indices,num_indices,1));
			xDelete<int>(indices);
		}
		xDelete<char>(autodiff_driver);
		/*}}}*/
		/*Deal with dependents first: {{{*/
		iomodel->FindConstant(&num_dependent_objects,"md.autodiff.num_dependent_objects");
		dependent_objects=new DataSet();
		num_dep=0;

		if(num_dependent_objects){
			iomodel->FindConstant(&names,&dummy,"md.autodiff.dependent_object_names");
			iomodel->FetchData(&types,&dummy,&dummy,"md.autodiff.dependent_object_types");
			iomodel->FetchData(&indices,&dummy,&dummy,"md.autodiff.dependent_object_indices");

			for(i=0;i<num_dependent_objects;i++){
				DependentObject* dep=new DependentObject(names[i],types[i],indices[i]);
				dependent_objects->AddObject(dep);
				num_dep+=dep->NumDependents();
			}

			/*Free ressources:*/
			for(i=0;i<num_dependent_objects;i++){
				char* string=names[i]; xDelete<char>(string);
			}
			xDelete<char*>(names);
			xDelete<int>(types);
			xDelete<int>(indices);
		}
		parameters->AddObject(new DataSetParam(AutodiffDependentObjectsEnum,dependent_objects));
		parameters->AddObject(new IntParam(AutodiffNumDependentsEnum,num_dep));

		delete dependent_objects;
		/*}}}*/
		/*Deal with independents: {{{*/

		/*Independents have already been recovered in iomodel->DeclareIndependents. Just do some more processing. 
		 *In particular, figure out num_independents, and create the state vector xp, or size num_independents x 1 :*/
		num_ind=iomodel->NumIndependents();
		parameters->AddObject(new IntParam(AutodiffNumIndependentsEnum,num_ind));

		if(num_ind){
			xp=xNew<IssmDouble>(num_ind);
			iomodel->FillIndependents(xp);
			parameters->AddObject(new DoubleVecParam(AutodiffXpEnum,xp,num_ind));
			xDelete<IssmDouble>(xp);
		}
		/*}}}*/
	}
	#endif
}
