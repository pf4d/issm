/*!\file Misfit.h
 * \brief: header file for Misfit object
 */

#ifndef _MISFIT_H_
#define _MISFIT_H_

/*Headers:*/
/*{{{*/
#include "./Definition.h"
#include "../datastructures/datastructures.h"
#include "./Elements/Element.h"
#include "./Elements/Elements.h"
#include "./FemModel.h"
#include "../modules/SurfaceAreax/SurfaceAreax.h"
#include "../classes/Params/Parameters.h"
#include "../classes/Inputs/Input.h"
#include "../classes/gauss/Gauss.h"
/*}}}*/
IssmDouble OutputDefinitionsResponsex(FemModel* femmodel,int output_enum);

class Misfit: public Object, public Definition{

	public: 

		int         definitionenum;
		bool        local;     
		int         model_enum;
		char*       name;
		int         observation_enum;
		char*       timeinterpolation;
		int         weights_enum;
		
		int         lock; // if lock is on, we just return the value stored in "misfit".  this is used so we don't compute misfit past the final_time
		IssmDouble  misfit; //value carried over in time.
		
		/*Misfit constructors, destructors :*/
		Misfit(){/*{{{*/

			this->definitionenum = -1;
			this->name = NULL;
			this->model_enum = UNDEF;
			this->observation_enum = UNDEF;
			this->weights_enum = UNDEF;
			this->timeinterpolation=NULL;
			this->local=true;
			this->misfit=0;
			this->lock=0;

		}
		/*}}}*/
		Misfit(char* in_name, int in_definitionenum, int in_model_enum, int in_observation_enum, char* in_timeinterpolation, bool in_local, int in_weights_enum){/*{{{*/

			this->definitionenum=in_definitionenum;
			
			this->name		= xNew<char>(strlen(in_name)+1);
			xMemCpy<char>(this->name,in_name,strlen(in_name)+1);

			this->timeinterpolation = xNew<char>(strlen(in_timeinterpolation)+1);
			xMemCpy<char>(this->timeinterpolation,in_timeinterpolation,strlen(in_timeinterpolation)+1);
						
			this->model_enum=in_model_enum;
			this->observation_enum=in_observation_enum;
			this->weights_enum=in_weights_enum;
			this->local=in_local;
			
			this->misfit=0;
			this->lock=0;
		}
		/*}}}*/
		~Misfit(){/*{{{*/
			if(this->name)xDelete(this->name);
			if(this->timeinterpolation)xDelete(this->timeinterpolation);
			this->misfit=0;
			this->lock=0;
		}
		/*}}}*/
		/*Object virtual function resolutoin: */
		Object* copy() {/*{{{*/
			Misfit* mf = new Misfit(this->name,this->definitionenum, this->model_enum,this->observation_enum,this->timeinterpolation,this->local,this->weights_enum);
			mf->misfit=this->misfit;
			mf->lock=this->lock;
			return (Object*) mf;
		}
		/*}}}*/
		void DeepEcho(void){/*{{{*/
			this->Echo();
		}
		/*}}}*/
		void Echo(void){/*{{{*/
			_printf_(" Misfit: " << name << " " << this->definitionenum << "\n");
			_printf_("    model_enum: " << model_enum << " " << EnumToStringx(model_enum) << "\n");
			_printf_("    observation_enum: " << observation_enum << " " << EnumToStringx(observation_enum) << "\n");
			_printf_("    weights_enum: " << weights_enum << " " << EnumToStringx(weights_enum) << "\n");
			_printf_("    timeinterpolation: " << timeinterpolation << "\n");
			_printf_("    local: " << local << "\n");
		}
		/*}}}*/
		int Id(void){/*{{{*/
			return -1;
		}
		/*}}}*/
		void Marshall(char** pmarshalled_data,int* pmarshalled_data_size, int marshall_direction){/*{{{*/
			_error_("not implemented yet!"); 
		} 
		/*}}}*/
		int ObjectEnum(void){/*{{{*/
			return MisfitEnum;
		}
		/*}}}*/
		/*Definition virtual function resolutoin: */
		int DefinitionEnum(){/*{{{*/
			return this->definitionenum;
		}
		/*}}}*/
		char* Name(){/*{{{*/
			char* name2=xNew<char>(strlen(this->name)+1);
			xMemCpy(name2,this->name,strlen(this->name)+1);

			return name2;
		}
		/*}}}*/
		 IssmDouble Response(FemModel* femmodel){/*{{{*/
				 
			 /*diverse: */
			 IssmDouble time,starttime,finaltime;
			 IssmDouble dt;
			 
			 /*recover time parameters: */
			 femmodel->parameters->FindParam(&starttime,TimesteppingStartTimeEnum);
			 femmodel->parameters->FindParam(&finaltime,TimesteppingFinalTimeEnum);
			 femmodel->parameters->FindParam(&time,TimeEnum);
			 femmodel->parameters->FindParam(&dt,TimesteppingTimeStepEnum);

			 if (this->local){ /*local computation: {{{*/

				 int i;
				 IssmDouble misfit_t=0.;
				 IssmDouble all_misfit_t=0.;
				 IssmDouble area_t=0.;
				 IssmDouble all_area_t;

			
				 /*If we are locked, return time average: */
				 if(this->lock)return misfit/(time-starttime);

				 for(i=0;i<femmodel->elements->Size();i++){
					 Element* element=(Element*)femmodel->elements->GetObjectByOffset(i);
					 misfit_t+=element->Misfit(model_enum,observation_enum,weights_enum);
					 area_t+=element->MisfitArea(weights_enum);
				 }

				 ISSM_MPI_Allreduce ( (void*)&misfit_t,(void*)&all_misfit_t,1,ISSM_MPI_DOUBLE,ISSM_MPI_SUM,IssmComm::GetComm());
				 ISSM_MPI_Allreduce ( (void*)&area_t,(void*)&all_area_t,1,ISSM_MPI_DOUBLE,ISSM_MPI_SUM,IssmComm::GetComm());
				 area_t=all_area_t;
				 misfit_t=all_misfit_t;
				 
				 /*Divide by surface area if not nill!: */
				 if (area_t!=0) misfit_t=misfit_t/area_t;

				 /*Add this time's contribution to curent misfit: */
				 misfit+=dt*misfit_t;

				 /*Do we lock? i.e. are we at final_time? :*/
				 if(time==finaltime)this->lock=1;

				 /*What we return is the value of misfit / time: */
				 return misfit/(time-starttime);
			 } /*}}}*/
			 else{ /*global computation: {{{ */
				 
				 IssmDouble model, observation;
				 
				 /*If we are locked, return time average: */
				 if(this->lock)return misfit/(time-starttime);

				 /*First, the global  model response: */
				 model=OutputDefinitionsResponsex(femmodel,this->model_enum);
				 /*Now, the observation is buried inside the elements, go fish it in the first element (cludgy, needs fixing): */
				 Element* element=(Element*)femmodel->elements->GetObjectByOffset(0); _assert_(element);
				 Input* input = element->GetInput(observation_enum); _assert_(input);
				 input->GetInputAverage(&observation);
				 
				 /*Add this time's contribution to curent misfit: */
				 misfit+=dt*(model-observation);
				 
				 /*Do we lock? i.e. are we at final_time? :*/
				 if(time==finaltime)this->lock=1;
				 
				 /*What we return is the value of misfit / time: */
				 return misfit/(time-starttime);
			 } /*}}}*/

		 }
			/*}}}*/
};

#endif  /* _MISFIT_H_ */
