/*!\file: sealevelrise_core_eustatic.cpp
 * \brief: eustatic core of the SLR solution (terms that are constant with respect to sea-level)
 */ 

#include "./cores.h"
#include "../toolkits/toolkits.h"
#include "../classes/classes.h"
#include "../shared/shared.h"
#include "../modules/modules.h"
#include "../solutionsequences/solutionsequences.h"

Vector<IssmDouble>* sealevelrise_core_eustatic(FemModel* femmodel){

	Vector<IssmDouble> *Sgi    = NULL;
	IssmDouble          Sgi_oceanaverage   = 0;

	/*parameters: */
	int  configuration_type;
	int  gsize;
	bool spherical=true;
	IssmDouble *latitude  = NULL;
	IssmDouble *longitude = NULL;
	IssmDouble *radius    = NULL;

	/*outputs:*/
	IssmDouble eustatic;

	/*recover parameters:*/
	femmodel->parameters->FindParam(&configuration_type,ConfigurationTypeEnum);

	/*first, recover lat,long and radius vectors from vertices: */
	VertexCoordinatesx(&latitude,&longitude,&radius,femmodel->vertices,spherical); 

	/*Figure out size of g-set deflection vector and allocate solution vector: */
	gsize = femmodel->nodes->NumberOfDofs(configuration_type,GsetEnum);
	
	/*Initialize:*/
	Sgi = new Vector<IssmDouble>(gsize);

	/*call the eustatic main module: */
	femmodel->SealevelriseEustatic(Sgi,&eustatic, latitude, longitude, radius); //this computes 

	/*we need to average Sgi over the ocean: RHS term  4 in Eq.4 of Farrel and clarke. Only the elements can do that: */
	Sgi_oceanaverage=femmodel->SealevelriseOceanAverage(Sgi);

	/*Sg is the sum of the pure eustatic component (term 3) and the contribution from the perturbation to the graviation potential due to the 
	 * presence of ice (terms 1 and 4 in Eq.4 of Farrel and Clarke):*/
	Sgi->Shift(-eustatic-Sgi_oceanaverage);

	/*save eustatic value for results: */
	femmodel->results->AddResult(new GenericExternalResult<IssmDouble>(femmodel->results->Size()+1,SealevelEustaticEnum,-eustatic));

	/*clean up and return:*/
	xDelete<IssmDouble>(latitude);
	xDelete<IssmDouble>(longitude);
	xDelete<IssmDouble>(radius);
	return Sgi;
}

