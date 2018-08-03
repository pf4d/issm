/*!\file Mergesolutionfromftogx
 * \brief merge solution back from f set into g set
 */

#include "../VecMergex/VecMergex.h"
#include "../../shared/io/io.h"
#include "./Mergesolutionfromftogx.h"

void	Mergesolutionfromftogx( Vector<IssmDouble>** pug, Vector<IssmDouble>* uf, Vector<IssmDouble>* ys, Nodes* nodes, Parameters* parameters, bool flag_ys0){

	/*output: */
	Vector<IssmDouble>* ug=NULL;

	/*intermediary: */
	int configuration_type;
	int gsize,fsize,ssize;

	/*Display message*/
	if(VerboseModule()) _printf0_("   Merging solution vector from fset to gset\n");

	/*first, get gsize, fsize and ssize: */
	parameters->FindParam(&configuration_type,ConfigurationTypeEnum);
	gsize=nodes->NumberOfDofs(configuration_type,GsetEnum);
	fsize=nodes->NumberOfDofs(configuration_type,FsetEnum);
	ssize=nodes->NumberOfDofs(configuration_type,SsetEnum);

	/*serialize uf and ys: those two vectors will be indexed by the nodes, who are the only ones 
	 *that know which values should be plugged into ug and where: */
	if(ssize){
		if(flag_ys0){
			ys->Set(0.0);
		}
	}

	/*initialize ug: */
	ug=new Vector<IssmDouble>(gsize);

	/*Merge f set back into g set: */
	if(fsize){
		VecMergex(ug,uf,nodes,parameters,FsetEnum);
	}

	/*Merge s set back into g set: */
	if(ssize){
		VecMergex(ug,ys,nodes,parameters,SsetEnum);
	}

	/*Assign correct pointer*/
	*pug=ug;
}
