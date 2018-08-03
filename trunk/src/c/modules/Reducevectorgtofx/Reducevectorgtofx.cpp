/*!\file Reducevectorgtofx
 * \brief reduce petsc vector from g set to s set (free dofs), using the nodeset partitioning 
 * vectors.
 */

#include "./Reducevectorgtofx.h"

void Reducevectorgtofx(Vector<IssmDouble>** puf, Vector<IssmDouble>* ug, Nodes* nodes,Parameters* parameters){

	/*output: */
	Vector<IssmDouble>* uf=NULL;

	/*variables: */
	int         configuration_type;
	int         fsize;
	IssmDouble *ug_serial = NULL;
	bool        oldalloc  = false;

	if(VerboseModule()) _printf0_("   Reduce vector from g to f set\n");

	/*first figure out fsize: */
	parameters->FindParam(&configuration_type,ConfigurationTypeEnum);
	fsize=nodes->NumberOfDofs(configuration_type,FsetEnum);

	int    analysis_type;
	parameters->FindParam(&analysis_type,AnalysisTypeEnum);
	int flocalsize = nodes->NumberOfDofsLocal(analysis_type,FsetEnum);

	if(fsize==0){
		uf=NULL;
	}
	else{
		/*allocate: */
		if(oldalloc)
		 uf=new Vector<IssmDouble>(fsize);
		else
		 uf=new Vector<IssmDouble>(flocalsize,fsize);

		if(nodes->NumberOfNodes(configuration_type)){ 

			/*serialize ug, so nodes can index into it: */
			ug_serial=ug->ToMPISerial();

			/*Go through all nodes, and ask them to retrieve values from ug, and plug them into uf: */
			for(int i=0;i<nodes->Size();i++){

				Node* node=(Node*)nodes->GetObjectByOffset(i);

				/*Check that this node corresponds to our analysis currently being carried out: */
				if (node->InAnalysis(configuration_type)){

					/*For this object, reduce values for enum set Fset: */
					node->VecReduce(uf,ug_serial,FsetEnum);
				}
			}
		}
		/*Assemble vector: */
		uf->Assemble();
	}

	/*Free ressources:*/
	xDelete<IssmDouble>(ug_serial);

	/*Assign output pointers:*/
	*puf=uf;
}
