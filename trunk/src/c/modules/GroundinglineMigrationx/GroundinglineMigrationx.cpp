/*!\file GroundinglineMigrationx
 * \brief: migration grounding line position.
 */

#include "../../shared/shared.h"
#include "../../toolkits/toolkits.h"
#include "../../classes/classes.h"
#include "./GroundinglineMigrationx.h"

void GroundinglineMigrationx(Elements* elements,Nodes* nodes, Vertices* vertices,Loads* loads,Materials* materials, Parameters* parameters){

	int                 migration_style,analysis_type;
	IssmDouble         *vertices_potentially_ungrounding = NULL;
	IssmDouble         *phi_ungrounding                  = NULL;
	Element            *element                          = NULL;


	/*retrieve parameters: */
	parameters->FindParam(&migration_style,GroundinglineMigrationEnum);
	parameters->FindParam(&analysis_type,AnalysisTypeEnum);

	if(migration_style==NoneEnum) return;
	
	if(VerboseModule()) _printf0_("   Migrating grounding line\n");

	/*Set toolkit to default*/
	ToolkitsOptionsFromAnalysis(parameters,DefaultAnalysisEnum);

	switch(migration_style){
		case SoftMigrationEnum:
			/*Create flag for grounded vertices above the hydrostatic equilibrium: */
			vertices_potentially_ungrounding=PotentialUngrounding(elements,vertices,parameters);
			/*propagate ice shelf into connex areas of the ice sheet that potentially want to unground: */
			phi_ungrounding=PropagateFloatingiceToGroundedNeighbors(elements,nodes,vertices,parameters,vertices_potentially_ungrounding);
			break;
		case ContactEnum:
			phi_ungrounding=ContactFSLevelset(elements,vertices);
			break;
		case AggressiveMigrationEnum:
		case SubelementMigrationEnum:
		case SubelementMigration2Enum:
		case GroundingOnlyEnum:
			/*Nothing additional to do here, MigrateGroundingLine takes care of everything*/
			break;
		default:
			_error_("Grounding line migration "<<EnumToStringx(migration_style) << " not supported yet!");
	}

	/*Migrate grounding line : */
	for(int i=0;i<elements->Size();i++){
		element=xDynamicCast<Element*>(elements->GetObjectByOffset(i));
		element->MigrateGroundingLine(phi_ungrounding);
	}

	/*free ressouces: */
	xDelete<IssmDouble>(vertices_potentially_ungrounding);
	xDelete<IssmDouble>(phi_ungrounding);
}

IssmDouble*    ContactFSLevelset(Elements* elements,Vertices* vertices){ /*{{{*/

	Vector<IssmDouble>* vertexgrounded = NULL;
	Vector<IssmDouble>* vertexfloating = NULL;
	IssmDouble*  serial_vertexgrounded = NULL;
	IssmDouble*  serial_vertexfloating = NULL;
	IssmDouble*  phi                   = NULL;

	/*Initialize vector with number of vertices*/
	int numberofvertices = vertices->NumberOfVertices();
	vertexgrounded = new Vector<IssmDouble>(numberofvertices);
	vertexfloating = new Vector<IssmDouble>(numberofvertices);
	phi            = xNew<IssmDouble>(numberofvertices);

	/*Fill vector vertices_potentially_floating: */
	for(int i=0;i<elements->Size();i++){
		Element* element=xDynamicCast<Element*>(elements->GetObjectByOffset(i));
		element->FSContactMigration(vertexgrounded,vertexfloating);
	}

	/*Assemble vector and serialize */
	vertexgrounded->Assemble();
	vertexfloating->Assemble();
	serial_vertexgrounded=vertexgrounded->ToMPISerial();
	serial_vertexfloating=vertexfloating->ToMPISerial();
	for(int i=0;i<numberofvertices;i++){
			if (serial_vertexgrounded[i]==1. && serial_vertexfloating[i]==1.) phi[i]=0.;
			else if (serial_vertexgrounded[i]==1) phi[i]=1;
			else if (serial_vertexfloating[i]==1) phi[i]=-1;
			else if (serial_vertexgrounded[i]>10) phi[i]=9999;
			else phi[i]=-9999;
	}

	/*free ressouces and return: */
	delete vertexgrounded;
	delete vertexfloating;
	xDelete<IssmDouble>(serial_vertexgrounded);
	xDelete<IssmDouble>(serial_vertexfloating);

	return phi;
}
/*}}}*/
IssmDouble*    PotentialUngrounding(Elements* elements,Vertices* vertices,Parameters* parameters){ /*{{{*/

	int                 i,numberofvertices;
	IssmDouble*         vertices_potentially_ungrounding      = NULL;
	Vector<IssmDouble>* vec_vertices_potentially_ungrounding  = NULL;
	Element*            element                               = NULL;

	/*Initialize vector with number of vertices*/
	numberofvertices=vertices->NumberOfVertices();
	vec_vertices_potentially_ungrounding=new Vector<IssmDouble>(numberofvertices); //grounded vertex that could start floating

	/*Fill vector vertices_potentially_floating: */
	for(i=0;i<elements->Size();i++){
		element=xDynamicCast<Element*>(elements->GetObjectByOffset(i));
		element->PotentialUngrounding(vec_vertices_potentially_ungrounding);
	}

	/*Assemble vector and serialize */
	vec_vertices_potentially_ungrounding->Assemble();
	vertices_potentially_ungrounding=vec_vertices_potentially_ungrounding->ToMPISerial();

	/*free ressouces and return: */
	delete vec_vertices_potentially_ungrounding;
	return vertices_potentially_ungrounding;
}
/*}}}*/
IssmDouble*    PropagateFloatingiceToGroundedNeighbors(Elements* elements,Nodes* nodes,Vertices* vertices,Parameters* parameters,IssmDouble* vertices_potentially_ungrounding){ /*{{{*/
	int                 i,analysis_type;
	int                 nflipped,local_nflipped;
	IssmDouble*         phi                                  = NULL;
	IssmDouble*         elements_neighboring_floatingce      = NULL;
	Vector<IssmDouble>* vec_elements_neighboring_floatingice = NULL;
	Vector<IssmDouble>* vec_phi                              = NULL;
	Element*            element                               = NULL;

	/*recover parameters: */
	parameters->FindParam(&analysis_type,AnalysisTypeEnum);

	/*recover vec_phi*/
	vec_phi=new Vector<IssmDouble>(vertices->NumberOfVertices());
	for(i=0;i<elements->Size();i++){
		Element* element=xDynamicCast<Element*>(elements->GetObjectByOffset(i));
		element->GetVectorFromInputs(vec_phi,MaskGroundediceLevelsetEnum,VertexPIdEnum);
	}
	vec_phi->Assemble();
	phi=vec_phi->ToMPISerial();

	nflipped=1; //bootstrap
	while(nflipped){

		/*Vector of size number of elements*/
		vec_elements_neighboring_floatingice=new Vector<IssmDouble>(elements->NumberOfElements(),true);

		/*Figure out if any of the nodes of the element will be floating -> elements neighbouting the floating ice*/
		for(i=0;i<elements->Size();i++){
			element=xDynamicCast<Element*>(elements->GetObjectByOffset(i));
			vec_elements_neighboring_floatingice->SetValue(element->Sid(),element->IsNodeOnShelfFromFlags(phi)?1.0:0.0,INS_VAL);
		}

		/*Assemble vector and serialize: */
		vec_elements_neighboring_floatingice->Assemble();
		elements_neighboring_floatingce=vec_elements_neighboring_floatingice->ToMPISerial();

		/*Go through elements_neighboring_floatingce, and update vector of the nodes that will start floating*/
		local_nflipped=0;
		for(i=0;i<elements->Size();i++){
			element=xDynamicCast<Element*>(elements->GetObjectByOffset(i));
			if(reCast<int,IssmDouble>(elements_neighboring_floatingce[element->Sid()])){
				local_nflipped+=element->UpdatePotentialUngrounding(vertices_potentially_ungrounding,vec_phi,phi);
			}
		}
		vec_phi->Assemble();

		ISSM_MPI_Allreduce(&local_nflipped,&nflipped,1,ISSM_MPI_INT,ISSM_MPI_SUM,IssmComm::GetComm());
		if(VerboseConvergence()) _printf0_("   Additional number of vertices allowed to unground: " << nflipped << "\n");

		/*Avoid leaks: */
		xDelete<IssmDouble>(elements_neighboring_floatingce);
		xDelete<IssmDouble>(phi);

		/*Assemble and serialize:*/
		delete vec_elements_neighboring_floatingice;
		phi=vec_phi->ToMPISerial();
	}

	/*Free ressources:*/
	delete vec_phi;
	xDelete<IssmDouble>(elements_neighboring_floatingce);

	return phi;
}
/*}}}*/
