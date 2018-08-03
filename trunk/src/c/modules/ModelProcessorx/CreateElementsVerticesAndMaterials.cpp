/*
 * CreateElementsNodesAndMaterialsStressbalanceHoriz.c:
 */

#include "../../toolkits/toolkits.h"
#include "../../classes/classes.h"
#include "../../shared/shared.h"
#include "./ModelProcessorx.h"

void CreateElementsVerticesAndMaterials(Elements* elements,Vertices* vertices,Materials* materials,IoModel* iomodel,const int nummodels){

	/*Intermediary*/
	int i;
	int materials_type;
	bool control_analysis;
	bool dakota_analysis;

	/*Fetch parameters: */
	iomodel->FindConstant(&control_analysis,"md.inversion.iscontrol");
	iomodel->FindConstant(&dakota_analysis,"md.qmu.isdakota");
	iomodel->FindConstant(&materials_type,"md.materials.type");

	/*Did we already create the elements? : */
	_assert_(elements->Size()==0);

	/*Create elements*/
	if(control_analysis)iomodel->FetchData(2,"md.inversion.min_parameters","md.inversion.max_parameters");
	switch(iomodel->meshelementtype){
		case TriaEnum:
			for(i=0;i<iomodel->numberofelements;i++){
				if(iomodel->my_elements[i]) elements->AddObject(new Tria(i+1,i,i,iomodel,nummodels));
			}
			break;
		case TetraEnum:
			for(i=0;i<iomodel->numberofelements;i++){
				if(iomodel->my_elements[i]) elements->AddObject(new Tetra(i+1,i,i,iomodel,nummodels));
			}
			break;
		case PentaEnum:
			iomodel->FetchData(2,"md.mesh.upperelements","md.mesh.lowerelements");
			for(i=0;i<iomodel->numberofelements;i++){
				if(iomodel->my_elements[i]) elements->AddObject(new Penta(i+1,i,i,iomodel,nummodels));
			}
			break;
		default:
			_error_("Mesh not supported yet");
	}

	/*Create materials*/
	switch(materials_type){
		case MaticeEnum:
			iomodel->FetchDataToInput(elements,"md.materials.rheology_B",MaterialsRheologyBEnum);
			iomodel->FetchDataToInput(elements,"md.materials.rheology_n",MaterialsRheologyNEnum);
			for (i=0;i<iomodel->numberofelements;i++) if(iomodel->my_elements[i]) materials->AddObject(new Matice(i+1,i,iomodel));
			switch(iomodel->domaindim){
				case 2:
					elements->InputDuplicate(MaterialsRheologyBEnum,MaterialsRheologyBbarEnum);
					break;
				case 3:
					break;
				default:
					_error_("Mesh not supported yet");
			}
			break;
		case MatenhancediceEnum:
			iomodel->FetchDataToInput(elements,"md.materials.rheology_B",MaterialsRheologyBEnum);
			iomodel->FetchDataToInput(elements,"md.materials.rheology_n",MaterialsRheologyNEnum);
			iomodel->FetchDataToInput(elements,"md.materials.rheology_E",MaterialsRheologyEEnum);
			for (i=0;i<iomodel->numberofelements;i++) if(iomodel->my_elements[i]) materials->AddObject(new Matice(i+1,i,iomodel));
			switch(iomodel->domaindim){
				case 2:
					elements->InputDuplicate(MaterialsRheologyBEnum,MaterialsRheologyBbarEnum);
					elements->InputDuplicate(MaterialsRheologyEEnum,MaterialsRheologyEbarEnum);
					break;
				case 3:
					break;
				default:
					_error_("Mesh not supported yet");
			}
			break;
		case MatdamageiceEnum:
			iomodel->FetchDataToInput(elements,"md.materials.rheology_B",MaterialsRheologyBEnum);
			iomodel->FetchDataToInput(elements,"md.materials.rheology_n",MaterialsRheologyNEnum);
			iomodel->FetchDataToInput(elements,"md.damage.D",DamageDEnum);
			for (i=0;i<iomodel->numberofelements;i++) if(iomodel->my_elements[i]) materials->AddObject(new Matice(i+1,i,iomodel));
			switch(iomodel->domaindim){
				case 2:
					elements->InputDuplicate(MaterialsRheologyBEnum,MaterialsRheologyBbarEnum);
					elements->InputDuplicate(DamageDEnum,DamageDbarEnum);
					break;
				case 3:
					break;
				default:
					_error_("Mesh not supported yet");
			}
			break;
		case MatestarEnum:
			iomodel->FetchDataToInput(elements,"md.materials.rheology_B",MaterialsRheologyBEnum);
			iomodel->FetchDataToInput(elements,"md.materials.rheology_Ec",MaterialsRheologyEcEnum);
			iomodel->FetchDataToInput(elements,"md.materials.rheology_Es",MaterialsRheologyEsEnum);
			for(i=0;i<iomodel->numberofelements;i++) if(iomodel->my_elements[i]) materials->AddObject(new Matestar(i+1,i,iomodel));
			switch(iomodel->domaindim){
				case 2:
					elements->InputDuplicate(MaterialsRheologyBEnum,MaterialsRheologyBbarEnum);
					elements->InputDuplicate(MaterialsRheologyEcEnum,MaterialsRheologyEcbarEnum);
					elements->InputDuplicate(MaterialsRheologyEsEnum,MaterialsRheologyEsbarEnum);
					break;
				case 3:
					break;
				default:
					_error_("Mesh not supported yet");
			}
			break;
		default:
			_error_("Materials "<<EnumToStringx(materials_type)<<" not supported");
	}

	/*Free data: */
	iomodel->DeleteData(7,"md.mesh.upperelements","md.mesh.lowerelements","md.material.rheology_B",
				"md.material.rheology_n","md.damage.D","md.inversion.min_parameters","md.inversion.max_parameters");

	/*Add new constant material property to materials, at the end: */
	materials->AddObject(new Matpar(iomodel->numberofelements+1,iomodel));//put it at the end of the materials

	/*Create vertices: */

	/*Fetch data:*/
	iomodel->FetchData(6,"md.mesh.x","md.mesh.y","md.mesh.z","md.geometry.base","md.geometry.thickness","md.mask.ice_levelset");
	if (iomodel->domaintype == Domain3DsurfaceEnum) iomodel->FetchData(3,"md.mesh.lat","md.mesh.long","md.mesh.r");
	
	CreateNumberNodeToElementConnectivity(iomodel);

	for(i=0;i<iomodel->numberofvertices;i++){
		if(iomodel->my_vertices[i]) vertices->AddObject(new Vertex(i+1,i,i,iomodel));
	}

	/*Free data: */
	iomodel->DeleteData(6,"md.mesh.x","md.mesh.y","md.mesh.z","md.geometry.base","md.geometry.thickness","md.mask.ice_levelset");
	if (iomodel->domaintype == Domain3DsurfaceEnum) iomodel->DeleteData(3,"md.mesh.lat","md.mesh.long","md.mesh.r");
}
