/*!\file:  FacesPartitioning.cpp
 * \brief: partition elements and nodes and vertices
 */ 

#include <string.h>
#include "../../classes/classes.h"
#include "../../shared/shared.h"
#include "./ModelProcessorx.h"

void FacesPartitioning(bool** pmy_faces,IoModel* iomodel){

	/*Intermediaries*/
	int elementnbf;

	/*Get faces and elements*/
	CreateFaces(iomodel);
	_assert_(iomodel->elementtofaceconnectivity);

	/*Mesh dependent variables*/
	if(iomodel->domaintype==Domain2DhorizontalEnum){
		elementnbf = 3;
	}
	else if(iomodel->domaintype==Domain2DverticalEnum){
		elementnbf = 3;
	}
	else if(iomodel->domaintype==Domain3DEnum){
		elementnbf = 5;
	}
	else{
		_error_("mesh dimension not supported yet");
	}
	/*output: */
	bool* my_faces=xNewZeroInit<bool>(iomodel->numberoffaces);

	for(int i=0;i<iomodel->numberofelements;i++){
		if(iomodel->my_elements[i]){
			for(int j=0;j<elementnbf;j++){
				_assert_(iomodel->elementtofaceconnectivity[i*elementnbf+j] >= 0);
				_assert_(iomodel->elementtofaceconnectivity[i*elementnbf+j] <  iomodel->numberoffaces);
				my_faces[iomodel->elementtofaceconnectivity[i*elementnbf+j]] = true;
			}
		}
	}

	/*Free data and assign output pointers */
	*pmy_faces=my_faces;
}
