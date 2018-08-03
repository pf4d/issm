/*!\file:  SetActiveNodesLSMx.h
 * \brief header file for updating single point constraints  for next time step
 */ 

#ifndef _SETACTIVENODESLSMX_H
#define _SETACTIVENODESLSMX_H

#include "../../classes/classes.h"

void SetActiveNodesLSMx(FemModel* femmodel);
void GetMaskOfIceVerticesLSMx(FemModel* femmodel);
void SetMaskOfIceElement(Vector<IssmDouble>* vec_mask_ice, Element* element);
#endif  /* _UPDATESPCSX_H */
