/*!\file:  InputUpdateFromConstantx.h
 * \brief header file for updating datasets from inputs
 */ 

#ifndef _UPDATEINPUTSFROMCONSTANTXX_H
#define _UPDATEINPUTSFROMCONSTANTXX_H

#include "../../classes/classes.h"

/* local prototypes: */
void InputUpdateFromConstantx(FemModel* femmodel,bool       constant,int name);
void InputUpdateFromConstantx(FemModel* femmodel,int        constant,int name);
void InputUpdateFromConstantx(FemModel* femmodel,IssmDouble constant,int name);
void InputUpdateFromConstantx(Elements* elements,IssmDouble constant,int name);

#endif  /* _UPDATEINPUTSFROMCONSTANTXX_H */
