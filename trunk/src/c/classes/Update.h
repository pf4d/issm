/*!\file: Update.h: abstract class used by some objects to update their properties
 * \brief prototypes for Update.h
 */ 

#ifndef _UPDATE_H_
#define  _UPDATE_H_

/*Headers:*/
#include "../shared/shared.h"
class IoModel;

class Update{

	public:

		virtual void  InputUpdateFromConstant(IssmDouble constant, int name)=0;
		virtual void  InputUpdateFromConstant(int constant, int name)=0;
		virtual void  InputUpdateFromConstant(bool constant, int name)=0;
		#ifdef _HAVE_DAKOTA_
		virtual void  InputUpdateFromMatrixDakota(IssmDouble* matrix, int rows, int ncols, int name, int type)=0;
		virtual void  InputUpdateFromVectorDakota(IssmDouble* vector, int name, int type)=0;
		#endif
		virtual void  InputUpdateFromIoModel(int index, IoModel* iomodel)=0;
		virtual void  InputUpdateFromVector(IssmDouble* vector, int name, int type)=0;

};

#endif //ifndef _UPDATE_H_
