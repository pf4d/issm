#ifndef _CONTAINER_PARAMETERS_H_
#define  _CONTAINER_PARAMETERS_H_
#include <stdio.h>

/*forward declarations */
class Param;
class DataSet;
template <class doublematrix> class Matrix;
template <class doubletype> class Vector;
#include "../../shared/shared.h"

#define NUMPARAMS ParametersENDEnum - ParametersSTARTEnum -1

/*!\brief Declaration of Parameters class.  
 *
 * Declaration of Parameters class.  Parameters are a static array of Parameter objects.
 */ 
class Parameters{

	private:
		Param* params[NUMPARAMS];

	public:

		/*constructors, destructors*/ 
		Parameters();
		~Parameters();

		/*numerics*/
		void  AddObject(Param* newparam);
		Parameters* Copy(void);
		void  DeepEcho();
		void  Echo();
		void  Delete(int enum_type);
		bool  Exist(int enum_type);
		void  Marshall(char** pmarshalled_data, int* pmarshalled_data_size, int marshall_direction);

		void  FindParam(bool* pinteger,int enum_type);
		void  FindParam(int* pinteger,int enum_type);
		void  FindParam(IssmDouble* pscalar, int enum_type);
		void  FindParam(IssmDouble* pscalar, int enum_type,IssmDouble time);
		void  FindParam(char** pstring,int enum_type);
		void  FindParam(char*** pstringarray,int* pM,int enum_type);
		void  FindParam(int** pintarray,int* pM,int enum_type);
		void  FindParam(int** pintarray,int* pM,int* PN,int enum_type);
		void  FindParam(IssmDouble** pIssmDoublearray,int* pM,int enum_type);
		void  FindParam(IssmDouble** pIssmDoublearray,int* pM,int* pN,int enum_type);
		void  FindParam(IssmDouble*** parray,int* pM, int** pmdims_array,int** pndims_array,int enum_type);
		void  FindParam(Vector<IssmDouble>** pvec,int enum_type);
		void  FindParam(Matrix<IssmDouble>** pmat,int enum_type);
		void  FindParam(FILE** pfid,int enum_type);
		void  FindParam(DataSet** pdataset, int enum_type);

		void  SetParam(bool boolean,int enum_type);
		void  SetParam(int integer,int enum_type);
		void  SetParam(IssmDouble scalar, int enum_type);
		void  SetParam(char* string,int enum_type);
		void  SetParam(char** stringarray,int M,int enum_type);
		void  SetParam(IssmDouble* IssmDoublearray,int M,int enum_type);
		void  SetParam(IssmDouble* IssmDoublearray,int M,int N,int enum_type);
		void  SetParam(int* intarray,int M,int enum_type);
		void  SetParam(int* intarray,int M,int N,int enum_type);
		void  SetParam(Vector<IssmDouble>* vec,int enum_type);
		void  SetParam(Matrix<IssmDouble>* mat,int enum_type);
		void  SetParam(FILE* fid,int enum_type);
		void  SetParam(DataSet* dataset,int enum_type);

		Param* FindParamObject(int enum_type);

};

/*Methods relating to parameters: */
char *OptionsFromAnalysis(Parameters *parameters,int analysis_type);
void  ToolkitsOptionsFromAnalysis(Parameters* parameters,int analysis_type);

#endif //ifndef _PARAMETERS_H_
