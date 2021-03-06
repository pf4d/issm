/*!\file:  ElementMatrix.h
 * \brief container for information needed to plug element matrix generated by elements 
 * into the Kff and Kfs global matrices. 
 * This object will hold the element matrix on the g-set, the local as well as global 
 * dof lists in the f and s sets.
 */ 

#ifndef _ELEMENT_MATRIX_H_
#define _ELEMENT_MATRIX_H_

/*Headers:*/
#include "../../datastructures/datastructures.h"
#include "../../toolkits/toolkits.h"
#include "../../shared/shared.h"
class Node;
template <class doublematrix> class Matrix;
class Parameters;

class ElementMatrix{

	public:

		int      nrows;
		int      ncols;
		bool     dofsymmetrical;
		IssmDouble*  values;

		//gset
		int*     gglobaldoflist;

		/*row wise: */
		//fset
		int      row_fsize;
		int*     row_flocaldoflist;
		int*     row_fglobaldoflist;
		//sset
		int      row_ssize;
		int*     row_slocaldoflist;
		int*     row_sglobaldoflist;

		/*column wise: */
		//fset
		int      col_fsize;
		int*     col_flocaldoflist;
		int*     col_fglobaldoflist;
		//sset
		int      col_ssize;
		int*     col_slocaldoflist;
		int*     col_sglobaldoflist;

		/*ElementMatrix constructors, destructors*/
		ElementMatrix();
		ElementMatrix(ElementMatrix* Ke);
		ElementMatrix(ElementMatrix* Ke1,ElementMatrix* Ke2);
		ElementMatrix(ElementMatrix* Ke1,ElementMatrix* Ke2,ElementMatrix* Ke3);
		ElementMatrix(Node** nodes,int numnodes,Parameters* parameters,int approximation=NoneApproximationEnum);
		~ElementMatrix();

		/*ElementMatrix specific routines*/
		void AddDiagonalToGlobal(Vector<IssmDouble>* pf);
		void AddToGlobal(Matrix<IssmDouble>* Kff, Matrix<IssmDouble>* Kfs);
		void AddToGlobal(Matrix<IssmDouble>* Jff);
		void CheckConsistency(void);
		void Echo(void);
		void Init(ElementMatrix* Ke);
		void Lump(void);
		void StaticCondensation(int numindices,int* indices);
		void Transpose(void);
};
#endif //#ifndef _ELEMENT_MATRIX_H_
