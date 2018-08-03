/*!\file Pengrid.h
 * \brief: header file for pengrid object */

#ifndef _PENGRID_H_
#define _PENGRID_H_

/*Headers:*/
/*{{{*/
#ifdef HAVE_CONFIG_H
#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif
#include "./Load.h"
class Hook;
class Inputs;
class Parameters;
class IoModel;
/*}}}*/

class Pengrid: public Load{

	private: 

		int        id;
		int        analysis_type;

		/*Hooks*/
		Hook* hnode;  //hook to 1 node
		Hook* helement;  //hook to 1 element
		Hook* hmatpar; //hook to 1 matpar

		/*Corresponding fields*/
		Node    *node;
		Element *element;
		Matpar  *matpar;

		Parameters* parameters; //pointer to solution parameters

		/*internals: */
		int active;
		int zigzag_counter;

	public:

		/*Pengrid constructors, destructors {{{*/
		Pengrid();
		Pengrid(int index, int id, IoModel* iomodel,int analysis_type);
		~Pengrid();
		/*}}}*/
		/*Object virtual functions definitions:{{{ */
		Object* copy();
		void  DeepEcho();
		void  Echo();
		int   Id(); 
		void Marshall(char** pmarshalled_data,int* pmarshalled_data_size, int marshall_direction);
		int   ObjectEnum();
		/*}}}*/
		/*Update virtual functions resolution: {{{*/
		void  InputUpdateFromConstant(IssmDouble constant, int name);
		void  InputUpdateFromConstant(int constant, int name);
		void  InputUpdateFromConstant(bool constant, int name);
		void  InputUpdateFromIoModel(int index, IoModel* iomodel){_error_("not implemented yet");};
		void  InputUpdateFromMatrixDakota(IssmDouble* matrix ,int nrows, int ncols, int name, int type);
		void  InputUpdateFromVector(IssmDouble* vector, int name, int type);
		void  InputUpdateFromVectorDakota(IssmDouble* vector, int name, int type);
		/*}}}*/
		/*Load virtual functions definitions: {{{*/
		void  Configure(Elements* elements,Loads* loads,Nodes* nodes,Vertices* vertices,Materials* materials,Parameters* parameters);
		void  CreateJacobianMatrix(Matrix<IssmDouble>* Jff){_error_("Not implemented yet");};
		void  CreateKMatrix(Matrix<IssmDouble>* Kff, Matrix<IssmDouble>* Kfs);
		void  CreatePVector(Vector<IssmDouble>* pf);
		void  GetNodesLidList(int* lidlist);
		void  GetNodesSidList(int* sidlist);
		int   GetNumberOfNodes(void);
		bool  InAnalysis(int analysis_type);
		bool  IsPenalty(void);
		void  PenaltyCreateJacobianMatrix(Matrix<IssmDouble>* Jff,IssmDouble kmax){_error_("Not implemented yet");};
		void  PenaltyCreateKMatrix(Matrix<IssmDouble>* Kff, Matrix<IssmDouble>* kfs, IssmDouble kmax);
		void  PenaltyCreatePVector(Vector<IssmDouble>* pf, IssmDouble kmax);
		void  ResetHooks();
		void  SetCurrentConfiguration(Elements* elements,Loads* loads,Nodes* nodes,Vertices* vertices,Materials* materials,Parameters* parameters);
		void  SetwiseNodeConnectivity(int* d_nz,int* o_nz,Node* node,bool* flags,int* flagsindices,int set1_enum,int set2_enum);
		/*}}}*/
		/*Pengrid management {{{*/
		void				ConstraintActivate(int* punstable);
		void           ConstraintActivateHydrologyDCInefficient(int* punstable);
		void           ConstraintActivateThermal(int* punstable);
		ElementMatrix* PenaltyCreateKMatrixHydrologyDCInefficient(IssmDouble kmax);
		ElementMatrix* PenaltyCreateKMatrixMelting(IssmDouble kmax);
		ElementMatrix* PenaltyCreateKMatrixThermal(IssmDouble kmax);
		ElementVector* PenaltyCreatePVectorHydrologyDCInefficient(IssmDouble kmax);
		ElementVector* PenaltyCreatePVectorMelting(IssmDouble kmax);
		ElementVector* PenaltyCreatePVectorThermal(IssmDouble kmax);
		void  ResetConstraint(void);
		/*}}}*/

};

#endif  /* _PENGRID_H_ */
