/*!\file Matice.h
 * \brief: header file for matice object
 */

#ifndef MATICE_H_
#define MATICE_H_

/*Headers:*/
/*{{{*/
#include "./Material.h"
#include "../Hook.h"
class IoModel;
class Elements;
class Element;
class Loads;
class Nodes;
class Vertices;
class Materials;
class Parameters;
class Gauss;
class Input;
/*}}}*/

class Matice: public Material{

	private: 
		int      mid;
		bool     isdamaged;
		bool     isenhanced;
		Hook    *helement;
		Element *element;

	public:
		/*Matice constructors, destructors: {{{*/
		Matice();
		Matice(int mid,int i, IoModel* iomodel);
		Matice(int mid,int i, int materialtype);
		~Matice();
		void Init(int mid,int i, int materialtype);
		/*}}}*/
		/*Object virtual functions definitions:{{{ */
		Object* copy();
		void  DeepEcho();
		void  Echo();
		int   Id(); 
		void Marshall(char** pmarshalled_data,int* pmarshalled_data_size, int marshall_direction);
		int   ObjectEnum();
		/*}}}*/
		/*Update virtual funictions definitions: {{{*/
		void  InputUpdateFromConstant(IssmDouble constant, int name);
		void  InputUpdateFromConstant(int constant, int name);
		void  InputUpdateFromConstant(bool constant, int name);
		void  InputUpdateFromIoModel(int index, IoModel* iomodel){_error_("not implemented");};
		void  InputUpdateFromMatrixDakota(IssmDouble* matrix, int nrow, int ncols, int name, int type);
		void  InputUpdateFromVector(IssmDouble* vector, int name, int type);
		void  InputUpdateFromVectorDakota(IssmDouble* vector, int name, int type);
		/*}}}*/
		/*Material virtual functions resolution: {{{*/
		void   Configure(Elements* elements);
		Material*  copy2(Element* element);
		void       GetViscosity(IssmDouble* pviscosity, IssmDouble eps_eff);
		void       GetViscosityBar(IssmDouble* pviscosity, IssmDouble eps_eff);
		void       GetViscosityComplement(IssmDouble* pviscosity_complement, IssmDouble* pepsilon);
		void       GetViscosityDComplement(IssmDouble*, IssmDouble*);
		void       GetViscosityDerivativeEpsSquare(IssmDouble* pmu_prime, IssmDouble* pepsilon);
		void       GetViscosity_B(IssmDouble* pviscosity, IssmDouble eps_eff);
		void       GetViscosity_D(IssmDouble* pviscosity, IssmDouble eps_eff);
		void       GetViscosity2dDerivativeEpsSquare(IssmDouble* pmu_prime, IssmDouble* pepsilon);
		IssmDouble GetA();
		IssmDouble GetAbar();
		IssmDouble GetB();
		IssmDouble GetBbar();
		IssmDouble GetD();
		IssmDouble GetDbar();
		IssmDouble GetE();
		IssmDouble GetEbar();
		IssmDouble GetN();
		bool       IsDamage();
		bool       IsEnhanced();
		void       ResetHooks();
		void       SetCurrentConfiguration(Elements* elements,Loads* loads,Nodes* nodes,Vertices* vertices,Materials* materials,Parameters* parameters);

		void       ViscosityFS(IssmDouble* pviscosity,int dim,IssmDouble* xyz_list,Gauss* gauss,Input* vx_input,Input* vy_input,Input* vz_input);
		void       ViscosityFSDerivativeEpsSquare(IssmDouble* pmu_prime,IssmDouble* epsilon);
		void       ViscosityHO(IssmDouble* pviscosity,int dim,IssmDouble* xyz_list,Gauss* gauss,Input* vx_input,Input* vy_input);
		void       ViscosityHODerivativeEpsSquare(IssmDouble* pmu_prime,IssmDouble* epsilon);
		void       ViscosityL1L2(IssmDouble* pviscosity,IssmDouble* xyz_list,Gauss* gauss,Input* vx_input,Input* vy_input,Input* surf);
		void       ViscositySSA(IssmDouble* pviscosity,int dim,IssmDouble* xyz_list,Gauss* gauss,Input* vx_input,Input* vy_input);
		void       ViscositySSADerivativeEpsSquare(IssmDouble* pmu_prime,IssmDouble* epsilon);
		void       ViscosityBFS(IssmDouble* pmudB,int dim,IssmDouble* xyz_list,Gauss* gauss,Input* vx_input,Input* vy_input,Input* vz_input){_error_("not supported");};
		void       ViscosityBHO(IssmDouble* pmudB,int dim,IssmDouble* xyz_list,Gauss* gauss,Input* vx_input,Input* vy_input){_error_("not supported");};
		void       ViscosityBSSA(IssmDouble* pmudB,int dim,IssmDouble* xyz_list,Gauss* gauss,Input* vx_input,Input* vy_input){_error_("not supported");};
		/*}}}*/
};

#endif  /* _MATICE_H_ */
