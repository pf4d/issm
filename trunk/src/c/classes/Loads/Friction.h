/*!\file Friction.h
 * \brief: header file for friction object
 */

#ifndef _FRICTION_H_
#define _FRICTION_H_

/*Headers:*/
/*{{{*/
class Inputs;
class Matpar;
class GaussPenta;
class GaussTria;
/*}}}*/

class Friction{

	public:
		int analysis_type;

		Element* element;
		int      dim;
		int      law;

		/*methods: */
		Friction();
		Friction(Element* element_in,int dim_in);
		~Friction();

		void  Echo(void);
		void  GetAlphaComplement(IssmDouble* alpha_complement,Gauss* gauss);
		void  GetAlphaHydroComplement(IssmDouble* alpha_complement,Gauss* gauss);
		void  GetAlphaTempComplement(IssmDouble* alpha_complement,Gauss* gauss);
		void  GetAlphaViscousComplement(IssmDouble* alpha_complement,Gauss* gauss);
		void  GetAlpha2(IssmDouble* palpha2,Gauss* gauss);
		void  GetAlpha2Coulomb(IssmDouble* palpha2,Gauss* gauss);
		void  GetAlpha2Hydro(IssmDouble* palpha2,Gauss* gauss);
		void  GetAlpha2Josh(IssmDouble* palpha2,Gauss* gauss);
		void  GetAlpha2Sommers(IssmDouble* palpha2,Gauss* gauss);
		void  GetAlpha2Temp(IssmDouble* palpha2,Gauss* gauss);
		void  GetAlpha2Viscous(IssmDouble* palpha2,Gauss* gauss);
		void  GetAlpha2WaterLayer(IssmDouble* palpha2,Gauss* gauss);
		void  GetAlpha2Weertman(IssmDouble* palpha2,Gauss* gauss);
		void  GetAlpha2WeertmanTemp(IssmDouble* palpha2,Gauss* gauss);
};

#endif  /* _FRICTION_H_ */
