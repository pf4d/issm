/*!\file GaussSeg.c
 * \brief: implementation of the GaussSeg object
 */

#include "./GaussSeg.h"
#include "../../shared/io/Print/Print.h"
#include "../../shared/Exceptions/exceptions.h"
#include "../../shared/MemOps/MemOps.h"
#include "../../shared/Enum/Enum.h"
#include "../../shared/Numerics/GaussPoints.h"
#include "../../shared/Numerics/constants.h"

/*GaussSeg constructors and destructors:*/
GaussSeg::GaussSeg(){/*{{{*/

	numgauss=-1;

	weights=NULL;
	coords1=NULL;

	weight=UNDEF;
	coord1=UNDEF;
}
/*}}}*/
GaussSeg::GaussSeg(int order){/*{{{*/

	IssmPDouble* pcoords1=NULL;
	IssmPDouble* pweights=NULL;

	/*Get gauss points*/
	this->numgauss = order;
	GaussLegendreLinear(&pcoords1,&pweights,order);
	
	this->coords1=xNew<IssmDouble>(numgauss);
	this->weights=xNew<IssmDouble>(numgauss);

	/*cast : */
	for(int i=0;i<numgauss;i++){
		this->coords1[i]=pcoords1[i];
		this->weights[i]=pweights[i];
	}

	/*Free ressources: */
	xDelete<IssmPDouble>(pcoords1);
	xDelete<IssmPDouble>(pweights);

	/*Initialize static fields as undefinite*/
	weight=UNDEF;
	coord1=UNDEF;
}
/*}}}*/
GaussSeg::GaussSeg(IssmDouble position){/*{{{*/

	/*Get gauss points*/
	this->numgauss = 1;
	this->coords1=xNew<IssmDouble>(numgauss);
	this->weights=xNew<IssmDouble>(numgauss);

	/*cast : */
	_assert_(position>=-1. && position<=+1.);
	this->coords1[0]=position;
	this->weights[0]=1.;

	/*Initialize static fields as undefinite*/
	weight=UNDEF;
	coord1=UNDEF;
}
/*}}}*/
GaussSeg::~GaussSeg(){/*{{{*/
	xDelete<IssmDouble>(weights);
	xDelete<IssmDouble>(coords1);
}
/*}}}*/

/*Methods*/
int GaussSeg::begin(void){/*{{{*/

	/*Check that this has been initialized*/
	_assert_(numgauss>0);
	_assert_(weights);
	_assert_(coords1);

	/*return first gauss index*/
	return 0;
}
/*}}}*/
void GaussSeg::Echo(void){/*{{{*/

	_printf_("GaussSeg:\n");
	_printf_("   numgauss: " << numgauss << "\n");

	if (weights){
	 _printf_("   weights = ["); 
	 for(int i=0;i<numgauss;i++) _printf_(" " << weights[i] << "\n");
	 _printf_("]\n");
	}
	else _printf_("weights = NULL\n");
	if (coords1){
	 _printf_("   coords1 = ["); 
	 for(int i=0;i<numgauss;i++) _printf_(" " << coords1[i] << "\n");
	 _printf_("]\n");
	}
	_printf_("   weight = " << weight << "\n");
	_printf_("   coord1 = " << coord1 << "\n");

}
/*}}}*/
int GaussSeg::end(void){/*{{{*/

	/*Check that this has been initialized*/
	_assert_(numgauss>0);
	_assert_(weights);
	_assert_(coords1);

	/*return last gauss index +1*/
	return numgauss;
}
/*}}}*/
int GaussSeg::Enum(void){/*{{{*/
	return GaussSegEnum;
}
/*}}}*/
void GaussSeg::GaussPoint(int ig){/*{{{*/

	/*Check input in debugging mode*/
	 _assert_(ig>=0 && ig< numgauss);

	 /*update static arrays*/
	 weight=weights[ig];
	 coord1=coords1[ig];
}
/*}}}*/
void GaussSeg::GaussNode(int finiteelement,int iv){/*{{{*/

	/*in debugging mode: check that the default constructor has been called*/
	_assert_(numgauss==-1);

	/*update static arrays*/
	switch(finiteelement){
		case P1Enum: case P1DGEnum:
			switch(iv){
				case 0: coord1=-1.; break;
				case 1: coord1=+1.; break;
				default: _error_("node index should be in [0 1]");
			}
			break;
		default: _error_("Finite element "<<EnumToStringx(finiteelement)<<" not supported");
	}

}
/*}}}*/
void GaussSeg::GaussVertex(int iv){/*{{{*/

	/*in debugging mode: check that the default constructor has been called*/
	_assert_(numgauss==-1);

	/*update static arrays*/
	switch(iv){
		case 0: coord1=-1.; break;
		case 1: coord1=+1.; break;
		default: _error_("vertex index should be in [0 1]");
	}
}
/*}}}*/
void GaussSeg::SynchronizeGaussBase(Gauss* gauss){/*{{{*/

	_error_("not supported");
}
/*}}}*/
