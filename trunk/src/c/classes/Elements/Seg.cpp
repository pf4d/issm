/*!\file Seg.cpp
 * \brief: implementation of the Segment object
 */
/*Headers:*/
/*{{{*/
#ifdef HAVE_CONFIG_H
#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#include <stdio.h>
#include <string.h>
#include "../classes.h"
#include "../../shared/shared.h"
/*}}}*/

/*Element macros*/
#define NUMVERTICES 2
/*Constructors/destructor/copy*/
Seg::Seg(int seg_id, int seg_sid, int index, IoModel* iomodel,int nummodels)/*{{{*/
		:ElementHook(nummodels,index+1,NUMVERTICES,iomodel){

			/*id: */
			this->id  = seg_id;
			this->sid = seg_sid;

			//this->parameters: we still can't point to it, it may not even exist. Configure will handle this.
			this->parameters = NULL;

			/*intialize inputs: */
			this->inputs  = new Inputs();

			/*initialize pointers:*/
			this->nodes    = NULL;
			this->vertices = NULL;
			this->material = NULL;
			this->matpar   = NULL;

			/*Only allocate pointer*/
			this->element_type_list=xNew<int>(nummodels);
		}
/*}}}*/
Seg::~Seg(){/*{{{*/
	this->parameters=NULL;
}
/*}}}*/
Object* Seg::copy(){/*{{{*/

	int i;
	Seg* seg=NULL;

	seg=new Seg();

	//deal with TriaRef mother class
	int nanalyses = this->numanalyses;
	if(nanalyses > 0){
		seg->element_type_list=xNew<int>(nanalyses);
		for(i=0;i<nanalyses;i++){
			if (this->element_type_list[i]) seg->element_type_list[i]=this->element_type_list[i];
			else seg->element_type_list[i] = 0;
		}
	}
	else seg->element_type_list = NULL;
	seg->element_type=this->element_type;
	seg->numanalyses=nanalyses;

	//deal with ElementHook mother class
	if (this->hnodes){
		seg->hnodes=xNew<Hook*>(seg->numanalyses);
		for(i=0;i<seg->numanalyses;i++){
			if (this->hnodes[i]) seg->hnodes[i] = (Hook*)(this->hnodes[i]->copy());
			else seg->hnodes[i] = NULL;
		}
	}
	else seg->hnodes = NULL;

	seg->hvertices = (Hook*)this->hvertices->copy();
	seg->hmaterial = (Hook*)this->hmaterial->copy();
	seg->hmatpar   = (Hook*)this->hmatpar->copy();
	seg->hneighbors = NULL;

	/*deal with Element fields: */
	seg->id  = this->id;
	seg->sid = this->sid;
	if(this->inputs) seg->inputs = (Inputs*)(this->inputs->Copy());
	else seg->inputs=new Inputs();

	/*point parameters: */
	seg->parameters=this->parameters;

	/*recover objects: */
	if (this->nodes){
		unsigned int num_nodes = 3;
		seg->nodes = xNew<Node*>(num_nodes); //we cannot rely on an analysis_counter to tell us which analysis_type we are running, so we just copy the nodes.
		for(i=0;i<num_nodes;i++) if(this->nodes[i]) seg->nodes[i]=this->nodes[i]; else seg->nodes[i] = NULL;
	}
	else seg->nodes = NULL;

	seg->vertices = (Vertex**)this->hvertices->deliverp();
	seg->material = (Material*)this->hmaterial->delivers();
	seg->matpar   = (Matpar*)this->hmatpar->delivers();

	return seg;


}
/*}}}*/
void Seg::Marshall(char** pmarshalled_data,int* pmarshalled_data_size, int marshall_direction){ /*{{{*/

	MARSHALLING_ENUM(SegEnum);

	/*Call parent classes: */
	ElementHook::Marshall(pmarshalled_data,pmarshalled_data_size,marshall_direction);
	Element::MarshallElement(pmarshalled_data,pmarshalled_data_size,marshall_direction,this->numanalyses);
	SegRef::Marshall(pmarshalled_data,pmarshalled_data_size,marshall_direction);

	vertices = (Vertex**)this->hvertices->deliverp();
	material = (Material*)this->hmaterial->delivers();
	matpar   = (Matpar*)this->hmatpar->delivers();

}
/*}}}*/

IssmDouble Seg::CharacteristicLength(void){/*{{{*/

	IssmDouble xyz_list[NUMVERTICES][3];
	IssmDouble x1,y1,x2,y2;

	/*Get xyz list: */
	::GetVerticesCoordinates(&xyz_list[0][0],vertices,NUMVERTICES);
	x1=xyz_list[0][0]; y1=xyz_list[0][1];
	x2=xyz_list[1][0]; y2=xyz_list[1][1];

	return sqrt((x2-x1)*(x2-x1) + (y2-y1)*(y2-y1));
}
/*}}}*/
int        Seg::FiniteElement(void){/*{{{*/
	return this->element_type;
}
/*}}}*/
void       Seg::GetIcefrontCoordinates(IssmDouble** pxyz_front,IssmDouble* xyz_list,int levelsetenum){/*{{{*/
	
	/* Intermediaries */
	int nrfrontnodes,index;
	IssmDouble  levelset[NUMVERTICES];

	/*Recover parameters and values*/
	GetInputListOnVertices(&levelset[0],levelsetenum);
	/* Get nodes where there is no ice */
	nrfrontnodes=0;
	for(int i=0;i<NUMVERTICES;i++){
		if(levelset[i]>=0.){
			index=i;
			nrfrontnodes++;
		}
	}

	_assert_(nrfrontnodes==1);

	IssmDouble* xyz_front = xNew<IssmDouble>(3);

	/* Return nodes */
	for(int dir=0;dir<3;dir++){
		xyz_front[dir]=xyz_list[3*index+dir];
	}

	*pxyz_front=xyz_front;
}/*}}}*/
int        Seg::GetNumberOfNodes(void){/*{{{*/
	return this->NumberofNodes(this->element_type);
}
/*}}}*/
int        Seg::GetNumberOfVertices(void){/*{{{*/
	return NUMVERTICES;
}
/*}}}*/
void       Seg::GetVerticesCoordinates(IssmDouble** pxyz_list){/*{{{*/

	IssmDouble* xyz_list = xNew<IssmDouble>(NUMVERTICES*3);
	::GetVerticesCoordinates(xyz_list,this->vertices,NUMVERTICES);

	/*Assign output pointer*/
	*pxyz_list = xyz_list;

}/*}}}*/
bool       Seg::IsIcefront(void){/*{{{*/

	bool isicefront;
	int i,nrice;
	IssmDouble ls[NUMVERTICES];

	/*Retrieve all inputs and parameters*/
	GetInputListOnVertices(&ls[0],MaskIceLevelsetEnum);

	/* If only one vertex has ice, there is an ice front here */
	isicefront=false;
	if(IsIceInElement()){
		nrice=0;       
		for(i=0;i<NUMVERTICES;i++)
		 if(ls[i]<0.) nrice++;
		if(nrice==1) isicefront= true;
	}
	return isicefront;
}/*}}}*/
void       Seg::JacobianDeterminant(IssmDouble* pJdet,IssmDouble* xyz_list,Gauss* gauss){/*{{{*/

	_assert_(gauss->Enum()==GaussSegEnum);
	this->GetJacobianDeterminant(pJdet,xyz_list,(GaussSeg*)gauss);

}
/*}}}*/
void       Seg::JacobianDeterminantSurface(IssmDouble* pJdet,IssmDouble* xyz_list,Gauss* gauss){/*{{{*/

	*pJdet = 1.;

}
/*}}}*/
Gauss*     Seg::NewGauss(void){/*{{{*/
	return new GaussSeg();
}
/*}}}*/
Gauss*     Seg::NewGauss(int order){/*{{{*/
	return new GaussSeg(order);
}
/*}}}*/
Gauss*     Seg::NewGauss(IssmDouble* xyz_list, IssmDouble* xyz_list_front,int order){/*{{{*/

	/*Output*/
	Gauss* gauss = NULL;

	if(xyz_list_front[0] == xyz_list[0]){
		gauss = new GaussSeg(-1.);
	}
	else if(xyz_list_front[0] == xyz_list[3*1+0]){
		gauss = new GaussSeg(+1.);
	}
	else{
		_error_("front is not located on element edge");
	}

	return gauss;
}
/*}}}*/
void       Seg::NodalFunctions(IssmDouble* basis, Gauss* gauss){/*{{{*/

	_assert_(gauss->Enum()==GaussSegEnum);
	this->GetNodalFunctions(basis,(GaussSeg*)gauss,this->element_type);

}
/*}}}*/
void       Seg::NodalFunctionsDerivatives(IssmDouble* dbasis,IssmDouble* xyz_list,Gauss* gauss){/*{{{*/

	_assert_(gauss->Enum()==GaussSegEnum);
	this->GetNodalFunctionsDerivatives(dbasis,xyz_list,(GaussSeg*)gauss,this->element_type);

}
/*}}}*/
void       Seg::NodalFunctionsP1(IssmDouble* basis, Gauss* gauss){/*{{{*/

	_assert_(gauss->Enum()==GaussSegEnum);
	this->GetNodalFunctions(basis,(GaussSeg*)gauss,P1Enum);

}
/*}}}*/
void       Seg::NodalFunctionsP2(IssmDouble* basis, Gauss* gauss){/*{{{*/

	_assert_(gauss->Enum()==GaussSegEnum);
	this->GetNodalFunctions(basis,(GaussSeg*)gauss,P2Enum);

}
/*}}}*/
void       Seg::NormalSection(IssmDouble* normal,IssmDouble* xyz_list_front){/*{{{*/

	IssmDouble* xyz_list = xNew<IssmDouble>(NUMVERTICES*3);
	::GetVerticesCoordinates(xyz_list,this->vertices,NUMVERTICES);

	if(xyz_list_front[0]>xyz_list[0])
	 normal[0]= + 1.;
	else
	 normal[0]= - 1.;

	xDelete<IssmDouble>(xyz_list);
}
/*}}}*/
int        Seg::ObjectEnum(void){/*{{{*/

	return SegEnum;

}
/*}}}*/
