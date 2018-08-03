/*
 * \file DataSet.cpp
 * \brief: Implementation of DataSet class
 */

/*Headers: {{{*/
#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#include <cstring>
#include <vector>
#include <functional>
#include <algorithm>
#include <iostream>

#include "../datastructures/datastructures.h"
#include "../shared/shared.h"
#include "../classes/classes.h"

using namespace std;
/*}}}*/

/*Constructors/Destructors*/
DataSet::DataSet(){/*{{{*/

	sorted=0;
	numsorted=0;
	presorted=0;
	enum_type=-1;
	sorted_ids=NULL;
	id_offsets=NULL;

}
/*}}}*/
DataSet::DataSet(int dataset_enum){/*{{{*/
	enum_type=dataset_enum;

	sorted=0;
	numsorted=0;
	presorted=0;
	sorted_ids=NULL;
	id_offsets=NULL;

}
/*}}}*/
DataSet* DataSet::Copy(void){/*{{{*/

	vector<Object*>::iterator obj;
	Object* object_copy=NULL;

	DataSet* copy=new DataSet(this->enum_type);

	copy->sorted=this->sorted;
	copy->numsorted=this->numsorted;
	copy->presorted=this->presorted;

	/*Now we need to deep copy the objects: */
	for ( obj=this->objects.begin() ; obj < this->objects.end(); obj++ ){
		/*Call copy on object: */
		object_copy = (*obj)->copy();
		copy->AddObject(object_copy);
	}

	/*Build id_offsets and sorted_ids*/
	int objsize = this->numsorted;
	if(this->sorted && objsize>0 && this->id_offsets){	
		/*Allocate new ids*/
		copy->id_offsets=xNew<int>(objsize);
		xMemCpy<int>(copy->id_offsets,this->id_offsets,objsize);
	}
	else copy->id_offsets=NULL;
	if(this->sorted && objsize>0 && this->sorted_ids){
		/*Allocate new ids*/
		copy->sorted_ids=xNew<int>(objsize);
		xMemCpy<int>(copy->sorted_ids,this->sorted_ids,objsize);
	}
	else copy->sorted_ids=NULL;

	return copy;
}
/*}}}*/
DataSet::~DataSet(){/*{{{*/
	clear();
	xDelete<int>(sorted_ids);
	xDelete<int>(id_offsets);
}
/*}}}*/

/*Specific methods*/
void  DataSet::Marshall(char** pmarshalled_data, int* pmarshalled_data_size, int marshall_direction){ /*{{{*/
	
	vector<Object*>::iterator obj;
	int obj_size=0;
	int obj_enum=0;
	int i;

	if(marshall_direction==MARSHALLING_FORWARD || marshall_direction==MARSHALLING_SIZE){
		obj_size=objects.size();
	}
	else{
		clear();
	}

	MARSHALLING_ENUM(DataSetEnum);
	MARSHALLING(enum_type);
	MARSHALLING(sorted);
	MARSHALLING(presorted);
	MARSHALLING(numsorted);

	/*Now branch according to direction of marshalling: */
	if(marshall_direction==MARSHALLING_FORWARD || marshall_direction==MARSHALLING_SIZE){
		if(!(this->sorted && numsorted>0 && this->id_offsets)){
			sorted_ids=NULL;
			id_offsets=NULL;
		  }
		MARSHALLING_DYNAMIC(sorted_ids,int,numsorted);
		MARSHALLING_DYNAMIC(id_offsets,int,numsorted);
		MARSHALLING(obj_size);

		/*Go through our objects, and marshall them into the buffer: */
		for( obj=this->objects.begin() ; obj < this->objects.end(); obj++ ){
			obj_enum=(*obj)->ObjectEnum();
			MARSHALLING(obj_enum);
			(*obj)->Marshall(pmarshalled_data,pmarshalled_data_size,marshall_direction);
		}
	}
	else{

		MARSHALLING_DYNAMIC(sorted_ids,int,numsorted);
		MARSHALLING_DYNAMIC(id_offsets,int,numsorted);
		if (!(this->sorted && numsorted>0)){
		 sorted_ids=NULL;
		 id_offsets=NULL;
		}
		MARSHALLING(obj_size);

		/*This is the heart of the demashalling method. We have a buffer coming
		 in, and we are supposed to create a dataset out of it. No such thing
		 as class orientation for buffers, we need to key off the enum of each
		 object stored in the buffer. */
		for(i=0;i<obj_size;i++){

			/*Recover enum of object first: */
			MARSHALLING(obj_enum); 

			/*Giant case statement to spin-up the right object, and demarshall into it the information 
			 *stored in the buffer: */
			if(obj_enum==NodeEnum){
				Node* node=NULL;
				node=new Node();
				node->Marshall(pmarshalled_data,pmarshalled_data_size,marshall_direction);
				this->AddObject(node);
			}
			else if(obj_enum==VertexEnum){
				Vertex* vertex=NULL;
				vertex=new Vertex();
				vertex->Marshall(pmarshalled_data,pmarshalled_data_size,marshall_direction);
				this->AddObject(vertex);
			}
			else if(obj_enum==MaticeEnum){
				Matice* matice=NULL;
				matice=new Matice();
				matice->Marshall(pmarshalled_data,pmarshalled_data_size,marshall_direction);
				this->AddObject(matice);
			}
			else if(obj_enum==MatestarEnum){
				Matestar* matestar=NULL;
				matestar=new Matestar();
				matestar->Marshall(pmarshalled_data,pmarshalled_data_size,marshall_direction);
				this->AddObject(matestar);
			}
			else if(obj_enum==MatparEnum){
				Matpar* matpar=NULL;
				matpar=new Matpar();
				matpar->Marshall(pmarshalled_data,pmarshalled_data_size,marshall_direction);
				this->AddObject(matpar);
			}
			else if(obj_enum==SpcStaticEnum){
				SpcStatic* spcstatic=NULL;
				spcstatic=new SpcStatic();
				spcstatic->Marshall(pmarshalled_data,pmarshalled_data_size,marshall_direction);
				this->AddObject(spcstatic);
			}
			else if(obj_enum==SpcDynamicEnum){
				SpcDynamic* spcdynamic=NULL;
				spcdynamic=new SpcDynamic();
				spcdynamic->Marshall(pmarshalled_data,pmarshalled_data_size,marshall_direction);
				this->AddObject(spcdynamic);
			}
			else if(obj_enum==SpcTransientEnum){
				SpcTransient* spctransient=NULL;
				spctransient=new SpcTransient();
				spctransient->Marshall(pmarshalled_data,pmarshalled_data_size,marshall_direction);
				this->AddObject(spctransient);
			}
			else if(obj_enum==TriaEnum){
				Tria* tria=NULL;
				tria=new Tria();
				tria->Marshall(pmarshalled_data,pmarshalled_data_size,marshall_direction);
				this->AddObject(tria);
			}
			else if(obj_enum==PentaEnum){
				Penta* penta=NULL;
				penta=new Penta();
				penta->Marshall(pmarshalled_data,pmarshalled_data_size,marshall_direction);
				this->AddObject(penta);
			}
			else if(obj_enum==TetraEnum){
				Tetra* tetra=NULL;
				tetra=new Tetra();
				tetra->Marshall(pmarshalled_data,pmarshalled_data_size,marshall_direction);
				this->AddObject(tetra);
			}
			else if(obj_enum==SegEnum){
				Seg* seg=NULL;
				seg=new Seg();
				seg->Marshall(pmarshalled_data,pmarshalled_data_size,marshall_direction);
				this->AddObject(seg);
			}
			else if(obj_enum==BoolInputEnum){
				BoolInput* boolinput=NULL;
				boolinput=new BoolInput();
				boolinput->Marshall(pmarshalled_data,pmarshalled_data_size,marshall_direction);
				this->AddObject(boolinput);
			}
			else if(obj_enum==DoubleInputEnum){
				DoubleInput* doubleinput=NULL;
				doubleinput=new DoubleInput();
				doubleinput->Marshall(pmarshalled_data,pmarshalled_data_size,marshall_direction);
				this->AddObject(doubleinput);
			}
			else if(obj_enum==IntInputEnum){
				IntInput* intinput=NULL;
				intinput=new IntInput();
				intinput->Marshall(pmarshalled_data,pmarshalled_data_size,marshall_direction);
				this->AddObject(intinput);
			}
			else if(obj_enum==ControlInputEnum){
				ControlInput* cinput=NULL;
				cinput=new ControlInput();
				cinput->Marshall(pmarshalled_data,pmarshalled_data_size,marshall_direction);
				this->AddObject(cinput);
			}
			else if(obj_enum==TransientInputEnum){
				TransientInput* transinput=NULL;
				transinput=new TransientInput();
				transinput->Marshall(pmarshalled_data,pmarshalled_data_size,marshall_direction);
				this->AddObject(transinput);
			}
			else if(obj_enum==TriaInputEnum){
				TriaInput* triainput=NULL;
				triainput=new TriaInput();
				triainput->Marshall(pmarshalled_data,pmarshalled_data_size,marshall_direction);
				this->AddObject(triainput);
			}
			else if(obj_enum==PentaInputEnum){
				PentaInput* pentainput=NULL;
				pentainput=new PentaInput();
				pentainput->Marshall(pmarshalled_data,pmarshalled_data_size,marshall_direction);
				this->AddObject(pentainput);
			}
			else if(obj_enum==TetraInputEnum){
				TetraInput* tetrainput=NULL;
				tetrainput=new TetraInput();
				tetrainput->Marshall(pmarshalled_data,pmarshalled_data_size,marshall_direction);
				this->AddObject(tetrainput);
			}
			else if(obj_enum==SegInputEnum){
				SegInput* seginput=NULL;
				seginput=new SegInput();
				seginput->Marshall(pmarshalled_data,pmarshalled_data_size,marshall_direction);
				this->AddObject(seginput);
			}
			else if(obj_enum==RiftfrontEnum){
				Riftfront* rift=NULL;
				rift=new Riftfront();
				rift->Marshall(pmarshalled_data,pmarshalled_data_size,marshall_direction);
				this->AddObject(rift);
			}
			else if(obj_enum==NumericalfluxEnum){
				Numericalflux* numflux=NULL;
				numflux=new Numericalflux();
				numflux->Marshall(pmarshalled_data,pmarshalled_data_size,marshall_direction);
				this->AddObject(numflux);
			}
			else if(obj_enum==PengridEnum){
				Pengrid* pengrid=NULL;
				pengrid=new Pengrid();
				pengrid->Marshall(pmarshalled_data,pmarshalled_data_size,marshall_direction);
				this->AddObject(pengrid);
			}
			else if(obj_enum==PenpairEnum){
				Penpair* penpair=NULL;
				penpair=new Penpair();
				penpair->Marshall(pmarshalled_data,pmarshalled_data_size,marshall_direction);
				this->AddObject(penpair);
			}
			else _error_("could not recognize enum type: " << obj_enum << ": " << EnumToStringx(obj_enum) ); 
		}
	}
}
/*}}}*/
int   DataSet::AddObject(Object* object){/*{{{*/

	_assert_(this);
	objects.push_back(object);

	return 1;
}
/*}}}*/
void  DataSet::clear(){/*{{{*/

/*  use reverse_iterator for efficiency in matlab memory manager
	(keeping old code in case it needs to revert back)  */

//	vector<Object*>::iterator object;
	vector<Object*>::reverse_iterator object;

//	for ( object=objects.begin() ; object < objects.end(); object++ ){
//		delete (*object);
//	}
	for ( object=objects.rbegin() ; object < objects.rend(); object++ ){
		delete (*object);
	}
	objects.clear();
}
/*}}}*/
int   DataSet::DeleteObject(Object* object){/*{{{*/

	vector<Object*>::iterator iterator;

	if(object){
		iterator = find(objects.begin(), objects.end(),object);
		delete *iterator;
		objects.erase(iterator);
	}

	return 1;

}
/*}}}*/
void  DataSet::DeepEcho(){/*{{{*/

	vector<Object*>::iterator object;

	_assert_(this);

	_printf0_("DataSet echo: " << objects.size() << " objects\n");

	for ( object=objects.begin() ; object < objects.end(); object++ ){

		/*Call deep echo on object: */
		(*object)->DeepEcho();

	}
}
/*}}}*/
void  DataSet::Echo(){/*{{{*/

	vector<Object*>::iterator object;

	_assert_(this);

	_printf0_("DataSet echo: " << objects.size() << " objects\n");

	for ( object=objects.begin() ; object < objects.end(); object++ ){

		/*Call echo on object: */
		(*object)->Echo();

	}
	return;
}
/*}}}*/
int   DataSet::GetEnum(){/*{{{*/
	return enum_type;
}
/*}}}*/
int   DataSet::GetEnum(int offset){/*{{{*/

	return objects[offset]->ObjectEnum();

}
/*}}}*/
Object* DataSet::GetObjectByOffset(int offset){/*{{{*/

	/*Check index in debugging mode*/
	_assert_(this!=NULL);
	_assert_(offset>=0);
	_assert_(offset<this->Size());

	return objects[offset];

}
/*}}}*/
Object* DataSet::GetObjectById(int* poffset,int eid){/*{{{*/

	int id_offset;
	int offset;

	_assert_(this);
	if(!sorted || objects.size()>numsorted)_error_("trying to binary search on a non-sorted dataset!");

	/*Carry out a binary search on the sorted_ids: */
	if(!binary_search(&id_offset,eid,sorted_ids,objects.size())){
		_error_("could not find object with id " << eid << " in DataSet " << EnumToStringx(enum_type));
	}

	/*Convert  the id offset into sorted offset: */
	offset=id_offsets[id_offset];

	/*Assign output pointers if requested:*/
	if(poffset)*poffset=offset;

	/*Return object at offset position in objects :*/
	return objects[offset];
}
/*}}}*/
void  DataSet::Presort(){/*{{{*/

	/*vector of objects is already sorted, just allocate the sorted ids and their
	 * offsets:*/
	if(objects.size()){

		/*Delete existing ids*/
		if(sorted_ids) xDelete<int>(sorted_ids);
		if(id_offsets) xDelete<int>(id_offsets);

		/*Allocate new ids*/
		sorted_ids=xNew<int>(objects.size());
		id_offsets=xNew<int>(objects.size());

		/*Build id_offsets and sorted_ids*/
		for(int i=0;i<objects.size();i++){
			id_offsets[i]=i;
			sorted_ids[i]=objects[i]->Id();
		}
	}

	/*set sorted flag: */
	numsorted=objects.size();
	sorted=1;
}
/*}}}*/
int   DataSet::Size(void){/*{{{*/
	_assert_(this!=NULL);

	return objects.size();
}
/*}}}*/
void  DataSet::Sort(){/*{{{*/

	/*Only sort if we are not already sorted: */
	if(!sorted){
		_error_("not implemented yet!");
	}
}
/*}}}*/
