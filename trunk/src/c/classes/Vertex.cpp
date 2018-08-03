/*!\file Vertex.c
 * \brief: implementation of the Vertex object
 */

/*Include files: {{{*/
#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#include <string.h>
#include "classes.h"
#include "shared/shared.h"
/*}}}*/

/*Vertex constructors and destructor:*/
Vertex::Vertex(){/*{{{*/
	return;
}
/*}}}*/
Vertex::Vertex(int vertex_id, int vertex_sid,int i, IoModel* iomodel){/*{{{*/

	this->id           = vertex_id;
	this->sid          = vertex_sid;
	this->pid          = UNDEF;

	_assert_(iomodel->Data("md.mesh.x") && iomodel->Data("md.mesh.y") && iomodel->Data("md.mesh.z"));
	this->x            = iomodel->Data("md.mesh.x")[i];
	this->y            = iomodel->Data("md.mesh.y")[i];
	this->z            = iomodel->Data("md.mesh.z")[i];
	this->domaintype     = iomodel->domaintype;

	switch(iomodel->domaintype){
		case Domain3DEnum:
			_assert_(iomodel->Data("md.geometry.base") && iomodel->Data("md.geometry.thickness"));
			this->sigma = (iomodel->Data("md.mesh.z")[i]-iomodel->Data("md.geometry.base")[i])/(iomodel->Data("md.geometry.thickness")[i]);
			break;
		case Domain3DsurfaceEnum:
			this->latitute     = iomodel->Data("md.mesh.lat")[i];
			this->longitude    = iomodel->Data("md.mesh.long")[i];
			this->R            = iomodel->Data("md.mesh.r")[i];
			break;
		case Domain2DhorizontalEnum:
			this->sigma = 0.;
			break;
		case Domain2DverticalEnum:
			_assert_(iomodel->Data("md.geometry.base") && iomodel->Data("md.geometry.thickness"));
			this->sigma = (iomodel->Data("md.mesh.y")[i]-iomodel->Data("md.geometry.base")[i])/(iomodel->Data("md.geometry.thickness")[i]);
			break;
	}

	_assert_(iomodel->numbernodetoelementconnectivity);
	this->connectivity = iomodel->numbernodetoelementconnectivity[i];

}
/*}}}*/
Vertex::~Vertex(){/*{{{*/
	return;
}
/*}}}*/

/*Object virtual functions definitions:*/
Object* Vertex::copy() {/*{{{*/

	return new Vertex(*this); 

}
/*}}}*/
void Vertex::DeepEcho(void){/*{{{*/
	this->Echo();
}
/*}}}*/
void Vertex::Echo(void){/*{{{*/

	_printf_("Vertex:\n");
	_printf_("   id: " << id << "\n");
	_printf_("   sid: " << sid << "\n");
	_printf_("   pid: " << pid << "\n");
	_printf_("   x: " << x << "\n");
	_printf_("   y: " << y << "\n");
	_printf_("   z: " << z << "\n");
	_printf_("   sigma: " << sigma << "\n");
	_printf_("   connectivity: " << connectivity << "\n");
	_printf_("   clone: " << clone << "\n");

	return;
}
/*}}}*/
int Vertex::Id(void){ return id; }/*{{{*/
/*}}}*/
void Vertex::Marshall(char** pmarshalled_data,int* pmarshalled_data_size, int marshall_direction){ /*{{{*/

	MARSHALLING_ENUM(VertexEnum);
	MARSHALLING(clone);
	MARSHALLING(domaintype);
	MARSHALLING(id);
	MARSHALLING(sid);
	MARSHALLING(pid);
	MARSHALLING(x);
	MARSHALLING(y);
	MARSHALLING(z);
	MARSHALLING(sigma);
	MARSHALLING(connectivity);

}
/*}}}*/
int Vertex::ObjectEnum(void){/*{{{*/

	return VertexEnum;

}
/*}}}*/

/*Vertex management: */
int        Vertex::Connectivity(void){return connectivity;}/*{{{*/
/*}}}*/
void       Vertex::DistributePids(int* ppidcount){/*{{{*/

	/*retrieve current pid*/
	int pidcount=*ppidcount;

	/*This vertex is a clone! Don't distribute pids, it will get them from another cpu!*/
	if(this->clone) return;

	/*This vertex should distribute its pid*/
	this->pid=pidcount;
	pidcount++;

	/*Assign output pointers: */
	*ppidcount=pidcount;
}
/*}}}*/
IssmDouble Vertex::GetLatitude(){/*{{{*/
	return this->latitute;
}
/*}}}*/
IssmDouble Vertex::GetLongitude(){/*{{{*/
	return this->longitude;
}
/*}}}*/
IssmDouble Vertex::GetRadius(){/*{{{*/
	return this->R;
}
/*}}}*/
IssmDouble Vertex::GetX(){/*{{{*/
	return this->x;
}
/*}}}*/
IssmDouble Vertex::GetY(){/*{{{*/
	return this->y;
}
/*}}}*/
IssmDouble Vertex::GetZ(){/*{{{*/
	return this->z;
}
/*}}}*/
void       Vertex::OffsetPids(int pidcount){/*{{{*/

	/*This vertex is a clone, don't offset the pids*/
	if(this->clone) return;

	/*This vertex should offset his pid, go ahead: */
	this->pid+=pidcount;
}
/*}}}*/
int        Vertex::Pid(void){ return pid; }/*{{{*/
/*}}}*/
void       Vertex::SetClone(int* minranks){/*{{{*/

	int my_rank;

	/*recover my_rank:*/
	my_rank=IssmComm::GetRank();

	if (minranks[this->sid]==my_rank){
		this->clone=false;
	}
	else{
		/*!there is a cpu with lower rank that has the same vertex, 
		therefore, I am a clone*/
		this->clone=true;
	}

}
/*}}}*/
void       Vertex::ShowTruePids(int* truepids){/*{{{*/

	/*Are we a clone? : */
	if(this->clone)return;

	/*Ok, we are not a clone, just plug our pid into truepids: */
	truepids[this->sid]=this->pid;
}
/*}}}*/
int        Vertex::Sid(void){ return sid; }/*{{{*/
/*}}}*/
void       Vertex::ToXYZ(Matrix<IssmDouble>* matrix){/*{{{*/

	IssmDouble xyz[3];
	int        indices[3];

	if (this->clone==true) return;

	xyz[0]=x;
	xyz[1]=y; 
	xyz[2]=z;
	indices[0]=0;
	indices[1]=1; 
	indices[2]=2;

	matrix->SetValues(1,&sid,3,&indices[0],&xyz[0],INS_VAL);
}
/*}}}*/
void       Vertex::UpdateClonePids(int* alltruepids){/*{{{*/

	/*If we are not a clone, don't update, we already have pids: */
	if(!this->clone)return;

	/*Ok, we are a clone node, but we did not create the pid for this vertex 
	 * Therefore, our pid is garbage right now. Go pick it up in the alltruepids: */
	this->pid=alltruepids[this->sid];
}
/*}}}*/
void       Vertex::UpdatePosition(Vector<IssmDouble>* vx,Vector<IssmDouble>* vy,Vector<IssmDouble>* vz,Parameters* parameters,IssmDouble* surface,IssmDouble* bed){/*{{{*/

	IssmDouble oldy,newy,vely;
	IssmDouble oldz,newz,velz;
	IssmDouble dt;

	/*Get time stepping*/
	parameters->FindParam(&dt,TimesteppingTimeStepEnum);

	/*sigma remains constant. z=bed+sigma*thickness*/
	switch(this->domaintype){
		case Domain2DhorizontalEnum:
			/*Nothing*/
			return;
		case Domain2DverticalEnum:
			oldy = this->y;
			newy = bed[this->pid]+sigma*(surface[this->pid] - bed[this->pid]);
			vely = (newy-oldy)/dt;
			this->y = newy;
			vy->SetValue(this->pid,vely,INS_VAL);
			_assert_(!xIsNan<IssmDouble>(vely));
			return;
		case Domain3DEnum:
			oldz = this->z;
			newz = bed[this->pid]+sigma*(surface[this->pid] - bed[this->pid]);
			velz = (newz-oldz)/dt;
			this->z = newz;
			vz->SetValue(this->pid,velz,INS_VAL);
			_assert_(!xIsNan<IssmDouble>(velz));
			return;
		default:
			_error_("not implemented");
	}
}
/*}}}*/
void       Vertex::VertexCoordinates(Vector<IssmDouble>* vx,Vector<IssmDouble>* vy,Vector<IssmDouble>* vz, bool spherical){/*{{{*/

	if (this->clone==true) return;

	if(!spherical){
		vx->SetValue(this->sid,this->x,INS_VAL);
		vy->SetValue(this->sid,this->y,INS_VAL);
		vz->SetValue(this->sid,this->z,INS_VAL);
	}
	else{
		vx->SetValue(this->sid,this->latitute,INS_VAL);
		vy->SetValue(this->sid,this->longitude,INS_VAL);
		vz->SetValue(this->sid,this->R,INS_VAL);
	}

	return;
}
/*}}}*/

/*Methods relating to Vertex, but not internal methods: */
void GetVerticesCoordinates(IssmDouble* xyz,Vertex** vertices, int numvertices,bool spherical){ /*{{{*/

	_assert_(vertices);
	_assert_(xyz);

	if(!spherical){
		for(int i=0;i<numvertices;i++) {
			xyz[i*3+0]=vertices[i]->GetX();
			xyz[i*3+1]=vertices[i]->GetY();
			xyz[i*3+2]=vertices[i]->GetZ();
		}
	}
	else{
		for(int i=0;i<numvertices;i++) {
			xyz[i*3+0]=vertices[i]->GetLatitude();
			xyz[i*3+1]=vertices[i]->GetLongitude();
			xyz[i*3+2]=vertices[i]->GetRadius();
		}
	}
}/*}}}*/
