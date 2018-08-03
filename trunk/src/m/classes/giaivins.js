//GIA class definition for Ivins and James model 
//
//   Usage:
//      giaivins=new giaivins();

function giaivins (){
	//methods
	this.setdefaultparameters = function(){// {{{

		this.cross_section_shape=1; //square as default (see iedge in GiaDeflectionCorex)
	
	}// }}}
	this.disp= function(){// {{{

		console.log(sprintf('   giaivins parameters:'));

		fielddisplay(this,'mantle_viscosity','mantle viscosity[Pa s]');
		fielddisplay(this,'lithosphere_thickness','lithosphere thickness (km)');
		fielddisplay(this,'cross_section_shape','1: square-edged (default). 2: elliptical.  See iedge in GiaDeflectionCore');

	}// }}}
	this.classname= function(){// {{{
		return "giaivins";
	}// }}}
	this.checkconsistency = function(md,solution,analyses) { // {{{

		if(!ArrayAnyEqual(ArrayIsMember('GiaAnalysis',analyses),1))return;

		checkfield(md,'fieldname','gia.mantle_viscosity','NaN',1,'Inf',1,'size',[md.mesh.numberofvertices,1],'>',0);
		checkfield(md,'fieldname','gia.lithosphere_thickness','NaN',1,'Inf',1,'size',[md.mesh.numberofvertices,1],'>',0);
		checkfield(md,'fieldname','gia.cross_section_shape','numel',[1],'values',[1,2]);

		//be sure that if we are running a masstransport ice flow model coupled with giaivins, that thickness forcings 
		//are not provided into the future.
		if (solution=='TransientSolution' & md.trans.ismasstransport & md.trans.isgia){
			//figure out if thickness is a transient forcing: 
			if (md.geometry.thickness.length == (md.mesh.numberofvertices+1)){
				//recover the furthest time "in time": 
				t=md.geometry.thickness[0].length;
				if(md.geometry.thickness[md.geometry.thickness.length-1][t-1]!=md.timestepping.start_time){
					md.checkmessage('if ismasstransport is on, transient thickness forcing for the giaivins model should not be provided in the future. Synchronize your start_time to correspond to the most recent transient thickness forcing timestep');
				}
			}
		}
	} // }}}
		this.marshall=function(md,prefix,fid) { //{{{
			WriteData(fid,prefix,'object',this,'fieldname','mantle_viscosity','format','DoubleMat','mattype',1);
			WriteData(fid,prefix,'object',this,'fieldname','lithosphere_thickness','format','DoubleMat','mattype',1,'scale',Math.pow(10,3)); //from km to m
			WriteData(fid,prefix,'object',this,'fieldname','cross_section_shape','format','Integer');
		}//}}}
		this.fix=function() { //{{{
			this.mantle_viscosity=NullFix(this.mantle_viscosity,NaN);
			this.lithosphere_thickness=NullFix(this.lithosphere_thickness,NaN);
		}//}}}
	//properties 
	// {{{

	this.mantle_viscosity              = NaN;
	this.lithosphere_thickness         = NaN;
	this.cross_section_shape           = 0;

	this.setdefaultparameters();
	//}}}
}
