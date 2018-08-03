//MASSTRANSPORT class definition
//
//   Usage:
//      masstransport=new masstransport();

function masstransport (){
	//methods
	this.setdefaultparameters = function(){// {{{

		//Type of stabilization to use 0:nothing 1:artificial_diffusivity 3:Discontinuous Galerkin
		this.stabilization=1;

		//Factor applied to compute the penalties kappa=max(stiffness matrix)*10^penalty_factor
		this.penalty_factor=3;

		//Minimum ice thickness that can be used
		this.min_thickness=1;

		//Hydrostatic adjustment
		this.hydrostatic_adjustment='Absolute';

		//default output
		this.requested_outputs=['default'];

	}// }}}
	this.disp= function(){// {{{
		console.log(sprintf('   Masstransport solution parameters:'));
		fielddisplay(this,'spcthickness','thickness constraints (NaN means no constraint) [m]');
		fielddisplay(this,'isfreesurface','do we use free surfaces (FS only) are mass conservation');
		fielddisplay(this,'min_thickness','minimum ice thickness allowed [m]');
		fielddisplay(this,'hydrostatic_adjustment',"adjustment of ice shelves surface and bed elevations: 'Incremental' or 'Absolute' ");
		fielddisplay(this,'stabilization','0: no, 1: artificial_diffusivity, 2: streamline upwinding, 3: discontinuous Galerkin, 4: Flux Correction Transport');

		console.log(sprintf('\n      %s','Penalty options:'));
		fielddisplay(this,'penalty_factor','offset used by penalties: penalty = Kmax*10^offset');
		fielddisplay(this,'vertex_pairing','pairs of vertices that are penalized');
		fielddisplay(this,'requested_outputs','additional outputs requested');

	}// }}}
	this.classname= function(){// {{{
		return "masstransport";
	}// }}}
    this.extrude = function(md) {//{{{
        this.spcthickness=project3d(md,'vector',this.spcthickness,'type','node');
        return this;
    }//}}}
		this.checkconsistency = function (md,solution,analyses){  // {{{

			//Early return
			if(!ArrayAnyEqual(ArrayIsMember('HydrologyShreveAnalysis',analyses),1) | (solution=='TransientSolution' & md.trans.ismasstransport==0)) return; 

			checkfield(md,'fieldname','masstransport.spcthickness','Inf',1,'timeseries',1);
			checkfield(md,'fieldname','masstransport.isfreesurface','values',[0 ,1]);
			checkfield(md,'fieldname','masstransport.hydrostatic_adjustment','values',['Absolute', 'Incremental']);
			checkfield(md,'fieldname','masstransport.stabilization','values',[0,1,2,3,4]);
			checkfield(md,'fieldname','masstransport.min_thickness','>',0);
			checkfield(md,'fieldname','masstransport.requested_outputs','stringrow',1);

		} // }}}
		this.marshall=function(md,prefix,fid) { //{{{

			var yts=md.constants.yts;

			WriteData(fid,prefix,'object',this,'fieldname','spcthickness','format','DoubleMat','mattype',1,'timeserieslength',md.mesh.numberofvertices+1,'yts',md.constants.yts);
			WriteData(fid,prefix,'object',this,'fieldname','isfreesurface','format','Boolean');
			WriteData(fid,prefix,'object',this,'fieldname','min_thickness','format','Double');
			WriteData(fid,prefix,'data',this.hydrostatic_adjustment,'format','String','name','md.masstransport.hydrostatic_adjustment');
			WriteData(fid,prefix,'object',this,'fieldname','stabilization','format','Integer');
			WriteData(fid,prefix,'object',this,'fieldname','vertex_pairing','format','DoubleMat','mattype',3);
			WriteData(fid,prefix,'object',this,'fieldname','penalty_factor','format','Double');

			//process requested outputs
			var outputs = this.requested_outputs;
			for (var i=0;i<outputs.length;i++){
				if (outputs[i] == 'default') {
					outputs.splice(i,1);
					var newoutputs=this.defaultoutputs(md);
					for (var j=0;j<newoutputs.length;j++) outputs.push(newoutputs[j]);
				}
			}
			WriteData(fid,prefix,'data',outputs,'name','md.masstransport.requested_outputs','format','StringArray');
		}//}}}
		this.defaultoutputs = function(md) { //{{{
			return ['Thickness','Surface','Base'];
		}//}}}
		this.fix=function() { //{{{
			this.spcthickness=NullFix(this.spcthickness,NaN);
			this.vertex_pairing=NullFix(this.vertex_pairing,NaN);
		}//}}}
	//properties 
	// {{{

	this.spcthickness           = NaN;
	this.isfreesurface          = 0;
	this.min_thickness          = 0;
	this.hydrostatic_adjustment = 0;
	this.stabilization          = 0;
	this.vertex_pairing         = NaN;
	this.penalty_factor         = 0;
	this.requested_outputs      = [];

	this.setdefaultparameters();
	//}}}
}
