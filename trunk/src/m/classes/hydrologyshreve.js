//HYDROLOGYSHREVE class definition
//
//   Usage:
//      hydrologyshreve=new hydrologyshreve();

function hydrologyshreve (){
	//methods
	this.setdefaultparameters = function(){// {{{

		//Type of stabilization to use 0:nothing 1:artificial_diffusivity
		this.stabilization=1;

	}// }}}
		this.disp= function(){// {{{

		console.log(sprintf('   hydrologyshreve solution parameters:'));
		fielddisplay(this,'spcwatercolumn','water thickness constraints (NaN means no constraint) [m]');
		fielddisplay(this,'stabilization','artificial diffusivity (default is 1). can be more than 1 to increase diffusivity.');

	}// }}}
    this.extrude = function(md) {//{{{
        return this;
    };//}}}
		this.classname= function(){// {{{
			return "hydrologyshreve";

		}// }}}
	this.type= function(){// {{{

		return "hydrologyshreve";
	}// }}}
		this.checkconsistency = function(md,solution,analyses) { //{{{

			//Early return
			if(!ArrayAnyEqual(ArrayIsMember('HydrologyShreveAnalysis',analyses),1)) return;

			checkfield(md,'fieldname','hydrology.spcwatercolumn','Inf',1,'timeseries',1);
			checkfield(md,'fieldname','hydrology.stabilization','>=',0);

		} // }}}
		this.marshall=function(md,prefix,fid) { //{{{
			WriteData(fid,prefix,'name','md.hydrology.model','data',2,'format','Integer');
			WriteData(fid,prefix,'object',this,'fieldname','spcwatercolumn','format','DoubleMat','mattype',1,'timeserieslength',md.mesh.numberofvertices+1,'yts',md.constants.yts);
			WriteData(fid,prefix,'object',this,'fieldname','stabilization','format','Double');
		}//}}}
		this.fix=function() { //{{{
			this.spcwatercolumn=NullFix(this.spcwatercolumn,NaN);
		}//}}}
	//properties 
	// {{{
	this.spcwatercolumn = NaN;
	this.stabilization  = 0;
	this.setdefaultparameters();
	//}}}
}
