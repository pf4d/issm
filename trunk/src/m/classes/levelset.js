//LEVELSET class definition
//
//   Usage:
//      levelset=new levelset();

function levelset (){
	//methods
	this.setdefaultparameters = function(){// {{{

		//stabilization = 2 by default
		this.stabilization		= 2;
		this.reinit_frequency	= 5;
	
	}// }}}
	this.disp= function(){// {{{

		console.log(sprintf('   Level-set parameters:'));
		fielddisplay(this,'stabilization','0: no, 1: artificial_diffusivity, 2: streamline upwinding');
		fielddisplay(this,'spclevelset','Levelset constraints (NaN means no constraint)');
		fielddisplay(this,'reinit_frequency','Amount of time steps after which the levelset function in re-initialized (NaN: no re-initialization).');

	}// }}}
    this.extrude = function(md) {//{{{
        this.spclevelset=project3d(md,'vector',this.spclevelset,'type','node');
        return this;
    }//}}}
	this.classname= function(){// {{{
		return "levelset";
	}// }}}
	this.checkconsistency = function(md,solution,analyses) { // {{{
		//Early return
		if (solution!='TransientSolution' | md.trans.ismovingfront==0) return;

		checkfield(md,'fieldname','levelset.spclevelset','Inf',1,'timeseries',1);
		checkfield(md,'fieldname','levelset.stabilization','values',[0,1,2]);
	} //}}}
	this.marshall=function(md,prefix,fid) { //{{{
		WriteData(fid,prefix,'object',this,'fieldname','stabilization','format','Integer');
		WriteData(fid,prefix,'object',this,'fieldname','spclevelset','format','DoubleMat','mattype',1,'timeserieslength',md.mesh.numberofvertices+1,'yts',md.constants.yts);
		WriteData(fid,prefix,'object',this,'fieldname','reinit_frequency','format','Integer');

	}//}}}
		this.fix=function() { //{{{
			this.spclevelset=NullFix(this.spclevelset,NaN);
		}//}}}
	//properties 
	// {{{

	this.stabilization		= 0;
	this.spclevelset			= NaN;
	this.reinit_frequency	= NaN;

	this.setdefaultparameters();
	//}}}
}
