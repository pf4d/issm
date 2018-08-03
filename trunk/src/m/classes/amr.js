//AMR class definition
//
//   Usage:
//      amr=new amr();

function amr (){
	//methods
	this.setdefaultparameters = function(){// {{{
		//level_max: 2 to 4
		this.level_max=2;

		//region_level_1: region around (m) the discontinuity (grounding line or ice front) where the mesh will be refined once (h=1).
		this.region_level_1=20000.;

		//region_level_max: region around (m) the discontinuity (grounding line or ice front) where the mesh will be refined with max level of refinement (h=level_max).
		this.region_level_max=15000.;
	}// }}}
	this.disp= function(){// {{{

		console.log(sprintf('   amr parameters:'));
		fielddisplay(this,'level_max','maximum refinement level (1, 2, 3 or 4)');
		fielddisplay(this,'region_level_1','region which will be refined once (level 1) [ m ]');
		fielddisplay(this,'region_level_max','region which will be refined with level_max [ m ]');

	}// }}}
	this.classname= function(){// {{{
		return "amr";

	}// }}}
		this.checkconsistency = function(md,solution,analyses) { //{{{
			
			checkfield(md,'fieldname','amr.level_max','numel',[1],'>=',0,'<=',4);
			checkfield(md,'fieldname','amr.region_level_1','numel',[1],'>',0,'NaN',1,'Inf',1);
			checkfield(md,'fieldname','amr.region_level_max','numel',[1],'>',0,'NaN',1,'Inf',1);
			if (this.region_level_1-this.region_level_max<0.2*this.region_level_1){
				md.checkmessage('region_level_max should be lower than 80% of region_level_1');
			}
		} // }}}
		this.marshall=function(md,prefix,fid) { //{{{

			WriteData(fid,prefix,'object',this,'fieldname','level_max','format','Integer');
			WriteData(fid,prefix,'object',this,'fieldname','region_level_1','format','Double');
			WriteData(fid,prefix,'object',this,'fieldname','region_level_max','format','Double');

		}//}}}
		this.fix=function() { //{{{
		}//}}}
	//properties 
	// {{{
	this.level_max				= 0;
	this.region_level_1     = 0.;
	this.region_level_max   = 0.;

	this.setdefaultparameters();
	//}}}
}
