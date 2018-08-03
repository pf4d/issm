//GROUNDINGLINE class definition
//
//   Usage:
//      groundingline=new groundingline();

function groundingline (){
	//methods
	this.setdefaultparameters = function(){// {{{
		//Type of migration
		this.migration='None';

	}// }}}
	this.disp= function(){// {{{
		console.log(sprintf('   grounding line migration parameters:'));
		fielddisplay(this,'migration',"type of grounding line migration: 'SoftMigration','AggressiveMigration','SubelementMigration','SubelementMigration2' or 'None'");


	}// }}}
	this.classname= function(){// {{{
		return "groundingline";
	}// }}}
		this.checkconsistency = function(md,solution,analyses) {// {{{

			checkfield(md,'fieldname','groundingline.migration','values',['None', 'AggressiveMigration', 'SoftMigration', 'SubelementMigration', 'SubelementMigration2', 'Contact', 'None', 'GroundingOnly']);

			if (this.migration !='None'){
				if (isNaN(md.geometry.bed)){
					md.checkmessage('requesting grounding line migration, but bathymetry is absent!');
				}
				for (var i=0;i<md.mesh.numberofvertices;i++){
					if(md.mask.groundedice_levelset[i]>0){
						md.checkmessage('base not equal to bed on grounded ice!');
						break;
					}
					if(md.geometry.bed[i] - md.geometry.base[i] > Math.pow(10,-9)){
						md = checkmessage(md,'bed superior to base on floating ice!');
						break;
					}
				}
			}
		} // }}}
		this.marshall=function(md,prefix,fid) { //{{{
			WriteData(fid,prefix,'data',this.migration,'name','md.groundingline.migration','format','String');
		}//}}}
		this.fix=function() { //{{{
		}//}}}
	//properties 
	// {{{
	this.migration    = '';
	this.setdefaultparameters();
	//}}}
}
