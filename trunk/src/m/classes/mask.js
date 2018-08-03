//MASK class definition
//
//   Usage:
//      mask= new mask();

function mask () {
	//properties 
	// {{{
		this.groundedice_levelset                           = NaN;
		this.ice_levelset                           = NaN;
		//}}}
	//methods 
		this.setdefaultparameters = function (){ //{{{
		} // }}}
		this.disp = function () { //{{{
			console.log(sprintf("   mask:")); 

			fielddisplay(this,"groundedice_levelset","is ice grounded ? grounded ice if > 0, grounding line position if = 0, floating ice if < 0");
			fielddisplay(this,"ice_levelset","presence of ice if < 0, icefront position if = 0, no ice if > 0");
		} //}}}
		this.extrude = function(md) {//{{{
			this.groundedice_levelset=project3d(md,'vector',this.groundedice_levelset,'type','node');
			this.ice_levelset=project3d(md,'vector',this.ice_levelset,'type','node');
			return this;
		}//}}}
		this.classname = function () { //{{{
			return "mask";
		} //}}}
		this.checkconsistency = function(md,solution,analyses){ //{{{

			checkfield(md,'fieldname','mask.groundedice_levelset','size',[md.mesh.numberofvertices, 1]);
			checkfield(md,'fieldname','mask.ice_levelset'        ,'size',[md.mesh.numberofvertices, 1]);
			var isice=NewArrayFill(md.mesh.numberofvertices,0); 
			for(var i=0;i<md.mesh.numberofvertices;i++)if(md.mask.ice_levelset[i]<=0)isice[i]=1;
			if (ArraySum(isice)==0){
				console.log('no ice present in the domain');
			}
			if (ArrayMax(md.mask.ice_levelset)<0){
				console.log('no ice front provided');
			}
		} // }}}
		this.marshall=function(md,prefix,fid) { //{{{
            WriteData(fid,prefix,'object',this,'fieldname','groundedice_levelset','format','DoubleMat','mattype',1);
            WriteData(fid,prefix,'object',this,'fieldname','ice_levelset','format','DoubleMat','mattype',1);
		}//}}}
		this.fix=function() { //{{{
		}//}}}

}
