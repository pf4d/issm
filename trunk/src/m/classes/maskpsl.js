//MASKPSL class definition
//
//   Usage:
//      maskpsl= new maskpsl();

function maskpsl () {
	//properties 
	// {{{
		this.groundedice_levelset                           = NaN;
		this.ice_levelset                           = NaN;
		this.land_levelset                           = NaN;
		this.ocean_levelset                           = NaN;
		//}}}
	//methods 
		this.setdefaultparameters = function (){ //{{{
		} // }}}
		this.disp = function () { //{{{
			console.log(sprintf("   mask:")); 

			fielddisplay(this,"groundedice_levelset","is ice grounded ? grounded ice if > 0, grounding line position if = 0, floating ice if < 0");
			fielddisplay(this,"ice_levelset","presence of ice if < 0, icefront position if = 0, no ice if > 0");
			fielddisplay(this,"ocean_levelset","is the vertex on the ocean? yes if = 1, no if = 0");
			fielddisplay(this,"land_levelset","is the vertex on land? yes if = 1, no if = 0");
		} //}}}
		this.classname = function () { //{{{
			return "maskpsl";
		} //}}}
		this.checkconsistency = function(md,solution,analyses){ //{{{

			checkfield(md,'fieldname','mask.groundedice_levelset','size',[md.mesh.numberofvertices, 1]);
			checkfield(md,'fieldname','mask.ice_levelset'        ,'size',[md.mesh.numberofvertices, 1]);
			checkfield(md,'fieldname','mask.ocean_levelset'        ,'size',[md.mesh.numberofvertices, 1]);
			checkfield(md,'fieldname','mask.land_levelset'        ,'size',[md.mesh.numberofvertices, 1]);
			
			var isice=NewArrayFill(md.mesh.numberofvertices,0); 
			for(var i=0;i<md.mesh.numberofvertices;i++)if(md.mask.ice_levelset[i]<=0)isice[i]=1;
			if (ArraySum(isice)==0){
				console.log('no ice present in the domain');
			}
			if (ArrayMax(md.mask.ice_levelset)<0){
				console.log('no ice front provided');
			}
				
			var icefront=NewArrayFill(md.mesh.numberofelements,0);
			for(var i=0;i<md.mesh.numberofelements;i++){
				for(var j=0;j<md.mesh.elements[0].length;j++){
					icefront[i]+=(md.mask.ice_levelset[md.mesh.elements[i][j]-1]==0);
				}
			}
			if ((ArrayMax(icefront)==3 & (md.mesh.elementtype() == 'Tria')) | (ArrayMax(icefront)==6 & md.mesh.elementtype() == 'Penta')){
				if (md.mesh.elementtype()=='Tria'){
					var pos=ArrayFindEqual(icefront,3); numberemptyelements=pos.length;
				}
				else if (md.mesh.elementtype() == 'Penta'){
					var pos=ArrayFindEqual(icefront,6); numberemptyelements=pos.length;
				}
				throw Error(sprintf(" %i have all nodes on ice front, change md.mask.ice_levelset to fix it",numberemptyelements));
			}
		} // }}}
		this.marshall=function(md,prefix,fid) { //{{{
			WriteData(fid,prefix,'object',this,'class','mask','fieldname','groundedice_levelset','format','DoubleMat','mattype',1);
			WriteData(fid,prefix,'object',this,'class','mask','fieldname','ice_levelset','format','DoubleMat','mattype',1);
			WriteData(fid,prefix,'object',this,'class','mask','fieldname','ocean_levelset','format','DoubleMat','mattype',1);
			WriteData(fid,prefix,'object',this,'class','mask','fieldname','land_levelset','format','DoubleMat','mattype',1);
		}//}}}
		this.fix=function() { //{{{
		}//}}}

}
