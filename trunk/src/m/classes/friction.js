//FRICTION class definition
//
//   Usage:
//      friction=friction();

function friction (){
	//methods
	this.setdefaultparameters = function(){ // {{{

	} // }}}
	this.disp= function (){// {{{
		console.log(sprintf('Basal shear stress parameters: Sigma_b = coefficient^2 * Neff ^r * |u_b|^(s-1) * u_b\n(effective stress Neff=rho_ice*g*thickness+rho_water*g*bed, r=q/p and s=1/p)'));
		fielddisplay(this,'coefficient','friction coefficient [SI]');
		fielddisplay(this,'p','p exponent');
		fielddisplay(this,'q','q exponent');
	} // }}}
    this.extrude = function(md) {//{{{
        this.coefficient = project3d(md, 'vector', this.coefficient, 'type', 'node', 'layer', 1);
        this.p = project3d(md, 'vector', this.p, 'type', 'element');
        this.q = project3d(md, 'vector', this.q, 'type', 'element');
        return this;
    }//}}}
	this.classname= function (){// {{{
		return "friction";
	} // }}}
		this.checkconsistency = function(md,solution,analyses){ //{{{

			//Early return
			if ((!ArrayAnyEqual(ArrayIsMember('StressbalanceAnalysis',analyses),1)) & (!ArrayAnyEqual(ArrayIsMember('StressbalanceAnalysis',analyses),1))){
				return; 
			}
			md = checkfield(md,'fieldname','friction.coefficient','timeseries',1,'NaN',1,'Inf',1);
			md = checkfield(md,'fieldname','friction.q','NaN',1,'Inf',1,'size',[md.mesh.numberofelements ,1]);
			md = checkfield(md,'fieldname','friction.p','NaN',1,'Inf',1,'size',[md.mesh.numberofelements ,1]);

		} // }}}
		this.marshall=function(md,prefix,fid) { //{{{
			var yts=md.constants.yts;

			WriteData(fid,prefix,'name','md.friction.law','data',1,'format','Integer');
			WriteData(fid,prefix,'object',this,'fieldname','coefficient','format','DoubleMat','mattype',1,'timeserieslength',md.mesh.numberofvertices+1,'yts',md.constants.yts);
			//WriteData(fid,prefix,'object',this,'fieldname','coefficient','format','DoubleMat','mattype',1);
			WriteData(fid,prefix,'object',this,'fieldname','p','format','DoubleMat','mattype',2);
			WriteData(fid,prefix,'object',this,'fieldname','q','format','DoubleMat','mattype',2);
			

		}//}}}
		this.fix=function() { //{{{
		}//}}}
	//properties 
	//{{{
	this.coefficient = NaN;
	this.p           = NaN;
	this.q           = NaN;
	this.setdefaultparameters();
	//}}}
}
