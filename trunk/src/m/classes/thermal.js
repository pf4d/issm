//THERMAL class definition
//
//   Usage:
//      thermal=new thermal();

function thermal (){
	//methods
	this.setdefaultparameters = function(){// {{{

		//Number of unstable constraints acceptable
		this.penalty_threshold=0;

		//Type of stabilization used
		this.stabilization=1;

		//Relative tolerance for the enthalpy convergence
		this.reltol=0.01;

		//Maximum number of iterations
		this.maxiter=100;

		//factor used to compute the values of the penalties: kappa=max(stiffness matrix)*10^penalty_factor
		this.penalty_factor=3;

		//Should we use cold ice (default) or enthalpy formulation
		this.isenthalpy=0;

		//will basal boundary conditions be set dynamically
		this.isdynamicbasalspc=0;

		//default output
		this.requested_outputs=['default'];

	}// }}}
	this.disp= function(){// {{{

		console.log(sprintf('   Thermal solution parameters:'));

		fielddisplay(this,'spctemperature','temperature constraints (NaN means no constraint) [K]');
		fielddisplay(this,'stabilization','0: no, 1: artificial_diffusivity, 2: SUPG');
		fielddisplay(this,'reltol','relative tolerance convergence criterion for enthalpy');
		fielddisplay(this,'maxiter','maximum number of non linear iterations');
		fielddisplay(this,'penalty_lock','stabilize unstable thermal constraints that keep zigzagging after n iteration (default is 0, no stabilization)');
		fielddisplay(this,'penalty_threshold','threshold to declare convergence of thermal solution (default is 0)');
		fielddisplay(this,'penalty_factor','scaling exponent (default is 3)');
		fielddisplay(this,'isenthalpy','use an enthalpy formulation to include temperate ice (default is 0)');
		fielddisplay(this,'isdynamicbasalspc','enable dynamic setting of basal forcing. required for enthalpy formulation (default is 0)');
		fielddisplay(this,'requested_outputs','additional outputs requested');

	}// }}}
	this.classname= function(){// {{{
		return "thermal";
	}// }}}
    this.extrude = function(md) {//{{{
        this.spctemperature=project3d(md,'vector',this.spctemperature,'type','node','layer',md.mesh.numberoflayers,'padding',NaN);
        if (md.initialization.temperature.length===md.mesh.numberofvertices) {
            this.spctemperature = ArrayFill2D(Create2DArray(md.mesh.numberofvertices, 1), NaN);
            var pos=ArrayFindNot2D(md.mesh.vertexonsurface, 0);
            pos.forEach(function(x) { // impose observed temperature on surface
                this.spctemperature[x] = md.initialization.temperature[x];
            });
        }

        return this;
    }//}}}
	this.checkconsistency = function(md,solution,analyses){ // {{{

		//Early return
		if(!ArrayAnyEqual(ArrayIsMember('ThermalAnalysis',analyses),1) & !ArrayAnyEqual(ArrayIsMember('EnthalpyAnalysis',analyses),1)  | (solution == 'TransientSolution' & md.trans.isthermal==0)) return;

		checkfield(md,'fieldname','thermal.stabilization','numel',[1],'values',[0 ,1, 2]);
		checkfield(md,'fieldname','thermal.spctemperature','Inf',1,'timeseries',1);
		if(ArrayAnyEqual(ArrayIsMember('EnthalpyAnalysis',analyses),1) & md.thermal.isenthalpy & md.mesh.dimension() == 3){
			
			for(var i=0;i<md.mesh.numberofvertices;i++){
				for(var j=0;j<md.thermal.spctemperature[0].length;j++){
					if (!isNaN(md.thermal.spctemperature[i][j])){
						var rep=md.geometry.surface[i]-md.mesh.z[i];
						if (md.thermal.spctemperature[i][j] <= md.materials.melting-md.materials.beta*md.materials.rho_ice*md.constants.g*rep+Math.pow(10,-5)){

							md.checkmessage('spctemperature should be less or equal than the adjusted melting point');
							break;
						}
					}
				}
			}
			checkfield(md,'fieldname','thermal.isenthalpy','numel',[1],'values',[0, 1]);
			checkfield(md,'fieldname','thermal.isdynamicbasalspc','numel', [1],'values',[0, 1]);
			if(md.thermal.isenthalpy){
				if (isNan(md.stressbalance.reltol)){
					md.checkmessage('for a steadystate computation, thermal.reltol (relative convergence criterion) must be defined!');
				}
			}
			checkfield(md,'fieldname','thermal.reltol','>',0.,'message','reltol must be larger than zero');
		}
		checkfield(md,'fieldname','thermal.requested_outputs','stringrow',1);
	} // }}} 
		this.marshall=function(md,prefix,fid) { //{{{
			WriteData(fid,prefix,'object',this,'fieldname','spctemperature','format','DoubleMat','mattype',1,'timeserieslength',md.mesh.numberofvertices+1,'yts',md.constants.yts);
			WriteData(fid,prefix,'object',this,'fieldname','penalty_threshold','format','Integer');
			WriteData(fid,prefix,'object',this,'fieldname','stabilization','format','Integer');
			WriteData(fid,prefix,'object',this,'fieldname','reltol','format','Double');
			WriteData(fid,prefix,'object',this,'fieldname','maxiter','format','Integer');
			WriteData(fid,prefix,'object',this,'fieldname','penalty_lock','format','Integer');
			WriteData(fid,prefix,'object',this,'fieldname','penalty_factor','format','Double');
			WriteData(fid,prefix,'object',this,'fieldname','isenthalpy','format','Boolean');
			WriteData(fid,prefix,'object',this,'fieldname','isdynamicbasalspc','format','Boolean');

			//process requested outputs
			var outputs = this.requested_outputs;
			for (var i=0;i<outputs.length;i++){
				if (outputs[i] == 'default') {
					outputs.splice(i,1);
					var newoutputs=this.defaultoutputs(md);
					for (var j=0;j<newoutputs.length;j++) outputs.push(newoutputs[j]);
				}
			}
			WriteData(fid,prefix,'data',outputs,'name','md.thermal.requested_outputs','format','StringArray');
        	}//}}}
		this.defaultoutputs = function(md) { //{{{

			if (this.isenthalpy) return ['Enthalpy','Temperature','Waterfraction','Watercolumn','BasalforcingsGroundediceMeltingRate'];
			else return ['Temperature','BasalforcingsGroundediceMeltingRate'];
		}//}}}
		this.fix=function() { //{{{
			this.spctemperature=NullFix(this.spctemperature,NaN);
		}//}}}
	//properties 
	// {{{

	this.spctemperature    = NaN;
	this.penalty_threshold = 0;
	this.stabilization     = 0;
	this.reltol				= 0;
	this.maxiter           = 0;
	this.penalty_lock      = 0;
	this.penalty_factor    = 0;
	this.isenthalpy        = 0;
	this.isdynamicbasalspc = 0;
	this.requested_outputs = [];

	this.setdefaultparameters();
	//}}}
}
