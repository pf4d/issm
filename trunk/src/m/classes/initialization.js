//INITIALIZATION class definition
//
//   Usage:
//      initialization=new initialization();

function initialization (){
	//methods
	this.setdefaultparameters = function(){// {{{
	}// }}}
	this.disp= function(){// {{{
		console.log(sprintf('   initial field values:'));

		fielddisplay(this,'vx','x component of velocity [m/yr]');
		fielddisplay(this,'vy','y component of velocity [m/yr]');
		fielddisplay(this,'vz','z component of velocity [m/yr]');
		fielddisplay(this,'vel','velocity norm [m/yr]');
		fielddisplay(this,'pressure','pressure field [Pa]');
		fielddisplay(this,'temperature','temperature [K]');
		fielddisplay(this,'waterfraction','fraction of water in the ice');
		fielddisplay(this,'sediment_head','sediment water head of subglacial system [m]');
		fielddisplay(this,'epl_head','epl water head of subglacial system [m]');
		fielddisplay(this,'epl_thickness','epl layer thickness [m]');
		fielddisplay(this,'watercolumn','thickness of subglacial water [m]');

	}// }}}
    this.extrude = function(md) {//{{{
        this.vx=project3d(md,'vector',this.vx,'type','node');
        this.vy=project3d(md,'vector',this.vy,'type','node');
        this.vz=project3d(md,'vector',this.vz,'type','node');
        this.vel=project3d(md,'vector',this.vel,'type','node');
        this.temperature=project3d(md,'vector',this.temperature,'type','node');
        this.waterfraction=project3d(md,'vector',this.waterfraction,'type','node');
        this.watercolumn=project3d(md,'vector',this.watercolumn,'type','node','layer',1);
        this.sediment_head=project3d(md,'vector',this.sediment_head,'type','node','layer',1);
        this.epl_head=project3d(md,'vector',this.epl_head,'type','node','layer',1);
        this.epl_thickness=project3d(md,'vector',this.epl_thickness,'type','node','layer',1);

        //Lithostatic pressure by default
        this.pressure=md.constants.g*md.materials.rho_ice*(md.geometry.surface-md.mesh.z);
        return this;
    }//}}}
		this.checkconsistency = function(md,solution,analyses) { //{{{
			if(ArrayAnyEqual(ArrayIsMember('StressbalanceAnalysis',analyses),1)){
				if (!(isNaN(md.initialization.vx) | isNaN(md.initialization.vy))){
					checkfield(md,'fieldname','initialization.vx','NaN',1,'Inf',1,'size',[md.mesh.numberofvertices, 1]);
					checkfield(md,'fieldname','initialization.vy','NaN',1,'Inf',1,'size',[md.mesh.numberofvertices, 1]);
				}
			}
			if(ArrayAnyEqual(ArrayIsMember('MasstransportAnalysis',analyses),1)){
				checkfield(md,'fieldname','initialization.vx','NaN',1,'Inf',1,'size',[md.mesh.numberofvertices, 1]);
				checkfield(md,'fieldname','initialization.vy','NaN',1,'Inf',1,'size',[md.mesh.numberofvertices, 1]);
			}
			if(ArrayAnyEqual(ArrayIsMember('BalancethicknessSolution',analyses),1) & (solution=='BalancethicknessSolution')){
				checkfield(md,'fieldname','initialization.vx','NaN',1,'Inf',1,'size',[md.mesh.numberofvertices, 1]);
				checkfield(md,'fieldname','initialization.vy','NaN',1,'Inf',1,'size',[md.mesh.numberofvertices, 1]);
				//Triangle with zero velocity
				for(var i=0;i<md.mesh.numberofelements;i++){
					var sum=0;
					for(var j=0;j<md.mesh.elements[0].length;j++){
						if  ((md.initialization.vx[md.mesh.elements[i][j]-1]==0) & (md.initialization.vy[md.mesh.elements[i][j]-1]==0)) sum+=1;
					}
					if (sum==md.mesh.elements[0].length){
						md.checkmessage('at least one triangle has all its vertices with a zero velocity');
					}
				}
			}
			if(ArrayAnyEqual(ArrayIsMember('ThermalAnalysis',analyses),1)){
				checkfield(md,'fieldname','initialization.vx','NaN',1,'Inf',1,'size',[md.mesh.numberofvertices, 1]);
				checkfield(md,'fieldname','initialization.vy','NaN',1,'Inf',1,'size',[md.mesh.numberofvertices, 1]);
				if (md.mesh.dimension() == 3){
					checkfield(md,'fieldname','initialization.vz','NaN',1,'Inf',1,'size',[md.mesh.numberofvertices ,1]);
				}
				checkfield(md,'fieldname','initialization.pressure','NaN',1,'Inf',1,'size',[md.mesh.numberofvertices ,1]);
				checkfield(md,'fieldname','initialization.temperature','NaN',1,'Inf',1,'size',[md.mesh.numberofvertices ,1]);
			}
			if( ArrayAnyEqual(ArrayIsMember('EnthalpyAnalysis',analyses),1) & md.thermal.isenthalpy){
				checkfield(md,'fieldname','initialization.waterfraction','>=',0,'size',[md.mesh.numberofvertices, 1]);
				checkfield(md,'fieldname','initialization.watercolumn'  ,'>=',0,'size',[md.mesh.numberofvertices, 1]);
			}
			if(ArrayAnyEqual(ArrayIsMember('HydrologyShreveAnalysis',analyses),1)){
				if (md.hydrology.type() == 'hydrologyshreve'){
					checkfield(md,'fieldname','initialization.watercolumn','NaN',1,'Inf',1,'size',[md.mesh.numberofvertices ,1]);
				}
			}
			if(ArrayAnyEqual(ArrayIsMember('HydrologyDCInefficientAnalysis',analyses),1)){
				if (md.hydrology.type() == 'hydrologydc'){
					checkfield(md,'fieldname','initialization.sediment_head','NaN',1,'Inf',1,'size',[md.mesh.numberofvertices, 1]);
				}
			}
			if(ArrayAnyEqual(ArrayIsMember('HydrologyDCEfficientAnalysis',analyses),1)){
				if (md.hydrology.type() == 'hydrologydc'){
					if (md.hydrology.isefficientlayer==1){
						checkfield(md,'fieldname','initialization.epl_head','NaN',1,'Inf',1,'size',[md.mesh.numberofvertices ,1]);
						checkfield(md,'fieldname','initialization.epl_thickness','NaN',1,'Inf',1,'size',[md.mesh.numberofvertices ,1]);
					}
				}
			}
		} //}}}
		this.marshall=function(md,prefix,fid) { //{{{

			var yts=md.constants.yts;

			WriteData(fid,prefix,'object',this,'fieldname','vx','format','DoubleMat','mattype',1,'scale',1./yts);
			WriteData(fid,prefix,'object',this,'fieldname','vy','format','DoubleMat','mattype',1,'scale',1./yts);
			WriteData(fid,prefix,'object',this,'fieldname','vz','format','DoubleMat','mattype',1,'scale',1./yts);
			WriteData(fid,prefix,'object',this,'fieldname','pressure','format','DoubleMat','mattype',1);
			WriteData(fid,prefix,'object',this,'fieldname','temperature','format','DoubleMat','mattype',1);
			WriteData(fid,prefix,'object',this,'fieldname','waterfraction','format','DoubleMat','mattype',1);
			WriteData(fid,prefix,'object',this,'fieldname','sediment_head','format','DoubleMat','mattype',1);
			WriteData(fid,prefix,'object',this,'fieldname','epl_head','format','DoubleMat','mattype',1);
			WriteData(fid,prefix,'object',this,'fieldname','epl_thickness','format','DoubleMat','mattype',1);
			WriteData(fid,prefix,'object',this,'fieldname','watercolumn','format','DoubleMat','mattype',1);

			if (md.thermal.isenthalpy){
				tpmp=NewArrayFill(md.mesh.numberofvertices,0);
				for (var i=0;i<md.mesh.numberofvertices;i++) tpmp[i]= md.materials.meltingpoint - md.materials.beta*md.initialization.pressure[i];
				enthalpy=NewArrayFill(md.mesh.numberofvertices,0); 
				for (var i=0;i<md.mesh.numberofvertices;i++)enthalpy[i] = md.materials.heatcapacity*(md.initialization.temperature[i]-md.constants.referencetemperature);
				
				for (var i=0;i<md.mesh.numberofvertices;i++)if(md.initialization.temperature[i]>=tpmp[i]){
					enthalpy[i] = md.materials.heatcapacity*(tpmp[i] - md.constants.referencetemperature) + md.materials.latentheat*md.initialization.waterfraction[i];
				}
				WriteData(fid,prefix,'data',enthalpy,'format','DoubleMat','mattype',1,'name','md.initialization.enthalpy');
			}
		}//}}}
		this.fix=function(md) { //{{{
			this.vx=FloatFix(this.vx,md.mesh.numberofvertices); 
			this.vy=FloatFix(this.vx,md.mesh.numberofvertices); 
			this.vy=FloatFix(this.vx,md.mesh.numberofvertices); 
			this.waterfraction=NullFix(this.waterfraction,NaN);
			this.sediment_head=NullFix(this.sediment_head,NaN);
			this.epl_head=NullFix(this.epl_head,NaN);
			this.epl_thickness=NullFix(this.epl_thickness,NaN);
			this.watercolumn=NullFix(this.watercolumn,NaN);
		}//}}}
	//properties 
	// {{{
	this.vx            = NaN;
	this.vy            = NaN;
	this.vz            = NaN;
	this.vel           = NaN;
	this.pressure      = NaN;
	this.temperature   = NaN;
	this.waterfraction = NaN;
	this.sediment_head = NaN;
	this.epl_head      = NaN;
	this.epl_thickness = NaN;
	this.watercolumn   = NaN;
	this.setdefaultparameters();

	//}}}
}
