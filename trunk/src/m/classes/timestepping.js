//TIMESTEPPING class definition
//
//   Usage:
//      timestepping=new timestepping();

function timestepping (){
	//methods
	this.setdefaultparameters = function(){// {{{
		//time between 2 time steps
		this.time_step=1./2.;

		//final time
		this.final_time=10.*this.time_step;

		//time adaptation? 
		this.time_adapt=0;
		this.cfl_coefficient=0.5;

		//should we interpolate forcings between timesteps?
		this.interp_forcings=1;
	}// }}}
	this.disp= function(){// {{{

		var unit;
		console.log(sprintf('   timestepping parameters:'));
		unit = 'yr';
		fielddisplay(this,'start_time','simulation starting time ['+ unit + ']');
		fielddisplay(this,'final_time','final time to stop the simulation ['+ unit + ']');
		fielddisplay(this,'time_step','length of time steps [' +unit+ ']');
		fielddisplay(this,'time_adapt','use cfl condition to define time step ? (0 or 1) ');
		fielddisplay(this,'cfl_coefficient','coefficient applied to cfl condition');
		fielddisplay(this,'interp_forcings','interpolate in time between requested forcing values ? (0 or 1)');

	}// }}}
	this.classname= function(){// {{{
		return "timestepping";

	}// }}}
		this.checkconsistency = function(md,solution,analyses) { //{{{
			
			checkfield(md,'fieldname','timestepping.start_time','numel',[1],'NaN',1,'Inf',1);
			checkfield(md,'fieldname','timestepping.final_time','numel',[1],'NaN',1,'Inf',1);
			checkfield(md,'fieldname','timestepping.time_step','numel',[1],'>=',0,'NaN',1,'Inf',1);
			checkfield(md,'fieldname','timestepping.time_adapt','numel',[1],'values',[0,1]);
			checkfield(md,'fieldname','timestepping.cfl_coefficient','numel',[1],'>',0,'<=',1);
			checkfield(md,'fieldname','timestepping.interp_forcings','numel',[1],'values',[0,1]);
			if (this.final_time-this.start_time<0){
				md.checkmessage('timestepping.final_time should be larger than timestepping.start_time');
			}
		} // }}}
		this.marshall=function(md,prefix,fid) { //{{{

			var scale;
			scale = md.constants.yts;
			
			WriteData(fid,prefix,'object',this,'fieldname','start_time','format','Double','scale',scale);
			WriteData(fid,prefix,'object',this,'fieldname','final_time','format','Double','scale',scale);
			WriteData(fid,prefix,'object',this,'fieldname','time_step','format','Double','scale',scale);
			WriteData(fid,prefix,'object',this,'fieldname','time_adapt','format','Boolean');
			WriteData(fid,prefix,'object',this,'fieldname','cfl_coefficient','format','Double');
			WriteData(fid,prefix,'object',this,'fieldname','interp_forcings','format','Boolean');

		}//}}}
		this.fix=function() { //{{{
		}//}}}
	//properties 
	// {{{
	this.start_time      = 0.;
	this.final_time      = 0.;
	this.time_step       = 0.;
	this.time_adapt      = 0;
	this.cfl_coefficient = 0.;
	this.interp_forcings = 1;

	this.setdefaultparameters();
	//}}}
}
