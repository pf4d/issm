%TIMESTEPPING Class definition
%
%   Usage:
%      timestepping=timestepping();

classdef timestepping
	properties (SetAccess=public) 
		start_time      = 0.;
		final_time      = 0.;
		time_step       = 0.;
		time_adapt      = 0;
		cfl_coefficient = 0.;
		interp_forcings = 1;
	end
	methods
		function self = timestepping(varargin) % {{{
			switch nargin
				case 0
					self=setdefaultparameters(self);
				otherwise
					error('constructor not supported');
			end
		end % }}}
		function self = setdefaultparameters(self) % {{{

			%time between 2 time steps
			self.time_step=1./2.;

			%final time
			self.final_time=10.*self.time_step;

			%time adaptation? 
			self.time_adapt=0;
			self.cfl_coefficient=0.5;

			%should we interpolate forcings between timesteps?
			self.interp_forcings=1;
		end % }}}
		function md = checkconsistency(self,md,solution,analyses) % {{{

			md = checkfield(md,'fieldname','timestepping.start_time','numel',[1],'NaN',1,'Inf',1);
			md = checkfield(md,'fieldname','timestepping.final_time','numel',[1],'NaN',1,'Inf',1);
			md = checkfield(md,'fieldname','timestepping.time_step','numel',[1],'>=',0,'NaN',1,'Inf',1);
			md = checkfield(md,'fieldname','timestepping.time_adapt','numel',[1],'values',[0 1]);
			md = checkfield(md,'fieldname','timestepping.cfl_coefficient','numel',[1],'>',0,'<=',1);
			md = checkfield(md,'fieldname','timestepping.interp_forcings','numel',[1],'values',[0 1]);
			if self.final_time-self.start_time<0,
				md = checkmessage(md,'timestepping.final_time should be larger than timestepping.start_time');
			end 
		end % }}}
		function disp(self) % {{{
			disp(sprintf('   timestepping parameters:'));

			unit = 'yr';
			fielddisplay(self,'start_time',['simulation starting time [' unit ']']);
			fielddisplay(self,'final_time',['final time to stop the simulation [' unit ']']);
			fielddisplay(self,'time_step',['length of time steps [' unit ']']);
			fielddisplay(self,'time_adapt','use cfl condition to define time step ? (0 or 1) ');
			fielddisplay(self,'cfl_coefficient','coefficient applied to cfl condition');
			fielddisplay(self,'interp_forcings','interpolate in time between requested forcing values ? (0 or 1)');

		end % }}}
		function marshall(self,prefix,md,fid) % {{{

			scale = md.constants.yts;
			WriteData(fid,prefix,'object',self,'fieldname','start_time','format','Double','scale',scale);
			WriteData(fid,prefix,'object',self,'fieldname','final_time','format','Double','scale',scale);
			WriteData(fid,prefix,'object',self,'fieldname','time_step','format','Double','scale',scale);
			WriteData(fid,prefix,'object',self,'fieldname','time_adapt','format','Boolean');
			WriteData(fid,prefix,'object',self,'fieldname','cfl_coefficient','format','Double');
			WriteData(fid,prefix,'object',self,'fieldname','interp_forcings','format','Boolean');
		end % }}}
		function savemodeljs(self,fid,modelname) % {{{
		
			writejsdouble(fid,[modelname '.timestepping.start_time'],self.start_time);
			writejsdouble(fid,[modelname '.timestepping.final_time'],self.final_time);
			writejsdouble(fid,[modelname '.timestepping.time_step'],self.time_step);
			writejsdouble(fid,[modelname '.timestepping.time_adapt'],self.time_adapt);
			writejsdouble(fid,[modelname '.timestepping.cfl_coefficient'],self.cfl_coefficient);
			writejsdouble(fid,[modelname '.timestepping.interp_forcings'],self.interp_forcings);

		end % }}}
	end
end
