%GIA class definition for Ivins and James model 
%
%   Usage:
%      giaivins=giaivins();

classdef giaivins
	properties (SetAccess=public) 
		mantle_viscosity              = NaN;
		lithosphere_thickness         = NaN;
		cross_section_shape           = 0;
	end
	methods
		function self = extrude(self,md) % {{{
			self.mantle_viscosity=project3d(md,'vector',self.mantle_viscosity,'type','node');
			self.lithosphere_thickness=project3d(md,'vector',self.lithosphere_thickness,'type','node');
		end % }}}
		function self = giaivins(varargin) % {{{
			switch nargin
				case 0
					self=setdefaultparameters(self);
				otherwise
					error('constructor not supported');
			end
		end % }}}
		function self = setdefaultparameters(self) % {{{
		self.cross_section_shape=1; %square as default (see iedge in GiaDeflectionCorex)
		end % }}}
		function md = checkconsistency(self,md,solution,analyses) % {{{

			if ~ismember('GiaAnalysis',analyses), return; end
			md = checkfield(md,'fieldname','gia.mantle_viscosity','NaN',1,'Inf',1,'size',[md.mesh.numberofvertices 1],'>',0);
			md = checkfield(md,'fieldname','gia.lithosphere_thickness','NaN',1,'Inf',1,'size',[md.mesh.numberofvertices 1],'>',0);
			md = checkfield(md,'fieldname','gia.cross_section_shape','numel',[1],'values',[1,2]);

			%be sure that if we are running a masstransport ice flow model coupled with giaivins, that thickness forcings 
			%are not provided into the future.
			if strcmp(solution,'TransientSolution') & md.transient.ismasstransport & md.transient.isgia,
				%figure out if thickness is a transient forcing: 
				if size(md.geometry.thickness,1)==md.mesh.numberofvertices+1,
					%recover the furthest time "in time": 
					if(thickness(end,end)~=md.timestepping.start_time),
						md = checkmessage(md,['if ismasstransport is on, transient thickness forcing'...
							' for the giaivins model should not be provided in the future.'...
							' Synchronize your start_time to correspond to the most recent transient'...
							' thickness forcing timestep']);
					end
				end
			end

		end % }}}
		function disp(self) % {{{
			disp(sprintf('   giaivins parameters:'));

			fielddisplay(self,'mantle_viscosity','mantle viscosity[Pa s]');
			fielddisplay(self,'lithosphere_thickness','lithosphere thickness (km)');
			fielddisplay(self,'cross_section_shape','1: square-edged (default). 2: elliptical.  See iedge in GiaDeflectionCore');

		end % }}}
		function marshall(self,prefix,md,fid) % {{{
			WriteData(fid,prefix,'object',self,'fieldname','mantle_viscosity','format','DoubleMat','mattype',1);
			WriteData(fid,prefix,'object',self,'fieldname','lithosphere_thickness','format','DoubleMat','mattype',1,'scale',10^3); %from km to m
			WriteData(fid,prefix,'object',self,'fieldname','cross_section_shape','format','Integer');
		end % }}}
		function savemodeljs(self,fid,modelname) % {{{
		
			writejsdouble(fid,[modelname '.gia.mantle_viscosity'],self.mantle_viscosity);
			writejsdouble(fid,[modelname '.gia.lithosphere_thickness'],self.lithosphere_thickness);
			writejsdouble(fid,[modelname '.gia.cross_section_shape'],self.cross_section_shape);

		end % }}}
	end
end
