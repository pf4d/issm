%CALVINGDEV class definition
%
%   Usage:
%      calvingdev=calvingdev();

classdef calvingdev
	properties (SetAccess=public) 
		stress_threshold_groundedice = 0.;
		stress_threshold_floatingice = 0.;
		meltingrate   = NaN;
	end
	methods
		function self = calvingdev(varargin) % {{{
			switch nargin
				case 0
					self=setdefaultparameters(self);
				case 1
					inputstruct=varargin{1};
					list1 = properties('calvingdev');
					list2 = fieldnames(inputstruct);
					for i=1:length(list1)
						fieldname = list1{i};
						if ismember(fieldname,list2),
							self.(fieldname) = inputstruct.(fieldname);
						end
					end
				otherwise
					error('constructor not supported');
			end
		end % }}}
		function self = extrude(self,md) % {{{
			self.meltingrate=project3d(md,'vector',self.meltingrate,'type','node');
		end % }}}
		function self = setdefaultparameters(self) % {{{

			%Default sigma max
			self.stress_threshold_groundedice = 1e6;
			self.stress_threshold_floatingice = 150e3;
		end % }}}
		function md = checkconsistency(self,md,solution,analyses) % {{{
			%Early return
			if (~strcmp(solution,'TransientSolution') | md.transient.ismovingfront==0), return; end

			md = checkfield(md,'fieldname','calving.stress_threshold_groundedice','>',0,'numel',1,'nan',1,'Inf',1);
			md = checkfield(md,'fieldname','calving.stress_threshold_floatingice','>',0,'numel',1,'nan',1,'Inf',1);
			md = checkfield(md,'fieldname','calving.meltingrate','NaN',1,'Inf',1,'timeseries',1,'>=',0);
		end % }}}
		function disp(self) % {{{
			disp(sprintf('   Calving Pi parameters:'));
			fielddisplay(self,'stress_threshold_groundedice','sigma_max applied to grounded ice only [Pa]');
			fielddisplay(self,'stress_threshold_floatingice','sigma_max applied to floating ice only [Pa]');
			fielddisplay(self,'meltingrate','melting rate at given location [m/a]');

		end % }}}
		function marshall(self,prefix,md,fid) % {{{
			yts=md.constants.yts;
			WriteData(fid,prefix,'name','md.calving.law','data',2,'format','Integer');
			WriteData(fid,prefix,'object',self,'fieldname','stress_threshold_groundedice','format','Double');
			WriteData(fid,prefix,'object',self,'fieldname','stress_threshold_floatingice','format','Double');
			WriteData(fid,prefix,'object',self,'fieldname','meltingrate','format','DoubleMat','mattype',1,'timeserieslength',md.mesh.numberofvertices+1,'yts',md.constants.yts,'scale',1./yts);
		end % }}}
	end
end
