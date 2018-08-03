%AMR Class definition
%
%   Usage:
%      amr=amr();

classdef amr
	properties (SetAccess=public) 
		level_max			= 0; 
		region_level_1		= 0;
		region_level_max	= 0;
	end
	methods
		function self = amr(varargin) % {{{
			switch nargin
				case 0
					self=setdefaultparameters(self);
				otherwise
					error('constructor not supported');
			end
		end % }}}
		function self = setdefaultparameters(self) % {{{

			%level_max: 2 to 4
			self.level_max=2;

			%region_level_1: region around (m) the discontinuity (grounding line or ice front) where the mesh will be refined once (h=1).
			self.region_level_1=20000.;
			
			%region_level_max: region around (m) the discontinuity (grounding line or ice front) where the mesh will be refined with max level of refinement (h=level_max). 
			self.region_level_max=15000.;

		end % }}}
		function md = checkconsistency(self,md,solution,analyses) % {{{

			md = checkfield(md,'fieldname','amr.level_max','numel',[1],'>=',0,'<=',4);
			md = checkfield(md,'fieldname','amr.region_level_1','numel',[1],'>',0,'NaN',1,'Inf',1);
			md = checkfield(md,'fieldname','amr.region_level_max','numel',[1],'>',0,'NaN',1,'Inf',1);
			if self.region_level_1-self.region_level_max<0.2*self.region_level_1, %it was adopted 20% of the region_level_1
				md = checkmessage(md,'region_level_max should be lower than 80% of region_level_1');
			end 
		end % }}}
		function disp(self) % {{{
			disp(sprintf('   amr parameters:'));

			fielddisplay(self,'level_max',['maximum refinement level (1, 2, 3 or 4)']);
			fielddisplay(self,'region_level_1',['region which will be refined once (level 1) [ m ]']);
			fielddisplay(self,'region_level_max',['region which will be refined with level_max [ m ]']);

		end % }}}
		function marshall(self,prefix,md,fid) % {{{

			scale = md.constants.yts;
			WriteData(fid,prefix,'object',self,'fieldname','level_max','format','Integer');
			WriteData(fid,prefix,'object',self,'fieldname','region_level_1','format','Double');
			WriteData(fid,prefix,'object',self,'fieldname','region_level_max','format','Double');
		end % }}}
		function savemodeljs(self,fid,modelname) % {{{
		
			writejsdouble(fid,[modelname '.amr.level_max'],self.level_max);
			writejsdouble(fid,[modelname '.amr.region_level_1'],self.region_level_1);
			writejsdouble(fid,[modelname '.amr.region_level_max'],self.region_level_max);

		end % }}}
	end
end
