%MASK class definition
%
%   Usage:
%      mask=mask();

classdef mask
	properties (SetAccess=public) 
		groundedice_levelset = NaN;
		ice_levelset         = NaN;
	end
	methods (Static)
		function self = loadobj(self) % {{{
			% This function is directly called by matlab when a model object is
			% loaded. Update old properties here

			%2014 February 5th
			if numel(self.ice_levelset)>1 & all(self.ice_levelset>=0),
				disp('WARNING: md.mask.ice_levelset>=0, you probably need to change the sign of this levelset');
			end

		end% }}}
	end
	methods
		function self = extrude(self,md) % {{{
			self.groundedice_levelset=project3d(md,'vector',self.groundedice_levelset,'type','node');
			self.ice_levelset=project3d(md,'vector',self.ice_levelset,'type','node');
		end % }}}
		function self = mask(varargin) % {{{
			switch nargin
				case 0
					self=setdefaultparameters(self);
				otherwise
					error('constructor not supported');
			end
		end % }}}
		function self = setdefaultparameters(self) % {{{

		end % }}}
		function md = checkconsistency(self,md,solution,analyses) % {{{

			md = checkfield(md,'fieldname','mask.groundedice_levelset','forcing',1,'NaN',1);
			md = checkfield(md,'fieldname','mask.ice_levelset','NaN',1,'size',[md.mesh.numberofvertices 1]);
			isice=(md.mask.ice_levelset<=0);
			if sum(isice)==0,
				warning('no ice present in the domain');
			end
			if max(md.mask.ice_levelset)<0,
				disp('WARNING: no ice front provided');
			end
		end % }}}
		function disp(self) % {{{
			disp(sprintf('   masks:'));

			fielddisplay(self,'groundedice_levelset','is ice grounded ? grounded ice if > 0, grounding line position if = 0, floating ice if < 0');
			fielddisplay(self,'ice_levelset','presence of ice if < 0, icefront position if = 0, no ice if > 0');
		end % }}}
		function marshall(self,prefix,md,fid) % {{{
			WriteData(fid,prefix,'object',self,'fieldname','groundedice_levelset','format','DoubleMat','mattype',1,'timeserieslength',md.mesh.numberofvertices+1,'yts',md.constants.yts);
			WriteData(fid,prefix,'object',self,'fieldname','ice_levelset','format','DoubleMat','mattype',1);
		end % }}}
		function savemodeljs(self,fid,modelname) % {{{
		
			writejs1Darray(fid,[modelname '.mask.groundedice_levelset'],self.groundedice_levelset);
			writejs1Darray(fid,[modelname '.mask.ice_levelset'],self.ice_levelset);

		end % }}}
	end
end
