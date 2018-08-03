%HYDROLOGYSHREVE class definition
%
%   Usage:
%      hydrologyshreve=hydrologyshreve();

classdef hydrologyshreve
	properties (SetAccess=public) 
		spcwatercolumn = NaN;
		stabilization  = 0;
	end
	methods
		function self = extrude(self,md) % {{{
		end % }}}
		function self = hydrologyshreve(varargin) % {{{
			switch nargin
				case 0
					self=setdefaultparameters(self);
				case 1
					self=structtoobj(self,varargin{1});
				otherwise
					error('constructor not supported');
			end
		end % }}}
		function self = setdefaultparameters(self) % {{{

			%Type of stabilization to use 0:nothing 1:artificial_diffusivity
			self.stabilization=1;
		end % }}}
		function md = checkconsistency(self,md,solution,analyses) % {{{

			%Early return
			if ~ismember('HydrologyShreveAnalysis',analyses)
				return;
			end

			md = checkfield(md,'fieldname','hydrology.spcwatercolumn','Inf',1,'timeseries',1);
			md = checkfield(md,'fieldname','hydrology.stabilization','>=',0);
		end % }}}
		function disp(self) % {{{
			disp(sprintf('   hydrologyshreve solution parameters:'));
			fielddisplay(self,'spcwatercolumn','water thickness constraints (NaN means no constraint) [m]');
			fielddisplay(self,'stabilization','artificial diffusivity (default is 1). can be more than 1 to increase diffusivity.');

		end % }}}
		function marshall(self,prefix,md,fid) % {{{
			WriteData(fid,prefix,'name','md.hydrology.model','data',2,'format','Integer');
			WriteData(fid,prefix,'object',self,'fieldname','spcwatercolumn','format','DoubleMat','mattype',1,'timeserieslength',md.mesh.numberofvertices+1,'yts',md.constants.yts);
			WriteData(fid,prefix,'object',self,'fieldname','stabilization','format','Double');
		end % }}}
		function savemodeljs(self,fid,modelname) % {{{
		
			writejs1Darray(fid,[modelname '.hydrology.spcwatercolumn'],self.spcwatercolumn);
			writejsdouble(fid,[modelname '.hydrology.stabilization'],self.stabilization);

		end % }}}
	end
end

