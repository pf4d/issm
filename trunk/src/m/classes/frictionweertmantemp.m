%FRICTIONWEERTMAN class definition
%
%   Usage:
%      frictionweertmantemp=frictionweertmantemp();

classdef frictionweertmantemp
	properties (SetAccess=public) 
		gamma = 0;
		C = NaN;
		m = NaN;
	end
	methods
		function self = frictionweertmantemp(varargin) % {{{
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

			%Early return
			if ~ismember('StressbalanceAnalysis',analyses) & ~ismember('ThermalAnalysis',analyses), return; end
			md = checkfield(md,'fieldname','friction.C','timeseries',1,'NaN',1,'Inf',1);
			md = checkfield(md,'fieldname','friction.m','NaN',1,'Inf',1,'size',[md.mesh.numberofelements 1]);
		end % }}}
		function disp(self) % {{{
			disp('Weertman sliding law parameters:');
			disp('      Sigma_b = C^(-1/m) * |u_b|^(1/m-1)  u_b * 1/f(T)');
			disp(' ');
			fielddisplay(self,'gamma','submelt sliding parameter f(T) = exp((T-Tpmp)/gamma)');
			fielddisplay(self,'C','friction coefficient [SI]');
			fielddisplay(self,'m','m exponent');
		end % }}}
		function marshall(self,prefix,md,fid) % {{{
			yts=md.constants.yts;

			WriteData(fid,prefix,'name','md.friction.law','data',6,'format','Integer');
			WriteData(fid,prefix,'class','friction','object',self,'fieldname','gamma','format','Double');
			WriteData(fid,prefix,'class','friction','object',self,'fieldname','C','format','DoubleMat','mattype',1,'timeserieslength',md.mesh.numberofvertices+1,'yts',md.constants.yts);
			WriteData(fid,prefix,'class','friction','object',self,'fieldname','m','format','DoubleMat','mattype',2);
			

		end % }}}
	end
end
