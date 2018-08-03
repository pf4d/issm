%FRICTION class definition
%
%   Usage:
%      friction=friction();

classdef friction
	properties (SetAccess=public) 
		coefficient = NaN;
		p           = NaN;
		q           = NaN;
	end
	methods
		function self = extrude(self,md) % {{{
			self.coefficient=project3d(md,'vector',self.coefficient,'type','node','layer',1);
			self.p=project3d(md,'vector',self.p,'type','element');
			self.q=project3d(md,'vector',self.q,'type','element');
		end % }}}
		function self = friction(varargin) % {{{
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
			if (strcmp(solution,'TransientSolution') &  md.transient.isstressbalance ==0 & md.transient.isthermal == 0), return; end

			md = checkfield(md,'fieldname','friction.coefficient','timeseries',1,'NaN',1,'Inf',1);
			md = checkfield(md,'fieldname','friction.q','NaN',1,'Inf',1,'size',[md.mesh.numberofelements 1]);
			md = checkfield(md,'fieldname','friction.p','NaN',1,'Inf',1,'size',[md.mesh.numberofelements 1]);
		end % }}}
		function disp(self) % {{{
			disp(sprintf('Basal shear stress parameters: Sigma_b = coefficient^2 * Neff ^r * |u_b|^(s-1) * u_b\n(effective stress Neff=rho_ice*g*thickness+rho_water*g*bed, r=q/p and s=1/p)'));
			fielddisplay(self,'coefficient','friction coefficient [SI]');
			fielddisplay(self,'p','p exponent');
			fielddisplay(self,'q','q exponent');
		end % }}}
		function marshall(self,prefix,md,fid) % {{{
			yts=md.constants.yts;

			WriteData(fid,prefix,'name','md.friction.law','data',1,'format','Integer');
			WriteData(fid,prefix,'object',self,'fieldname','coefficient','format','DoubleMat','mattype',1,'timeserieslength',md.mesh.numberofvertices+1,'yts',md.constants.yts);
			WriteData(fid,prefix,'object',self,'fieldname','p','format','DoubleMat','mattype',2);
			WriteData(fid,prefix,'object',self,'fieldname','q','format','DoubleMat','mattype',2);
			

		end % }}}
		function savemodeljs(self,fid,modelname) % {{{
		
			writejs1Darray(fid,[modelname '.friction.coefficient'],self.coefficient);
			writejs1Darray(fid,[modelname '.friction.p'],self.p);
			writejs1Darray(fid,[modelname '.friction.q'],self.q);

		end % }}}
	end
end
