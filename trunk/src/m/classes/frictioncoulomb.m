%FRICTIONCOULOMB class definition
%
%   Usage:
%      frictioncoulomb=frictioncoulomb();

classdef frictioncoulomb
	properties (SetAccess=public) 
		coefficient        = NaN;
		p                  = NaN;
		q                  = NaN;
		coefficientcoulomb = NaN;
	end
	methods
		function self = extrude(self,md) % {{{
			self.coefficient=project3d(md,'vector',self.coefficient,'type','node','layer',1);
			self.coefficientcoulomb=project3d(md,'vector',self.coefficientcoulomb,'type','node','layer',1);
			self.p=project3d(md,'vector',self.p,'type','element');
			self.q=project3d(md,'vector',self.q,'type','element');
		end % }}}
		function self = frictioncoulomb(varargin) % {{{
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
			md = checkfield(md,'fieldname','friction.coefficient','timeseries',1,'NaN',1,'Inf',1);
			md = checkfield(md,'fieldname','friction.coefficientcoulomb','timeseries',1,'NaN',1,'Inf',1);
			md = checkfield(md,'fieldname','friction.q','NaN',1,'Inf',1,'size',[md.mesh.numberofelements 1]);
			md = checkfield(md,'fieldname','friction.p','NaN',1,'Inf',1,'size',[md.mesh.numberofelements 1]);
		end % }}}
		function disp(self) % {{{
			disp(sprintf('Basal shear stress parameters: Sigma_b = min( coefficient^2 * Neff ^r * |u_b|^(s-1) * u_b\n, coefficientcoulomb^2 * rho_i * g * (h-h_f)) (effective stress Neff=rho_ice*g*thickness+rho_water*g*bed, r=q/p and s=1/p, floatation thickness h_f=max(0,-rho_sw / rho_i * bed))'));
			fielddisplay(self,'coefficient','power law (Weertman) friction coefficient [SI]');
			fielddisplay(self,'coefficientcoulomb','Coulomb friction coefficient [SI]');
			fielddisplay(self,'p','p exponent');
			fielddisplay(self,'q','q exponent');
		end % }}}
		function marshall(self,prefix,md,fid) % {{{

			WriteData(fid,prefix,'name','md.friction.law','data',7,'format','Integer');
			WriteData(fid,prefix,'object',self,'fieldname','coefficient','format','DoubleMat','mattype',1,'timeserieslength',md.mesh.numberofvertices+1,'yts',md.constants.yts);
			WriteData(fid,prefix,'object',self,'fieldname','coefficientcoulomb','format','DoubleMat','mattype',1);
			WriteData(fid,prefix,'object',self,'fieldname','p','format','DoubleMat','mattype',2);
			WriteData(fid,prefix,'object',self,'fieldname','q','format','DoubleMat','mattype',2);

		end % }}}
	end
end
