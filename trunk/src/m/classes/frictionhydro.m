%FRICTIONWEERTMAN class definition
%
%   Usage:
%      friction=frictionhydro();

classdef frictionhydro
	properties (SetAccess=public) 
		Coupling           = 0;
		q                  = NaN;
		C                  = NaN;
		As                 = NaN;
		effective_pressure = NaN; 
	end
	methods
		function self = frictionhydro(varargin) % {{{
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
			md = checkfield(md,'fieldname','friction.Coupling','numel',[1],'values',[0 1]);
			md = checkfield(md,'fieldname','friction.q','NaN',1,'Inf',1,'size',[md.mesh.numberofelements 1]);
			md = checkfield(md,'fieldname','friction.C','NaN',1,'Inf',1,'size',[md.mesh.numberofelements 1]);
			md = checkfield(md,'fieldname','friction.As','NaN',1,'Inf',1,'size',[md.mesh.numberofelements 1]);
			if self.Coupling==0,
				md = checkfield(md,'fieldname','friction.effective_pressure','NaN',1,'Inf',1,'timeseries',1);
	    end
		end % }}}
		function self = extrude(self,md) % {{{
			self.q=project3d(md,'vector',self.q,'type','element');
			self.C=project3d(md,'vector',self.C,'type','element');
			self.As=project3d(md,'vector',self.As,'type','element');
			if self.Coupling==0,
				self.effective_pressure=project3d(md,'vector',self.effective_pressure,'type','node','layer',1);
			end
	  end % }}}
		function disp(self) % {{{
			disp(sprintf('Effective Pressure based friction law described in Gagliardini 2007'));
			fielddisplay(self,'Coupling','Coupling flag, 1 for coupling and 0 for forcing');
			fielddisplay(self,'q','friction law exponent q>=1');
			fielddisplay(self,'C','friction law max value [SI]');
			fielddisplay(self,'As','Sliding Parameter without cavitation [m Pa^-n s^-1]');
			fielddisplay(self,'effective_pressure','Effective Pressure for the forcing if not coupled [Pa]');
		end % }}}
		function marshall(self,prefix,md,fid) % {{{
			WriteData(fid,prefix,'name','md.friction.law','data',3,'format','Integer');
			WriteData(fid,prefix,'class','friction','object',self,'fieldname','Coupling','format','Integer');
			WriteData(fid,prefix,'class','friction','object',self,'fieldname','q','format','DoubleMat','mattype',2);
			WriteData(fid,prefix,'class','friction','object',self,'fieldname','C','format','DoubleMat','mattype',2);
			WriteData(fid,prefix,'class','friction','object',self,'fieldname','As','format','DoubleMat','mattype',2);
			if self.Coupling==0,
				WriteData(fid,prefix,'class','friction','object',self,'fieldname','effective_pressure','format','DoubleMat','mattype',1,'timeserieslength',md.mesh.numberofvertices+1,'yts',md.constants.yts);
			end
	  end % }}}
	end
end
