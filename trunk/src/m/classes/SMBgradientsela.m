%SMBgradientsela Class definition
%
%   Usage:
%      SMBgradientsela=SMBgradientsela();

classdef SMBgradientsela
	properties (SetAccess=public) 
		ela    = NaN;
		b_pos  = NaN;
		b_neg  = NaN;
		b_max  = NaN;
		b_min  = NaN;
		requested_outputs      = {};
	end
	methods
		function self = SMBgradientsela(varargin) % {{{
			switch nargin
				case 0
					self=setdefaultparameters(self);
				otherwise
					error('constructor not supported');
			end
		end % }}}
		function self = extrude(self,md) % {{{

			%Nothing for now

		end % }}}
		function list = defaultoutputs(self,md) % {{{
			list = {''};
		end % }}}
		function self = initialize(self,md) % {{{

			%Nothing done for now

		end % }}}
		function self = setdefaultparameters(self) % {{{

			self.b_max=9999;
			self.b_min=-9999;

		end % }}}
		function md = checkconsistency(self,md,solution,analyses) % {{{

			if ismember('MasstransportAnalysis',analyses),
				md = checkfield(md,'fieldname','smb.ela','timeseries',1,'NaN',1,'Inf',1);
				md = checkfield(md,'fieldname','smb.b_pos','timeseries',1,'NaN',1,'Inf',1);
				md = checkfield(md,'fieldname','smb.b_neg','timeseries',1,'NaN',1,'Inf',1);
				md = checkfield(md,'fieldname','smb.b_max','timeseries',1,'NaN',1,'Inf',1);
				md = checkfield(md,'fieldname','smb.b_min','timeseries',1,'NaN',1,'Inf',1);
			end
			md = checkfield(md,'fieldname','smb.requested_outputs','stringrow',1);
		end % }}}
		function disp(self) % {{{
			disp(sprintf('   surface forcings parameters:'));

			disp(sprintf('\n   SMB gradients ela parameters:'));
			fielddisplay(self,'ela',' equilibrium line altitude from which deviation is used to calculate smb using the smb gradients ela method [m a.s.l.]');
			fielddisplay(self,'b_pos',' vertical smb gradient (dB/dz) above ela');
			fielddisplay(self,'b_neg',' vertical smb gradient (dB/dz) below ela');
			fielddisplay(self,'b_max',' upper cap on smb rate, default: 9999 (no cap) [m ice eq./yr] ');
			fielddisplay(self,'b_min',' lower cap on smb rate, default: -9999 (no cap) [m ice eq./yr]');
			fielddisplay(self,'requested_outputs','additional outputs requested');

		end % }}}
		function marshall(self,prefix,md,fid) % {{{

			yts=md.constants.yts;

			WriteData(fid,prefix,'name','md.smb.model','data',9,'format','Integer');
			WriteData(fid,prefix,'object',self,'class','smb','fieldname','ela','format','DoubleMat','mattype',1,'timeserieslength',md.mesh.numberofvertices+1,'yts',md.constants.yts);
			WriteData(fid,prefix,'object',self,'class','smb','fieldname','b_pos','format','DoubleMat','mattype',1,'scale',1./yts,'timeserieslength',md.mesh.numberofvertices+1,'yts',md.constants.yts);
			WriteData(fid,prefix,'object',self,'class','smb','fieldname','b_neg','format','DoubleMat','mattype',1,'scale',1./yts,'timeserieslength',md.mesh.numberofvertices+1,'yts',md.constants.yts);
			WriteData(fid,prefix,'object',self,'class','smb','fieldname','b_max','format','DoubleMat','mattype',1,'scale',1./yts,'timeserieslength',md.mesh.numberofvertices+1,'yts',md.constants.yts);
			WriteData(fid,prefix,'object',self,'class','smb','fieldname','b_min','format','DoubleMat','mattype',1,'scale',1./yts,'timeserieslength',md.mesh.numberofvertices+1,'yts',md.constants.yts);
			
			%process requested outputs
			outputs = self.requested_outputs;
			pos  = find(ismember(outputs,'default'));
			if ~isempty(pos),
				outputs(pos) = [];                         %remove 'default' from outputs
				outputs      = [outputs defaultoutputs(self,md)]; %add defaults
			end
			WriteData(fid,prefix,'data',outputs,'name','md.smb.requested_outputs','format','StringArray');

		end % }}}
	end
end
