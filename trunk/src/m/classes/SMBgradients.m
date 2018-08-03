%SMBgradients Class definition
%
%   Usage:
%      SMBgradients=SMBgradients();

classdef SMBgradients
	properties (SetAccess=public) 
		href   = NaN;
		smbref = NaN;
		b_pos  = NaN;
		b_neg  = NaN;
		requested_outputs      = {};
	end
	methods
		function self = SMBgradients(varargin) % {{{
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

			%Nothing for now

		end % }}}
		function md = checkconsistency(self,md,solution,analyses) % {{{

			if ismember('MasstransportAnalysis',analyses),
				md = checkfield(md,'fieldname','smb.href','timeseries',1,'NaN',1,'Inf',1);
				md = checkfield(md,'fieldname','smb.smbref','timeseries',1,'NaN',1,'Inf',1);
				md = checkfield(md,'fieldname','smb.b_pos','timeseries',1,'NaN',1,'Inf',1);
				md = checkfield(md,'fieldname','smb.b_neg','timeseries',1,'NaN',1,'Inf',1);
			end
			md = checkfield(md,'fieldname','smb.requested_outputs','stringrow',1);
		end % }}}
		function disp(self) % {{{
			disp(sprintf('   surface forcings parameters:'));

			disp(sprintf('\n   SMB gradients parameters:'));
			fielddisplay(self,'href',' reference elevation from which deviation is used to calculate SMB adjustment in smb gradients method [m]');
			fielddisplay(self,'smbref',' reference smb from which deviation is calculated in smb gradients method [mm/yr water equiv]');
			fielddisplay(self,'b_pos',' slope of hs - smb regression line for accumulation regime required if smb gradients is activated');
			fielddisplay(self,'b_neg',' slope of hs - smb regression line for ablation regime required if smb gradients is activated');
			fielddisplay(self,'requested_outputs','additional outputs requested');

		end % }}}
		function marshall(self,prefix,md,fid) % {{{

			yts=md.constants.yts;

			WriteData(fid,prefix,'name','md.smb.model','data',6,'format','Integer');
			WriteData(fid,prefix,'object',self,'class','smb','fieldname','href','format','DoubleMat','mattype',1,'timeserieslength',md.mesh.numberofvertices+1,'yts',md.constants.yts);
			WriteData(fid,prefix,'object',self,'class','smb','fieldname','smbref','format','DoubleMat','mattype',1,'scale',1./yts,'timeserieslength',md.mesh.numberofvertices+1,'yts',md.constants.yts);
			WriteData(fid,prefix,'object',self,'class','smb','fieldname','b_pos','format','DoubleMat','mattype',1,'scale',1./yts,'timeserieslength',md.mesh.numberofvertices+1,'yts',md.constants.yts);
			WriteData(fid,prefix,'object',self,'class','smb','fieldname','b_neg','format','DoubleMat','mattype',1,'scale',1./yts,'timeserieslength',md.mesh.numberofvertices+1,'yts',md.constants.yts);
			
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
