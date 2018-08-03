%SMBcomponents Class definition
%
%   Usage:
%      SMBcomponents=SMBcomponents();

classdef SMBcomponents
	properties (SetAccess=public) 
		accumulation = NaN;
		runoff = NaN;
		evaporation = NaN;
		requested_outputs      = {};
	end
	methods
		function self = SMBcomponents(varargin) % {{{
			switch nargin
				case 0
				otherwise
					error('constructor not supported');
			end
		end % }}}
		function self = extrude(self,md) % {{{

			self.accumulation=project3d(md,'vector',self.accumulation,'type','node');
			self.runoff=project3d(md,'vector',self.runoff,'type','node');
			self.evaporation=project3d(md,'vector',self.evaporation,'type','node');

		end % }}}
		function list = defaultoutputs(self,md) % {{{

			list = {''};

		end % }}}
		function self = initialize(self,md) % {{{

			if isnan(self.accumulation)
				self.accumulation=zeros(md.mesh.numberofvertices,1);
				disp('      no smb.accumulation specified: values set as zero');
			end
			if isnan(self.evaporation)
				self.evaporation=zeros(md.mesh.numberofvertices,1);
				disp('      no smb.evaporation specified: values set as zero');
			end
			if isnan(self.runoff)
				self.runoff=zeros(md.mesh.numberofvertices,1);
				disp('      no smb.runoff specified: values set as zero');
			end

		end % }}}
		function md = checkconsistency(self,md,solution,analyses) % {{{

			if ismember('MasstransportAnalysis',analyses),
				md = checkfield(md,'fieldname','smb.accumulation','timeseries',1,'NaN',1,'Inf',1);
			end
			if ismember('BalancethicknessAnalysis',analyses),
				md = checkfield(md,'fieldname','smb.accumulation','size',[md.mesh.numberofvertices 1],'NaN',1,'Inf',1);
			end
			if ismember('MasstransportAnalysis',analyses),
				md = checkfield(md,'fieldname','smb.runoff','timeseries',1,'NaN',1,'Inf',1);
			end
			if ismember('BalancethicknessAnalysis',analyses),
				md = checkfield(md,'fieldname','smb.runoff','size',[md.mesh.numberofvertices 1],'NaN',1,'Inf',1);
			end
			if ismember('MasstransportAnalysis',analyses),
				md = checkfield(md,'fieldname','smb.evaporation','timeseries',1,'NaN',1,'Inf',1);
			end
			if ismember('BalancethicknessAnalysis',analyses),
				md = checkfield(md,'fieldname','smb.evaporation','size',[md.mesh.numberofvertices 1],'NaN',1,'Inf',1);
			end
			md = checkfield(md,'fieldname','smb.requested_outputs','stringrow',1);
		end % }}}
		function disp(self) % {{{
			disp(sprintf('   surface forcings parameters (SMB=accumulation-runoff-evaporation) :'));
			fielddisplay(self,'accumulation','accumulated snow [m/yr ice eq]');
			fielddisplay(self,'runoff','amount of ice melt lost from the ice column [m/yr ice eq]');
			fielddisplay(self,'evaporation','amount of ice lost to evaporative processes [m/yr ice eq]');
			fielddisplay(self,'requested_outputs','additional outputs requested');
		end % }}}
		function marshall(self,prefix,md,fid) % {{{

			yts=md.constants.yts;

			WriteData(fid,prefix,'name','md.smb.model','data',2,'format','Integer');
			WriteData(fid,prefix,'object',self,'class','smb','fieldname','accumulation','format','DoubleMat','mattype',1,'scale',1./yts,'timeserieslength',md.mesh.numberofvertices+1,'yts',md.constants.yts);
			WriteData(fid,prefix,'object',self,'class','smb','fieldname','runoff','format','DoubleMat','mattype',1,'scale',1./yts,'timeserieslength',md.mesh.numberofvertices+1,'yts',md.constants.yts);
			WriteData(fid,prefix,'object',self,'class','smb','fieldname','evaporation','format','DoubleMat','mattype',1,'scale',1./yts,'timeserieslength',md.mesh.numberofvertices+1,'yts',md.constants.yts);
			
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
