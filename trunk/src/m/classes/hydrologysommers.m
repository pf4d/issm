%HYDROLOGYSOMMERS class definition
%
%   Usage:
%      hydrologysommers=hydrologysommers();

classdef hydrologysommers
	properties (SetAccess=public) 
		head            = NaN;
		gap_height      = NaN;
		bump_spacing    = NaN;
		bump_height     = NaN;
		englacial_input = NaN;
		moulin_input    = NaN;
		reynolds        = NaN;
		spchead         = NaN;
		neumannflux     = NaN;
		relaxation      = 0;
		storage         = 0;
	end
	methods
		function self = extrude(self,md) % {{{
		end % }}}
		function self = hydrologysommers(varargin) % {{{
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
	      % Set under-relaxation parameter to be 1 (no under-relaxation of nonlinear iteration)	
			self.relaxation=1;
			self.storage=0;
		end % }}}
		function md = checkconsistency(self,md,solution,analyses) % {{{

			%Early return
			if ~ismember('HydrologySommersAnalysis',analyses)
				return;
			end

			md = checkfield(md,'fieldname','hydrology.head','size',[md.mesh.numberofvertices 1],'NaN',1,'Inf',1);
			md = checkfield(md,'fieldname','hydrology.gap_height','>=',0,'size',[md.mesh.numberofelements 1],'NaN',1,'Inf',1);
			md = checkfield(md,'fieldname','hydrology.bump_spacing','>',0,'size',[md.mesh.numberofelements 1],'NaN',1,'Inf',1);
			md = checkfield(md,'fieldname','hydrology.bump_height','>=',0,'size',[md.mesh.numberofelements 1],'NaN',1,'Inf',1);
			md = checkfield(md,'fieldname','hydrology.englacial_input','>=',0,'NaN',1,'Inf',1,'timeseries',1);
			md = checkfield(md,'fieldname','hydrology.moulin_input','>=',0,'NaN',1,'Inf',1,'timeseries',1);
			md = checkfield(md,'fieldname','hydrology.reynolds','>',0,'size',[md.mesh.numberofelements 1],'NaN',1,'Inf',1);
			md = checkfield(md,'fieldname','hydrology.neumannflux','timeseries',1,'NaN',1,'Inf',1);
			md = checkfield(md,'fieldname','hydrology.spchead','size',[md.mesh.numberofvertices 1]);	
         md = checkfield(md,'fieldname','hydrology.relaxation','>=',0);	
			md = checkfield(md,'fieldname','hydrology.storage','>=',0);
		end % }}}
		function disp(self) % {{{
			disp(sprintf('   hydrologysommers solution parameters:'));
			fielddisplay(self,'head','subglacial hydrology water head (m)');
			fielddisplay(self,'gap_height','height of gap separating ice to bed (m)');
			fielddisplay(self,'bump_spacing','characteristic bedrock bump spacing (m)');
			fielddisplay(self,'bump_height','characteristic bedrock bump height (m)');
			fielddisplay(self,'englacial_input','liquid water input from englacial to subglacial system (m/yr)');
			fielddisplay(self,'moulin_input','liquid water input from moulins (at the vertices) to subglacial system (m^3/s)');
			fielddisplay(self,'reynolds','Reynolds'' number');
			fielddisplay(self,'neumannflux','water flux applied along the model boundary (m^2/s)');
			fielddisplay(self,'spchead','water head constraints (NaN means no constraint) (m)');
			fielddisplay(self,'relaxation','under-relaxation coefficient for nonlinear iteration');
			fielddisplay(self,'storage','englacial storage coefficient (void ratio)');
		end % }}}
		function marshall(self,prefix,md,fid) % {{{

			yts=md.constants.yts;

			WriteData(fid,prefix,'name','md.hydrology.model','data',3,'format','Integer');
			WriteData(fid,prefix,'object',self,'class','hydrology','fieldname','head','format','DoubleMat','mattype',1);
			WriteData(fid,prefix,'object',self,'class','hydrology','fieldname','gap_height','format','DoubleMat','mattype',2);
			WriteData(fid,prefix,'object',self,'class','hydrology','fieldname','bump_spacing','format','DoubleMat','mattype',2);
			WriteData(fid,prefix,'object',self,'class','hydrology','fieldname','bump_height','format','DoubleMat','mattype',2);
			WriteData(fid,prefix,'object',self,'class','hydrology','fieldname','englacial_input','format','DoubleMat','mattype',1,'scale',1./yts,'timeserieslength',md.mesh.numberofvertices+1,'yts',md.constants.yts);
			WriteData(fid,prefix,'object',self,'class','hydrology','fieldname','moulin_input','format','DoubleMat','mattype',1,'timeserieslength',md.mesh.numberofvertices+1,'yts',md.constants.yts);
			WriteData(fid,prefix,'object',self,'class','hydrology','fieldname','reynolds','format','DoubleMat','mattype',2);
			WriteData(fid,prefix,'object',self,'class','hydrology','fieldname','neumannflux','format','DoubleMat','mattype',2,'timeserieslength',md.mesh.numberofelements+1,'yts',md.constants.yts);
			WriteData(fid,prefix,'object',self,'class','hydrology','fieldname','spchead','format','DoubleMat','mattype',1);
         WriteData(fid,prefix,'object',self,'class','hydrology','fieldname','relaxation','format','Double');
			WriteData(fid,prefix,'object',self,'class','hydrology','fieldname','storage','format','Double');
		end % }}}
	end
end

