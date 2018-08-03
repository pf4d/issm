%SLR class definition
%
%   Usage:
%      slr=slr();

classdef slr
	properties (SetAccess=public) 
		deltathickness = NaN;
		sealevel       = NaN; 
		maxiter        = 0;
		reltol         = 0;
		abstol         = 0;
		love_h         = 0; %provided by PREM model
		love_k         = 0; %ideam
		love_l         = 0; %ideam
		tide_love_k    = 0; %ideam
		tide_love_h    = 0; %ideam
		fluid_love     = 0; 
		equatorial_moi = 0; 
		polar_moi		= 0; 
		angular_velocity = 0;
		rigid          = 0;
		elastic        = 0;
		rotation       = 0;
		ocean_area_scaling = 0;
		degacc         = 0;
		requested_outputs      = {};
		transitions    = {};
	end
	methods
		function self = slr(varargin) % {{{
			switch nargin
				case 0
					self=setdefaultparameters(self);
				otherwise
					error('constructor not supported');
			end
		end % }}}
		function self = setdefaultparameters(self) % {{{
		
		%Convergence criterion: absolute, relative and residual
		self.reltol=NaN; %default
		self.abstol=0.001; %1 mm of sea level rise

		%maximum of non-linear iterations.
		self.maxiter=10;

		%computational flags: 
		self.rigid=1;
		self.elastic=1;
		self.rotation=0;
		self.ocean_area_scaling=0;

		%tidal love numbers: 
		self.tide_love_h=0.6149; %degree 2
		self.tide_love_k=0.3055; % degree 2
	
		%secular fluid love number: 
		self.fluid_love=0.942; 
		
		%moment of inertia: 
		self.equatorial_moi=8.0077*10^37; % [kg m^2] 
		self.polar_moi		 =8.0345*10^37; % [kg m^2] 

		% mean rotational velocity of earth 
		self.angular_velocity=7.2921*10^-5; % [s^-1] 

		%numerical discretization accuracy
		self.degacc=.01;
		
		%output default:
		self.requested_outputs={'default'};

		%transitions should be a cell array of vectors: 
		self.transitions={};
		
		end % }}}
		function md = checkconsistency(self,md,solution,analyses) % {{{

			if ~ismember('SealevelriseAnalysis',analyses), return; end
			md = checkfield(md,'fieldname','slr.deltathickness','NaN',1,'Inf',1,'size',[md.mesh.numberofelements 1]);
			md = checkfield(md,'fieldname','slr.sealevel','NaN',1,'Inf',1,'size',[md.mesh.numberofvertices 1]);
			md = checkfield(md,'fieldname','slr.love_h','NaN',1,'Inf',1);
			md = checkfield(md,'fieldname','slr.love_k','NaN',1,'Inf',1);
			md = checkfield(md,'fieldname','slr.love_l','NaN',1,'Inf',1);
			md = checkfield(md,'fieldname','slr.tide_love_h','NaN',1,'Inf',1);
			md = checkfield(md,'fieldname','slr.tide_love_k','NaN',1,'Inf',1);
			md = checkfield(md,'fieldname','slr.fluid_love','NaN',1,'Inf',1);
			md = checkfield(md,'fieldname','slr.equatorial_moi','NaN',1,'Inf',1);
			md = checkfield(md,'fieldname','slr.polar_moi','NaN',1,'Inf',1);
			md = checkfield(md,'fieldname','slr.angular_velocity','NaN',1,'Inf',1);
			md = checkfield(md,'fieldname','slr.reltol','size',[1 1]);
			md = checkfield(md,'fieldname','slr.abstol','size',[1 1]);
			md = checkfield(md,'fieldname','slr.maxiter','size',[1 1],'>=',1);
			md = checkfield(md,'fieldname','slr.degacc','size',[1 1],'>=',1e-10);
			md = checkfield(md,'fieldname','slr.requested_outputs','stringrow',1);

			%check that love numbers are provided at the same level of accuracy: 
			if (size(self.love_h,1)~=size(self.love_k,1) | size(self.love_h,1)~=size(self.love_l,1)),
				error('slr error message: love numbers should be provided at the same level of accuracy');
			end

			%cross check that whereever we have an ice load, the mask is <0 on each vertex: 
			pos=find(self.deltathickness);
			maskpos=md.mask.ice_levelset(md.mesh.elements(pos,:)); 
			[els,vertices]=find(maskpos>0);
			if length(els),
				error('slr checkconsistency fail: there are elements with ice loads where some vertices are not on the ice!');
			end

		end % }}}
		function list=defaultoutputs(self,md) % {{{
			list = {'Sealevel'};
		end % }}}
		function disp(self) % {{{
			disp(sprintf('   slr parameters:'));

			fielddisplay(self,'deltathickness','thickness change (main loading of the slr solution core [m]');
			fielddisplay(self,'sealevel','current sea level (prior to computation) [m]');
			fielddisplay(self,'reltol','sea level rise relative convergence criterion, (default, NaN: not applied)');
			fielddisplay(self,'abstol','sea level rise absolute convergence criterion, NaN: not applied');
			fielddisplay(self,'maxiter','maximum number of nonlinear iterations');
			fielddisplay(self,'love_h','load Love number for radial displacement');
			fielddisplay(self,'love_k','load Love number for gravitational potential perturbation');
			fielddisplay(self,'love_l','load Love number for horizontal displacements');
			fielddisplay(self,'tide_love_k','tidal load Love number (deg 2)');
			fielddisplay(self,'tide_love_h','tidal load Love number (deg 2)');
			fielddisplay(self,'fluid_love','secular fluid Love number');
			fielddisplay(self,'equatorial_moi','mean equatorial moment of inertia [kg m^2]');
			fielddisplay(self,'polar_moi','polar moment of inertia [kg m^2]');
			fielddisplay(self,'angular_velocity','mean rotational velocity of earth [per second]'); 
			fielddisplay(self,'rigid','rigid earth graviational potential perturbation');
			fielddisplay(self,'elastic','elastic earth graviational potential perturbation');
			fielddisplay(self,'rotation','earth rotational potential perturbation');
			fielddisplay(self,'ocean_area_scaling','correction for model representation of ocean area [default: No correction]'); 
			fielddisplay(self,'degacc','accuracy (default .01 deg) for numerical discretization of the Green''s functions');
			fielddisplay(self,'transitions','indices into parts of the mesh that will be icecaps');
			fielddisplay(self,'requested_outputs','additional outputs requested');

		end % }}}
		function marshall(self,prefix,md,fid) % {{{
			WriteData(fid,prefix,'object',self,'fieldname','deltathickness','format','DoubleMat','mattype',2);
			%WriteData(fid,prefix,'object',self,'fieldname','deltathickness','format','DoubleMat','mattype',1,'timeserieslength',md.mesh.numberofelements+1);
			WriteData(fid,prefix,'object',self,'fieldname','sealevel','mattype',1,'format','DoubleMat','timeserieslength',md.mesh.numberofvertices+1,'yts',md.constants.yts);
			WriteData(fid,prefix,'object',self,'fieldname','reltol','format','Double');
			WriteData(fid,prefix,'object',self,'fieldname','abstol','format','Double');
			WriteData(fid,prefix,'object',self,'fieldname','maxiter','format','Integer');
			WriteData(fid,prefix,'object',self,'fieldname','love_h','format','DoubleMat','mattype',1);
			WriteData(fid,prefix,'object',self,'fieldname','love_k','format','DoubleMat','mattype',1);
			WriteData(fid,prefix,'object',self,'fieldname','love_l','format','DoubleMat','mattype',1);
			WriteData(fid,prefix,'object',self,'fieldname','tide_love_h','format','Double');
			WriteData(fid,prefix,'object',self,'fieldname','tide_love_k','format','Double');
			WriteData(fid,prefix,'object',self,'fieldname','fluid_love','format','Double');
			WriteData(fid,prefix,'object',self,'fieldname','equatorial_moi','format','Double');
			WriteData(fid,prefix,'object',self,'fieldname','polar_moi','format','Double');
			WriteData(fid,prefix,'object',self,'fieldname','angular_velocity','format','Double');
			WriteData(fid,prefix,'object',self,'fieldname','rigid','format','Boolean');
			WriteData(fid,prefix,'object',self,'fieldname','elastic','format','Boolean');
			WriteData(fid,prefix,'object',self,'fieldname','rotation','format','Boolean');
			WriteData(fid,prefix,'object',self,'fieldname','ocean_area_scaling','format','Boolean');
			WriteData(fid,prefix,'object',self,'fieldname','degacc','format','Double');
			WriteData(fid,prefix,'object',self,'fieldname','transitions','format','MatArray');
			
			%process requested outputs
			outputs = self.requested_outputs;
			pos  = find(ismember(outputs,'default'));
			if ~isempty(pos),
				outputs(pos) = [];                         %remove 'default' from outputs
				outputs      = [outputs defaultoutputs(self,md)]; %add defaults
			end
			WriteData(fid,prefix,'data',outputs,'name','md.slr.requested_outputs','format','StringArray');

		end % }}}
		function savemodeljs(self,fid,modelname) % {{{
		
			writejs1Darray(fid,[modelname '.slr.deltathickness'],self.deltathickness);
			writejs1Darray(fid,[modelname '.slr.sealevel'],self.sealevel);
			writejsdouble(fid,[modelname '.slr.maxiter'],self.maxiter);
			writejsdouble(fid,[modelname '.slr.reltol'],self.reltol);
			writejsdouble(fid,[modelname '.slr.abstol'],self.abstol);
			writejs1Darray(fid,[modelname '.slr.love_h'],self.love_h);
			writejs1Darray(fid,[modelname '.slr.love_k'],self.love_k);
			writejs1Darray(fid,[modelname '.slr.love_l'],self.love_l);
			writejsdouble(fid,[modelname '.slr.tide_love_k'],self.tide_love_k);
			writejsdouble(fid,[modelname '.slr.tide_love_h'],self.tide_love_h);
			writejsdouble(fid,[modelname '.slr.fluid_love'],self.fluid_love);
			writejsdouble(fid,[modelname '.slr.equatorial_moi'],self.equatorial_moi);
			writejsdouble(fid,[modelname '.slr.polar_moi'],self.polar_moi);
			writejsdouble(fid,[modelname '.slr.angular_velocity'],self.angular_velocity);
			writejsdouble(fid,[modelname '.slr.rigid'],self.rigid);
			writejsdouble(fid,[modelname '.slr.elastic'],self.elastic);
			writejsdouble(fid,[modelname '.slr.rotation'],self.rotation);
			writejsdouble(fid,[modelname '.slr.ocean_area_scaling'],self.ocean_area_scaling);
			writejsdouble(fid,[modelname '.slr.degacc'],self.degacc);
			writejscellstring(fid,[modelname '.slr.requested_outputs'],self.requested_outputs);
			writejscellarray(fid,[modelname '.slr.transitions'],self.transitions);
		end % }}}
	end
end
