%MATERIALS class definition
%
%   Usage:
%      materials=materials();

classdef materials < dynamicprops
	properties (SetAccess=public) 
		nature={};
		%all properties are dynamic.
	end
	methods
		function self = materials(varargin) % {{{
			if nargin==0
				self.nature={'ice'};
			else 
				self.nature=varargin;
			end
			
			%check this is acceptable: 
			for i=1:length(self.nature),
				if ~(strcmpi(self.nature{i},'litho') | strcmpi(self.nature{i},'ice')), 
					error('materials constructor error message: nature of the material not supported yet! (''ice'' or ''litho'')');
				end
			end
			
			%start filling in the dynamic fields: 
			for i=1:length(self.nature),
				nat=self.nature{i}; 
				switch nat
				case 'ice'
					self.addprop('rho_ice');
					self.addprop('rho_water');
					self.addprop('rho_freshwater');
					self.addprop('mu_water');
					self.addprop('heatcapacity');
					self.addprop('latentheat');
					self.addprop('thermalconductivity');
					self.addprop('temperateiceconductivity');
				    self.addprop('meltingpoint');
					self.addprop('beta');
					self.addprop('mixed_layer_capacity');
					self.addprop('thermal_exchange_velocity');
					self.addprop('rheology_B');
					self.addprop('rheology_n');
					self.addprop('rheology_law');
				case 'litho'
					self.addprop('numlayers');
				    self.addprop('radius');
					self.addprop('viscosity');
					self.addprop('lame_lambda');
					self.addprop('lame_mu');
					self.addprop('burgers_viscosity');
					self.addprop('burgers_mu');
					self.addprop('isburgers');
					self.addprop('density');
					self.addprop('issolid');
				otherwise
					error('materials constructor error message: nature of the material not supported yet! (''ice'' or ''litho'')');
				end
			end
			%set default parameters:
			self.setdefaultparameters();

		end % }}}
		function self = setdefaultparameters(self) % {{{

			for i=1:length(self.nature),
				nat=self.nature{i}; 
				switch nat
				case 'ice'
					%ice density (kg/m^3)
					self.rho_ice=917.;

					%ocean water density (kg/m^3)
					self.rho_water=1023.;

					%fresh water density (kg/m^3)
					self.rho_freshwater=1000.;

					%water viscosity (N.s/m^2)
					self.mu_water=0.001787;  

					%ice heat capacity cp (J/kg/K)
					self.heatcapacity=2093.;

					%ice latent heat of fusion L (J/kg)
					self.latentheat=3.34*10^5;

					%ice thermal conductivity (W/m/K)
					self.thermalconductivity=2.4;
					
					%wet ice thermal conductivity (W/m/K)
					self.temperateiceconductivity=.24;

					%the melting point of ice at 1 atmosphere of pressure in K
					self.meltingpoint=273.15;

					%rate of change of melting point with pressure (K/Pa)
					self.beta=9.8*10^-8;

					%mixed layer (ice-water interface) heat capacity (J/kg/K)
					self.mixed_layer_capacity=3974.;

					%thermal exchange velocity (ice-water interface) (m/s)
					self.thermal_exchange_velocity=1.00*10^-4;

					%Rheology law: what is the temperature dependence of B with T
					%available: none, paterson and arrhenius
					self.rheology_law='Paterson';

				case 'litho'
					%we default to a configuration that enables running GIA solutions using giacaron and/or giaivins. 
					self.numlayers=2;

					%surface, then the lab (lithosphere/asthenosphere boundary) then the center of the earth 
					%(with 1d3 to avoid numerical singularities) 
					self.radius=[6378*1e3;6278*1e3;1e3]; 

					self.viscosity=[1e40;1e21]; %lithosphere and mantle viscosity (respectively) [Pa.s]
					self.lame_mu=[6.7*10^10;1.45*1e11];  % (Pa) %lithosphere and mantle shear modulus (respectively) [Pa]
					self.lame_lambda=[6.7*10^10;1.45*1e11];  % (Pa) %lithosphere and mantle lamba parameter (respectively) [Pa]
					self.burgers_viscosity=[NaN;NaN];
					self.burgers_mu=[NaN;NaN];
					self.isburgers=[false;false];
					self.density=[3.32*1e3;3.34*1e3];  % (Pa) %lithosphere and mantle density [kg/m^3]
					self.issolid=[true;true]; % is layer solid or liquid.

				otherwise
					error('materials constructor error message: nature of the material not supported yet! (''ice'' or ''litho'')');
				end
			end
		end % }}}
		function disp(self) % {{{
			disp(sprintf('   Materials:'));

			for i=1:length(self.nature),
				nat=self.nature{i}; 
				switch nat
				case 'ice'
					disp(sprintf('   \nIce:'));
					fielddisplay(self,'rho_ice','ice density [kg/m^3]');
					fielddisplay(self,'rho_water','ocean water density [kg/m^3]');
					fielddisplay(self,'rho_freshwater','fresh water density [kg/m^3]');
					fielddisplay(self,'mu_water','water viscosity [N s/m^2]');
					fielddisplay(self,'heatcapacity','heat capacity [J/kg/K]');
					fielddisplay(self,'thermalconductivity',['ice thermal conductivity [W/m/K]']);
					fielddisplay(self,'temperateiceconductivity','temperate ice thermal conductivity [W/m/K]');
					fielddisplay(self,'meltingpoint','melting point of ice at 1atm in K');
					fielddisplay(self,'latentheat','latent heat of fusion [J/kg]');
					fielddisplay(self,'beta','rate of change of melting point with pressure [K/Pa]');
					fielddisplay(self,'mixed_layer_capacity','mixed layer capacity [W/kg/K]');
					fielddisplay(self,'thermal_exchange_velocity','thermal exchange velocity [m/s]');
					fielddisplay(self,'rheology_B','flow law parameter [Pa/s^(1/n)]');
					fielddisplay(self,'rheology_n','Glen''s flow law exponent');
					fielddisplay(self,'rheology_law',['law for the temperature dependance of the rheology: ''None'', ''BuddJacka'', Cuffey'', ''CuffeyTemperate'', ''Paterson'', ''Arrhenius'' or ''LliboutryDuval''']);
				case 'litho'
					disp(sprintf('   \nLitho:'));
					fielddisplay(self,'numlayers','number of layers (default 2)');
					fielddisplay(self,'radius','array describing the radius for each interface (numlayers+1) [m]');
					fielddisplay(self,'viscosity','array describing each layer''s viscosity (numlayers) [Pa.s]');
					fielddisplay(self,'lame_lambda','array describing the lame lambda parameter (numlayers) [Pa]');
					fielddisplay(self,'lame_mu','array describing the shear modulus for each layers (numlayers) [Pa]');
					fielddisplay(self,'burgers_viscosity','array describing each layer''s transient viscosity, only for Burgers rheologies  (numlayers) [Pa.s]');
					fielddisplay(self,'burgers_mu','array describing each layer''s transient shear modulus, only for Burgers rheologies  (numlayers) [Pa]');
					fielddisplay(self,'isburgers','array describing whether we adopt a MaxWell (0) or Burgers (1) rheology (default 0)');
					fielddisplay(self,'density','array describing each layer''s density (numlayers) [kg/m^3]');
					fielddisplay(self,'issolid','array describing whether the layer is solid or liquid (default 1) (numlayers)');
				otherwise
					error('materials constructor error message: nature of the material not supported yet! (''ice'' or ''litho'')');
				end
			end
		end % }}}
		function md = checkconsistency(self,md,solution,analyses) % {{{

			for i=1:length(self.nature),
				nat=self.nature{i}; 
				switch nat
				case 'ice'
					md = checkfield(md,'fieldname','materials.rho_ice','>',0);
					md = checkfield(md,'fieldname','materials.rho_water','>',0);
					md = checkfield(md,'fieldname','materials.rho_freshwater','>',0);
					md = checkfield(md,'fieldname','materials.mu_water','>',0);
					md = checkfield(md,'fieldname','materials.rheology_B','>',0,'timeseries',1,'NaN',1,'Inf',1);
					md = checkfield(md,'fieldname','materials.rheology_n','>',0,'size',[md.mesh.numberofelements 1]);
					md = checkfield(md,'fieldname','materials.rheology_law','values',{'None' 'BuddJacka' 'Cuffey' 'CuffeyTemperate' 'Paterson' 'Arrhenius' 'LliboutryDuval'});
				case 'litho'
					if ~ismember('LoveAnalysis',analyses), return; end
					md = checkfield(md,'fieldname','materials.numlayers','NaN',1,'Inf',1,'>',0,'numel',1);
					md = checkfield(md,'fieldname','materials.radius','NaN',1,'Inf',1,'size',[md.materials.numlayers+1 1],'>',0);
					md = checkfield(md,'fieldname','materials.lame_mu','NaN',1,'Inf',1,'size',[md.materials.numlayers 1],'>',0);
					md = checkfield(md,'fieldname','materials.lame_lambda','NaN',1,'Inf',1,'size',[md.materials.numlayers 1],'>',0);
					md = checkfield(md,'fieldname','materials.issolid','NaN',1,'Inf',1,'size',[md.materials.numlayers 1],'>',0,'<',2);
					md = checkfield(md,'fieldname','materials.density','NaN',1,'Inf',1,'size',[md.materials.numlayers 1],'>',0);
					md = checkfield(md,'fieldname','materials.viscosity','NaN',1,'Inf',1,'size',[md.materials.numlayers 1],'>',0);
					md = checkfield(md,'fieldname','materials.isburgers','NaN',1,'Inf',1,'size',[md.materials.numlayers 1],'>',0,'<',2);
					md = checkfield(md,'fieldname','materials.burgers_viscosity','NaN',1,'Inf',1,'size',[md.materials.numlayers 1],'>',0);
					md = checkfield(md,'fieldname','materials.burgers_mu','NaN',1,'Inf',1,'size',[md.materials.numlayers 1],'>',0);

				otherwise
					error('materials constructor error message: nature of the material not supported yet! (''ice'' or ''litho'')');
				end
			end

		end % }}}
		function intnat = naturetointeger(strnat) % {{{
			intnat=zeros(length(strnat),1);
			for i=1:length(strnat),
				switch strnat{i},
				case 'damageice'
					intnat(i)=1; 
				case 'estar'
					intnat(i)=2; 
				case 'ice'
					intnat(i)=3; 
				case 'enhancedice'
					intnat(i)=4; 
				case 'litho'
					intnat(i)=5;
				otherwise
					error('materials constructor error message: nature of the material not supported yet! (''ice'' or ''litho'')');
				end
			end
		end % }}}
		function marshall(self,prefix,md,fid) % {{{

			%1: MatdamageiceEnum 2: MatestarEnum 3: MaticeEnum 4: MatenhancediceEnum 5: MaterialsEnum 
			WriteData(fid,prefix,'name','md.materials.type','data',6,'format','Integer');
			WriteData(fid,prefix,'name','md.materials.nature','data',naturetointeger(self.nature),'format','IntMat','mattype',1);

			for i=1:length(self.nature),
				nat=self.nature{i}; 
				switch nat
				case 'ice'
					WriteData(fid,prefix,'object',self,'class','materials','fieldname','rho_ice','format','Double');
					WriteData(fid,prefix,'object',self,'class','materials','fieldname','rho_water','format','Double');
					WriteData(fid,prefix,'object',self,'class','materials','fieldname','rho_freshwater','format','Double');
					WriteData(fid,prefix,'object',self,'class','materials','fieldname','mu_water','format','Double');
					WriteData(fid,prefix,'object',self,'class','materials','fieldname','heatcapacity','format','Double');
					WriteData(fid,prefix,'object',self,'class','materials','fieldname','latentheat','format','Double');
					WriteData(fid,prefix,'object',self,'class','materials','fieldname','thermalconductivity','format','Double');
					WriteData(fid,prefix,'object',self,'class','materials','fieldname','temperateiceconductivity','format','Double');
					WriteData(fid,prefix,'object',self,'class','materials','fieldname','meltingpoint','format','Double');
					WriteData(fid,prefix,'object',self,'class','materials','fieldname','beta','format','Double');
					WriteData(fid,prefix,'object',self,'class','materials','fieldname','mixed_layer_capacity','format','Double');
					WriteData(fid,prefix,'object',self,'class','materials','fieldname','thermal_exchange_velocity','format','Double');
					WriteData(fid,prefix,'object',self,'class','materials','fieldname','rheology_B','format','DoubleMat','mattype',1,'timeserieslength',md.mesh.numberofvertices+1,'yts',md.constants.yts);
					WriteData(fid,prefix,'object',self,'class','materials','fieldname','rheology_n','format','DoubleMat','mattype',2);
					WriteData(fid,prefix,'data',self.rheology_law,'name','md.materials.rheology_law','format','String');
				case 'litho'
					WriteData(fid,prefix,'object',self,'class','materials','fieldname','numlayers','format','Double');
					WriteData(fid,prefix,'object',self,'class','materials','fieldname','radius','format','DoubleMat','mattype',1);
					WriteData(fid,prefix,'object',self,'class','materials','fieldname','lame_mu','format','DoubleMat','mattype',1);
					WriteData(fid,prefix,'object',self,'class','materials','fieldname','lame_lambda ','format','DoubleMat','mattype',1);
					WriteData(fid,prefix,'object',self,'class','materials','fieldname','issolid','format','BooleanMat','mattype',1);
					WriteData(fid,prefix,'object',self,'class','materials','fieldname','density','format','DoubleMat','mattype',1); 
					WriteData(fid,prefix,'object',self,'class','materials','fieldname','viscosity','format','DoubleMat','mattype',1); 
					WriteData(fid,prefix,'object',self,'class','materials','fieldname','isburgers','format','BooleanMat','mattype',1); 
					WriteData(fid,prefix,'object',self,'class','materials','fieldname','burgers_viscosity','format','DoubleMat','mattype',1); 
					WriteData(fid,prefix,'object',self,'class','materials','fieldname','burgers_mu','format','DoubleMat','mattype',1); 
				otherwise
					error('materials constructor error message: nature of the material not supported yet! (''ice'' or ''litho'')');
				end
			end
		end % }}}
		function self = extrude(self,md) % {{{
			for i=1:length(self.nature),
				nat=self.nature{i}; 
				switch nat
				case 'ice'
					self.rheology_B=project3d(md,'vector',self.rheology_B,'type','node');
					self.rheology_n=project3d(md,'vector',self.rheology_n,'type','element');
				end
			end
		end % }}}
		function savemodeljs(self,fid,modelname) % {{{
	
			for i=1:length(self.nature),
				nat=self.nature{i}; 
				switch nat
				case 'ice'
					writejsdouble(fid,[modelname '.materials.rho_ice'],self.rho_ice);
					writejsdouble(fid,[modelname '.materials.rho_water'],self.rho_water);
					writejsdouble(fid,[modelname '.materials.rho_freshwater'],self.rho_freshwater);
					writejsdouble(fid,[modelname '.materials.mu_water'],self.mu_water);
					writejsdouble(fid,[modelname '.materials.heatcapacity'],self.heatcapacity);
					writejsdouble(fid,[modelname '.materials.latentheat'],self.latentheat);
					writejsdouble(fid,[modelname '.materials.thermalconductivity'],self.thermalconductivity);
					writejsdouble(fid,[modelname '.materials.temperateiceconductivity'],self.temperateiceconductivity);
					writejsdouble(fid,[modelname '.materials.meltingpoint'],self.meltingpoint);
					writejsdouble(fid,[modelname '.materials.beta'],self.beta);
					writejsdouble(fid,[modelname '.materials.mixed_layer_capacity'],self.mixed_layer_capacity);
					writejsdouble(fid,[modelname '.materials.thermal_exchange_velocity'],self.thermal_exchange_velocity);
					writejsdouble(fid,[modelname '.materials.mixed_layer_capacity'],self.mixed_layer_capacity);
					writejs1Darray(fid,[modelname '.materials.rheology_B'],self.rheology_B);
					writejs1Darray(fid,[modelname '.materials.rheology_n'],self.rheology_n);
					writejsstring(fid,[modelname '.materials.rheology_law'],self.rheology_law);
				case 'litho'
					writejsdouble(fid,[modelname '.materials.numlayers'],self.numlayers);
					writejsdouble(fid,[modelname '.materials.radius'],self.radius);
					writejsdouble(fid,[modelname '.materials.lame_mu'],self.lame_mu);
					writejsdouble(fid,[modelname '.materials.lame_lambda'],self.lame_lambda);
					writejsdouble(fid,[modelname '.materials.issolid'],self.issolid);
					writejsdouble(fid,[modelname '.materials.density'],self.density); 
					writejsdouble(fid,[modelname '.materials.viscosity'],self.viscosity); 
					writejsdouble(fid,[modelname '.materials.isburgers'],self.isburgers); 
					writejsdouble(fid,[modelname '.materials.burgers_viscosity'],self.burgers_viscosity); 
					writejsdouble(fid,[modelname '.materials.burgers_mu'],self.burgers_mu); 
				otherwise
					error('materials constructor error message: nature of the material not supported yet! (''ice'' or ''litho'')');
				end
			end



		end % }}}
	end
end


