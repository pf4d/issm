%SMBpdd Class definition
%
%   Usage:
%      SMBpdd=SMBpdd();

classdef SMBpdd
	properties (SetAccess=public) 
		precipitation             = NaN;
		monthlytemperatures       = NaN;
		desfac                    = 0;
		s0p                       = NaN;
		s0t                       = NaN;
		rlaps                     = 0;
		rlapslgm                  = 0; 
		Pfac                      = NaN;
		Tdiff                     = NaN;
		sealev                    = NaN;
		isdelta18o                = 0;
		ismungsm                  = 0;
		delta18o                  = NaN;
		delta18o_surface          = NaN;
		temperatures_presentday   = NaN;
		temperatures_lgm          = NaN;
		precipitations_presentday = NaN;
		precipitations_lgm        = NaN;
		requested_outputs      = {};
	end
	methods
		function self = SMBpdd(varargin) % {{{
			switch nargin
				case 0
					self=setdefaultparameters(self);
				otherwise
					error('constructor not supported');
			end
		end % }}}
		function self = extrude(self,md) % {{{
			if(self.isdelta18o==0 & self.ismungsm==0),self.precipitation=project3d(md,'vector',self.precipitation,'type','node');end
			if(self.isdelta18o==0 & self.ismungsm==0),self.monthlytemperatures=project3d(md,'vector',self.monthlytemperatures,'type','node');end
			if(self.isdelta18o),self.temperatures_lgm=project3d(md,'vector',self.temperatures_lgm,'type','node');end
			if(self.isdelta18o),self.temperatures_presentday=project3d(md,'vector',self.temperatures_presentday,'type','node');end
			if(self.isdelta18o),self.precipitations_presentday=project3d(md,'vector',self.precipitations_presentday,'type','node');end
			if(self.isdelta18o),self.precipitations_lgm=project3d(md,'vector',self.precipitations_lgm,'type','node');end
			if(self.ismungsm),self.temperatures_lgm=project3d(md,'vector',self.temperatures_lgm,'type','node');end
			if(self.ismungsm),self.temperatures_presentday=project3d(md,'vector',self.temperatures_presentday,'type','node');end
			if(self.ismungsm),self.precipitations_presentday=project3d(md,'vector',self.precipitations_presentday,'type','node');end
			if(self.ismungsm),self.precipitations_lgm=project3d(md,'vector',self.precipitations_lgm,'type','node');end
			self.s0p=project3d(md,'vector',self.s0p,'type','node');
			self.s0t=project3d(md,'vector',self.s0t,'type','node');

		end % }}}
		function list = defaultoutputs(self,md) % {{{
			list = {''};
		end % }}}
		function self = initialize(self,md) % {{{
                    
			if isnan(self.s0p),
				self.s0p=zeros(md.mesh.numberofvertices,1);
				disp('      no SMBpdd.s0p specified: values set as zero');
			end
			if isnan(self.s0t),
				self.s0t=zeros(md.mesh.numberofvertices,1);
				disp('      no SMBpdd.s0t specified: values set as zero');
			end

		end % }}}
		function self = setdefaultparameters(self) % {{{

		  self.isdelta18o = 0;
		  self.ismungsm   = 0;
		  self.desfac     = 0.5;
		  self.rlaps      = 6.5;
		  self.rlapslgm   = 6.5;
                  
		end % }}}
		function md = checkconsistency(self,md,solution,analyses) % {{{

			if ismember('MasstransportAnalysis',analyses),
				md = checkfield(md,'fieldname','smb.desfac','<=',1,'numel',1);
				md = checkfield(md,'fieldname','smb.s0p','>=',0,'NaN',1,'Inf',1,'size',[md.mesh.numberofvertices 1]);
				md = checkfield(md,'fieldname','smb.s0t','>=',0,'NaN',1,'Inf',1,'size',[md.mesh.numberofvertices 1]);
				md = checkfield(md,'fieldname','smb.rlaps','>=',0,'numel',1);
				md = checkfield(md,'fieldname','smb.rlapslgm','>=',0,'numel',1);
				if(self.isdelta18o==0 & self.ismungsm==0)
					md = checkfield(md,'fieldname','smb.monthlytemperatures','timeseries',1,'NaN',1,'Inf',1);
					md = checkfield(md,'fieldname','smb.precipitation','timeseries',1,'NaN',1,'Inf',1);
				elseif(self.isdelta18o==1) 
					md = checkfield(md,'fieldname','smb.delta18o','NaN',1,'Inf',1,'size',[2,NaN],'singletimeseries',1);
					md = checkfield(md,'fieldname','smb.delta18o_surface','NaN',1,'Inf',1,'size',[2,NaN],'singletimeseries',1);
					md = checkfield(md,'fieldname','smb.temperatures_presentday','size',[md.mesh.numberofvertices+1 12],'NaN',1,'Inf',1,'timeseries',1);
					md = checkfield(md,'fieldname','smb.temperatures_lgm','size',[md.mesh.numberofvertices+1 12],'NaN',1,'Inf',1,'timeseries',1);
					md = checkfield(md,'fieldname','smb.precipitations_presentday','size',[md.mesh.numberofvertices+1 12],'NaN',1,'Inf',1,'timeseries',1);
					md = checkfield(md,'fieldname','smb.precipitations_lgm','size',[md.mesh.numberofvertices+1 12],'NaN',1,'Inf',1,'timeseries',1);                                       
					md = checkfield(md,'fieldname','smb.Tdiff','NaN',1,'Inf',1,'size',[2,NaN],'singletimeseries',1);
					md = checkfield(md,'fieldname','smb.sealev','NaN',1,'Inf',1,'size',[2,NaN],'singletimeseries',1);
				elseif(self.ismungsm==1) 
					md = checkfield(md,'fieldname','smb.temperatures_presentday','size',[md.mesh.numberofvertices+1 12],'NaN',1,'Inf',1,'timeseries',1);
					md = checkfield(md,'fieldname','smb.temperatures_lgm','size',[md.mesh.numberofvertices+1 12],'NaN',1,'Inf',1,'timeseries',1);
					md = checkfield(md,'fieldname','smb.precipitations_presentday','size',[md.mesh.numberofvertices+1 12],'NaN',1,'Inf',1,'timeseries',1);
					md = checkfield(md,'fieldname','smb.precipitations_lgm','size',[md.mesh.numberofvertices+1 12],'NaN',1,'Inf',1,'timeseries',1);                                       
					md = checkfield(md,'fieldname','smb.Pfac','NaN',1,'Inf',1,'size',[2,NaN],'singletimeseries',1);
					md = checkfield(md,'fieldname','smb.Tdiff','NaN',1,'Inf',1,'size',[2,NaN],'singletimeseries',1);
					md = checkfield(md,'fieldname','smb.sealev','NaN',1,'Inf',1,'size',[2,NaN],'singletimeseries',1);
				end
			end
			md = checkfield(md,'fieldname','smb.requested_outputs','stringrow',1);
			
		end % }}}
		function disp(self) % {{{
			disp(sprintf('   surface forcings parameters:'));

			disp(sprintf('\n   PDD and deltaO18 parameters:'));
			fielddisplay(self,'isdelta18o','is temperature and precipitation delta18o parametrisation activated (0 or 1, default is 0)');
			fielddisplay(self,'ismungsm','is temperature and precipitation mungsm parametrisation activated (0 or 1, default is 0)');
			fielddisplay(self,'desfac','desertification elevation factor (between 0 and 1, default is 0.5) [m]');
			fielddisplay(self,'s0p','should be set to elevation from precip source (between 0 and a few 1000s m, default is 0) [m]');
			fielddisplay(self,'s0t','should be set to elevation from temperature source (between 0 and a few 1000s m, default is 0) [m]');
			fielddisplay(self,'rlaps','present day lapse rate [degree/km]');
			fielddisplay(self,'rlapslgm','LGM lapse rate [degree/km]');
			if(self.isdelta18o==0 & self.ismungsm==0)
				fielddisplay(self,'monthlytemperatures',['monthly surface temperatures [K], required if pdd is activated and delta18o not activated']);
				fielddisplay(self,'precipitation',['monthly surface precipitation [m/yr water eq], required if pdd is activated and delta18o or mungsm not activated']);
			elseif(self.isdelta18o==1)
				fielddisplay(self,'delta18o','delta18o [per mil], required if pdd is activated and delta18o activated');
				fielddisplay(self,'delta18o_surface','surface elevation of the delta18o site, required if pdd is activated and delta18o activated');
				fielddisplay(self,'temperatures_presentday','monthly present day surface temperatures [K], required if delta18o/mungsm is activated');
				fielddisplay(self,'temperatures_lgm','monthly LGM surface temperatures [K], required if delta18o or mungsm is activated');
				fielddisplay(self,'precipitations_presentday','monthly surface precipitation [m/yr water eq], required if delta18o/mungsm is activated');
				fielddisplay(self,'precipitations_lgm','monthly surface precipitation [m/yr water eq], required if delta18o/mungsm is activated');
				fielddisplay(self,'Tdiff','time interpolation parameter for temperature, 1D(year), required if mungsm is activated');
				fielddisplay(self,'sealev','sea level [m], 1D(year), required if mungsm is activated');
			elseif(self.ismungsm==1)
				fielddisplay(self,'temperatures_presentday','monthly present day surface temperatures [K], required if delta18o/mungsm is activated');
				fielddisplay(self,'temperatures_lgm','monthly LGM surface temperatures [K], required if delta18o or mungsm is activated');
				fielddisplay(self,'precipitations_presentday','monthly surface precipitation [m/yr water eq], required if delta18o/mungsm is activated');
				fielddisplay(self,'precipitations_lgm','monthly surface precipitation [m/yr water eq], required if delta18o/mungsm is activated');
				fielddisplay(self,'Pfac','time interpolation parameter for precipitation, 1D(year), required if mungsm is activated');
				fielddisplay(self,'Tdiff','time interpolation parameter for temperature, 1D(year), required if mungsm is activated');
				fielddisplay(self,'sealev','sea level [m], 1D(year), required if mungsm is activated');
			end
			fielddisplay(self,'requested_outputs','additional outputs requested');
		end % }}}
		function marshall(self,prefix,md,fid) % {{{

			yts=md.constants.yts;

			WriteData(fid,prefix,'name','md.smb.model','data',4,'format','Integer');

			WriteData(fid,prefix,'object',self,'class','smb','fieldname','isdelta18o','format','Boolean');
			WriteData(fid,prefix,'object',self,'class','smb','fieldname','ismungsm','format','Boolean');
			WriteData(fid,prefix,'object',self,'class','smb','fieldname','desfac','format','Double');
			WriteData(fid,prefix,'object',self,'class','smb','fieldname','s0p','format','DoubleMat','mattype',1);
			WriteData(fid,prefix,'object',self,'class','smb','fieldname','s0t','format','DoubleMat','mattype',1);
			WriteData(fid,prefix,'object',self,'class','smb','fieldname','rlaps','format','Double');
			WriteData(fid,prefix,'object',self,'class','smb','fieldname','rlapslgm','format','Double');

			if(self.isdelta18o==0 & self.ismungsm==0)
				%WriteData(fid,prefix,'object',self,'class','smb','fieldname','monthlytemperatures','format','DoubleMat','mattype',1);
				WriteData(fid,prefix,'object',self,'class','smb','fieldname','monthlytemperatures','format','DoubleMat','mattype',1,'timeserieslength',md.mesh.numberofvertices+1,'yts',md.constants.yts);
				WriteData(fid,prefix,'object',self,'class','smb','fieldname','precipitation','format','DoubleMat','mattype',1,'scale',1./yts,'timeserieslength',md.mesh.numberofvertices+1,'yts',md.constants.yts);
			elseif self.isdelta18o
				WriteData(fid,prefix,'object',self,'class','smb','fieldname','temperatures_presentday','format','DoubleMat','mattype',1,'timeserieslength',md.mesh.numberofvertices+1,'yts',md.constants.yts);
				WriteData(fid,prefix,'object',self,'class','smb','fieldname','temperatures_lgm','format','DoubleMat','mattype',1,'timeserieslength',md.mesh.numberofvertices+1,'yts',md.constants.yts);
				WriteData(fid,prefix,'object',self,'class','smb','fieldname','precipitations_presentday','format','DoubleMat','mattype',1,'scale',1./yts,'timeserieslength',md.mesh.numberofvertices+1,'yts',md.constants.yts);
				WriteData(fid,prefix,'object',self,'class','smb','fieldname','precipitations_lgm','format','DoubleMat','mattype',1,'scale',1./yts,'timeserieslength',md.mesh.numberofvertices+1,'yts',md.constants.yts);
				WriteData(fid,prefix,'object',self,'class','smb','fieldname','delta18o_surface','format','DoubleMat','mattype',1,'timeserieslength',2,'yts',md.constants.yts);
				WriteData(fid,prefix,'object',self,'class','smb','fieldname','delta18o','format','DoubleMat','mattype',1,'timeserieslength',2,'yts',md.constants.yts);
				WriteData(fid,prefix,'object',self,'class','smb','fieldname','Tdiff','format','DoubleMat','mattype',1,'timeserieslength',2,'yts',md.constants.yts);
				WriteData(fid,prefix,'object',self,'class','smb','fieldname','sealev','format','DoubleMat','mattype',1,'timeserieslength',2,'yts',md.constants.yts);
			elseif self.ismungsm
				WriteData(fid,prefix,'object',self,'class','smb','fieldname','temperatures_presentday','format','DoubleMat','mattype',1,'timeserieslength',md.mesh.numberofvertices+1,'yts',md.constants.yts);
				WriteData(fid,prefix,'object',self,'class','smb','fieldname','temperatures_lgm','format','DoubleMat','mattype',1,'timeserieslength',md.mesh.numberofvertices+1,'yts',md.constants.yts);
				WriteData(fid,prefix,'object',self,'class','smb','fieldname','precipitations_presentday','format','DoubleMat','mattype',1,'scale',1./yts,'timeserieslength',md.mesh.numberofvertices+1,'yts',md.constants.yts);
				WriteData(fid,prefix,'object',self,'class','smb','fieldname','precipitations_lgm','format','DoubleMat','mattype',1,'scale',1./yts,'timeserieslength',md.mesh.numberofvertices+1,'yts',md.constants.yts);
				WriteData(fid,prefix,'object',self,'class','smb','fieldname','Pfac','format','DoubleMat','mattype',1,'timeserieslength',2,'yts',md.constants.yts);
				WriteData(fid,prefix,'object',self,'class','smb','fieldname','Tdiff','format','DoubleMat','mattype',1,'timeserieslength',2,'yts',md.constants.yts);
				WriteData(fid,prefix,'object',self,'class','smb','fieldname','sealev','format','DoubleMat','mattype',1,'timeserieslength',2,'yts',md.constants.yts);
			end
			
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
