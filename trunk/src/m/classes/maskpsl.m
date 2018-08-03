%MASKPSL class definition
%
%   Usage:
%      maskpsl=maskpsl();

classdef maskpsl 
	properties (SetAccess=public) 
		groundedice_levelset = NaN;
		ice_levelset         = NaN;
		ocean_levelset = NaN;
		land_levelset = NaN;
	end
	methods (Static)
		function self = loadobj(self) % {{{
			% This function is directly called by matlab when a model object is
			% loaded. Update old properties here

			%2014 February 5th
			if numel(self.ice_levelset)>1 & all(self.ice_levelset>=0),
				disp('WARNING: md.mask.ice_levelset>=0, you probably need to change the sign of this levelset');
			end

		end% }}}
	end
	methods
		function self = extrude(self,md) % {{{
			self.groundedice_levelset=project3d(md,'vector',self.groundedice_levelset,'type','node');
			self.ice_levelset=project3d(md,'vector',self.ice_levelset,'type','node');
			self.ocean_levelset=project3d(md,'vector',self.ocean_levelset,'type','node');
			self.land_levelset=project3d(md,'vector',self.land_levelset,'type','node');
		end % }}}
		function self = mask(varargin) % {{{
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

			md = checkfield(md,'fieldname','mask.groundedice_levelset','size',[md.mesh.numberofvertices 1]);
			md = checkfield(md,'fieldname','mask.ice_levelset'        ,'size',[md.mesh.numberofvertices 1]);
			md = checkfield(md,'fieldname','mask.ocean_levelset','size',[md.mesh.numberofvertices 1]);
			md = checkfield(md,'fieldname','mask.land_levelset','size',[md.mesh.numberofvertices 1]);
			isice=(md.mask.ice_levelset<=0);
			if sum(isice)==0,
				warning('no ice present in the domain');
			end
			if max(md.mask.ice_levelset)<0,
				warning('no ice front provided');
			end
			icefront=sum(md.mask.ice_levelset(md.mesh.elements)==0,2);
			if (max(icefront)==3 & strcmp(elementtype(md.mesh),'Tria')) | (max(icefront==6) & strcmp(elementtype(md.mesh),'Penta')),
				error('At least one element has all nodes on ice front, change md.mask.ice_levelset to fix it')
			end
		end % }}}
		function disp(self) % {{{
			disp(sprintf('   masks:'));

			fielddisplay(self,'groundedice_levelset','is ice grounded ? grounded ice if > 0, grounding line position if = 0, floating ice if < 0');
			fielddisplay(self,'ice_levelset','presence of ice if < 0, icefront position if = 0, no ice if > 0');
			fielddisplay(self,'ocean_levelset','is the vertex on the ocean ? yes if = 1, no if = 0');
			fielddisplay(self,'land_levelset','is the vertex on the land ? yes if = 1, no if = 0');
		end % }}}
		function marshall(self,prefix,md,fid) % {{{
			WriteData(fid,prefix,'object',self,'class','mask','fieldname','groundedice_levelset','format','DoubleMat','mattype',1);
			WriteData(fid,prefix,'object',self,'class','mask','fieldname','ice_levelset','format','DoubleMat','mattype',1);
			WriteData(fid,prefix,'object',self,'class','mask','fieldname','ocean_levelset','format','DoubleMat','mattype',1);
			WriteData(fid,prefix,'object',self,'class','mask','fieldname','land_levelset','format','DoubleMat','mattype',1);
		end % }}}
		function savemodeljs(self,fid,modelname) % {{{
		
			fprintf(fid,'%s.mask=new maskpsl();\n',modelname);
			writejs1Darray(fid,[modelname '.mask.groundedice_levelset'],self.groundedice_levelset);
			writejs1Darray(fid,[modelname '.mask.ice_levelset'],self.ice_levelset);
			writejs1Darray(fid,[modelname '.mask.ocean_levelset'],self.ocean_levelset);
			writejs1Darray(fid,[modelname '.mask.land_levelset'],self.land_levelset);

		end % }}}
	end
end
