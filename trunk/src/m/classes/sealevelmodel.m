%SEALEVELMODEL class definition
%
%   Usage:
%      slm = sealevelmodel(varargin)
%
%      where varargin is a variable list of options: 
%
%   Example: 
%      slm = sealevel('icecap',md_greenland,'icecap',md_antarctica,'earth',md_earth);

classdef sealevelmodel
	properties (SetAccess=public) %Model fields
		% {{{
		icecaps          = {}; % list of ice cap models
		earth            = 0;  % model for the whole earth
		cluster          = 0;
		miscellaneous    = 0;
		settings         = 0;
		private          = 0;
		%}}}
	end
	methods
		function slm = sealevelmodel(varargin) % {{{

			if nargin==0, 
				slm=setdefaultparameters(slm);
			else 
				slm=setdefaultparameters(slm);

				options=pairoptions(varargin{:}); 
			
				%recover all the icecap models: 
				slm.icecaps=getfieldvalues(options,'ice_cap',{}); 
				
				%recover the earth model:
				slm.earth = getfieldvalue(options,'earth');
			end
		end
		%}}}
		function checkconsistency(slm,solutiontype) % {{{

			%is the coupler turned on? 
			for i=1:length(slm.icecaps),
				if slm.icecaps{i}.transient.iscoupler==0,
					error(sprintf('sealevelmodel checkconsistenty error:  icecap model %s should have the transient coupler option turned on!',slm.icecaps{i}.miscellaneous.name));
				end
			end
				
			if slm.earth.transient.iscoupler==0,
				error('sealevelmodel checkconsistenty error:  earth model should have the transient coupler option turned on!');
			end

			%check that the transition vectors have the right size: 
			for i=1:length(slm.icecaps),
				if slm.icecaps{i}.mesh.numberofvertices ~= length(slm.earth.slr.transitions{i}),
					error('sealevelmodel checkconsistenty issue with size of transition vectors!');
				end
			end


		end
		%}}}
		function slm = setdefaultparameters(slm) % {{{

			%initialize subclasses
			slm.icecaps           = {};
			slm.earth             = {};
			slm.miscellaneous     = miscellaneous();
			slm.settings          = settings();
			slm.private           = private();
			slm.cluster           = generic();
		end
		%}}}
		function disp(self) % {{{
			disp(sprintf('%19s: %-22s -- %s','icecaps'         ,['[' num2str(length(self.icecaps)) 'x1 ' class(self.icecaps) ']'],'ice caps'));
			disp(sprintf('%19s: %-22s -- %s','earth'           ,['[1x1 ' class(self.earth) ']'],'earth'));
			disp(sprintf('%19s: %-22s -- %s','settings'        ,['[1x1 ' class(self.settings) ']'],'settings properties'));
			disp(sprintf('%19s: %-22s -- %s','cluster'         ,['[1x1 ' class(self.cluster) ']'],'cluster parameters (number of cpus...)'));
			disp(sprintf('%19s: %-22s -- %s','miscellaneous'   ,['[1x1 ' class(self.miscellaneous) ']'],'miscellaneous fields'));
		end % }}}
	end
end
