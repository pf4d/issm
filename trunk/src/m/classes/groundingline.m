%GROUNDINGLINE class definition
%
%   Usage:
%      groundingline=groundingline();

classdef groundingline
	properties (SetAccess=public) 
		migration    = '';
	end
	methods
		function self = groundingline(varargin) % {{{
			switch nargin
				case 0
					self=setdefaultparameters(self);
				otherwise
					error('constructor not supported');
			end
		end % }}}
		function self = setdefaultparameters(self) % {{{

			%Type of migration
			self.migration='None';

		end % }}}
		function md = checkconsistency(self,md,solution,analyses) % {{{

			md = checkfield(md,'fieldname','groundingline.migration','values',{'None' 'AggressiveMigration' 'SoftMigration' 'SubelementMigration' 'SubelementMigration2' 'Contact' 'None' 'GroundingOnly'});

			if ~strcmp(self.migration,'None'),
				if isnan(md.geometry.bed),
					md = checkmessage(md,['requesting grounding line migration, but bathymetry is absent!']);
				end
				pos=find(md.mask.groundedice_levelset>0. & md.mask.ice_levelset<=0);
				if any(abs(md.geometry.base(pos)-md.geometry.bed(pos))>10^-10),
					md = checkmessage(md,['base not equal to bed on grounded ice!']);
				end
				pos=find(md.mask.groundedice_levelset<=0. & md.mask.ice_levelset<=0);
				if any(md.geometry.bed(pos) - md.geometry.base(pos) > 10^-9),
					md = checkmessage(md,['bed superior to base on floating ice!']);
				end
			end

		end % }}}
		function disp(self) % {{{
			disp(sprintf('   grounding line migration parameters:'));
			fielddisplay(self,'migration','type of grounding line migration: ''SoftMigration'',''AggressiveMigration'',''SubelementMigration'',''SubelementMigration2'' or ''None''');

		end % }}}
		function marshall(self,prefix,md,fid) % {{{
			WriteData(fid,prefix,'data',self.migration,'name','md.groundingline.migration','format','String');
		end % }}}
		function savemodeljs(self,fid,modelname) % {{{
		
			writejsstring(fid,[modelname '.groundingline.migration'],self.migration);

		end % }}}
	end
end
