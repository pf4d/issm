%TRANSIENT class definition
%
%   Usage:
%      transient=transient();

classdef transient
	properties (SetAccess=public) 
		issmb             = 0;
		ismasstransport   = 0;
		isstressbalance   = 0;
		isthermal         = 0;
		isgroundingline   = 0;
		isgia             = 0;
		isesa             = 0;
		isdamageevolution = 0;
		ismovingfront     = 0;
		ishydrology       = 0;
		isslr             = 0;
		iscoupler         = 0;
		amr_frequency     = 0;
		isoceancoupling   = 0;
		requested_outputs = {};
	end
	methods
		function self = transient(varargin) % {{{
			switch nargin
				case 0
					self = setdefaultparameters(self);
				otherwise
					error('constructor not supported');
			end
		end % }}}
		function self = deactivateall(self) % {{{

			%full analysis: Stressbalance, Masstransport and Thermal but no groundingline migration for now
			self.issmb           = 0;
			self.ismasstransport = 0;
			self.isstressbalance = 0;
			self.isthermal       = 0;
			self.isgroundingline = 0;
			self.isgia           = 0;
			self.isesa           = 0;
			self.isdamageevolution = 0;
			self.ismovingfront   =0;
			self.ishydrology     = 0;
			self.isslr           = 0;
			self.isoceancoupling = 0;
			self.iscoupler       = 0;
			self.amr_frequency	= 0;

			%default output
			self.requested_outputs={};
		end % }}}
		function self = setdefaultparameters(self) % {{{

			%full analysis: Stressbalance, Masstransport and Thermal but no groundingline migration for now
			self.issmb           = 1;
			self.ismasstransport = 1;
			self.isstressbalance = 1;
			self.isthermal       = 1;
			self.isgroundingline = 0;
			self.isgia           = 0;
			self.isesa           = 0;
			self.isdamageevolution = 0;
			self.ismovingfront   = 0;
			self.ishydrology     = 0;
			self.isslr           = 0;
			self.isoceancoupling = 0;
			self.iscoupler       = 0;
			self.amr_frequency	= 1;

			%default output
			self.requested_outputs={'default'};
		end % }}}
		function list = defaultoutputs(self,md) % {{{
			if(self.issmb)
				list = {'SmbMassBalance'};
			else
				list = {};
			end
		end % }}}
		function md = checkconsistency(self,md,solution,analyses) % {{{

			%Early return
			if ~strcmp(solution,'TransientSolution'), return; end

			md = checkfield(md,'fieldname','transient.issmb','numel',[1],'values',[0 1]);
			md = checkfield(md,'fieldname','transient.ismasstransport','numel',[1],'values',[0 1]);
			md = checkfield(md,'fieldname','transient.isstressbalance','numel',[1],'values',[0 1]);
			md = checkfield(md,'fieldname','transient.isthermal','numel',[1],'values',[0 1]);
			md = checkfield(md,'fieldname','transient.isgroundingline','numel',[1],'values',[0 1]);
			md = checkfield(md,'fieldname','transient.isgia','numel',[1],'values',[0 1]);
			md = checkfield(md,'fieldname','transient.isesa','numel',[1],'values',[0 1]);
			md = checkfield(md,'fieldname','transient.isdamageevolution','numel',[1],'values',[0 1]);
			md = checkfield(md,'fieldname','transient.ismovingfront','numel',[1],'values',[0 1]);
			md = checkfield(md,'fieldname','transient.ishydrology','numel',[1],'values',[0 1]);
			md = checkfield(md,'fieldname','transient.requested_outputs','stringrow',1);
			md = checkfield(md,'fieldname','transient.isslr','numel',[1],'values',[0 1]);
			md = checkfield(md,'fieldname','transient.isoceancoupling','numel',[1],'values',[0 1]);
			md = checkfield(md,'fieldname','transient.iscoupler','numel',[1],'values',[0 1]);
			md = checkfield(md,'fieldname','transient.amr_frequency','numel',[1],'>=',0,'NaN',1,'Inf',1);

			if (~strcmp(solution,'TransientSolution') & md.transient.iscoupling==1), 
				md = checkmessage(md,['Coupling with ocean model can only be performed for transient simulations!']);
			end
		end % }}}
		function disp(self) % {{{
			disp(sprintf('   transient solution parameters:'));

			fielddisplay(self,'issmb','indicates whether a surface mass balance solution is used in the transient');
			fielddisplay(self,'ismasstransport','indicates whether a masstransport solution is used in the transient');
			fielddisplay(self,'isstressbalance','indicates whether a stressbalance solution is used in the transient');
			fielddisplay(self,'isthermal','indicates whether a thermal solution is used in the transient');
			fielddisplay(self,'isgroundingline','indicates whether a groundingline migration is used in the transient');
			fielddisplay(self,'isgia','indicates whether a postglacial rebound model is used in the transient');
			fielddisplay(self,'isesa','indicates whether an elastic adjustment model is used in the transient');
			fielddisplay(self,'isdamageevolution','indicates whether damage evolution is used in the transient');
			fielddisplay(self,'ismovingfront','indicates whether a moving front capability is used in the transient');
			fielddisplay(self,'ishydrology','indicates whether an hydrology model is used');
			fielddisplay(self,'isslr','indicates whether a sea-level rise solution is used in the transient');
			fielddisplay(self,'isoceancoupling','indicates whether a coupling with an ocean model is used in the transient');
			fielddisplay(self,'iscoupler','indicates whether different models are being run with need for coupling');
			fielddisplay(self,'amr_frequency','frequency at which mesh is refined in simulations with multiple time_steps');
			fielddisplay(self,'requested_outputs','list of additional outputs requested');

		end % }}}
		function marshall(self,prefix,md,fid) % {{{
			WriteData(fid,prefix,'object',self,'fieldname','issmb','format','Boolean');
			WriteData(fid,prefix,'object',self,'fieldname','ismasstransport','format','Boolean');
			WriteData(fid,prefix,'object',self,'fieldname','isstressbalance','format','Boolean');
			WriteData(fid,prefix,'object',self,'fieldname','isthermal','format','Boolean');
			WriteData(fid,prefix,'object',self,'fieldname','isgroundingline','format','Boolean');
			WriteData(fid,prefix,'object',self,'fieldname','isgia','format','Boolean');
			WriteData(fid,prefix,'object',self,'fieldname','isesa','format','Boolean');
			WriteData(fid,prefix,'object',self,'fieldname','isdamageevolution','format','Boolean');
			WriteData(fid,prefix,'object',self,'fieldname','ishydrology','format','Boolean');
			WriteData(fid,prefix,'object',self,'fieldname','ismovingfront','format','Boolean');
			WriteData(fid,prefix,'object',self,'fieldname','isslr','format','Boolean');
			WriteData(fid,prefix,'object',self,'fieldname','isoceancoupling','format','Boolean');
			WriteData(fid,prefix,'object',self,'fieldname','iscoupler','format','Boolean');
			WriteData(fid,prefix,'object',self,'fieldname','amr_frequency','format','Integer');

			%process requested outputs
			outputs = self.requested_outputs;
			pos  = find(ismember(outputs,'default'));
			if ~isempty(pos),
				outputs(pos) = [];                         %remove 'default' from outputs
				outputs      = [outputs defaultoutputs(self,md)]; %add defaults
			end
			WriteData(fid,prefix,'data',outputs,'name','md.transient.requested_outputs','format','StringArray');
		end % }}}
		function savemodeljs(self,fid,modelname) % {{{
		
			writejsdouble(fid,[modelname '.trans.issmb'],self.issmb);
			writejsdouble(fid,[modelname '.trans.ismasstransport'],self.ismasstransport);
			writejsdouble(fid,[modelname '.trans.isstressbalance'],self.isstressbalance);
			writejsdouble(fid,[modelname '.trans.isthermal'],self.isthermal);
			writejsdouble(fid,[modelname '.trans.isgroundingline'],self.isgroundingline);
			writejsdouble(fid,[modelname '.trans.isgia'],self.isgia);
			writejsdouble(fid,[modelname '.trans.isesa'],self.isesa);
			writejsdouble(fid,[modelname '.trans.isdamageevolution'],self.isdamageevolution);
			writejsdouble(fid,[modelname '.trans.ismovingfront'],self.ismovingfront);
			writejsdouble(fid,[modelname '.trans.ishydrology'],self.ishydrology);
			writejsdouble(fid,[modelname '.trans.isslr'],self.isslr);
			writejsdouble(fid,[modelname '.trans.isoceancoupling'],self.isoceancoupling);
			writejsdouble(fid,[modelname '.trans.iscoupler'],self.iscoupler);
			writejsdouble(fid,[modelname '.trans.amr_frequency'],self.amr_frequency);
			writejscellstring(fid,[modelname '.trans.requested_outputs'],self.requested_outputs);

		end % }}}
	end
end
