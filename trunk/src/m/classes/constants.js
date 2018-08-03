//CONSTANTS class definition
//
//   Usage:
//      constants=constants();

function constants() {
	//methods 
		this.setdefaultparameters = function (){ //{{{

			//acceleration due to gravity (m/s^2)
			this.g=9.81;

			//Earth's rotation speed 
			this.omega = 7.292*1e-5;

			//converstion from year to seconds
			this.yts=365.*24.*3600.;

			//the reference temperature for enthalpy model (cf Aschwanden)
			this.referencetemperature=223.15;
		}// }}}
		this.disp = function () { //{{{
			console.log(sprintf("   Constants parameters:")); 
			
			fielddisplay(this,'g','gravitational acceleration [m/s^2]');
			fielddisplay(this,'omega','angular velocity of Earth [rad/s]');
			fielddisplay(this,'yts','number of seconds in a year [s/yr]');
			fielddisplay(this,'referencetemperature','reference temperature used in the enthalpy model [K]');

		} //}}}
		this.classname = function () { //{{{
			return "constants";

		} //}}}
		this.checkconsistency = function(md,solution,analyses) {//% {{{

			checkfield(md,'fieldname','constants.g','>=',0,'size',[1,1]); //We allow 0 for validation tests
			checkfield(md,'fieldname','constants.omega','>=',0,'size',[1,1]);
			checkfield(md,'fieldname','constants.yts','>',0,'size',[1,1]);
			checkfield(md,'fieldname','constants.referencetemperature','size',[1,1]);

		} // }}}
		this.marshall=function(md,prefix,fid) { //{{{
			WriteData(fid,prefix,'object',this,'fieldname','g','format','Double');
			WriteData(fid,prefix,'object',this,'fieldname','yts','format','Double');
			WriteData(fid,prefix,'object',this,'fieldname','referencetemperature','format','Double');
		}//}}}
		this.fix=function() { //{{{
		}//}}}
	//properties 
	// {{{
		this.g                    = 0.;
		this.omega                = 0.;
		this.yts                  = 0.;
		this.referencetemperature = 0.;
		this.setdefaultparameters();
		//}}}
}
