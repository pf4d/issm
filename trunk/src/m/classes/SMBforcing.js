//SMBforcing Class definition
//
//   Usage:
//      SMB=SMBforcing();

function SMBforcing(){
	//methods
	this.setdefaultparameters = function(){// {{{
		this.requested_outputs=['default'];
	} // }}}
	this.disp = function(){ // {{{
		console.log(sprintf('   surface forcings parameters:'));
		fielddisplay(this,'mass_balance','surface mass balance [m/yr ice eq]');
		fielddisplay(this,'requested_outputs','additional outputs requested');
	} // }}}
	this.defaultoutputs = function(){ // {{{
		return '';
	}//}}}
    this.classname = function(){ // {{{
        return "SMBforcing";
    } // }}}
    this.extrude = function(md) {//{{{
        this.mass_balance=project3d(md,'vector',this.mass_balance,'type','node');
        return this;
    }//}}}
    this.initialize = function(md) {// {{{

        if (isNaN(this.mass_balance)){
            this.mass_balance=NewArrayFill(md.mesh.numberofvertices,0);
            console.log('      no smb.mass_balance specified: values set as zero');
        }

    } // }}}
    this.checkconsistency = function(md,solution,analyses) { //{{{

        if(ArrayAnyEqual(ArrayIsMember('MasstransportAnalysis',analyses),1)){
            checkfield(md,'fieldname','smb.mass_balance','timeseries',1,'NaN',1,'Inf',1);
        }
        if(ArrayAnyEqual(ArrayIsMember('BalancethicknessAnalysis',analyses),1)){
            checkfield(md,'fieldname','smb.mass_balance','size',[md.mesh.numberofvertices,1],'NaN',1,'Inf',1);
        }
        checkfield(md,'fieldname','smb.requested_outputs','stringrow',1);

    } // }}}
    this.marshall=function(md,prefix,fid) { //{{{

        var yts=md.constants.yts;

        WriteData(fid,prefix,'name','md.smb.model','data',1,'format','Integer');
        WriteData(fid,prefix,'object',this,'class','smb','fieldname','mass_balance','format','DoubleMat','mattype',1,'scale',1./yts,'timeserieslength',md.mesh.numberofvertices+1,'yts',md.constants.yts);

        //process requested outputs
        var outputs = this.requested_outputs.slice();
        for (var i=0;i<outputs.length;i++){
            if (outputs[i] == 'default') {
                outputs.splice(i,1);
                var newoutputs=this.defaultoutputs(md);
                for (var j=0;j<newoutputs.length;j++) outputs.push(newoutputs[j]);
            }
        }
        WriteData(fid,prefix,'data',outputs,'name','md.smb.requested_outputs','format','StringArray');

    }//}}}
    this.fix=function() { //{{{
    }//}}}
	//properties 
    // {{{
	this.mass_balance = NaN;
	this.requested_outputs      = [];
	this.setdefaultparameters();
    // }}}
}
