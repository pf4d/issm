from fielddisplay import fielddisplay
from EnumDefinitions import *
from checkfield import checkfield
from WriteData import WriteData
import numpy

class mismipbasalforcings(object):
    """
    MISMIP Basal Forcings class definition

        Usage:
	    mismipbasalforcings=mismipbasalforcings()
    """

    def __init__(self,md): # {{{

        self.groundedice_melting_rate = float('NaN')
        self.meltrate_factor = float('NaN')
        self.threshold_thickness = float('NaN')
        self.upperdepth_melt = float('NaN')
        self.geothermalflux = float('NaN')

        if numpy.all(numpy.isnan(self.groundedice_melting_rate)):
            self.groundedice_melting_rate=numpy.zeros(md.mesh.numberofvertices)
            print(' no basalforcings.groundedice_melting_rate specified: values set as zero')

	self.setdefaultparameters()

    #}}}
    def __repr__(self): # {{{
        string=" MISMIP+ basal melt parameterization\n"
        string="%s\n%s"%(string,fielddisplay(self,"groundedice_melting_rate","basal melting rate (positive if melting) [m/yr]"))
        string="%s\n%s"%(string,fielddisplay(self,"meltrate_factor","Melt-rate rate factor [1/yr] (sign is opposite to MISMIP+ benchmark to remain consistent with ISSM convention of positive values for melting)"))
        string="%s\n%s"%(string,fielddisplay(self,"threshold_thickness","Threshold thickness for saturation of basal melting [m]"))
        string="%s\n%s"%(string,fielddisplay(self,"upperdepth_melt","Depth above which melt rate is zero [m]"))
        string="%s\n%s"%(string,fielddisplay(self,"geothermalflux","Geothermal heat flux [W/m^2]"))

	return string
    #}}}
    def extrude(self,md): # {{{
	self.coefficient=project3d(md,'vector',self.coefficient,'type','node','layer',1)
	self.p=project3d(md,'vector',self.p,'type','element')
	self.q=project3d(md,'vector',self.q,'type','element')
	return self
    #}}}
    def setdefaultparameters(self): # {{{

        # default values for melting parameterization
        self.meltrate_factor = 0.2
        self.threshold_thickness = 75.
        self.upperdepth_melt = -100.

	return self
    #}}}
    def checkconsistency(self,md,solution,analyses):    # {{{

	#Early return
        if MasstransportAnalysisEnum() in analyses and not (solution==TransientSolutionEnum() and md.transient.ismasstransport==0):

	    md = checkfield(md,'fieldname','basalforcings.groundedice_melting_rate','NaN',1,'Inf',1,'timeseries',1)
	    md = checkfield(md,'fieldname','basalforcings.meltrate_factor','>=',0,'numel',[1])
	    md = checkfield(md,'fieldname','basalforcings.threshold_thickness','>=',0,'numel',[1])
	    md = checkfield(md,'fieldname','basalforcings.upperdepth_melt','<=',0,'numel',[1])

        if BalancethicknessAnalysisEnum() in analyses:

	    md = checkfield(md,'fieldname','basalforcings.groundedice_melting_rate','NaN',1,'Inf',1,'size',[md.mesh.numberofvertices])
	    md = checkfield(md,'fieldname','basalforcings.meltrate_factor','>=',0,'numel',[1])
	    md = checkfield(md,'fieldname','basalforcings.threshold_thickness','>=',0,'numel',[1])
	    md = checkfield(md,'fieldname','basalforcings.upperdepth_melt','<=',0,'numel',[1])

        if ThermalAnalysisEnum() in analyses and not (solution==TransientSolutionEnum() and md.transient.isthermal==0):

	    md = checkfield(md,'fieldname','basalforcings.groundedice_melting_rate','NaN',1,'Inf',1,'timeseries',1)
	    md = checkfield(md,'fieldname','basalforcings.meltrate_factor','>=',0,'numel',[1])
	    md = checkfield(md,'fieldname','basalforcings.threshold_thickness','>=',0,'numel',[1])
	    md = checkfield(md,'fieldname','basalforcings.upperdepth_melt','<=',0,'numel',[1])
	    md = checkfield(md,'fieldname','basalforcings.geothermalflux','NaN',1,'Inf',1,'timeseries',1,'>=',0)
	return md
    # }}}
    def marshall(self,md,fid):    # {{{

        yts=md.constants.yts
        if yts!=365.2422*24.*3600.:
            print('WARNING: value of yts for MISMIP+ runs different from ISSM default!')

        floatingice_melting_rate = numpy.zeros((md.mesh.numberofvertices,1))
        floatingice_melting_rate = md.basalforcings.meltrate_factor*numpy.tanh((md.geometry.base-md.geometry.bed)/md.basalforcings.threshold_thickness)*numpy.amax(md.basalforcings.upperdepth_melt-md.geometry.base,0)

	WriteData(fid,'enum',BasalforcingsEnum(),'data',MismipFloatingMeltRateEnum(),'format','Integer')
	WriteData(fid,'data',floatingice_melting_rate,'format','DoubleMat','enum',BasalforcingsFloatingiceMeltingRateEnum(),'mattype',1,'scale',1./yts,'timeserieslength',md.mesh.numberofvertices+1)
	WriteData(fid,'object',self,'fieldname','groundedice_melting_rate','format','DoubleMat','enum',BasalforcingsGroundediceMeltingRateEnum(),'mattype',1,'scale',1./yts,'timeserieslength',md.mesh.numberofvertices+1)
	WriteData(fid,'object',self,'fieldname','geothermalflux','enum',BasalforcingsGeothermalfluxEnum(),'format','DoubleMat','mattype',1,'timeserieslength',md.mesh.numberofvertices+1)
	WriteData(fid,'object',self,'fieldname','meltrate_factor','format','Double','enum',BasalforcingsMeltrateFactorEnum(),'scale',1./yts)
	WriteData(fid,'object',self,'fieldname','threshold_thickness','format','Double','enum',BasalforcingsThresholdThicknessEnum())
	WriteData(fid,'object',self,'fieldname','upperdepth_melt','format','Double','enum',BasalforcingsUpperdepthMeltEnum())

    # }}}
