sed -i 's/(md,'fieldname','friction.coefficient','timeseries',1,'NaN',1)/(md,fieldname='friction.coefficient',timeseries=1,NaN=1)/g' ./classes/frictioncoulomb.py	
sed -i 's/(md,'fieldname','friction.coefficientcoulomb','timeseries',1,'NaN',1)/(md,fieldname='friction.coefficientcoulomb',timeseries=1,NaN=1)/g' ./classes/frictioncoulomb.py
sed -i 's/(md,'fieldname','friction.q','NaN',1,'size',[md.mesh.numberofelements])/(md,fieldname='friction.q',NaN=1,size=[md.mesh.numberofelements])/g' ./classes/frictioncoulomb.py	
sed -i 's/(md,'fieldname','friction.p','NaN',1,'size',[md.mesh.numberofelements])/(md,fieldname='friction.p',NaN=1,size=[md.mesh.numberofelements])/g' ./classes/frictioncoulomb.py
sed -i 's/(md,'fieldname','autodiff.obufsize','>=',524288)/(md,fieldname='autodiff.obufsize',ge=524288)/g' ./classes/autodiff.py
sed -i 's/(md,'fieldname','autodiff.lbufsize','>=',524288)/(md,fieldname='autodiff.lbufsize',ge=524288)/g' ./classes/autodiff.py	
sed -i 's/(md,'fieldname','autodiff.cbufsize','>=',524288)/(md,fieldname='autodiff.cbufsize',ge=524288)/g' ./classes/autodiff.py		
sed -i 's/(md,'fieldname','autodiff.tbufsize','>=',524288)/(md,fieldname='autodiff.tbufsize',ge=524288)/g' ./classes/autodiff.py		
sed -i 's/(md,'fieldname','autodiff.gcTriggerRatio','>=',2.0)/(md,fieldname='autodiff.gcTriggerRatio',ge=2.0)/g' ./classes/autodiff.py
sed -i 's/(md,'fieldname','autodiff.gcTriggerMaxSize','>=',2000000)/(md,fieldname='autodiff.gcTriggerMaxSize',ge=2000000)/g' ./classes/autodiff.py
sed -i 's/(md,'fieldname','autodiff.driver','values',['fos_forward','fov_forward','fov_forward_all','fos_reverse','fov_reverse','fov_reverse_all'])/(md,fieldname='autodiff.driver',values=['fos_forward','fov_forward','fov_forward_all','fos_reverse','fov_reverse','fov_reverse_all'])/g' ./classes/autodiff.py
sed -i 's/(md,'fieldname','materials.rho_ice','>',0)/(md,fieldname='materials.rho_ice',gt=0)/g' ./classes/matice.py
sed -i 's/(md,'fieldname','materials.rho_water','>',0)/(md,fieldname='materials.rho_water',gt=0)/g' ./classes/matice.py
sed -i 's/(md,'fieldname','materials.rho_freshwater','>',0)/(md,fieldname='materials.rho_freshwater',gt=0)/g' ./classes/matice.py
sed -i 's/(md,'fieldname','materials.mu_water','>',0)/(md,fieldname='materials.mu_water',gt=0)/g' ./classes/matice.py		
sed -i 's/(md,'fieldname','materials.rheology_B','>',0,'timeseries',1,'NaN',1)/(md,fieldname='materials.rheology_B',gt=0,timeseries=1,NaN=1)/g' ./classes/matice.py
sed -i 's/(md,'fieldname','materials.rheology_n','>',0,'size',[md.mesh.numberofelements])/(md,fieldname='materials.rheology_n',gt=0,size=[md.mesh.numberofelements])/g' ./classes/matice.py
sed -i 's/(md,'fieldname','materials.rheology_law','values',['None','Cuffey','Paterson','Arrhenius','LliboutryDuval'])/(md,fieldname='materials.rheology_law',values=['None','Cuffey','Paterson','Arrhenius','LliboutryDuval'])/g' ./classes/matice.py		
sed -i 's/(md,'fieldname','materials.lithosphere_shear_modulus','>',0,'numel',[1]);/(md,fieldname='materials.lithosphere_shear_modulus',gt=0,numel=[1]);/g' ./classes/matice.py
sed -i 's/(md,'fieldname','materials.lithosphere_density','>',0,'numel',[1]);/(md,fieldname='materials.lithosphere_density',gt=0,numel=[1]);/g' ./classes/matice.py
sed -i 's/(md,'fieldname','materials.mantle_shear_modulus','>',0,'numel',[1]);/(md,fieldname='materials.mantle_shear_modulus',gt=0,numel=[1]);/g' ./classes/matice.py
sed -i 's/(md,'fieldname','materials.mantle_density','>',0,'numel',[1]);/(md,fieldname='materials.mantle_density',gt=0,numel=[1]);/g' ./classes/matice.py:		
sed -i 's/(md,'fieldname','smb.desfac','<=',1,'numel',[1])/(md,fieldname='smb.desfac',le=1,numel=[1])/g' ./classes/SMBd18opdd.py:			
sed -i 's/(md,'fieldname','smb.s0p','>=',0,'NaN',1,'size',[md.mesh.numberofvertices,1])/(md,fieldname='smb.s0p',ge=0,NaN=1,size=[md.mesh.numberofvertices,1])/g' ./classes/SMBd18opdd.py:			
sed -i 's/(md,'fieldname','smb.s0t','>=',0,'NaN',1,'size',[md.mesh.numberofvertices,1])/(md,fieldname='smb.s0t',ge=0,NaN=1,size=[md.mesh.numberofvertices,1])/g' ./classes/SMBd18opdd.py:			
sed -i 's/(md,'fieldname','smb.rlaps','>=',0,'numel',[1])/(md,fieldname='smb.rlaps',ge=0,numel=[1])/g' ./classes/SMBd18opdd.py:			
sed -i 's/(md,'fieldname','smb.rlapslgm','>=',0,'numel',[1])/(md,fieldname='smb.rlapslgm',ge=0,numel=[1])/g' ./classes/SMBd18opdd.py:			
sed -i 's/(md,'fieldname','smb.temperatures_presentday','size',[md.mesh.numberofvertices+1,12],'NaN',1,'timeseries',1)/(md,fieldname='smb.temperatures_presentday',size=[md.mesh.numberofvertices+1,12],NaN=1,timeseries=1)/g' ./classes/SMBd18opdd.py:				
sed -i 's/(md,'fieldname','smb.precipitations_presentday','size',[md.mesh.numberofvertices+1,12],'NaN',1,'timeseries',1)/(md,fieldname='smb.precipitations_presentday',size=[md.mesh.numberofvertices+1,12],NaN=1,timeseries=1)/g' ./classes/SMBd18opdd.py:				
sed -i 's/(md,'fieldname','smb.delta18o',NaN=1,'size',[2,numpy.nan],'singletimeseries',1)/(md,fieldname='smb.delta18o',NaN=1,size=[2,numpy.nan],'singletimesed -i 's/(md,'fieldname','smb.dpermil','>=',0,'numel',[1])/(md,fieldname='smb.dpermil',ge=0,numel=[1])/g' ./classes/SMBd18opdd.py:				
sed -i 's/(md,'fieldname','masstransport.requested_outputs','stringrow',1)/(md,fieldname='masstransport.requested_outputs','stringrow',1)/g' ./classes/SMBd18opdd.py:		
sed -i 's/(md,'fieldname','inversion.iscontrol','values',[0,1])/(md,fieldname='inversion.iscontrol',values=[0,1])/g' ./classes/m1qn3inversion.py:		
sed -i 's/(md,'fieldname','inversion.incomplete_adjoint','values',[0,1])/(md,fieldname='inversion.incomplete_adjoint',values=[0,1])/g' ./classes/m1qn3inversion.py:		
sed -i 's/(md,'fieldname','inversion.control_parameters','cell',1,'values',supportedcontrols())/(md,fieldname='inversion.control_parameters','cell',1,values=supportedcontrols())/g' ./classes/m1qn3inversion.py:		
sed -i 's/(md,'fieldname','inversion.control_scaling_factors','size',[num_controls],'>',0,'NaN',1)/(md,fieldname='inversion.control_scaling_factors',size=[num_controls],gt=0,NaN=1)/g' ./classes/m1qn3inversion.py:		
sed -i 's/(md,'fieldname','inversion.maxsteps','numel',[1],'>=',0)/(md,fieldname='inversion.maxsteps',numel=[1],ge=0)/g' ./classes/m1qn3inversion.py:		
sed -i 's/(md,'fieldname','inversion.maxiter','numel',[1],'>=',0)/(md,fieldname='inversion.maxiter',numel=[1],ge=0)/g' ./classes/m1qn3inversion.py:		
sed -i 's/(md,'fieldname','inversion.dxmin','numel',[1],'>',0.)/(md,fieldname='inversion.dxmin',numel=[1],gt=0.)/g' ./classes/m1qn3inversion.py:		
sed -i 's/(md,'fieldname','inversion.gttol','numel',[1],'>',0.)/(md,fieldname='inversion.gttol',numel=[1],gt=0.)/g' ./classes/m1qn3inversion.py:		
sed -i 's/(md,'fieldname','inversion.cost_functions','size',[num_costfunc],'values',supportedcostfunctions())/(md,fieldname='inversion.cost_functions',size=[num_costfunc],values=supportedcostfunctions())/g' ./classes/m1qn3inversion.py:		
sed -i 's/(md,'fieldname','inversion.cost_functions_coefficients','size',[md.mesh.numberofvertices,num_costfunc],'>=',0)/(md,fieldname='inversion.cost_functions_coefficients',size=[md.mesh.numberofvertices,num_costfunc],ge=0)/g' ./classes/m1qn3inversion.py:		
sed -i 's/(md,'fieldname','inversion.min_parameters','size',[md.mesh.numberofvertices,num_controls])/(md,fieldname='inversion.min_parameters',size=[md.mesh.numberofvertices,num_controls])/g' ./classes/m1qn3inversion.py:		
sed -i 's/(md,'fieldname','inversion.max_parameters','size',[md.mesh.numberofvertices,num_controls])/(md,fieldname='inversion.max_parameters',size=[md.mesh.numberofvertices,num_controls])/g' ./classes/m1qn3inversion.py:		
sed -i 's/(md,'fieldname','inversion.thickness_obs','size',[md.mesh.numberofvertices],'NaN',1)/(md,fieldname='inversion.thickness_obs',size=[md.mesh.numberofvertices],NaN=1)/g' ./classes/m1qn3inversion.py:			
sed -i 's/(md,'fieldname','inversion.vx_obs','size',[md.mesh.numberofvertices],'NaN',1)/(md,fieldname='inversion.vx_obs',size=[md.mesh.numberofvertices],NaN=1)/g' ./classes/m1qn3inversion.py:			
sed -i 's/(md,'fieldname','inversion.vy_obs','size',[md.mesh.numberofvertices],'NaN',1)/(md,fieldname='inversion.vy_obs',size=[md.mesh.numberofvertices],NaN=1)/g' ./classes/m1qn3inversion.py:			
sed -i 's///g' ./classes/masstransport.py:from checkfield import checkfield
sed -i 's/(md,'fieldname','masstransport.spcthickness','timeseries',1)/(md,fieldname='masstransport.spcthickness',timeseries=1)/g' ./classes/masstransport.py:		
sed -i 's/(md,'fieldname','masstransport.isfreesurface','values',[0,1])/(md,fieldname='masstransport.isfreesurface',values=[0,1])/g' ./classes/masstransport.py:		
sed -i 's/(md,'fieldname','masstransport.hydrostatic_adjustment','values',['Absolute','Incremental'])/(md,fieldname='masstransport.hydrostatic_adjustment',values=['Absolute','Incremental'])/g' ./classes/masstransport.py:		
sed -i 's/(md,'fieldname','masstransport.stabilization','values',[0,1,2,3,4])/(md,fieldname='masstransport.stabilization',values=[0,1,2,3,4])/g' ./classes/masstransport.py:		
sed -i 's/(md,'fieldname','masstransport.min_thickness','>',0)/(md,fieldname='masstransport.min_thickness',gt=0)/g' ./classes/masstransport.py:		
sed -i 's/(md,'fieldname','masstransport.requested_outputs','stringrow',1)/(md,fieldname='masstransport.requested_outputs','stringrow',1)/g' ./classes/masstransport.py:		
sed -i 's///g' ./classes/mismipbasalforcings.py:from checkfield import checkfield
sed -i 's/(md,'fieldname','basalforcings.groundedice_melting_rate','NaN',1,'timeseries',1)/(md,fieldname='basalforcings.groundedice_melting_rate',NaN=1,timeseries=1)/g' ./classes/mismipbasalforcings.py:	    
sed -i 's/(md,'fieldname','basalforcings.meltrate_factor','>=',0,'numel',[1])/(md,fieldname='basalforcings.meltrate_factor',ge=0,numel=[1])/g' ./classes/mismipbasalforcings.py:	    
sed -i 's/(md,'fieldname','basalforcings.threshold_thickness','>=',0,'numel',[1])/(md,fieldname='basalforcings.threshold_thickness',ge=0,numel=[1])/g' ./classes/mismipbasalforcings.py:	    
sed -i 's/(md,'fieldname','basalforcings.upperdepth_melt','<=',0,'numel',[1])/(md,fieldname='basalforcings.upperdepth_melt',le=0,numel=[1])/g' ./classes/mismipbasalforcings.py:	    
sed -i 's/(md,'fieldname','basalforcings.groundedice_melting_rate','NaN',1,'size',[md.mesh.numberofvertices])/(md,fieldname='basalforcings.groundedice_melting_rate',NaN=1,size=[md.mesh.numberofvertices])/g' ./classes/mismipbasalforcings.py:	    
sed -i 's/(md,'fieldname','basalforcings.meltrate_factor','>=',0,'numel',[1])/(md,fieldname='basalforcings.meltrate_factor',ge=0,numel=[1])/g' ./classes/mismipbasalforcings.py:	    
sed -i 's/(md,'fieldname','basalforcings.threshold_thickness','>=',0,'numel',[1])/(md,fieldname='basalforcings.threshold_thickness',ge=0,numel=[1])/g' ./classes/mismipbasalforcings.py:	    
sed -i 's/(md,'fieldname','basalforcings.upperdepth_melt','<=',0,'numel',[1])/(md,fieldname='basalforcings.upperdepth_melt',le=0,numel=[1])/g' ./classes/mismipbasalforcings.py:	    
sed -i 's/(md,'fieldname','basalforcings.groundedice_melting_rate','NaN',1,'timeseries',1)/(md,fieldname='basalforcings.groundedice_melting_rate',NaN=1,timeseries=1)/g' ./classes/mismipbasalforcings.py:	    
sed -i 's/(md,'fieldname','basalforcings.meltrate_factor','>=',0,'numel',[1])/(md,fieldname='basalforcings.meltrate_factor',ge=0,numel=[1])/g' ./classes/mismipbasalforcings.py:	    
sed -i 's/(md,'fieldname','basalforcings.threshold_thickness','>=',0,'numel',[1])/(md,fieldname='basalforcings.threshold_thickness',ge=0,numel=[1])/g' ./classes/mismipbasalforcings.py:	    
sed -i 's/(md,'fieldname','basalforcings.upperdepth_melt','<=',0,'numel',[1])/(md,fieldname='basalforcings.upperdepth_melt',le=0,numel=[1])/g' ./classes/mismipbasalforcings.py:	    
sed -i 's/(md,'fieldname','basalforcings.geothermalflux','NaN',1,'timeseries',1,'>=',0)/(md,fieldname='basalforcings.geothermalflux',NaN=1,timeseries=1,ge=0)/g' ./classes/mismipbasalforcings.py:	    
sed -i 's///g' ./classes/timestepping.py:from checkfield import checkfield
sed -i 's/(md,'fieldname','timestepping.start_time','numel',[1],'NaN',1)/(md,fieldname='timestepping.start_time',numel=[1],NaN=1)/g' ./classes/timestepping.py:		
sed -i 's/(md,'fieldname','timestepping.final_time','numel',[1],'NaN',1)/(md,fieldname='timestepping.final_time',numel=[1],NaN=1)/g' ./classes/timestepping.py:		
sed -i 's/(md,'fieldname','timestepping.time_step','numel',[1],'>=',0,'NaN',1)/(md,fieldname='timestepping.time_step',numel=[1],ge=0,NaN=1)/g' ./classes/timestepping.py:		
sed -i 's/(md,'fieldname','timestepping.time_adapt','numel',[1],'values',[0,1])/(md,fieldname='timestepping.time_adapt',numel=[1],values=[0,1])/g' ./classes/timestepping.py:		
sed -i 's/(md,'fieldname','timestepping.cfl_coefficient','numel',[1],'>',0,'<=',1)/(md,fieldname='timestepping.cfl_coefficient',numel=[1],gt=0,le=1)/g' ./classes/timestepping.py:		
sed -i 's/(md,'fieldname','timestepping.interp_forcings','numel',[1],'values',[0,1])/(md,fieldname='timestepping.interp_forcings',numel=[1],values=[0,1])/g' ./classes/timestepping.py:		
sed -i 's///g' ./classes/calving.py:from checkfield import checkfield
sed -i 's/(md,'fieldname','calving.spclevelset','timeseries',1)/(md,fieldname='calving.spclevelset',timeseries=1)/g' ./classes/calving.py:		
sed -i 's/(md,'fieldname','calving.stabilization','values',[0,1,2]);/(md,fieldname='calving.stabilization',values=[0,1,2]);/g' ./classes/calving.py:		
sed -i 's/(md,'fieldname','calving.calvingrate','>=',0,'timeseries',1,'NaN',1);/(md,fieldname='calving.calvingrate',ge=0,timeseries=1,NaN=1);/g' ./classes/calving.py:		
sed -i 's/(md,'fieldname','calving.meltingrate','>=',0,'timeseries',1,'NaN',1);/(md,fieldname='calving.meltingrate',ge=0,timeseries=1,NaN=1);/g' ./classes/calving.py:		
sed -i 's///g' ./classes/SMBcomponents.py:from checkfield import *
sed -i 's/(md,'fieldname','smb.accumulation','timeseries',1,'NaN',1)/(md,fieldname='smb.accumulation',timeseries=1,NaN=1)/g' ./classes/SMBcomponents.py:			
sed -i 's/(md,'fieldname','smb.accumulation','size',[md.mesh.numberofvertices],'NaN',1)/(md,fieldname='smb.accumulation',size=[md.mesh.numberofvertices],NaN=1)/g' ./classes/SMBcomponents.py:			
sed -i 's/(md,'fieldname','smb.runoff','timeseries',1,'NaN',1)/(md,fieldname='smb.runoff',timeseries=1,NaN=1)/g' ./classes/SMBcomponents.py:			
sed -i 's/(md,'fieldname','smb.runoff','size',[md.mesh.numberofvertices],'NaN',1)/(md,fieldname='smb.runoff',size=[md.mesh.numberofvertices],NaN=1)/g' ./classes/SMBcomponents.py:			
sed -i 's/(md,'fieldname','smb.evaporation','timeseries',1,'NaN',1)/(md,fieldname='smb.evaporation',timeseries=1,NaN=1)/g' ./classes/SMBcomponents.py:			
sed -i 's/(md,'fieldname','smb.evaporation','size',[md.mesh.numberofvertices],'NaN',1)/(md,fieldname='smb.evaporation',size=[md.mesh.numberofvertices],NaN=1)/g' ./classes/SMBcomponents.py:			
sed -i 's/(md,'fieldname','masstransport.requested_outputs','stringrow',1)/(md,fieldname='masstransport.requested_outputs','stringrow',1)/g' ./classes/SMBcomponents.py:		
sed -i 's///g' ./classes/flaim.py:from checkfield import checkfield
sed -i 's/(md,'fieldname','flaim.tracks','file',1)/(md,fieldname='flaim.tracks','file',1)/g' ./classes/flaim.py:		
sed -i 's/(md,'fieldname','flaim.targets','file',1)/(md,fieldname='flaim.targets','file',1)/g' ./classes/flaim.py:			
sed -i 's/(md,'fieldname','flaim.criterion',numel=[md.mesh.numberofvertices,md.mesh.numberofelements])/(md,fieldname='flaim.criterion','numel',[md.mesh.numberofvertices,md.mesh.numberofelements])/g' ./classes/flaim.py:			
sed -i 's/(md,'fieldname',"autodiff.independents[%d].fov_forward_indices" % i,'>=',1,'<=',self.nods,'size',[float('NaN'),1])/(md,fieldname="autodiff.independents[%d].fov_forward_indices" % i,ge=1,le=self.nods,size=[float('NaN'),1])/g' ./classes/independent.py:			
sed -i 's///g' ./classes/calvinglevermann.py:from checkfield import checkfield
sed -i 's/(md,'fieldname','calving.spclevelset','timeseries',1)/(md,fieldname='calving.spclevelset',timeseries=1)/g' ./classes/calvinglevermann.py:		
sed -i 's/(md,'fieldname','calving.stabilization','values',[0,1,2]);/(md,fieldname='calving.stabilization',values=[0,1,2]);/g' ./classes/calvinglevermann.py:		
sed -i 's/(md,'fieldname','calving.coeff','size',[md.mesh.numberofvertices],'>',0)/(md,fieldname='calving.coeff',size=[md.mesh.numberofvertices],gt=0)/g' ./classes/calvinglevermann.py:		
sed -i 's/(md,'fieldname','calving.meltingrate','NaN',1,'size',[md.mesh.numberofvertices],'>=',0)/(md,fieldname='calving.meltingrate',NaN=1,size=[md.mesh.numberofvertices],ge=0)/g' ./classes/calvinglevermann.py:		
sed -i 's///g' ./classes/steadystate.py:from checkfield import checkfield
sed -i 's/(md,'fieldname','steadystate.requested_outputs','stringrow',1)/(md,fieldname='steadystate.requested_outputs','stringrow',1)/g' ./classes/steadystate.py:		
sed -i 's///g' ./classes/SMBmeltcomponents.py:from checkfield import *
sed -i 's/(md,'fieldname','smb.accumulation','timeseries',1,'NaN',1)/(md,fieldname='smb.accumulation',timeseries=1,NaN=1)/g' ./classes/SMBmeltcomponents.py:			
sed -i 's/(md,'fieldname','smb.accumulation','size',[md.mesh.numberofvertices],'NaN',1)/(md,fieldname='smb.accumulation',size=[md.mesh.numberofvertices],NaN=1)/g' ./classes/SMBmeltcomponents.py:			
sed -i 's/(md,'fieldname','smb.melt','timeseries',1,'NaN',1)/(md,fieldname='smb.melt',timeseries=1,NaN=1)/g' ./classes/SMBmeltcomponents.py:			
sed -i 's/(md,'fieldname','smb.melt','size',[md.mesh.numberofvertices],'NaN',1)/(md,fieldname='smb.melt',size=[md.mesh.numberofvertices],NaN=1)/g' ./classes/SMBmeltcomponents.py:			
sed -i 's/(md,'fieldname','smb.refreeze','timeseries',1,'NaN',1)/(md,fieldname='smb.refreeze',timeseries=1,NaN=1)/g' ./classes/SMBmeltcomponents.py:			
sed -i 's/(md,'fieldname','smb.refreeze','size',[md.mesh.numberofvertices],'NaN',1)/(md,fieldname='smb.refreeze',size=[md.mesh.numberofvertices],NaN=1)/g' ./classes/SMBmeltcomponents.py:			
sed -i 's/(md,'fieldname','smb.evaporation','timeseries',1,'NaN',1)/(md,fieldname='smb.evaporation',timeseries=1,NaN=1)/g' ./classes/SMBmeltcomponents.py:			
sed -i 's/(md,'fieldname','smb.evaporation','size',[md.mesh.numberofvertices],'NaN',1)/(md,fieldname='smb.evaporation',size=[md.mesh.numberofvertices],NaN=1)/g' ./classes/SMBmeltcomponents.py:			
sed -i 's/(md,'fieldname','masstransport.requested_outputs','stringrow',1)/(md,fieldname='masstransport.requested_outputs','stringrow',1)/g' ./classes/SMBmeltcomponents.py:		
sed -i 's///g' ./classes/matdamageice.py:from checkfield import checkfield
sed -i 's/(md,'fieldname','materials.rho_ice','>',0)/(md,fieldname='materials.rho_ice',gt=0)/g' ./classes/matdamageice.py:		
sed -i 's/(md,'fieldname','materials.rho_water','>',0)/(md,fieldname='materials.rho_water',gt=0)/g' ./classes/matdamageice.py:		
sed -i 's/(md,'fieldname','materials.rho_freshwater','>',0)/(md,fieldname='materials.rho_freshwater',gt=0)/g' ./classes/matdamageice.py:		
sed -i 's/(md,'fieldname','materials.mu_water','>',0)/(md,fieldname='materials.mu_water',gt=0)/g' ./classes/matdamageice.py:		
sed -i 's/(md,'fieldname','materials.rheology_B','>',0,'size',[md.mesh.numberofvertices])/(md,fieldname='materials.rheology_B',gt=0,size=[md.mesh.numberofvertices])/g' ./classes/matdamageice.py:		
sed -i 's/(md,'fieldname','materials.rheology_n','>',0,'size',[md.mesh.numberofelements])/(md,fieldname='materials.rheology_n',gt=0,size=[md.mesh.numberofelements])/g' ./classes/matdamageice.py:		
sed -i 's/(md,'fieldname','materials.rheology_law','values',['None','Cuffey','Paterson','Arrhenius','LliboutryDuval'])/(md,fieldname='materials.rheology_law',values=['None','Cuffey','Paterson','Arrhenius','LliboutryDuval'])/g' ./classes/matdamageice.py:		
sed -i 's/(md,'fieldname','materials.lithosphere_shear_modulus','>',0,numel=[1]);/(md,fieldname='materials.lithosphere_shear_modulus',gt=0,'numel',[1]);/g' ./classes/matdamageice.py:		
sed -i 's/(md,'fieldname','materials.lithosphere_density','>',0,numel=[1]);/(md,fieldname='materials.lithosphere_density',gt=0,'numel',[1]);/g' ./classes/matdamageice.py:		
sed -i 's/(md,'fieldname','materials.mantle_shear_modulus','>',0,numel=[1]);/(md,fieldname='materials.mantle_shear_modulus',gt=0,'numel',[1]);/g' ./classes/matdamageice.py:		
sed -i 's/(md,'fieldname','materials.mantle_density','>',0,'numel',[1]);/(md,fieldname='materials.mantle_density',gt=0,numel=[1]);/g' ./classes/matdamageice.py:		
sed -i 's///g' ./classes/massfluxatgate.py:from checkfield import checkfield
sed -i 's/(md,'field',self.definitionenum,'values',[Outputdefinition1Enum(),Outputdefinition2Enum(),Outputdefinition3Enum(),Outputdefinition4Enum(),Outputdefinition5Enum(),Outputdefinition6Enum(),Outputdefinition7Enum(),Outputdefinition8Enum(),Outputdefinition9Enum(),Outputdefinition10Enum()])/(md,'field',self.definitionenum,values=[Outputdefinition1Enum(),Outputdefinition2Enum(),Outputdefinition3Enum(),Outputdefinition4Enum(),Outputdefinition5Enum(),Outputdefinition6Enum(),Outputdefinition7Enum(),Outputdefinition8Enum(),Outputdefinition9Enum(),Outputdefinition10Enum()])/g' ./classes/massfluxatgate.py:			
sed -i 's///g' ./classes/gia.py:from checkfield import checkfield
sed -i 's/(md,'fieldname','gia.mantle_viscosity','NaN',1,'size',[md.mesh.numberofvertices,1],'>',0)/(md,fieldname='gia.mantle_viscosity',NaN=1,size=[md.mesh.numberofvertices,1],gt=0)/g' ./classes/gia.py:		
sed -i 's/(md,'fieldname','gia.lithosphere_thickness','NaN',1,'size',[md.mesh.numberofvertices,1],'>',0)/(md,fieldname='gia.lithosphere_thickness',NaN=1,size=[md.mesh.numberofvertices,1],gt=0)/g' ./classes/gia.py:		
sed -i 's/(md,'fieldname','gia.cross_section_shape','numel',[1],'values',[1,2])/(md,fieldname='gia.cross_section_shape',numel=[1],values=[1,2])/g' ./classes/gia.py:		
sed -i 's///g' ./classes/balancethickness.py:from checkfield import checkfield
sed -i 's/(md,'fieldname','balancethickness.spcthickness')/(md,fieldname='balancethickness.spcthickness')/g' ./classes/balancethickness.py:		
sed -i 's/(md,'fieldname','balancethickness.thickening_rate','size',[md.mesh.numberofvertices],'NaN',1)/(md,fieldname='balancethickness.thickening_rate',size=[md.mesh.numberofvertices],NaN=1)/g' ./classes/balancethickness.py:		
sed -i 's/(md,'fieldname','balancethickness.stabilization','size',[1],'values',[0,1,2,3])/(md,fieldname='balancethickness.stabilization',size=[1],values=[0,1,2,3])/g' ./classes/balancethickness.py:		
sed -i 's///g' ./classes/SMBgradients.py:from checkfield import checkfield
sed -i 's/(md,'fieldname','smb.href','timeseries',1,'NaN',1)/(md,fieldname='smb.href',timeseries=1,NaN=1)/g' ./classes/SMBgradients.py:			
sed -i 's/(md,'fieldname','smb.smbref','timeseries',1,'NaN',1)/(md,fieldname='smb.smbref',timeseries=1,NaN=1)/g' ./classes/SMBgradients.py:			
sed -i 's/(md,'fieldname','smb.b_pos','timeseries',1,'NaN',1)/(md,fieldname='smb.b_pos',timeseries=1,NaN=1)/g' ./classes/SMBgradients.py:			
sed -i 's/(md,'fieldname','smb.b_neg','timeseries',1,'NaN',1)/(md,fieldname='smb.b_neg',timeseries=1,NaN=1)/g' ./classes/SMBgradients.py:			
sed -i 's/(md,'fieldname','masstransport.requested_outputs','stringrow',1)/(md,fieldname='masstransport.requested_outputs','stringrow',1)/g' ./classes/SMBgradients.py:		
sed -i 's///g' ./classes/mesh2d.py:from checkfield import checkfield
sed -i 's/(md,'fieldname','mesh.x','NaN',1,'size',[md.mesh.numberofvertices])/(md,fieldname='mesh.x',NaN=1,size=[md.mesh.numberofvertices])/g' ./classes/mesh2d.py:		
sed -i 's/(md,'fieldname','mesh.y','NaN',1,'size',[md.mesh.numberofvertices])/(md,fieldname='mesh.y',NaN=1,size=[md.mesh.numberofvertices])/g' ./classes/mesh2d.py:		
sed -i 's/(md,'fieldname','mesh.elements','NaN',1,'>',0,'values',numpy.arange(1,md.mesh.numberofvertices+1))/(md,fieldname='mesh.elements',NaN=1,gt=0,values=numpy.arange(1,md.mesh.numberofvertices+1))/g' ./classes/mesh2d.py:		
sed -i 's/(md,'fieldname','mesh.elements','size',[md.mesh.numberofelements,3])/(md,fieldname='mesh.elements',size=[md.mesh.numberofelements,3])/g' ./classes/mesh2d.py:		
sed -i 's/(md,'fieldname','mesh.numberofelements','>',0)/(md,fieldname='mesh.numberofelements',gt=0)/g' ./classes/mesh2d.py:		
sed -i 's/(md,'fieldname','mesh.numberofvertices','>',0)/(md,fieldname='mesh.numberofvertices',gt=0)/g' ./classes/mesh2d.py:		
sed -i 's/(md,'fieldname','mesh.average_vertex_connectivity','>=',9,'message',"'mesh.average_vertex_connectivity' should be at least 9 in 2d")/(md,fieldname='mesh.average_vertex_connectivity',ge=9,'message',"'mesh.average_vertex_connectivity' should be at least 9 in 2d")/g' ./classes/mesh2d.py:		
sed -i 's/(md,'fieldname','inversion.iscontrol','values',[0, 1])/(md,fieldname='inversion.iscontrol',values=[0, 1])/g' ./classes/adinversion.py:		
sed -i 's/(md,'fieldname','inversion.control_parameters','cell',1,'values',\/(md,fieldname='inversion.control_parameters','cell',1,values=\/g' ./classes/adinversion.py:		
sed -i 's/(md,'fieldname','inversion.control_scaling_factors','size',[1, num_controls],'>',0,float('Nan'),1)/(md,fieldname='inversion.control_scaling_factors',size=[1, num_controls],gt=0,float('Nan'),1)/g' ./classes/adinversion.py:		
sed -i 's/(md,'fieldname','inversion.maxsteps','numel',1,'>=',0)/(md,fieldname='inversion.maxsteps',numel=1,ge=0)/g' ./classes/adinversion.py:		
sed -i 's/(md,'fieldname','inversion.maxiter','numel',1,'>=',0)/(md,fieldname='inversion.maxiter',numel=1,ge=0)/g' ./classes/adinversion.py:		
sed -i 's/(md,'fieldname','inversion.dxmin','numel',1,'>',0)/(md,fieldname='inversion.dxmin',numel=1,gt=0)/g' ./classes/adinversion.py:		
sed -i 's/(md,'fieldname','inversion.gttol','numel',1,'>',0)/(md,fieldname='inversion.gttol',numel=1,gt=0)/g' ./classes/adinversion.py:		
sed -i 's/(md,'fieldname','inversion.cost_functions','size',[1, num_costfunc],'values', [i for i in range(101,106)]+[201]+[i for i in range(501,507)]+[i for i in range(601,605)]+[i for i in range(1001, 1011)])/(md,fieldname='inversion.cost_functions',size=[1, num_costfunc],values= [i for i in range(101,106)]+[201]+[i for i in range(501,507)]+[i for i in range(601,605)]+[i for i in range(1001, 1011)])/g' ./classes/adinversion.py:		
sed -i 's/(md,'fieldname','inversion.cost_functions_coefficients','size',[md.mesh.numberofvertices, num_costfunc],'>=',0)/(md,fieldname='inversion.cost_functions_coefficients',size=[md.mesh.numberofvertices, num_costfunc],ge=0)/g' ./classes/adinversion.py:		
sed -i 's/(md,'fieldname','inversion.min_parameters','size',[md.mesh.numberofvertices, num_controls])/(md,fieldname='inversion.min_parameters',size=[md.mesh.numberofvertices, num_controls])/g' ./classes/adinversion.py:		
sed -i 's/(md,'fieldname','inversion.max_parameters','size',[md.mesh.numberofvertices, num_controls])/(md,fieldname='inversion.max_parameters',size=[md.mesh.numberofvertices, num_controls])/g' ./classes/adinversion.py:		
sed -i 's/(md,'fieldname','inversion.thickness_obs','size',[md.mesh.numberofvertices, 1],float('Nan'),1)/(md,fieldname='inversion.thickness_obs',size=[md.mesh.numberofvertices, 1],float('Nan'),1)/g' ./classes/adinversion.py:			
sed -i 's/(md,'fieldname','inversion.surface_obs','size',[md.mesh.numberofvertices, 1], float('Nan'),1)/(md,fieldname='inversion.surface_obs',size=[md.mesh.numberofvertices, 1], float('Nan'),1)/g' ./classes/adinversion.py:			
sed -i 's/(md,'fieldname','inversion.thickness_obs','size',[md.mesh.numberofvertices, 1],float('Nan'),1)/(md,fieldname='inversion.thickness_obs',size=[md.mesh.numberofvertices, 1],float('Nan'),1)/g' ./classes/adinversion.py:			
sed -i 's/(md,'fieldname','inversion.vx_obs','size',[md.mesh.numberofvertices, 1],float('Nan'),1)/(md,fieldname='inversion.vx_obs',size=[md.mesh.numberofvertices, 1],float('Nan'),1)/g' ./classes/adinversion.py:			
sed -i 's/(md,'fieldname','inversion.vy_obs','size',[md.mesh.numberofvertices, 1],float('Nan'),1)/(md,fieldname='inversion.vy_obs',size=[md.mesh.numberofvertices, 1],float('Nan'),1)/g' ./classes/adinversion.py:				
sed -i 's///g' ./classes/damage.py:from checkfield import checkfield
sed -i 's/(md,'fieldname','damage.isdamage','numel',[1],'values',[0,1])/(md,fieldname='damage.isdamage',numel=[1],values=[0,1])/g' ./classes/damage.py:		
sed -i 's/(md,'fieldname','damage.D','>=',0,'<=',self.max_damage,'size',[md.mesh.numberofvertices])/(md,fieldname='damage.D',ge=0,le=self.max_damage,size=[md.mesh.numberofvertices])/g' ./classes/damage.py:			
sed -i 's/(md,'fieldname','damage.max_damage','<',1,'>=',0)/(md,fieldname='damage.max_damage','<',1,ge=0)/g' ./classes/damage.py:			
sed -i 's/(md,'fieldname','damage.law','numel',[1],'values',[0,1,2,3])/(md,fieldname='damage.law',numel=[1],values=[0,1,2,3])/g' ./classes/damage.py:			
sed -i 's/(md,'fieldname','damage.spcdamage','timeseries',1)/(md,fieldname='damage.spcdamage',timeseries=1)/g' ./classes/damage.py:			
sed -i 's/(md,'fieldname','damage.stabilization','numel',[1],'values',[0,1,2,4])/(md,fieldname='damage.stabilization',numel=[1],values=[0,1,2,4])/g' ./classes/damage.py:			
sed -i 's/(md,'fieldname','damage.maxiter','>=0',0)/(md,fieldname='damage.maxiter','>=0',0)/g' ./classes/damage.py:			
sed -i 's/(md,'fieldname','damage.elementinterp','values',['P1','P2'])/(md,fieldname='damage.elementinterp',values=['P1','P2'])/g' ./classes/damage.py:			
sed -i 's/(md,'fieldname','damage.stress_threshold','>=',0)/(md,fieldname='damage.stress_threshold',ge=0)/g' ./classes/damage.py:			
sed -i 's/(md,'fieldname','damage.kappa','>',1)/(md,fieldname='damage.kappa',gt=1)/g' ./classes/damage.py:			
sed -i 's/(md,'fieldname','damage.healing','>=',0)/(md,fieldname='damage.healing',ge=0)/g' ./classes/damage.py:			
sed -i 's/(md,'fieldname','damage.c1','>=',0)/(md,fieldname='damage.c1',ge=0)/g' ./classes/damage.py:			
sed -i 's/(md,'fieldname','damage.c2','>=',0)/(md,fieldname='damage.c2',ge=0)/g' ./classes/damage.py:			
sed -i 's/(md,'fieldname','damage.c3','>=',0)/(md,fieldname='damage.c3',ge=0)/g' ./classes/damage.py:			
sed -i 's/(md,'fieldname','damage.c4','>=',0)/(md,fieldname='damage.c4',ge=0)/g' ./classes/damage.py:			
sed -i 's/(md,'fieldname','damage.healing','>=',0)/(md,fieldname='damage.healing',ge=0)/g' ./classes/damage.py:			
sed -i 's/(md,'fieldname','damage.equiv_stress','numel',[1],'values',[0,1])/(md,fieldname='damage.equiv_stress',numel=[1],values=[0,1])/g' ./classes/damage.py:			
sed -i 's/(md,'fieldname','damage.requested_outputs','stringrow',1)/(md,fieldname='damage.requested_outputs','stringrow',1)/g' ./classes/damage.py:			
sed -i 's///g' ./classes/friction.py:from checkfield import checkfield
sed -i 's/(md,'fieldname','friction.coefficient','timeseries',1,'NaN',1)/(md,fieldname='friction.coefficient',timeseries=1,NaN=1)/g' ./classes/friction.py:		
sed -i 's/(md,'fieldname','friction.q','NaN',1,'size',[md.mesh.numberofelements])/(md,fieldname='friction.q',NaN=1,size=[md.mesh.numberofelements])/g' ./classes/friction.py:		
sed -i 's/(md,'fieldname','friction.p','NaN',1,'size',[md.mesh.numberofelements])/(md,fieldname='friction.p',NaN=1,size=[md.mesh.numberofelements])/g' ./classes/friction.py:		
sed -i 's///g' ./classes/thermal.py:from checkfield import checkfield
sed -i 's/(md,'fieldname','thermal.stabilization','numel',[1],'values',[0,1,2])/(md,fieldname='thermal.stabilization',numel=[1],values=[0,1,2])/g' ./classes/thermal.py:		
sed -i 's/(md,'fieldname','thermal.spctemperature','timeseries',1)/(md,fieldname='thermal.spctemperature',timeseries=1)/g' ./classes/thermal.py:		
sed -i 's/(md,'fieldname','thermal.spctemperature[numpy.nonzero(numpy.logical_not(numpy.isnan(md.thermal.spctemperature[0:md.mesh.numberofvertices,:])))]','<',md.materials.meltingpoint-md.materials.beta*md.materials.rho_ice*md.constants.g*replicate[pos],'message',"spctemperature should be below the adjusted melting point")/(md,fieldname='thermal.spctemperature[numpy.nonzero(numpy.logical_not(numpy.isnan(md.thermal.spctemperature[0:md.mesh.numberofvertices,:])))]','<',md.materials.meltingpoint-md.materials.beta*md.materials.rho_ice*md.constants.g*replicate[pos],'message',"spctemperature should be below the adjusted melting point")/g' ./classes/thermal.py:			
sed -i 's/(md,'fieldname','thermal.isenthalpy','numel',[1],'values',[0,1])/(md,fieldname='thermal.isenthalpy',numel=[1],values=[0,1])/g' ./classes/thermal.py:			
sed -i 's/(md,'fieldname','thermal.isdynamicbasalspc','numel',[1],'values',[0,1]);/(md,fieldname='thermal.isdynamicbasalspc',numel=[1],values=[0,1]);/g' ./classes/thermal.py:			
sed -i 's/(md,'fieldname','thermal.reltol','>',0.,'message',"reltol must be larger than zero");/(md,fieldname='thermal.reltol',gt=0.,'message',"reltol must be larger than zero");/g' ./classes/thermal.py:				
sed -i 's/(md,'fieldname','thermal.requested_outputs','stringrow',1)/(md,fieldname='thermal.requested_outputs','stringrow',1)/g' ./classes/thermal.py:		
sed -i 's///g' ./classes/constants.py:from checkfield import checkfield
sed -i 's/(md,'fieldname','constants.g','>',0,'size',[1])/(md,fieldname='constants.g',gt=0,size=[1])/g' ./classes/constants.py:		
sed -i 's/(md,'fieldname','constants.yts','>',0,'size',[1])/(md,fieldname='constants.yts',gt=0,size=[1])/g' ./classes/constants.py:		
sed -i 's/(md,'fieldname','constants.referencetemperature','size',[1])/(md,fieldname='constants.referencetemperature',size=[1])/g' ./classes/constants.py:		
sed -i 's///g' ./classes/SMBforcing.py:from checkfield import checkfield
sed -i 's/(md,'fieldname','smb.mass_balance','timeseries',1,'NaN',1)/(md,fieldname='smb.mass_balance',timeseries=1,NaN=1)/g' ./classes/SMBforcing.py:			
sed -i 's/(md,'fieldname','smb.mass_balance','size',[md.mesh.numberofvertices],'NaN',1)/(md,fieldname='smb.mass_balance',size=[md.mesh.numberofvertices],NaN=1)/g' ./classes/SMBforcing.py:			
sed -i 's/(md,'fieldname','masstransport.requested_outputs','stringrow',1)/(md,fieldname='masstransport.requested_outputs','stringrow',1)/g' ./classes/SMBforcing.py:		
sed -i 's///g' ./classes/settings.py:from checkfield import checkfield
sed -i 's/(md,'fieldname','settings.results_on_nodes','numel',[1],'values',[0,1])/(md,fieldname='settings.results_on_nodes',numel=[1],values=[0,1])/g' ./classes/settings.py:		
sed -i 's/(md,'fieldname','settings.io_gather','numel',[1],'values',[0,1])/(md,fieldname='settings.io_gather',numel=[1],values=[0,1])/g' ./classes/settings.py:		
sed -i 's/(md,'fieldname','settings.lowmem','numel',[1],'values',[0,1])/(md,fieldname='settings.lowmem',numel=[1],values=[0,1])/g' ./classes/settings.py:		
sed -i 's/(md,'fieldname','settings.output_frequency','numel',[1],'>=',1)/(md,fieldname='settings.output_frequency',numel=[1],ge=1)/g' ./classes/settings.py:		
sed -i 's/(md,'fieldname','settings.recording_frequency','numel',[1],'>=',0)/(md,fieldname='settings.recording_frequency',numel=[1],ge=0)/g' ./classes/settings.py:		
sed -i 's/(md,'fieldname','settings.waitonlock','numel',[1])/(md,fieldname='settings.waitonlock',numel=[1])/g' ./classes/settings.py:		
sed -i 's///g' ./classes/SMBpdd.py:from checkfield import checkfield
sed -i 's/(md,'fieldname','smb.desfac','<=',1,'numel',[1])/(md,fieldname='smb.desfac',le=1,numel=[1])/g' ./classes/SMBpdd.py:			
sed -i 's/(md,'fieldname','smb.s0p','>=',0,'NaN',1,'size',[md.mesh.numberofvertices,1])/(md,fieldname='smb.s0p',ge=0,NaN=1,size=[md.mesh.numberofvertices,1])/g' ./classes/SMBpdd.py:			
sed -i 's/(md,'fieldname','smb.s0t','>=',0,'NaN',1,'size',[md.mesh.numberofvertices,1])/(md,fieldname='smb.s0t',ge=0,NaN=1,size=[md.mesh.numberofvertices,1])/g' ./classes/SMBpdd.py:			
sed -i 's/(md,'fieldname','smb.rlaps','>=',0,'numel',[1])/(md,fieldname='smb.rlaps',ge=0,numel=[1])/g' ./classes/SMBpdd.py:			
sed -i 's/(md,'fieldname','smb.rlapslgm','>=',0,'numel',[1])/(md,fieldname='smb.rlapslgm',ge=0,numel=[1])/g' ./classes/SMBpdd.py:			
sed -i 's/(md,'fieldname','smb.monthlytemperatures','NaN',1,'timeseries',1)/(md,fieldname='smb.monthlytemperatures',NaN=1,timeseries=1)/g' ./classes/SMBpdd.py:				
sed -i 's/(md,'fieldname','smb.precipitation','NaN',1,'timeseries',1)/(md,fieldname='smb.precipitation',NaN=1,timeseries=1)/g' ./classes/SMBpdd.py:				
sed -i 's/(md,'fieldname','smb.delta18o','NaN',1,'size',[2,numpy.nan],'singletimeseries',1)/(md,fieldname='smb.delta18o',NaN=1,size=[2,numpy.nan],'singletimeseries',1)/g' ./classes/SMBpdd.py:				
sed -i 's/(md,'fieldname','smb.delta18o_surface','NaN',1,'size',[2,numpy.nan],'singletimeseries',1)/(md,fieldname='smb.delta18o_surface',NaN=1,size=[2,numpy.nan],'singletimeseries',1)/g' ./classes/SMBpdd.py:				
sed -i 's/(md,'fieldname','smb.temperatures_presentday','size',[md.mesh.numberofvertices+1,12],'NaN',1,'timeseries',1)/(md,fieldname='smb.temperatures_presentday',size=[md.mesh.numberofvertices+1,12],NaN=1,timeseries=1)/g' ./classes/SMBpdd.py:				
sed -i 's/(md,'fieldname','smb.temperatures_lgm','size',[md.mesh.numberofvertices+1,12],'NaN',1,'timeseries',1)/(md,fieldname='smb.temperatures_lgm',size=[md.mesh.numberofvertices+1,12],NaN=1,timeseries=1)/g' ./classes/SMBpdd.py:				
sed -i 's/(md,'fieldname','smb.precipitations_presentday','size',[md.mesh.numberofvertices+1,12],'NaN',1,'timeseries',1)/(md,fieldname='smb.precipitations_presentday',size=[md.mesh.numberofvertices+1,12],NaN=1,timeseries=1)/g' ./classes/SMBpdd.py:				
sed -i 's/(md,'fieldname','smb.precipitations_lgm','size',[md.mesh.numberofvertices+1,12],'NaN',1,'timeseries',1)                                       /(md,fieldname='smb.precipitations_lgm','size',[md.mesh.numberofvertices+1,12],NaN=1,timeseries=1)                                       /g' ./classes/SMBpdd.py:				
sed -i 's/(md,'fieldname','smb.Tdiff','NaN',1,size=[2,numpy.nan],'singletimeseries',1)/(md,fieldname='smb.Tdiff',NaN=1,'size',[2,numpy.nan],'singletimeseries',1)/g' ./classes/SMBpdd.py:				
sed -i 's/(md,'fieldname','smb.sealev','NaN',1,size=[2,numpy.nan],'singletimeseries',1)/(md,fieldname='smb.sealev',NaN=1,'size',[2,numpy.nan],'singletimeseries',1)/g' ./classes/SMBpdd.py:				
sed -i 's/(md,'fieldname','smb.temperatures_presentday',size=[md.mesh.numberofvertices+1,12],'NaN',1,'timeseries',1)/(md,fieldname='smb.temperatures_presentday','size',[md.mesh.numberofvertices+1,12],NaN=1,timeseries=1)/g' ./classes/SMBpdd.py:				
sed -i 's/(md,'fieldname','smb.temperatures_lgm',size=[md.mesh.numberofvertices+1,12],'NaN',1,timeseries=1)/(md,fieldname='smb.temperatures_lgm','size',[md.mesh.numberofvertices+1,12],NaN=1,'timeseries',1)/g' ./classes/SMBpdd.py:				
sed -i 's/(md,'fieldname','smb.precipitations_presentday',size=[md.mesh.numberofvertices+1,12],'NaN',1,timeseries=1)/(md,fieldname='smb.precipitations_presentday','size',[md.mesh.numberofvertices+1,12],NaN=1,'timeseries',1)/g' ./classes/SMBpdd.py:				
sed -i 's/(md,'fieldname','smb.precipitations_lgm',size=[md.mesh.numberofvertices+1,12],'NaN',1,timeseries=1)                                       /(md,fieldname='smb.precipitations_lgm','size',[md.mesh.numberofvertices+1,12],NaN=1,'timeseries',1)                                       /g' ./classes/SMBpdd.py:				
sed -i 's/(md,'fieldname','smb.Pfac','NaN',1,size=[2,numpy.nan],'singletimeseries',1)/(md,fieldname='smb.Pfac',NaN=1,'size',[2,numpy.nan],'singletimeseries',1)/g' ./classes/SMBpdd.py:				
sed -i 's/(md,'fieldname','smb.Tdiff','NaN',1,size=[2,numpy.nan],'singletimeseries',1)/(md,fieldname='smb.Tdiff',NaN=1,'size',[2,numpy.nan],'singletimeseries',1)/g' ./classes/SMBpdd.py:				
sed -i 's/(md,'fieldname','smb.sealev','NaN',1,size=[2,numpy.nan],'singletimeseries',1)/(md,fieldname='smb.sealev',NaN=1,'size',[2,numpy.nan],'singletimeseries',1)/g' ./classes/SMBpdd.py:				
sed -i 's/(md,'fieldname','masstransport.requested_outputs','stringrow',1)/(md,fieldname='masstransport.requested_outputs','stringrow',1)/g' ./classes/SMBpdd.py:		
sed -i 's///g' ./classes/toolkits.py:from checkfield import checkfield
sed -i 's///g' ./classes/transient.py:from checkfield import checkfield
sed -i 's/(md,'fieldname','transient.issmb','numel',[1],'values',[0,1])/(md,fieldname='transient.issmb',numel=[1],values=[0,1])/g' ./classes/transient.py:		
sed -i 's/(md,'fieldname','transient.ismasstransport','numel',[1],'values',[0,1])/(md,fieldname='transient.ismasstransport',numel=[1],values=[0,1])/g' ./classes/transient.py:		
sed -i 's/(md,'fieldname','transient.isstressbalance','numel',[1],'values',[0,1])/(md,fieldname='transient.isstressbalance',numel=[1],values=[0,1])/g' ./classes/transient.py:		
sed -i 's/(md,'fieldname','transient.isthermal','numel',[1],'values',[0,1])/(md,fieldname='transient.isthermal',numel=[1],values=[0,1])/g' ./classes/transient.py:		
sed -i 's/(md,'fieldname','transient.isgroundingline','numel',[1],'values',[0,1])/(md,fieldname='transient.isgroundingline',numel=[1],values=[0,1])/g' ./classes/transient.py:		
sed -i 's/(md,'fieldname','transient.isgia','numel',[1],'values',[0,1])/(md,fieldname='transient.isgia',numel=[1],values=[0,1])/g' ./classes/transient.py:		
sed -i 's/(md,'fieldname','transient.isdamageevolution','numel',[1],'values',[0,1])/(md,fieldname='transient.isdamageevolution',numel=[1],values=[0,1])/g' ./classes/transient.py:		
sed -i 's/(md,'fieldname','transient.islevelset','numel',[1],'values',[0,1])/(md,fieldname='transient.islevelset',numel=[1],values=[0,1])/g' ./classes/transient.py:		
sed -i 's/(md,'fieldname','transient.ishydrology','numel',[1],'values',[0,1])/(md,fieldname='transient.ishydrology',numel=[1],values=[0,1])/g' ./classes/transient.py:		
sed -i 's/(md,'fieldname','transient.iscalving','numel',[1],'values',[0,1]);/(md,'fieldname','transient.iscalving',numel=[1],values=[0,1]);/g' ./classes/transient.py:		
sed -i 's/(md,'fieldname','transient.requested_outputs','stringrow',1)/(md,fieldname='transient.requested_outputs','stringrow',1)/g' ./classes/transient.py:		
sed -i 's///g' ./classes/basalforcings.py:from checkfield import checkfield
sed -i 's/(md,'fieldname','basalforcings.groundedice_melting_rate','NaN',1,'timeseries',1)/(md,fieldname='basalforcings.groundedice_melting_rate',NaN=1,timeseries=1)/g' ./classes/basalforcings.py:			
sed -i 's/(md,'fieldname','basalforcings.floatingice_melting_rate','NaN',1,'timeseries',1)/(md,fieldname='basalforcings.floatingice_melting_rate',NaN=1,timeseries=1)/g' ./classes/basalforcings.py:			
sed -i 's/(md,'fieldname','basalforcings.groundedice_melting_rate','NaN',1,size=[md.mesh.numberofvertices])/(md,fieldname='basalforcings.groundedice_melting_rate',NaN=1,'size',[md.mesh.numberofvertices])/g' ./classes/basalforcings.py:			
sed -i 's/(md,'fieldname','basalforcings.floatingice_melting_rate','NaN',1,size=[md.mesh.numberofvertices])/(md,fieldname='basalforcings.floatingice_melting_rate',NaN=1,'size',[md.mesh.numberofvertices])/g' ./classes/basalforcings.py:			
sed -i 's/(md,'fieldname','basalforcings.groundedice_melting_rate','NaN',1,'timeseries',1)/(md,fieldname='basalforcings.groundedice_melting_rate',NaN=1,timeseries=1)/g' ./classes/basalforcings.py:			
sed -i 's/(md,'fieldname','basalforcings.floatingice_melting_rate','NaN',1,'timeseries',1)/(md,fieldname='basalforcings.floatingice_melting_rate',NaN=1,timeseries=1)/g' ./classes/basalforcings.py:			
sed -i 's/(md,'fieldname','basalforcings.geothermalflux','NaN',1,'timeseries',1,'>=',0)/(md,fieldname='basalforcings.geothermalflux',NaN=1,timeseries=1,ge=0)/g' ./classes/basalforcings.py:			
sed -i 's///g' ./classes/stressbalance.py:from checkfield import checkfield
sed -i 's/(md,'fieldname','stressbalance.spcvx','timeseries',1)/(md,fieldname='stressbalance.spcvx',timeseries=1)/g' ./classes/stressbalance.py:		
sed -i 's/(md,'fieldname','stressbalance.spcvy','timeseries',1)/(md,fieldname='stressbalance.spcvy',timeseries=1)/g' ./classes/stressbalance.py:		
sed -i 's/(md,'fieldname','stressbalance.spcvz','timeseries',1)/(md,fieldname='stressbalance.spcvz',timeseries=1)/g' ./classes/stressbalance.py:			
sed -i 's/(md,'fieldname','stressbalance.restol',size=[1],'>',0)/(md,fieldname='stressbalance.restol','size',[1],gt=0)/g' ./classes/stressbalance.py:		
sed -i 's/(md,'fieldname','stressbalance.reltol',size=[1])/(md,fieldname='stressbalance.reltol','size',[1])/g' ./classes/stressbalance.py:		
sed -i 's/(md,'fieldname','stressbalance.abstol',size=[1])/(md,fieldname='stressbalance.abstol','size',[1])/g' ./classes/stressbalance.py:		
sed -i 's/(md,'fieldname','stressbalance.isnewton','numel',[1],'values',[0,1,2])/(md,fieldname='stressbalance.isnewton',numel=[1],values=[0,1,2])/g' ./classes/stressbalance.py:		
sed -i 's/(md,'fieldname','stressbalance.FSreconditioning',size=[1],'NaN',1)/(md,fieldname='stressbalance.FSreconditioning','size',[1],NaN=1)/g' ./classes/stressbalance.py:		
sed -i 's/(md,'fieldname','stressbalance.viscosity_overshoot',size=[1],'NaN',1)/(md,fieldname='stressbalance.viscosity_overshoot','size',[1],NaN=1)/g' ./classes/stressbalance.py:		
sed -i 's/(md,'fieldname','stressbalance.maxiter',size=[1],'>=',1)/(md,fieldname='stressbalance.maxiter','size',[1],ge=1)/g' ./classes/stressbalance.py:		
sed -i 's/(md,'fieldname','stressbalance.referential',size=[md.mesh.numberofvertices,6])/(md,fieldname='stressbalance.referential','size',[md.mesh.numberofvertices,6])/g' ./classes/stressbalance.py:		
sed -i 's/(md,'fieldname','stressbalance.loadingforce',size=[md.mesh.numberofvertices,3])/(md,fieldname='stressbalance.loadingforce','size',[md.mesh.numberofvertices,3])/g' ./classes/stressbalance.py:		
sed -i 's/(md,'fieldname','stressbalance.requested_outputs','stringrow',1);/(md,'fieldname','stressbalance.requested_outputs','stringrow',1);/g' ./classes/stressbalance.py:		
sed -i 's///g' ./classes/initialization.py:from checkfield import checkfield
sed -i 's/(md,'fieldname','initialization.vx','NaN',1,size=[md.mesh.numberofvertices])/(md,fieldname='initialization.vx',NaN=1,'size',[md.mesh.numberofvertices])/g' ./classes/initialization.py:				
sed -i 's/(md,'fieldname','initialization.vy','NaN',1,size=[md.mesh.numberofvertices])/(md,fieldname='initialization.vy',NaN=1,'size',[md.mesh.numberofvertices])/g' ./classes/initialization.py:				
sed -i 's/(md,'fieldname','initialization.vx','NaN',1,size=[md.mesh.numberofvertices])/(md,fieldname='initialization.vx',NaN=1,'size',[md.mesh.numberofvertices])/g' ./classes/initialization.py:			
sed -i 's/(md,'fieldname','initialization.vy','NaN',1,size=[md.mesh.numberofvertices])/(md,fieldname='initialization.vy',NaN=1,'size',[md.mesh.numberofvertices])/g' ./classes/initialization.py:			
sed -i 's/(md,'fieldname','initialization.vx','NaN',1,size=[md.mesh.numberofvertices])/(md,fieldname='initialization.vx',NaN=1,'size',[md.mesh.numberofvertices])/g' ./classes/initialization.py:			
sed -i 's/(md,'fieldname','initialization.vy','NaN',1,size=[md.mesh.numberofvertices])/(md,fieldname='initialization.vy',NaN=1,'size',[md.mesh.numberofvertices])/g' ./classes/initialization.py:			
sed -i 's/(md,'fieldname','initialization.vx','NaN',1,size=[md.mesh.numberofvertices])/(md,fieldname='initialization.vx',NaN=1,'size',[md.mesh.numberofvertices])/g' ./classes/initialization.py:			
sed -i 's/(md,'fieldname','initialization.vy','NaN',1,size=[md.mesh.numberofvertices])/(md,fieldname='initialization.vy',NaN=1,'size',[md.mesh.numberofvertices])/g' ./classes/initialization.py:			
sed -i 's/(md,'fieldname','initialization.temperature','NaN',1,size=[md.mesh.numberofvertices])/(md,fieldname='initialization.temperature',NaN=1,'size',[md.mesh.numberofvertices])/g' ./classes/initialization.py:			
sed -i 's/(md,'fieldname','initialization.vz','NaN',1,size=[md.mesh.numberofvertices])/(md,fieldname='initialization.vz',NaN=1,'size',[md.mesh.numberofvertices])/g' ./classes/initialization.py:				
sed -i 's/(md,'fieldname','initialization.pressure','NaN',1,'size',[md.mesh.numberofvertices])/(md,fieldname='initialization.pressure',NaN=1,size=[md.mesh.numberofvertices])/g' ./classes/initialization.py:			
sed -i 's/(md,'fieldname','initialization.waterfraction','>=',0,'size',[md.mesh.numberofvertices])/(md,fieldname='initialization.waterfraction',ge=0,size=[md.mesh.numberofvertices])/g' ./classes/initialization.py:				
sed -i 's/(md,'fieldname','initialization.watercolumn'  ,'>=',0,'size',[md.mesh.numberofvertices])/(md,fieldname='initialization.watercolumn'  ,ge=0,size=[md.mesh.numberofvertices])/g' ./classes/initialization.py:				
sed -i 's/(md,'fieldname','initialization.watercolumn','NaN',1,'size',[md.mesh.numberofvertices])/(md,fieldname='initialization.watercolumn',NaN=1,size=[md.mesh.numberofvertices])/g' ./classes/initialization.py:				
sed -i 's/(md,'fieldname','initialization.sediment_head','NaN',1,'size',[md.mesh.numberofvertices,1])/(md,fieldname='initialization.sediment_head',NaN=1,size=[md.mesh.numberofvertices,1])/g' ./classes/initialization.py:				
sed -i 's/(md,'fieldname','initialization.epl_head','NaN',1,'size',[md.mesh.numberofvertices,1])/(md,fieldname='initialization.epl_head',NaN=1,size=[md.mesh.numberofvertices,1])/g' ./classes/initialization.py:					
sed -i 's/(md,'fieldname','initialization.epl_thickness','NaN',1,'size',[md.mesh.numberofvertices,1])/(md,fieldname='initialization.epl_thickness',NaN=1,size=[md.mesh.numberofvertices,1])/g' ./classes/initialization.py:					
sed -i 's///g' ./classes/hydrologydc.py:from checkfield import checkfield
sed -i 's/(md,'fieldname','hydrology.water_compressibility','numel',[1],'>',0.)/(md,fieldname='hydrology.water_compressibility',numel=[1],gt=0.)/g' ./classes/hydrologydc.py:		
sed -i 's/(md,'fieldname','hydrology.isefficientlayer','numel',[1],'values',[0,1])/(md,fieldname='hydrology.isefficientlayer',numel=[1],values=[0,1])/g' ./classes/hydrologydc.py:		
sed -i 's/(md,'fieldname','hydrology.penalty_factor','>',0.,'numel',[1])/(md,fieldname='hydrology.penalty_factor',gt=0.,numel=[1])/g' ./classes/hydrologydc.py:		
sed -i 's/(md,'fieldname','hydrology.penalty_lock','>=',0.,'numel',[1])/(md,fieldname='hydrology.penalty_lock',ge=0.,numel=[1])/g' ./classes/hydrologydc.py:		
sed -i 's/(md,'fieldname','hydrology.rel_tol','>',0.,'numel',[1])/(md,fieldname='hydrology.rel_tol',gt=0.,numel=[1])/g' ./classes/hydrologydc.py:		
sed -i 's/(md,'fieldname','hydrology.max_iter','>',0.,'numel',[1])/(md,fieldname='hydrology.max_iter',gt=0.,numel=[1])/g' ./classes/hydrologydc.py:		
sed -i 's/(md,'fieldname','hydrology.sedimentlimit_flag','numel',[1],'values',[0,1,2,3])/(md,fieldname='hydrology.sedimentlimit_flag',numel=[1],values=[0,1,2,3])/g' ./classes/hydrologydc.py:		
sed -i 's/(md,'fieldname','hydrology.transfer_flag','numel',[1],'values',[0,1])/(md,fieldname='hydrology.transfer_flag',numel=[1],values=[0,1])/g' ./classes/hydrologydc.py:		
sed -i 's/(md,'fieldname','hydrology.sedimentlimit','>',0.,'numel',[1])/(md,fieldname='hydrology.sedimentlimit',gt=0.,numel=[1])/g' ./classes/hydrologydc.py:			
sed -i 's/(md,'fieldname','hydrology.leakage_factor','>',0.,'numel',[1])/(md,fieldname='hydrology.leakage_factor',gt=0.,numel=[1])/g' ./classes/hydrologydc.py:			
sed -i 's/(md,'fieldname','hydrology.basal_moulin_input','NaN',1,'timeseries',1)/(md,fieldname='hydrology.basal_moulin_input',NaN=1,timeseries=1)/g' ./classes/hydrologydc.py:		
sed -i 's/(md,'fieldname','hydrology.spcsediment_head','timeseries',1)/(md,fieldname='hydrology.spcsediment_head',timeseries=1)/g' ./classes/hydrologydc.py:		
sed -i 's/(md,'fieldname','hydrology.sediment_compressibility','>',0.,'numel',[1])/(md,fieldname='hydrology.sediment_compressibility','>',0.,numel=[1])/g' ./classes/hydrologydc.py:		
sed -i 's/(md,'fieldname','hydrology.sediment_porosity',gt=0.,'numel',[1])/(md,fieldname='hydrology.sediment_porosity','>',0.,numel=[1])/g' ./classes/hydrologydc.py:		
sed -i 's/(md,'fieldname','hydrology.sediment_thickness',gt=0.,'numel',[1])/(md,fieldname='hydrology.sediment_thickness','>',0.,numel=[1])/g' ./classes/hydrologydc.py:		
sed -i 's/(md,'fieldname','hydrology.sediment_transmitivity','>=',0,'size',[md.mesh.numberofvertices,1])/(md,fieldname='hydrology.sediment_transmitivity',ge=0,size=[md.mesh.numberofvertices,1])/g' ./classes/hydrologydc.py:		
sed -i 's/(md,'fieldname','hydrology.spcepl_head','timeseries',1)/(md,fieldname='hydrology.spcepl_head',timeseries=1)/g' ./classes/hydrologydc.py:			
sed -i 's/(md,'fieldname','hydrology.mask_eplactive_node','size',[md.mesh.numberofvertices,1],'values',[0,1])/(md,fieldname='hydrology.mask_eplactive_node',size=[md.mesh.numberofvertices,1],values=[0,1])/g' ./classes/hydrologydc.py:			
sed -i 's/(md,'fieldname','hydrology.epl_compressibility','>',0.,'numel',[1])/(md,fieldname='hydrology.epl_compressibility',gt=0.,numel=[1])/g' ./classes/hydrologydc.py:			
sed -i 's/(md,'fieldname','hydrology.epl_porosity','>',0.,'numel',[1])/(md,fieldname='hydrology.epl_porosity',gt=0.,numel=[1])/g' ./classes/hydrologydc.py:			
sed -i 's/(md,'fieldname','hydrology.epl_max_thickness','numel',[1],'>',0.)/(md,fieldname='hydrology.epl_max_thickness',numel=[1],gt=0.)/g' ./classes/hydrologydc.py:			
sed -i 's/(md,'fieldname','hydrology.epl_initial_thickness','numel',[1],'>',0.)/(md,fieldname='hydrology.epl_initial_thickness',numel=[1],gt=0.)/g' ./classes/hydrologydc.py:			
sed -i 's/(md,'fieldname','hydrology.epl_colapse_thickness','numel',[1],'>',0.)/(md,fieldname='hydrology.epl_colapse_thickness',numel=[1],gt=0.)/g' ./classes/hydrologydc.py:			
sed -i 's/(md,'fieldname','hydrology.epl_thick_comp','numel',[1],'values',[0,1])/(md,fieldname='hydrology.epl_thick_comp',numel=[1],values=[0,1])/g' ./classes/hydrologydc.py:			
sed -i 's/(md,'fieldname','hydrology.eplflip_lock','>=',0.,'numel',[1])/(md,fieldname='hydrology.eplflip_lock',ge=0.,numel=[1])/g' ./classes/hydrologydc.py:			
sed -i 's/(md,'fieldname','hydrology.epl_conductivity','numel',[1],'>',0.)/(md,fieldname='hydrology.epl_conductivity',numel=[1],gt=0.)/g' ./classes/hydrologydc.py:			
sed -i 's///g' ./classes/linearbasalforcings.py:from checkfield import checkfield
sed -i 's/(md,'fieldname','basalforcings.groundedice_melting_rate','NaN',1,'timeseries',1)/(md,fieldname='basalforcings.groundedice_melting_rate',NaN=1,timeseries=1)/g' ./classes/linearbasalforcings.py:			
sed -i 's/(md,'fieldname','basalforcings.deepwater_melting_rate','>=',0);/(md,'fieldname','basalforcings.deepwater_melting_rate',ge=0);/g' ./classes/linearbasalforcings.py:			
sed -i 's/(md,'fieldname','basalforcings.deepwater_elevation','<',md.basalforcings.upperwater_elevation);/(md,'fieldname','basalforcings.deepwater_elevation','<',md.basalforcings.upperwater_elevation);/g' ./classes/linearbasalforcings.py:			
sed -i 's/(md,'fieldname','basalforcings.upperwater_elevation','<',0);/(md,'fieldname','basalforcings.upperwater_elevation','<',0);/g' ./classes/linearbasalforcings.py:			
sed -i 's/(md,'fieldname','basalforcings.groundedice_melting_rate','NaN',1,'size',[md.mesh.numberofvertices])/(md,fieldname='basalforcings.groundedice_melting_rate',NaN=1,size=[md.mesh.numberofvertices])/g' ./classes/linearbasalforcings.py:			
sed -i 's/(md,'fieldname','basalforcings.deepwater_melting_rate','>=',0);/(md,'fieldname','basalforcings.deepwater_melting_rate',ge=0);/g' ./classes/linearbasalforcings.py:			
sed -i 's/(md,'fieldname','basalforcings.deepwater_elevation','<',md.basalforcings.upperwater_elevation);/(md,'fieldname','basalforcings.deepwater_elevation','<',md.basalforcings.upperwater_elevation);/g' ./classes/linearbasalforcings.py:			
sed -i 's/(md,'fieldname','basalforcings.upperwater_elevation','<',0);/(md,'fieldname','basalforcings.upperwater_elevation','<',0);/g' ./classes/linearbasalforcings.py:			
sed -i 's/(md,'fieldname','basalforcings.groundedice_melting_rate','NaN',1,'timeseries',1)/(md,fieldname='basalforcings.groundedice_melting_rate',NaN=1,timeseries=1)/g' ./classes/linearbasalforcings.py:			
sed -i 's/(md,'fieldname','basalforcings.deepwater_melting_rate','>=',0);/(md,'fieldname','basalforcings.deepwater_melting_rate',ge=0);/g' ./classes/linearbasalforcings.py:			
sed -i 's/(md,'fieldname','basalforcings.deepwater_elevation','<',md.basalforcings.upperwater_elevation);/(md,'fieldname','basalforcings.deepwater_elevation','<',md.basalforcings.upperwater_elevation);/g' ./classes/linearbasalforcings.py:			
sed -i 's/(md,'fieldname','basalforcings.upperwater_elevation','<',0);/(md,'fieldname','basalforcings.upperwater_elevation','<',0);/g' ./classes/linearbasalforcings.py:			
sed -i 's/(md,'fieldname','basalforcings.geothermalflux','NaN',1,'timeseries',1,'>=',0)/(md,fieldname='basalforcings.geothermalflux',NaN=1,timeseries=1,ge=0)/g' ./classes/linearbasalforcings.py:			
sed -i 's///g' ./classes/outputdefinition.py:from checkfield import checkfield
sed -i 's/(md,'fieldname','outputdefinition.definitions','cell',1)/(md,fieldname='outputdefinition.definitions','cell',1)/g' ./classes/outputdefinition.py:		
sed -i 's///g' ./classes/frictionweertman.py:from checkfield import checkfield
sed -i 's/(md,'fieldname','friction.C','timeseries',1,'NaN',1)/(md,fieldname='friction.C',timeseries=1,NaN=1)/g' ./classes/frictionweertman.py:		
sed -i 's/(md,'fieldname','friction.m','NaN',1,'size',[md.mesh.numberofelements])/(md,fieldname='friction.m',NaN=1,size=[md.mesh.numberofelements])/g' ./classes/frictionweertman.py:		
sed -i 's///g' ./classes/miscellaneous.py:from checkfield import checkfield
sed -i 's/(md,'fieldname','miscellaneous.name','empty',1)/(md,fieldname='miscellaneous.name','empty',1)/g' ./classes/miscellaneous.py:		
sed -i 's///g' ./classes/mask.py:from checkfield import checkfield
sed -i 's/(md,'fieldname','mask.ice_levelset'        ,'size',[md.mesh.numberofvertices])/(md,fieldname='mask.ice_levelset'        ,size=[md.mesh.numberofvertices])/g' ./classes/mask.py:		
sed -i 's///g' ./classes/hydrologyshreve.py:from checkfield import checkfield
sed -i 's/(md,'fieldname','hydrology.spcwatercolumn','timeseries',1)/(md,fieldname='hydrology.spcwatercolumn',timeseries=1)/g' ./classes/hydrologyshreve.py:		
sed -i 's/(md,'fieldname','hydrology.stabilization','>=',0)/(md,fieldname='hydrology.stabilization',ge=0)/g' ./classes/hydrologyshreve.py:		
sed -i 's///g' ./classes/private.py:from checkfield import checkfield
sed -i 's///g' ./classes/rifts.py:from checkfield import checkfield
sed -i 's/(md,'fieldname',"rifts.riftstruct[%d]['fill']" % i,'values',[WaterEnum(),AirEnum(),IceEnum(),MelangeEnum()])/(md,fieldname="rifts.riftstruct[%d]['fill']" % i,values=[WaterEnum(),AirEnum(),IceEnum(),MelangeEnum()])/g' ./classes/rifts.py:				
sed -i 's///g' ./classes/groundingline.py:from checkfield import checkfield
sed -i 's/(md,fieldname='groundingline.migration',values=['None','AggressiveMigration','SoftMigration','SubelementMigration','SubelementMigration2','Contact','GroundingOnly'])/(md,fieldname='groundingline.migration',values=['None','AggressiveMigration','SoftMigration','SubelementMigration','SubelementMigration2','Contact','GroundingOnly'])/g' ./classes/groundingline.py:		
sed -i 's///g' ./classes/taoinversion.py:from checkfield import checkfield
sed -i 's/(md,'fieldname','inversion.iscontrol','values',[0, 1])/(md,fieldname='inversion.iscontrol',values=[0, 1])/g' ./classes/taoinversion.py:		
sed -i 's/(md,'fieldname','inversion.incomplete_adjoint','values',[0, 1])/(md,fieldname='inversion.incomplete_adjoint',values=[0, 1])/g' ./classes/taoinversion.py:		
sed -i 's/(md,'fieldname','inversion.control_parameters','cell',1,'values',supportedcontrols())/(md,fieldname='inversion.control_parameters','cell',1,values=supportedcontrols())/g' ./classes/taoinversion.py:		
sed -i 's/(md,'fieldname','inversion.maxsteps','numel',1,'>=',0)/(md,fieldname='inversion.maxsteps',numel=1,ge=0)/g' ./classes/taoinversion.py:		
sed -i 's/(md,'fieldname','inversion.maxiter','numel',1,'>=',0)/(md,fieldname='inversion.maxiter',numel=1,ge=0)/g' ./classes/taoinversion.py:		
sed -i 's/(md,'fieldname','inversion.fatol','numel',1,'>=',0)/(md,fieldname='inversion.fatol',numel=1,ge=0)/g' ./classes/taoinversion.py:		
sed -i 's/(md,'fieldname','inversion.frtol','numel',1,'>=',0)/(md,fieldname='inversion.frtol',numel=1,ge=0)/g' ./classes/taoinversion.py:		
sed -i 's/(md,'fieldname','inversion.gatol','numel',1,'>=',0)/(md,fieldname='inversion.gatol',numel=1,ge=0)/g' ./classes/taoinversion.py:		
sed -i 's/(md,'fieldname','inversion.grtol','numel',1,'>=',0)/(md,fieldname='inversion.grtol',numel=1,ge=0)/g' ./classes/taoinversion.py:		
sed -i 's/(md,'fieldname','inversion.gttol','numel',1,'>=',0)/(md,fieldname='inversion.gttol',numel=1,ge=0)/g' ./classes/taoinversion.py:		
sed -i 's/(md,'fieldname','inversion.algorithm','values',{'blmvm','cg','lmvm'})/(md,fieldname='inversion.algorithm',values={'blmvm','cg','lmvm'})/g' ./classes/taoinversion.py:			
sed -i 's/(md,'fieldname','inversion.algorithm','values',{'tao_blmvm','tao_cg','tao_lmvm'})/(md,fieldname='inversion.algorithm',values={'tao_blmvm','tao_cg','tao_lmvm'})/g' ./classes/taoinversion.py:			
sed -i 's/(md,'fieldname','inversion.cost_functions','size',[1, num_costfunc],'values',supportedcostfunctions())/(md,fieldname='inversion.cost_functions',size=[1, num_costfunc],values=supportedcostfunctions())/g' ./classes/taoinversion.py:		
sed -i 's/(md,'fieldname','inversion.cost_functions_coefficients','size',[md.mesh.numberofvertices, num_costfunc],'>=',0)/(md,fieldname='inversion.cost_functions_coefficients',size=[md.mesh.numberofvertices, num_costfunc],ge=0)/g' ./classes/taoinversion.py:		
sed -i 's/(md,'fieldname','inversion.min_parameters','size',[md.mesh.numberofvertices, num_controls])/(md,fieldname='inversion.min_parameters',size=[md.mesh.numberofvertices, num_controls])/g' ./classes/taoinversion.py:		
sed -i 's/(md,'fieldname','inversion.max_parameters','size',[md.mesh.numberofvertices, num_controls])/(md,fieldname='inversion.max_parameters',size=[md.mesh.numberofvertices, num_controls])/g' ./classes/taoinversion.py:		
sed -i 's/(md,'fieldname','inversion.thickness_obs','size',[md.mesh.numberofvertices, 1],'NaN',1)/(md,fieldname='inversion.thickness_obs',size=[md.mesh.numberofvertices, 1],NaN=1)/g' ./classes/taoinversion.py:			
sed -i 's/(md,'fieldname','inversion.thickness_obs','size',[md.mesh.numberofvertices, 1],'NaN',1)/(md,fieldname='inversion.thickness_obs',size=[md.mesh.numberofvertices, 1],NaN=1)/g' ./classes/taoinversion.py:			
sed -i 's/(md,'fieldname','inversion.vx_obs','size',[md.mesh.numberofvertices, 1],'NaN',1)/(md,fieldname='inversion.vx_obs',size=[md.mesh.numberofvertices, 1],NaN=1)/g' ./classes/taoinversion.py:			
sed -i 's/(md,'fieldname','inversion.vy_obs','size',[md.mesh.numberofvertices, 1],'NaN',1)/(md,fieldname='inversion.vy_obs',size=[md.mesh.numberofvertices, 1],NaN=1)/g' ./classes/taoinversion.py:			
sed -i 's///g' ./classes/flowequation.py:from checkfield import checkfield
sed -i 's/(md,'fieldname','flowequation.isSIA','numel',[1],'values',[0,1])/(md,fieldname='flowequation.isSIA',numel=[1],values=[0,1])/g' ./classes/flowequation.py:		
sed -i 's/(md,'fieldname','flowequation.isSSA','numel',[1],'values',[0,1])/(md,fieldname='flowequation.isSSA',numel=[1],values=[0,1])/g' ./classes/flowequation.py:		
sed -i 's/(md,'fieldname','flowequation.isL1L2','numel',[1],'values',[0,1])/(md,fieldname='flowequation.isL1L2',numel=[1],values=[0,1])/g' ./classes/flowequation.py:		
sed -i 's/(md,'fieldname','flowequation.isHO','numel',[1],'values',[0,1])/(md,fieldname='flowequation.isHO',numel=[1],values=[0,1])/g' ./classes/flowequation.py:		
sed -i 's/(md,'fieldname','flowequation.isFS','numel',[1],'values',[0,1])/(md,fieldname='flowequation.isFS',numel=[1],values=[0,1])/g' ./classes/flowequation.py:		
sed -i 's/(md,'fieldname','flowequation.fe_SSA','values',['P1','P1bubble','P1bubblecondensed','P2','P2bubble'])/(md,fieldname='flowequation.fe_SSA',values=['P1','P1bubble','P1bubblecondensed','P2','P2bubble'])/g' ./classes/flowequation.py:		
sed -i 's/(md,'fieldname','flowequation.fe_HO' ,'values',['P1','P1bubble','P1bubblecondensed','P1xP2','P2xP1','P2','P2bubble','P1xP3','P2xP4'])/(md,fieldname='flowequation.fe_HO' ,values=['P1','P1bubble','P1bubblecondensed','P1xP2','P2xP1','P2','P2bubble','P1xP3','P2xP4'])/g' ./classes/flowequation.py:		
sed -i 's/(md,'fieldname','flowequation.fe_FS' ,'values',['P1P1','P1P1GLS','MINIcondensed','MINI','TaylorHood','XTaylorHood','OneLayerP4z','CrouzeixRaviart'])/(md,fieldname='flowequation.fe_FS' ,values=['P1P1','P1P1GLS','MINIcondensed','MINI','TaylorHood','XTaylorHood','OneLayerP4z','CrouzeixRaviart'])/g' ./classes/flowequation.py:		
sed -i 's/(md,'fieldname','flowequation.borderSSA','size',[md.mesh.numberofvertices],'values',[0,1])/(md,fieldname='flowequation.borderSSA',size=[md.mesh.numberofvertices],values=[0,1])/g' ./classes/flowequation.py:		
sed -i 's/(md,'fieldname','flowequation.borderHO','size',[md.mesh.numberofvertices],'values',[0,1])/(md,fieldname='flowequation.borderHO',size=[md.mesh.numberofvertices],values=[0,1])/g' ./classes/flowequation.py:		
sed -i 's/(md,'fieldname','flowequation.borderFS','size',[md.mesh.numberofvertices],'values',[0,1])/(md,fieldname='flowequation.borderFS',size=[md.mesh.numberofvertices],values=[0,1])/g' ./classes/flowequation.py:		
sed -i 's/(md,'fieldname','flowequation.augmented_lagrangian_r','numel',[1],'>',0.)/(md,fieldname='flowequation.augmented_lagrangian_r',numel=[1],gt=0.)/g' ./classes/flowequation.py:		
sed -i 's/(md,'fieldname','flowequation.augmented_lagrangian_rhop','numel',[1],'>',0.)/(md,fieldname='flowequation.augmented_lagrangian_rhop',numel=[1],gt=0.)/g' ./classes/flowequation.py:		
sed -i 's/(md,'fieldname','flowequation.augmented_lagrangian_rlambda','numel',[1],'>',0.)/(md,fieldname='flowequation.augmented_lagrangian_rlambda',numel=[1],gt=0.)/g' ./classes/flowequation.py:		
sed -i 's/(md,'fieldname','flowequation.augmented_lagrangian_rholambda','numel',[1],'>',0.)/(md,fieldname='flowequation.augmented_lagrangian_rholambda',numel=[1],gt=0.)/g' ./classes/flowequation.py:		
sed -i 's/(md,'fieldname','flowequation.XTH_theta','numel',[1],'>=',0.,'<',.5)/(md,fieldname='flowequation.XTH_theta',numel=[1],ge=0.,'<',.5)/g' ./classes/flowequation.py:		
sed -i 's/(md,'fieldname','flowequation.vertex_equation','size',[md.mesh.numberofvertices],'values',[1,2])/(md,fieldname='flowequation.vertex_equation',size=[md.mesh.numberofvertices],values=[1,2])/g' ./classes/flowequation.py:			
sed -i 's/(md,'fieldname','flowequation.element_equation','size',[md.mesh.numberofelements],'values',[1,2])/(md,fieldname='flowequation.element_equation',size=[md.mesh.numberofelements],'values',[1,2])/g' ./classes/flowequation.py:			
sed -i 's/(md,'fieldname','flowequation.vertex_equation','size',[md.mesh.numberofvertices],values=numpy.arange(0,8+1))/(md,fieldname='flowequation.vertex_equation',size=[md.mesh.numberofvertices],'values',numpy.arange(0,8+1))/g' ./classes/flowequation.py:			
sed -i 's/(md,'fieldname','flowequation.element_equation','size',[md.mesh.numberofelements],values=numpy.arange(0,8+1))/(md,fieldname='flowequation.element_equation',size=[md.mesh.numberofelements],'values',numpy.arange(0,8+1))/g' ./classes/flowequation.py:			
sed -i 's///g' ./classes/geometry.py:from checkfield import checkfield
sed -i 's/(md,'fieldname','geometry.surface'  ,'NaN',1,'size',[md.mesh.numberofvertices])/(md,fieldname='geometry.surface'  ,NaN=1,size=[md.mesh.numberofvertices])/g' ./classes/geometry.py:		
sed -i 's/(md,'fieldname','geometry.base'      ,'NaN',1,'size',[md.mesh.numberofvertices])/(md,fieldname='geometry.base'      ,NaN=1,size=[md.mesh.numberofvertices])/g' ./classes/geometry.py:		
sed -i 's/(md,'fieldname','geometry.thickness','NaN',1,'size',[md.mesh.numberofvertices],'>',0,'timeseries',1)/(md,fieldname='geometry.thickness',NaN=1,size=[md.mesh.numberofvertices],gt=0,timeseries=1)/g' ./classes/geometry.py:		
sed -i 's/(md,'fieldname','geometry.bed','NaN',1,'size',[md.mesh.numberofvertices])/(md,fieldname='geometry.bed',NaN=1,size=[md.mesh.numberofvertices])/g' ./classes/geometry.py:			
sed -i 's///g' ./classes/inversion.py:from checkfield import checkfield
sed -i 's/(md,'fieldname','inversion.iscontrol','values',[0,1])/(md,fieldname='inversion.iscontrol',values=[0,1])/g' ./classes/inversion.py:		
sed -i 's/(md,'fieldname','inversion.incomplete_adjoint','values',[0,1])/(md,fieldname='inversion.incomplete_adjoint',values=[0,1])/g' ./classes/inversion.py:		
sed -i 's/(md,'fieldname','inversion.control_parameters','cell',1,'values',supportedcontrols())/(md,fieldname='inversion.control_parameters','cell',1,values=supportedcontrols())/g' ./classes/inversion.py:		
sed -i 's/(md,'fieldname','inversion.nsteps','numel',[1],'>=',0)/(md,fieldname='inversion.nsteps',numel=[1],ge=0)/g' ./classes/inversion.py:		
sed -i 's/(md,'fieldname','inversion.maxiter_per_step','size',[md.inversion.nsteps],'>=',0)/(md,fieldname='inversion.maxiter_per_step',size=[md.inversion.nsteps],ge=0)/g' ./classes/inversion.py:		
sed -i 's/(md,'fieldname','inversion.step_threshold','size',[md.inversion.nsteps])/(md,fieldname='inversion.step_threshold',size=[md.inversion.nsteps])/g' ./classes/inversion.py:		
sed -i 's/(md,'fieldname','inversion.cost_functions','size',[num_costfunc],'values',supportedcostfunctions())/(md,fieldname='inversion.cost_functions',size=[num_costfunc],values=supportedcostfunctions())/g' ./classes/inversion.py:		
sed -i 's/(md,'fieldname','inversion.cost_functions_coefficients','size',[md.mesh.numberofvertices,num_costfunc],'>=',0)/(md,fieldname='inversion.cost_functions_coefficients',size=[md.mesh.numberofvertices,num_costfunc],ge=0)/g' ./classes/inversion.py:		
sed -i 's/(md,'fieldname','inversion.gradient_scaling','size',[md.inversion.nsteps,num_controls])/(md,fieldname='inversion.gradient_scaling',size=[md.inversion.nsteps,num_controls])/g' ./classes/inversion.py:		
sed -i 's/(md,'fieldname','inversion.min_parameters','size',[md.mesh.numberofvertices,num_controls])/(md,fieldname='inversion.min_parameters',size=[md.mesh.numberofvertices,num_controls])/g' ./classes/inversion.py:		
sed -i 's/(md,'fieldname','inversion.max_parameters','size',[md.mesh.numberofvertices,num_controls])/(md,fieldname='inversion.max_parameters',size=[md.mesh.numberofvertices,num_controls])/g' ./classes/inversion.py:		
sed -i 's/(md,'fieldname','inversion.thickness_obs','size',[md.mesh.numberofvertices],'NaN',1)/(md,fieldname='inversion.thickness_obs',size=[md.mesh.numberofvertices],NaN=1)/g' ./classes/inversion.py:			
sed -i 's/(md,'fieldname','inversion.vx_obs','size',[md.mesh.numberofvertices],'NaN',1)/(md,fieldname='inversion.vx_obs',size=[md.mesh.numberofvertices],NaN=1)/g' ./classes/inversion.py:			
sed -i 's/(md,'fieldname','inversion.vy_obs','size',[md.mesh.numberofvertices],'NaN',1)/(md,fieldname='inversion.vy_obs',size=[md.mesh.numberofvertices],NaN=1)/g' ./classes/inversion.py:			
sed -i 's///g' ./classes/qmu.py:from checkfield import checkfield
sed -i 's///g' ./classes/mesh3dprisms.py:from checkfield import *
sed -i 's/(md,'fieldname','mesh.x','NaN',1,'size',[md.mesh.numberofvertices])/(md,fieldname='mesh.x',NaN=1,size=[md.mesh.numberofvertices])/g' ./classes/mesh3dprisms.py:		
sed -i 's/(md,'fieldname','mesh.y','NaN',1,'size',[md.mesh.numberofvertices])/(md,fieldname='mesh.y',NaN=1,size=[md.mesh.numberofvertices])/g' ./classes/mesh3dprisms.py:		
sed -i 's/(md,'fieldname','mesh.z','NaN',1,'size',[md.mesh.numberofvertices])/(md,fieldname='mesh.z',NaN=1,size=[md.mesh.numberofvertices])/g' ./classes/mesh3dprisms.py:		
sed -i 's/(md,'fieldname','mesh.elements','NaN',1,'>',0,'values',numpy.arange(1,md.mesh.numberofvertices+1))/(md,fieldname='mesh.elements',NaN=1,gt=0,values=numpy.arange(1,md.mesh.numberofvertices+1))/g' ./classes/mesh3dprisms.py:		
sed -i 's/(md,'fieldname','mesh.elements','size',[md.mesh.numberofelements,6])/(md,fieldname='mesh.elements',size=[md.mesh.numberofelements,6])/g' ./classes/mesh3dprisms.py:		
sed -i 's/(md,'fieldname','mesh.numberoflayers','>=',0)/(md,fieldname='mesh.numberoflayers',ge=0)/g' ./classes/mesh3dprisms.py:		
sed -i 's/(md,'fieldname','mesh.numberofelements','>',0)/(md,fieldname='mesh.numberofelements',gt=0)/g' ./classes/mesh3dprisms.py:		
sed -i 's/(md,'fieldname','mesh.numberofvertices','>',0)/(md,fieldname='mesh.numberofvertices',gt=0)/g' ./classes/mesh3dprisms.py:		
sed -i 's/(md,'fieldname','mesh.vertexonbase','size',[md.mesh.numberofvertices],'values',[0,1])/(md,fieldname='mesh.vertexonbase',size=[md.mesh.numberofvertices],values=[0,1])/g' ./classes/mesh3dprisms.py:		
sed -i 's/(md,'fieldname','mesh.vertexonsurface','size',[md.mesh.numberofvertices],'values',[0,1])/(md,fieldname='mesh.vertexonsurface',size=[md.mesh.numberofvertices],values=[0,1])/g' ./classes/mesh3dprisms.py:		
sed -i 's/(md,'fieldname','mesh.average_vertex_connectivity','>=',24,'message',"'mesh.average_vertex_connectivity' should be at least 24 in 3d")/(md,fieldname='mesh.average_vertex_connectivity',ge=24,'message',"'mesh.average_vertex_connectivity' should be at least 24 in 3d")/g' ./classes/mesh3dprisms.py:		
sed -i 's///g' ./consistency/checkfield.py:def checkfield(md,**kwargs):

