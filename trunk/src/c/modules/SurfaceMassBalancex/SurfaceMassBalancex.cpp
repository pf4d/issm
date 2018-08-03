/*!\file SurfaceMassBalancex
 * \brief: calculates SMB 
 */

#include "./SurfaceMassBalancex.h"
#include "../../shared/shared.h"
#include "../../toolkits/toolkits.h"


void SmbGradientsx(FemModel* femmodel){/*{{{*/

	// void SurfaceMassBalancex(hd,agd,ni){
	//    INPUT parameters: ni: working size of arrays
	//    INPUT: surface elevation (m): hd(NA)
	//    OUTPUT: mass-balance (m/yr ice): agd(NA)
	int v;
	IssmDouble rho_water;                   // density of fresh water
	IssmDouble rho_ice;                     // density of ice
	IssmDouble yts;								// conversion factor year to second

	/*Loop over all the elements of this partition*/
	for(int i=0;i<femmodel->elements->Size();i++){
		Element* element=xDynamicCast<Element*>(femmodel->elements->GetObjectByOffset(i));

		/*Allocate all arrays*/
		int         numvertices = element->GetNumberOfVertices();
		IssmDouble* Href        = xNew<IssmDouble>(numvertices); // reference elevation from which deviations are used to calculate the SMB adjustment
		IssmDouble* Smbref      = xNew<IssmDouble>(numvertices); // reference SMB to which deviations are added
		IssmDouble* b_pos       = xNew<IssmDouble>(numvertices); // Hs-SMB relation parameter
		IssmDouble* b_neg       = xNew<IssmDouble>(numvertices); // Hs-SMB relation paremeter
		IssmDouble* s           = xNew<IssmDouble>(numvertices); // surface elevation (m)
		IssmDouble* smb         = xNew<IssmDouble>(numvertices);

		/*Recover SmbGradients*/
		element->GetInputListOnVertices(Href,SmbHrefEnum);
		element->GetInputListOnVertices(Smbref,SmbSmbrefEnum);
		element->GetInputListOnVertices(b_pos,SmbBPosEnum);
		element->GetInputListOnVertices(b_neg,SmbBNegEnum);

		/*Recover surface elevation at vertices: */
		element->GetInputListOnVertices(s,SurfaceEnum);

		/*Get material parameters :*/
		rho_ice=element->matpar->GetMaterialParameter(MaterialsRhoIceEnum);
		rho_water=element->matpar->GetMaterialParameter(MaterialsRhoFreshwaterEnum);

		/* Get constants */
		femmodel->parameters->FindParam(&yts,ConstantsYtsEnum);

		// loop over all vertices
		for(v=0;v<numvertices;v++){
			if(Smbref[v]>0){
				smb[v]=Smbref[v]+b_pos[v]*(s[v]-Href[v]);
			}
			else{
				smb[v]=Smbref[v]+b_neg[v]*(s[v]-Href[v]);
			}

			smb[v]=smb[v]/1000*rho_water/rho_ice;      // SMB in m/y ice
		}  //end of the loop over the vertices

		/*Add input to element and Free memory*/
		element->AddInput(SmbMassBalanceEnum,smb,P1Enum);
		xDelete<IssmDouble>(Href);
		xDelete<IssmDouble>(Smbref);
		xDelete<IssmDouble>(b_pos);
		xDelete<IssmDouble>(b_neg);
		xDelete<IssmDouble>(s);
		xDelete<IssmDouble>(smb);
	}

}/*}}}*/
void SmbGradientsElax(FemModel* femmodel){/*{{{*/

	// void SurfaceMassBalancex(hd,agd,ni){
	//    INPUT parameters: ni: working size of arrays
	//    INPUT: surface elevation (m): hd(NA)
	//    OUTPUT: surface mass-balance (m/yr ice): agd(NA)
	int v;

	/*Loop over all the elements of this partition*/
	for(int i=0;i<femmodel->elements->Size();i++){
		Element* element=xDynamicCast<Element*>(femmodel->elements->GetObjectByOffset(i));

		/*Allocate all arrays*/
		int         numvertices = element->GetNumberOfVertices();
		IssmDouble* ela       = xNew<IssmDouble>(numvertices); // Equilibrium Line Altitude (m a.s.l) to which deviations are used to calculate the SMB
		IssmDouble* b_pos       = xNew<IssmDouble>(numvertices); // SMB gradient above ELA (m ice eq. per m elevation change)
		IssmDouble* b_neg       = xNew<IssmDouble>(numvertices); // SMB gradient below ELA (m ice eq. per m elevation change)
		IssmDouble* b_max       = xNew<IssmDouble>(numvertices); // Upper cap on SMB rate (m/y ice eq.)
		IssmDouble* b_min       = xNew<IssmDouble>(numvertices); // Lower cap on SMB rate (m/y ice eq.)
		IssmDouble* s           = xNew<IssmDouble>(numvertices); // Surface elevation (m a.s.l.)
		IssmDouble* smb         = xNew<IssmDouble>(numvertices); // SMB (m/y ice eq.)

		/*Recover ELA, SMB gradients, and caps*/
		element->GetInputListOnVertices(ela,SmbElaEnum);
		element->GetInputListOnVertices(b_pos,SmbBPosEnum);
		element->GetInputListOnVertices(b_neg,SmbBNegEnum);
		element->GetInputListOnVertices(b_max,SmbBMaxEnum);
		element->GetInputListOnVertices(b_min,SmbBMinEnum);

		/*Recover surface elevation at vertices: */
		element->GetInputListOnVertices(s,SurfaceEnum);

		/*Loop over all vertices, calculate SMB*/
		for(v=0;v<numvertices;v++){
			// if surface is above the ELA
			if(s[v]>ela[v]){		
				smb[v]=b_pos[v]*(s[v]-ela[v]);
			}
			// if surface is below or equal to the ELA
			else{
				smb[v]=b_neg[v]*(s[v]-ela[v]);
			}

			// if SMB is larger than upper cap, set SMB to upper cap
			if(smb[v]>b_max[v]){
				smb[v]=b_max[v];
			}
			// if SMB is smaller than lower cap, set SMB to lower cap
			if(smb[v]<b_min[v]){
				smb[v]=b_min[v];
			}
		}  //end of the loop over the vertices

		/*Add input to element and Free memory*/
		element->AddInput(SmbMassBalanceEnum,smb,P1Enum);
		xDelete<IssmDouble>(ela);
		xDelete<IssmDouble>(b_pos);
		xDelete<IssmDouble>(b_neg);
		xDelete<IssmDouble>(b_max);
		xDelete<IssmDouble>(b_min);
		xDelete<IssmDouble>(s);
		xDelete<IssmDouble>(smb);

	}

}/*}}}*/
void Delta18oParameterizationx(FemModel* femmodel){/*{{{*/

	for(int i=0;i<femmodel->elements->Size();i++){
		Element* element=xDynamicCast<Element*>(femmodel->elements->GetObjectByOffset(i));
		element->Delta18oParameterization();
	}

}/*}}}*/
void MungsmtpParameterizationx(FemModel* femmodel){/*{{{*/

	for(int i=0;i<femmodel->elements->Size();i++){
		Element* element=xDynamicCast<Element*>(femmodel->elements->GetObjectByOffset(i));
		element->MungsmtpParameterization();
	}

}/*}}}*/
void Delta18opdParameterizationx(FemModel* femmodel){/*{{{*/

	for(int i=0;i<femmodel->elements->Size();i++){
		Element* element=xDynamicCast<Element*>(femmodel->elements->GetObjectByOffset(i));
		element->Delta18opdParameterization();
	}

}/*}}}*/
void PositiveDegreeDayx(FemModel* femmodel){/*{{{*/

	// void PositiveDegreeDayx(hd,vTempsea,vPrec,agd,Tsurf,ni){
	//    note "v" prefix means 12 monthly means, ie time dimension
	//    INPUT parameters: ni: working size of arrays
	//    INPUT: surface elevation (m): hd(NA)
	//    monthly mean surface sealevel temperature (degrees C): vTempsea(NA
	//    ,NTIME) 
	//    monthly mean precip rate (m/yr water equivalent): vPrec(NA,NTIME)
	//    OUTPUT: mass-balance (m/yr ice): agd(NA)
	//    mean annual surface temperature (degrees C): Tsurf(NA)

	int    i, it, jj, itm;
	IssmDouble DT = 0.02, sigfac, snormfac;
	IssmDouble signorm = 5.5;      // signorm : sigma of the temperature distribution for a normal day 
	IssmDouble siglim;       // sigma limit for the integration which is equal to 2.5 sigmanorm
	IssmDouble signormc = signorm - 0.5;     // sigma of the temperature distribution for cloudy day
	IssmDouble siglimc, siglim0, siglim0c;
	IssmDouble tstep, tsint, tint, tstepc;
	int    NPDMAX = 1504, NPDCMAX = 1454;
	//IssmDouble pdds[NPDMAX]={0}; 
	//IssmDouble pds[NPDCMAX]={0};
	IssmDouble pddt, pd ; // pd : snow/precip fraction, precipitation falling as snow
	IssmDouble PDup, PDCUT = 2.0;    // PDcut: rain/snow cutoff temperature (C)
	IssmDouble tstar; // monthly mean surface temp

	bool ismungsm;

	IssmDouble *pdds    = NULL;
	IssmDouble *pds     = NULL;
	Element    *element = NULL;

	pdds=xNew<IssmDouble>(NPDMAX+1); 
	pds=xNew<IssmDouble>(NPDCMAX+1); 

	// Get ismungsm parameter
	femmodel->parameters->FindParam(&ismungsm,SmbIsmungsmEnum);

	/* initialize PDD (creation of a lookup table)*/
	tstep    = 0.1;
	tsint    = tstep*0.5;
	sigfac   = -1.0/(2.0*pow(signorm,2));
	snormfac = 1.0/(signorm*sqrt(2.0*acos(-1.0)));
	siglim   = 2.5*signorm;
	siglimc  = 2.5*signormc;
	siglim0  = siglim/DT + 0.5;
	siglim0c = siglimc/DT + 0.5;
	PDup     = siglimc+PDCUT;

	itm = reCast<int,IssmDouble>((2*siglim/DT + 1.5));

	if(itm >= NPDMAX) _error_("increase NPDMAX in massBalance.cpp");
	for(it = 0; it < itm; it++){  
		//    tstar = REAL(it)*DT-siglim;
		tstar = it*DT-siglim;
		tint = tsint;
		pddt = 0.;
		for ( jj = 0; jj < 600; jj++){
			if (tint > (tstar+siglim)){break;}
			pddt = pddt + tint*exp(sigfac*(pow((tint-tstar),2)))*tstep;
			tint = tint+tstep;
		}
		pdds[it] = pddt*snormfac;
	}
	pdds[itm+1] = siglim + DT;

	//*********compute PD(T) : snow/precip fraction. precipitation falling as snow
	tstepc   = 0.1;
	tsint    = PDCUT-tstepc*0.5;
	signormc = signorm - 0.5;
	sigfac   = -1.0/(2.0*pow(signormc,2));
	snormfac = 1.0/(signormc*sqrt(2.0*acos(-1.0)));
	siglimc  = 2.5*signormc ;
	itm = reCast<int,IssmDouble>((PDCUT+2.*siglimc)/DT + 1.5);
	if(itm >= NPDCMAX) _error_("increase NPDCMAX in p35com");
	for(it = 0; it < itm; it++ ){
		tstar = it*DT-siglimc;
		//    tstar = REAL(it)*DT-siglimc;
		tint = tsint;          // start against upper bound
		pd = 0.;
		for (jj = 0; jj < 600; jj++){
			if (tint<(tstar-siglimc)) {break;}
			pd = pd + exp(sigfac*(pow((tint-tstar),2)))*tstepc;
			tint = tint-tstepc;
		}
		pds[it] = pd*snormfac;  // gaussian integral lookup table for snow fraction
	}
	pds[itm+1] = 0.;
	//     *******END initialize PDD

	for(i=0;i<femmodel->elements->Size();i++){
		element=xDynamicCast<Element*>(femmodel->elements->GetObjectByOffset(i));
		element->PositiveDegreeDay(pdds,pds,signorm,ismungsm);
	}
	/*free ressouces: */
	xDelete<IssmDouble>(pdds);
	xDelete<IssmDouble>(pds);
}/*}}}*/
void SmbHenningx(FemModel* femmodel){/*{{{*/

	/*Intermediaries*/
	IssmDouble  z_critical = 1675.;
	IssmDouble  dz = 0;
	IssmDouble  a = -15.86;
	IssmDouble  b = 0.00969;
	IssmDouble  c = -0.235;
	IssmDouble  f = 1.;
	IssmDouble  g = -0.0011;
	IssmDouble  h = -1.54e-5;
	IssmDouble  smb,smbref,anomaly,yts,z;
    
    /* Get constants */
    femmodel->parameters->FindParam(&yts,ConstantsYtsEnum);
    /*iomodel->FindConstant(&yts,"md.constants.yts");*/
    /*this->parameters->FindParam(&yts,ConstantsYtsEnum);*/
    /*Mathieu original*/
    /*IssmDouble  smb,smbref,z;*/
    
	/*Loop over all the elements of this partition*/
	for(int i=0;i<femmodel->elements->Size();i++){
		Element* element=xDynamicCast<Element*>(femmodel->elements->GetObjectByOffset(i));

		/*Get reference SMB (uncorrected) and allocate all arrays*/
		int         numvertices = element->GetNumberOfVertices();
		IssmDouble* surfacelist = xNew<IssmDouble>(numvertices);
		IssmDouble* smblistref  = xNew<IssmDouble>(numvertices);
		IssmDouble* smblist     = xNew<IssmDouble>(numvertices);
		element->GetInputListOnVertices(surfacelist,SurfaceEnum);
		element->GetInputListOnVertices(smblistref,SmbSmbrefEnum);

		/*Loop over all vertices of element and correct SMB as a function of altitude z*/
		for(int v=0;v<numvertices;v++){

			/*Get vertex elevation, anoma smb*/
			z      = surfacelist[v];
			anomaly = smblistref[v];

            /* Henning edited acc. to Riannes equations*/
            /* Set SMB maximum elevation, if dz = 0 -> z_critical = 1675 */
            z_critical = z_critical + dz;
            
            /* Calculate smb acc. to the surface elevation z */
            if(z<z_critical){
				smb = a + b*z + c;
			}
			else{
			  smb = (a + b*z)*(f + g*(z-z_critical) + h*(z-z_critical)*(z-z_critical)) + c;
			}
            
			/* Compute smb including anomaly,
				correct for number of seconds in a year [s/yr]*/
			smb = smb/yts + anomaly;


			/*Update array accordingly*/
			smblist[v] = smb;

		}

		/*Add input to element and Free memory*/
		element->AddInput(SmbMassBalanceEnum,smblist,P1Enum);
		xDelete<IssmDouble>(surfacelist);
		xDelete<IssmDouble>(smblistref);
		xDelete<IssmDouble>(smblist);
	}

}/*}}}*/
void SmbComponentsx(FemModel* femmodel){/*{{{*/

	// void SmbComponentsx(acc,evap,runoff,ni){
	//    INPUT parameters: ni: working size of arrays
	//    INPUT: surface accumulation (m/yr water equivalent): acc
	//    surface evaporation (m/yr water equivalent): evap
	//    surface runoff (m/yr water equivalent): runoff
	//    OUTPUT: mass-balance (m/yr ice): agd(NA)
	int v;

	/*Loop over all the elements of this partition*/
	for(int i=0;i<femmodel->elements->Size();i++){
		Element* element=xDynamicCast<Element*>(femmodel->elements->GetObjectByOffset(i));

		/*Allocate all arrays*/
		int         numvertices = element->GetNumberOfVertices();
		IssmDouble* acc         = xNew<IssmDouble>(numvertices); 
		IssmDouble* evap        = xNew<IssmDouble>(numvertices);
		IssmDouble* runoff      = xNew<IssmDouble>(numvertices); 
		IssmDouble* smb         = xNew<IssmDouble>(numvertices);

		/*Recover Smb Components*/
		element->GetInputListOnVertices(acc,SmbAccumulationEnum);
		element->GetInputListOnVertices(evap,SmbEvaporationEnum);
		element->GetInputListOnVertices(runoff,SmbRunoffEnum);

		// loop over all vertices
		for(v=0;v<numvertices;v++){
			smb[v]=acc[v]-evap[v]-runoff[v];
		}  //end of the loop over the vertices

		/*Add input to element and Free memory*/
		element->AddInput(SmbMassBalanceEnum,smb,P1Enum);
		xDelete<IssmDouble>(acc);
		xDelete<IssmDouble>(evap);
		xDelete<IssmDouble>(runoff);
		xDelete<IssmDouble>(smb);
	}

}/*}}}*/
void SmbMeltComponentsx(FemModel* femmodel){/*{{{*/

	// void SmbMeltComponentsx(acc,evap,melt,refreeze,ni){
	//    INPUT parameters: ni: working size of arrays
	//    INPUT: surface accumulation (m/yr water equivalent): acc
	//    surface evaporation (m/yr water equivalent): evap
	//    surface melt (m/yr water equivalent): melt
	//    refreeze of surface melt (m/yr water equivalent): refreeze
	//    OUTPUT: mass-balance (m/yr ice): agd(NA)
	int v;

	/*Loop over all the elements of this partition*/
	for(int i=0;i<femmodel->elements->Size();i++){
		Element* element=xDynamicCast<Element*>(femmodel->elements->GetObjectByOffset(i));

		/*Allocate all arrays*/
		int         numvertices = element->GetNumberOfVertices();
		IssmDouble* acc         = xNew<IssmDouble>(numvertices);
		IssmDouble* evap        = xNew<IssmDouble>(numvertices); 
		IssmDouble* melt        = xNew<IssmDouble>(numvertices);
		IssmDouble* refreeze    = xNew<IssmDouble>(numvertices);
		IssmDouble* smb         = xNew<IssmDouble>(numvertices);

		/*Recover Smb Components*/
		element->GetInputListOnVertices(acc,SmbAccumulationEnum);
		element->GetInputListOnVertices(evap,SmbEvaporationEnum);
		element->GetInputListOnVertices(melt,SmbMeltEnum);
		element->GetInputListOnVertices(refreeze,SmbRefreezeEnum);

		// loop over all vertices
		for(v=0;v<numvertices;v++){
			smb[v]=acc[v]-evap[v]-melt[v]+refreeze[v];
		}  //end of the loop over the vertices

		/*Add input to element and Free memory*/
		element->AddInput(SmbMassBalanceEnum,smb,P1Enum);
		xDelete<IssmDouble>(acc);
		xDelete<IssmDouble>(evap);
		xDelete<IssmDouble>(melt);
		xDelete<IssmDouble>(refreeze);
		xDelete<IssmDouble>(smb);
	}

}/*}}}*/
