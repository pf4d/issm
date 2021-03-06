import numpy

def cuffey(temperature):
	"""
	CUFFEY - calculates ice rigidity as a function of temperature

	   rigidity (in s^(1/3)Pa) is the flow law parameter in the flow law sigma=B*e(1/3)
		(Cuffey and Paterson, p75). 
	   temperature is in Kelvin degrees

	   Usage:
	      rigidity=cuffey(temperature)
	"""
	
	if numpy.any(temperature<0.):
		raise RuntimeError("input temperature should be in Kelvin (positive)")
	
	T = temperature.reshape(-1,)-273.15

	rigidity=numpy.zeros_like(T)
	pos=numpy.nonzero(T<=-45)
	rigidity[pos]=10**8*(-0.000396645116301*(T[pos]+50)**3+ 0.013345579471334*(T[pos]+50)**2  -0.356868703259105*(T[pos]+50)+7.272363035371383)
	pos=numpy.nonzero(numpy.logical_and(-45<=T,T<-40))
	rigidity[pos]=10**8*(-0.000396645116301*(T[pos]+45)**3+ 0.007395902726819*(T[pos]+45)**2  -0.253161292268336*(T[pos]+45)+5.772078366321591)
	pos=numpy.nonzero(numpy.logical_and(-40<=T,T<-35))
	rigidity[pos]=10**8*(0.000408322072669*(T[pos]+40)**3+  0.001446225982305*(T[pos]+40)**2  -0.208950648722716*(T[pos]+40)+4.641588833612773)
	pos=numpy.nonzero(numpy.logical_and(-35<=T,T<-30))
	rigidity[pos]=10**8*(-0.000423888728124*(T[pos]+35)**3+ 0.007571057072334*(T[pos]+35)**2  -0.163864233449525*(T[pos]+35)+3.684031498640382)
	pos=numpy.nonzero(numpy.logical_and(-30<=T,T<-25))
	rigidity[pos]=10**8*(0.000147154327025*(T[pos]+30)**3+ 0.001212726150476*(T[pos]+30)**2  -0.119945317335478*(T[pos]+30)+3.001000667185614)
	pos=numpy.nonzero(numpy.logical_and(-25<=T,T<-20))
	rigidity[pos]=10**8*(-0.000193435838672*(T[pos]+25)**3+ 0.003420041055847*(T[pos]+25)**2  -0.096781481303861*(T[pos]+25)+2.449986525148220)
	pos=numpy.nonzero(numpy.logical_and(-20<=T,T<-15))
	rigidity[pos]=10**8*(0.000219771255067*(T[pos]+20)**3+  0.000518503475772*(T[pos]+20)**2  -0.077088758645767*(T[pos]+20)+2.027400665191131)
	pos=numpy.nonzero(numpy.logical_and(-15<=T,T<-10))
	rigidity[pos]=10**8*(-0.000653438900191*(T[pos]+15)**3+ 0.003815072301777*(T[pos]+15)**2  -0.055420879758021*(T[pos]+15)+1.682390865739973)
	pos=numpy.nonzero(numpy.logical_and(-10<=T,T<-5))
	rigidity[pos]=10**8*(0.000692439419762*(T[pos]+10)**3 -0.005986511201093 *(T[pos]+10)**2 -0.066278074254598*(T[pos]+10)+1.418983411970382)
	pos=numpy.nonzero(numpy.logical_and(-5<=T,T<-2))
	rigidity[pos]=10**8*(-0.000132282004110*(T[pos]+5)**3 +0.004400080095332*(T[pos]+5)**2    -0.074210229783403*(T[pos]+5)+ 1.024485188140279)
	pos=numpy.nonzero(-2<=T)
	rigidity[pos]=10**8*(-0.000132282004110*(T[pos]+2)**3 +0.003209542058346*(T[pos]+2)**2    -0.051381363322371*(T[pos]+2)+ 0.837883605537096)

	#Now make sure that rigidity is positive
	pos=numpy.nonzero(rigidity<0)
	rigidity[pos]=1**6 

	return rigidity

