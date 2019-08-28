#pragma rtGlobals=3		// Use modern global access method and strict wave access.

//**********************************************************************************
Function allatonce_FitCoa(pw,yw,xw) : FitFunc
	Wave pw,yw,xw

	//CurveFitDialog/ These comments were created by the Curve Fitting dialog. Altering them will
	//CurveFitDialog/ make the function less convenient to work with in the Curve Fitting dialog.
	//CurveFitDialog/ Equation:
	//CurveFitDialog/ f(p(w)) = SOM_v1
	//CurveFitDialog/ End of Equation
	//CurveFitDialog/ Independent Variables 1
	// CurveFitDialog/ pw
	//CurveFitDialog/ Coefficients 6
	//CurveFitDialog/ pw[0] = Frag
	//CurveFitDialog/ pw[1] = dLVP
	//CurveFitDialog/ pw[2] = Oxy1
	//CurveFitDialog/ pw[3] = Oxy2
	//CurveFitDialog/ pw[4] = Oxy3
	//CurveFitDialog/ pw[5] = Oxy4

	SOM_v1(fitcoefs=pw,fittime=xw,quiet=1)
	wave coa_time_interp
	yw = coa_time_interp
End

//**********************************************************************************
Function allatonce_FitLogOffsetCoa(pw,yw,xw) : FitFunc
	Wave pw,yw,xw

	//CurveFitDialog/ These comments were created by the Curve Fitting dialog. Altering them will
	//CurveFitDialog/ make the function less convenient to work with in the Curve Fitting dialog.
	//CurveFitDialog/ Equation:
	//CurveFitDialog/ f(p(w)) = SOM_v1
	//CurveFitDialog/ End of Equation
	//CurveFitDialog/ Independent Variables 1
	// CurveFitDialog/ pw
	//CurveFitDialog/ Coefficients 6
	//CurveFitDialog/ pw[0] = Frag
	//CurveFitDialog/ pw[1] = dLVP
	//CurveFitDialog/ pw[2] = Oxy1
	//CurveFitDialog/ pw[3] = Oxy2
	//CurveFitDialog/ pw[4] = Oxy3
	//CurveFitDialog/ pw[5] = Oxy4

	SOM_v1(fitcoefs=pw,fittime=xw,quiet=1)
	wave coa_time_interp
	yw = log(2+coa_time_interp)
End

//**********************************************************************************
// Fit the particle-deposition wall-loss corrected data
Function FitCoaWallCorr(pw,yw,xw) : FitFunc
	Wave pw,yw,xw

	//CurveFitDialog/ These comments were created by the Curve Fitting dialog. Altering them will
	//CurveFitDialog/ make the function less convenient to work with in the Curve Fitting dialog.
	//CurveFitDialog/ Equation:
	//CurveFitDialog/ f(p(w)) = SOM_v1
	//CurveFitDialog/ End of Equation
	//CurveFitDialog/ Independent Variables 1
	// CurveFitDialog/ pw
	//CurveFitDialog/ Coefficients 6
	//CurveFitDialog/ pw[0] = Frag
	//CurveFitDialog/ pw[1] = dLVP
	//CurveFitDialog/ pw[2] = Oxy1
	//CurveFitDialog/ pw[3] = Oxy2
	//CurveFitDialog/ pw[4] = Oxy3
	//CurveFitDialog/ pw[5] = Oxy4

	SOM_v1(fitcoefs=pw,fittime=xw,quiet=1)
	wave coa_time_WallCorr_interp
	yw = coa_time_WallCorr_interp
End

//**********************************************************************************
Function FitCoa_GECKO(pw,yw,xw) : FitFunc
	Wave pw,yw,xw

	//CurveFitDialog/ These comments were created by the Curve Fitting dialog. Altering them will
	//CurveFitDialog/ make the function less convenient to work with in the Curve Fitting dialog.
	// CurveFitDialog/ Equation:
	//CurveFitDialog/ f(p(w)) = SOM_v1
	//CurveFitDialog/ End of Equation
	//CurveFitDialog/ Independent Variables 1
	// CurveFitDialog/ pw
	//CurveFitDialog/ Coefficients 6
	//CurveFitDialog/ pw[0] = Frag
	//CurveFitDialog/ pw[1] = dLVP
	//CurveFitDialog/ pw[2] = Oxy1
	//CurveFitDialog/ pw[3] = Oxy2
	//CurveFitDialog/ pw[4] = Oxy3
	//CurveFitDialog/ pw[5] = Oxy4

	SOM_v1(fitcoefs=pw,fittime=xw,quiet=1)
	wave coa_time
	yw = coa_time
End



//**********************************************************************************
Function allatonce_FitO2C(pw,yw,xw) : FitFunc
	Wave pw,yw,xw

	//CurveFitDialog/ These comments were created by the Curve Fitting dialog. Altering them will
	//CurveFitDialog/ make the function less convenient to work with in the Curve Fitting dialog.
	//CurveFitDialog/ Coefficients 6
	//CurveFitDialog/ pw[0] = Frag
	//CurveFitDialog/ pw[1] = dLVP
	//CurveFitDialog/ pw[2] = Oxy1
	//CurveFitDialog/ pw[3] = Oxy2
	//CurveFitDialog/ pw[4] = Oxy3
	//CurveFitDialog/ pw[5] = Oxy4

	//RunSOM(xw,pw)
	//SOM_FitData(pw,xw)
	SOM_v1(fitcoefs=pw,fittime=xw,quiet=1)
	wave O2C_time_interp
	yw = O2C_time_interp
End

//**********************************************************************************
Function allatonce_FitAll(pw,yw,xw) : FitFunc
	Wave pw,yw,xw

	//CurveFitDialog/ These comments were created by the Curve Fitting dialog. Altering them will
	//CurveFitDialog/ make the function less convenient to work with in the Curve Fitting dialog.
	//CurveFitDialog/ Equation:
	//CurveFitDialog/ f(x) = funny stuff
	//CurveFitDialog/ End of Equation
	//CurveFitDialog/ Independent Variables 1
	//CurveFitDialog/ xw
	//CurveFitDialog/ Coefficients 6
	//CurveFitDialog/ pw[0] = Frag
	//CurveFitDialog/ pw[1] = dLVP
	//CurveFitDialog/ pw[2] = Oxy1
	//CurveFitDialog/ pw[3] = Oxy2
	//CurveFitDialog/ pw[4] = Oxy3
	//CurveFitDialog/ pw[5] = Oxy4
	//CurveFitDialog/ pw[6] = Coa_weight

	wave Coa_experiment, O2C_experiment, Coa_O2C_weighted
	Coa_O2C_weighted = Coa_experiment*pw[6] + O2C_experiment*(1-pw[6])
	
	SOM_v1(fitcoefs=pw,fittime=xw,quiet=1)
	wave coa_time
	wave O2C_time
	yw = coa_time*pw[6] + O2C_time*(1-pw[6])
End

//**********************************************************************************
Function allatonce_FitCxO1(pw,yw,xw) : FitFunc
	Wave pw,yw,xw

	//CurveFitDialog/ These comments were created by the Curve Fitting dialog. Altering them will
	//CurveFitDialog/ make the function less convenient to work with in the Curve Fitting dialog.
	// CurveFitDialog/ Equation:
	//CurveFitDialog/ f(p(w)) = SOM_v1
	//CurveFitDialog/ End of Equation
	//CurveFitDialog/ Independent Variables 1
	// CurveFitDialog/ pw
	//CurveFitDialog/ Coefficients 6
	//CurveFitDialog/ pw[0] = Frag
	//CurveFitDialog/ pw[1] = dLVP
	//CurveFitDialog/ pw[2] = Oxy1
	//CurveFitDialog/ pw[3] = Oxy2
	//CurveFitDialog/ pw[4] = Oxy3
	//CurveFitDialog/ pw[5] = Oxy4

	SOM_v1(fitcoefs=pw,fittime=xw,quiet=1)
	wave ParticleMass_time
	yw = ParticleMass_Time[0][1][p]
End

//
//**********************************************************************************
Function allatonce_FitCxO2(pw,yw,xw) : FitFunc
	Wave pw,yw,xw

	//CurveFitDialog/ These comments were created by the Curve Fitting dialog. Altering them will
	//CurveFitDialog/ make the function less convenient to work with in the Curve Fitting dialog.
	// CurveFitDialog/ Equation:
	//CurveFitDialog/ f(p(w)) = SOM_v1
	//CurveFitDialog/ End of Equation
	//CurveFitDialog/ Independent Variables 1
	// CurveFitDialog/ pw
	//CurveFitDialog/ Coefficients 6
	//CurveFitDialog/ pw[0] = Frag
	//CurveFitDialog/ pw[1] = dLVP
	//CurveFitDialog/ pw[2] = Oxy1
	//CurveFitDialog/ pw[3] = Oxy2
	//CurveFitDialog/ pw[4] = Oxy3
	//CurveFitDialog/ pw[5] = Oxy4

	SOM_v1(fitcoefs=pw,fittime=xw,quiet=1)
	wave ParticleMass_time
	yw = ParticleMass_Time[0][2][p]
End

//**********************************************************************************
Function allatonce_FitCxO3(pw,yw,xw) : FitFunc
	Wave pw,yw,xw

	//CurveFitDialog/ These comments were created by the Curve Fitting dialog. Altering them will
	//CurveFitDialog/ make the function less convenient to work with in the Curve Fitting dialog.
	// CurveFitDialog/ Equation:
	//CurveFitDialog/ f(p(w)) = SOM_v1
	//CurveFitDialog/ End of Equation
	//CurveFitDialog/ Independent Variables 1
	// CurveFitDialog/ pw
	//CurveFitDialog/ Coefficients 6
	//CurveFitDialog/ pw[0] = Frag
	//CurveFitDialog/ pw[1] = dLVP
	//CurveFitDialog/ pw[2] = Oxy1
	//CurveFitDialog/ pw[3] = Oxy2
	//CurveFitDialog/ pw[4] = Oxy3
	//CurveFitDialog/ pw[5] = Oxy4

	SOM_v1(fitcoefs=pw,fittime=xw,quiet=1)
	wave ParticleMass_time
	variable nLayers = dimsize(ParticleMass_time,2)
	make/o/d/n=(nLayers) myConc
	myConc = ParticleMass_time[0][3][p]
	yw = myConc
End

//**********************************************************************************
Function allatonce_FitCxO4(pw,yw,xw) : FitFunc
	Wave pw,yw,xw

	//CurveFitDialog/ These comments were created by the Curve Fitting dialog. Altering them will
	//CurveFitDialog/ make the function less convenient to work with in the Curve Fitting dialog.
	// CurveFitDialog/ Equation:
	//CurveFitDialog/ f(p(w)) = SOM_v1
	//CurveFitDialog/ End of Equation
	//CurveFitDialog/ Independent Variables 1
	// CurveFitDialog/ pw
	//CurveFitDialog/ Coefficients 6
	//CurveFitDialog/ pw[0] = Frag
	//CurveFitDialog/ pw[1] = dLVP
	//CurveFitDialog/ pw[2] = Oxy1
	//CurveFitDialog/ pw[3] = Oxy2
	//CurveFitDialog/ pw[4] = Oxy3
	//CurveFitDialog/ pw[5] = Oxy4

	SOM_v1(fitcoefs=pw,fittime=xw,quiet=1)
	wave ParticleMass_time
	yw = ParticleMass_Time[0][4][p]
End

//**********************************************************************************
Function allatonce_FitCxO5(pw,yw,xw) : FitFunc
	Wave pw,yw,xw

	//CurveFitDialog/ These comments were created by the Curve Fitting dialog. Altering them will
	//CurveFitDialog/ make the function less convenient to work with in the Curve Fitting dialog.
	// CurveFitDialog/ Equation:
	//CurveFitDialog/ f(p(w)) = SOM_v1
	//CurveFitDialog/ End of Equation
	//CurveFitDialog/ Independent Variables 1
	// CurveFitDialog/ pw
	//CurveFitDialog/ Coefficients 6
	//CurveFitDialog/ pw[0] = Frag
	//CurveFitDialog/ pw[1] = dLVP
	//CurveFitDialog/ pw[2] = Oxy1
	//CurveFitDialog/ pw[3] = Oxy2
	//CurveFitDialog/ pw[4] = Oxy3
	//CurveFitDialog/ pw[5] = Oxy4

	SOM_v1(fitcoefs=pw,fittime=xw,quiet=1)
	wave ParticleMass_time
	yw = ParticleMass_Time[0][5][p]
End

//**********************************************************************************
Function allatonce_FitCxO6(pw,yw,xw) : FitFunc
	Wave pw,yw,xw

	//CurveFitDialog/ These comments were created by the Curve Fitting dialog. Altering them will
	//CurveFitDialog/ make the function less convenient to work with in the Curve Fitting dialog.
	// CurveFitDialog/ Equation:
	//CurveFitDialog/ f(p(w)) = SOM_v1
	//CurveFitDialog/ End of Equation
	//CurveFitDialog/ Independent Variables 1
	// CurveFitDialog/ pw
	//CurveFitDialog/ Coefficients 6
	//CurveFitDialog/ pw[0] = Frag
	//CurveFitDialog/ pw[1] = dLVP
	//CurveFitDialog/ pw[2] = Oxy1
	//CurveFitDialog/ pw[3] = Oxy2
	//CurveFitDialog/ pw[4] = Oxy3
	//CurveFitDialog/ pw[5] = Oxy4

	SOM_v1(fitcoefs=pw,fittime=xw,quiet=1)
	wave ParticleMass_time
	yw = ParticleMass_Time[0][6][p]
End

//**********************************************************************************
Function allatonce_FitCoa_kw_alpha(pw,yw,xw) : FitFunc
	Wave pw,yw,xw

	//CurveFitDialog/ These comments were created by the Curve Fitting dialog. Altering them will
	//CurveFitDialog/ make the function less convenient to work with in the Curve Fitting dialog.
	//CurveFitDialog/ Equation:
	//CurveFitDialog/ f(p(w)) = SOM_v1
	//CurveFitDialog/ End of Equation
	//CurveFitDialog/ Independent Variables 1
	// CurveFitDialog/ pw
	//CurveFitDialog/ Coefficients 8
	//CurveFitDialog/ pw[0] = Frag
	//CurveFitDialog/ pw[1] = dLVP
	//CurveFitDialog/ pw[2] = Oxy1
	//CurveFitDialog/ pw[3] = Oxy2
	//CurveFitDialog/ pw[4] = Oxy3
	//CurveFitDialog/ pw[5] = Oxy4
	//CurveFitDialog/ pw[6] = kwall
	//CurveFitDialog/ pw[7] = alpha
	
	make/o/d/n=(6) root:pw_short
	wave pw_short = root:pw_short
	pw_short = pw[x]
	NVAR WLRgas = root:gasWLR
	WLRgas = pw[6]
	NVAR alpha = root:alpha
	alpha = pw[7]
	
	SOM_v1(fitcoefs=pw_short,fittime=xw,quiet=1)
	wave coa_time_interp
	yw = coa_time_interp
End

//**********************************************************************************
Function FitCoa_kw_alpha_low1(pw,yw,xw) : FitFunc
	Wave pw,yw,xw

	//CurveFitDialog/ These comments were created by the Curve Fitting dialog. Altering them will
	//CurveFitDialog/ make the function less convenient to work with in the Curve Fitting dialog.
	//CurveFitDialog/ Equation:
	//CurveFitDialog/ f(p(w)) = SOM_v1
	//CurveFitDialog/ End of Equation
	//CurveFitDialog/ Independent Variables 1
	// CurveFitDialog/ pw
	//CurveFitDialog/ Coefficients 8
	//CurveFitDialog/ pw[0] = Frag
	//CurveFitDialog/ pw[1] = dLVP
	//CurveFitDialog/ pw[2] = Oxy1
	//CurveFitDialog/ pw[3] = Oxy2
	//CurveFitDialog/ pw[4] = Oxy3
	//CurveFitDialog/ pw[5] = Oxy4
	//CurveFitDialog/ pw[6] = kwall
	//CurveFitDialog/ pw[7] = alpha
	
	make/o/d/n=(6) root:pw_short
	wave pw_short = root:pw_short
	pw_short = pw[x]
	NVAR WLRgas = root:gasWLR
	WLRgas = pw[6]
	NVAR alpha = root:alpha
	alpha = pw[7]
	setfromcaltech("Toluene2013","low1")
	
	SOM_v1(fitcoefs=pw_short,fittime=xw,quiet=1)
	wave coa_time_interp
	yw = coa_time_interp
End

//**********************************************************************************
Function FitCoa_kw_alpha_low2(pw,yw,xw) : FitFunc
	Wave pw,yw,xw

	//CurveFitDialog/ These comments were created by the Curve Fitting dialog. Altering them will
	//CurveFitDialog/ make the function less convenient to work with in the Curve Fitting dialog.
	//CurveFitDialog/ Equation:
	//CurveFitDialog/ f(p(w)) = SOM_v1
	//CurveFitDialog/ End of Equation
	//CurveFitDialog/ Independent Variables 1
	// CurveFitDialog/ pw
	//CurveFitDialog/ Coefficients 8
	//CurveFitDialog/ pw[0] = Frag
	//CurveFitDialog/ pw[1] = dLVP
	//CurveFitDialog/ pw[2] = Oxy1
	//CurveFitDialog/ pw[3] = Oxy2
	//CurveFitDialog/ pw[4] = Oxy3
	//CurveFitDialog/ pw[5] = Oxy4
	//CurveFitDialog/ pw[6] = kwall
	//CurveFitDialog/ pw[7] = alpha
	
	make/o/d/n=(6) root:pw_short
	wave pw_short = root:pw_short
	pw_short = pw[x]
	NVAR WLRgas = root:gasWLR
	WLRgas = pw[6]
	NVAR alpha = root:alpha
	alpha = pw[7]
	setfromcaltech("Toluene2013","low2")
	
	SOM_v1(fitcoefs=pw_short,fittime=xw,quiet=1)
	wave coa_time_interp
	yw = coa_time_interp
End

//**********************************************************************************
Function FitCoa_kw_alpha_low3(pw,yw,xw) : FitFunc
	Wave pw,yw,xw

	//CurveFitDialog/ These comments were created by the Curve Fitting dialog. Altering them will
	//CurveFitDialog/ make the function less convenient to work with in the Curve Fitting dialog.
	//CurveFitDialog/ Equation:
	//CurveFitDialog/ f(p(w)) = SOM_v1
	//CurveFitDialog/ End of Equation
	//CurveFitDialog/ Independent Variables 1
	// CurveFitDialog/ pw
	//CurveFitDialog/ Coefficients 8
	//CurveFitDialog/ pw[0] = Frag
	//CurveFitDialog/ pw[1] = dLVP
	//CurveFitDialog/ pw[2] = Oxy1
	//CurveFitDialog/ pw[3] = Oxy2
	//CurveFitDialog/ pw[4] = Oxy3
	//CurveFitDialog/ pw[5] = Oxy4
	//CurveFitDialog/ pw[6] = kwall
	//CurveFitDialog/ pw[7] = alpha
	
	make/o/d/n=(6) root:pw_short
	wave pw_short = root:pw_short
	pw_short = pw[x]
	NVAR WLRgas = root:gasWLR
	WLRgas = pw[6]
	NVAR alpha = root:alpha
	alpha = pw[7]
	setfromcaltech("Toluene2013","low3")
	
	SOM_v1(fitcoefs=pw_short,fittime=xw,quiet=1)
	wave coa_time_interp
	yw = coa_time_interp
End

//**********************************************************************************
Function FitCoa_kw_alpha_low4(pw,yw,xw) : FitFunc
	Wave pw,yw,xw

	//CurveFitDialog/ These comments were created by the Curve Fitting dialog. Altering them will
	//CurveFitDialog/ make the function less convenient to work with in the Curve Fitting dialog.
	//CurveFitDialog/ Equation:
	//CurveFitDialog/ f(p(w)) = SOM_v1
	//CurveFitDialog/ End of Equation
	//CurveFitDialog/ Independent Variables 1
	// CurveFitDialog/ pw
	//CurveFitDialog/ Coefficients 8
	//CurveFitDialog/ pw[0] = Frag
	//CurveFitDialog/ pw[1] = dLVP
	//CurveFitDialog/ pw[2] = Oxy1
	//CurveFitDialog/ pw[3] = Oxy2
	//CurveFitDialog/ pw[4] = Oxy3
	//CurveFitDialog/ pw[5] = Oxy4
	//CurveFitDialog/ pw[6] = kwall
	//CurveFitDialog/ pw[7] = alpha
	
	make/o/d/n=(6) root:pw_short
	wave pw_short = root:pw_short
	pw_short = pw[x]
	NVAR WLRgas = root:gasWLR
	WLRgas = pw[6]
	NVAR alpha = root:alpha
	alpha = pw[7]
	setfromcaltech("Toluene2013","low4")
	
	SOM_v1(fitcoefs=pw_short,fittime=xw,quiet=1)
	wave coa_time_interp
	yw = coa_time_interp
End

//**********************************************************************************
Function FitCoa_kw_alpha_low5(pw,yw,xw) : FitFunc
	Wave pw,yw,xw

	//CurveFitDialog/ These comments were created by the Curve Fitting dialog. Altering them will
	//CurveFitDialog/ make the function less convenient to work with in the Curve Fitting dialog.
	//CurveFitDialog/ Equation:
	//CurveFitDialog/ f(p(w)) = SOM_v1
	//CurveFitDialog/ End of Equation
	//CurveFitDialog/ Independent Variables 1
	// CurveFitDialog/ pw
	//CurveFitDialog/ Coefficients 8
	//CurveFitDialog/ pw[0] = Frag
	//CurveFitDialog/ pw[1] = dLVP
	//CurveFitDialog/ pw[2] = Oxy1
	//CurveFitDialog/ pw[3] = Oxy2
	//CurveFitDialog/ pw[4] = Oxy3
	//CurveFitDialog/ pw[5] = Oxy4
	//CurveFitDialog/ pw[6] = kwall
	//CurveFitDialog/ pw[7] = alpha
	
	make/o/d/n=(6) root:pw_short
	wave pw_short = root:pw_short
	pw_short = pw[x]
	NVAR WLRgas = root:gasWLR
	WLRgas = pw[6]
	NVAR alpha = root:alpha
	alpha = pw[7]
	setfromcaltech("Toluene2013","low5")
	
	SOM_v1(fitcoefs=pw_short,fittime=xw,quiet=1)
	wave coa_time_interp
	yw = coa_time_interp
End