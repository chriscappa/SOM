#pragma rtGlobals=1		// Use modern global access method and strict wave access.

strconstant ksDFroot= "root:"
strconstant ksDFVOCbase = "root:VOCs"

//******************************************************************************************
// The Statistical Oxidation Model
// Actually started keeping notes on updates on 07/07/2013
// 10/16/2013: 1. Saved as v6 (from v5)
//10/13/2013: Added capability of turning off equilibrium partitioning to particles so that only gas-phase chemistry and gas wall-loss can be considered
// 08/20/2013: 1. fixed gas-phase wall-loss...it was in the wrong spot and didn't actually do anything. Values around 1e-6/s affect the results. 
//			  2. Added particle wall loss option, based on wall-losses in Loza et al. Atmos. Chem. Phys., 12, 151–167, 2012
// 07/20/2013: Added option to perform eqm. calculations based on mass concentrations, as opposed to molecular concentrations
//			it is not clear which is "better" to use, but some (e.g. Donahue) prefer mass, with the argument that it better captures the actual physicallity of the situation, which is that bigger molecules take up more space.
//			search for string "units" to update and select which method to use
// 07/11/2013: MAJOR update
//			1. it was determined that there was an error in how the fragmentation was being applied. this has now been corrected so that
//			fragmentation probabilities depend on the oxygen content of the PARENT species. Oops. Corrected both in SOM_v1 and in SOM_Heterogeneous
//			2. added multithreading to SOM_v1 and heterogeneous that cut off serious time. Hooray.
// 07/07/2013: added an option to account for wall-losses of gas-phase species
//			assumes a first order loss, with gasWLR in units of per second.

FUNCTION SOM_MC([SOMparams,InitialConcMatrix,fitcoefs,fittime,quiet])
	wave SOMparams// 2D wave; rows = parameters for a given compound; columns = different compounds
					// [0] = Ncarbons for each species
					// [1] = FragSlope
				 	// [2] = delta_logCstar_perO
				 	// [3] = Pfunc[0]
				 	// [4] = Pfunc[1]
				 	// [5] = Pfunc[2]
				 	// [6] = Pfunc[3]
				 	// [7] = krxn (leave 0 if default should be used)
	wave InitialConcMatrix // 2D wave; each column contains concentrations for a given VOC precursor class (in ppm)
						// each row contains the conc. for a given Nc species, starting from the top
				 	
	wave fitcoefs		// [0] = FragSlope
				 	// [1] = delta_logCstar_perO
				 	// [2] = Pfunc[0]
				 	// [3] = Pfunc[1]
				 	// [4] = Pfunc[2]
				 	// [5] = Pfunc[3]
				 	// [6] = GasWLR	// gas-phase wall loss rate
	wave fittime		// time wave to use, if fitting data
	variable quiet		// 0 (default) = print, 1 = limited printing
	
	// everything should happen in root folder...yes, this is annoying but more trouble to rewrite everything
	string df_home = ksDFroot // "root:"
	string df_base = ksDFVOCbase // "root:VOCs" // Folder in which waves for individual VOCs are stored (in subfolders)
	string df_VOC // TBD; subfolders for storing individual waves from individual VOCs
	string cdf // current data folder

	setdatafolder $df_home
	cdf = getdatafolder(1)
//	setdatafolder root:
	
// Deal with number of precursor VOCs (more than one is a fundamental shift...need to check calculations carefully)
	variable nVOCs // # of precursor VOCs
	if(ParamIsDefault(SOMparams))
		nVOCs = 1 // 1 is the default
	else
		nVOCs = dimsize(SOMparams,1) // # of columns = # of precursor VOCs
	endif
	variable idex_VOCs // index for looping around different VOCs
	
	make/o/d/n=(nVOCs)/FREE nCforVOCs
	nCforVOCs = SOMparams[0][p]
	wavestats/q nCforVOCs
	variable nCmax = V_max // maximum number of VOCs...all waves will be based on this, but populated according to actual nC for each compound

	if(ParamIsDefault(quiet))
		quiet = 0		
	endif

//***//***// This bit is if 1 compound only //**//**//**
	// specify SOM parameters
	NVAR FragSlope = $(df_home + "FragSlope") // either mfrag or cfrag, depending on fragmentation method
	NVAR delta_logCstar_perO = $(df_home + "delta_logCstar_perO") // change in logCstar per oxygen added
	NVAR ProbOx1 = $(df_home + "ProbOx1")
	NVAR ProbOx2 = $(df_home + "ProbOx2")
	NVAR ProbOx3 = $(df_home + "ProbOx3")
	NVAR ProbOx4 = $(df_home + "ProbOx4")	// probability of adding some number of oxygen atoms
	NVAR gasWLR = $(df_home + "gasWLR")			// gas phase wall loss first order rate coefficient (1/s)
	variable numFitCoefs = numpnts(FitCoefs)	// number of elements in FitCoefs; needed to deal with more than standard 6 FitCoefs
	if(ParamIsDefault(fitcoefs))
		// use values from panel
		if(ProbOx1<0)
			  ProbOx1 = 1e-4 
		endif
		if(ProbOx2<0)
		 	 ProbOx2 = 1e-4 
		endif
		 if(ProbOx3<0)
		 	 ProbOx3 = 1e-4 
		 endif
		 if(ProbOx4<0)
		 	 ProbOx4 = 1e-4 
		 endif
	else
		// set values from FitCoefs
		FragSlope = fitcoefs[0]
		delta_logCstar_perO = fitcoefs[1]
		ProbOx1 = fitcoefs[2]
		 if(ProbOx1<0)
			  ProbOx1 = 1e-4 
	  	 endif
		ProbOx2 = fitcoefs[3]
		 if(ProbOx2<0)
		 	 ProbOx2 = 1e-4 
		 endif
		ProbOx3 = fitcoefs[4]
		 if(ProbOx3<0)
		 	 ProbOx3 = 1e-4 
		 endif
		ProbOx4 = fitcoefs[5]
		 if(ProbOx4<0)
		 	 ProbOx4 = 1e-4 
		 endif
		 if(numFitCoefs > 6)
		 	gasWLR = fitCoefs[6]
		 endif
	endif
	// make a matrix to hold oxygen addition probabilities for later use
	// normalize so that total probability sums to 1
	variable ProbOxSum = ProbOx1 + ProbOx2 + ProbOx3 + ProbOx4
	ProbOx1 /= ProbOxSum
	ProbOx2 /= ProbOxSum
	ProbOx3 /= ProbOxSum
	ProbOx4 /= ProbOxSum
//***//***// Above bit is if 1 compound only //**//**//**
	
	// For multiple VOCs
		variable mfrag
		variable DLVP
		variable POx1
		variable POx2
		variable POx3
		variable POx4
		ProbOxSum = POx1 + POx2 + POx3 + POx4
		POx1 /= ProbOxSum
		POx2 /= ProbOxSum
		POx3 /= ProbOxSum
		POx4 /= ProbOxSum
	
	// Specify time step
	NVAR nsteps = $(df_home + "nSteps") // number of time steps --> this is not an input, but will be calculated later
	NVAR timestep = $(df_home + "timestep")// timestep in seconds --> this is an input, may be adjusted below
	if(gasWLR > 0)	// deal with issues of numerical stability
		if(timestep > 0.05/gasWLR)
			timestep = 0.05/gasWLR
		endif
	endif
	NVAR MaxTime_hours = $(df_home + "MaxTime_hours") // length of simulation, in hours --> this is an output
	nsteps = ceil(MaxTime_hours*60*60/timestep)+1	// total number of time steps
	
	// Set this variable to control number of iterations during fitting
	// This doesn't always seem to do anything...
	variable V_FitMaxIters = 10 // max number of iterations
	
	if(quiet==0)
		print "//"
		print "SOM-MP Run at " + time() + " on " + date()
	endif
	
// Access variables from panel, and create a few more variables for use
// Oligomerization
	NVAR OligomerizationIsOn = $(df_home + "OligomerizationIsON")// oligomerization is on (1) or off (0)
	NVAR krxn_base_olig	= $(df_home + "krxn_base_olig")	// oligomerization bimolecular rate coefficient
// Gas-phase wall loss
	variable gasWLmethod = 1 // 0 = irreversible uptake, 1 = reversible uptake
	NVAR Cwall	= $(df_home + "Cwall")			// "effective" wall concentration; 10 mg/m3 is default choice
	NVAR kwall_gas_scaling = $(df_home + "kwall_gas_scaling")	// composition dependent wall loss rate coefficient; not really used 
	kwall_gas_scaling = 0 // set default to zero on 12/4/16
	NVAR noSOA = $(df_home + "noSOA")			// Set this to 1 if just the gas-phase chemistry is to be run w/o SOA formation
	NVAR Krechmer_Cw = $(df_home + "Krechmer_Cw")	// Set this to 0 to use a constant Cw value; Set to 1 to use Krechmer (2016) relationship

	NVAR Ncarbons = $(df_home + "Ncarbons")// number of carbon atoms in parent molecule --> input for 1 precursor
	variable Noxygens_parent = 0	// adjust this to change the number of oxygens on the parent molecule
	variable Nox_precursor = 0 // number of oxygens in precursor species...oops, created this twice. Same as Noxygens_parent
	variable nC
	variable nO 
	nC = nCarbons
//	nO = ceil(20/(nC^0.33)) //+ 2 // Max number of oxygen atoms that can be added
	nO = 8 // 8 seems to provide robust results. Seven also works reasonably well, but can limit O:C
	variable Nspecies = nCmax*nO		// added 08/26/14 for dynamic partitioning; updated to nCmax (from nC) on 10/24/2015

	if(OligomerizationIsOn==1) // Oligomerization is not ready for primetime yet (note from 10/24/15)
		nC = nC*2
		nO = nO*2 // max number of oxygen atoms that can be added
	endif
	
	// concentrations
	NVAR ctot_ppm = $(df_home + "ctot_ppm") // concentration of parent hydrocarbon (gas + particle) in ppm
	variable Ctot_init // initial mass loading (in ug/m^3) of the total organic (gas + particle) material
	// info about particles
	NVAR Nparticles = $(df_home + "Nparticles") // #/cm^3 --> input
	variable Np = Nparticles*1e6 // #/cm^3 --> #/m^3
	NVAR Density_Base = $(df_home + "Density_Base") // g/cm^3 --> input
	variable density = density_base*1000 // g/cm^3 --> kg/m^3
	// parameters related to kinetics
	NVAR krxn_parent = $(df_home + "krxn_parent") // cm3/molecules.s; set to zero to use default values
	NVAR krxn_method = $(df_home + "krxn_method") // 4 is the method of choice, now (10/24/15), so hardwired below
	krxn_method = 4
	NVAR O3_yn = $(df_home + "O3_yn") // 0 = OH reactions, 1 = O3 reaction (no secondary chemistry)
	NVAR O3_conc = $(df_home + "O3_conc") // O3 concentration, if O3_yn = 1
	NVAR usescalingforoxidant = $(df_home + "usescalingforoxidant") // 0 = constant; 1 = exponential decay; 2 = interpolate observed [OH] or [O3] to model timebase (requires experimental data)
	// Parameters for Van Krevelen
	NVAR H_per_O = $(df_home + "H_per_O") // user selectable; has minimal influence on anything except H:C ratio
	NVAR Hadjustment = $(df_home + "Hadjustment") // user selectable; meant to account for difference between MW calculation for hydrocarbon (H = 2C+2) and actual VOC
						// however, we now just set this to 0 all the time, as it doesn't really matter (i.e. gets accounted for in fit parameters and makes life simple)
	// some constants
	variable m, n, o, k, i, j
	variable Na = 6.022e23	 // molecules/mol...Avagadros Number
	NVAR Temp_G = $(df_home + "Temp_G")			// temperature, in Kelvin, as global variable
	NVAR Pressure_G = $(df_home + "Pressure_G")		// pressure, in atm, as global variable
	variable Tvar = Temp_G	// Kelvin
	variable Pressure = Pressure_G * 760 * 133.322 // Pa
	variable NumberDensityAir = Pressure / (1.381e-23 * Tvar) // molecules per m^3
	NVAR EqmMethod = $(df_home + "EqmMethod")
	string units
	if(EqmMethod == 0) // "molecules" or "mass"; affects how eqm. partitioning is performed. See EqmCalc().
		units = "molecules"	
	else
		units = "mass"
	endif

// 2. Specify various properties.

  	variable Diff_CO2 = 0.138e-4 // m^2/s, diffusion coefficient for CO2 (will scale other values from this)
  	variable MW_CO2 = 44	// g/mol, MW of CO2

	// Parameters for size distribution; added 10/16/13; updated 092914
	NVAR SizeSpread = $(df_home + "SizeSpread")		// log-normal standard deviation for seed particle distribution; added 10/16/13
	NVAR SeedDiameter = $(df_home + "SeedDiameter")	// seed particle diameter in nm; added 10/16/13 as an alternative to the seed concentration
  	NVAR SeedVolConc = $(df_home + "SeedVolConc")		// um^3/cm^3, will be calculated
	NVAR SeedSurfaceArea = $(df_home + "SeedSurfaceArea")	// to be calculated
	NVAR polydisperse = $(df_home + "polydisperse")	// 0 = monodisperse, 1 = polydisperse
	NVAR nSizeBinsG = $(df_home + "nSizeBinsG") 	// Global variable (not on panel); set to 7; so nSizeBins just hardwired in next line
	variable nSizeBins = nSizeBinsG 
	variable index_size
	variable dlogDp			// dlogDp associated with generated size distribution
	variable DpStart = SeedDiameter	// initial median particle diameter
  	variable VolumeSeed = (pi/6)*((DpStart*1e-9)^3)	// m^3/particle, seed volume per particle
  	variable VolumeOrgPerParticle		// volume of organic per particle, to be calculated
	
	if(polydisperse==0)
		nSizeBins = 1	// monodisperse
	endif
	
	setdatafolder $df_home // just making sure that we are really set to the appropriate folder
	// make size distributions and calculate some initial size things
	dlogDp = makelognormaldistn(SizeSpread,SeedDiameter,Nparticles,nbins = nsizebins) // (log(Dp2)-log(Dp1)) = constant
	wave Diameter = $(df_home + "Diameter")	// nm, by size
	wave dNdlogDp = $(df_home + "dNdlogDp")	// p/cm^3, conc, by size
	wave dSdlogDp = $(df_home + "dSdlogDp")	// m2/m3 (I think)
	wave dVdlogDp = $(df_home + "dVdlogDp")	// um3/m3 (I think)
	variable NpCalc	// particle number concentration
	make/o/d/n=(nsizebins)/FREE logDp	// log(Dp)
	logDp = log(diameter)
	Integrate/METH=1 dNdlogDp/X=logDp/D=INT_Result
	NpCalc = INT_Result[nsizebins-1]	// integrated number concentration, p/cm3
	Integrate/METH=1 dVdlogDp/X=logDp/D=INT_Result
	SeedVolConc = INT_Result[nsizebins-1]/1e9 // um^3/cm^3
	Integrate/METH=1 dSdlogDp/X=logDp/D=INT_Result
	SeedSurfaceArea = INT_Result[nsizebins-1]/1e6 // um^2/cm^3
  	make/o/d/n=(nSizeBins) VolumeSeedPerBin = (pi/6)*((Diameter[p]*1e-9)^3)	// m^3/particle
	variable SeedMW = 250 // g/mol; DOS = 427
  	make/o/d/n=(nSizeBins) MoleculesSeedPerBin = VolumeSeedPerBin * 1e6 * 1 * (1/SeedMW) * Na * (dlogDp*dNdlogDp) // m^3/p * cm3/m3 * g/cm3 * mol/g * molecules/mol * p/cm^3
  	make/o/d/n=(nSizeBins) VolumeOrgPerBin // m^3/particle
//  	make/o/d/n=(nSizeBins) dDpdt_on, dDpdt_off
  	
  	if(Polydisperse==0) // monodisperse
  		SeedVolConc = VolumeSeed*NParticles*1e18	// um^3/cm^3
		SeedSurfaceArea = (4*pi*(DpStart*1e-3/2)^2)*(Nparticles) // um^2/cm^3
	endif

	// Parameters for Particle Wall Loss; Added 08/20/13
	// empirical expression determined from Loza et al., ACP, 2012
	setdatafolder $df_home
	NVAR ParticleWallLoss = $(df_home + "ParticleWallLoss")	// 0 = no loss, 1 = loss
	variable particleWLR	// size-dependent wall loss rate in per second; calculated at each model time-step below
	make/o/d/n=(numpnts(dndlogdp)) particleWLR_wave
	note particleWLR_wave "units = 1/s; determined from Loza et al, 2012 for Caltech chamber"
	
	 // Parameters for dynamic partitioning, added 11/12/13
  	NVAR KineticMassTransfer 	// 0 = equilibrium calculation; 1 = dynamic partitioning
	NVAR alpha 					// accommodation coefficient, timestep will scale with this
	variable DynamicPartitioningMethod = 0 // set to 0 to use Euler-like method (preferred) 	

 //2c.  //Create some parameters and waves FOR ZAVERI METHOD (first implementation, v7.4)
 	// THIS IS JUST THE START FOR MULTICOMPONENT SOM -- DOES NOT WORK
//	make/d/o/n=(nC,nO) Diffusivity_matrix = Diff_CO2*(MW_CO2/MW_matrix) // m^2/s; assume that Diffusivity scales with MW and use CO2 as reference case
//	make/o/d/n=(nC,nO) MeanFreePath_Matrix = 3*Diffusivity_Matrix/sqrt((8*1.381e-23*Tvar)/(pi*(MW_matrix/(1000*Na))))	// meters; D/(lambda*c_bar) = 1/3, for Fuchs expression
//	//
//	make/o/d/n=(nC,nO,nSizeBins) SatRatio_3D = 0
//	make/o/d/n=(nC,nO,nSizeBins) tmatrix_3D=0
//	make/o/d/n=(nC,nO,nSizeBins) Phi_GPP = 0
//	make/o/d/n=(nC,nO) tarray=0
//	make/o/d/n=(nC,nO,nSizeBins) Knudsen_3D = nan, Fuchs_3D = nan
//	make/o/d/n=(nC,nO,nSizeBins) deltaDynamicPartitioningPre = nan, deltaDynamicPartitioning=nan
//	make/o/d/n=(nSizeBins) Kelvin
	variable deltat_min, deltat_min_new
	variable ntime_GPP, ntime_GPP_new // the rounded number of points that can be allowed for a GPP timestep

	NVAR MaxTime_hours // length of simulation, in hours // v7.4 moved from above allow for timestep to change, as necessary for dynamic partitioning
	nsteps = ceil(MaxTime_hours*60*60/timestep)//+1	// total number of time steps


	 // Method for dynamic partitioning
	variable MassTransferMaxTime = 0.3	// seconds, for operator split, added 11/12/13
  		// specify time step for dynamic partitioning based on alpha
  		// Usual timestep is 60 seconds, so useful to keep as divisors of 60
  	if(alpha < 1e-4)
  		MassTransferMaxTime = 30
  	elseif(alpha < 0.5e-3)
  		MassTransferMaxTime = 30
  	elseif(alpha <= 1e-2)
  		MassTransferMaxTime = 10
  	elseif(alpha <= 1e-1)
  		MassTransferMaxTime = 10//1
  	elseif(alpha <= 1)
  		MassTransferMaxTime = 0.3
  	endif
   	variable MassTransferTimeScaling = round(TimeStep/MassTransferMaxTime)	// number of timesteps in operator split, added 11/12/13
   	MassTransferTimeScaling = MassTransferTimeScaling < 1 ? 1 : MassTransferTimeScaling

	// absorbing seed
	NVAR AbsorbingSeed		// 0 = no, 1 = yes
	variable SeedMass = 0.0000000005 // ug/m^3
	if(AbsorbingSeed == 0)
		SeedMass = 1e-6		// ug/m^3
	else
		SeedMass = SeedVolConc * 1.0	// ug/m^3, and where 1.0 is the assumed density
	endif
	variable SeedMolecules = SeedMass * 1e-6 * 1e-6 * Na / SeedMW // molecules/cm^3
	NVAR DilutionVar //variable DilutionVar = 0//6.12 // % per hour (6.12% from Dzepina et al., 2011)
	NVAR logCstar_adjustmentfactor //variable logCstar_adjustmentfactor = 0.9 // multiplicative factor to adjust logCstar values. Important for isomer simulations

// Decide how to represent volatility...legacy and not currently in use (5/27/13)	
	string Cstar_evolution = "no" // yes to allow for time varying functional group evolution
	variable delta_logCstar_perO_final = 2.2 // relative to delta_logCstar_perO
	variable delta_logCstar_perO_perStep = (delta_logCstar_perO_final-delta_logCstar_perO)/nsteps
// Decide how to implement fragmentation	
	NVAR Pfrag_type // 0 = cfrag*Nox; 1 = (O:C)^mfrag
	NVAR Frag_Method = small_fragments // 0 for random number generator, 1 for equal probs., 2 to only produce fragments with 1 carbon (HCHO, CO2, CH4)
  // some other silly things
	variable steady_state = 0, SSvalue = 0 // 0 for base case, any number for something else keep parent gas-phase abundance at initial abundance
// OH concentration...many of these are not in use (5/27/13)
	NVAR OHconc // [OH] in molecules/cm^3
	NVAR OH_scale // scaling factor to have [OH] decrease exponentially with time; OHconc_t = OHconc_0*exp(-timewv[idex-1]*oh_scale)
	variable OH_profile = 2 // 1 to have a diurnal profile of OH concentrations; 2 for anything else
	variable OH_counter = 0 // related to diurnal profile
	variable day_counter = 0 // related to diurnal profile
	variable day_iseven = 1 // related to diurnal profile
	variable OH_max = 4e6 // related to diurnal profile
	variable OH_min = 1e5 // related to diurnal profile
// Parameters for heterogeneous chemistry
	variable MW_OH = 17/6.022e23	// molecular weight in g/molecule
	variable Diff_OH = 0.3e-5			// m^2/s
	variable MFP_OH = 3*Diff_OH/sqrt((8*1.381e-23*Tvar)/(pi*MW_OH/1000))	// mean free path in meters, for Fuchs-Sutugin correction
//	variable OH_Flux = (OHconc*1e6)*1.381e-23*Tvar/sqrt(2*pi* (17/1000/6.022e23)*1.381e-23*Tvar) // molecules m^-2 s^-1, will be calculated on the fly
//	variable OH_Rate // molecules/s; will be calculated from OH_Flux and particle surface area
	NVAR gammaOH // OH reactive uptake coefficient
	NVAR hetchem // 0 for no heterogeneous chemistry; 1 to include
// First generation products?
	NVAR FirstGenProductsOnly
// Parameters for sequential partitioning model (SPM)
	NVAR SPM // 0 to ignore, 1 to include
// Parameters for timing
	variable t11, t22
			
//3. Write a bit of info to a wave
	string RunInfoStr = ""
	make/O/T RunInfo = ""
	make/o/d RunValues

	RunInfoStr = "nC;Ctot_init_ug_m3;Ctot_ppm;OHconc_molec_cm3;OH_scale;FragSlope;delta_logCstar_perO;ProbOx1;ProbOx2;ProbOx3;ProbOx4;"
	RunInfoStr += "PfragMethod;SeedMass;SeedMW;AbsorbingSeed;RunTime;nsteps;timestep;DilutionVar;Small_Fragments;H_per_O;Hadjustment;hetchem;gammaOH;SPM;"
	RunInfoStr += "krxn_parent;O3_yn;O3_conc_ppb;logCstar_adjustmentfactor;kwall_gas;kwall_aer;Oligomerization;krxn_base_olig;EqmMethod"
	i= ItemsInList(RunInfoStr)
	Make/O/T/N=(i) RunInfo= StringFromList(p,RunInfoStr)
	RunValues = {nC,Ctot_init,Ctot_ppm,OHconc,OH_scale,FragSlope,delta_logCstar_perO,ProbOx1,ProbOx2,ProbOx3,ProbOx4}
	wavestats/Q RunValues
	RunValues[x+V_npnts] = {Pfrag_type,seedmass,seedMW,AbsorbingSeed,MaxTime_Hours,nsteps,timestep,DilutionVar,Frag_Method,H_per_O,Hadjustment,hetchem,gammaOH,SPM}
	wavestats/Q RunValues
	RunValues[x+V_npnts] = {krxn_parent,O3_yn,O3_conc,logCstar_adjustmentfactor,gasWLR,ParticleWallLoss,OligomerizationIsOn,krxn_base_olig,EqmMethod}

// 4. Create a bunch of waves (multi-dimensional) for use in calculations	
	
	SOM_MakeWaves(SOMparams,nSteps,nSizeBins,nCmax,nO,timestep) // creates a bunch of VOC-specific waves in folders
	// The number of waves created depends on the number of columns (precursor classes) in the wave "SOMparams"
	// Everything is created in the "VOC" folder and subfolders therein; things are not created in root:
	
	// now make general waves (non-VOC specific) in root folder
	setdatafolder $df_home // create in root
	
	make/d/o/n=(nsteps) SeedConc_time = SeedVolConc
	setscale/P x, 0, (timestep/60/60), "hours", SeedConc_time
	//
	make/d/o/n=(nsteps) TimeW = 0
	note TimeW "Reaction time [hrs]"
	//
	make/d/o/n=(nsteps) O2C_time = nan
	note O2C_time "Oxygen-to-Carbon ratio of total SOA"
	setscale/P x, 0, (timestep/60/60), "hours", O2C_time
	setscale/P d, 0, 0, "O:C", O2C_time
	make/d/o/n=(nsteps) H2C_time = nan
	note H2C_time "hydrogen-to-carbon ratio of total SOA"
	setscale/P x, 0, (timestep/60/60), "hours", H2C_time
	setscale/P d, 0, 0, "H:C", H2C_time
	make/d/o/n=(nsteps) OCseed_time = 0
	//
	make/d/o/n=(nsteps) Coa_time = 0
	note Coa_time "Total SOA mass concentration [ug/m^3]"
	note Coa_time "nVOCs used = " + num2str(nVOCs) 
	note Coa_time "Pfrag = using type " + num2str(Pfrag_type) + " and method " + num2str(Frag_Method)
	setscale/P x, 0, (timestep/60/60), "hours", Coa_time
	setscale/P d, 0, 0, "ug/m\S3\M", Coa_time
	//
	make/d/o/n=(nsteps) Coa_time_wallcorr = 0
	note Coa_time_wallcorr "Particle wall-loss corrected total SOA mass concentration [ug/m^3]"
	note Coa_time_wallcorr "Coa_time_wallcorr = Coa_time + WallParticleMolec_Matrix"
	setscale/P x, 0, (timestep/60/60), "hours", Coa_time_wallcorr
	setscale/P d, 0, 0, "ug/m\S3\M", Coa_time_wallcorr
	//
	make/d/o/n=(nsteps) deltaHC_time = 0
	setscale/P x, 0, (timestep/60/60), "hours" deltaHC_time
//	make/d/o/n=(nsteps) HC_ppb_time = 0
//	note HC_ppb_time "Total Parent HC concentration [ppb]"
//	setscale/P x, 0, (timestep/60/60), "hours" HC_ppb_time
//	setscale/P d, 0, 0, "ppb" HC_ppb_time
//	make/d/o/n=(nsteps) deltaHCfrac_time = 0
//	setscale/P x, 0, (timestep/60/60), "hours", deltaHCfrac_time
	make/d/o/n=(nsteps) Yield_time = 0
	note Yield_time "Total aerosol mass yield"
	setscale/P x, 0, (timestep/60/60), "hours", Yield_time
	//
	make/d/o/n=(nsteps) OH_wave = OHconc
	note OH_wave "[OH] in molecules/cm^3"
	setscale/P x, 0, (timestep/60/60), "hours", OH_wave
	setscale/P d, 0, 0, "molecules/cm\S3\M", OH_wave
	make/d/o/n=(nsteps) OH_exposure = 0
	make/d/o/n=(nsteps) Lifetimes = 0
	note Lifetimes "Number of oxidation lifetimes"
	setscale/P x, 0, (timestep/60/60), "hours", Lifetimes
	//
	make/o/d/n=(nsteps) Dp_time = 0; 
	Dp_time = DpStart // particle diameter in nm
	note Dp_time, "Particle diameter, including seed [nm]"
	setscale/P x, 0, (timestep/60/60), "hours", Dp_time
	setscale/P d, 0, 0, "nm", Dp_time
	//
	make/o/d/n=(nsteps) Np_time = Np/1e6	// particle number concentration
	note Np_time, "Particle Number Concentration [p/cc]"
	setscale/P x, 0, (timestep/60/60), "hours", Np_Time
	//
	make/o/d/n=(nsteps,nSizeBins) dNdlogDp_time=0	// updated 9/29/14; added 10/16/13
	dNdlogDp_time[][] = dNdlogDp[q]
	setscale/P x, 0, (timestep/60/60), "hours", dNdlogDp_time
	note dNdlogDp_time "dNdlogDp, concentration is p/m^3"
	make/o/d/n=(nsteps,nSizeBins) dVdlogDp_time=0
	setscale/P x, 0, (timestep/60/60), "hours", dVdlogDp_time
	//
	make/o/d/n=(nSizeBins) dNdlogDp_final=0	// updated 9/29/14; added 10/16/13
	dNdlogDp_final[][] = dNdlogDp[q]
	note dNdlogDp_time "dNdlogDp, concentration is p/m^3"
	//
	make/o/d/n=(nsteps,nSizeBins) NpPerBin_time=0	// added 9/29/14
	NpPerBin_time[][] = dNdlogDp[q]*dlogDp
	setscale/P x, 0, (timestep/60/60), "hours", NpPerBin_time
	note NpPerBin_time "concentration per bin is p/m^3"
	//
	make/o/d/n=(nsteps,nSizeBins) Diameter_time=0	// updated 9/29/14; added 10/16/13
	Diameter_time[][] = Diameter[q]
	setscale/P x, 0, (timestep/60/60), "hours", Diameter_time
	note Diameter_time "Changing Diameter associated with size distribution, in nm"	
	//
	make/o/d/n=(nSizeBins) Diameter_final=0	// added 02/27/16
	Diameter_final[][] = Diameter[q]
	note Diameter_final "Final Diameter associated with size distribution, in nm"	
	
	// TOTAL molecules in particle phase (sum over all VOCs)
	make/d/o/n=(nCmax,nO) ParticleMolecules_Matrix = 0 // molecules/cm^3
	setscale/P x, nCmax, -1, "nCarbons" ParticleMolecules_Matrix
	setscale/P y, 0, 1, "nOxygens" ParticleMolecules_Matrix
	note ParticleMolecules_Matrix "particle phase species concentration in molecules/cm^3"
	note ParticleMolecules_Matrix "Ncarbons = Rows; Noxygens = columns"
	wave ParticleMoleculesSum_Matrix = ParticleMolecules_Matrix // this lives in root:, as do all of these waves
	// TOTAL mass in particle phase (sum over all VOCs)
	make/d/o/n=(nCmax,nO) ParticleMass_Matrix = 0
	setscale/P x, nCmax, -1, "nCarbons" ParticleMass_Matrix
	setscale/P y, 0, 1, "nOxygens" ParticleMass_Matrix
	wave ParticleMassSum_matrix = ParticleMass_Matrix
	// TOTAL molecules in gas phase (sum over all VOCs)
	make/d/o/n=(nCmax,nO) GasMolecules_Matrix = 0
	setscale/P x, nCmax, -1, "nCarbons" GasMolecules_Matrix
	setscale/P y, 0, 1, "nOxygens" GasMolecules_Matrix
	wave GasMoleculesSum_matrix = GasMolecules_Matrix
	// TOTAL mass in gas phase (sum over all VOCs)
	make/d/o/n=(nCmax,nO) GasMass_Matrix = 0
	setscale/P x, nCmax, -1, "nCarbons" GasMass_Matrix
	setscale/P y, 0, 1, "nOxygens" GasMass_Matrix
	wave GasMassSum_matrix = GasMass_Matrix
	// particle molecules, time
	make/d/o/n=(nCmax,nO,nsteps) GasMolecules_Time = 0
	note GasMolecules_Time "gas-phase species concentration in molecules/cm^3"
	note GasMolecules_Time "Ncarbons = Rows; Noxygens = columns; Layers = time [hrs]"
	setscale/P z, 0, (timestep/60/60), "hours", GasMolecules_Time
	setscale/P x, nCmax, -1, "nCarbons" GasMolecules_Time
	setscale/P y, 0, 1, "nOxygens" GasMolecules_Time
	wave GasMoleculesSum_Time = GasMolecules_Time
	// particle molecules, time
	make/d/o/n=(nCmax,nO,nsteps) ParticleMolecules_Time = 0
	note ParticleMolecules_Time "particle-phase species concentration in molecules/cm^3"
	note ParticleMolecules_Time "Ncarbons = Rows; Noxygens = columns; Layers = time [hrs]"
	setscale/P z, 0, (timestep/60/60), "hours", ParticleMolecules_Time
	setscale/P x, nCmax, -1, "nCarbons" ParticleMolecules_Time
	setscale/P y, 0, 1, "nOxygens" ParticleMolecules_Time
	wave ParticleMoleculesSum_Time = ParticleMolecules_Time
	// particle mass, time
	make/d/o/n=(nCmax,nO,nsteps) ParticleMass_Time = 0
	note ParticleMass_Time "particle-phase species concentration in molecules/cm^3"
	note ParticleMass_Time "Ncarbons = Rows; Noxygens = columns; Layers = time [hrs]"
	setscale/P z, 0, (timestep/60/60), "hours", ParticleMass_Time
	setscale/P x, nCmax, -1, "nCarbons" ParticleMass_Time
	setscale/P y, 0, 1, "nOxygens" ParticleMass_Time
	wave ParticleMassSum_time = ParticleMass_time
	// 3D version of ParticleMolecules_matrix
	make/o/d/n=(nCmax,nO,nSizeBins) PM3D_matrix = 0 
	setscale/P x, nCmax, -1, "nCarbons" PM3D_matrix
	setscale/P y, 0, 1, "nOxygens" PM3D_matrix
	wave PM3Dsum_matrix = PM3D_matrix
	// 3D version of ParticleMass_matrix
	make/o/d/n=(nCmax,nO,nSizeBins) PMM3D_matrix = 0 
	setscale/P x, nCmax, -1, "nCarbons" PMM3D_matrix
	setscale/P y, 0, 1, "nOxygens" PMM3D_matrix
	wave PMM3Dsum_matrix = PMM3D_matrix
	// 3D version of ParticleMoleFraction_matrix
	make/o/d/n=(nCmax,nO,nSizeBins) PMF3D_matrix = 0 
	setscale/P x, nCmax, -1, "nCarbons" PMF3D_matrix
	setscale/P y, 0, 1, "nOxygens" PMF3D_matrix
	wave PMF3Dsum_matrix = PMF3D_matrix
	// SOA associated with this VOC
	make/d/o/n=(nsteps) Coa_time = 0
	note Coa_time "SOA mass concentration [ug/m^3]"
	setscale/P x, 0, (timestep/60/60), "hours", Coa_time
	setscale/P d, 0, 0, "ug/m\S3\M", Coa_time
	wave CoaSum_Time = Coa_time
	// wall corrected value
	make/d/o/n=(nsteps) Coa_time_wallcorr = 0
	note Coa_time_wallcorr "Particle wall-loss corrected SOA mass concentration [ug/m^3]"
	note Coa_time_wallcorr "Coa_time_wallcorr = Coa_time + WallParticleMolec_Matrix"
	setscale/P x, 0, (timestep/60/60), "hours", Coa_time_wallcorr
	setscale/P d, 0, 0, "ug/m\S3\M", Coa_time_wallcorr
	wave CoaSum_Time_wallcorr = Coa_time_wallcorr
	// Amount reacted for all compounds
	make/d/o/n=(nsteps) deltaHC_time = 0
	setscale/P x, 0, (timestep/60/60), "hours" deltaHC_time
	note deltaHC_time "total amount reacted, in ug/m3"
	setscale/P d, 0, 0, "ug/m\S3\M", deltaHC_time
	wave deltaHCsum_time = deltaHC_time
	make/d/o/n=(nsteps) HC_ugm3_time = 0
	note HC_ugm3_time "Sum of parent HC concentrations [ug/m3]"
	setscale/P x, 0, (timestep/60/60), "hours" HC_ugm3_time
	setscale/P d, 0, 0, "ug/m\S3\M" HC_ugm3_time
	wave HCsum_ugm3_time = HC_ugm3_time
	// Effective yield for all compounds
	make/d/o/n=(nsteps) Yield_time = 0
	note Yield_time "Aerosol mass yield over all compounds"
	setscale/P x, 0, (timestep/60/60), "hours", Yield_time
	wave YieldSum_time = Yield_time
	//
	make/d/o/n=(nCmax) GasMolecules_Parent = 0
	variable HCstart, HCfinish
	//
	make/o/d/n=(nCmax,nO) ParticleVolume_matrix
	//
	// O2C for this compound class
	make/d/o/n=(nsteps) O2C_time = nan
	note O2C_time "Oxygen-to-Carbon ratio of SOA"
	setscale/P x, 0, (timestep/60/60), "hours", O2C_time
	setscale/P d, 0, 0, "O:C", O2C_time
	wave O2Csum_time = O2C_time
	make/d/o/n=(nsteps) H2C_time = nan
	note H2C_time "hydrogen-to-carbon ratio of SOA"
	setscale/P x, 0, (timestep/60/60), "hours", H2C_time
	setscale/P d, 0, 0, "H:C", H2C_time
	wave H2Csum_time = H2C_time
	make/d/o/n=(nCmax,nO) TotalCarbonAtoms = 0
	make/d/o/n=(nCmax,nO) TotalOxygenAtoms = 0
	make/d/o/n=(nCmax,nO) TotalHydrogenAtoms = 0	
	//
	variable Ave_nC
	variable Ave_nO
	variable Ave_nH
	variable molecules_PM_previous = 0
	
	TimeW = (x)*timestep/60/60 // hours

// Loop around for different VOCs
	cdf = getdatafolder(1) // get current data folder
	setdatafolder $df_home // this is likely root:
	for(idex_VOCs=0;idex_VOCs<nVOCs;idex_VOCs+=1)
		df_VOC = df_base + ":VOC" + num2istr(idex_VOCs+1) // set to VOC specific data folder
		setdatafolder $df_VOC
	// 2. Populate O:C, H:C, Oxidation State and MW Matrices for each VOC
		SOM_CreateAtomMatrices(nCmax,nO,H_per_O,Hadjustment) // changed from nC to nCmax
		//$$$ Will need to access these later
		//		wave C_matrix, O_matrix, H_matrix, OC_matrix, HC_matrix, Ox_matrix, MW_matrix
	
	// 2. Populate Diffusivity and mean free path matrices; these exist in the VOC subfolders
		wave Diffusivity_Matrix, MW_matrix, MeanFreePath_Matrix
		Diffusivity_Matrix = Diff_CO2*(MW_CO2/MW_matrix)	// m^2/s; assume that Diffusivity scales with MW and use CO2 as reference case
		MeanFreePath_Matrix = 3*Diffusivity_Matrix/sqrt((8*1.381e-23*Tvar)/(pi*(MW_matrix/(1000*Na))))	// meters; D/(lambda*c_bar) = 1/3, for Fuchs expression	
	// 3. Populate initial logC* matrix
		wave logCstar_Matrix, MW_matrix, SatConc_matrix
		delta_logCstar_perO = SOMparams[2][idex_VOCs] // get VOC specific value
		logCstar_Matrix[][] = -0.0337*MW_Matrix[p] + 11.56 // estimated from Lide, CRC data for saturated hydrocarbons
		logCstar_Matrix[][] = logCstar_Matrix[p][0] - y*delta_logCstar_perO // adjust C* values based on number of oxygens per molecule
		logCstar_matrix += logCstar_adjustmentfactor // one can "adjust" the logCstar matrix from the base case assumptions; best to have this be 0; FIXED from * to + on 8/4/18
		SatConc_Matrix = (10^(logCstar_Matrix)) * 1e-6 * Na / MW_Matrix // molecules/m^3; C* --> saturation concentration

	// 4. Enthalpy of vaporization; added 09/03/13
		SOM_deltaHvap(logCstar_Matrix=logCstar_matrix,ConstantDHvap=0)	// ConstantDHvap = 0 means not constant, otherwise enter values in kJ/mol; operates in current folder
		//$$$ wave deltaHvap_matrix
	
	// 5. OH reaction rate coefficients
		wave C_matrix, O_matrix // these may be be "oversized" depending on max number of carbons in the biggest VOC
		if(OligomerizationIsOn==0)	// no oligomerization
			SOM_RateCoefficients(C_matrix,O_matrix,0,O3_yn,TempK=Tvar,method=krxn_method) // operates in current VOC folder
			wave krxn_matrix
			// Adjust parent rate coefficient if desired
			if(SOMparams[7][idex_VOCs] == 0)
				// do nothing
			else
				krxn_matrix[nCmax - SOMparams[0][idex_VOCs]][0] = SOMparams[7][idex_VOCs]// // need to adjust appropriate element for oversized matrices
			endif
		else	// oligomerization $$$ Likely broken for multiple compounds
			SOM_RateCoefficients(C_matrix,O_matrix,0,O3_yn,TempK=Tvar,method=krxn_method)
			wave krxn_matrix
//			krxn_matrix[nCarbons][Noxygens_parent] = krxn_parent
			if(SOMparams[7][idex_VOCs] == 0)
				// do nothing
			else
				krxn_matrix[SOMparams[0][idex_VOCs]][0] = SOMparams[7][idex_VOCs]// // need to adjust appropriate element for oversized matrices
			endif
//			krxn_matrix[SOMparams[0][idex_VOCs]][0] = SOMparams[7][idex_VOCs] // set to VOC specific value
		endif
		setscale/P x, nC, -1, "nCarbons" krxn_matrix
		setscale/P y, 0, 1, "nOxygens" krxn_matrix

		if(FirstGenProductsOnly==1)
			krxn_matrix[][1,] = 0
		endif
	
		if(OligomerizationIsOn==1)
			SOM_OligomerRateCoef(nCmax*2,nO,krxn_base=krxn_base_olig) // NEED TO CHECK THIS AFTER UPDATING FOR MULTIPLE VOCS $$$
			wave krxn_matrix_olig
		endif
	
	// 6. Populate fragmentation probability matrix, which is a 3D matrix with rows = #carbons, columns = #oxygens, layers = probability
		FragSlope = SOMparams[1][idex_VOCs] // update for each VOC; at this point this should be a global variable that lives in root
//		SOM_Fragmentation(Frag_Method,FragSlope,Pfrag_type,C_matrix,O_matrix,OC_matrix)
		wave C_matrix_small, O_matrix_small, OC_matrix_small // subsets of the "big" matrices to deal with varying carbon numbers between precursors
		SOM_Fragmentation_MP(Frag_Method,FragSlope,Pfrag_type,nCmax,C_matrix_small,O_matrix_small,OC_matrix_small) // operates in current folder
		// wave Prob1_Matrix, Prob2_Matrix, Prob_Matrix, Frag1_Matrix, Frag1_Array // $$$ WILL NEED TO ACCESS LATER

	// 7. Decide here how to input initial concentration	
		// Completely updated for dealing with multiple compounds
		// InitialConcMatrix wave holds the initial concentrations (in ppm)
		wave GasMass_matrix, GasMolecules_Matrix, TotalMolecules_Matrix, MW_matrix, H_matrix

// CDC 07/05/17		GasMass_matrix[nCmax-SOMparams[0][idex_VOCs],][0] = InitialConcMatrix[p-(nCmax-SOMparams[0][idex_VOCs])][idex_VOCs]
		GasMass_matrix[][0] = InitialConcMatrix[p][idex_VOCs] // updated 07/05/17 by CDC
		GasMass_matrix = GasMass_matrix*(NumberDensityAir*1e-6/Na)*(MW_matrix*1e6)
		GasMass_Matrix = (numtype(H_Matrix) == 2 ? NaN : GasMass_Matrix) // take care of pesky NaN values
		GasMolecules_Matrix = GasMass_Matrix * 1e-6 * 1e-6 * Na / MW_Matrix // convert ug/m^3 to molecules/cm^3
		TotalMolecules_Matrix = GasMolecules_Matrix // (molecules/cm^3)
		GasMoleculesSum_matrix += GasMolecules_matrix // added 08/04/18
		GasMassSum_matrix += GasMass_matrix // added 08/04/18
		setdatafolder $df_home
	endfor // end of VOC precursor loop


	// Calculate initial Equilibrium Distribution (for multiple VOCs)
	EqmCalc_MultipleVOC(SOMparams,SeedConc=SeedMolecules,SeedMW=SeedMW,units="mass") // changed to mass on 08/04/18
  	
	cdf = getdatafolder(1)
	// Populate concentration matrices for different VOCs in their subfolders
	for(idex_VOCs=0;idex_VOCs<nVOCs;idex_VOCs+=1) // Start VOC looping
		df_VOC = df_base + ":VOC" + num2istr(idex_VOCs+1) // set to VOC specific data folder
		setdatafolder $df_VOC
		// Access the waves you want
		wave GasMolecules_Matrix, ParticleMolecules_matrix, TotalMolecules_matrix
		wave PM3D_matrix
		wave Coa_eqm = Coa	// Does not include seed particle mass
		// Populate the 2D waves
		GasMolecules_matrix = TotalMolecules_Matrix - Coa_eqm
		ParticleMolecules_matrix = Coa_eqm
		PM3D_matrix = ParticleMolecules_matrix[p][q]/nSizeBins	// divide mass among all size bins to start
		TotalMolecules_Matrix = GasMolecules_matrix + ParticleMolecules_Matrix
		// add to total
		ParticleMoleculesSum_matrix += ParticleMolecules_matrix // This is necessary to allow calculation of mole fractions	; adds individual VOCs to total sum, which lives in root
	
	// Initialize Kinetics by populating 3D waves for each VOC
		wave GasMolecules_Time, ParticleMolecules_time, TotalMolecules_time
		GasMolecules_Time[][][0] = GasMolecules_matrix[p][q]
		ParticleMolecules_Time[][][0] = ParticleMolecules_Matrix[p][q]
		TotalMolecules_Time[][][0] = TotalMolecules_Matrix[p][q]
		
	// Determine parameters associated with gas-phase wall loss for each VOC
		wave logCstar_matrix, O_matrix
		SOM_GasPhaseWallLoss(logCstar_matrix,O_matrix,gasWLR,Cw_base=Cwall,kwg_increase_perO=kwall_gas_scaling,Krechmer_Cw=Krechmer_Cw)
		// wave kwg_off		// 2D matrix with composition dependent rate of desorption (1/s) from walls	
		setdatafolder $df_home
	endfor	// end loop over VOCs
	
	// Get mole fractions
	cdf = getdatafolder(1)
	for(idex_VOCs=0;idex_VOCs<nVOCs;idex_VOCs+=1) // start loop over VOCs
		df_VOC = df_base + ":VOC" + num2istr(idex_VOCs+1) // set to VOC specific data folder
		setdatafolder $df_VOC

		wave ParticleMoleFraction_matrix, ParticleMolecules_matrix, PMF3D_matrix, ParticleMoleFraction_time
		wavestats/q ParticleMoleculesSum_matrix // this contains molecules summed over all species; this wave lives in root
		ParticleMoleFraction_matrix = ParticleMolecules_Matrix/V_sum
		PMF3D_matrix[][][] = ParticleMoleFraction_matrix[p][q] 
		ParticleMoleFraction_Time[][][0] = ParticleMoleFraction_Matrix[p][q]
		setdatafolder $df_home
	endfor	// end loop over VOCs
	
	// initial particle diameter (seemingly in a random location...should probably move somewhere better)
   	Dp_time[0] = DpStart
   	// some other variables...
	variable TimeVar
	variable FragVar
	variable StepVar

// 10b. Decide whether to utilize OH wave that has been fit to VOC decay (using linear interpolations)
//	SetOH_from_experiment(krxn_matrix[0][0])
	wave OH_exp	// molecules/cm^3
	wave OH_exp_time	// hours
	if(UseScalingForOxidant==2)
		Interpolate2/T=1/N=200/I=3/Y=OH_exp_Interp/X=TimeW OH_exp_time, OH_exp
		wave OH_exp_Interp
		wavestats/q OH_exp
		OH_exp_Interp = OH_exp_Interp < V_min ? V_min : OH_exp_Interp // do not let interpolated values be < minimum observed value
	endif
	// if O3...
	if(O3_yn==1)
		if(UseScalingForOxidant==2)
			OH_wave = OH_exp_interp*1e-9 * (1e-6*NumberDensityAir)
		else
			OH_wave = O3_conc * 1e-9 * (1e-6*NumberDensityAir)
		endif	
	endif
	
	if(quiet==0)
//		print "Ctot = " + num2str(ctot_init) + " in ug/m^3 for Nc = " + num2str(nC) + " with FragSlope = " + num2str(FragSlope) + " and dlVP = " + num2str(delta_logCstar_perO)+ " and ProbOx = " + num2str(ProbOx1) + ";"+ num2str(ProbOx2) + ";"+ num2str(ProbOx3) + ";"+ num2str(ProbOx4) + ";"
	endif
// End 1. Decision Finish
// 10d. Some timing stuff		
	t11 = ticks
	variable FirstTime = 1
	variable counter = 0

// 11. START KINETICS
	for(m=1;m<=nsteps-1;m+=1)
		// 1: Get current OH or O3 concentration (since it may not be constant)	
		if(OH_profile == 1)
			if(timestep*OH_counter > 43200)
				OH_counter = 0
				day_counter += 1
				// check for is even
				day_iseven = day_counter/2
				if(day_iseven == floor(day_iseven))
					day_iseven = 1 // 1 = true = daytime; 0 = false = nighttime
				else
					day_iseven = 0
				endif
			endif
			
			OHConc = (-1*(timestep*OH_counter-21600)^2 + 21600^2) / (21600^2)
			OHConc = (OH_max*(1/((day_counter+1))) - OH_min)*OHConc + OH_min
			OH_counter+=1
			OH_wave[m-1] = OHconc
			if(m>250)
				OH_wave[m-1] = 0
			endif
		elseif(UseScalingForOxidant == 1)
			OH_wave[m-1] = SOM_SetOHconc(timew,OHconc,OH_scale,m)
		elseif(UseScalingForOxidant == 2 && O3_yn!=1)
			OH_wave[m-1] = OH_exp_Interp[m-1]//OH_exp[m-1]
		elseif(O3_yn!=1)	// O3!!
			// do nothing
		endif
		// run differently if O3 is reactant...in this case OH_conc is actually [O3]
		if(O3_yn==5 && UseScalingForOxidant!=2)
			SOM_setO3conc(OH_wave,timew,Tvar,m,O3_conc)
		elseif(O3_yn==1 && UseScalingForOxidant==2)
			 // good to go // OH_wave[m-1] = OH_exp_interp[m-1]*1e-9 * (1e-6*NumberDensityAir)	// molecules/cm^3
		elseif(O3_yn==1 && UseScalingForOxidant==1)
			OH_wave[m-1] =  O3_conc * 1e-9 * (1e-6*NumberDensityAir)*(OH_scale)^(m-1)
		endif

	// 2: Particle loss to the walls; split for multiple compound simulations
		if(ParticleWallLoss!=0)
			Diameter = Diameter_time[m-1][p]
			particleWLR_wave[] = SOM_ParticleWallLossRate(Diameter,particlewallloss)	// get size dependent wall-loss rate (1/s)
			for(index_size=0;index_size<nSizeBins;index_size+=1)
				particleWLR_wave[index_size] = SOM_ParticleWallLossRate(Diameter_time[m-1][index_size],particlewallloss)	// get size dependent wall-loss rate
			endfor
			NpPerBin_time[m][] = NpPerBin_time[m-1][q] - NpPerBin_time[m-1][q]*particleWLR_wave[q]*timestep				// loss of particles
			ParticleMolecules_Matrix[][][] *= (NpPerBin_Time[m][r] / NpPerBin_time[m-1][r]) // added 8/4/18
		endif	
		// calculate total species concentration after (1) gas-phase reaction, (2) gas-phase wall loss and (3) particle-phase wall loss
		TotalMolecules_Matrix = GasMolecules_matrix + ParticleMolecules_Matrix // this is based on previous iteration

		cdf = getdatafolder(1)
		for(idex_VOCs=0;idex_VOCs<nVOCs;idex_VOCs+=1) // start loop over VOCs
			df_VOC = df_base + ":VOC" + num2istr(idex_VOCs+1) // the full path to the VOC specific folder
			setdatafolder $df_VOC // set to VOC specific data folder
			// access a bunch of waves for the current VOC in the VOC-specific data folder
			wave TotalMolecules_matrix, GasMolecules_matrix, ParticleMolecules_matrix
			wave RxnStep_Matrix_minus, krxn_matrix
			wave RxnStep_Matrix_plus, RxnStep_matrix_plus1, RxnStep_Matrix_plus2, RxnStep_Matrix_plus3, RxnStep_Matrix_plus4
			wave Frag1_Matrix, RxnStep_Matrix_Frag
			wave Prob_matrix, Prob1_matrix, Prob2_matrix
			wave RxnStep_array_minus, Frag1_array, RxnStep_array_minus
			wave WPM3D_matrix, PM3D_matrix, WallParticleMolec_Matrix
			wave WallMolecules_matrix, WallMolecules_time, kwg_off
			
	// 3: Run gas phase reaction to determine delta[HC] for a given VOC
			// 3a: Chemical LOSS
			if(O3_yn==0)	// reaction with OH
				RxnStep_Matrix_minus[][] = GasMolecules_Matrix[p][q] * krxn_matrix[p][q] * timestep * OH_wave[m-1] // -d[HC] = [CnOm]*k_nm*delta_t*[OH]
				RxnStep_Matrix_minus[][nO-1] = 0 // needed as limit on rxn...CO2 is unreactive
				RxnStep_Matrix_minus = (numtype(RxnStep_Matrix_minus)==2 ? 0 : RxnStep_Matrix_minus) // Turn NaN's to zero's' for later math.
				GasMolecules_matrix[][] -= RxnStep_Matrix_minus // account for LOSS due to reaction
			elseif(O3_yn==1)	// reaction with O3; added on 08/30/13 to try and speed things up
				RxnStep_Matrix_minus[0][] = GasMolecules_Matrix[0][q] * krxn_matrix[0][q] * timestep * OH_wave[m-1] // -d[HC] = [CnOm]*k_nm*delta_t*[OH]
				RxnStep_Matrix_minus = (numtype(RxnStep_Matrix_minus)==2 ? 0 : RxnStep_Matrix_minus) // Turn NaN's to zero's' for later math.
				GasMolecules_matrix[0][] -= RxnStep_Matrix_minus[0][q] // account for LOSS due to reaction
			endif
			// 3b: Chemical FORMATION - no fragmentation yet
			variable FragMethTest = 1
			// Access VOC-specific parameters (note: SOMparams lives in root:)
			mfrag = SOMparams[1][idex_VOCs]
			POx1 = SOMparams[3][idex_VOCs]
			POx2 = SOMparams[4][idex_VOCs]
			POx3 = SOMparams[5][idex_VOCs]
			POx4 = SOMparams[6][idex_VOCs]
			ProbOxSum = POx1 + POx2 + POx3 + POx4
			POx1 /= ProbOxSum
			POx2 /= ProbOxSum
			POx3 /= ProbOxSum
			POx4 /= ProbOxSum
			// this is the correct implementation of fragmentation. changed 07/11/13
			if(mfrag!=0 || O3_yn==0)	// if fragmentation is on and this is a reaction with OH
				RxnStep_Matrix_plus1=0; RxnStep_Matrix_plus2=0; RxnStep_Matrix_plus3=0; RxnStep_Matrix_plus4=0
				RxnStep_Matrix_plus1[][1,] = RxnStep_Matrix_minus[p][q-1]*POx1*(1-Frag1_Matrix[p][q-1])
				RxnStep_Matrix_plus2[][2,] = RxnStep_Matrix_minus[p][q-2]*POx2*(1-Frag1_Matrix[p][q-2])
				RxnStep_Matrix_plus2[][nO-1] += (RxnStep_Matrix_minus[p][nO-2])*POx2*(1-Frag1_Matrix[p][nO-2])
				RxnStep_Matrix_plus3[][3,] = RxnStep_Matrix_minus[p][q-3]*POx3*(1-Frag1_Matrix[p][q-3])
				RxnStep_Matrix_plus3[][nO-1] += (RxnStep_Matrix_minus[p][nO-2]*(1-Frag1_matrix[p][nO-2])+RxnStep_Matrix_minus[p][nO-3]*(1-Frag1_matrix[p][nO-3]))*POx3
				RxnStep_Matrix_plus4[][4,] = RxnStep_Matrix_minus[p][q-4]*POx4*(1-Frag1_Matrix[p][q-4])
				RxnStep_Matrix_Plus4[][nO-1] += (RxnStep_Matrix_minus[p][nO-2]*(1-Frag1_matrix[p][nO-2])+RxnStep_Matrix_minus[p][nO-3]*(1-Frag1_matrix[p][nO-3])+RxnStep_Matrix_minus[p][nO-4]*(1-Frag1_matrix[p][nO-4]))*POx4
		 		RxnStep_Matrix_plus[][] = RxnStep_Matrix_plus1 + RxnStep_matrix_plus2 + RxnStep_Matrix_Plus3 + RxnStep_matrix_Plus4
		 	else	// for O3 reactions, added 08/29/13
				RxnStep_Matrix_plus = 0;RxnStep_Matrix_plus1=0; RxnStep_Matrix_plus2=0; RxnStep_Matrix_plus3=0; RxnStep_Matrix_plus4=0
				RxnStep_Matrix_plus1[0][1,] = RxnStep_Matrix_minus[0][q-1]*POx1
				RxnStep_Matrix_plus2[0][2,] = RxnStep_Matrix_minus[0][q-2]*POx2
				RxnStep_Matrix_plus2[0][nO-1] += (RxnStep_Matrix_minus[0][nO-2])*POx2
				RxnStep_Matrix_plus3[0][3,] = RxnStep_Matrix_minus[0][q-3]*POx3
				RxnStep_Matrix_plus3[0][nO-1] += (RxnStep_Matrix_minus[0][nO-2]+RxnStep_Matrix_minus[0][nO-3])*POx3
				RxnStep_Matrix_plus4[0][4,] = RxnStep_Matrix_minus[0][q-4]*POx4
				RxnStep_Matrix_Plus4[0][nO-1] += (RxnStep_Matrix_minus[p][nO-2]+RxnStep_Matrix_minus[0][nO-3]+RxnStep_Matrix_minus[0][nO-4])*POx4
		 		RxnStep_Matrix_plus[0][] = RxnStep_Matrix_plus1[0][q] + RxnStep_matrix_plus2[0][q] + RxnStep_Matrix_Plus3[0][q] + RxnStep_matrix_Plus4[0][q]
			endif
			// 3c: Determine what molecules are formed as a result of fragmentation
			if(mFrag != 0) // only run this if fragmentation indeed occurs.
				RxnStep_Matrix_Frag = 0			
				for(i=0;i<nCmax;i+=1)
					RxnStep_array_Minus[i*nO,(i+1)*nO-1] = RxnStep_Matrix_Minus[i][p-i*nO]
				endfor
				multithread prob_matrix = (Prob1_Matrix+Prob2_Matrix)*Frag1_Array[r]*RxnStep_Array_Minus[r]
				MatrixOp/O RxnStep_Matrix_Frag = sumBeams(prob_matrix)
				GasMolecules_Matrix[][] += RxnStep_Matrix_Frag + RxnStep_Matrix_Plus // New version as of 07/11/13
			else // No fragmentation
				GasMolecules_matrix[][] += RxnStep_Matrix_plus // molecules formed per step, accounting for different numbers of oxygens added per step
			endif
			
	// 3d: Wall loss of gas-phase species; added 7/14/2013, fixed on 8/20/2013
			if(gasWLR != 0)
				if(gasWLmethod == 0) // irreversible uptake
					GasMolecules_matrix[][] -=  GasMolecules_matrix[p][q] * gasWLR * timestep 
				else // reversible partitioning
					// wall deposition
					RxnStep_Matrix_minus[][] = GasMolecules_matrix[p][q] * gasWLR * timestep	// gas-phase wall loss
					RxnStep_Matrix_minus = RxnStep_Matrix_minus > GasMolecules_matrix ? GasMolecules_matrix : RxnStep_Matrix_minus
					GasMolecules_matrix[][] -= RxnStep_Matrix_minus
					WallMolecules_Matrix[][] += RxnStep_Matrix_minus
					// wall desorption
					kwg_off = kwg_off*timestep > 1 ? 1/timestep : kwg_off
					RxnStep_Matrix_plus[][] = WallMolecules_matrix[p][q] * kwg_off * timestep	// desorption from walls
					GasMolecules_matrix[][] += RxnStep_matrix_Plus[p][q] 
					WallMolecules_matrix[][] -= RxnStep_Matrix_Plus[p][q]
					WallMolecules_Time[][][m] = WallMolecules_Matrix[p][q]
				endif
			endif
			
	// 3e: Particle loss to the walls; added 08/20/2013; updated 09/29/14
			// This may not be working correctly for multi-component simulations
			if(ParticleWallLoss!=0)
				WPM3D_Matrix[][][] += ((NpPerBin_time[m-1][r]-NpPerBin_time[m][r])/NpPerBin_time[m-1][r])*PM3D_Matrix[p][q][r]	// wall-bound particles, binned
				PM3D_Matrix[][][] -= ((NpPerBin_time[m-1][r]-NpPerBin_time[m][r])/NpPerBin_time[m-1][r])*PM3D_Matrix[p][q][r]		// suspended particles, binned
				MatrixOp/O WallParticleMolec_Matrix = sumBeams(WPM3D_Matrix)	// total wall-bound material
				MatrixOp/O ParticleMolecules_Matrix = sumBeams(PM3D_Matrix)		// total suspended material
			endif	
			// calculate total species concentration after (1) gas-phase reaction, (2) gas-phase wall loss and (3) particle-phase wall loss
			TotalMolecules_Matrix = GasMolecules_matrix + ParticleMolecules_Matrix // Gas is new; Particle is from previous step, adjusted for wall losses

			setdatafolder $df_home // returns to the default folder (aka root:)
		endfor // end loop over VOCs for gas-phase kinetics and wall loss
		
	// 4: Mass Transfer to Particle Phase - Part 1; need to split calculations for multiple VOCs
		if(KineticMassTransfer==0)
			// Eqm Partitioning Calculation
			if(SPM == 0) // as normal
				wavestats/q ParticleMoleculesSum_Matrix // this lives in root:
				// This function accesses data folders for you
				EqmCalc_MultipleVOC(SOMparams,Coa_guess=V_sum,SeedConc=SeedMolecules,SeedMW=SeedMW,units="molecules",Tvar=Tvar)
			else
				// $$$ sequential partitioning model...need to rewrite for multiple VOCs
				wavestats/q ParticleMoleculesSum_Matrix
				EqmCalc_MultipleVOC(SOMparams,Coa_guess=V_sum/m,SeedConc=SeedMolecules,SeedMW=SeedMW,units="molecules",Tvar=Tvar,SPM=1)
			endif
		else
			 // Dynamic partitioning, assuming absorptive partitioning
			 // 1. Do some "precalculations" of things that don't change substantially over a time-step
			for(idex_VOCs=0;idex_VOCs<nVOCs;idex_VOCs+=1) // start loop over VOCs
				df_VOC = df_base + ":VOC" + num2istr(idex_VOCs+1) // set to VOC specific data folder
				setdatafolder $df_VOC
				// access VOC-specific waves
				wave RxnStepKinetic3D_Pre2, SatConc_Matrix
				
				if(idex_VOCs==0)
					wave Knudsen_Matrix, MeanFreePath_Matrix, Beta_Matrix, RxnStepKinetic3D_Pre1, Diffusivity_Matrix // access these from the first VOC folder b/c they are the same for all species
					for(index_size=0;index_size<nSizeBins;index_size+=1)
		 				Knudsen_Matrix = 2*MeanFreePath_Matrix/(Diameter_time[m-1][index_size]*1e-9)	// Fuchs, assuming Delta=lambda; see S&P p. 602
						Beta_Matrix = (0.75*alpha*(1+Knudsen_Matrix))/(Knudsen_Matrix^2+Knudsen_Matrix+0.283*Knudsen_Matrix*alpha+0.75*alpha) // Fuchs-Sutugin
						RxnStepKinetic3D_Pre1[][][index_size] = 4*pi*(Diameter_time[m-1][index_size]*1e-9/2)*Diffusivity_Matrix[p][q]*(timestep/MassTransferTimeScaling)*NpPerBin_time[m][index_size]*Beta_Matrix[p][q]*1e6
					endfor
				endif
				RxnStepKinetic3D_Pre2[][][] = RxnStepKinetic3D_Pre1[p][q][r]*SatConc_Matrix[p][q]*1e-6 // the saturation concentration wave is specific to each VOC
				setdatafolder $df_home
			endfor // end loop over VOCs
				
			// Now, transfer some mass for each VOC; Separate Condensation and Evaporation steps	 			
			for(k=0;k<MassTransferTimeScaling;k+=1)
				PM3Dsum_matrix = 0 // reset to zero
				// first loop around VOCs to get total mass summed over all VOCs for each size bin for calculation of mole or mass fraction
				// Also calculate change due to CONDENSATION
				// (This should probably be checked over...)
				for(idex_VOCs=0;idex_VOCs<nVOCs;idex_VOCs+=1) // start loop over VOCs
					df_VOC = df_base + ":VOC" + num2istr(idex_VOCs+1) // set to VOC specific data folder
					setdatafolder $df_VOC
					wave RxnStepKinetic3D_Cond, GasMolecules_Matrix, SumKinetics, PM3D_Matrix
	 				// Allow condensation, then calculate mass fraction
	 				RxnStepKinetic3D_Cond[][][] = RxnStepKinetic3D_Pre1[p][q][r]*GasMolecules_Matrix[p][q] // determine condensation amount; molecules/bin (delta t already accounted for)
	 				MatrixOp/O SumKinetics = sumBeams(RxnStepKinetic3D_Cond) // sum over particle size bins
	 				RxnStepKinetic3D_Cond[][][] = SumKinetics[p][q] > GasMolecules_Matrix[p][q] ? 0 : RxnStepKinetic3D_Cond[p][q][r] // conditional to set things to zero if too much mass is transferred
	 				PM3D_Matrix[][][] += RxnStepKinetic3D_Cond[p][q][r] // add condensed molecules to particle phase across size bins; molecules/bin
	 				GasMolecules_matrix -= SumKinetics[p][q] // loss of gas-phase mass for each VOC due to condensation
	 				PM3Dsum_matrix += PM3D_Matrix // summing over all VOCs (this is wave PM3D_matrix in root)
					setdatafolder $df_home
				endfor // end loop over VOCs
				// second loop around VOCs to get mole fractions, etc., and calculate EVAPORATION
				for(idex_VOCs=0;idex_VOCs<nVOCs;idex_VOCs+=1) // start loop over VOCs
					df_VOC = df_base + ":VOC" + num2istr(idex_VOCs+1) // set to VOC specific data folder
					setdatafolder $df_VOC
					wave MoleFraction3D, PM3D_matrix, MW_matrix, PMM3D_matrix, ParticleMass_matrix
					wave RxnStepKinetic3D_Evap, RxnStepKinetic3D_Pre2, SumKinetics
					wave GasMolecules_matrix, ParticleMolecules_matrix
					if(idex_VOCs==0)
						PMM3Dsum_matrix[][][] = PM3Dsum_matrix[p][q][r]*MW_matrix[p][q] // This is the total, mass weighted molecules in each size bin; weight by molecular weight to use mass fraction; this lives in root:
					endif
					PMM3D_matrix[][][] = PM3D_Matrix[p][q][r]*MW_matrix[p][q] // This is the VOC-specific mass-weighted molecules in each size bin; weight by molecular weight to use mass fraction
	 				// loop over size bins
	 				for(index_size=0;index_size<nSizeBins;index_size+=1) // need to normalize each size bin individually, so loop through bins
	 					ParticleMass_Matrix = PMM3D_matrix[p][q][index_size] // get mass matrix for each VOC for a given size bin
	 					ParticleMassSum_Matrix = PMM3Dsum_matrix[p][q][index_size] // get total mass matrix (summed over VOCs) for a given size bin
	 					wavestats/q ParticleMassSum_Matrix // total mass in a given size bin, summed over all SOM species
	 					if(V_sum==0) // there are problems when there are zero molecules in the condensed phase
							MoleFraction3D[][][index_size] = numtype(ParticleMass_Matrix[p][q])==2 ? nan : 1 // deal with pesky nans (mostly important for first step)
						else 
							if(absorbingseed==0) // non-absorbing seed
								MoleFraction3D[][][index_size] = ParticleMass_Matrix[p][q]/(V_sum) // Divide each species concentration by the total concentration for each size bin
							else // absorbing seed
								MoleFraction3D[][][index_size] = ParticleMass_Matrix[p][q]/(V_sum+MoleculesSeedPerBin[index_size]*seedMW) // Divide each species concentration by the total concentration (including seed) for each size bin
							endif
						endif
					endfor // end loop over size bins
					// Calculate change due to EVAPORATION
	 				RxnStepKinetic3D_Evap[][][] = RxnStepKinetic3D_Pre2[p][q][r]*MoleFraction3D[p][q][r] // evaporation amount; molecules/particle
	 				RxnStepKinetic3D_Evap[][][] = RxnStepKinetic3D_Evap[p][q][r] > PM3D_Matrix[p][q][r] ? PM3D_Matrix[p][q][r] : RxnStepKinetic3D_Evap[p][q][r] // conditional
	 				MatrixOp/O SumKinetics = SumBeams(RxnStepKinetic3D_Evap) // for adding to the gas phase; sum across all particle bins
	 				GasMolecules_matrix[][] += SumKinetics // update gas-phase matrix for each VOC due to evaporation
	 				PM3D_matrix[][][] -= RxnStepKinetic3D_Evap[p][q][r] // update total particle mass for each VOC in each size bin
	 			endfor // end loop over VOCs
			endfor // end loop over mass transfer time
			
			// HETEROGENEOUS CHEMISTRY
			if(hetchem==1)
				for(idex_VOCs=0;idex_VOCs<nVOCs;idex_VOCs+=1) // start loop over VOCs
					df_VOC = df_base + ":VOC" + num2istr(idex_VOCs+1) // set to VOC specific data folder
					setdatafolder $df_VOC
					// access VOC-specific waves
					wave PM3D_matrix, MoleFraction3D, Frag1_matrix
					wave HetChem_matrix_minus, HetChem_matrix_plus
					
					if(idex_VOCs==0)
						wave Knudsen_Matrix, MeanFreePath_Matrix, Beta_Matrix, RxnStepKinetic3D_Pre1, Diffusivity_Matrix // access these from the first VOC folder b/c they are the same for all species
						for(index_size=0;index_size<nSizeBins;index_size+=1)
			 				Knudsen_Matrix = 2*MFP_OH/(Diameter_time[m-1][index_size]*1e-9)	// Fuchs, assuming Delta=lambda; see S&P p. 602
							Beta_Matrix = (0.75*(1+Knudsen_Matrix))/(Knudsen_Matrix^2+Knudsen_Matrix+0.283*Knudsen_Matrix+0.75) // Fuchs-Sutugin
							RxnStepKinetic3D_Pre1[][][index_size] = (gammaOH*NpPerBin_time[m][index_size]*pi*(Diameter_time[m-1][index_size]*1e-9)^2*Beta_Matrix*OH_wave[m-1]*1e6*1.381e-23*Tvar)/sqrt(2*pi* (17/1000/6.022e23)*1.381e-23*Tvar)
						endfor
					endif
					//for(index_size=0;index_size<nSizeBins;index_size+=1)
					HetChem_matrix_minus = RxnStepKinetic3D_Pre1[p][q][r]*MoleFraction3D[p][q][r]*timestep // the mole fractions are calculated in the previous for loop for evaporation
					HetChem_matrix_minus[][nO-1][] = 0 // limits reaction to max oxygens
					HetChem_matrix_minus[nCmax-1][2] = 0 // CO2 is unreactive	
					HetChem_Matrix_minus = (numtype(HetChem_Matrix_minus)==2 ? 0 : HetChem_Matrix_minus) // Turn NaN's to zero's' for later math.
					HetChem_Matrix_minus = HetChem_matrix_minus > PM3D_matrix ? PM3D_matrix : HetChem_Matrix_Minus
					HetChem_matrix_plus = 0
					HetChem_Matrix_Plus[][1,][] = HetChem_Matrix_Minus[p][q-1][r]*(1-Frag1_Matrix[p][q-1])
					HetChem_Matrix_plus = (numtype(HetChem_Matrix_plus)==2 ? 0 : HetChem_Matrix_plus)
					PM3D_matrix += (HetChem_matrix_plus - HetChem_matrix_minus)
					// for now, just assume that all fragments are highly volatile and do not track (mass is not conserved)
					setdatafolder $df_home
				endfor // end loop over VOCs
			endif
			
			// Sum over all bins for each VOC to get total condensed mass for that VOC by 'species' and the total overall mass
			ParticleMoleculesSum_Matrix = 0
			PM3Dsum_matrix = 0 // added 8/4/18
			for(idex_VOCs=0;idex_VOCs<nVOCs;idex_VOCs+=1) // start loop over VOCs
				df_VOC = df_base + ":VOC" + num2istr(idex_VOCs+1) // set to VOC specific data folder
				setdatafolder $df_VOC
				wave ParticleMolecules_matrix, PM3D_matrix		
	 			MatrixOp/O ParticleMolecules_Matrix = sumBeams(PM3D_Matrix)	// sum over size bins (layers) to get total condensed mass for a given VOC
	 			ParticleMoleculesSum_Matrix += ParticleMolecules_Matrix // sum over all VOCs to get total condensed mass
	 			PM3Dsum_matrix += PM3D_Matrix // added 8/4/18
	 		endfor
		endif

		// 11f: Mass Transfer to Particle Phase - Part 2
		// This is only important for instantaneous equilibrium (and a bit of an awkward way to do this)
		setdatafolder $df_home
		ParticleMoleculesSum_Matrix = 0 // reset value to zero
		for(idex_VOCs=0;idex_VOCs<nVOCs;idex_VOCs+=1) // start loop over VOCs
			df_VOC = df_base + ":VOC" + num2istr(idex_VOCs+1) // set to VOC specific data folder
			setdatafolder $df_VOC
			if(KineticMassTransfer==0)
				// Eqm Partitioning Calculation --> move molecules around
				if(SPM == 0) // as normal
					wave Coa_eqm = Coa // does not include seed particle mass
					wave TotalMolecules_Matrix, GasMolecules_Matrix, ParticleMolecules_Matrix
					if(NoSOA==0)	// only operate if you want to include particle-gas partitioning...which is the default case
						GasMolecules_matrix = TotalMolecules_Matrix - Coa_eqm
						ParticleMolecules_matrix = Coa_eqm
					endif
				else
					// $$$ sequential partitioning model...need to rewrite for multiple VOCs
					wave Coa_eqm = Coa // does not include seed particle mass
					wave TotalMolecules_Matrix, GasMolecules_Matrix, ParticleMolecules_Matrix
					if(NoSOA==0)	// only operate if you want to include particle-gas partitioning...which is the default case
						GasMolecules_matrix = TotalMolecules_Matrix - Coa_eqm
						ParticleMolecules_matrix = Coa_eqm
					endif
				endif
			else
				 // Dynamic partitioning, assuming absorptive partitioning; this seems to be repetitive with above
			endif
//			ParticleMoleculesSum_matrix += ParticleMolecules_matrix
			setdatafolder $df_home
		endfor // end loop over VOCs

// 11ff: Oligomerization reactions
// Oligomerization is not ready.
//		if(OligomerizationIsOn==1)
//			SOM_Oligomerization(ParticleMolecules_Matrix,krxn_matrix_Olig,timestep,Np,VolumeOrgPerParticle,FirstTime=FirstTime,method=2)
//			wave ParticleMolecules_Matrix
//		endif

// 11fff: Final concentrations from this iteration
		TotalCarbonAtoms = 0
		TotalOxygenAtoms = 0
		TotalHydrogenAtoms = 0

		// dilution, added by CDC on 07.05.17
		if(DilutionVar != 0)
			SeedMolecules *= 1-(DilutionVar/100)*(TimeStep/60/60)
			MoleculesSeedPerBin *= 1-(DilutionVar/100)*(TimeStep/60/60)
			SeedConc_time[m] = SeedConc_time[m-1]*(1-(DilutionVar/100)*(TimeStep/60/60))
			
			for(idex_VOCs=0;idex_VOCs<nVOCs;idex_VOCs+=1) // start loop over VOCs
				df_VOC = df_base + ":VOC" + num2istr(idex_VOCs+1) // set to VOC specific data folder
				setdatafolder $df_VOC
				wave TotalMolecules_matrix, GasMolecules_matrix, ParticleMolecules_matrix, ParticleMass_matrix, WallParticleMolec_Matrix
				wave TotalMolecules_time, GasMolecules_time, ParticleMolecules_time
				wave Coa_time, Coa_time_wallcorr, MW_matrix
				wave C_matrix, CarbonAtoms, O_matrix, OxygenAtoms, H_matrix, HydrogenAtoms
				wave O2C_time, H2C_time
				wave deltaHC_time, Yield_time, HC_ugm3_time
				
				// Populate time-dependent waves
				GasMolecules_matrix *= 1-(DilutionVar/100)*(TimeStep/60/60)
				ParticleMolecules_Matrix *= 1-(DilutionVar/100)*(TimeStep/60/60)
				setdatafolder $df_home
			endfor // end loop over VOCs	
		endif
	
		for(idex_VOCs=0;idex_VOCs<nVOCs;idex_VOCs+=1) // start loop over VOCs
			df_VOC = df_base + ":VOC" + num2istr(idex_VOCs+1) // set to VOC specific data folder
			setdatafolder $df_VOC
			wave TotalMolecules_matrix, GasMolecules_matrix, ParticleMolecules_matrix, ParticleMass_matrix, WallParticleMolec_Matrix
			wave TotalMolecules_time, GasMolecules_time, ParticleMolecules_time
			wave Coa_time, Coa_time_wallcorr, MW_matrix
			wave C_matrix, CarbonAtoms, O_matrix, OxygenAtoms, H_matrix, HydrogenAtoms
			wave O2C_time, H2C_time
			wave deltaHC_time, Yield_time, HC_ugm3_time
			
			// Populate time-dependent waves
			TotalMolecules_Matrix = GasMolecules_matrix + ParticleMolecules_Matrix
			GasMolecules_Time[][][m] = GasMolecules_Matrix[p][q]
			ParticleMolecules_Time[][][m] = ParticleMolecules_Matrix[p][q]
			TotalMolecules_Time[][][m] = TotalMolecules_Matrix[p][q]
			ParticleMoleculesSum_matrix[][] += ParticleMolecules_Matrix[p][q]
			
			// Suspended SOA concentration
			ParticleMass_Matrix = ParticleMolecules_Matrix[p][q] * 1e6 * 1e6 * MW_Matrix[p][q]/Na
			wavestats/q ParticleMass_matrix
			Coa_time[m] = V_sum // Current SOA concentration from a given VOC (ug/m^3)
			CoaSum_time[m] += V_sum // this was initialized to zero earlier

			// Correct for particle deposition
			ParticleMass_Matrix = (ParticleMolecules_Matrix[p][q] + WallParticleMolec_Matrix[p][q]) * 1e6 * 1e6 * MW_Matrix[p][q]/Na
			wavestats/q ParticleMass_Matrix
			Coa_time_wallcorr[m] = V_sum // Current OA concentration, corrected for particle wall-loss (ug/m^3)

			// 11g: O:C calculation
			// VOC-specific	
			wavestats/Q ParticleMolecules_matrix
			CarbonAtoms = C_matrix*ParticleMolecules_matrix//V_Sum
			OxygenAtoms = O_matrix*ParticleMolecules_matrix//V_Sum
			HydrogenAtoms = H_matrix*ParticleMolecules_matrix//V_sum
			wavestats/Q CarbonAtoms
			Ave_nC = V_Sum
			wavestats/Q OxygenAtoms
			Ave_nO = V_sum
			wavestats/Q HydrogenAtoms
			Ave_nH = V_sum
			O2C_time[m] = Ave_No/Ave_Nc
			H2C_time[m] = Ave_Nh/Ave_Nc	
			// Sum over all VOCs		
			TotalCarbonAtoms += CarbonAtoms
			TotalOxygenAtoms += OxygenAtoms
			TotalHydrogenAtoms += HydrogenAtoms
			
			// Yields and VOC reacted
			GasMolecules_Parent = GasMolecules_time[p][0][0] * 1e6 * 1e6 * MW_Matrix[p][0]/Na
			wavestats/q GasMolecules_Parent
			HCstart = V_sum
			if(m==1)
				HC_ugm3_time[0] = V_sum
				HCsum_ugm3_time[0] += V_sum
			endif
			GasMolecules_Parent = GasMolecules_time[p][0][m] * 1e6 * 1e6 * MW_Matrix[p][0]/Na
			wavestats/q GasMolecules_Parent
			HCfinish = V_sum
			HC_ugm3_time[m] = V_sum
			HCsum_ugm3_time[m] += V_sum
			deltaHC_time[m] = HCstart - HCfinish // ug/m3
			Yield_time[m] = Coa_time[m]/deltaHC_time[m]
			deltaHCsum_time[m] += deltaHC_time[m]
			setdatafolder $df_home
		endfor // end loop over VOCs	

		ParticleMoleculesSum_Time[][][m] = ParticleMoleculesSum_Matrix[p][q]
		wavestats/Q TotalCarbonAtoms
		Ave_nC = V_Sum
		wavestats/Q TotalOxygenAtoms
		Ave_nO = V_sum
		wavestats/Q TotalHydrogenAtoms
		Ave_nH = V_sum
		O2Csum_time[m] = Ave_nO/Ave_nC
		H2Csum_time[m] = Ave_nH/Ave_nC		
		
		// loop over sizes to get new diameters...this may not be quite perfect yet (09/29/14)
		variable VolOrgTot = 0
		if(KineticMassTransfer==0) // eqm. partitioning)
			PM3D_matrix = 0
			PM3D_matrix[][][0] = ParticleMolecules_Matrix[p][q]
		endif
		for(index_size=0;index_size<nSizeBins;index_size+=1)
			ParticleVolume_matrix = (PM3D_matrix[p][q][index_size]/NpPerBin_time[m][index_size])*(MW_Matrix[p][q]/(Na*Density_Base*1e6))
			wavestats/Q ParticleVolume_Matrix
			VolumeOrgPerParticle = V_sum // m^3/p
			Diameter_time[m][index_size] = 1e9*((6/pi)*(VolumeOrgPerParticle+VolumeSeedPerBin[index_size]))^(1/3)
			VolumeOrgPerBin[index_size] = VolumeOrgPerParticle
			VolOrgTot += VolumeOrgPerParticle
		endfor	
		MatrixOp/O Np_time = sumRows(NpPerBin_time)	// total number concentration of particles, summed over all bins
		Dp_time[m] = 1e9*((6/pi)*(VolOrgTot+VolumeSeed))^(1/3)	// for polydisperse simulations, this is an overall "average" diameter

		if(m==100 && quiet == 0)
			t22 = ticks
			print "Calculations will likely take ~ " + num2str(((t22-t11)/60)*nsteps/100) + " seconds"
		endif
		FirstTime = 0

	endfor // end big kinetics loop
	
	GasMassSum_Matrix = 0
	for(idex_VOCs=0;idex_VOCs<nVOCs;idex_VOCs+=1) // start loop over VOCs
		df_VOC = df_base + ":VOC" + num2istr(idex_VOCs+1) // set to VOC specific data folder
		setdatafolder $df_VOC
		wave GasMolecules_matrix, MW_matrix		
		GasMassSum_Matrix += GasMolecules_Matrix  * 1e6 * 1e6 * MW_Matrix[p][0]/Na// sum over all VOCs to get total condensed mass
 		endfor

	
	// YIELDS
	YieldSum_time = CoaSum_time/deltaHCsum_time
	Diameter_final = Diameter_time[m-1][p]
	dNdlogDp_final = dNdlogDp_time[m-1][p]
	
	// INTERPOLATE
	// Select whether to interpolate the results to an experimental time base
	// This is important when data fitting
	// The wave Time_experiment must exist and this only works if the simulation was run for longer than the experiment
	variable InterpToExperiment = 0 // 0 = do nothing; 1 = interpolate
	if(InterpToExperiment==1)
		wave Time_experiment
		Interpolate2/T=2/E=1/I=3/Y=Coa_time_interp/X=Time_experiment Coa_time
	endif
	
	// Print results
	if(quiet==0)
		print "[SOA] =  " + num2str(CoaSum_time[nsteps-1]) + " with Yield = " + num2str(YieldSum_time[nsteps-1]) + " and average O:C = " + num2str(O2Csum_time[nsteps-1])
	endif
	
	SOM_KillWaves()
	SOM_KillWaves_MultiComponent()
	
	Graph_MC_yield()
	Graph_MC_GandP()

END

//*************************************************************************************************************************
// Calculate gas-particle equilbirum based on "Pankow" theory (i.e. Raoult's Law)
// Written by CDC on 10/24/15
// This is based on EqmCalc(), but extended to work for multiple VOCs
Function EqmCalc_MultipleVOC(SOMparams,[Coa_guess,SeedConc,SeedMW,units,Tvar,SPM])
	wave SOMparams // wave containing info on VOCs
	variable Coa_guess	// initial guess, in molecules/cm^3 or ug/m^3
	variable SeedConc // seed concentration in molecules/cm^3 or ug/m^3
	variable SeedMW // seed MW, if absorbing
	string units	// either "molecules" or "mass". "molecules" is default
	variable Tvar	// temperature, in K --> added 09/03/13
	variable SPM // sequential partitioning model --> added 02/28/16

	if(ParamIsDefault(units))
		units = "molecules"
	endif
	if(ParamIsDefault(SeedConc))
		seedconc = 0
	endif
	if(ParamIsDefault(SeedMW))
		SeedMW = 400
	endif
	if(ParamIsDefault(Coa_Guess))
//		wavestats/Q ConcMatrix
		Coa_Guess = 1 //V_sum/4
	endif
	variable Tadjust
	if(ParamIsDefault(Tvar))	// added 09/03/13
		Tadjust = 0
	else
		Tadjust = 1
	endif
	if(ParamIsDefault(SPM))
		SPM = 0
	endif
	
	// prep and deal with seed
	variable Coa_tol	// tolerance on result
	variable Seed
	NVAR AbsorbingSeed = $(ksDFroot + "absorbingseed")
	if(stringmatch(units,"molecules"))
		Coa_tol = 1	// tolerance is 1 molecule
		Seed = SeedConc * 1e6 // convert molecules/cm^3 to molecules/m^3
	elseif(stringmatch(units,"mass"))
		Coa_tol = 1e-6	// tolerance is 1e-6 ug/m^3
		Seed = SeedConc * 1e6 * 1e6 * SeedMW / 6.022e23
	endif
	if(AbsorbingSeed==1)	
		Coa_Guess = seed
	endif
	Variable Coa_calc
	Variable Coa_dif
	variable counter=0	// added 08/31/13

	// deal with mulitple VOCs
	string df_home = ksDFroot // "root:"
	string df_base = ksDFVOCbase // "root:VOCs" // Folder in which waves for individual VOCs are stored (in subfolders)
	string df_VOC // TBD; subfolders for storing individual waves from individual VOCs
	setdatafolder $df_home // root:
	string cdf = getdatafolder(1)
	
	variable nVOCs = dimsize(SOMparams,1) // # of columns = # of precursor VOCs
	variable idex_VOCs // index for looping around different VOCs
	variable nC, nO
	string str_ConcMatrix, str_ConcSat25CforTD, str_MWwave, str_deltaHvap_matrix
	
	for(idex_VOCs=0;idex_VOCs<nVOCs;idex_VOCs+=1) // loop over individual VOCs
		df_VOC = df_base + ":VOC" + num2istr(idex_VOCs+1)
		setdatafolder $df_VOC
		// access waves
		if(SPM==0) // all material can partition
			wave ConcMatrix = TotalMolecules_Matrix // in molecules/cm^3 
		else // sequential partitioning model
			wave ConcMatrix = GasMolecules_matrix // in molecules/cm^3
		endif
		wave ConcSat25CForTD = SatConc_matrix // saturation vapor pressure in molecules/m^3
		wave MWwave = MW_matrix // molecular weight matrix in g/mol
		wave deltaHvap_matrix	// matrix of deltaHvap values, in kJ/mol --> added 09/03/13
		nC = dimsize(ConcMatrix,0)
		nO = dimsize(ConcMatrix,1)
		// create some new waves in the VOC folders
		make/o/d/n=(nC,nO) ConcWave=0, FracAero=0, Coa=0
		make/o/d/n=(nC,nO) SatConc_temp=0

		if(idex_VOCs==0) 
			duplicate/o/d Coa $(df_home + "Coa") // root:Coa // make a "total" wave in the root folder
			wave CoaSum = $(df_home + "Coa") // root:Coa
			CoaSum = 0
		endif
		
		if(stringmatch(units,"molecules"))
			ConcWave = ConcMatrix * 1e6		// convert molecules/cm^3 to molecules/m^3
			ConcWave = numtype(ConcWave)==2 ? 0 : ConcWave
			SatConc_temp = ConcSat25CforTD	// no units change
		elseif(stringmatch(units,"mass"))
			ConcWave = ConcMatrix * 1e6 * 1e6 * MWwave / 6.022e23	// convert molecules/cm^3 to ug/m^3
			ConcWave = numtype(ConcWave)==2 ? 0 : ConcWave
			SatConc_temp = ConcSat25CforTD * 1e6 * MWwave / 6.022e23	// convert molecules/m^3 to ug/m^3
		endif	
	
		if(Tadjust==1)
			// adjust vapor pressure based on Tvar and deltaHvap_matrix --> added 09/03/13
			// assumes ConcSat25CForTD at 25C
			variable Tref = 298
			if(stringmatch(units,"molecules"))
				SatConc_temp[][] = SatConc_temp[p][q]*exp(-1*(deltaHvap_matrix[p][q]*1000/8.314)*(1/Tvar - 1/Tref))
			elseif(stringmatch(units,"mass"))
				SatConc_temp[][] = SatConc_temp[p][q]*(Tref/Tvar)*exp(-1*(deltaHvap_matrix[p][q]*1000/8.314)*(1/Tvar - 1/Tref))
			endif
		endif

		WaveStats/Q ConcWave
		Coa_guess += V_sum/4 // keep adding for each VOC

		setdatafolder $df_home // back to start folder
	endfor // end multiple VOCs (for now)
	
	do // loop until tolerance is met
		if(AbsorbingSeed==1)
			Coa_calc = seed // reset to seed value every time, if absorbing seed
		endif
		CoaSum = 0
		for(idex_VOCs=0;idex_VOCs<nVOCs;idex_VOCs+=1) // loop over VOCs
			// go to a specific VOC folder
			df_VOC = df_base + ":VOC" + num2istr(idex_VOCs+1)
			setdatafolder $df_VOC
			
			wave ConcWave
			wave FracAero
			wave SatConc_temp
			wave Coa
			FracAero = 1/(1+(SatConc_temp)/Coa_guess) // particle-phase fraction for a given VOC (partitioning equation)
			Coa = FracAero * ConcWave // Coa for a given VOC
			wavestats/Q Coa // get total amount
			Coa_calc += V_Sum // accumulate for each VOC
			CoaSum += Coa // accumulate for each VOC
		endfor
		Coa_dif = abs(Coa_calc - Coa_guess) // compare actual with guess value
		if(Coa_dif < Coa_tol) // compare to tolerance
			break
		endif
		Coa_guess = Coa_calc // reset "guess" for next iteration
		counter += 1 // update counter
		if(counter==5)	// reset tolerance based on current guess; this speeds things up since the desired tolerance depends on the amount of OA there is.
			Coa_tol = max(Coa_tol,0.01*(Coa_calc-Seed))
		endif
	while(1)

	CoaSum = 0 // reset 
	for(idex_VOCs=0;idex_VOCs<nVOCs;idex_VOCs+=1) // loop over VOCs
		// go to a specific VOC folder
		df_VOC = df_base + ":VOC" + num2istr(idex_VOCs+1)
		setdatafolder $df_VOC
		wave Coa
		wave MWwave = MW_matrix
		if(stringmatch(units,"molecules"))
			Coa = Coa * 1e-6 // convert molecules/m^3 to molecules/cm^3...does not include seed mass
			CoaSum += Coa
		elseif(stringmatch(units,"mass"))
			Coa = Coa * 6.022e23 / (1e6 * 1e6 * MWwave)	// convert ug/m^3 to molecules/cm^3
			CoaSum += Coa
		endif
	endfor
	setdatafolder $cdf
End	

//*****************************************************************************************************************
Function SOM_KillWaves_MultiComponent()

	string df_base = ksDFVOCbase // "root:VOCs"
	string df_VOC
	variable idex_VOCs
	variable nVOCs
	wave SOMparams = root:SOMparams
	nVOCs = dimsize(SOMparams,1)
	
		for(idex_VOCs=0;idex_VOCs<nVOCs;idex_VOCs+=1) // start loop over VOCs
			df_VOC = df_base + ":VOC" + num2istr(idex_VOCs+1) // set to VOC specific data folder
			setdatafolder $df_VOC
			killwaves/z Beta_matrix
			killwaves/z CarbonAtoms
			killwaves/z ConcWave
			killwaves/z C_matrix
			killwaves/z C_Matrix_small
			killwaves/z C_matrix_temp
//			killwaves/z deltaHC_time
			killwaves/z deltaHCfrac_time
			killwaves/z deltaHvap_matrix
			killwaves/z Diffusivity_matrix
			killwaves/z FracAero
			killwaves/z Frag1_array
			killwaves/z Frag1_array_small
			killwaves/z Frag1_matrix
			killwaves/z Frag1_matrix_small
			killwaves/z GasMass_matrix
			killwaves/z GasMass_time
			killwaves/z GasMolecules_matrix
			killwaves/z GasMolecules_matrix_log
			killwaves/z GasMolecules_time
			killwaves/z HC_matrix
			killwaves/z HC_ppb_time
//			killwaves/z HC_ugm3_time
			killwaves/z HetChem_matrix_minus
			killwaves/z HetChem_matrix_plus
			killwaves/z HydrogenAtoms
			killwaves/z H_matrix
			killwaves/z Knudsen_matrix
			killwaves/z krxn_matrix
			killwaves/z kwg_off
			killwaves/z lifetimes
			killwaves/z logCstar_matrix
			killwaves/z logCstar_mean_time
			killwaves/z MeanFreePath_matrix
			killwaves/z Molecules_time
			killwaves/z MoleFraction3D
			killwaves/z MW_matrix
			killwaves/z Ncarbons_wave
			killwaves/z Nc_ave_time
			killwaves/z Noxygens_wave
			killwaves/z No_ave_time
			killwaves/z OCseed_time
			killwaves/z OC_matrix
			killwaves/z OC_matrix_small
			killwaves/z O_matrix_small
			killwaves/z OxygenAtoms
			killwaves/z Ox_matrix
			killwaves/z O_matrix
			killwaves/z ParticleFrac_matrix
			killwaves/z ParticleMass_matrix
//			killwaves/z ParticleMass_time
			killwaves/z ParticleMolecules_matrix
			killwaves/z ParticleMolecules_time
			killwaves/z ParticleMoleFraction_matrix
			killwaves/z ParticleMoleFraction_time
			killwaves/z ParticleVolume_matrix
			killwaves/z Particle_div_total
			killwaves/z Particle_fraction
			killwaves/z Particle_Fraction_log
			killwaves/z PM3D_matrix
			killwaves/z PMF3D_matrix
			killwaves/z PMM3D_matrix
			killwaves/z Prob1_matrix
			killwaves/z Prob1_matrix_small
			killwaves/z Prob2_matrix
			killwaves Prob2_matrix_small
			killwaves/z Prob_matrix
			killwaves/z Prob_matrix_small
			killwaves/z RxnStepKinetic3D_cond
			killwaves/z RxnStepKinetic3D_evap
			killwaves/z RxnStepKinetic3D_Pre1
			killwaves/z RxnStepKinetic3D_Pre2
			killwaves/z RxnStep_array_minus
			killwaves/z RxnStep_matrix_Frag
			killwaves/z RxnStep_matrix_kinetic_off
			killwaves/z RxnStep_matrix_kinetic_on
			killwaves/z RxnStep_matrix_kinetic_pre1
			killwaves/z RxnStep_matrix_kinetic_pre2
			killwaves/z RxnStep_matrix_kinetic_time
			killwaves/z RxnStep_matrix_minus
			killwaves/z RxnStep_matrix_plus
			killwaves/z RxnStep_matrix_plus1
			killwaves/z RxnStep_matrix_plus2
			killwaves/z RxnStep_matrix_plus3
			killwaves/z RxnStep_matrix_plus4
			killwaves/z SatConc_matrix
			killwaves/z SatConc_temp
			killwaves/z SumKinetics
			killwaves/z TimeW
			killwaves/z TotalMass_matrix
			killwaves/z TotalMass_time
			killwaves/z TotalMolecules_matrix
			killwaves/z TotalMolecules_matrix_log
			killwaves/z TotalMolecules_time
			killwaves/z wallmolecules_matrix
			killwaves/z wallmolecules_time
			killwaves/z wallparticlemolec_matrix
			killwaves/z wpm3d_matrix
			killwaves/z yield2_time
			setdatafolder root:
		endfor // end loop over VOCs

End

//************************************************************************************************************
Function SOM_MakeWaves(SOMparams,nSteps,nSizeBins, nCmax, nO, timestep)
	wave SOMparams// 2D wave; rows = parameters for a given compound; columns = different compounds
					// [0] = Ncarbons for each species
					// [1] = FragSlope
				 	// [2] = delta_logCstar_perO
				 	// [3] = Pfunc[0]
				 	// [4] = Pfunc[1]
				 	// [5] = Pfunc[2]
				 	// [6] = Pfunc[3]
				 	// [7] = krxn (leave 0 if default should be used)
	variable nSteps	// number of time steps
	variable nSizeBins // number of size bins (for polydisperse)
	variable nCmax // max number of carbon atoms based on all VOCs
	variable nO // max number of oxygen atoms
	variable timestep // timestep in seconds
//	variable SeedVolConc // seed volume concentration (PROBABLY NEED TO REMOVE THIS)
	
	string df_root = ksDFroot // "root:"
	string df_base = ksDFVOCbase // "root:VOCs" // base VOC folder
	string df_VOC // VOC-specific folder

	setdatafolder $df_root
	if(datafolderexists(df_base)==0)
		print "yes"
		newdatafolder $df_base
	endif	
	variable nVOCs = dimsize(SOMparams,1) // number of precursor compounds (= # of columns)
	variable i
	variable nC // # of carbon atoms of precursor
	
	variable maxNumberOfSpecies = 10 
	
//	for(i=1;i<=nVOCs;i+=1)
	for(i=1;i<=maxNumberOfSpecies;i+=1)
	
	// Set folders
		setdatafolder $df_base // root:VOCs
		df_VOC = df_base + ":VOC" + num2istr(i)// "root:VOCs:VOCX"
		if(datafolderexists(df_VOC)==0)
			newdatafolder $df_VOC // make new folder as necessary
		endif
		setdatafolder $df_VOC // go to folder for a given VOC
		
	// Get Key parameters for a given VOC
		// create all matrices based on maximum number of carbon atoms...will make things slower, but easier
		nC = nCmax //SOMparams[0][i-1] // first row of this matrix contains number of carbons for a given precursor VOC
		
		if(i <= nVOCs)
			variable mFrag = SOMparams[1][i-1]
			variable DLVP = SOMparams[2][i-1]
			variable POx1 = SOMparams[3][i-1]
			variable POx2 = SOMparams[4][i-1]
			variable POx3 = SOMparams[5][i-1]
			variable POx4 = SOMparams[6][i-1]
			variable ProbOxSum = POx1 + POx2 + POx3 + POx4
			POx1 /= ProbOxSum
			POx2 /= ProbOxSum
			POx3 /= ProbOxSum
			POx4 /= ProbOxSum
			SOMparams[3][i-1] = POx1
			SOMparams[4][i-1] = POx2
			SOMparams[5][i-1] = POx3
			SOMparams[6][i-1] = POx4
		else
			// do nothing, as these compounds don't exist and all we want to do
			// is set them to zero in case they were to exist otherwise
			nc = 1
			no = 1
		endif
		
	// Make some waves
		make/o/d/n=(nC,nO) SatConc_Matrix = 0
		make/o/d/n=(nC,nO) logCstar_matrix = 0
		//
		make/d/o/n=(nC,nO) GasMolecules_Matrix = 0
		note GasMolecules_Matrix "gas-phase species concentration in molecules/cm^3"
		note GasMolecules_Matrix "Ncarbons = Rows; Noxygens = columns"
		make/d/o/n=(nC,nO) GasMolecules_Matrix_Log = 0
		make/d/o/n=(nC,nO) GasMass_Matrix = 0
		// 
		make/d/o/n=(nC,nO) ParticleMolecules_Matrix = 0 // molecules/cm^3
		note ParticleMolecules_Matrix "particle phase species concentration in molecules/cm^3"
		note ParticleMolecules_Matrix "Ncarbons = Rows; Noxygens = columns"
		// create 3D version to deal with size-dependent wall losses; added 9/29/14
		make/d/o/n=(nC,nO,nSizeBins) PM3D_Matrix = 0 // molecules/cm^3
		note PM3D_Matrix "particle phase species concentration in molecules/cm^3"
		note PM3D_Matrix "Ncarbons = Rows; Noxygens = columns; nSizeBins = layers"
		//
		make/d/o/n=(nC,nO,nSizeBins) PMM3D_Matrix = 0 // g/cm^3
		note PMM3D_Matrix "particle phase species concentration in mass units"
		note PMM3D_Matrix "Ncarbons = Rows; Noxygens = columns; nSizeBins = layers"
		//
		make/d/o/n=(nC,nO) ParticleMass_Matrix = 0
		setscale/P x, nC, -1, "nCarbons" ParticleMass_Matrix
		setscale/P y, 0, 1, "nOxygens" ParticleMass_Matrix
		make/o/d/n=(nC,nO) TotalMolecules_Matrix = 0
		make/o/d/n=(nC,nO) TotalMass_Matrix = 0
		make/d/o/n=(nC,nO) TotalMolecules_Matrix_Log = 0
		make/d/o/n=(nC,nO) Particle_Fraction = 0
		make/d/o/n=(nC,nO) Particle_Fraction_Log = 0
		make/o/d/n=(nC,nO) ParticleMoleFraction_matrix = 0
		make/o/d/n=(nC,nO,nSizeBins) PMF3D_matrix = 0 // 3D version of ParticleMoleFraction_matrix; added 9/29/14
		make/o/d/n=(nC,nO) ParticleVolume_matrix
		make/d/o/n=(nC,nO) Particle_div_Total = 0
		//
		make/d/o/n=(nC,nO) CarbonAtoms = 0
		make/d/o/n=(nC,nO) OxygenAtoms = 0
		make/d/o/n=(nC,nO) HydrogenAtoms = 0
		//
		make/d/o/n=(nC,nO) Diffusivity_matrix = 0
		make/o/d/n=(nC,nO) MeanFreePath_Matrix = 0
		make/o/d/n=(nC,nO) Knudsen_Matrix = 0
		make/o/d/n=(nC,nO) Beta_Matrix = 0
		//
		make/d/o/n=(nC,nO) RxnStep_Matrix_plus = 0
		make/d/o/n=(nC,nO) RxnStep_Matrix_plus1 = 0
		make/d/o/n=(nC,nO) RxnStep_Matrix_plus2 = 0
		make/d/o/n=(nC,nO) RxnStep_Matrix_plus3 = 0
		make/d/o/n=(nC,nO) RxnStep_Matrix_plus4 = 0
		make/d/o/n=(nC,nO) RxnStep_Matrix_minus = 0
		make/o/d/n=(nC,nO) SumKinetics = 0 // added 8/4/18 for v7.3.8 to fix a bug in SOM_MC
		//
		make/o/d/n=(nC,nO) RxnStep_Matrix_Kinetic_On = 0
		make/o/d/n=(nC,nO) RxnStep_Matrix_Kinetic_Off = 0
		make/o/d/n=(nC,nO) RxnStep_Matrix_Kinetic_Pre1 = 0
		make/o/d/n=(nC,nO) RxnStep_Matrix_Kinetic_Pre2 = 0
		make/o/d/n=(nC,nO,nSteps) RxnStep_Matrix_Kinetic_Time = 0
		make/o/d/n=(nC*nO) RxnStep_array_minus
		make/d/o/n=(nC,nO) RxnStep_Matrix_frag = 0
		//
		make/d/o/n=(nC,nO) HetChem_Matrix_plus = 0
		make/d/o/n=(nC,nO) HetChem_Matrix_minus = 0
		// gas molecules
		make/d/o/n=(nC,nO,nsteps) GasMolecules_Time = 0
		note GasMolecules_Time "gas-phase species concentration in molecules/cm^3"
		note GasMolecules_Time "Ncarbons = Rows; Noxygens = columns; Layers = time [hrs]"
		setscale/P z, 0, (timestep/60/60), GasMolecules_Time
		setscale/P x, nC, -1, "nCarbons", GasMolecules_Time
		setscale/P y, 0, 1, "nOxygens", GasMolecules_time
		// gas mass
		make/d/o/n=(nC,nO,nsteps) GasMass_Time = 0
		note GasMass_Time "gas-phase species concentration in molecules/cm^3"
		note GasMass_Time "Ncarbons = Rows; Noxygens = columns; Layers = time [hrs]"
		setscale/P z, 0, (timestep/60/60), GasMass_Time
		setscale/P x, nC, -1, "nCarbons", GasMass_Time
		setscale/P y, 0, 1, "nOxygens", GasMass_Time
		// particle molecules
		make/d/o/n=(nC,nO,nsteps) ParticleMolecules_Time = 0
		note ParticleMolecules_Time "particle-phase species concentration in molecules/cm^3"
		note ParticleMolecules_Time "Ncarbons = Rows; Noxygens = columns; Layers = time [hrs]"
		setscale/P z, 0, (timestep/60/60), "hours", ParticleMolecules_Time
		setscale/P x, nC, -1, "nCarbons" ParticleMolecules_Time
		setscale/P y, 0, 1, "nOxygens" ParticleMolecules_Time
		// particle mass
		make/d/o/n=(nC,nO,nsteps) ParticleMass_Time = 0
		note ParticleMass_Time "particle-phase species concentration in molecules/cm^3"
		note ParticleMass_Time "Ncarbons = Rows; Noxygens = columns; Layers = time [hrs]"
		setscale/P z, 0, (timestep/60/60), "hours", ParticleMass_Time
		setscale/P x, nC, -1, "nCarbons" ParticleMass_Time
		setscale/P y, 0, 1, "nOxygens" ParticleMass_Time
		// total molecules
		make/d/o/n=(nC,nO,nsteps) TotalMolecules_Time = 0
		note TotalMolecules_Time "total species (gas+particle) concentration in molecules/cm^3"
		note TotalMolecules_Time "Ncarbons = Rows; Noxygens = columns; Layers = time [hrs]"
		setscale/P z, 0, (timestep/60/60), TotalMolecules_Time
		setscale/P x, nC, -1, "nCarbons" TotalMolecules_Time
		setscale/P y, 0, 1, "nOxygens" TotalMolecules_Time
		setscale/P d, 0, 0, "molecules/cm\S3", TotalMolecules_Time
		make/d/o/n=(nC,nO,nsteps) TotalMass_Time = 0
		make/o/d/n=(nC,nO,nsteps) ParticleMoleFraction_time = 0
		make/d/o/n=(nC,nO) ParticleFrac_Matrix = 0
		// vapor and particle wall loss
		make/d/o/n=(nC,nO,nsteps) WallMolecules_Time = 0
		note WallMolecules_Time "wall concentration in molecules/cm^3"
		note WallMolecules_Time "Ncarbons = Rows; Noxygens = columns; Layers = time [hrs]"
		setscale/P z, 0, (timestep/60/60), WallMolecules_Time
		setscale/P x, nC, -1, "nCarbons" WallMolecules_Time
		setscale/P y, 0, 1, "nOxygens" WallMolecules_Time
		setscale/P d, 0, 0, "molecules/cm\S3", WallMolecules_Time
		make/d/o/n=(nC,nO) WallMolecules_Matrix = 0 // absorbed wall molecules
		note WallMolecules_Matrix "wall species concentration in molecules/cm^3"
		note WallMolecules_Matrix "Ncarbons = Rows; Noxygens = columns"
		make/d/o/n=(nC,nO) WallParticleMolec_Matrix = 0 // deposited particles; added 9/29/14
		note WallParticleMolec_Matrix "Material associated with particle wall deposition, in molecules/cm^3"
		note WallParticleMolec_Matrix "Ncarbons = Rows; Noxygens = columns"
		// 3D version to deal with size-dependent wall losses; added 9/29/14
		make/d/o/n=(nC,nO,nsizebins) WPM3D_Matrix = 0 // molecules/cm^3
		note WPM3D_Matrix "particle phase species concentration in wall-deposited particles in molecules/cm^3"
		note WPM3D_Matrix "Ncarbons = Rows; Noxygens = columns; nSizeBins = layers"
		// Heterogeneous Chemistry
		make/o/d/n=(nC,nO,nsizebins) HetChem_Matrix_Minus // loss due to reaction
		make/o/d/n=(nC,nO,nsizebins) HetChem_Matrix_Plus // formation due to reaction
	
		make/d/o/n=(nsteps) TimeW = 0
		note TimeW "Reaction time [hrs]"
		// O2C for this compound class
		make/d/o/n=(nsteps) O2C_time = nan
		note O2C_time "Oxygen-to-Carbon ratio of SOA"
		setscale/P x, 0, (timestep/60/60), "hours", O2C_time
		setscale/P d, 0, 0, "O:C", O2C_time
		make/d/o/n=(nsteps) H2C_time = nan
		note H2C_time "hydrogen-to-carbon ratio of SOA"
		setscale/P x, 0, (timestep/60/60), "hours", H2C_time
		setscale/P d, 0, 0, "H:C", H2C_time
		// mean volatility for this compound class
		make/d/o/n=(nsteps) logCstar_mean_time = 0
		setscale/P x, 0, (timestep/60/60), "hours", logCstar_mean_time
		make/d/o/n=(nsteps) OCseed_time = 0
		// SOA associated with this VOC
		make/d/o/n=(nsteps) Coa_time = 0
		note Coa_time "SOA mass concentration [ug/m^3]"
		note Coa_time "Pfunc = " +num2str(POx1)+";"+num2str(POx2)+";"+num2str(POx3)+";"+num2str(POx4)
		note Coa_time "Pfrag = " + num2str(mfrag)
		note Coa_time "DLVP = " + num2str(DLVP)
		setscale/P x, 0, (timestep/60/60), "hours", Coa_time
		setscale/P d, 0, 0, "ug/m\S3\M", Coa_time
		// wall corrected value
		make/d/o/n=(nsteps) Coa_time_wallcorr = 0
		note Coa_time_wallcorr "Particle wall-loss corrected SOA mass concentration [ug/m^3]"
		note Coa_time_wallcorr "Coa_time_wallcorr = Coa_time + WallParticleMolec_Matrix"
		setscale/P x, 0, (timestep/60/60), "hours", Coa_time_wallcorr
		setscale/P d, 0, 0, "ug/m\S3\M", Coa_time_wallcorr
		// VOC reacted for this compound
		make/d/o/n=(nsteps) deltaHC_time = 0
		setscale/P x, 0, (timestep/60/60), "hours" deltaHC_time
		make/d/o/n=(nsteps) HC_ugm3_time = 0
		note HC_ugm3_time "Parent HC concentration [ugm3]"
		setscale/P x, 0, (timestep/60/60), "hours" HC_ugm3_time
		setscale/P d, 0, 0, "ug/m\S3\M" HC_ugm3_time
		make/d/o/n=(nsteps) deltaHCfrac_time = 0
		setscale/P x, 0, (timestep/60/60), "hours", deltaHCfrac_time
		// Yield for this compound
		make/d/o/n=(nsteps) Yield_time = 0
		note Yield_time "Aerosol mass yield"
		setscale/P x, 0, (timestep/60/60), "hours", Yield_time
		make/d/o/n=(nsteps) Yield2_time = 0
		make/d/o/n=(nsteps) Nc_ave_time = 0
		make/d/o/n=(nsteps) No_ave_time = 0
		// Lifetimes are VOC specific
		make/d/o/n=(nsteps) Lifetimes = 0
		note Lifetimes "Number of oxidation lifetimes"
		setscale/P x, 0, (timestep/60/60), "hours", Lifetimes
		// 
		make/o/d/n=(nsteps) Molecules_time = 0
		make/o/d/n=(nC+1) Ncarbons_Wave = nC-x
		make/o/d/n=(nO+1) Noxygens_Wave = x
	
	 	make/o/d/n=(nC,nO,nSizeBins) RxnStepKinetic3D_Cond, RxnStepKinetic3D_Evap
	 	make/o/d/n=(nC,nO,nSizeBins) RxnStepKinetic3D_Pre1, RxnStepKinetic3D_Pre2
		// mole fraction for this compound	
	 	make/o/d/n=(nC,nO,nSizeBins) MoleFraction3D
	 	
	 	// Make smaller waves to deal with fragmentation
	 	make/o/d/n=(SOMparams[0][i-1],nO) C_matrix_small, O_matrix_small, OC_matrix_small
	 	C_matrix_small = SOMparams[0][i-1]-x
	 	O_matrix_small = y
	 	OC_matrix_small = O_matrix_small/C_matrix_small
 	
 	endfor	
 	setdatafolder $df_root
End



//*************************************************************************************************************************
// Calculate gas-particle equilbirum based on "Pankow" theory (i.e. Raoult's Law)
// and "sequential partitioning model" of Cappa and Wilson (2011)
// Written by CDC on 10/24/15
// This is based on EqmCalc(), but extended to work for multiple VOCs
Function SPMcalc_MultipleVOC(SOMparams,[Coa_guess,SeedConc,SeedMW,units,Tvar])
	wave SOMparams // wave containing info on VOCs
	variable Coa_guess	// initial guess, in molecules/cm^3 or ug/m^3
	variable SeedConc // seed concentration in molecules/cm^3 or ug/m^3
	variable SeedMW // seed MW, if absorbing
	string units	// either "molecules" or "mass". "molecules" is default
	variable Tvar	// temperature, in K --> added 09/03/13

	if(ParamIsDefault(units))
		units = "molecules"
	endif
	if(ParamIsDefault(SeedConc))
		seedconc = 0
	endif
	if(ParamIsDefault(SeedMW))
		SeedMW = 400
	endif
	if(ParamIsDefault(Coa_Guess))
//		wavestats/Q ConcMatrix
		Coa_Guess = 1 //V_sum/4
	endif
	variable Tadjust
	if(ParamIsDefault(Tvar))	// added 09/03/13
		Tadjust = 0
	else
		Tadjust = 1
	endif
	
	// prep and deal with seed
	variable Coa_tol	// tolerance on result
	variable Seed
	if(stringmatch(units,"molecules"))
		Coa_tol = 1	// tolerance is 1 molecule
		Seed = SeedConc * 1e6 // convert molecules/cm^3 to molecules/m^3
	elseif(stringmatch(units,"mass"))
		Coa_tol = 1e-6	// tolerance is 1e-6 ug/m^3
		Seed = SeedConc * 1e6 * 1e6 * SeedMW / 6.022e23
	endif	
	Coa_Guess = seed
	Variable Coa_calc
	Variable Coa_dif
	variable counter=0	// added 08/31/13

	// deal with mulitple VOCs
	variable nVOCs = dimsize(SOMparams,1) // # of columns = # of precursor VOCs
	setdatafolder root:
	string cdf = getdatafolder(1)
	string df_home = "root:"
	string df_base = "root:VOCs" // Folder in which waves for individual VOCs are stored (in subfolders)
	string df_VOC // TBD; subfolders for storing individual waves from individual VOCs
	variable idex_VOCs // index for looping around different VOCs
	variable nC, nO
	string str_ConcMatrix, str_ConcSat25CforTD, str_MWwave, str_deltaHvap_matrix
	
	for(idex_VOCs=0;idex_VOCs<nVOCs;idex_VOCs+=1)
		df_VOC = df_base + ":VOC" + num2istr(idex_VOCs+1)
		setdatafolder $df_VOC
		// access waves
		// The only difference between this and the eqm model is this line, I think.
		wave ConcMatrix = GasMolecules_Matrix  // in molecules/cm^3 
		wave ConcSat25CForTD = SatConc_matrix // saturation vapor pressure in molecules/m^3
		wave MWwave = MW_matrix // molecular weight matrix in g/mol
		wave deltaHvap_matrix	// matrix of deltaHvap values, in kJ/mol --> added 09/03/13
		nC = dimsize(ConcMatrix,0)
		nO = dimsize(ConcMatrix,1)
		// create some new waves in the VOC folders
		make/o/d/n=(nC,nO) ConcWave=0, FracAero=0, Coa=0
		make/o/d/n=(nC,nO) SatConc_temp=0

		if(idex_VOCs==0) 
			duplicate/o/d Coa root:Coa // make a "total" wave in the root folder
			wave CoaSum = root:Coa
			CoaSum = 0
		endif
		
		if(stringmatch(units,"molecules"))
			ConcWave = ConcMatrix * 1e6		// convert molecules/cm^3 to molecules/m^3
			ConcWave = numtype(ConcWave)==2 ? 0 : ConcWave
			SatConc_temp = ConcSat25CforTD	// no units change
		elseif(stringmatch(units,"mass"))
			ConcWave = ConcMatrix * 1e6 * 1e6 * MWwave / 6.022e23	// convert molecules/cm^3 to ug/m^3
			ConcWave = numtype(ConcWave)==2 ? 0 : ConcWave
			SatConc_temp = ConcSat25CforTD * 1e6 * MWwave / 6.022e23	// convert molecules/m^3 to ug/m^3
		endif	
	
		if(Tadjust==1)
			// adjust vapor pressure based on Tvar and deltaHvap_matrix --> added 09/03/13
			// assumes ConcSat25CForTD at 25C
			variable Tref = 298
			if(stringmatch(units,"molecules"))
				SatConc_temp[][] = SatConc_temp[p][q]*exp(-1*(deltaHvap_matrix[p][q]*1000/8.314)*(1/Tvar - 1/Tref))
			elseif(stringmatch(units,"mass"))
				SatConc_temp[][] = SatConc_temp[p][q]*(Tref/Tvar)*exp(-1*(deltaHvap_matrix[p][q]*1000/8.314)*(1/Tvar - 1/Tref))
			endif
		endif

		WaveStats/Q ConcWave
		Coa_guess += V_sum/4 // keep adding for each VOC

		setdatafolder $df_home // back to start folder
	endfor // end multiple VOCs (for now)
	
	do // loop until tolerance is met
		Coa_calc = seed // reset to seed value every time
		CoaSum = 0
		for(idex_VOCs=0;idex_VOCs<nVOCs;idex_VOCs+=1) // loop over VOCs
			// go to a specific VOC folder
			df_VOC = df_base + ":VOC" + num2istr(idex_VOCs+1)
			setdatafolder $df_VOC
			
			wave ConcWave
			wave FracAero
			wave SatConc_temp
			wave Coa
			FracAero = 1/(1+(SatConc_temp)/Coa_guess) // particle-phase fraction for a given VOC (partitioning equation)
			Coa = FracAero * ConcWave // Coa for a given VOC
			wavestats/Q Coa // get total amount
			Coa_calc += V_Sum // accumulate for each VOC
			CoaSum += Coa // accumulate for each VOC
		endfor
		Coa_dif = abs(Coa_calc - Coa_guess) // compare actual with guess value
		if(Coa_dif < Coa_tol) // compare to tolerance
			break
		endif
		Coa_guess = Coa_calc // reset "guess" for next iteration
		counter += 1 // update counter
		if(counter==5)	// reset tolerance based on current guess; this speeds things up since the desired tolerance depends on the amount of OA there is.
			Coa_tol = max(Coa_tol,0.01*(Coa_calc-Seed))
		endif
	while(1)

	CoaSum = 0 // reset 
	for(idex_VOCs=0;idex_VOCs<nVOCs;idex_VOCs+=1) // loop over VOCs
		// go to a specific VOC folder
		df_VOC = df_base + ":VOC" + num2istr(idex_VOCs+1)
		setdatafolder $df_VOC
		wave Coa
		wave MWwave = MW_matrix
		if(stringmatch(units,"molecules"))
			Coa = Coa * 1e-6 // convert molecules/m^3 to molecules/cm^3...does not include seed mass
			CoaSum += Coa
		elseif(stringmatch(units,"mass"))
			Coa = Coa * 6.022e23 / (1e6 * 1e6 * MWwave)	// convert ug/m^3 to molecules/cm^3
			CoaSum += Coa
		endif
	endfor
	setdatafolder $cdf
End	

//****************************************************************************************************************

Function SetCalNex(VWL,NOx)
	variable VWL
	string NOx
	
	setdatafolder root:
	wave SOMparams = root:SOMparams
	NVAR gasWLR

	setdatafolder root:BestFitParams
	wave mfrag
	wave deltaLVP
	wave pOx1
	wave pOx2
	wave pOx3
	wave pOx4
	
	make/o/d/n=6 idex_wave
	variable i
	
	if(VWL==0.001 && stringmatch(NOx,"low"))
		idex_wave = {0,1,3,7,5,9}
	elseif(VWL==0.001 && stringmatch(NOx,"high"))
		idex_wave = {12,13,15,19,17,20}
	elseif(VWL==0.0025 && stringmatch(NOx,"low"))
		idex_wave = {23,24,26,30,28,32}
	elseif(VWL==0.0025 && stringmatch(NOx,"high"))
		idex_wave = {35,36,38,42,40,44}
	elseif(VWL==0 && stringmatch(NOx,"low"))
		idex_wave = {47,48,50,54,52,56}
	elseif(VWL==0 && stringmatch(NOx,"high"))
		idex_wave = {59,60,62,66,64,67}
	else
		abort "bad input"
	endif
	
	for(i=0;i<6;i+=1)
		SOMparams[1][i] = mfrag[idex_wave[i]]
		SOMparams[2][i] = deltaLVP[idex_wave[i]]
		SOMparams[3][i] = pOx1[idex_wave[i]]
		SOMparams[4][i] = pOx2[idex_wave[i]]
		SOMparams[5][i] = pOx3[idex_wave[i]]
		SOMparams[6][i] = pOx4[idex_wave[i]]
	endfor
	
	killwaves/Z idex_wave
	setdatafolder root:
End


Function GraphSize()

	wave dndlogdp_time
	wave diameter_time
	variable npnts = dimsize(dndlogdp_time,0)
	display dndlogdp_time[0][] vs diameter_time[0][]
	modifygraph log(bottom)=1
	variable i
	for(i=1;i<npnts;i+=100)
		appendtograph dndlogdp_time[i][] vs diameter_time[i][]
	endfor
end

//*****************************************************************************
Function SOM_SetupMulticomponent()

	setdatafolder root:
	NVAR nCompoundClasses // # of compound classes
	string TheCompoundClasses = "isoprene;aPinene;sesq;benzene;toluene;xylenes;naphthalene;n-alkane;b-alkane;c-alkane"
	string C1, C2, C3, C4, C5, C6, C7, C8, C9, C10
	string selectedCompounds
	if(nCompoundClasses==1)
		Prompt C1, "Choose C#1", popup, TheCompoundClasses
		DoPrompt "Choose Compounds", C1
		selectedCompounds = C1
	elseif(nCompoundClasses==2)
		Prompt C1, "Choose C#1", popup, TheCompoundClasses
		Prompt C2, "Choose C#2", popup, TheCompoundClasses
		DoPrompt "Choose Compounds", C1, C2
		selectedCompounds = C1 + ";" + C2
	elseif(nCompoundClasses==3)
		Prompt C1, "Choose C#1", popup, TheCompoundClasses
		Prompt C2, "Choose C#2", popup, TheCompoundClasses
		Prompt C3, "Choose C#3", popup, TheCompoundClasses
		DoPrompt "Choose Compounds", C1, C2, C3
		selectedCompounds = C1 + ";" + C2 + ";" + C3
	elseif(nCompoundClasses==4)
		Prompt C1, "Choose C#1", popup, TheCompoundClasses
		Prompt C2, "Choose C#2", popup, TheCompoundClasses
		Prompt C3, "Choose C#3", popup, TheCompoundClasses
		Prompt C4, "Choose C#4", popup, TheCompoundClasses
		DoPrompt "Choose Compounds", C1, C2, C3, C4
		selectedCompounds = C1 + ";" + C2 + ";" + C3 + ";" + C4
	elseif(nCompoundClasses==5)
		Prompt C1, "Choose C#1", popup, TheCompoundClasses
		Prompt C2, "Choose C#2", popup, TheCompoundClasses
		Prompt C3, "Choose C#3", popup, TheCompoundClasses
		Prompt C4, "Choose C#4", popup, TheCompoundClasses
		Prompt C5, "Choose C#5", popup, TheCompoundClasses
		DoPrompt "Choose Compounds", C1, C2, C3, C4, C5
		selectedCompounds = C1 + ";" + C2 + ";" + C3 + ";" + C4 + ";" + C5
	elseif(nCompoundClasses==6)
		Prompt C1, "Choose C#1", popup, TheCompoundClasses
		Prompt C2, "Choose C#2", popup, TheCompoundClasses
		Prompt C3, "Choose C#3", popup, TheCompoundClasses
		Prompt C4, "Choose C#4", popup, TheCompoundClasses
		Prompt C5, "Choose C#5", popup, TheCompoundClasses
		Prompt C6, "Choose C#6", popup, TheCompoundClasses
		DoPrompt "Choose Compounds", C1, C2, C3, C4, C5, C6
		selectedCompounds = C1 + ";" + C2 + ";" + C3 + ";" + C4 + ";" + C5 + ";" + C6
	elseif(nCompoundClasses==7)
		Prompt C1, "Choose C#1", popup, TheCompoundClasses
		Prompt C2, "Choose C#2", popup, TheCompoundClasses
		Prompt C3, "Choose C#3", popup, TheCompoundClasses
		Prompt C4, "Choose C#4", popup, TheCompoundClasses
		Prompt C5, "Choose C#5", popup, TheCompoundClasses
		Prompt C6, "Choose C#6", popup, TheCompoundClasses
		Prompt C7, "Choose C#7", popup, TheCompoundClasses
		DoPrompt "Choose Compounds", C1, C2, C3, C4, C5, C6, C7
		selectedCompounds = C1 + ";" + C2 + ";" + C3 + ";" + C4 + ";" + C5 + ";" + C6 + ";" + C7
	elseif(nCompoundClasses==8)
		Prompt C1, "Choose C#1", popup, TheCompoundClasses
		Prompt C2, "Choose C#2", popup, TheCompoundClasses
		Prompt C3, "Choose C#3", popup, TheCompoundClasses
		Prompt C4, "Choose C#4", popup, TheCompoundClasses
		Prompt C5, "Choose C#5", popup, TheCompoundClasses
		Prompt C6, "Choose C#6", popup, TheCompoundClasses
		Prompt C7, "Choose C#7", popup, TheCompoundClasses
		Prompt C8, "Choose C#8", popup, TheCompoundClasses
		DoPrompt "Choose Compounds", C1, C2, C3, C4, C5, C6, C7, C8
		selectedCompounds = C1 + ";" + C2 + ";" + C3 + ";" + C4 + ";" + C5 + ";" + C6 + ";" + C7 + ";" + C8
	elseif(nCompoundClasses==9)
		Prompt C1, "Choose C#1", popup, TheCompoundClasses
		Prompt C2, "Choose C#2", popup, TheCompoundClasses
		Prompt C3, "Choose C#3", popup, TheCompoundClasses
		Prompt C4, "Choose C#4", popup, TheCompoundClasses
		Prompt C5, "Choose C#5", popup, TheCompoundClasses
		Prompt C6, "Choose C#6", popup, TheCompoundClasses
		Prompt C7, "Choose C#7", popup, TheCompoundClasses
		Prompt C8, "Choose C#8", popup, TheCompoundClasses
		Prompt C9, "Choose C#9", popup, TheCompoundClasses
		DoPrompt "Choose Compounds", C1, C2, C3, C4, C5, C6, C7, C8, C9
		selectedCompounds = C1 + ";" + C2 + ";" + C3 + ";" + C4 + ";" + C5 + ";" + C6 + ";" + C7 + ";" + C8 + ";" + C9
	elseif(nCompoundClasses==10)
		Prompt C1, "Choose C#1", popup, TheCompoundClasses
		Prompt C2, "Choose C#2", popup, TheCompoundClasses
		Prompt C3, "Choose C#3", popup, TheCompoundClasses
		Prompt C4, "Choose C#4", popup, TheCompoundClasses
		Prompt C5, "Choose C#5", popup, TheCompoundClasses
		Prompt C6, "Choose C#6", popup, TheCompoundClasses
		Prompt C7, "Choose C#7", popup, TheCompoundClasses
		Prompt C8, "Choose C#8", popup, TheCompoundClasses
		Prompt C9, "Choose C#9", popup, TheCompoundClasses
		Prompt C10, "Choose C#10", popup, TheCompoundClasses
		DoPrompt "Choose Compounds", C1, C2, C3, C4, C5, C6, C7, C8, C9, C10
		selectedCompounds = C1 + ";" + C2 + ";" + C3 + ";" + C4 + ";" + C5 + ";" + C6 + ";" + C7 + ";" + C8 + ";" + C9 + ";" + C10
	else
		abort "Bad Choice in nCompoundClasses"		
	endif
	
	if(V_flag==1)
		abort "Cancelled update. No parameters changed."
	endif
	
	print "You have just set the following compounds:"
	print selectedCompounds

	make/o/d/n=(8,nCompoundClasses) SOMparams=0 // the SOM parameters for all compounds being condsidered
	note/K SOMparams, "Compounds are "

//	print selectedCompounds
	
	variable i, j
	string myCompound
	SVAR LowOrHighNOx
	string lowOrHighNOx_temp
	SVAR VWLcondition
	NVAR maxCforAlkanes
	string df_fits = "root:BestFitParams"
	wave kwall = $(df_fits+":kwall")
	wave/T VOC = $(df_fits+":VOC")
	wave/T NOx = $(df_fits + ":NOx")
	wave mfrag_wave = $(df_fits + ":mfrag")
	wave deltaLVP_wave = $(df_fits + ":deltaLVP")
	wave pOx1_wave = $(df_fits + ":pOx1")
	wave pOx2_wave = $(df_fits + ":pOx2")
	wave pOx3_wave = $(df_fits + ":pOx3")
	wave pOx4_wave = $(df_fits + ":pOx4")
	wave nCarbons_wave = $(df_fits + ":Ncarbons")
	string myVOC, myNOx, mykwall
	variable maxNumCarbons = 0
	variable flag = 0
	
	// Select appropriate fit conditions to use
	for(i=0;i<nCompoundClasses;i+=1)
		myCompound = stringfromlist(i,selectedCompounds,";")
		if(stringmatch(myCompound,"xylenes")==1)
			myCompound = "mxylene"
		elseif(stringmatch(myCompound,"sesq")==1)
			myCompound = "aPinene"
			flag = 1
		elseif(stringmatch(myCompound,"n-alkane")==1)
			myCompound = "dodecane"
			flag = 2
		elseif(stringmatch(myCompound,"b-alkane")==1)
			myCompound = "methylundecane" 
			flag = 2
		elseif(stringmatch(myCompound,"c-alkane")==1)
			myCompound = "hexylcyclohexane"
			flag = 2
		elseif(stringmatch(myCompound,"toluene")==1)
			myCompound = "toluene2013"
		endif
		// loop through all fits and select the right one for the compound/nox/kwall condition of interest
		for(j=0;j<100;j+=1)
			lowOrHighNOx_temp = lowOrHighNOx
			myVOC = VOC[j]
			myNOx = NOx[j]
			if(kwall[j]==0)
				mykwall = "zero"
			elseif(kwall[j]==1e-4)
				mykwall = "low"
			elseif(kwall[j]==2.5e-4)
				mykwall = "high"
			endif
			if(stringmatch(myCompound,"mxylene")==1 && stringmatch(lowOrHighNOx,"high"))
				lowOrHighNOx_temp = "high2"
			elseif(stringmatch(myCompound,"mxylene")==1 && stringmatch(lowOrHighNOx,"low"))
				lowOrHighNOx_temp = "low2"
			elseif(stringmatch(myCompound,"toluene2013")==1 && stringmatch(lowOrHighNOx,"high"))
				lowOrHighNOx_temp = "high2"
			elseif(stringmatch(myCompound,"toluene2013")==1 && stringmatch(lowOrHighNOx,"low"))
				lowOrHighNOx_temp = "low3"
			endif
			if(stringmatch(myVOC,myCompound)==1 && stringmatch(myNOx,LowOrHighNOx_temp)==1 && stringmatch(mykwall,VWLcondition)==1)
				break
			endif
		endfor
		
		// set SOM parameters for each compound
		if(flag==0)
			SOMparams[0][i] = Ncarbons_wave[j]
			if(Ncarbons_wave[j] > maxNumCarbons)
				maxNumCarbons = Ncarbons_wave[j]
			endif
		elseif(flag==1)
			SOMparams[0][i] = 15 // sesquiterepenes are special
			if(15 > maxNumCarbons)
				maxNumCarbons = 15
			endif
		elseif(flag==2)
			SOMparams[0][i] = maxCforAlkanes
			if(maxCforAlkanes > maxNumCarbons)
				maxNumCarbons = maxCforAlkanes
			endif
		endif
		SOMparams[1][i] = mfrag_wave[j]
		SOMparams[2][i] = deltaLVP_wave[j]
		SOMparams[3][i] =pOx1_wave[j]
		SOMparams[4][i] =pOx2_wave[j]
		SOMparams[5][i] =pOx3_wave[j]
		SOMparams[6][i] =pOx4_wave[j]
		if(flag==2)
			SOMparams[7][i] = 0 // use defaul alkane parameterization
		else
			SOMparams[7][i] = Set_kOH(myCompound)
		endif
		flag = 0 // reset flag to zero
		note SOMparams, myCompound+";"
		SetDimLabel 1, i, $(stringfromlist(i,selectedCompounds,";")), SOMparams
endfor

	// set initial concentrations
	// for alkanes, assume a distribution of IVOCs and specify total concentration
	string InitialConcentrations
	Prompt InitialConcentrations, "Initial Concentrations in ppb as X;Y;Z..."
	DoPrompt "Initial Concentrations as X;Y;Z...", InitialConcentrations
	
	string df_IVOC = "root:IVOCs"
	make/o/d/n=(maxNumCarbons,nCompoundClasses) InitialConcMatrix = 0
	DoWindow/F SOM_Setup
//	wave InitialConcMatrix // initial concentrations in ppm
		// for non-alkanes, the concentration goes in the first row, and the SOM code will parse this
		// for alkanes, the concentrations are spread throughout the rows and use a weighted distribution

	for(i=0;i<nCompoundClasses;i+=1)	
		myCompound = stringfromlist(i,selectedCompounds,";")
		if(stringmatch(myCompound,"xylenes")==1)
			myCompound = "mxylene"
		elseif(stringmatch(myCompound,"sesq")==1)
			myCompound = "aPinene"
			flag = 1
		elseif(stringmatch(myCompound,"n-alkane")==1)
			myCompound = "dodecane"
			flag = 2
		elseif(stringmatch(myCompound,"b-alkane")==1)
			myCompound = "methylundecane" 
			flag = 2
		elseif(stringmatch(myCompound,"c-alkane")==1)
			myCompound = "hexylcyclohexane"
			flag = 2
		elseif(stringmatch(myCompound,"toluene")==1)
			myCompound = "toluene2013"
		endif

		// now set concentrations
		if(stringmatch(myCompound,"dodecane")==1)
			wave weightingfactors = $(df_IVOC+":n_alkanes_norm")
			InitialConcMatrix[][i] = weightingfactors[p]*str2num(stringfromlist(i,InitialConcentrations,";"))/1000
		elseif(stringmatch(myCompound,"methylundecane")==1)
			wave weightingfactors = $(df_IVOC+":b_alkanes_norm")
			InitialConcMatrix[][i] = weightingfactors[p]*str2num(stringfromlist(i,InitialConcentrations,";"))/1000
		elseif(stringmatch(myCompound,"hexylcyclohexane")==1)
			wave weightingfactors = $(df_IVOC+":c_alkanes_norm")
			InitialConcMatrix[][i] = weightingfactors[p]*str2num(stringfromlist(i,InitialConcentrations,";"))/1000
		else // all other compounds
// CDC 07/05/17	InitialConcMatrix[0][i] = str2num(stringfromlist(i,InitialConcentrations,";"))/1000 // ppb --> ppm
			InitialConcMatrix[maxNumCarbons - SOMparams[0][i]][i] = str2num(stringfromlist(i,InitialConcentrations,";"))/1000 // ppb --> ppm // changed on 07/05/17 by CDC
		endif
		SetDimLabel 1, i, $(stringfromlist(i,selectedCompounds,";")), InitialConcMatrix
	endfor
	
End

//*****************************************************************************
Function SOM_SetupSingleComponent()

	setdatafolder root:
	
	variable i, j
	string myCompound
	SVAR LowOrHighNOx
	string lowOrHighNOx_temp
	SVAR VWLcondition
	SVAR TheCompound
	
	NVAR maxCforAlkanes
	NVAR isopreneisspecial
	
	string df_fits = "root:BestFitParams"
	wave kwall = $(df_fits+":kwall")
	wave/T VOC = $(df_fits+":VOC")
	wave/T NOx = $(df_fits + ":NOx")
	wave mfrag_wave = $(df_fits + ":mfrag")
	wave deltaLVP_wave = $(df_fits + ":deltaLVP")
	wave pOx1_wave = $(df_fits + ":pOx1")
	wave pOx2_wave = $(df_fits + ":pOx2")
	wave pOx3_wave = $(df_fits + ":pOx3")
	wave pOx4_wave = $(df_fits + ":pOx4")
	wave nCarbons_wave = $(df_fits + ":Ncarbons")
	string myVOC, myNOx, mykwall
	variable maxNumCarbons = 0
	variable flag = 0
	
	// Select appropriate fit conditions to use
	for(i=0;i<1;i+=1)
		myCompound = TheCompound
		if(stringmatch(myCompound,"xylenes")==1)
			myCompound = "mxylene"
		elseif(stringmatch(myCompound,"sesq")==1)
			myCompound = "aPinene"
			flag = 1
		elseif(stringmatch(myCompound,"n-alkane")==1)
			myCompound = "dodecane"
			flag = 2
		elseif(stringmatch(myCompound,"b-alkane")==1)
			myCompound = "methylundecane" 
			flag = 2
		elseif(stringmatch(myCompound,"c-alkane")==1)
			myCompound = "hexylcyclohexane"
			flag = 2
		elseif(stringmatch(myCompound,"toluene")==1)
			myCompound = "toluene2013"
		endif
		// loop through all fits and select the right one for the compound/nox/kwall condition of interest
		for(j=0;j<150;j+=1)
			lowOrHighNOx_temp = lowOrHighNOx
			myVOC = VOC[j]
			myNOx = NOx[j]
			if(kwall[j]==0)
				mykwall = "zero"
			elseif(kwall[j]==1e-4)
				mykwall = "low"
			elseif(kwall[j]==2.5e-4)
				mykwall = "high"
			endif
			if(stringmatch(myCompound,"mxylene")==1 && stringmatch(lowOrHighNOx,"high"))
				lowOrHighNOx_temp = "high2"
			elseif(stringmatch(myCompound,"mxylene")==1 && stringmatch(lowOrHighNOx,"low"))
				lowOrHighNOx_temp = "low2"
			elseif(stringmatch(myCompound,"toluene2013")==1 && stringmatch(lowOrHighNOx,"high"))
				lowOrHighNOx_temp = "high2"
			elseif(stringmatch(myCompound,"toluene2013")==1 && stringmatch(lowOrHighNOx,"low"))
				lowOrHighNOx_temp = "low3"
			endif
			if(stringmatch(myVOC,myCompound)==1 && stringmatch(myNOx,LowOrHighNOx_temp)==1 && stringmatch(mykwall,VWLcondition)==1)
				break 
			endif
			if(j==102)
				abort  "The combination of VOC, NOx and VWL does not exist"
			endif
		endfor
	endfor
	
	NVAR delta_logCstar_perO
	NVAR fragslope
	NVAR probOx1, probOx2, ProbOx3, ProbOx4
	NVAR krxn_parent
	NVAR gasWLR
	NVAR Ncarbons
	print flag
			// set SOM parameters for each compound
		if(flag==0)
			Ncarbons = Ncarbons_wave[j]
		elseif(flag==1)
			Ncarbons = 15 // sesquiterepenes are special
		elseif(flag==2)
			// do nothing, v7.3.7
			//Ncarbons = maxCforAlkanes
		endif
		fragslope = mfrag_wave[j]
		delta_logCstar_perO = deltaLVP_wave[j]
		probOx1 =pOx1_wave[j]
		probOx2 =pOx2_wave[j]
		probOx3 =pOx3_wave[j]
		probOx4 =pOx4_wave[j]
		if(flag==2)
			krxn_parent = 0 // use defaul alkane parameterization
		else
			krxn_parent = Set_kOH(myCompound)
		endif	
		if(stringmatch(myCompound,"isoprene_alt")) // added for v7.3.5
			isopreneisspecial = 1
		else
			isopreneisspecial = 0
		endif
End

//************************************************************************
Function Graph_MC_soa()

	wave SOMparams
	variable nCompounds = dimsize(SOMparams,1)
	string df_base = "root:VOCs"
	string df_VOC
	string current
	variable i
	
	DoWindow/F SOA
	if(V_flag==1)
		KillWindow SOA
	endif
	display/N=SOA as "SOA"
	for(i=0;i<nCompounds;i+=1)
		df_VOC = df_base+":VOC"+num2istr(i+1)
		wave Coa = $(df_VOC+":Coa_time")
		appendtograph Coa
		if(i==0)
			current = "Coa_time"
			textbox/C/A=LT/N=text0 "\\Z08"

		else
			current = "Coa_time#"+num2istr(i)
		endif
		if(stringmatch(GetDimLabel(SOMparams,1,i),"isoprene"))
			ModifyGraph rgb($current)=(65280,32768,58880) // pink
			AppendText/N=text0 "\\s("+current+") isoprene"
		elseif(stringmatch(GetDimLabel(SOMparams,1,i),"aPinene"))
			ModifyGraph rgb($current)=(65280,54528,32768) // orange
			AppendText/N=text0 "\\s("+current+") aPinene"
		elseif(stringmatch(GetDimLabel(SOMparams,1,i),"sesq"))
			ModifyGraph rgb($current)=(51456,44032,58880) // purple
			AppendText/N=text0 "\\s("+current+") sequiterpene"
		elseif(stringmatch(GetDimLabel(SOMparams,1,i),"benzene"))
			ModifyGraph rgb($current)=(52224,17408,0) // auburn
			AppendText/N=text0 "\\s("+current+") benzene"
		elseif(stringmatch(GetDimLabel(SOMparams,1,i),"toluene"))
			ModifyGraph rgb($current)=(52224,34816,0) // brown
			AppendText/N=text0 "\\s("+current+") toluene"
		elseif(stringmatch(GetDimLabel(SOMparams,1,i),"xylenes"))
			ModifyGraph rgb($current)=(16384,48896,65280) // light blue
			AppendText/N=text0 "\\s("+current+") xylenes"
		elseif(stringmatch(GetDimLabel(SOMparams,1,i),"n-alkane"))
			ModifyGraph rgb($current)=(0,0,65280) // blue
			AppendText/N=text0 "\\s("+current+") n-alkanes"
		elseif(stringmatch(GetDimLabel(SOMparams,1,i),"b-alkane"))
			ModifyGraph rgb($current)=(0,52224,0) // green
			AppendText/N=text0 "\\s("+current+") b-alkanes"
		elseif(stringmatch(GetDimLabel(SOMparams,1,i),"c-alkane"))
			ModifyGraph rgb($current)=(30464,30464,30464) // gray
			AppendText/N=text0 "\\s("+current+") c-alkanes"
		elseif(stringmatch(GetDimLabel(SOMparams,1,i),"naphthalene"))
			ModifyGraph rgb($current)=(39168,0,0) // deep red
			AppendText/N=text0 "\\s("+current+") naphthalene"
		endif
 	endfor
 	
 	ModifyGraph mode=7,hbFill=2,toMode=3
	ModifyGraph mode=7,hbFill=2,toMode=3
	ModifyGraph fSize=14,axThick=1.5,standoff=0
	Label left "[SOA] (ug/m3)\\e"
	Label bottom "Time (hrs)\\e"

End

//************************************************************************
Function Graph_MC_GandP()

	wave SOMparams
	variable nCompounds = dimsize(SOMparams,1)
	string df_base = "root:VOCs"
	string df_VOC
	string current
	variable i
	
	DoWindow/F GasAndParticle
	if(V_flag==1)
		KillWindow GasAndParticle
	endif
	display/N=GasAndParticle as "GasAndParticle"
	for(i=0;i<nCompounds;i+=1)
		df_VOC = df_base+":VOC"+num2istr(i+1)
		wave Coa = $(df_VOC+":Coa_time")
		appendtograph Coa
		if(i==0)
			current = "Coa_time"
			textbox/C/A=LT/N=text0 "\\Z08"
		else
			current = "Coa_time#"+num2istr(i)
		endif
		if(stringmatch(GetDimLabel(SOMparams,1,i),"isoprene"))
			ModifyGraph rgb($current)=(65280,32768,58880) // pink
			AppendText/N=text0 "\\s("+current+") isoprene"
		elseif(stringmatch(GetDimLabel(SOMparams,1,i),"aPinene"))
			ModifyGraph rgb($current)=(65280,54528,32768) // orange
			AppendText/N=text0 "\\s("+current+") aPinene"
		elseif(stringmatch(GetDimLabel(SOMparams,1,i),"sesq"))
			ModifyGraph rgb($current)=(51456,44032,58880) // purple
			AppendText/N=text0 "\\s("+current+") sequiterpene"
		elseif(stringmatch(GetDimLabel(SOMparams,1,i),"benzene"))
			ModifyGraph rgb($current)=(52224,17408,0) // auburn
			AppendText/N=text0 "\\s("+current+") benzene"
		elseif(stringmatch(GetDimLabel(SOMparams,1,i),"toluene"))
			ModifyGraph rgb($current)=(52224,34816,0) // brown
			AppendText/N=text0 "\\s("+current+") toluene"
		elseif(stringmatch(GetDimLabel(SOMparams,1,i),"xylenes"))
			ModifyGraph rgb($current)=(16384,48896,65280) // light blue
			AppendText/N=text0 "\\s("+current+") xylenes"
		elseif(stringmatch(GetDimLabel(SOMparams,1,i),"n-alkane"))
			ModifyGraph rgb($current)=(0,0,65280) // blue
			AppendText/N=text0 "\\s("+current+") n-alkanes"
		elseif(stringmatch(GetDimLabel(SOMparams,1,i),"b-alkane"))
			ModifyGraph rgb($current)=(0,52224,0) // green
			AppendText/N=text0 "\\s("+current+") b-alkanes"
		elseif(stringmatch(GetDimLabel(SOMparams,1,i),"c-alkane"))
			ModifyGraph rgb($current)=(30464,30464,30464) // gray
			AppendText/N=text0 "\\s("+current+") c-alkanes"
		elseif(stringmatch(GetDimLabel(SOMparams,1,i),"naphthalene"))
			ModifyGraph rgb($current)=(39168,0,0) // deep red
			AppendText/N=text0 "\\s("+current+") naphthalene"
		endif
 	endfor
 	
 	for(i=0;i<nCompounds;i+=1)
		df_VOC = df_base+":VOC"+num2istr(i+1)
		wave HC = $(df_VOC+":HC_ugm3_time")
		appendtograph/l=L_Gas /b=B_Gas HC
		if(i==0)
			current = "HC_ugm3_time"
		else
			current = "HC_ugm3_time#"+num2istr(i)
		endif
		if(stringmatch(GetDimLabel(SOMparams,1,i),"isoprene"))
			ModifyGraph rgb($current)=(65280,32768,58880) // pink
		elseif(stringmatch(GetDimLabel(SOMparams,1,i),"aPinene"))
			ModifyGraph rgb($current)=(65280,54528,32768) // orange
		elseif(stringmatch(GetDimLabel(SOMparams,1,i),"sesq"))
			ModifyGraph rgb($current)=(51456,44032,58880) // purple
		elseif(stringmatch(GetDimLabel(SOMparams,1,i),"benzene"))
			ModifyGraph rgb($current)=(52224,17408,0) // auburn
		elseif(stringmatch(GetDimLabel(SOMparams,1,i),"toluene"))
			ModifyGraph rgb($current)=(52224,34816,0) // brown
		elseif(stringmatch(GetDimLabel(SOMparams,1,i),"xylenes"))
			ModifyGraph rgb($current)=(16384,48896,65280) // light blue
		elseif(stringmatch(GetDimLabel(SOMparams,1,i),"n-alkane"))
			ModifyGraph rgb($current)=(0,0,65280) // blue
		elseif(stringmatch(GetDimLabel(SOMparams,1,i),"b-alkane"))
			ModifyGraph rgb($current)=(0,52224,0) // green
		elseif(stringmatch(GetDimLabel(SOMparams,1,i),"c-alkane"))
			ModifyGraph rgb($current)=(30464,30464,30464) // gray
		elseif(stringmatch(GetDimLabel(SOMparams,1,i),"naphthalene"))
			ModifyGraph rgb($current)=(39168,0,0) // deep red
		endif
 	endfor
 	
 	ModifyGraph mode=7,hbFill=2,toMode=3
	ModifyGraph mode=7,hbFill=2,toMode=3
	ModifyGraph fSize=14,axThick=1.5,standoff=0
	ModifyGraph freePos(L_Gas)={0,bottom},freePos(B_Gas)={0,L_Gas}
	ModifyGraph axisEnab(left)={0,0.46},axisEnab(L_Gas)={0.54,1}
	SetAxis L_Gas 0,*
	Label left "[SOA] (ug/m3)\\e"
	Label bottom "Time (hrs)\\e"
	Label L_Gas "[VOC] (ug/m3)\\e"
	ModifyGraph axOffset(L_Gas)=-5
	ModifyGraph axOffset=0,lblPos(L_Gas)=60
	ModifyGraph width=290,height=400

End

//************************************************************************
Function Graph_MC_Yield()

	setdatafolder root:
	wave SOMparams
	wave Coa_time
	variable nCompounds = dimsize(SOMparams,1)
	string df_base = "root:VOCs"
	string df_VOC
	string current
	variable i
	
	DoWindow/F Yield
	if(V_flag==1)
		KillWindow Yield
	endif
	display/N=Yield as "Yield"
	for(i=0;i<nCompounds;i+=1)
		df_VOC = df_base+":VOC"+num2istr(i+1)
		wave Yield = $(df_VOC+":Yield_time")
		appendtograph Yield vs Coa_time
		if(i==0)
			current = "Yield_time"
			textbox/C/A=LT/N=text0 "\\Z08"
		else
			current = "Yield_time#"+num2istr(i)
		endif
		if(stringmatch(GetDimLabel(SOMparams,1,i),"isoprene"))
			ModifyGraph rgb($current)=(65280,32768,58880) // pink
			AppendText/N=text0 "\\s("+current+") isoprene"
		elseif(stringmatch(GetDimLabel(SOMparams,1,i),"aPinene"))
			ModifyGraph rgb($current)=(65280,54528,32768) // orange
			AppendText/N=text0 "\\s("+current+") aPinene"
		elseif(stringmatch(GetDimLabel(SOMparams,1,i),"sesq"))
			ModifyGraph rgb($current)=(51456,44032,58880) // purple
			AppendText/N=text0 "\\s("+current+") sequiterpene"
		elseif(stringmatch(GetDimLabel(SOMparams,1,i),"benzene"))
			ModifyGraph rgb($current)=(52224,17408,0) // auburn
			AppendText/N=text0 "\\s("+current+") benzene"
		elseif(stringmatch(GetDimLabel(SOMparams,1,i),"toluene"))
			ModifyGraph rgb($current)=(52224,34816,0) // brown
			AppendText/N=text0 "\\s("+current+") toluene"
		elseif(stringmatch(GetDimLabel(SOMparams,1,i),"xylenes"))
			ModifyGraph rgb($current)=(16384,48896,65280) // light blue
			AppendText/N=text0 "\\s("+current+") xylenes"
		elseif(stringmatch(GetDimLabel(SOMparams,1,i),"n-alkane"))
			ModifyGraph rgb($current)=(0,0,65280) // blue
			AppendText/N=text0 "\\s("+current+") n-alkanes"
		elseif(stringmatch(GetDimLabel(SOMparams,1,i),"b-alkane"))
			ModifyGraph rgb($current)=(0,52224,0) // green
			AppendText/N=text0 "\\s("+current+") b-alkanes"
		elseif(stringmatch(GetDimLabel(SOMparams,1,i),"c-alkane"))
			ModifyGraph rgb($current)=(30464,30464,30464) // gray
			AppendText/N=text0 "\\s("+current+") c-alkanes"
		elseif(stringmatch(GetDimLabel(SOMparams,1,i),"naphthalene"))
			ModifyGraph rgb($current)=(39168,0,0) // deep red
			AppendText/N=text0 "\\s("+current+") naphthalene"
		endif
 	endfor
 	
	ModifyGraph lsize=2
	ModifyGraph fSize=14,axThick=1.5,standoff=0
	Label left "Mass Yield"
	Label bottom "[SOA] (ug/m3)\e"
	SetAxis bottom 0.1,*
	ModifyGraph log(bottom)=1
	ModifyGraph width=290,height=250

End

//************************************************************************************************
// Pop up functions and buttons
Function PopMenuProc_NOx(pa) : PopupMenuControl
	STRUCT WMPopupAction &pa

	SVAR LowOrHighNOx
	
	switch( pa.eventCode )
		case 2: // mouse up
			Variable popNum = pa.popNum
			String popStr = pa.popStr
			LowOrHighNOx = popStr
			break
		case -1: // control being killed
			break
	endswitch

	return 0
End

Function PopMenuProc_VWL(pa) : PopupMenuControl
	STRUCT WMPopupAction &pa

	SVAR VWLcondition
	switch( pa.eventCode )
		case 2: // mouse up
			Variable popNum = pa.popNum
			String popStr = pa.popStr
			VWLcondition = popStr
			break
		case -1: // control being killed
			break
	endswitch

	return 0
End

Function PopMenuProc_SingleCompound(pa) : PopupMenuControl
	STRUCT WMPopupAction &pa

	SVAR TheCompound
	switch( pa.eventCode )
		case 2: // mouse up
			Variable popNum = pa.popNum
			String popStr = pa.popStr
			TheCompound = popStr
			break
		case -1: // control being killed
			break
	endswitch

	return 0
End

Function SetParams_MPSOM(ba) : ButtonControl
	STRUCT WMButtonAction &ba

	switch( ba.eventCode )
		case 2: // mouse up
			SOM_SetupMulticomponent()
			break
		case -1: // control being killed
			break
	endswitch

	return 0
End

Function SetParams_SOM(ba) : ButtonControl
	STRUCT WMButtonAction &ba

	switch( ba.eventCode )
		case 2: // mouse up
			SOM_SetupSingleComponent()
			break
		case -1: // control being killed
			break
	endswitch

	return 0
End

Function Run_SOM_MP(ctrlName) : ButtonControl
	String ctrlName
	
	setdatafolder root:
	SOM_MC(SOMparams=SOMparams,InitialConcMatrix=InitialConcMatrix)

End
// END POPUP FUNCTIONS AND BUTTONS

//********************************************************************************************
// Recreate the SOM panel
Window sompanel_022716() : Panel
	PauseUpdate; Silent 1		// building window...
	NewPanel /W=(687,127,1194,793) as "sompanel_022716"
	ShowTools/A
	SetDrawLayer UserBack
	DrawRect 151,568,316,647
	DrawText 167,585,"\\f01Heterogeneous Chem"
	SetDrawEnv fillfgc= (65280,65280,48896)
	DrawRect 14,3,141,170
	DrawText 28,23,"\\f01General Operation"
	DrawRect 15,471,142,521
	DrawRect 14,173,141,317
	DrawText 16,190,"\\f01Kinetics and Oxidants"
	DrawText 172,485,"Model Parameters"
	SetDrawEnv fsize= 8
	DrawText 106,201,"0 = OH"
	SetDrawEnv fsize= 8
	DrawText 106,209,"1 = O3"
	SetDrawEnv fsize= 7
	DrawText 93,217,"Constant=0"
	SetDrawEnv fsize= 7
	DrawText 93,223,"Scaling=1"
	SetDrawEnv fsize= 7
	DrawText 93,229,"Experiment=2"
	SetDrawEnv fsize= 8
	DrawText 247,612,"0=1oxygen"
	SetDrawEnv fsize= 8
	DrawText 247,623,"1=multiple oxygens"
	SetDrawEnv fillfgc= (65280,65280,48896),fillbgc= (65280,65280,48896)
	DrawRect 150,4,315,182
	SetDrawEnv fillfgc= (57344,65280,48896)
	DrawRect 14,321,142,421
	DrawText 30,339,"\\f01Wall Loss Terms"
	DrawLine 21,272,135,272
	SetDrawEnv fillfgc= (57344,65280,48896)
	DrawRect 151,185,316,240
	DrawText 175,201,"\\f01Dynamic Partitioning"
	SetDrawEnv fillfgc= (65280,65280,48896),fillbgc= (65280,65280,48896)
	DrawRect 14,426,142,466
	DrawText 174,21,"\\f01Size Distribution Info"
	SetDrawEnv fillfgc= (65280,54528,48896),fillbgc= (65280,54528,48896)
	DrawRect 15,527,142,627
	DrawLine 156,140,308,140
	SetDrawEnv fillfgc= (65280,48896,48896)
	DrawRect 152,245,316,564
	DrawText 181,262,"\\f01SOM Parameters"
	DrawText 177,277,"(single component)"
	DrawLine 158,427,309,427
	SetDrawEnv fsize= 8
	DrawText 231,521,"0 = Epstein relationship"
	DrawText 58,545,"\\f01 CFSTR"
	SetDrawEnv fillfgc= (65280,48896,48896)
	DrawRect 321,7,490,129
	DrawText 340,25,"\\f01Multicomponent Setup"
	SetDrawEnv fillfgc= (65280,54528,48896),fillbgc= (65280,54528,48896)
	DrawRect 321,184,490,340
	DrawText 335,203,"\\f01Single Component Setup"
	SetDrawEnv fillfgc= (65280,65280,48896),fillbgc= (65280,65280,48896)
	DrawRect 321,133,490,179
	DrawText 370,150,"\\f01Fit Condition"
	SetVariable setvar0,pos={154,41},size={152,16},title="# Particles (p/cm3)"
	SetVariable setvar0,value= Nparticles
	SetVariable setvar1,pos={186,82},size={128,16},title="Density (g/cm3)"
	SetVariable setvar1,value= Density_Base
	SetVariable setvar2,pos={49,61},size={82,16},limits={-inf,inf,0},value= Nsteps
	SetVariable setvar3,pos={27,24},size={111,16},title="Timestep (s)"
	SetVariable setvar3,value= timestep
	SetVariable setvar4,pos={22,79},size={116,16},title="Dilution (%/hr)"
	SetVariable setvar4,value= DilutionVar
	SetVariable setvar5,pos={20,97},size={118,16},title="VP adjustment"
	SetVariable setvar5,value= logCstar_adjustmentfactor
	CheckBox check0,pos={193,587},size={59,14},title="ON/OFF",variable= hetchem
	SetVariable setvar6,pos={156,626},size={110,16},title="Gamma OH",value= gammaOH
	CheckBox check1,pos={17,475},size={123,14},title="\\Z10Sequential Partitioning"
	CheckBox check1,fSize=11,variable= SPM
	SetVariable setvar7,pos={192,342},size={121,16},value= ProbOx1
	SetVariable setvar8,pos={192,363},size={121,16},value= ProbOx2
	SetVariable setvar9,pos={192,384},size={121,16},value= ProbOx3
	SetVariable setvar08,pos={192,406},size={121,16},value= ProbOx4
	SetVariable setvar10,pos={44,277},size={95,16},title="[OH]",value= OHconc
	SetVariable setvar11,pos={16,231},size={123,16},title="ScalingFactor"
	SetVariable setvar11,value= OH_scale
	SetVariable setvar12,pos={158,278},size={155,16},title="# carbon atoms"
	SetVariable setvar12,value= Ncarbons
	SetVariable setvar13,pos={206,318},size={107,20},title="m\\Bfrag\\M"
	SetVariable setvar13,limits={0,inf,0.01},value= FragSlope
	SetVariable setvar14,pos={208,298},size={105,18},title="\\F'symbol'D\\F'arial'lVP"
	SetVariable setvar14,limits={0,3,0.01},value= delta_logCstar_perO
	SetVariable setvar15,pos={164,433},size={149,16},title="[VOC] (ppm)"
	SetVariable setvar15,value= ctot_ppm
	Button button0,pos={353,261},size={104,25},proc=Run_SOM,title="Run SOM",fSize=16
	Button button0,fStyle=1,fColor=(48896,52992,65280),valueColor=(65280,0,0)
	SetVariable setvar16,pos={231,453},size={82,16},title="H-per-O",value= H_per_O
	SetVariable setvar21,pos={206,473},size={107,16},title="H adjustment"
	SetVariable setvar21,value= Hadjustment
	SetVariable setvar22,pos={18,251},size={121,16},title="krxn(parent)"
	SetVariable setvar22,value= krxn_parent
	SetVariable setvar23,pos={19,192},size={82,16},title="OH or O3?",value= O3_yn
	SetVariable setvar25,pos={19,296},size={120,16},title="[O3] (ppb)"
	SetVariable setvar25,value= O3_conc
	SetVariable setvar26,pos={20,212},size={71,16},title="Method",fSize=12
	SetVariable setvar26,value= UseScalingForOxidant
	SetVariable setvar06,pos={18,42},size={120,16},title="RunTime (hrs)"
	SetVariable setvar06,value= MaxTime_Hours
	SetVariable setvar07,pos={160,605},size={82,16},title="Oxygens?"
	SetVariable setvar07,value= AddMultipleOxygens
	SetVariable setvar01,pos={152,162},size={160,16},title="Seed Conc (um3/cm3)"
	SetVariable setvar01,limits={-inf,inf,0},value= SeedVolConc
	Button button3,pos={326,289},size={156,25},proc=Run_SOM_Exp,title="Run SOM & Interp"
	Button button3,fSize=16,fStyle=1,fColor=(48896,52992,65280)
	Button button3,valueColor=(65280,0,0)
	CheckBox check3,pos={34,431},size={89,14},title="Oligomerization"
	CheckBox check3,variable= OligomerizationIsOn
	SetVariable setvar24,pos={21,447},size={115,16},title="krxn,base\\M"
	SetVariable setvar24,value= krxn_base_olig
	SetVariable setvar27,pos={19,340},size={117,16},title="WLR(g) (1/s)"
	SetVariable setvar27,value= gasWLR
	SetVariable setvar28,pos={20,400},size={117,16},title="WLR(p) (scale)"
	SetVariable setvar28,value= ParticleWallLoss
	CheckBox check5,pos={154,97},size={65,26},title="Absorbing\r seed?   "
	CheckBox check5,variable= AbsorbingSeed,side= 1
	SetVariable setvar29,pos={19,358},size={117,20},title="C\\Bwall\\M (mg/m3)"
	SetVariable setvar29,value= Cwall
	SetVariable setvar09,pos={18,115},size={120,16},title="Temperature (K)"
	SetVariable setvar09,value= Temp_G
	SetVariable setvar18,pos={20,133},size={118,16},title="Pressure (atm)"
	SetVariable setvar18,value= Pressure_G
	SetVariable setvar32,pos={193,493},size={120,16},title="DHvap (kJ/mol)"
	SetVariable setvar32,limits={0,inf,1},value= DHvap_G
	CheckBox check6,pos={17,491},size={80,14},title="Turn off SOA",variable= NoSOA
	SetVariable setvar02,pos={260,60},size={53,20},title="\\F'Symbol's\\F'Arial'\\Bg"
	SetVariable setvar02,value= SizeSpread
	SetVariable setvar03,pos={153,60},size={106,20},title="\\Z08D\\Bp,seed\\M\\Z08(nm)"
	SetVariable setvar03,value= SeedDiameter
	SetVariable setvar33,pos={167,203},size={127,16},title="Eqm=0; Dyn=1"
	SetVariable setvar33,value= KineticMassTransfer
	SetVariable setvar34,pos={185,221},size={109,16},title="Acc. coeff."
	SetVariable setvar34,value= alpha
	SetVariable setvar04,pos={163,144},size={149,16},title="Seed SA (um2/cm3)"
	SetVariable setvar04,limits={-inf,inf,0},value= SeedSurfaceArea
	CheckBox check7,pos={155,24},size={83,14},title="Polydisperse?"
	CheckBox check7,variable= polydisperse,side= 1
	SetVariable setvar35,pos={19,547},size={119,16},title="Bag Vol. (m3)"
	SetVariable setvar35,limits={0,inf,1},value= CSTR_Volume
	SetVariable setvar05,pos={43,604},size={93,16},title="Tau (h)"
	SetVariable setvar05,limits={-inf,inf,0},value= CSTR_ResidenceTime_hr
	SetVariable setvar37,pos={18,566},size={120,16},title="Flow Rate (lpm)"
	SetVariable setvar37,limits={0,inf,1},value= CSTR_FlowRate
	SetVariable setvar38,pos={41,585},size={97,16},title="# Times"
	SetVariable setvar38,limits={0,inf,1},value= CSTR_NumTimes
	Button button2,pos={349,97},size={115,25},proc=Run_SOM_MP,title="Run MC-SOM"
	Button button2,fSize=16,fStyle=1,fColor=(48896,52992,65280)
	Button button2,valueColor=(65280,0,0)
	SetVariable setvar19,pos={360,28},size={121,16},title="# Compounds"
	SetVariable setvar19,limits={1,10,1},value= nCompoundClasses
	PopupMenu popup0,pos={411,153},size={69,21},proc=PopMenuProc_NOx,title="NOx"
	PopupMenu popup0,mode=1,popvalue="low",value= #"\"low;high\""
	PopupMenu VWLcondition,pos={329,153},size={76,21},proc=PopMenuProc_VWL,title="VWL"
	PopupMenu VWLcondition,mode=3,popvalue="high",value= #"\"zero;low;high\""
	SetVariable setvar20,pos={354,47},size={127,16},title="Max alkane nC"
	SetVariable setvar20,limits={1,35,1},value= maxCforAlkanes
	Button button4,pos={373,69},size={65,25},proc=SetParams_MPSOM,title="Setup"
	Button button4,fSize=16,fStyle=1,fColor=(48896,52992,65280)
	Button button4,valueColor=(65280,0,0)
	PopupMenu Compound,pos={338,206},size={118,21},proc=PopMenuProc_SingleCompound,title="Compound"
	PopupMenu Compound,mode=5,popvalue="toluene",value= #"\"isoprene;aPinene;sesq;benzene;toluene;xylenes;naphthalene;n-alkane;b-alkane;c-alkane;isoprene_Alt\""
	Button button5,pos={353,232},size={104,25},proc=SetParams_SOM,title="Set Params"
	Button button5,fSize=16,fStyle=1,fColor=(48896,52992,65280)
	Button button5,valueColor=(65280,0,0)
	CheckBox check2,pos={17,506},size={84,14},title="First Gen Only"
	CheckBox check2,variable= FirstGenProductsOnly
	SetVariable setvar17,pos={243,22},size={65,16},title="# bins"
	SetVariable setvar17,limits={1,50,1},value= nSizeBinsG
	CheckBox check8,pos={348,320},size={114,14},title="Isoprene Is Special?"
	CheckBox check8,variable= IsopreneIsSpecial,side= 1
	CheckBox check9,pos={33,381},size={87,17},title="Krechmer C\\Bw\\M?"
	CheckBox check9,variable= Krechmer_Cw,side= 1
	SetVariable setvar30,pos={183,526},size={130,16},title="# O atoms in parent"
	SetVariable setvar30,value= Nox_precursor
	CheckBox check4,pos={20,152},size={112,14},title="Yield exclude POA?"
	CheckBox check4,variable= YieldCalc,side= 1
	SetVariable setvar31,pos={202,545},size={111,16},title="Cut off small nC"
	SetVariable setvar31,value= CutOffSmallStuff
	SetVariable setvar36,pos={167,122},size={130,16},title="NucleationTime (h)"
	SetVariable setvar36,value= NucleationTimeEmpirical
	SetVariable setvar39,pos={220,103},size={95,16},title=" Min Seed",value= minSeed
EndMacro

Window SizeDistribution() : Graph
	PauseUpdate; Silent 1		// building window...
	Display /W=(908.25,77.75,1276.5,382.25) dNdlogDp_time[0][*] vs Diameter_time[0][*] as "Size Distribution"
	AppendToGraph dNdlogDp_final vs Diameter_final
	ModifyGraph margin(top)=7,margin(right)=7
	ModifyGraph lSize=2
	ModifyGraph lStyle(dNdlogDp_final)=3
	ModifyGraph rgb(dNdlogDp_time)=(0,0,0),rgb(dNdlogDp_final)=(65280,0,0)
	ModifyGraph log(bottom)=1
	ModifyGraph tick=2
	ModifyGraph mirror=1
	ModifyGraph fSize=14
	ModifyGraph standoff=0
	ModifyGraph axOffset(left)=-2,axOffset(bottom)=-0.333333
	ModifyGraph axThick=1.5
	Label left "dN/dlogD\\Bp\\M (p cm\\S-3\\M)"
	Label bottom "Diameter (nm)"
EndMacro


function displayit()

	wave thewave = particlemolecules_time
	wave timew
	variable nrows = dimsize(thewave,0)
	variable ncols = dimsize(thewave,1)
	variable i, j
	
	make/o/d/n=(numpnts(timew)) thetestwave
	variable testvar
	variable counter=0
	string thename
	
	display
	for(i=0;i<ncols;i+=1)
		for(j=0;j<nrows;j+=1)
			testvar = thewave[j][i][0]
			
			if(numtype(testvar)!=2)
				appendtograph thewave[j][i][] vs timew
				counter+=1
			endif
		endfor
	endfor

ModifyGraph mode=7,hbFill=2
•ModifyGraph rgb=(26112,52224,0)
	
	wave thewave = gasmolecules_time
	counter=0
	for(i=0;i<ncols;i+=1)
		for(j=0;j<nrows;j+=1)
		testvar = thewave[j][i][0]
			
			if(numtype(testvar)!=2)
			appendtograph thewave[j][i][] vs timew
			if(i!=0 & j !=0) 
				thename = "GasMolecules_Time#"+num2istr(counter)
//				print thename
				ModifyGraph rgb($thename)=(16384/i,168384/i,65280/i)
			endif
			counter+=1
			endif
		endfor
	endfor
	ModifyGraph hbFill=2,toMode=2

end

Function doit()
	
	wave thewave = totalcarbon_time_ppb
	string wvstr = "totalcarbon_time_ppb"
	variable nrows = dimsize(thewave,0)
	variable ncols = dimsize(thewave,1)
	variable i, j
	make/o/d/n=(dimsize(thewave,2)) tempwave

	Make/O/N=26/FREE pc_red = {  65280, 50000, 39168, 52224, 65280, 0, 16384, 0, 0, 0, 0, 0, 0, 32768, 39168, 30464, 21760, 0, 34816, 43520, 52224, 39168, 29440, 65280, 26112, 0}
	Make/O/N=26/FREE pc_green={   0, 17480, 0, 0, 32512, 0, 48896, 17408, 0, 65280, 52224, 26112, 65280, 65280, 39168, 30464, 21760, 0, 34816, 43520, 34816, 0, 0, 43520, 17408, 0}
	Make/O/N=26/FREE pc_blue={    17408, 0, 15616, 20736, 16384, 65280, 65280, 26112, 52224, 65280, 0, 0, 0, 65280, 0, 30464, 21760, 0, 34816, 43520, 0, 31232, 58880, 0, 0, 0}

	pc_red = {0, 25535, 2, 0, 39321, 48059, 65535, 0, 16385, 65535, 0, 65535, 2, 0, 39321, 48059, 65535, 0, 16385, 65535}
	pc_green = {0, 16385, 39321, 0, 1, 48059, 32768, 65535, 65535, 32768, 0, 16385, 39321, 0, 1, 48059, 32768, 65535, 65535, 32768}
	pc_blue = {0, 16385, 1, 65535, 31457, 48059, 32768, 0, 65535, 58981, 0, 16385, 1, 65535, 31457, 48059, 32768, 0, 65535, 58981}

	string this_trace
	variable counter = 0	
	display
	for(i=nrows-1;i>=0;i-=1)
		for(j=ncols-1;j>=0;j-=1)
			tempwave = thewave[i][j][p]
			wavestats/q tempwave
			if(numtype(V_avg)!=2)
				appendtograph thewave[i][j][]
				if(counter==0)
					this_trace = wvstr
				else
					this_trace = wvstr + "#" + num2istr(counter)
				endif
				ModifyGraph rgb($this_trace)=(0,0,0)
				ModifyGraph plusrgb($this_trace)=(pc_red[i],pc_green[i],pc_blue[i])
				ModifyGraph negrgb($this_trace)=(pc_red[i],pc_green[i],pc_blue[i])
				counter += 1
			endif
		endfor
	endfor
	
	ModifyGraph mode=7,useNegRGB=1,usePlusRGB=1,hbFill=4,toMode=2,rgb=(0,0,0)
	ModifyGraph toMode=3

ENd


