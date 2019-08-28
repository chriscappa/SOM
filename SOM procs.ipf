#pragma rtGlobals=1		// Use modern global access method and strict wave access.
#include <Global Fit 2>
#include <KBColorizeTraces>

//*************************************************************************************************************************************************
Function TopOfSOM_toolkit_Procs() 
	//The only purpose of this function is to be the reference point for opening this procedure window from the menu bar.
	//=> leave this function at the top of the procedure window
End

//*********************************************************************
menu  "SOM"
	
	"Show SOM procs"
	
	Submenu "SOM"
		"Run SOM", SOM_v1()
		"Recreate Fit Window", MakeGraph_FitWindow()
		"Create SOM Panel", MakePanel_SOMpanel()
		"Create Fit Table", MakeTable_FitSetup()
		"Show SOM v1 proc"
	end
	Submenu "MultiCompnent"

	end
end

Function ShowSOMprocs()
	DisplayProcedure "TopOfSOM_toolkit_Procs"
End

Function ShowSOMv1proc()
	DisplayProcedure "SOM_v1"
End

//************NOTES**********************//
// v7.4.3 08/26/19
//	1. Added a number of waves to SOM_KillWaves()...most noted with 7.4.3
//	2. Replace all instances of Noxygens_parent (local variable) with Nox_precursor (global)
//	3. Major change: vapor wall losses now can be treated using the Zaveri (MOSAIC) method by setting VWL_method == 1 (set to default)
//	4. Lots of comments added
//	5. Added dropdown menu, with various popups/recreation functions
// v7.4.2 05/26/19
// 	1. Updated RunInfo and RunParameters (in the associated sub-folders) for low-NOx toluene (2013) based on observed size distributions, number concentrations, and widths
//	2. Created function "LoadFromCaltech_seed" that will set the seed diameter, spread, and number concentration
// 	3. Created "RecalculateChi2_single()" function to calculate ChiSq for a given experiment; will discard points with signals much lower than the peak (as these cause problems sometimes)
//	4. Updated "RecalculateChi2()" function to use interpolated output from SOM
// V 7.4.1 - 1/16/19
//	1. Updated to allow for T-dependence, based on Epstein relationship, modified by Cappa/Jimenez relationship
// 	2. Minor updates to how kinetic (krxn) is calculated to account for temperature
// V 7.4 - 11.17.18
// 	1. Updated SOM_V1 to allow use of either Euler method (original) or Zaveri method (new) for gas-particle dynamic partitioning
//		MCSOM not updated yet
//	    - This is now the default for gas-particle kinetic partitioning; gives results (usually) very similar to equilibrium value if alpha = 1
//	    - For some systems, can run into trouble at early times in the simulation if the surface area is not sufficiently large
//	    - Can lead to steps in concentration, but typically end up with correct final answer (if simulation runs long enough)
//	    - Can manage this problem by setting a minimum absorbing seed concentration, b/c problems are with mole fraction calculation when concentrations are really small
// 	    - This doesn't fix everything. But it helps a lot
//	2. Update how mass transfer time step is dealt with for Euler method. Make this timestep vary as 1/(alpha*log(seedsurfacearea)). But,
//	    allow kinetic loop to reset if the initial guess is not good. Will break after four attempts.
// 	     - turns out this was unnecessary, as the MOSAIC method now works. But it illustrates some differences between the two
//	3. Update vapor wall loss to have it run in a loop, depending on the value of the kwall_on. When larger than 1e-4, will run with multiple steps
// V 7.3.8
//	1. added calculations of "HOM" concentrations. This is a starting point and not currently used. User specifies (hard-wired) volatility limit on HOM
//	2. added empirical method of treating nucleation. User specifies as a global variable the "NucleationTimeEmpirical" in hours. This is the time at which
//	   it is assumed that nucleation occurs. Prior to this it is assumed that the particle concentration is zero. At this point, 'seed' aerosol is turned on with the
//	  diameter, width and number concentration specified
// 	3. Fixed bug with logCstar_adjustmentfactor in SOM_MC; previously had logCstar_matrix *= logCstar_adjustmentfactor 
//	  but this should be logCstar_matrix += logCstar_adjustmentfactor; this was really messing things up
// V 7.3.7
//	- updated treatment of precursors with oxygens. this was not robust, and variable on panel added
//	- updated "VP adjustment" on panel to be +/- rather than *; see variable logCstar_adjustment
// 	- added ability to exclude "POA," that is initially partitioned mass from low-volatility parent compounds, from yield calc
// 		-- a button was added on the panel
//	- updated function "SOM_SetupSingleComponent()" so that for alkanes the # of carbon atoms is not changed
//		when one clicks the "Set Params" button
// V 7.3.4
// 	- added ability to use Krechmer (2016) parameterization for Cw, the effective absorbing wall massuc
//	- added ability to use alternate isoprene scheme, which assumes that krxn(C5O1-C5O4) = krxn(isoprene)
// V 7.3.0
// 	- added ability to use a distribution of initial VOC's of a given type (i.e. C5-C15)
// Version 7.2
// Update 09/29/14
//	- Improved ability to deal with particle wall deposition
//	- added ability to use a polydisperse (but log normal) particle size distribution
//		- the number of bins is hardwired in the code. More bins = longer calculation time.
//	- added ability to keep track of wall-corrected SOA from particle deposition
//		by calculating Coa_time_wallcorr in addition to Coa_time, which is now just the suspended
//		particle mass concentration
// Version 7.1
// Update 08/26/14
// 	- attempted to use built in IGOR ODE solver to speed up dynamic partitioning
// 		this was a waste of time, as it actually slowed things down with no noticable improvement
//		in terms of accuracy. This has been commented out.
//	- Modified "brute force" dynamic partitioning code to improve speed by "precalculating" a number of things
//		that do not change much over a given period
//		Nearly a factor of 2 increase in speed.
//	- Modified timesteps associated with "brute force" dynamic partitioning code to better account for different
//		needs when different alpha's are used (larger alphas require smaller timesteps)
//		This leads to substantial increases in speed.
// 04/18/14
// 	- v7 created. No major updates from v6c. Just cleaning things up for the most part
// 	- added ability to fit data twice (SOM_BatchFit_FitTwice). This takes the results from an initial fit and uses them
//		to perform a second fit. Meant to deal with situations where the final fit is not the "best" that the fitting algorithm can do

	variable/G Nparticles // number of particles (for heterogeneous chemistry)
	variable/G Density_Base // density of particles (for heterogeneous chemistry)
	variable/G Pfrag_type = 1 // Hardwired to 1 for v7.3.3 // 0 = cfrag*Nox; 1 = (O:C)^mfrag
	variable/G Nsteps // number of steps
	variable/G timestep // timestep in seconds
	variable/G DilutionVar // dilution in percent per hour
	variable/G small_fragments = 1 // Hardwired to 1 for v7.3.3 // 0 for random; 1 for equal probabilities; 2 for small fragments (loose 1 carbon)
	variable/G logCstar_adjustmentfactor // multiplicative factor to adjust base vapor pressure
	variable/G hetchem // 0 for no heterogeneous; 1 to include
	variable/G gammaOH // reactive uptake coefficient
	variable/G SPM // 0 = no; 1 = use sequential partitioning model
	Variable/G ProbOx1, ProbOx2, ProbOx3, ProbOx4 // probability of adding X oxygens--> is renormalized later
	variable/G OHconc // [OH] in molecules/cm^3
	variable/G OH_scale // >0 if OH decays with time --> OH = OHconc_initial*(1-(m*timestep/timemax))^oh_scale
	variable/G Ncarbons // number of carbon atoms in parent molecule
	variable/G FragSlope // either mfrag or cfrag, dependnig on fragmentation method
	variable/G delta_logCstar_perO // change in logCstar per oxygen added
	variable/G ctot_ppm // concentration of parent hydrocarbon (gas + particle) in ppm
	variable/G H_per_O // number of H atoms lost per O added (van Krevelen slope)
	variable/G Hadjustment // adjustment to account for cases where there are fewer H's than a saturated hydrocarbon
	variable/G SaveVar // 1 to save data into "SaveFolderStr"
	string/G SaveFolderStr // where to save data; appends onto "root:"; will create the folder
	variable/G krxn_parent // set specific value for parent compound
	variable/G krxn_method = 4 // Hardwired to 4 in v7.3.3 // method to use for specifying krxn matrix; 0 is legacy, 1 is newish, 2 is newer (09/02/13)...4 is the one to use, and detailed in Zhang et al. (2014), PNAS
	variable/G O3_yn // O3 reaction or OH reaction
	variable/G O3_conc // O3 concentration in ppb
	variable/G UseScalingForOxidant // decide how to treat OH or O3 concentration time dependence
	string/G rxnwavestr
	variable/G MaxTime_Hours
	variable/G AddMultipleOxygens
	variable/G SeedVolConc
	variable/G krxn_base_olig
	variable/G OligomerizationIsOn
	variable/G gasWLR = 0	// gas-phase wall loss in per second
	variable/G Cwall			// "effective" wall OA concentration in mg/m^3, added 8/27/2013
	variable/G kwall_gas_scaling	// scaling of gasWLR; 0 = composition independent; > 1 = increases with added oxygen. 0 < x < 1 = decreases with added oxygen --> added 8/27/2013
	variable/G ParticleWallLoss = 0	// 0 = no wall loss, 1 = calculate wall loss --> added 08/20/2013
	variable/G EqmMethod = 1 // Hardwired to 1 for v7.3.3 // 0 = molecules; 1 = mass  --> added 07/20/2013
	variable/G AbsorbingSeed 		// 0 = non-absorbing, 1 = absorbing
	variable/G Temp_G	// temperature in K		// added 09/03/13
	variable/G Pressure_G	// pressure in atm		// added 09/03/13
	variable/G DHvap_G		// Enthalpy of vaporization for SOM species; 0 = Epstein relationship; 1 = constant value in kJ/mol; added 09/03/13
	Variable/G NoSOA=0		// flag to turn off equilibrium partitioning (0 = partitioning as normal, 1 = no partitioning); added 10/13/13
	variable/G SizeSpread		// log-normal standard deviation for seed particle distribution; added 10/16/13
	variable/G SeedDiameter	// seed particle diameter in nm; added 10/16/13 as an alternative to the seed concentration
	variable/G KineticMassTransfer		// option to treat gas-particle partitioning dynamically, 0 = no, 1 = yes, added 11/12/13
	variable/G alpha	// mass accommodation coefficient for dynamic gas-particle partitioning, added 11/12/13
	variable/G SeedSurfaceArea	// seed surface area in um^2/cm^3, output, added 11/12/13
	variable/G Polydisperse // 0 = mono; 1 = poly
	variable/G CSTR_Volume // m^3; volume of the reactor for CSTR; added 6/23/15
	variable/G CSTR_FlowRate // lpm; flowrate of air through the reactor for CSTR; added 6/23/15
	variable/G CSTR_ResidenceTime_hr // hrs; (BagVolume/(FlowRate*1e-3/60))/3600; added 6/23/15
	variable/G CSTR_ResidenceTime_sec // sec; BagResidenceTime_hr*=3600
	variable/G CSTR_NumTimes // number of residence times to consider (# of times model will be run)
	variable/G nSizeBinsG //= 7 // number of size bins for polydisperse calculations; made a global variable 6/23/15
	// some new additions for v7.3.3
	variable/G nCompoundClasses // added for v7.3.3 // # of compounds in multicomponent run
	string/G LowOrHighNOx // added for v7.3.3 // low or high, for selecting from fit parameters
	string/G VWLcondition // added for v7.3.3 // zero, low (kwall = 1e-4) or high (kwall = 2.5e-4)
	variable/G maxCforAlkanes // added for v7.3.3 // maximum number of carbon atoms for multicomponent simulation of alkanes
	string/G TheCompound // added for v7.3.3 // the compound of interest to set for single component calcs
	variable/G FirstGenProductsOnly // added for v7.3.3 // to only have first generation products formed
	variable/G IsopreneIsSpecial
	variable/G Krechmer_Cw
	variable/G Nox_precursor // added as global variable for v7.3.7
	variable/G YieldCalc // added for v7.3.7
	variable/G CutOffSmallStuff
	variable/G NucleationTimeEmpirical
	variable/G DynamicPartitioningMethod // 0 = Euler (legacy for SOM_v1); 1 = MOSAIC
	variable/G/D minSeed
	variable/G Ea_krxn
	
Function SetGlobalVarsToDefault()

	NVAR Nparticles
	Nparticles = 40000 // number of particles (for heterogeneous chemistry)
	NVAR Density_Base // density of particles (for heterogeneous chemistry)
	Density_Base = 1 // g/cm^3
	NVAR Pfrag_type // 0 = cfrag*Nox; 1 = (O:C)^mfrag
	Pfrag_type = 1
	NVAR Nsteps // number of steps
	Nsteps = 500
	NVAR timestep // timestep in seconds
	timestep = 50
	NVAR DilutionVar // dilution in percent per hour
	DilutionVar = 0
	NVAR small_fragments // 0 for random; 1 for small fragments
	Small_Fragments = 0
	NVAR logCstar_adjustmentfactor // multiplicative factor to adjust base vapor pressure
	logCstar_adjustmentfactor = 1
	NVAR hetchem // 0 for no heterogeneous; 1 to include
	hetchem = 0
	NVAR gammaOH // reactive uptake coefficient
	gammaOH = 1
	NVAR SPM // 0 = no; 1 = use sequential partitioning model
	SPM = 0
	NVAR ProbOx1, ProbOx2, ProbOx3, ProbOx4 // probability of adding X oxygens--> is renormalized later
	ProbOx1 = 0.25
	ProbOx2 = 0.25
	ProbOx3 = 0.25
	ProbOx4 = 0.25
	NVAR OHconc // [OH] in molecules/cm^3
	OHconc = 2e6
	NVAR OH_scale // >0 if OH decays with time --> OH = OHconc_initial*(1-(m*timestep/timemax))^oh_scale
	OH_scale = 0
	NVAR Ncarbons // number of carbon atoms in parent molecule
	Ncarbons = 10
	NVAR FragSlope // either mfrag or cfrag, dependnig on fragmentation method
	FragSlope = 0.3
	NVAR delta_logCstar_perO // change in logCstar per oxygen added
	delta_logCstar_perO = 1.8
	NVAR ctot_ppm // concentration of parent hydrocarbon (gas + particle) in ppm
	ctot_ppm = 0.1
	NVAR H_per_O // number of H atoms lost per O added (van Krevelen slope)
	H_per_O = 1
	NVAR Hadjustment // adjustment to account for cases where there are fewer H's than a saturated hydrocarbon
	Hadjustment = 0
	NVAR SaveVar // 1 to save data into "SaveFolderStr"
	SaveVar = 0
	string/G SaveFolderStr // where to save data; appends onto "root:"; will create the folder
	SaveFolderStr = "myFolder"
	NVAR krxn_parent // set specific value for parent compound
	krxn_parent = 1e-10
	NVAR krxn_method // method to use for specifying krxn matrix; 0 is legacy, 1 is newish, 2 is newest (09/02/13)
	krxn_method = 4
	NVAR O3_yn // O3 reaction or OH reaction
	O3_yn = 0
	NVAR O3_conc // O3 concentration in ppb
	O3_conc = 0
	NVAR UseScalingForOxidant // decide how to treat OH or O3 concentration time dependence
	UseScalingForOxidant = 0
	string/G rxnwavestr
	rxnwavestr = ""
	NVAR MaxTime_Hours
	maxtime_hours = 20
	NVAR AddMultipleOxygens
	AddMultipleOxygens = 0
	NVAR SeedVolConc
	SeedVolConc = 10
	NVAR krxn_base_olig
	krxn_base_olig = 0
	NVAR OligomerizationIsOn
	OligomerizationIsOn = 0
	NVAR gasWLR 
	gasWLR = 0	// gas-phase wall loss in per second
	NVAR Cwall			// "effective" wall OA concentration in mg/m^3, added 8/27/2013
	Cwall = 10
	NVAR kwall_gas_scaling	// scaling of gasWLR; 0 = composition independent; > 1 = increases with added oxygen. 0 < x < 1 = decreases with added oxygen --> added 8/27/2013
	kwall_gas_scaling = 0
	NVAR EqmMethod 
	EqmMethod = 1 // 0 = molecules; 1 = mass  --> added 07/20/2013
	NVAR ParticleWallLoss 
	ParticleWallLoss = 0	// 0 = no wall loss, 1 = calculate wall loss --> added 08/20/2013
	NVAR AbsorbingSeed 		// 0 = non-absorbing, 1 = absorbing
	AbsorbingSeed = 0
	NVAR Temp_G	// temperature in K		// added 09/03/13
	Temp_G = 298
	NVAR Pressure_G	// pressure in atm		// added 09/03/13
	Pressure_G = 1
	NVAR DHvap_G		// Enthalpy of vaporization for SOM species; 0 = Epstein relationship; 1 = constant value in kJ/mol; added 09/03/13
	DHvap_G = 0
	NVAR NoSOA
	NoSOA = 0		// flag to turn off equilibrium partitioning (0 = partitioning as normal, 1 = no partitioning); added 10/13/13
	NVAR SizeSpread		// log-normal standard deviation for seed particle distribution; added 10/16/13
	SizeSpread = 1.1
	NVAR SeedDiameter	// seed particle diameter in nm; added 10/16/13 as an alternative to the seed concentration
	SeedDiameter = 70
	NVAR KineticMassTransfer		// option to treat gas-particle partitioning dynamically, 0 = no, 1 = yes, added 11/12/13
	KineticMassTransfer = 0
	NVAR alpha	// mass accommodation coefficient for dynamic gas-particle partitioning, added 11/12/13
	alpha = 1
	NVAR SeedSurfaceArea	// seed surface area in um^2/cm^3, output, added 11/12/13
	SeedSurfaceArea = 1000
End

//******************************************************************************************
// The Statistical Oxidation Model
// Actually started keeping notes on updates on 07/07/2013
// 10/16/2013: 1. Saved as v6 (from v5)
//10/13/2013: Added capability of turning off equilibrium partitioning to particles so that only gas-phase chemistry and gas wall-loss can be considered
// 08/20/2013: 1. fixed gas-phase wall-loss...it was in the wrong spot and didn't actually do anything. Values around 1e-6/s affect the results. 
//			  2. Added particle wall loss option, based on wall-losses in Loza et al. Atmos. Chem. Phys., 12, 151167, 2012
// 07/20/2013: Added option to perform eqm. calculations based on mass concentrations, as opposed to molecular concentrations
//			it is not clear which is "better" to use, but some (e.g. Donahue) prefer mass, with the argument that it better captures the actual physicallity of the situation, which is that bigger molecules take up more space.
//			search for string "units" to update and select which method to use
// 07/11/2013: MAJOR update
//			1. it was determined that there was an error in how the fragmentation was being applied. this has now been corrected so that
//			fragmentation probabilities depend on the oxygen content of the PARENT species. Oops. Corrected both in SOM_v1 and in SOM_Heterogeneous
//			2. added multithreading to SOM_v1 and heterogeneous that cut off serious time. Hooray.
// 07/07/2013: added an option to account for wall-losses of gas-phase species
//			assumes a first order loss, with gasWLR in units of per second.


//********************************************************************************
// This is the main SOM function. You can either input fit coefficients, or the panel values will be used by default
FUNCTION SOM_v1([fitcoefs,fittime,quiet,cleanup])
	wave fitcoefs		// [0] = FragSlope
				 	// [1] = delta_logCstar_perO
				 	// [2] = Pfunc[0]
				 	// [3] = Pfunc[1]
				 	// [4] = Pfunc[2]
				 	// [5] = Pfunc[3]
				 	// [6] = GasWLR	// gas-phase wall loss rate
	wave fittime		// time wave to use, if fitting data --> typically the experimental time (e.g. "Time_experiment")
	variable quiet		// 0 (default) = print, 1 = limited printing
	variable cleanup	// 0 = do not kill waves; 1 = kill a bunch of waves at the end (default)
	
	// everything should happen in root folder...yes, this is annoying but more trouble to rewrite everything
	setdatafolder root:

	if(ParamIsDefault(quiet)) // default is to print "full" information
		quiet = 0		
	endif
	
	if(ParamIsDefault(cleanup)) 
		cleanup = 1
	endif
	
	// specify SOM parameters
	NVAR FragSlope // either mfrag or cfrag, depending on fragmentation method --> from the panel
	NVAR delta_logCstar_perO // change in logCstar per oxygen added --> from the panel
	NVAR ProbOx1, ProbOx2, ProbOx3, ProbOx4	// probability of adding some number of oxygen atoms
	NVAR gasWLR			// gas phase wall loss first order rate coefficient (1/s)
	variable numFitCoefs = numpnts(FitCoefs)	// number of elements in FitCoefs; needed to deal with more than standard 6 FitCoefs

	if(ParamIsDefault(fitcoefs)) // If no fitcoefs wave is used in input, use values from panel. 
		// use values from panel; do not allow negative oxygen additiona probabilities
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
	
//	NVAR MaxTime_hours // length of simulation, in hours // v7.4 this is being moved below to allow for timestep to change, as necessary for dynamic partitioning
//	nsteps = ceil(MaxTime_hours*60*60/timestep)+1	// total number of time steps
	
	// Select whether to interpolate the results to an experimental time base
	variable InterpToExperiment
	if(ParamIsDefault(fittime))
		InterpToExperiment = 0	// do nothing and use default number of points
	else
		InterpToExperiment = 1	// do something, and interpolate at the end to n = numpnts(fittime) points
	endif

	// Set this variable to control number of iterations during fitting
	// This doesn't always seem to do anything...probably needs to be set somewhere else to be effective
	variable V_FitMaxIters = 10 // max number of iterations
	
	if(quiet==0)
		print "//"
		print "Model Run at " + time() + " on " + date()
	endif
	
// Access variables from panel, and create a few more variables for use
// Oligomerization
	NVAR OligomerizationIsOn// oligomerization is on (1) or off (0)
	NVAR krxn_base_olig		// oligomerization bimolecular rate coefficient
// Gas-phase wall loss
	variable gasWLmethod = 1 // 0 = irreversible uptake, 1 = reversible uptake
	NVAR Cwall				// "effective" wall concentration; typical value = 10 mg/m3
	NVAR kwall_gas_scaling	// composition dependent wall loss rate coefficient (1/s); value for Caltech chamber ~2.5x10^-4 1/s
	NVAR noSOA 			// optional value to set C* values very large, to follow only the gas-phase chemistry
	NVAR Krechmer_Cw		// optional value to use the Krechmer wall loss parameterization (see below)
// Precursor properties and defining range of "species" to be considered
	NVAR Ncarbons // number of carbon atoms in parent molecule
	NVAR Nox_precursor // = 0 // number of oxygens in precursor species
	variable nC
	variable nO 
	nC = nCarbons
	nO = ceil(20/(nC^0.33)) //+ 2 // Max number of oxygen atoms that can be added
	nO = 7 // 10 //8
	if(OligomerizationIsOn==1)
		nC = nC*2
		nO = nO*2 // max number of oxygen atoms that can be added
	endif
	variable Nspecies = nC*nO		// added 08/26/14 for dynamic partitioning
	
	// precursor concentrations (for single precursor)
	NVAR ctot_ppm // concentration of parent hydrocarbon (gas + particle) in ppm
	variable Ctot_init // initial mass loading (in ug/m^3) of the total organic (gas + particle) material

	// info about particles
	NVAR Nparticles // total particle concentration; #/cm^3
	variable Np = Nparticles*1e6 // total particle concentration; #/cm^3 --> #/m^3
	NVAR Density_Base // particle material density; g/cm^3
	variable density = density_base*1000 // particle material density; g/cm^3 --> kg/m^3
	// parameters related to kinetics
	NVAR krxn_parent // reaction rate coefficient for precursor species; if krxn_parent == 0, treat as hydrocarbon. Otherwise use value
	NVAR krxn_method // Hardwired to 4 in v7.3.3 // method to use for specifying krxn matrix; 0 is legacy, 1 is newish, 2 is newer (09/02/13)...4 is the one to use, and detailed in Zhang et al. (2014), PNAS
	NVAR IsopreneIsSpecial // assume that the kOH values for C5O1-C5O4 are the same as k_isoprene
	NVAR O3_yn // 0 = OH is oxidant; 1 = O3 is oxidant and there *should* be no 2nd generation chemistry (not fully tested throughout)
	NVAR O3_conc // O3 concentration (ppb)
	NVAR usescalingforoxidant // 0 = use constant [OH]; 1 = allow for time varying OH based on scaling factor; 2 = interpolate observed [OH] or [O3] to model timebase
		// note: usescalingforoxidant == 2 requires that both OH_exp (observed concentration; molec/cc) and OH_exp_time (experimental time, hrs) waves exist
	
	// some constants
	variable m, n, o, k, i, j // general variables
	variable Na = 6.022e23 // molecules/mol...Avagadros Number
	NVAR Temp_G // temperature, in Kelvin, as global variable; specified in panel (may be adjusted throughout with time-varying T)
	NVAR Pressure_G // pressure, in atm, as global variable ; specified in panel
	variable Tvar = Temp_G //298.15 // kelvin
	NVAR Ea_krxn // activation energy for reactions, kJ/mol
	variable Pressure = Pressure_G * 760 * 133.322 // Pa
	variable NumberDensityAir = Pressure / (1.381e-23 * Tvar) // molecules per m^3
	
	// Equilibrium method (mass or molecules)
	NVAR EqmMethod
	string units
	if(EqmMethod == 0) // "molecules" or "mass"; affects how eqm. partitioning is performed. See EqmCalc().
		units = "molecules"	
	else
		units = "mass" // this is the typical value
	endif
	
	NVAR DynamicPartitioningMethod // 0 = Euler, 1 = Zaveri --> specifies the method to use for dynamic partitioning
	DynamicPartitioningMethod = 1 // This is hardwired to use the (more robust) Zaveri method
	variable Cheat = 0 // temporary variable if you want to test stuff

//// 1. Decide here how to input initial concentration	

// 2. Specify various properties

  	variable Diff_CO2 = 0.138e-4 // m^2/s; diffusitivity of CO2 (the reference value for other diffusivities)
  	variable MW_CO2 = 44	// g/mol; MW of CO2
  	
// 2a. Populate O:C, H:C, Oxidation State and MW Matrices 
	// Parameters for Van Krevelen
	NVAR H_per_O // # of H atoms lost per O added --> means that van Krevelen slope is specified, not a real result
	NVAR Hadjustment // rarely used; additive scaling factor to account for e.g. double bonds, rings, in calculation of MW (+1 = MW are decreased by 1)
	//
	SOM_CreateAtomMatrices(nC,nO,H_per_O,Hadjustment) // create...
	wave C_matrix
	wave O_matrix
	wave H_matrix
	wave OC_matrix
	wave HC_matrix
	wave Ox_matrix
	wave MW_matrix

// 2b. Parameters for size distribution; added 10/16/13; updated 09/29/14
	NVAR SizeSpread	// log-normal standard deviation (>= 1.0) for seed particle distribution; added 10/16/13
	NVAR SeedDiameter	// seed particle diameter in nm; added 10/16/13 as an alternative to the seed concentration
  	NVAR SeedVolConc 	// um^3/cm^3 --> calculated from size distribution
	NVAR SeedSurfaceArea // (um^2/cm3 --> calculated from size distribution
	NVAR polydisperse // 0 = monodisperse particles; 1 = polydisperse
	NVAR nSizeBinsG // number of size bins
	variable nSizeBins = nSizeBinsG 
	variable index_size // counter for the size bins
	variable dlogDp	// dlogDp associated with generated size distribution
	variable DpStart = SeedDiameter
  	variable VolumeSeed = (pi/6)*((DpStart*1e-9)^3)	// m^3/particle; seed volume
  	variable VolumeOrgPerParticle	// the organic volume per particle
	
	if(polydisperse==0)
		nSizeBins = 1
	endif

	dlogDp = makelognormaldistn(SizeSpread,SeedDiameter,Nparticles,nbins = nsizebins) // create size distribution arrays and matrices
	wave Diameter	// nm; a wave with the midpoint diameters
	wave dNdlogDp	// p/cm^3
	wave dSdlogDp
	wave dVdlogDp
	variable NpCalc // number concentration
	make/o/d/n=(nsizebins) logDp // logDp = log(Diameter)
	logDp = log(diameter)
	make/o/d/n=(nSizeBins) dVdlogDp_norm // normalized volume concentration
	Integrate/METH=1 dNdlogDp/X=logDp/D=INT_Result // integrate number distribution
	NpCalc = INT_Result[nsizebins-1] // # concentration; this should be the same as the specified # of particles
	Integrate/METH=1 dVdlogDp/X=logDp/D=INT_Result // integrate volume distribution
	SeedVolConc = INT_Result[nsizebins-1]/1e9 // um^3/cm^3; seed volume concentration
	dVdlogDp_norm = dlogdp*1e-9*dVdlogDp / SeedVolConc // normalize the seed volume concentration
	wavestats/q dVdlogDp_norm
	dVdlogDp_norm /= V_sum // normalized volume, according to the number of particles and size
	Integrate/METH=1 dSdlogDp/X=logDp/D=INT_Result
	SeedSurfaceArea = INT_Result[nsizebins-1]/1e6 // um^2/cm^3
	
  	make/o/d/n=(nSizeBins) VolumeSeedPerBin = (pi/6)*((Diameter[p]*1e-9)^3)	// m^3/particle; volume of a single particle in each size bin
	variable SeedMW = 250 // g/mol; DOS = 427
	variable SeedDensity = 1 // g/cm3
  	make/o/d/n=(nSizeBins) MoleculesSeedPerBin = VolumeSeedPerBin * 1e6 * SeedDensity * (1/SeedMW) * Na * (dlogDp*dNdlogDp) // molecules/cm3 = m^3/p * cm3/m3 * g/cm3 * mol/g * molecules/mol * p/cm^3
	note/k MoleculesSeedPerBin "molecules/cm3-air"
  	make/o/d/n=(nSizeBins) MassSeedPerBin = VolumeSeedPerBin * 1e6 * SeedDensity * (dlogDp*dNdlogDp) * 1e6 // g/cm3 = m^3/p * cm3/m3 * g/cm3 * p/cm^3 * cm3/m3
  	note/k MassSeedPerBin "g/m3-air"
  	// set a minimum seed concentration for any calculation (from v7.4.0)
  	NVAR minSeed // ug/m3
	make/o/d/n=(nSizeBins) MassSeedPerBin_noseed = minSeed*1e-6*dVdlogDp_norm//(pi/6)*((10*1e-9)^3)* 1e6 * SeedDensity * (dlogDp*dNdlogDp) * 1e6 // use this if there is no seed to prevent problems with dynamic partitioning
  	make/o/d/n=(nSizeBins) VolumeOrgPerBin // m^3/particle
  	make/o/d/n=(nSizeBins) dDpdt_on, dDpdt_off
	make/o/d/n=(nSizeBins) NumParticlesByBin = dNdlogDp*dlogDp // added with v7.4 for Zaveri method; product of dNdlogDp and the width of each bin // p/cm3
  	variable/D Dk10 = 3.75 // added with v7.4, for calculating Kelvin term: will use Kelvin = 10^(Dk10/Dp), where Dk10 = 3.75 nm, from Trostl et al., Nature, 2015

  	if(Polydisperse==0)
  		SeedVolConc = VolumeSeed*NParticles*1e18	// um^3/cm^3
		SeedSurfaceArea = (4*pi*(DpStart*1e-3/2)^2)*(Nparticles) // um^2/cm^3
	endif

	// Parameters for Particle Wall Loss; Added 08/20/13
	// empirical expression determined from Loza et al., ACP, 2012
	NVAR ParticleWallLoss	// 0 = no loss, 1 = loss
	variable/D particleWLR	// size-dependent wall loss rate in per second; calculated at each model time-step below
	make/o/d/n=(numpnts(dndlogdp)) particleWLR_wave
	note particleWLR_wave "units = 1/s; determined from Loza et al, 2012 for Caltech chamber"
	
// Specify time step, and some conditionals to help with stability and setting the correct timestep
	NVAR nsteps // number of time steps
	NVAR timestep // timestep in seconds
	variable/D timestep_new // for updating the timestep (v7.4)
  	NVAR KineticMassTransfer 	// turned into global variable on 11/12/13
	variable/D min_timestep_gas = 30 // seconds; the minimum timestep to use for gas-phase oxidation (typically, this is less restrictive than gas-particle partitioning)
	variable/D timestep_desired = 180 // seconds
	if(gasWLR > 0)	// deal with issues of numerical stability
		if(timestep > 0.05/gasWLR)
			timestep = 0.05/gasWLR
		endif
	endif
	if(timestep < timestep_desired)
		timestep_desired = timestep
		min_timestep_gas = timestep
	else
		timestep = timestep_desired
	endif
   	
 //2c.  //Create some parameters and waves FOR ZAVERI METHOD (first implementation, v7.4)
	NVAR alpha 	// accommodation coefficient, added 11/12/13
	make/d/o/n=(nC,nO) Diffusivity_matrix = Diff_CO2*(MW_CO2/MW_matrix) // m^2/s; assume that Diffusivity scales with MW and use CO2 as reference case
	make/o/d/n=(nC,nO) MeanFreePath_Matrix = 3*Diffusivity_Matrix/sqrt((8*1.381e-23*Tvar)/(pi*(MW_matrix/(1000*Na))))	// meters; D/(lambda*c_bar) = 1/3, for Fuchs expression
	//
	make/o/d/n=(nC,nO,nSizeBins) SatRatio_3D = 0 // The saturation ratio (Eqn. 18 in Zaveri)
	make/o/d/n=(nC,nO,nSizeBins) tmatrix_3D=0 // adaptive timestep (Eqn. 72/73 in Zaveri)
	make/o/d/n=(nC,nO,nSizeBins) Phi_GPP = 0 // relative driving force (Eqn. 70,71 in Zaveri)
	make/o/d/n=(nC,nO) tarray=0
	make/o/d/n=(nC,nO,nSizeBins) Knudsen_3D = nan, Fuchs_3D = nan
	make/o/d/n=(nC,nO,nSizeBins) deltaDynamicPartitioningPre = nan, deltaDynamicPartitioning=nan // matrices used in calculations
	make/o/d/n=(nSizeBins) Kelvin // Kelvin term
	variable/D deltat_min, deltat_min_new // the minimum timestep
	variable/D ntime_GPP, ntime_GPP_new // the rounded number of points that can be allowed for a GPP timestep

	NVAR MaxTime_hours // length of simulation, in hours // v7.4 moved from above allow for timestep to change, as necessary for dynamic partitioning
	nsteps = ceil(MaxTime_hours*60*60/timestep)//+1	// total number of time steps
	
// 2d. Parameters for dynamic partitioning, when using EULER method --> LEGACY
  	variable/D MassTransferMaxTime = 0.3	// seconds, for operator split, added 11/12/13
  		// specify time step for dynamic partitioning based on alpha
  	if(alpha < 1e-4)
  		MassTransferMaxTime = 30
  	elseif(alpha < 0.5e-3)
  		MassTransferMaxTime = 30
  	elseif(alpha <= 1e-2)
  		MassTransferMaxTime = 10
  	elseif(alpha <= 1e-1)
  		MassTransferMaxTime = 1
  	elseif(alpha <= 1)
  		MassTransferMaxTime = 0.1
  	endif
  	// try something with v7.4 to have the timestep vary depending on the alpha and the seed surface area
  	// set a maximum mass transfer timestep to 30 seconds; for Euler method
  	// allow for restarting of a run if this is not a good guess. 
  	MassTransferMaxTime = 5/alpha/(log(SeedSurfaceArea)-1)^2
	MassTransferMaxTime = min(MassTransferMaxTime,30)

   	variable/D MassTransferTimeScaling = round(TimeStep/MassTransferMaxTime)	// number of timesteps in operator split, added 11/12/13
   	MassTransferTimeScaling = MassTransferTimeScaling < 1 ? 1 : MassTransferTimeScaling

//**	// absorbing seed...NEED TO CHECK ON THIS WITH RESPECT TO THE ABOVE SEED STUFF
	NVAR AbsorbingSeed		// 0 = no, 1 = yes
	variable/D SeedMass = 0.0000000005 // ug/m^3
	if(AbsorbingSeed == 0)
		SeedMass = 1e-6		// ug/m^3
	else
		SeedMass = SeedVolConc * 1.0	// ug/m^3, and where 1.0 is the assumed density
	endif
	variable/D SeedMolecules = SeedMass * 1e-6 * 1e-6 * Na / SeedMW // molecules/cm^3

// For dilution
	NVAR DilutionVar //variable DilutionVar = 0//6.12 // % per hour (6.12% from Dzepina et al., 2011)
// If you want to adjust the precursor volatilties for some reason...
	NVAR logCstar_adjustmentfactor //variable logCstar_adjustmentfactor = 0.9 // multiplicative factor to adjust logCstar values. Important for isomer simulations

// Decide how to represent volatility...legacy and not currently in use (5/27/13)	
	string Cstar_evolution = "no" // yes to allow for time varying functional group evolution
	variable/D delta_logCstar_perO_final = 2.2 // relative to delta_logCstar_perO
	variable/D delta_logCstar_perO_perStep = (delta_logCstar_perO_final-delta_logCstar_perO)/nsteps
// Decide how to implement fragmentation	
	NVAR Pfrag_type // 0 = cfrag*Nox; 1 = (O:C)^mfrag (now default)
	NVAR Frag_Method = small_fragments // 0 for random number generator, 1 for equal probs (default)., 2 to only produce fragments with 1 carbon (HCHO, CO2, CH4)
  // some other silly things
	variable steady_state = 0, SSvalue = 0 // 0 for base case, any number for something else keep parent gas-phase abundance at initial abundance
// OH concentration...many of these are not in general use (5/27/13)
	NVAR OHconc // [OH] in molecules/cm^3
	NVAR OH_scale // scaling factor to have [OH] decrease exponentially with time; OHconc_t = OHconc_0*exp(-timewv[idex-1]*oh_scale)
	variable/D OH_profile = 2 // 1 to have a diurnal profile of OH concentrations; 2 for anything else
	variable/D OH_counter = 0 // related to diurnal profile
	variable/D day_counter = 0 // related to diurnal profile
	variable/D day_iseven = 1 // related to diurnal profile
	variable/D OH_max = 4e6 // related to diurnal profile
	variable/D OH_min = 1e5 // related to diurnal profile
// Parameters for heterogeneous chemistry
	variable/D OH_Flux = (OHconc*1e6)*1.381e-23*Tvar/sqrt(2*pi* (17/1000/6.022e23)*1.381e-23*Tvar) // molecules m^-2 s^-1
	variable/D OH_Rate // molecules/s; will be calculated from OH_Flux and particle surface area
	NVAR gammaOH // OH reactive uptake coefficient
	NVAR hetchem // 0 for no heterogeneous chemistry; 1 to include
// Parameters for sequential partitioning model (SPM)
	NVAR SPM // 0 to ignore (typical), 1 to include
// Parameters for printing amount of time taken
	variable/D t11, t22
// making the matrix smaller, in essence (not used)
	NVAR CutOffSmallStuff
// HOM and nucleation
	// Can either specify a time that nucleation "turns on", or can try and have nucleation
	// depend on the concentration of "HOMs" in the model; an unknown is what to clasify as HOM
	// Currently (v7.4.3) the empirical "on" time method works, while the [HOM] method is not operational
	NVAR NucleationTimeEmpirical // when nucleation should occur; h after start of experiment
	variable NucWait = 0 // will change to 1 once nucleation occurs
	variable/D HOM_logCstar = -2
	variable/D HOM, NucleationRate
	variable/D nuc_a1 = 0.04001 // parameters used later for nucleation rate calculation
	variable/D nuc_a2 = 1.848
	variable/D nuc_a5 = 0.1863
	make/o/d/n=(nSteps) HOM_time=0, NucleationRate_time = 0
	string HOM_note = "Concentration of HOM (logC* < " + num2str(HOM_logCstar) + "); molecules/cm3"
	note/k HOM_time HOM_note // "Concentration of HOM (logC* < 0.1); molecules/cm3"
	note/k NucleationRate_time "Nucleation Rate; p/(cm3.s); from Trostl et al., Nature, 2015"
	setscale/P x, 0, (timestep/60/60), "hours", HOM_time
	setscale/P x, 0, (timestep/60/60), "hours", NucleationRate_time
			
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
	make/o/d/n=(nC,nO) SatConc_Matrix = 0
	note/k SatConc_matrix "Saturation concentration, in molecules/m3"
	make/o/d/n=(nC,nO,nSizeBins) SatConc_Matrix_3D = 0
	make/o/d/n=(nC,nO) logCstar_matrix = 0
	//
	make/o/d/n=(nC,nO) SatConc_matrix_adj = 0
	note/k SatConc_matrix_adj "2D matrix with saturation concentrations, adjusted by the mole fraction of that species"
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
	make/o/d/n=(nC,nO,nSizeBins) PM3D_Matrix_Current=0
	make/o/d/n=(nC,nO) ParticleMolecules_matrix_cur=0
	make/o/d/n=(nC,nO) TotalMolecules_matrix_current=0
	make/o/d/n=(nC,nO) GasMolecules_matrix_current=0, GasMolecules_matrix_current2=0
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
	make/d/o/n=(nC,nO)/FREE C_Matrix_temp = 0
	make/d/o/n=(nC,nO)/FREE O_Matrix_temp = 0
	make/d/o/n=(nC,nO)/FREE H_Matrix_temp = 0
	//
//	make/d/o/n=(nC,nO)/FREE Diffusivity_matrix = 0
//	make/o/d/n=(nC,nO)/FREE MeanFreePath_Matrix = 0
//	make/o/d/n=(nC,nO)/FREE Knudsen_Matrix = 0
//	make/o/d/n=(nC,nO)/FREE Beta_Matrix = 0
	// v7.4 moving these above to use Zaveri method for dynamic partitioning
//	make/d/o/n=(nC,nO) Diffusivity_matrix = 0
//	make/o/d/n=(nC,nO) MeanFreePath_Matrix = 0
	make/o/d/n=(nC,nO) Knudsen_Matrix = 0
	make/o/d/n=(nC,nO) Beta_Matrix = 0
		//
	make/d/o/n=(nC,nO) RxnStep_Matrix_plus = 0
	make/d/o/n=(nC,nO) RxnStep_Matrix_plus1 = 0
	make/d/o/n=(nC,nO) RxnStep_Matrix_plus2 = 0
	make/d/o/n=(nC,nO) RxnStep_Matrix_plus3 = 0
	make/d/o/n=(nC,nO) RxnStep_Matrix_plus4 = 0
	make/d/o/n=(nC,nO) RxnStep_Matrix_minus = 0
	//
	make/o/d/n=(nC,nO)/FREE RxnStep_Matrix_Kinetic_On = 0
	make/o/d/n=(nC,nO)/FREE RxnStep_Matrix_Kinetic_Off = 0
	make/o/d/n=(nC,nO)/FREE RxnStep_Matrix_Kinetic_Pre1 = 0
	make/o/d/n=(nC,nO)/FREE RxnStep_Matrix_Kinetic_Pre2 = 0
	make/o/d/n=(nC,nO,nSteps) RxnStep_Matrix_Kinetic_Time = 0
	make/o/d/n=(nC*nO) RxnStep_array_minus
	make/d/o/n=(nC,nO) RxnStep_Matrix_frag = 0
	//
	make/d/o/n=(nC,nO) HetChem_Matrix_plus = 0
	make/d/o/n=(nC,nO) HetChem_Matrix_minus = 0
	//
	make/d/o/n=(nC,nO,nsteps) GasMolecules_Time = 0
	note GasMolecules_Time "gas-phase species concentration in molecules/cm^3"
	note GasMolecules_Time "Ncarbons = Rows; Noxygens = columns; Layers = time [hrs]"
	setscale/P z, 0, (timestep/60/60), GasMolecules_Time
	setscale/P x, nC, -1, "nCarbons", GasMolecules_Time
	setscale/P y, 0, 1, "nOxygens", GasMolecules_time
	//
	make/d/o/n=(nC,nO,nsteps) GasMass_Time = 0
	note GasMass_Time "gas-phase species concentration in molecules/cm^3"
	note GasMass_Time "Ncarbons = Rows; Noxygens = columns; Layers = time [hrs]"
	setscale/P z, 0, (timestep/60/60), GasMass_Time
	setscale/P x, nC, -1, "nCarbons", GasMass_Time
	setscale/P y, 0, 1, "nOxygens", GasMass_Time
	//
	make/d/o/n=(nC,nO,nsteps) ParticleMolecules_Time = 0
	note ParticleMolecules_Time "particle-phase species concentration in molecules/cm^3"
	note ParticleMolecules_Time "Ncarbons = Rows; Noxygens = columns; Layers = time [hrs]"
	setscale/P z, 0, (timestep/60/60), "hours", ParticleMolecules_Time
	setscale/P x, nC, -1, "nCarbons" ParticleMolecules_Time
	setscale/P y, 0, 1, "nOxygens" ParticleMolecules_Time
	//
	make/d/o/n=(nC,nO,nsteps) ParticleMass_Time = 0
	note ParticleMass_Time "particle-phase species concentration in molecules/cm^3"
	note ParticleMass_Time "Ncarbons = Rows; Noxygens = columns; Layers = time [hrs]"
	setscale/P z, 0, (timestep/60/60), "hours", ParticleMass_Time
	setscale/P x, nC, -1, "nCarbons" ParticleMass_Time
	setscale/P y, 0, 1, "nOxygens" ParticleMass_Time
	//
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
	//
	make/d/o/n=(nsteps) SeedConc_time = SeedVolConc
	setscale/P x, 0, (timestep/60/60), "hours", SeedConc_time
	// testing wall loss method
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
	make/o/d/n=(nC,nO) WallMass_matrix = 0 // absorbed wall mass; v7.4.3
	note WallMass_matrix "wall species concentration in ug/m3"
	make/d/o/n=(nC,nO) WallParticleMolec_Matrix = 0 // deposited particles; added 9/29/14
	note WallParticleMolec_Matrix "Material associated with particle wall deposition, in molecules/cm^3"
	note WallParticleMolec_Matrix "Ncarbons = Rows; Noxygens = columns"
	// 3D version to deal with size-dependent wall losses; added 9/29/14
	make/d/o/n=(nC,nO,nsizebins) WPM3D_Matrix = 0 // molecules/cm^3
	note WPM3D_Matrix "particle phase species concentration in wall-deposited particles in molecules/cm^3"
	note WPM3D_Matrix "Ncarbons = Rows; Noxygens = columns; nSizeBins = layers"

	make/d/o/n=(nsteps) TimeW = 0
	note TimeW "Reaction time [hrs]"
	make/d/o/n=(nsteps) O2C_time = nan
	note O2C_time "Oxygen-to-Carbon ratio of SOA"
	setscale/P x, 0, (timestep/60/60), "hours", O2C_time
	setscale/P d, 0, 0, "O:C", O2C_time
	make/d/o/n=(nsteps) H2C_time = nan
	note H2C_time "hydrogen-to-carbon ratio of SOA"
	setscale/P x, 0, (timestep/60/60), "hours", H2C_time
	setscale/P d, 0, 0, "H:C", H2C_time
	make/d/o/n=(nsteps) logCstar_mean_time = 0
	setscale/P x, 0, (timestep/60/60), "hours", logCstar_mean_time
	make/d/o/n=(nsteps) OCseed_time = 0
	//
	make/d/o/n=(nsteps) Coa_time = 0
	note/k Coa_time "SOA mass concentration [ug/m^3]"
	note Coa_time "VOC concentration = " + num2str(Ctot_ppm*1000) + " ppb"
	note Coa_time "Pfunc = " +num2str(ProbOx1)+";"+num2str(ProbOx2)+";"+num2str(ProbOx3)+";"+num2str(ProbOx4)
	note Coa_time "Pfrag = " + num2str(FragSlope) + " using type " + num2str(Pfrag_type) + " and method " + num2str(Frag_Method)
	note Coa_time "DLVP = " + num2str(delta_logCstar_perO)
	note Coa_time "kwall = " + num2str(gasWLR) + " s^-1"
	if(Krechmer_Cw==0)
		note Coa_time "Cw = " + num2str (Cwall) + " mg/m3"
	else
		note Coa_time "Cw from Krechmer (2016)"
	endif
	if(kineticmasstransfer == 0)
		note Coa_time "equilibrium partitioning"
	else
		note Coa_time "kinetic partitioning with alpha = " + num2str(alpha)
		note Coa_time "Np = " + num2str(Nparticles) + " p/cm3 & Dp,seed = " + num2str(SeedDiameter) + " nm"
	endif
	note Coa_time "Timestep = " + num2str(timestep) + " s"
	setscale/P x, 0, (timestep/60/60), "hours", Coa_time
	setscale/P d, 0, 0, "ug/m\S3\M", Coa_time
	//
	make/d/o/n=(nsteps) Coa_time_wallcorr = 0
	note Coa_time_wallcorr "Particle wall-loss corrected SOA mass concentration [ug/m^3]"
	note Coa_time_wallcorr "Coa_time_wallcorr = Coa_time + WallParticleMolec_Matrix"
	setscale/P x, 0, (timestep/60/60), "hours", Coa_time_wallcorr
	setscale/P d, 0, 0, "ug/m\S3\M", Coa_time_wallcorr
	//
	make/d/o/n=(nsteps) deltaHC_time = 0
	setscale/P x, 0, (timestep/60/60), "hours" deltaHC_time
	make/d/o/n=(nsteps) HC_ppb_time = 0
	note HC_ppb_time "Parent HC concentration [ppb]"
	setscale/P x, 0, (timestep/60/60), "hours" HC_ppb_time
	make/d/o/n=(nsteps) HC_ugm3_time = 0
	note HC_ppb_time "Parent HC concentration [ug/m3]"
	setscale/P x, 0, (timestep/60/60), "hours" HC_ugm3_time
	setscale/P d, 0, 0, "ppb" HC_ppb_time
	make/d/o/n=(nsteps) deltaHCfrac_time = 0
	setscale/P x, 0, (timestep/60/60), "hours", deltaHCfrac_time
	make/d/o/n=(nsteps) Yield_time = 0
	note Yield_time "Aerosol mass yield"
	setscale/P x, 0, (timestep/60/60), "hours", Yield_time
	make/d/o/n=(nsteps) Yield2_time = 0
	make/d/o/n=(nsteps) Nc_ave_time = 0
	make/d/o/n=(nsteps) No_ave_time = 0
	make/d/o/n=(nsteps) OH_wave = OHconc
	note OH_wave "[OH] in molecules/cm^3"
	setscale/P x, 0, (timestep/60/60), "hours", OH_wave
	setscale/P d, 0, 0, "molecules/cm\S3\M", OH_wave
	make/d/o/n=(nsteps) OH_exposure = 0
	make/d/o/n=(nsteps) Lifetimes = 0
	note Lifetimes "Number of oxidation lifetimes"
	setscale/P x, 0, (timestep/60/60), "hours", Lifetimes
	make/o/d/n=(nsteps) Molecules_time = 0
	make/o/d/n=(nC+1) Ncarbons_Wave = nC-x
	make/o/d/n=(nO+1) Noxygens_Wave = x
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
	make/o/d/n=(nsteps) Temperature_time = Temp_G // v7.4.1
	setscale/P x, 0, (timestep/60/60), "hours", Temperature_time
	note temperature_time "the temperature, in K"

 	make/o/d/n=(nC,nO,nSizeBins) RxnStepKinetic3D_Cond, RxnStepKinetic3D_Evap
 	make/o/d/n=(nC,nO,nSizeBins) RxnStepKinetic3D_Pre1, RxnStepKinetic3D_Pre2
 	make/o/d/n=(nC,nO,nSizeBins) MoleFraction3D
 	make/o/d/n=(nC,nO) MoleFraction2D
 	// v7.4.3 --> wall loss
 	make/o/d/n=(nC,nO) tmatrix_2D
 	make/o/d/n=(nC,nO) SatRatio_2D
 	make/o/d/n=(nC,nO) sumDeltaDP
 	make/o/d/n=(nC,nO) TotalWallMolecules_current
 	
 	//
 	make/o/d/n=(nSteps) TimeStep_wave = timestep
	
	//
	variable/D Ave_Nc
	variable/D Ave_No
	variable/D Ave_Nh
	variable/D molecules_PM_previous = 0
	
	TimeW = (x)*timestep/60/60 // hours

// 5c. Auto-Generate Dynamic Partitioning ODE procedure (added 08/26/14)
//	DP_Setup(Nspecies) // creates DynamicPartitioningFuncs.ipf --> a failed attempt at using the Igor solver

// 6. Populate initial logC* matrix
	logCstar_Matrix[][] = -0.0337*MW_Matrix[p] + 11.56 // estimated from Lide, CRC data for saturated hydrocarbons
	logCstar_Matrix[][] = logCstar_Matrix[p][0] - y*delta_logCstar_perO // adjust C* values based on number of oxygens per molecule
	logCstar_matrix += logCstar_adjustmentfactor // one can "adjust" the logCstar matrix from the base case assumptions; best to have this be 0
	duplicate/o logCstar_matrix logCstar_matrix_ref
//	logcstar_matrix[0][6] = 0.00001
	SatConc_Matrix = (10^(logCstar_Matrix)) * 1e-6 * Na / MW_Matrix // molecules/m^3; C* --> saturation concentration
// 6b. Enthalpy of vaporization; added 09/03/13
	SOM_deltaHvap(logCstar_Matrix=logCstar_matrix_ref,ConstantDHvap=0)	// ConstantDHvap = 0 means not constant, otherwise enter values in kJ/mol
	wave deltaHvap_matrix
	variable Tref = 298.15 // reference temperature; the above formula is for the reference temperatur0^
	// temp_g is the specified temperature
	logCstar_matrix = log((10^logCstar_matrix_ref)*(Tref/Temp_G)*exp(-(1000*deltaHvap_matrix/8.314)*(1/Temp_G - 1/Tref))) // v7.4.1
	SatConc_Matrix = (10^(logCstar_Matrix)) * 1e-6 * Na / MW_Matrix // molecules/m^3; C* --> saturation concentration
	
// 7. OH reaction rate coefficients
	// 09/01/13 --> added new method to estimate rate coefficients that is more continuous w.r.t. Nc and No.
	if(OligomerizationIsOn==0)	// no oligomerization
		SOM_RateCoefficients(C_matrix,O_matrix,0,O3_yn,TempK=Tvar,method=krxn_method,Ea=Ea_krxn)
		wave krxn_matrix
	// Adjust parent rate coefficient if desired
		if(krxn_parent == 0)
			// do nothing
		else
			krxn_matrix[0][0] = krxn_parent
		endif
	// Decide whether to use alternative scheme for isoprene where krxn(C5O1-C5O4) = krxn(isoprene)
		if(IsopreneIsSpecial)
			krxn_matrix[0][1,4] = krxn_parent
		endif
	else	// oligomerization
		SOM_RateCoefficients(C_matrix,O_matrix,0,O3_yn,TempK=Tvar,method=krxn_method,Ea=Ea_krxn)
		wave krxn_matrix
		krxn_matrix[nCarbons][Nox_precursor] = krxn_parent
	endif
	setscale/P x, nC, -1, "nCarbons" krxn_matrix
	setscale/P y, 0, 1, "nOxygens" krxn_matrix
	
	// Uncomment to have 1st generation products only
	variable/D krxn_00 = krxn_matrix[0][0]
//	krxn_matrix=0
//	krxn_matrix[0][0] = krxn_00
	
	if(OligomerizationIsOn==1)
		SOM_OligomerRateCoef(nC,nO,krxn_base=krxn_base_olig)
		wave krxn_matrix_olig
	endif
	
// 8. Populate fragmentation probability matrix, which is a 3D matrix with rows = #carbons, columns = #oxygens, layers = probability
	SOM_Fragmentation(Frag_Method,FragSlope,Pfrag_type,C_matrix,O_matrix,OC_matrix)
	wave Prob1_Matrix
	wave Prob2_Matrix
	wave Prob_Matrix
	wave Frag1_Matrix
	wave Frag1_Array
	
	if(Nox_precursor > 0)
//		Frag1_matrix[0][Nox_precursor] = 0.7
//		Frag1_array[Nox_precursor] = 0.7
	endif

// 9. Populate initial mass/molecules waves	

// 1. Decide here how to input initial concentration	
	// can specify multiple compounds to hold initial concentrations
	variable multiple_compounds = 0
	make/o/d/n=(nC) Ctot_init_wave = 0	 // an array with the initial concentrations of all species
	setscale/P x, nC, -1, Ctot_init_wave
	variable/D MWstart
	if(OligomerizationIsOn==0) // no oligomerization
		MWstart = 12*nC + (nC-2)*2 + 2*3 - Hadjustment - Nox_precursor*H_per_O + Nox_precursor*16
		Ctot_init = Ctot_ppm*(NumberDensityAir*1e-6/Na)*(MWstart*1e6) // ppm to ug/m^3
		Ctot_init_wave[0] = Ctot_init
		if(multiple_compounds==0) // just one VOC
			Ctot_init_wave[0] = Ctot_init
		else // multiple compounds
			// totally empirical and up to user to define this; SOM is not really set up for this (yet)
//			Ctot_init_wave = Ctot_ppm*10^(-1*(nC-p)/10)//Ctot_ppm*(p+1)^4 // assume a particular relationship between carbon number and abundance; made up for now
			wave IVOC_conc_Zhao
			Ctot_init_wave[0,9] = IVOC_conc_zhao
			wavestats/q Ctot_init_wave // get sum to normalize
			Ctot_ppm = V_sum
//			Ctot_init_wave = (Ctot_init_wave/V_sum)*Ctot_ppm
			Ctot_init_wave *= (NumberDensityAir*1e-6/Na)*(MW_matrix[p][0]*1e6) // ppm to ug/m^3
			wavestats/q Ctot_init_wave
			Ctot_init = V_sum // total concentration of gas-phase precursors
		endif
	else
		MWstart = 12*nCarbons + (nCarbons-2)*2 + 2*3 - Hadjustment - Nox_precursor*H_per_O + Nox_precursor*16
		Ctot_init = Ctot_ppm*(NumberDensityAir*1e-6/Na)*(MWstart*1e6) // ppm to ug/m^3
		Ctot_init_wave[nCarbons] = Ctot_init
	endif

	if(quiet==0)
		print "Ctot = " + num2str(ctot_init) + " in ug/m^3 for Nc = " + num2str(nC) + " with FragSlope = " + num2str(FragSlope) + " and dlVP = " + num2str(delta_logCstar_perO)+ " and ProbOx = " + num2str(ProbOx1) + ";"+ num2str(ProbOx2) + ";"+ num2str(ProbOx3) + ";"+ num2str(ProbOx4) + ";"
	endif
// End 1. Decision Finish

	GasMass_matrix[][Nox_precursor] = Ctot_init_wave[p] // Set initial gas-phase mass concentration as total 
	
	if(cheat==1) // for testing purposes
		GasMass_Matrix[][] = Ctot_init_wave[0]/numpnts(GasMass_matrix)
	endif
	
	GasMass_Matrix = (numtype(H_Matrix) == 2 ? NaN : GasMass_Matrix) // take care of pesky NaN values
	GasMolecules_Matrix = GasMass_Matrix * 1e-6 * 1e-6 * Na / MW_Matrix // convert ug/m^3 to molecules/cm^3
	TotalMolecules_Matrix = GasMolecules_Matrix // (molecules/cm^3)
	variable noEqm = 0	
   // Populate initial Gas + particle = total waves/matrices	
	if(noEqm==0)
		// set initial condition as equilibrium --> most important when dealing with precursors having low volatility
		EqmCalc(GasMolecules_Matrix,SatConc_Matrix,MW_Matrix,SeedConc=SeedMolecules,SeedMW=SeedMW,units=units)
//		EqmCalc(GasMolecules_Matrix,SatConc_Matrix,MW_Matrix,Coa_guess = V_Sum,SeedConc=SeedMolecules,SeedMW=SeedMW,units=units)
		wave Coa_eqm = Coa	// Does not include seed particle mass, i.e. is just the mass concentration of the "SOA"
		GasMolecules_matrix = TotalMolecules_Matrix - Coa_eqm
		ParticleMolecules_matrix = Coa_eqm
		Coa_eqm *= MW_matrix[p][q]*1e12/6.022e23
		wavestats/q Coa_eqm
		Coa_time[0] = V_sum
	else
		// no eqm to start and everything is left in the gas phase
	endif
	
	if(DynamicPartitioningMethod==1) // Zaveri method, v7.4
		// do not allow things to initially equilibrate; they will quickly if alpha is large
		GasMolecules_matrix = TotalMolecules_matrix
		ParticleMolecules_matrix = 0
		Coa_time[0] = 0
	endif
	
	PM3D_matrix = ParticleMolecules_matrix[p][q] * dVdlogDp_norm[r] // /nSizeBins	// divide mass among all size bins to start; Changed for v7.4 
	TotalMolecules_Matrix = GasMolecules_matrix + ParticleMolecules_Matrix
	if(DynamicPartitioningMethod==0 || KineticMassTransfer==0)
		wavestats/q ParticleMolecules_matrix
		ParticleMoleFraction_matrix = ParticleMolecules_Matrix/V_sum
		ParticleMoleFraction_Time[][][m] = ParticleMoleFraction_Matrix[p][q]
		PMF3D_matrix[][][] = ParticleMoleFraction_matrix[p][q]
	else
		ParticleMoleFraction_matrix = 0
		ParticleMoleFraction_time[][][0] = 0
		PMF3D_matrix[][][] = 0
	endif
   // initial particle diameter
   	Dp_time[0] = DpStart
   	
// 10. Initialize Kinetics (populate time-dependent waves that store values based on the initial values)
	GasMolecules_Time[][][0] = GasMolecules_matrix[p][q]
	ParticleMolecules_Time[][][0] = ParticleMolecules_Matrix[p][q]
	TotalMolecules_Time[][][0] = TotalMolecules_Matrix[p][q]
	
	if(steady_state!=0)
		Ctot_init = GasMolecules_Matrix[0][0]
		Ctot_init_wave = GasMolecules_matrix[p][0]
	endif
	
	make/o/d/n=(nC,nO) TempWave
	variable/D TimeVar
	variable/D FragVar
	variable/D StepVar

// 10b. Decide whether to utilize OH wave that has been fit to VOC decay (using linear interpolations)
	// WILL NEED TO REVISIT WITH ADAPTIVE TIMESTEPPING V7.4
//	SetOH_from_experiment(krxn_matrix[0][0])
	wave OH_exp	// molecules/cm^3
	wave OH_exp_time	// hours
	if(UseScalingForOxidant==2)
		Interpolate2/T=1/N=200/I=3/Y=OH_exp_Interp/X=TimeW OH_exp_time, OH_exp // interpolate OH_exp from OH_exp_time to TimeW
		wave OH_exp_Interp // should have same number of points as TimeW
		wavestats/q OH_exp
		OH_exp_Interp = OH_exp_Interp < V_min ? V_min : OH_exp_Interp // do not let interpolated values be < minimum observed value (typically occurs if you interpolate outside of an appropriate range)
	endif
	// if O3...
	if(O3_yn==1)
		if(UseScalingForOxidant==2)
			OH_wave = OH_exp_interp*1e-9 * (1e-6*NumberDensityAir)
		else
			OH_wave = O3_conc * 1e-9 * (1e-6*NumberDensityAir)
		endif	
	endif
	
	variable/D initVOC = GasMolecules_Time[0][0][0]
	
// 10c. Gas-phase wall loss initialization
	SOM_GasPhaseWallLoss(logCstar_matrix,O_matrix,gasWLR,Cw_base=Cwall,kwg_increase_perO=kwall_gas_scaling,Krechmer_Cw=Krechmer_Cw)
	wave kwg_off		// 2D matrix with composition dependent rate of desorption (1/s) from walls
	wave Cw_wave // effective wall concentration (ug/m3); v7.4.3
	variable VWL_method = 1 // 0 = Euler, 1 = based on Zaveri (MOSAIC); v7.4.3

// 10e. Some timing stuff		
	t11 = ticks
	variable FirstTime = 1
	variable counter = 0
	variable count = 0

	timestep_new = timestep
	// Store initial values in case there is a need to restart (this sometimes is necessary)
	duplicate/o GasMolecules_matrix, GasMolecules_matrix_init
	duplicate/o TotalMolecules_matrix, TotalMolecules_matrix_init
	duplicate/o ParticleMolecules_matrix, PM_matrix_init
	duplicate/o PM3D_matrix, PM3D_Matrix_init
	
//	wave kOH_special	// testing!!!
	timestep_wave[0] = deltat_min

// 11. START KINETICS
	m = 1 // set the first index value
	do // v7.4: change to do loop
		timestep = timestep_new
		timew[m] = timew[m-1] + timestep/3600

//		timestep = timestep_new
//		deltat_min = deltat_min_new
//		ntime_GPP = ntime_GPP_new
		
//	gasmolecules_matrix[0][0] = initVOC
//		krxn_matrix = kOH_special[p][q][m-1]
// 11a: simulates continuous input of parent molecules (i.e. emissions) to keep constant total gas-phase amount.
//		if(steady_state!=0) 
//			if(Multiple_compounds==0)
//				GasMolecules_Matrix[0][0] += RxnStep_matrix_minus[0][0]// (Ctot_init_wave[0]*steady_state)
//			else
//				GasMolecules_Matrix[][0] += RxnStep_matrix_minus[p][0]//(Ctot_init_wave*steady_state)
//			endif
//		endif
// 11b: Get current OH or O3 concentration (since it may not be constant)	

		if(OH_profile == 1) // in case you want to try and use a diurnal profile
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
		elseif(UseScalingForOxidant == 1)
			// Allow the OH concentration to vary with time; exponential function
			OH_wave[m-1] = SOM_SetOHconc(timew,OHconc,OH_scale,m)
		elseif(UseScalingForOxidant == 2 && O3_yn!=1)
			// use the interpolated OH concentration, based on the observations
			OH_wave[m-1] = OH_exp_Interp[m-1] //OH_exp[m-1]
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
		
		if(cheat==1)
			OH_wave[m-1] = 0
			krxn_matrix = 0
		endif
		
		// update the logCstar matrix based on current temperature;
		// if you want a time-varying temperature, need to edit
		logCstar_matrix = log((10^logCstar_matrix_ref)*(Tref/Temperature_time[m-1])*exp(-(1000*deltaHvap_matrix/8.314)*(1/Temperature_time[m-1] - 1/Tref)))
		SatConc_Matrix = (10^(logCstar_Matrix)) * 1e-6 * Na / MW_Matrix // molecules/m^3; C* --> saturation concentration

// 11bc: this can be used to test "reactive uptake", ala Shiraiwa et al., 2013 PNAS	
//		if(m >= nsteps-nsteps/2 && counter==0)
//			gasmolecules_matrix[Nc/2][1] += 5e11
//			OH_wave[m-1] = 0	
//			counter = 1
//		elseif(m > nsteps/3)
//			OH_wave[m-1] = 0				
//		endif

// 11c: Run gas phase reaction to determine delta[HC]
		// Chemical LOSS
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

		// FORMATION - no fragmentation yet
		variable FragMethTest = 1
		// this is the correct implementation of fragmentation. changed 07/11/13
		if(FragSlope!=0 || O3_yn==0)	// if fragmentation is on and this is a reaction with OH
			RxnStep_Matrix_plus1=0; RxnStep_Matrix_plus2=0; RxnStep_Matrix_plus3=0; RxnStep_Matrix_plus4=0
			RxnStep_Matrix_plus1[][1,] = RxnStep_Matrix_minus[p][q-1]*ProbOx1*(1-Frag1_Matrix[p][q-1])
			RxnStep_Matrix_plus2[][2,] = RxnStep_Matrix_minus[p][q-2]*ProbOx2*(1-Frag1_Matrix[p][q-2])
			RxnStep_Matrix_plus2[][nO-1] += (RxnStep_Matrix_minus[p][nO-2])*ProbOx2*(1-Frag1_Matrix[p][nO-2])
			RxnStep_Matrix_plus3[][3,] = RxnStep_Matrix_minus[p][q-3]*ProbOx3*(1-Frag1_Matrix[p][q-3])
			RxnStep_Matrix_plus3[][nO-1] += (RxnStep_Matrix_minus[p][nO-2]*(1-Frag1_matrix[p][nO-2])+RxnStep_Matrix_minus[p][nO-3]*(1-Frag1_matrix[p][nO-3]))*ProbOx3
			RxnStep_Matrix_plus4[][4,] = RxnStep_Matrix_minus[p][q-4]*ProbOx4*(1-Frag1_Matrix[p][q-4])
			RxnStep_Matrix_Plus4[][nO-1] += (RxnStep_Matrix_minus[p][nO-2]*(1-Frag1_matrix[p][nO-2])+RxnStep_Matrix_minus[p][nO-3]*(1-Frag1_matrix[p][nO-3])+RxnStep_Matrix_minus[p][nO-4]*(1-Frag1_matrix[p][nO-4]))*ProbOx4
	 		RxnStep_Matrix_plus[][] = RxnStep_Matrix_plus1 + RxnStep_matrix_plus2 + RxnStep_Matrix_Plus3 + RxnStep_matrix_Plus4
	 	else	// for O3 reactions, added 08/29/13
			RxnStep_Matrix_plus = 0;RxnStep_Matrix_plus1=0; RxnStep_Matrix_plus2=0; RxnStep_Matrix_plus3=0; RxnStep_Matrix_plus4=0
			RxnStep_Matrix_plus1[0][1,] = RxnStep_Matrix_minus[0][q-1]*ProbOx1
			RxnStep_Matrix_plus2[0][2,] = RxnStep_Matrix_minus[0][q-2]*ProbOx2
			RxnStep_Matrix_plus2[0][nO-1] += (RxnStep_Matrix_minus[0][nO-2])*ProbOx2
			RxnStep_Matrix_plus3[0][3,] = RxnStep_Matrix_minus[0][q-3]*ProbOx3
			RxnStep_Matrix_plus3[0][nO-1] += (RxnStep_Matrix_minus[0][nO-2]+RxnStep_Matrix_minus[0][nO-3])*ProbOx3
			RxnStep_Matrix_plus4[0][4,] = RxnStep_Matrix_minus[0][q-4]*ProbOx4
			RxnStep_Matrix_Plus4[0][nO-1] += (RxnStep_Matrix_minus[p][nO-2]+RxnStep_Matrix_minus[0][nO-3]+RxnStep_Matrix_minus[0][nO-4])*ProbOx4
	 		RxnStep_Matrix_plus[0][] = RxnStep_Matrix_plus1[0][q] + RxnStep_matrix_plus2[0][q] + RxnStep_Matrix_Plus3[0][q] + RxnStep_matrix_Plus4[0][q]
		endif

// 11d: Determine what molecules are formed as a result of fragmentation
		if(FragSlope != 0) // only run this if fragmentation indeed occurs.
			RxnStep_Matrix_Frag = 0			
			for(i=0;i<nC;i+=1)
				RxnStep_array_Minus[i*nO,(i+1)*nO-1] = RxnStep_Matrix_Minus[i][p-i*nO]
			endfor
			multithread prob_matrix = (Prob1_Matrix+Prob2_Matrix)*Frag1_Array[r]*RxnStep_Array_Minus[r]
			MatrixOp/O RxnStep_Matrix_Frag = sumBeams(prob_matrix)
			GasMolecules_Matrix[][] += RxnStep_Matrix_Frag + RxnStep_Matrix_Plus // New version as of 07/11/13
		else // No fragmentation
			GasMolecules_matrix[][] += RxnStep_Matrix_plus // molecules formed per step, accounting for different numbers of oxygens added per step
		endif

// 11e: Wall loss of gas-phase species; added 7/14/2013, fixed on 8/20/2013
		if(gasWLR != 0)
			if(VWL_method==0)
			// the Euler method was used exclusively prior to v7.3.4;  step splitting added at v7.4
				variable nsteps_WLR // # of steps to take; v7.4
				nsteps_WLR = ceil((timestep/30)*(gasWLR/1e-4)) // must be integer; empirical expression here
				nsteps_WLR = max(1,nsteps_WLR) // must have at least one step
				if(m==1)
					print nsteps_WLR // print for diagnostic
				endif
				if(gasWLmethod == 0) // irreversible uptake
					GasMolecules_matrix[][] -=  GasMolecules_matrix[p][q] * gasWLR * timestep 
				else // reversible partitioning
					for(k=0;k<nsteps_WLR;k+=1) // split this into multiple steps in v7.4, just in case (could still probably be more robust)
						// wall deposition
						RxnStep_Matrix_minus[][] = GasMolecules_matrix[p][q] * gasWLR * timestep/nsteps_WLR	// gas-phase wall loss
						RxnStep_Matrix_minus = RxnStep_Matrix_minus > GasMolecules_matrix ? GasMolecules_matrix : RxnStep_Matrix_minus
						GasMolecules_matrix[][] -= RxnStep_Matrix_minus
						WallMolecules_Matrix[][] += RxnStep_Matrix_minus
						// wall desorption
						kwg_off = kwg_off*timestep > 1 ? 1/timestep : kwg_off
						RxnStep_Matrix_plus[][] = WallMolecules_matrix[p][q] * kwg_off * timestep/nsteps_WLR	// desorption from walls
						GasMolecules_matrix[][] += RxnStep_matrix_Plus[p][q] 
						WallMolecules_matrix[][] -= RxnStep_Matrix_Plus[p][q]
					endfor
					WallMolecules_Time[][][m] = WallMolecules_Matrix[p][q]
				endif
			elseif(VWL_method==1) // based on Zaveri (MOSAIC) method, assuming walls can be treated as "liquid"; added v7.4.3
		 		// First, get current total concentrations of each species
	 			GasMolecules_matrix_current = GasMolecules_matrix // gas concentration (molecules/cm3); this is from the current time step, after adjusting for gas-phase reactions
	 			GasMolecules_matrix_current2 = GasMolecules_matrix // repeat, to use later
	 			TotalWallMolecules_current = GasMolecules_matrix + WallMolecules_matrix // (molecules/cm3)
	 			WallMass_Matrix = WallMolecules_matrix * (MW_matrix/Na) * 1e12 // ug/m3 = (molecules/cm3)*(1e6 cm3/m3)*(mol/molec)*(g/mol)*(1e6 ug/g)
				//
			 // Next, get current mole fraction (really, mass fraction)
	 			MoleFraction2D = WallMass_matrix / (WallMass_matrix+Cw_wave) // divide number of molecules (really, mass) for each species by the effective wall concentration for each species
//		 	// Get saturation concentation based on previous time step
				SatConc_matrix_adj = MoleFraction2D*SatConc_matrix*1e-6 // Zaveri (2008), Eqn. 71, C*,t; units: molecules/cm3 = unitless * molecules/m3
				SatConc_matrix_adj = numtype(SatConc_matrix_adj)==2 ? 0 : SatConc_matrix_adj
				phi_GPP = 0
				phi_GPP[][][0] = GasMolecules_matrix_current[p][q] > SatConc_matrix_adj[p][q] ? GasMolecules_matrix_current[p][q] : SatConc_Matrix_adj[p][q] // determine which is larger, gas phase concentration or saturation concentration (phi_GPP is 3D, but this is okay)
				phi_GPP[][][0] =(GasMolecules_matrix_current[p][q]-SatConc_Matrix_adj[p][q])/phi_GPP[p][q][0]
				tmatrix_2D = gasWLR*abs(phi_GPP[p][q][0])
				tmatrix_2D = tmatrix_2D==0 ? nan : 1/tmatrix_2D // assume that alpha for uptake onto walls is unity (or change this if you want)
				wavestats/q tmatrix_2D
				deltat_min = V_min // choose the minimum time step as the answer
				MassTransferTimeScaling = max(1,ceil(timestep/deltat_min)) // lower limit of unity for number of steps for GPP; adopted from simpleSOM
				ntime_GPP = MassTransferTimeScaling // number of steps; this is just using different variable names
			
			// now that you have the appropriate timestep, calculate the mass transfer using MOSAIC method
				for(k=0;k<ntime_GPP;k+=1)
					// Get current mass fraction
					WallMass_Matrix = WallMolecules_matrix * (MW_matrix/Na) * 1e12 // ug/m3 = (molecules/cm3)*(1e6 cm3/m3)*(mol/molec)*(g/mol)*(1e6 ug/g)
					MoleFraction2D = WallMass_Matrix / (WallMass_Matrix+Cw_wave) // fraction of wall that is a given species
			 		// Get current saturation concentration and saturation ratio, based on mole fraction
					SatConc_matrix_adj = MoleFraction2D*SatConc_matrix*1e-6 // Zaveri (2008), Eqn. 71, C*,t; units: molecules/cm3 = unitless * molecules/m3
					SatConc_matrix_adj = numtype(SatConc_matrix_adj)==2 ? 0 : SatConc_matrix_adj
					SatRatio_2D = WallMolecules_matrix==0 ? 1 : (SatConc_matrix_adj* (MW_matrix/Na) * 1e12) / (WallMass_Matrix + Cw_wave) // Saturation ratio; SatConc / Wall Stuff; unitless: molecules/cm3 / molecules/cm3
					// calculate gas-phase concentration at t+1 (Eqn. 23 from Zaveri): implemented in steps
					deltaDynamicPartitioning[][][0] = 1*(WallMolecules_matrix[p][q]) / (1+(timestep/ntime_GPP)*gasWLR*SatRatio_2D[p][q]) // first part of sum in numerator
						// note: for the above equation, Zaveri has this with a negative sign. But it seems as if it should be positive. 
					sumDeltaDP = deltaDynamicPartitioning[p][q][0] // leftover from copying from particles
					GasMolecules_matrix_current2[][] = TotalWallMolecules_current - sumDeltaDP // numerator in Eqn. 23
					deltaDynamicPartitioning[][][0] = (gasWLR)/(1+(timestep/ntime_GPP)*gasWLR*SatRatio_2D[p][q]) //  denominator 2nd term (the term to be summed) in Eqn. 23
					sumDeltaDP = deltaDynamicPartitioning[p][q][0] // sum in denominator from Eqn. 23; sum over all sizes
					GasMolecules_matrix_current2 /= (1+(timestep/ntime_GPP)*sumDeltaDP) // this is the gas concentration at timestep t+1
					// calculate particle-phase concentration at t+1 (Eqn. 19 from Zaveri). Uses just-calculated gas-phase concentration
					WallMolecules_matrix = (WallMolecules_matrix[p][q] + (timestep/ntime_GPP)*gasWLR*GasMolecules_matrix_current2[p][q]) // numerator
					WallMolecules_matrix /= (1 + (timestep/ntime_GPP)*gasWLR*SatRatio_2D[p][q]) // Particle concentration at t+1; divide numerator by denominator
					// reset values for next iteration
					GasMolecules_matrix = GasMolecules_matrix_current2 // gas-phase concentrations
					TotalWallMolecules_current = WallMolecules_Matrix + GasMolecules_Matrix // recalculate, just to be sure			
//					endif
				endfor
			endif
		endif
		
// 11ee: Particle loss to the walls; added 08/20/2013; updated 09/29/14
		if(ParticleWallLoss!=0 && timew[m] >= NucleationTimeEmpirical) // added nucleation time conditional in v7.4.3
		// v7.4 NEED TO MAKE SURE THAT THIS ACCOUNTS FOR LOSS OF SEED MATERIAL TOO
			Diameter = Diameter_time[m-1][p] // extract array from the time matrix
			particleWLR_wave[] = SOM_ParticleWallLossRate(Diameter,particlewallloss)	// get size dependent wall-loss rate (1/s)
			for(index_size=0;index_size<nSizeBins;index_size+=1)
				particleWLR_wave[index_size] = SOM_ParticleWallLossRate(Diameter_time[m-1][index_size],particlewallloss)	// get size dependent wall-loss rate (1/s)
			endfor

			NpPerBin_time[m][] = NpPerBin_time[m-1][q] - NpPerBin_time[m-1][q]*particleWLR_wave[q]*timestep				// loss of particles
			WPM3D_Matrix[][][] += ((NpPerBin_time[m-1][r]-NpPerBin_time[m][r])/NpPerBin_time[m-1][r])*PM3D_Matrix[p][q][r]	// wall-bound particles, binned
			PM3D_Matrix[][][] -= ((NpPerBin_time[m-1][r]-NpPerBin_time[m][r])/NpPerBin_time[m-1][r])*PM3D_Matrix[p][q][r]		// suspended particles, binned
			MatrixOp/O WallParticleMolec_Matrix = sumBeams(WPM3D_Matrix)	// total wall-bound material
			MatrixOp/O ParticleMolecules_Matrix = sumBeams(PM3D_Matrix)		// total suspended material
		endif
		
		if(CutOffSmallStuff!=0)
			// cut off small things
			GasMolecules_matrix[nC-CutOffSmallStuff,][] = 0
			ParticleMolecules_matrix[nC-CutOffSmallStuff,][] = 0
		endif		
		// calculate total species concentration after (1) gas-phase reaction, (2) gas-phase wall loss and (3) particle-phase wall loss
		TotalMolecules_Matrix = GasMolecules_matrix + ParticleMolecules_Matrix // Gas is new; Particle is from previous step, adjusted for wall losses
				
// 11f: Mass Transfer to Particle Phase
		if(KineticMassTransfer==0)
			// Eqm Partitioning Calculation
			timestep_new = timestep
			if(SPM == 0) // as normal
				if(hetchem==1 && O3_yn!=1) // heterogeneous chemistry & OH
					SOM_Heterogeneous(m,ParticleMolecules_Matrix,Np,Tvar,gammaOH,timestep,FragSlope,ProbOx1,ProbOx2,ProbOx3,ProbOx4,SeedMolecules)
				endif
				wavestats/q ParticleMolecules_Matrix
				EqmCalc(TotalMolecules_Matrix,SatConc_Matrix,MW_Matrix,Coa_guess = V_Sum,SeedConc=SeedMolecules,SeedMW=SeedMW,units=units,Tvar=Tvar,deltaHvap_matrix=deltaHvap_matrix)
				wave Coa_eqm = Coa // does not include seed particle mass
				if(NoSOA==0)	// only operate if you want to include particle-gas partitioning...which is the defaul case
					GasMolecules_matrix = TotalMolecules_Matrix - Coa_eqm
					ParticleMolecules_matrix = Coa_eqm
				endif
			else // use Sequential Parititioning Model: See Cappa and Wilson, ACP, 2011
				//EqmCalc_Seed_SOM(GasMolecules_Matrix,SatConc_Matrix,MW_matrix,SeedMolecules,molecules_PM_previous) // Requires as input molecules/cm^3; molecules/m^3; g/mol for units of input waves
				EqmCalc(GasMolecules_Matrix,SatConc_Matrix,MW_Matrix,Coa_guess = molecules_PM_previous,SeedConc=SeedMolecules,SeedMW=SeedMW,units=units)
				wave Coa
				if(hetchem==1 && O3_yn!=1)
					SOM_Heterogeneous(m,ParticleMolecules_Matrix,Np,Tvar,gammaOH,timestep,FragSlope,ProbOx1,ProbOx2,ProbOx3,ProbOx4,SeedMolecules)
					wavestats/q Coa
					EqmCalc(GasMolecules_Matrix,SatConc_Matrix,MW_Matrix,Coa_guess = V_Sum,SeedConc=SeedMolecules,SeedMW=SeedMW,units=units)
				endif
				wave Coa_eqm = Coa // does not include seed particle mass
				wavestats/q Coa_eqm
				if(stringmatch(units,"mass"))
					molecules_PM_previous = V_sum * MW_matrix[0][0]*1e12/6.022e23
				else
					molecules_PM_previous = V_sum
				endif
				GasMolecules_matrix -= Coa_eqm // Total Gas from previous step - newly formed aerosol
				ParticleMolecules_matrix += Coa_eqm // add particle mass and then assume no longer partitions
			endif
		else // Dynamic partitioning, assuming absorptive partitioning
			if(DynamicPartitioningMethod==0) // Euler method; Now Legacy for SOM_V1 (eventually delete this)
			timestep_new = timestep
			// If first step, allow equilibration based on initial gas-phase reaction to deal with zero mass on particles (v7.4)
				if(m==1) // first step
					TotalMolecules_Matrix = TotalMolecules_matrix==0 ? 1 : TotalMolecules_matrix
					EqmCalc(TotalMolecules_Matrix,SatConc_Matrix,MW_Matrix,Coa_guess = V_Sum,SeedConc=SeedMolecules,SeedMW=SeedMW,units=units,Tvar=Tvar,deltaHvap_matrix=deltaHvap_matrix)
					GasMolecules_matrix = TotalMolecules_Matrix - Coa_eqm
					ParticleMolecules_matrix = Coa_eqm
					PM3D_Matrix = ParticleMolecules_matrix[p][q] * dVdlogDp_norm[r] // distribute across size bins
				endif
			// loop over all particle sizes
	 			for(index_size=0;index_size<nSizeBins;index_size+=1)
	 				Knudsen_Matrix = 2*MeanFreePath_Matrix/(Diameter_time[m-1][index_size]*1e-9)	// Fuchs, assuming Delta=lambda; see S&P p. 602
					Beta_Matrix = (0.75*alpha*(1+Knudsen_Matrix))/(Knudsen_Matrix^2+Knudsen_Matrix+0.283*Knudsen_Matrix*alpha+0.75*alpha) // Fuchs-Sutugin
					if(timew[m-1] < NucleationTimeEmpirical) // Assume zero particles, no condensation surface area, for nucleation
						RxnStepKinetic3D_Pre1[][][index_size] = 0
					else // nucleation! or seed aerosol present for condensation
						RxnStepKinetic3D_Pre1[][][index_size] = 4*pi*(Diameter_time[m-1][index_size]*1e-9/2)*Diffusivity_Matrix[p][q]*(timestep/MassTransferTimeScaling)*NpPerBin_time[m][index_size]*Beta_Matrix[p][q]*1e6
					endif
				endfor
					RxnStepKinetic3D_Pre2[][][] = RxnStepKinetic3D_Pre1[p][q][r]*SatConc_Matrix[p][q]*1e-6
	
		 		// Now, transfer some mass	 			
				variable testGPP=0

				for(k=0;k<MassTransferTimeScaling;k+=1)
				if(testGPP==1) // this is here for testing. default location is below
	 				MoleFraction3D[][][] = PM3D_Matrix[p][q][r]*MW_matrix[p][q] // weight by molecular weight to use mass fraction
		 				for(index_size=0;index_size<nSizeBins;index_size+=1) // need to normalize each size bin; so loop through bins
		 					ParticleMass_Matrix = MoleFraction3D[p][q][index_size]
		 					wavestats/q ParticleMass_matrix
		 					if(V_sum==0)
								MoleFraction3D[][][index_size] = numtype(ParticleMass_Matrix[p][q])==2 ? nan : 1
							else 
								if(absorbingseed==0)
									MoleFraction3D[][][index_size] = ParticleMass_Matrix[p][q]/(V_sum)
								else
									MoleFraction3D[][][index_size] = ParticleMass_Matrix[p][q]/(V_sum+MoleculesSeedPerBin[index_size]*seedMW)
								endif
							endif
						endfor
					endif
	if(k==1)
//	abort
	endif					// condensation
	 				RxnStepKinetic3D_Cond[][][] = RxnStepKinetic3D_Pre1[p][q][r]*GasMolecules_Matrix[p][q] // determine condensation amount; molecules/bin (delta t already accounted for)
	 				MatrixOp/O SumKinetics = sumBeams(RxnStepKinetic3D_Cond) // sum over particle size bins
	 				RxnStepKinetic3D_Cond[][][] = SumKinetics[p][q] > GasMolecules_Matrix[p][q] ? 0 : RxnStepKinetic3D_Cond[p][q][r] // conditional to set things to zero if too much mass is transferred
	 				PM3D_Matrix[][][] += RxnStepKinetic3D_Cond[p][q][r] // add to particle phase across size bins; molecules/bin
	 				GasMolecules_matrix -= SumKinetics[p][q]
					if(testGPP==0) // this is the default, (pre v7.4)
		 				MoleFraction3D[][][] = PM3D_Matrix[p][q][r]*MW_matrix[p][q] // weight by molecular weight to use mass fraction
		 				for(index_size=0;index_size<nSizeBins;index_size+=1) // need to normalize each size bin; so loop through bins
		 					ParticleMass_Matrix = MoleFraction3D[p][q][index_size]
		 					wavestats/q ParticleMass_matrix
		 					if(V_sum==0)
								MoleFraction3D[][][index_size] = numtype(ParticleMass_Matrix[p][q])==2 ? nan : 1
							else 
								if(absorbingseed==0)
									MoleFraction3D[][][index_size] = ParticleMass_Matrix[p][q]/(V_sum)
								else
									MoleFraction3D[][][index_size] = ParticleMass_Matrix[p][q]/(V_sum+MoleculesSeedPerBin[index_size]*seedMW)
								endif
							endif
						endfor
					endif

	 				// evaporation
	 				RxnStepKinetic3D_Evap[][][] = RxnStepKinetic3D_Pre2[p][q][r]*MoleFraction3D[p][q][r] // evaporation amount; molecules/particle
	 				RxnStepKinetic3D_Evap[][][] = RxnStepKinetic3D_Evap[p][q][r] > PM3D_Matrix[p][q][r] ? PM3D_Matrix[p][q][r] : RxnStepKinetic3D_Evap[p][q][r] // conditional
	 				MatrixOp/O SumKinetics = SumBeams(RxnStepKinetic3D_Evap) // for adding to the gas phase

	 				GasMolecules_matrix[][] += SumKinetics
	 				PM3D_matrix[][][] -= RxnStepKinetic3D_Evap[p][q][r]
	 			endfor
	 			MatrixOp/O ParticleMolecules_Matrix = sumBeams(PM3D_Matrix)	
	 			
	 			// Check to make sure that the timestep was small enough (added with v7.4)
	 			wavestats/q GasMolecules_matrix // this goes very negative when the mass transfer timestep is too large
	 			if(V_min < -1e5)
	 				// things are bad, so restart
	 				count += 1
	 				MassTransferTimeScaling *= round((count+1))
	 				print "restarting with mass transfer timestep = " + num2str((timestep/masstransfertimescaling)) + " seconds"
	 				GasMolecules_matrix = GasMolecules_matrix_init
	 				TotalMolecules_matrix = TotalMolecules_matrix_init
	 				ParticleMolecules_matrix = PM_matrix_init
	 				PM3D_matrix = PM3D_matrix_init
	 				m = 1 // reset m, to start over
	 				if(count > 4)
	 					abort "The mass transfer timestep must be way off...aborting"
	 				endif
	 			else
	 				// do nothing
	 			endif
	 			
	 		elseif(DynamicPartitioningMethod==1) // ZAVERI method.  New with v7.4
		 		if(timew[m-1] >= NucleationTimeEmpirical) // Assume zero particles if t < t_nuc, no condensation surface area, for nucleation
		 			// NOTE: THIS WILL NOT WORK WITH NUCLEATION AT THE MOMENT v7.4
					if(m==-1) // if first step, get initial values assuming partial equilibrium. 
						// if m<1, this will do nothing. this was an original hack that doesn't seem necesary any longer
						TotalMolecules_Matrix = TotalMolecules_matrix==0 ? 0.1 : TotalMolecules_matrix // this adds just the tiniest bit of mass
						EqmCalc(TotalMolecules_Matrix,SatConc_Matrix,MW_Matrix,Coa_guess = V_Sum,SeedConc=SeedMolecules,SeedMW=SeedMW,units=units,Tvar=Tvar,deltaHvap_matrix=deltaHvap_matrix)
						GasMolecules_matrix = TotalMolecules_Matrix - Coa_eqm/10 // Gas conc. molecules/cm3
						wavestats/q Coa_eqm
						ParticleMolecules_matrix = Coa_eqm/10 // particle conc., molecules/cm3
						PM3D_Matrix = ParticleMolecules_matrix[p][q] * dVdlogDp_norm[r] // distribute particle molecules across size bins
					endif
		 			//
		 		// First, get current total concentrations of each species accross all size bins
		 			PM3D_matrix_current = PM3D_Matrix // particle concentration, in each size bin, molecules/cm3-air; this is from the last time step
		 			matrixop/o ParticleMolecules_matrix_cur = sumBeams(PM3D_matrix_current) // particle concentration, summed over size bins
		 			GasMolecules_matrix_current = GasMolecules_matrix // gas concentration; this is from the current time step, after adjusting for gas-phase reactions
		 			GasMolecules_matrix_current2 = GasMolecules_matrix // repeat, to use later
		 			TotalMolecules_Matrix_current = TotalMolecules_matrix // total concentration, not size resolved
					//
				// calculate updated mass transfer coefficients, that will be assumed constant over this timestep	 			
				 	Knudsen_3D[][][] = 2*MeanFreePath_Matrix[p][q]/(Diameter_time[m-1][r]*1e-9) // Knudsen-number
					Fuchs_3D[][] = (0.75*alpha*(1+Knudsen_3D))/(Knudsen_3D^2+Knudsen_3D+0.283*Knudsen_3D*alpha+0.75*alpha) // Fuchs-Sutugin; Eqn. 6 in Zaveri (2008)
					Kelvin = 10^(Dk10/Diameter_time[m-1][p]) // Kelvin parameter
					deltaDynamicPartitioningPre[][][] = 4*pi*(Diameter_time[m-1][r]*1e-9/2)*Diffusivity_Matrix[p][q]*(NumParticlesByBin[r]*1e6)*Fuchs_3D[p][q][r] // units: 1/s; this is Eqn. 5 in Zaveri (2008), but with the timestep included
					//
				 // Next, get current mole fraction
			 		MoleFraction3D = (PM3D_matrix_current[p][q][r]/Na)*MW_matrix[p][q]*1e6 // *(1/NumParticlesByBin[r]), actually the mass fraction // g/m3 = molecules/cm3 * mol/molecule * g/mol * cm3/m3
			 		for(index_size=0;index_size<nSizeBins;index_size+=1)
			 			MoleFraction2D = MoleFraction3D[p][q][index_size] // this is the total mass for a given size bin, excluding the seed (units: g/m3)
			 			wavestats/q molefraction2D // add the total mass for this size bin (units: g/m3)
	 					if(V_sum==0)
							MoleFraction3D[][][index_size] = numtype(MoleFraction2D[p][q])==2 ? nan : 1
						else
				 			if(AbsorbingSeed==0)
				 				MoleFraction3D[][][index_size] = MoleFraction2D[p][q] / (V_sum + MassSeedPerBin_noSeed[index_size]) // divide the mass concentration in size bin i by the total mass concentration in size bin i
				 			else
				 				MoleFraction3D[][][index_size] = MoleFraction2D[p][q] / (V_sum + MassSeedPerBin[index_size])// divide the mass concentration in size bin i by the total mass concentration in size bin i, including the seed
				 			endif
				 		endif
			 		endfor
					// 
			 	// Get saturation concentation based on previous time step
			 		SatConc_Matrix_3D[][][] = MoleFraction3D[p][q][r]*Kelvin[r]*SatConc_Matrix[p][q]*1e-6 // // Zaveri (2008), Eqn. 71, C*,t; units: molecules/cm3 = unitless * unitless * molecules/m3
	//				SatConc_Matrix_3D[][][] = Kelvin[r]*SatConc_Matrix[p][q]*1e-6 // // Zaveri (2008), Eqn. 71, C*,t
			 		SatConc_Matrix_3D = numtype(SatConc_Matrix_3D)==2 ? 0 : SatConc_Matrix_3D
			 		//
			 	 // Now, get appropriate adaptive timestep
			 	 if(stringmatch(units,"molecules"))
					phi_GPP = GasMolecules_matrix_current[p][q] > SatConc_Matrix_3D[p][q][r] ? GasMolecules_matrix_current[p][q] : SatConc_Matrix_3D // determine which is larger, gas phase concentration or saturation concentration
					phi_GPP =(GasMolecules_matrix_current[p][q]-SatConc_Matrix_3D[p][q][r])/phi_GPP[p][q][r]
				elseif(stringmatch(units,"mass"))
					phi_GPP = GasMolecules_matrix_current[p][q] > SatConc_Matrix_3D[p][q][r] ? GasMolecules_matrix_current[p][q] : SatConc_Matrix_3D // determine which is larger, gas phase concentration or saturation concentration
					phi_GPP =(GasMolecules_matrix_current[p][q]-SatConc_Matrix_3D[p][q][r])/phi_GPP[p][q][r]
				endif
					phi_GPP = abs(phi_GPP)
					phi_GPP *= deltaDynamicPartitioningPre[p][q][r]
					matrixop/o tmatrix_2D = sumbeams(phi_GPP) // sum over all size bins
					tmatrix_2D = tmatrix_2D==0 ? nan : alpha/tmatrix_2D
					wavestats/q tmatrix_2D
					deltat_min = V_min
					MassTransferTimeScaling = max(1,ceil(timestep/deltat_min)) // lower limit of unity for number of steps for GPP; adopted from simpleSOM
					ntime_GPP = MassTransferTimeScaling // this is just using different variable names
					timestep_wave[m] = deltat_min // store the current value
				if(m==1)
					print "The GPP timestep from the first iteration = " + num2str(deltat_min)	
				endif
				// now that you have the appropriate timestep, calculate the mass transfer using MOSAIC method
					for(k=0;k<ntime_GPP;k+=1)
						// Get current mass fraction
						if(stringmatch(units,"molecules"))
	//						// don't need to do anything
							MoleFraction3D = (PM3D_matrix_current[p][q][r]/Na)*1e6 // *(1/NumParticlesByBin[r]), actually the mass fraction // g/m3 = molecules/cm3 * mol/molecule * g/mol * cm3/m3
					 		for(index_size=0;index_size<nSizeBins;index_size+=1)
					 			MoleFraction2D = MoleFraction3D[p][q][index_size] // this is the total mass for a given size bin, excluding the seed (units: g/m3)
					 			wavestats/q molefraction2D // get the total mass for this size bin (units: g/m3)
			 					if(V_sum==0)
									MoleFraction3D[][][index_size] = numtype(MoleFraction2D[p][q])==2 ? nan : 1
								else
						 			if(AbsorbingSeed==0)
						 				MoleFraction3D[][][index_size] = MoleFraction2D[p][q] / (V_sum + MassSeedPerBin_noSeed[index_size]/seedMW) // divide the mass concentration in size bin i by the total mass concentration in size bin i
						 			else
						 				MoleFraction3D[][][index_size] = MoleFraction2D[p][q] / (V_sum + MassSeedPerBin[index_size]/(seedMW))// divide the mass concentration in size bin i by the total mass concentration in size bin i, including the seed
						 			endif
						 		endif
					 		endfor
					 		// Get current saturation concentration and saturation ratio, based on mole fraction
				 			SatConc_Matrix_3D[][][] = MoleFraction3D[p][q][r]*Kelvin[r]*SatConc_Matrix[p][q]*1e-6 // // Zaveri (2008), Eqn. 71, C*,t; units: molecules/cm3 = unitless * unitless * molecules/m3
					 		SatConc_Matrix_3D = numtype(SatConc_Matrix_3D)==2 ? 0 : SatConc_Matrix_3D
							SatRatio_3D = PM3D_matrix_current==0 ? 1 : SatConc_Matrix_3D[p][q][r] / PM3D_matrix_current[p][q][r] // Saturation ratio; SatConc / [OA]; unitless: molecules/cm3 / molecules/cm3; very important to have the conditional here
							// calculate gas-phase concentration at t+1 (Eqn. 23 from Zaveri): implemented in steps
							deltaDynamicPartitioning[][][] = 1*(PM3D_matrix_current[p][q][r]) / (1+(timestep/ntime_GPP)*deltadynamicpartitioningpre[p][q][r]*SatRatio_3D[p][q][r]) // first part of sum in numerator
								// note: for the above equation, Zaveri has this with a negative sign. But it seems as if it should be positive. 
							matrixOp/o sumDeltaDP = sumbeams(deltaDynamicPartitioning) // sum in numerator, sums over all size bins
							GasMolecules_matrix_current2[][] = TotalMolecules_Matrix_current - sumDeltaDP // numerator in Eqn. 23
							deltaDynamicPartitioning[][][] = (deltaDynamicPartitioningPre[p][q][r])/(1+(timestep/ntime_GPP)*deltadynamicpartitioningpre[p][q][r]*SatRatio_3D[p][q][r]) //  denominator 2nd term (the term to be summed) in Eqn. 23
							matrixOp/o sumDeltaDP = sumbeams(deltaDynamicPartitioning) // sum in denominator from Eqn. 23; sum over all sizes
							GasMolecules_matrix_current2 /= (1+(timestep/ntime_GPP)*sumDeltaDP) // this is the gas concentration at timestep t+1
							// calculate particle-phase concentration at t+1 (Eqn. 19 from Zaveri). Uses just-calculated gas-phase concentration
							PM3D_matrix_current = (PM3D_matrix_current[p][q][r] + (timestep/ntime_GPP)*deltadynamicpartitioningpre[p][q][r]*GasMolecules_matrix_current2[p][q]) // numerator
							PM3D_matrix_current /= (1 + (timestep/ntime_GPP)*deltadynamicpartitioningpre[p][q][r]*SatRatio_3D[p][q][r]) // Particle concentration at t+1; divide numerator by denominator
							// reset values for next iteration
							GasMolecules_matrix = GasMolecules_matrix_current2 // gas-phase concentrations
							PM3D_matrix = PM3D_matrix_current // size bin specific particle concentrations
							MatrixOp/O ParticleMolecules_Matrix = sumBeams(PM3D_Matrix_current)	 // total SOA concentration across all size bins; molecules/cm3-air	
							TotalMolecules_Matrix_current = ParticleMolecules_Matrix + GasMolecules_Matrix // recalculate, just to be sure
						elseif(stringmatch(units,"mass"))
	//						// don't need to do anything
							MoleFraction3D = (PM3D_matrix_current[p][q][r]/Na)*MW_matrix[p][q]*1e6 // *(1/NumParticlesByBin[r]), actually the mass fraction // g/m3 = molecules/cm3 * mol/molecule * g/mol * cm3/m3
					 		for(index_size=0;index_size<nSizeBins;index_size+=1)
					 			MoleFraction2D = MoleFraction3D[p][q][index_size] // this is the total mass for a given size bin, excluding the seed (units: g/m3)
					 			wavestats/q molefraction2D // get the total mass for this size bin (units: g/m3)
			 					if(V_sum==0)
									MoleFraction3D[][][index_size] = numtype(MoleFraction2D[p][q])==2 ? nan : 1
								else
						 			if(AbsorbingSeed==0)
						 				MoleFraction3D[][][index_size] = MoleFraction2D[p][q] / (V_sum + MassSeedPerBin_noSeed[index_size]) // divide the mass concentration in size bin i by the total mass concentration in size bin i
						 			else
						 				MoleFraction3D[][][index_size] = MoleFraction2D[p][q] / (V_sum + MassSeedPerBin[index_size])// divide the mass concentration in size bin i by the total mass concentration in size bin i, including the seed
						 			endif
						 		endif
					 		endfor
					 		// Get current saturation concentration and saturation ratio, based on mole fraction
				 			SatConc_Matrix_3D[][][] = MoleFraction3D[p][q][r]*Kelvin[r]*SatConc_Matrix[p][q]*1e-6 // // Zaveri (2008), Eqn. 71, C*,t; units: molecules/cm3 = unitless * unitless * molecules/m3
					 		SatConc_Matrix_3D = numtype(SatConc_Matrix_3D)==2 ? 0 : SatConc_Matrix_3D
							SatRatio_3D = PM3D_matrix_current==0 ? 1 : SatConc_Matrix_3D[p][q][r] / PM3D_matrix_current[p][q][r] // Saturation ratio; SatConc / [OA]; unitless: molecules/cm3 / molecules/cm3
							// calculate gas-phase concentration at t+1 (Eqn. 23 from Zaveri): implemented in steps
							deltaDynamicPartitioning[][][] = 1*(PM3D_matrix_current[p][q][r]) / (1+(timestep/ntime_GPP)*deltadynamicpartitioningpre[p][q][r]*SatRatio_3D[p][q][r]) // first part of sum in numerator
								// note: for the above equation, Zaveri has this with a negative sign. But it seems as if it should be positive. 
							matrixOp/o sumDeltaDP = sumbeams(deltaDynamicPartitioning) // sum in numerator, sums over all size bins
							GasMolecules_matrix_current2[][] = TotalMolecules_Matrix_current - sumDeltaDP // numerator in Eqn. 23
							deltaDynamicPartitioning[][][] = (deltaDynamicPartitioningPre[p][q][r])/(1+(timestep/ntime_GPP)*deltadynamicpartitioningpre[p][q][r]*SatRatio_3D[p][q][r]) //  denominator 2nd term (the term to be summed) in Eqn. 23
							matrixOp/o sumDeltaDP = sumbeams(deltaDynamicPartitioning) // sum in denominator from Eqn. 23; sum over all sizes
							GasMolecules_matrix_current2 /= (1+(timestep/ntime_GPP)*sumDeltaDP) // this is the gas concentration at timestep t+1
							// calculate particle-phase concentration at t+1 (Eqn. 19 from Zaveri). Uses just-calculated gas-phase concentration
							PM3D_matrix_current = (PM3D_matrix_current[p][q][r] + (timestep/ntime_GPP)*deltadynamicpartitioningpre[p][q][r]*GasMolecules_matrix_current2[p][q]) // numerator
							PM3D_matrix_current /= (1 + (timestep/ntime_GPP)*deltadynamicpartitioningpre[p][q][r]*SatRatio_3D[p][q][r]) // Particle concentration at t+1; divide numerator by denominator
							// reset values for next iteration
							GasMolecules_matrix = GasMolecules_matrix_current2 // gas-phase concentrations
							PM3D_matrix = PM3D_matrix_current // size bin specific particle concentrations
							MatrixOp/O ParticleMolecules_Matrix = sumBeams(PM3D_Matrix_current)	 // total SOA concentration across all size bins; molecules/cm3-air	
							TotalMolecules_Matrix_current = ParticleMolecules_Matrix + GasMolecules_Matrix // recalculate, just to be sure			
						endif
					endfor
	//				GasMolecules_matrix = GasMolecules_matrix_current2 // gas-phase concentrations
	//				PM3D_matrix = PM3D_matrix_current[p][q][r] // size bin specific particle concentrations
	//				MatrixOp/O ParticleMolecules_Matrix = sumBeams(PM3D_Matrix)	 // total SOA concentration across all size bins; molecules/cm3-air		 		else
		 			// no partitioning at all
		 		endif // conditional for empirical nucleation
		 	endif // conditional for Zaveri method for GPP
		endif // conditional for instantaneous vs. dynamic partitioning
		
// 11ff: Oligomerization reactions (these work okay, but could be a lot better)
		if(OligomerizationIsOn==1)
			SOM_Oligomerization(ParticleMolecules_Matrix,krxn_matrix_Olig,timestep,Np,VolumeOrgPerParticle,FirstTime=FirstTime,method=2)
			wave ParticleMolecules_Matrix
		endif

// 11fff: Final concentrations from this iteration		
		TotalMolecules_Matrix = GasMolecules_matrix + ParticleMolecules_Matrix
		GasMolecules_Time[][][m] = GasMolecules_Matrix[p][q]
		ParticleMolecules_Time[][][m] = ParticleMolecules_Matrix[p][q]
		TotalMolecules_Time[][][m] = TotalMolecules_Matrix[p][q]
		
		// HOM concentration and nucleation (IN PROGRESS; does not do anything important yet)
		wavestats/q logcstar_matrix
		if(V_min <= HOM_logCstar)
			extract/FREE gasmolecules_matrix, HOMwave, logcstar_matrix <= HOM_logCstar // get concentration of "HOM" species
		else // must have something...so set the lowest volatility bin as the "HOM" bin
			extract/FREE gasmolecules_matrix, HOMwave, logcstar_matrix <= V_min  // get concentration of "HOM" species
		endif
		wavestats/q HOMwave
		HOM = V_sum
		HOM_time[m] = HOM
		NucleationRate = nuc_a1*(HOM/1e7)^(nuc_a2 + nuc_a5/(HOM/1e7))
		NucleationRate_time[m] = NucleationRate
	
		// Suspended SOA concentration
		TempWave = ParticleMolecules_Matrix[p][q] * 1e6 * 1e6 * MW_Matrix[p][q]/Na
		wavestats/q tempwave
		Coa_time[m] = V_sum // Current SOA concentration (ug/m^3)
		
		// Correct for particle deposition
		TempWave = (ParticleMolecules_Matrix[p][q] + WallParticleMolec_Matrix[p][q]) * 1e6 * 1e6 * MW_Matrix[p][q]/Na
		wavestats/q tempwave
		Coa_time_wallcorr[m] = V_sum // Current OA concentration, corrected for particle wall-loss (ug/m^3)
		
		// loop over sizes to get new diameters...this may not be quite perfect yet (09/29/14)
		variable/D VolOrgTot = 0
		if(KineticMassTransfer==0) // eqm. partitioning)
			PM3D_matrix = 0
			PM3D_matrix[][][0] = ParticleMolecules_Matrix[p][q]
		endif
		for(index_size=0;index_size<nSizeBins;index_size+=1)
			ParticleVolume_matrix = (PM3D_matrix[p][q][index_size]/NpPerBin_time[m][index_size])*(MW_Matrix[p][q]/(Na*Density_Base*1e6)) // volume of each component in a size bin, per particle
			wavestats/Q ParticleVolume_Matrix // Add all components in a size bin
			VolumeOrgPerParticle = V_sum // m^3/p 
			Diameter_time[m][index_size] = 1e9*((6/pi)*(VolumeOrgPerParticle+VolumeSeedPerBin[index_size]))^(1/3) // add the seed volume per particle, and convert to diameter
			VolumeOrgPerBin[index_size] = VolumeOrgPerParticle
			VolOrgTot += VolumeOrgPerParticle
		endfor	

		Dp_time[m] = 1e9*((6/pi)*(VolOrgTot+VolumeSeed))^(1/3)	// for polydisperse simulations, this is an overall "average" diameter
				
// 11g: O:C calculation		
			wavestats/Q ParticleMolecules_matrix
			C_matrix_temp = C_matrix*ParticleMolecules_matrix/V_Sum
			O_matrix_temp = O_matrix*ParticleMolecules_matrix/V_Sum
			H_matrix_temp = H_matrix*ParticleMolecules_matrix/V_sum
			wavestats/Q C_Matrix_temp
			Ave_Nc = V_Sum
			wavestats/Q O_matrix_temp
			Ave_No = V_sum
			wavestats/Q H_matrix_temp
			Ave_Nh = V_sum
			O2C_time[m] = Ave_No/Ave_Nc
			H2C_time[m] = Ave_Nh/Ave_Nc
		
// 11h: DILUTION
		if(DilutionVar != 0)
			Totalmolecules_Matrix *= 1-(DilutionVar/100)*(TimeStep/60/60)
			GasMolecules_Matrix *= 1-(DilutionVar/100)*(TimeStep/60/60)
			ParticleMolecules_Matrix *= 1-(DilutionVar/100)*(TimeStep/60/60)
			SeedMolecules *= 1-(DilutionVar/100)*(TimeStep/60/60)
			MoleculesSeedPerBin *= 1-(DilutionVar/100)*(TimeStep/60/60)
			SeedConc_time[m] = SeedConc_time[m-1]*(1-(DilutionVar/100)*(TimeStep/60/60))
			NpPerBin_time[m] = NpPerBin_time[m-1]*(1-(DilutionVar/100)*(TimeStep/60/60)) // v7.4.3
		endif
		
// 11i: OH Exposure & lifetimes
		OH_exposure[m-1] = OH_exposure[m-2]+OH_wave[m-1]*timestep // molecules cm^-3 hr^-1
		Lifetimes[m-1] = Lifetimes[m-2] + krxn_matrix[0][Nox_precursor] * OH_wave[m-1] * timestep // based on the precursor
		Lifetimes[m] = Lifetimes[m-2] + 2* krxn_matrix[0][Nox_precursor] * OH_wave[m-1] * timestep // based on the precursor
		TimeW[m] = TimeW[m-1]+timestep/3600
		if(m==100 && quiet == 0)
			t22 = ticks
			print "Calculations will likely take ~ " + num2str(((t22-t11)/60)*nsteps/100) + " seconds"
		endif
		FirstTime = 0
		m+=1
//	endfor
	while(Timew[m] < MaxTime_hours) // v7.4 change to do loop
// END KINETICS	
	
	t22 = ticks
	if(quiet==0)
		print "Kinetics took "+num2str((t22-t11)/60) + " seconds to execute"
	endif
	
	ParticleMass_Matrix = ParticleMolecules_Matrix[p][q][nsteps] * 1e6 * 1e6 * MW_Matrix[p][q]/Na
	GasMass_Matrix = GasMolecules_Matrix[p][q][nsteps] * 1e6 * 1e6 * MW_Matrix[p][q]/Na	
	TotalMass_matrix = particlemass_matrix + gasmass_matrix
	ParticleMass_Time = ParticleMolecules_Time[p][q][r]*1e12*MW_Matrix[p][q]/Na
	GasMass_Time = GasMolecules_Time[p][q][r]*1e12*MW_Matrix[p][q]/Na
	duplicate/o TimeW TimeW_minutes
	TimeW_minutes = TimeW*60
	
	Make/o/d/n=(nsteps,nSizeBins) logDp_time
	logDp_time = log(diameter_time)
	differentiate/DIM=1 logDp_time/D=dlogDp_time
	dNdlogDp_time[][] = NpPerBin_time/dlogDp_time
	dVdlogDp_time[][] = dNdlogDp_time[p][q]*(pi/6)*(diameter_time[p][q]*1e-9)^3/1e18
	
	MatrixOp/O Np_time = sumRows(NpPerBin_time)
	Np_time = timew < NucleationTimeEmpirical ? 0 : Np_time
	
// 12: Do some random things for graphing final distribution of products
	variable type = 1
	if(type == 0) // molecules
		wavestats/Q totalmolecules_matrix;
		GasMolecules_Matrix_Log = log(GasMolecules_Matrix/V_sum)
		TotalMolecules_Matrix_Log = log(TotalMolecules_Matrix/V_sum)
		wavestats/Q ParticleMolecules_Matrix
		Particle_Fraction_Log = log(ParticleMolecules_Matrix/V_sum)
		Particle_Fraction = ParticleMolecules_Matrix/V_Sum
		Particle_div_Total = (ParticleMolecules_Matrix/TotalMolecules_Matrix)
			// Scale the graph showing the composition results
		DoWindow/F Graph11
		ModifyImage TotalMolecules_Matrix_Log ctab= {-4,0,Rainbow,0}; 
		ModifyImage Particle_Fraction_Log ctab= {-4,0,Rainbow,0}; 
		ModifyImage GasMolecules_Matrix_Log ctab= {-4,0,Rainbow,0}; 
	else // mass
		wavestats/Q totalmass_matrix
		GasMolecules_Matrix_Log = log(GasMass_Matrix/V_sum)
//		Particle_Fraction_Log = log(ParticleMass_Matrix/V_sum)
		TotalMolecules_Matrix_Log = log(TotalMass_Matrix/V_sum)
		wavestats/Q ParticleMass_Matrix
		Particle_Fraction_Log = log(ParticleMass_Matrix/V_sum)
		Particle_Fraction = ParticleMass_Matrix/V_Sum
		Particle_div_Total = (ParticleMass_Matrix/TotalMass_Matrix)
	endif

	variable/D Coa_var = Coa_time[0] // Final OA mass concentration
	
//	deltahc_time = (gasmolecules_time[0][0][0] - GasMolecules_time[0][0][p])* 1e6 * 1e6 * MW_Matrix[0][0]/Na
	deltahc_time = (gasmolecules_time[0][Nox_precursor][0] - GasMolecules_time[0][Nox_precursor][p])* 1e6 * 1e6 * MW_Matrix[0][Nox_precursor]/Na
//	deltaHCfrac_time = deltaHC_time/(gasmolecules_time[0][0][0] * 1e6 * 1e6 * MW_Matrix[0][0]/Na)
	deltaHCfrac_time = deltaHC_time/(gasmolecules_time[0][Nox_precursor][0] * 1e6 * 1e6 * MW_Matrix[0][Nox_precursor]/Na)
	deltaHCfrac_time = abs(1-deltaHCfrac_time)
//	deltaHCfrac_time = 1-deltaHCfrac_time
//	HC_ppb_time = 1e9*GasMolecules_time[0][0][p]/(NumberDensityAir*1e-6)//1e3*(Ctot_ppm - deltahc_time/((2.46e25*1e-6/Na)*(MWstart*1e6)))
//	HC_ugm3_time = GasMolecules_time[0][0][p] * 1e6 * 1e6 * MW_Matrix[0][0]/Na
	HC_ppb_time = 1e9*GasMolecules_time[0][Nox_precursor][p]/(NumberDensityAir*1e-6)//1e3*(Ctot_ppm - deltahc_time/((2.46e25*1e-6/Na)*(MWstart*1e6)))
	HC_ugm3_time = GasMolecules_time[0][Nox_precursor][p] * 1e6 * 1e6 * MW_Matrix[0][Nox_precursor]/Na

	NVAR YieldCalc
	if(yieldcalc==0)
		Yield_time = Coa_time/deltaHC_time
	else
		Yield_time = (Coa_time-Coa_time[0])/deltaHC_time // subtract "primary" SOA from precursor
	endif
	Yield_time = (Coa_time < 1e-7 ? 0 : Yield_time)
	Yield2_time = (Coa_time - Coa_var)/(deltaHC_Time)
	Yield2_Time = (Coa_time < 1e-7 ? 0 : Yield2_time)
	Ncarbons_wave += 0.5
	Noxygens_wave -= 0.5
	
	wavestats/Q Coa_time
	variable/D tempvar = V_max
	wavestats/Q Yield_time
	variable/D YieldVar = Yield_time[nsteps]//V_Max
	
// Print results
	if(quiet==0)
		print "[SOA] =  " + num2str(tempvar)+ "ug/m^3 and Yield = " + num2str(YieldVar) + " for Nc = " + num2str(nC) + " and O:C = " + num2str(O2C_time[nsteps]) + " and ProbOx = " + num2str(ProbOx1) + ";"+ num2str(ProbOx2) + ";"+ num2str(ProbOx3) + ";"+ num2str(ProbOx4) + ";"
	endif
//	print "Chi Square values are " + num2str(ChiSq_Coa) + " for Coa and " + num2str(ChiSq_O2C) + " for O2C."
	
// Extract into VBS framework (1 for mass, 0 for molecules)
//	VBS(1) // 1 for mass, 0 for molecules
//	ParticleMass_Time = ParticleMolecules_Time[p][q][r] /(1e-6 * 1e-6 * Na / MW_Matrix[p][q])
//	TotalMass_Time = TotalMolecules_Time[p][q][r] /(1e-6 * 1e-6 * Na / MW_Matrix[p][q])
//	VBS_time(1)

	if(InterpToExperiment==1)
		setdatafolder root:
		// interpolate Coa_time, Coa_time_wallcorr and O2C_time to Experiment time (fittime)
		Interpolate2/T=1/I=3/Y=Coa_time_Interp/X=FitTime TimeW, Coa_time
		Interpolate2/T=1/I=3/Y=Coa_time_Wallcorr_Interp/X=FitTime TimeW, Coa_time_wallcorr
		Interpolate2/T=1/I=3/Y=deltaHC_time_Interp/X=FitTime TimeW, deltaHC_time

		wavestats/q O2C_time
		if(V_npnts > 10)
			Interpolate2/T=1/I=3/Y=O2C_time_Interp/X=FitTime TimeW, O2C_time
		else
			duplicate/o Coa_time_interp O2C_time_interp
			O2C_time_interp = nan
		endif
		wave Coa_time_Interp
		setscale/P x 0, (timestep/60/60), "hours", Coa_time_interp
		setscale/P d, 0, 0, "ug/m\S3\M", Coa_time_interp
		wave Coa_time_WallCorr_Interp
		setscale/P x 0, (timestep/60/60), "hours", Coa_time_WallCorr_interp
		setscale/P d, 0, 0, "ug/m\S3\M", Coa_time_interp
		setscale/P x 0, (timestep/60/60), "hours", O2C_time_interp
		setscale/P d, 0, 0, "O:C", O2C_time_interp
		Coa_time_interp = FitTime>TimeW[nSteps-1] ? nan : Coa_time_interp
		Coa_time_WallCorr_interp = FitTime>TimeW[nSteps-1] ? nan : Coa_time_WallCorr_interp
		O2C_time_interp = FitTime>TimeW[nSteps-1] ? nan : O2C_time_interp
		// Calculate ChiSquare Values for Coa and O:C
		variable/G ChiSq_Coa
		variable/G ChiSq_O2C
		wave Coa_experiment
		if(waveexists(Coa_experiment)==1)
	//		EstimateError(1,2,30)
			wave Coa_experiment
			wave O2C_experiment
			wave Coa_experiment_err
			wave O2C_experiment_err
			wave Coa_experiment_mask
			wave O2C_experiment_mask
			duplicate/O Coa_experiment ChiSqTest
//			Coa_experiment_err = max(0.5,0.15*coa_experiment)
			ChiSqTest = ((Coa_time_interp-Coa_experiment)/Coa_experiment_err)^2 
			ChiSqTest *= Coa_experiment_mask
			wavestats/Q ChiSqTest
			ChiSq_Coa = (1/(numpnts(Coa_time_interp)-6-1))*V_sum
			O2C_experiment_err = max(0.05,0.3*O2C_experiment)
			ChiSqTest = ((O2C_time_interp-O2C_experiment)/O2C_experiment_err)^2
			ChiSqTest *= O2C_experiment_mask
			wavestats/Q ChiSqTest
			ChiSq_O2C = V_sum
		else
			ChiSq_Coa = nan
			ChiSq_O2C = nan
		endif
	endif

	// Kill unneeded waves	
	if(cleanup==1)
		SOM_KillWaves()
	endif
	
END

//****************************************************************************************************************************
Function SOM_ParticleWallLossRate(dp,particlewallloss)
	// written as sub-function on 09/29/14
	// uses empirical fit to data from Loza et al., ACP, 2012
	variable dp // nm
	variable particlewallloss // scaling factor
	
	variable particleWLR
	if(dp < 40)	// this condition is made up to deal with UCR data at early times in an experiment
		particleWLR = 0
	elseif(dp < 280)	// nm
		particleWLR = ParticleWallLoss*10^(-0.50785 - 2.048*log(dp) + 0.99519*log(dp)^2 - 0.30528*log(dp)^3)/60	// per second, from Loza 2012
	else
		particleWLR = ParticleWallLoss*10^(68.58 - 87.421*log(dp) + 34.237*log(dp)^2 - 4.3463*log(dp)^3)/60	// per second, from Loza 2012
	endif
	
	return particleWLR
End

//********************************************************************************************************
Function SetO3_Chan(a,b,pt,scale)
	variable a,b,pt, scale

	wave O3_conc_wave
	wavestats/q O3_conc_wave
	
	o3_conc_wave = a*1e12*(1-exp(-3*x/b))
	
	variable m
	variable npnts = numpnts(O3_conc_wave)
	
	for(m=pt;m<=npnts;m+=1)
		O3_conc_wave[m] = O3_conc_wave[m-1]*((npnts-m)/npnts)^scale
	endfor
end

//********************************************************************************************************
Function FindTarget(value,tol, upper, lower)
	variable value, tol,upper, lower
	
	lower = upper/2
	variable dif
	wave Coa_time
	nvar ctot_ppm
	variable keepgoing = 1
	variable Coa_current
	variable guess = (upper + lower)/2
	
	do
	ctot_ppm = guess
	som_v1()
	coa_current = coa_time[numpnts(coa_time)]
	dif = coa_current - value
		if(abs(dif) < tol)
			keepgoing = 0
		elseif(dif > 0)
			guess = (ctot_ppm + lower)/2
			lower = lower
			upper = ctot_ppm
		else
			guess = (ctot_ppm + upper)/2
			lower = ctot_ppm
			upper = upper
		endif
		print guess
	while(keepgoing==1)
	//	print keepgoing
	
End

//*************************************************************************************
Function MakeFancyGraph(matrix,xwave)
	wave matrix	// a 3D matrix, with time as layers
	wave xwave
	
	variable nRows = dimsize(matrix,0)
	variable nCols = dimsize(matrix,1)
	variable nLayers = dimsize(matrix,2)
	
	variable spacing = 1/nRows
	variable i, j
	string axisname
	display
	for(i=0;i<nRows;i+=1)
		for(j=0;j<nCols;j+=1)
			if(i==0)
				appendtograph matrix[i][j][] vs xwave
			else
				axisname = "L"+num2str(i)
				appendtograph/L=$axisname matrix[i][j][] vs xwave
			endif
		endfor
	endfor
	stackallaxes("",0,0)	
End

//***********************************************************************************************************
Function SOM_GetAtomicRatios(ParticleMolecules_Time,C_matrix,O_matrix,H_matrix)
// LEGACY
	wave ParticleMolecules_Time	// 3D matrix with particle species concentration as function of time
	wave C_matrix	// number of carbon atoms per molecule
	wave O_Matrix	// number of oxygen atoms per molecule
	wave H_Matrix	// number of hydrogen atoms per molecule
	
	variable nRows = dimsize(ParticleMolecules_Time,0)
	variable nCols = dimsize(ParticleMolecules_Time,1)
	variable nLayers = dimsize(ParticleMolecules_Time,2)
	variable i
	variable C_var
	variable O_var
	variable H_var
	variable t1, t2
	t1 = ticks
	make/o/d/n=(nRows,nCols,nLayers)/FREE Carbon=0, Oxygen=0, Hydrogen=0
	make/o/d/n=(nRows,nCols)/FREE SingleLayer
	make/o/d/n=(nLayers) O2C_time, H2C_time
	
	Carbon = ParticleMolecules_Time[p][q][r]*C_matrix[p][q]
	Oxygen = ParticleMolecules_Time[p][q][r]*O_matrix[p][q]
	Hydrogen = ParticleMolecules_Time[p][q][r]*H_matrix[p][q]
	
	for(i=0;i<nLayers;i+=1)
		SingleLayer = Carbon[p][q][i]
		wavestats/q SingleLayer
		C_var = V_sum
		SingleLayer = Oxygen[p][q][i]
		wavestats/q SingleLayer
		O_var = V_sum
		SingleLayer = Hydrogen[p][q][i]
		wavestats/q SingleLayer
		H_var = V_sum
		O2C_time[i] = O_var/C_var
		H2C_time[i] = H_var/C_var
	endfor
	t2 = ticks
	print ((t2-t1)/60)
End

//************************************************************************************************************
Function SOM_RateCoefficients(C_matrix,O_matrix,krxn_parent,O3_yn,[TempK,method,NO3adjustment,Ea])
	// OH reaction rate coefficients: increase with Nc, increase with added oxygen (to a point)
	// based on structure activity relationship from Kwok and Atkinson (1995), Atmos Environ, 29, 1685.
	// --> updated based on GECKO-A
	// can set different krxn for the parent species
	// if cals are for O3, treat differently
	wave C_matrix
	wave O_matrix
	variable krxn_parent
	variable O3_yn
	variable TempK
	variable method	// how do you want to specify your kOH matrix? 0 = old method, 1 = "new" method, 2 = "newer" method that is more continuous, 3 = kOH(X,Y) = b*kOH(0,0)
	variable NO3adjustment	// there is some evidence (from GECKO) that CxO3 species react about as fast as the parent compound. 0 = default parameterization; 1 = CxO3 = CxO0
	variable Ea // activation energy, kJ/mol
	
	if(paramisdefault(TempK))
		TempK = 298
	endif
	if(paramisdefault(method))
		method = 2
	endif
	if(paramisdefault(NO3adjustment))
		NO3adjustment = 0
	endif
	
	variable Ncarbons = dimsize(C_matrix,0)
	variable Noxygens = dimsize(C_matrix,1)
	make/o/d/n=(Ncarbons,Noxygens) krxn_matrix = 0

	NVAR FragSlope = root:FragSlope // fixed on 8/4/18 for use with SOM-MC.
	variable kbase = 0
	variable kprim = 1.36e-13//1.43e-13
	variable ksec = 9.34e-13//8.38e-13
	variable ktert = 1.94e-12//1.82e-12
	variable Foh = 3.5
	variable Fch2 = 1.23
	variable Fco = 0.75
	variable Fch2co = 3.9
	variable Nc_var
	variable No_var
	variable f1 = 1
	variable f2 = 3.7
	variable a_kin, b
	variable m, n

//		variable Ea //activation energy, kJ/mol.K
		variable multiplier //= (TempK^2)*(exp(-1*Ea*1000/(8.314*TempK)))

	variable kinetics = method
	if(kinetics == 1) // use new formulation	
		for(m=1;m<=Ncarbons;m+=1)
			for(n=1;n<=Noxygens;n+=1)
				kbase = 0
				Fch2 = 1.23
				Nc_var = C_matrix[m-1][n-1]
				No_var = O_matrix[m-1][n-1]
				if(No_var == 0)
					fch2 = 1.23
					if(Nc_var == 1) // methane
						kbase = kprim
					else
					kbase += 2*kprim*Fch2
						if(Nc_var < 5)
							kbase += (Nc_var-2)*ksec*Fch2*1
						else
							kbase += (Nc_var-4)*ksec*(Fch2)^2
							kbase += (2)*ksec*(Fch2)*1
						endif
					endif
				elseif(Nc_var - No_var >= 2) // molecules with >0 oxygens, but less than or equal Nc-2 oxygens
					a_kin = No_var/(Nc_var-2)
					if(No_var > (Nc_var-2)/2)
						b = 1 - abs((Nc_var-2)/2-No_var)/((Nc_var-2)/2)
					else
						b = 1
					endif
					kbase += 2*kprim*((1-a_kin)*fch2+a_kin*((fco+fch2co+foh)/3))
						if(Nc_var < 5)
							fch2 = 1
						else
							fch2 = (1*2 + (Nc_var-2)*fch2)/Nc_var
						endif
					kbase += (Nc_var-2-No_var)*ksec*((1-a_kin)*fch2*fch2 + 1*a_kin*(b*fco*fch2co+fco*fch2+b*fch2*fch2co)/(2*b+1))
					kbase += 0.5*No_var*ktert*((1-a_kin)*foh*fch2*fch2 + a_kin*foh*((fco+b*fch2co+fch2)/(b+2))^2)
				elseif(Nc_var-No_var == 1) // molecules with Nc-1 oxygens
					kbase += kprim*((fco+fch2co+foh)/3)
						if(Nc_var < 5)
							fch2 = 1
						else
							fch2 = (1*2 + (Nc_var-2)*fch2)/Nc_var
						endif
					kbase += 0.5*ksec*(foh*((fco+fch2co+fch2)/3))
					kbase += 0.5*No_var*ktert*(foh*((fco+fch2)/2)^2)
				elseif(Nc_var-No_var == 0) // molecules with Nc = No oxygens
						if(Nc_var < 5)
							fch2 = 1
						else
							fch2 = (1*2 + (Nc_var-2)*1.23)/Nc_var
						endif
					kbase += ksec*(foh*((fco+fch2co+fch2)/3))
					kbase += 0.5*No_var*ktert*(foh*((fco+fch2)/2)^2)
				elseif(2*Nc_var-No_Var <= 0) // molecules with 2*Nc = No oxygens
					kbase = 1e-14
				else // for molecules with O:C > 1 but < 2
						if(Nc_var < 5)
							fch2 = 1
						else
							fch2 = (1*2 + (Nc_var-2)*fch2)/Nc_var
						endif
					kbase += 0.5*No_var*ktert*(foh*((fco+fch2)/2)^2) * (1 - No_var/(2*Nc_var))
				endif
				krxn_matrix[m-1][n-1] = kbase			
			endfor
		endfor
	elseif(kinetics == 0) // use old formulation	
		for(m=1;m<=Ncarbons;m+=1)
			for(n=1;n<=Noxygens;n+=1)
				kbase = 0
				Nc_var = C_matrix[m-1][n-1]
				No_var = O_matrix[m-1][n-1]
				if(No_var == 0)
					if(Nc_var == 1) // methane
						kbase = kprim
					else
					kbase += 2*kprim*Fch2
					kbase += (Nc_var-2)*ksec*(Fch2)^2
					endif
				elseif(Nc_var - No_var >= 2) // molecules with >0 oxygens, but less than or equal Nc-2 oxygens
					kbase += 2*kprim*Fch2 // added on 03/07/12; v1.1a
					kbase += (Nc_var-2-No_var)*ksec*Fch2^2 // added on 03/07/12; v1.1a
					kbase += No_var*0.5*ktert*Foh
				elseif(Nc_var-No_var == 1) // molecules with Nc-1 oxygens
					kbase += kprim
					kbase += No_var*0.5*ktert*Foh
				elseif(Nc_var-No_var == 0) // molecules with Nc = No oxygens
					kbase += No_var*0.5*ktert*Foh
				elseif(2*Nc_var-No_Var <= 0) // molecules with 2*Nc = No oxygens
					kbase = 1e-14
				else // for molecules with O:C > 1 but < 2
					kbase += (2*Nc_var-No_var)*0.5*ktert*Foh
				endif
				krxn_matrix[m-1][n-1] = kbase
			endfor
		endfor
	elseif(kinetics==2)	// formulation from 9/1/2013 based on further analysis
		krxn_matrix = 0
		make/o/d/n=(3) W_coef_Carbon = {-15.103,-3.9481,-0.79796}
		make/o/d/n=(4) W_coef_O2C = {0.73201,3.7216,0.56103,0.35511}
		krxn_matrix[][0] = 10^(W_coef_Carbon[0]+W_coef_Carbon[1]*C_matrix[p][0]^W_coef_Carbon[2])
		krxn_matrix[][1,] = krxn_matrix[p][0]*(W_coef_O2C[0]+W_coef_O2C[1]*exp(-(((O_matrix[p][q]/C_matrix[p][q])-W_coef_O2C[2])/W_coef_O2C[3])^2))
		Ea = 1 // activation energy, kJ/mol.K
		multiplier = (TempK^2)*(exp(-1*Ea*1000/(8.314*TempK)))
		krxn_matrix[][] *= multiplier
	elseif(kinetics==3)	// adjustment on kinetics = 2 with new parameters
		krxn_matrix = 0
		make/o/d/n=(3) W_coef_Carbon = {-15.103,-3.9481,-0.79796}
//		make/o/d/n=(4) W_coef_O2C = {0.73201,3.4,0.4,0.35511}
		make/o/d/n=(4) W_coef_O2C = {0.73201,3.4,0.3,0.35511}
		krxn_matrix[][0] = 10^(W_coef_Carbon[0]+W_coef_Carbon[1]*C_matrix[p][0]^W_coef_Carbon[2])
		krxn_matrix[][1,] = krxn_matrix[p][0]*(W_coef_O2C[0]+W_coef_O2C[1]*exp(-(((O_matrix[p][q]/C_matrix[p][q])-W_coef_O2C[2])/W_coef_O2C[3])^2))
		krxn_matrix[1,][1,] = krxn_matrix[0][q]*(C_matrix[p][q]/(Ncarbons+1))
		Ea = 1 // activation energy, kJ/mol.K
		multiplier = (TempK^2)*(exp(-1*Ea*1000/(8.314*TempK)))
		krxn_matrix[][] *= multiplier
	elseif(kinetics==4)	//based on analysis of GECKO binned kOH values for MILAGRO outflow and Nc >= 10 (no aromatics)
		krxn_matrix = 0
		make/o/d/n=(3)/FREE W_coef_Carbon = {-15.103,-3.9481,-0.79796}	// parameters to describe variation in kOH for alkanes (straight chain)
		make/o/d/n=(nCarbons)/FREE sigma, x0, scalingfactor	// parameters of log normal distribution with f(x) = (scalingfactor/(sigma*sqrt(2*pi))) * exp(-1*(ln(x)-ln(x0))^2/(2*sigma^2))
		make/o/d/n=(nCarbons,nOxygens)/FREE logNorm = 0
		sigma = C_matrix[p][0] <= 15 ? (0.0214*C_matrix[p][0]+0.5238) : (-0.115*C_matrix[p][0]+2.695)	// spread is function of carbon number
		x0 = C_matrix[p][0] <= 15 ? (0.0314*C_matrix[p][0] + 0.9871) : (0.25*C_matrix[p][0] - 2.183)		// peak diameter is f(nC)
		scalingfactor = -0.2583*C_matrix[p][0] + 5.8944												// scaling factor is f(nC)
		lognorm = 1+(scalingfactor[p]/(sigma[p]*sqrt(2*pi))) * exp(-1*(ln(O_matrix[p][q]+0.01)-ln(x0[p]))^2/(2*sigma[p]^2))	// overall scaling depends on nC and nO
		krxn_matrix[][] = 10^(W_coef_Carbon[0]+W_coef_Carbon[1]*C_matrix[p][q]^W_coef_Carbon[2])	// base alkane behavior, developed to account for T-dependence with Ea = 1
		krxn_matrix *= lognorm[p][q]				// adjust for oxygen addition	
//		Ea = 10 	// comes from input										// activation energy, kJ/mol.K
		multiplier = (298.15^2)*(exp(-1*1*1000/(8.314*298.15)))	// Temperature adjustment a la Arrhenius
		krxn_matrix[][] *= multiplier							// Final result
		multiplier = exp(-1*(Ea*1000/8.314)*(1/TempK-1/298.15)) // v7.4.1
		krxn_matrix[][] *= multiplier
	endif // end kinetics choice
	
// Adjust parent rate coefficient if desired
	if(krxn_parent == 0)
		// do nothing
	else
		krxn_matrix[0][0] = krxn_parent
	endif

// Adjust to deal with nitrates
	if(NO3adjustment==1)
		krxn_matrix[][3] = krxn_matrix[p][0]
	endif
	
// If O3 reaction, set up rxn matrix differently and set Pfrag to zero
	if(O3_yn==1)
		krxn_matrix[][] = 0
		krxn_matrix[0][0] = krxn_parent
		FragSlope = 0
	endif
End

//**************************************************************************************************
Function SOM_Fragmentation(Frag_Method,FragSlope,Pfrag_type,C_matrix,O_matrix,OC_matrix)
// 8. Populate fragmentation probability matrix, which is a 3D matrix with rows = #carbons, columns = #oxygens, chunks = probability
// 	Each molecule is assigned a probability of, if it falls apart, what fragment(s) will form.
// 	x = Ncarbons_product, y = Noxygens_product, z = Ncarbons_parent, t = Nox_precursor
//	SetRandomSeed(0.1) // constrains the noise generator to be "reproducable"
	// For z dimension, read across rows, then down columns
	// e.g. 00, 01, 02, 03, 04, 10, 11, 12, 13, 14, 20, 21, 22, 23
// 05/13/12 - Allow one additional oxygen atom to be added to one of the fragmented products (in Prob 2 Matrix) and assume that no products can have zero oxygens
// 08/24/13 - Added new Pfrag_type: Pfrag = (O:C)*qfrag
	variable Frag_Method
	variable FragSlope
	variable Pfrag_Type
	wave C_matrix
	wave O_matrix
	wave OC_matrix
	
	wave RefMatrix = OC_matrix
	variable Ncarbons = dimsize(RefMatrix,0)
	variable Noxygens = dimsize(RefMatrix,1)
	
	make/o/d/n=(Ncarbons,Noxygens,Ncarbons*Noxygens) Prob1_Matrix = 0
	make/o/d/n=(Ncarbons,Noxygens,Ncarbons*Noxygens) Prob2_Matrix = 0
	make/o/d/n=(Ncarbons,Noxygens) Frag1_Matrix = 0
	make/o/d/n=(Ncarbons*Noxygens) Frag1_Array = 0
	make/o/d/n=(Ncarbons,Noxygens,Ncarbons*Noxygens) Prob_Matrix = 0
	make/o/d/n=(Ncarbons,Noxygens)/FREE Prob1_Subset = 0
	variable m, n, i, k
	variable ProbIndex

   // Matrix #1 (1st fragment)
	if(Frag_Method == 3 || Frag_Method == 4) // 0 = random fragments; // 2 = weight by gaussian to produce more fragments that have 1/2 of parent
		setrandomseed/BETR 0.25
		Prob1_matrix = abs(enoise(0.25))
		for(m=1;m<=NCarbons;m+=1)
			for(n=1;n<=Noxygens;n+=1)
				ProbIndex = (m-1)*Noxygens + (n-1)
				Prob1_Matrix[][n+1,][ProbIndex] = (n < Noxygens-1 ? 0 : Prob1_Matrix[p][q][ProbIndex])//(n < Noxygens-1 ? 0 : Prob1_Matrix[p][q][ProbIndex])
				Prob1_Matrix[,m-1][][ProbIndex] = 0
				Prob1_Matrix[][][ProbIndex] = (numtype(RefMatrix[p][q]) == 2 ? 0 : Prob1_Matrix[p][q][ProbIndex])
				Prob1_Matrix[][0][ProbIndex] = 0 // 051312 --> cannot have fragments with zero oxygens
				Prob1_Subset = Prob1_Matrix[p][q][ProbIndex]
				WaveStats/Q Prob1_Subset
				Prob1_Matrix[][][ProbIndex] = (V_sum == 0 ? 0 : Prob1_Matrix[p][q][ProbIndex]/V_sum)	
			endfor
		endfor
		
		if(Frag_Method == 4) // use this to weight by gaussian to produce more fragments that have 1/2 of parent
			make/o/d/n=(Ncarbons) Weighting_wv
			variable A = 1
			variable sigma = 2*(12/Ncarbons)
			for(n=1;n<=Ncarbons;n+=1)
				sigma = 3*((Ncarbons-(n-1))/Ncarbons) // allow width to vary with Nc
				weighting_wv = A*exp(-(((Ncarbons+(n-1))/2-x)/sigma)^2) // gaussian weighting
				prob1_matrix[][][(n-1)*Noxygens,n*Noxygens-1] *= weighting_wv[p] // apply gaussian weighting
			endfor
			// renormalize
			for(m=1;m<=NCarbons;m+=1)
			for(n=1;n<=Noxygens;n+=1)
				ProbIndex = (m-1)*Noxygens + (n-1)
				Prob1_Subset = Prob1_Matrix[p][q][ProbIndex]
				WaveStats/Q Prob1_Subset
				Prob1_Matrix[][][ProbIndex] = (V_sum == 0 ? 0 : Prob1_Matrix[p][q][ProbIndex]/V_sum)	
			endfor
			endfor
		endif
	
	elseif(Frag_Method == 2) // Use this to only produce fragments with 1 carbon atom
		Prob1_Matrix[][][] = 0
		Prob1_Matrix[Ncarbons-1][0,2][] = 1/3
		Prob1_Matrix[Ncarbons-1][0,1][;Noxygens] = 0.5
		Prob1_Matrix[Ncarbons-1][2][;Noxygens] = 0
	elseif(Frag_Method == 1) // Equal probabilities
		for(m=1;m<=Ncarbons;m+=1)
			for(n=1;n<=Noxygens;n+=1)
				ProbIndex = (m-1)*Noxygens + (n-1)
				Prob1_Subset = nan
				for(i=1;i<=Ncarbons;i+=1)
					for(k=1;k<=Noxygens;k+=1)
						Prob1_Subset[i-1][k-1] = (i>m) && (k<=n+1 && k>1) && ((n-k+1)/(i-m))<2? 1 : nan
					endfor
				endfor
				Prob1_Subset = numtype(RefMatrix[p][q])==2 ? nan : Prob1_Subset
				wavestats/q Prob1_Subset
				Prob1_Subset *= 1/V_npnts
				Prob1_Subset = numtype(Prob1_Subset)==2 ? 0 : Prob1_Subset
				Prob1_Matrix[][][ProbIndex] = Prob1_Subset[p][q]				
			endfor
		endfor
	
	elseif(Frag_Method == 0) // Better Random Probabilities
		setrandomseed/BETR 0.25
		Prob1_matrix = abs(enoise(0.25))
		for(m=1;m<=Ncarbons;m+=1)
			for(n=1;n<=Noxygens;n+=1)
				ProbIndex = (m-1)*Noxygens + (n-1)
				Prob1_Subset = nan
				for(i=1;i<=Ncarbons;i+=1)
					for(k=1;k<=Noxygens;k+=1)
						Prob1_Subset[i-1][k-1] = (i>m) && (k<=n+1 && k>1) && ((n-k+1)/(i-m))<2 ? Prob1_Matrix[p][q][ProbIndex] : 0
					endfor
				endfor
				Prob1_Subset = numtype(RefMatrix[p][q])==2 ? 0 : Prob1_Subset
				wavestats/q Prob1_Subset
				Prob1_Subset = V_sum == 0 ? 0 : Prob1_Subset/V_sum
				//Prob1_Subset = numtype(Prob1_Subset)==2 ? 0 : Prob1_Subset
				Prob1_Matrix[][][ProbIndex] = Prob1_Subset[p][q]				
			endfor
		endfor
	else	// some stuff
		for(m=1;m<=Ncarbons;m+=1)
			for(n=1;n<=Noxygens;n+=1)
				ProbIndex = (m-1)*Noxygens + (n-1)
				Prob1_Subset = nan
				for(i=1;i<=Ncarbons;i+=1)
					for(k=1;k<=Noxygens;k+=1)
						Prob1_Subset[i-1][k-1] = (i>m) && (k<=n+1 && k>1) && ((n-k+1)/(i-m))<2? 1 : nan
					endfor
				endfor
				Prob1_Subset = numtype(RefMatrix[p][q])==2 ? nan : Prob1_Subset
				Prob1_Subset[][3,] = nan
				wavestats/q Prob1_Subset
				Prob1_Subset *= 1/V_npnts
				Prob1_Subset = numtype(Prob1_Subset)==2 ? 0 : Prob1_Subset
				Prob1_Matrix[][][ProbIndex] = Prob1_Subset[p][q]				
			endfor
		endfor
	endif
   // Matrix #2 (complementary fragment to have mass balance)	
	for(m=1;m<=Ncarbons;m+=1)
		for(n=1;n<=Noxygens;n+=1)
			ProbIndex = (m-1)*Noxygens + (n-1)
			for(i=1;i<=Ncarbons;i+=1)
				for(k=1;k<=Noxygens;k+=1)
					Prob2_Matrix[i-1][k][ProbIndex] = Prob1_Matrix[Ncarbons - ((i-1) - (m-1))][(n-1) - (k-2)][ProbIndex] // 051312 - changed "k-1" to "k" to effectively have 1 oxygen added to these molecules
				endfor
			endfor
			Prob2_Matrix[][n+1,][ProbIndex] = (n < Noxygens-1 ? 0 : Prob2_Matrix[p][q][ProbIndex])
			Prob2_Matrix[,m-1][][ProbIndex] = 0
			Prob2_Matrix[][][ProbIndex] = (numtype(RefMatrix[p][q]) == 2 ? 0 : Prob2_Matrix[p][q][ProbIndex])
		endfor
	endfor
  // Determine fragmentation operator
	if(Pfrag_type==0) // Pfrag = Noxygens*cfrag
		Frag1_Matrix = FragSlope *O_matrix * (C_matrix/12)^(0/3)
		Frag1_Matrix = Frag1_Matrix < 0 ? 0 : Frag1_Matrix
	elseif(Pfrag_type==1) // Pfrag = (O:C)^(mfrag)
		Frag1_Matrix = (OC_matrix)^(FragSlope)
	elseif(Pfrag_type==2) // Pfrag = (Noxygens+1)*cfrag
		Frag1_Matrix = FragSlope * (O_matrix+1)
	elseif(Pfrag_type==3) // constant
		Frag1_Matrix = FragSlope
	elseif(Pfrag_type==4) // Pfrag = (O:C)*qfrag
		Frag1_Matrix = FragSlope*(1-log(C_matrix/(O_matrix+1)))//log(O_matrix+1)//(OC_matrix)*FragSlope
		Frag1_Matrix = Frag1_Matrix > 1 ? 1 : Frag1_Matrix
		Frag1_matrix = Frag1_Matrix < 0 ? 0 : Frag1_matrix
	endif
	Frag1_Matrix = (Frag1_Matrix > 1 || numtype(Frag1_Matrix) == 2 ? 1 : Frag1_Matrix)
//	OxygensPerRxn_Matrix = 1
	for(i=0;i<Ncarbons;i+=1)
		Frag1_array[i*nOxygens,(i+1)*nOxygens-1] = Frag1_Matrix[i][p-i*nOxygens]
	endfor
End

//**************************************************************************************************
Function SOM_Fragmentation_MP(Frag_Method,FragSlope,Pfrag_type,nCmax,C_matrix,O_matrix,OC_matrix)
// 8. Populate fragmentation probability matrix, which is a 3D matrix with rows = #carbons, columns = #oxygens, chunks = probability
// 	Each molecule is assigned a probability of, if it falls apart, what fragment(s) will form.
// 	x = Ncarbons_product, y = Noxygens_product, z = Ncarbons_parent, t = Nox_precursor
//	SetRandomSeed(0.1) // constrains the noise generator to be "reproducable"
	// For z dimension, read across rows, then down columns
	// e.g. 00, 01, 02, 03, 04, 10, 11, 12, 13, 14, 20, 21, 22, 23
// 05/13/12 - Allow one additional oxygen atom to be added to one of the fragmented products (in Prob 2 Matrix) and assume that no products can have zero oxygens
// 08/24/13 - Added new Pfrag_type: Pfrag = (O:C)*qfrag
	variable Frag_Method
	variable FragSlope
	variable Pfrag_Type
	variable nCmax
	wave C_matrix
	wave O_matrix
	wave OC_matrix
	
	wave RefMatrix = OC_matrix
	variable Ncarbons = dimsize(RefMatrix,0)
	variable Noxygens = dimsize(RefMatrix,1)
	
	make/o/d/n=(Ncarbons,Noxygens,Ncarbons*Noxygens) Prob1_Matrix = 0
	make/o/d/n=(Ncarbons,Noxygens,Ncarbons*Noxygens) Prob2_Matrix = 0
	make/o/d/n=(Ncarbons,Noxygens) Frag1_Matrix = 0
	make/o/d/n=(Ncarbons*Noxygens) Frag1_Array = 0
	make/o/d/n=(Ncarbons,Noxygens,Ncarbons*Noxygens) Prob_Matrix = 0
	make/o/d/n=(Ncarbons,Noxygens)/FREE Prob1_Subset = 0
	variable m, n, i, k
	variable ProbIndex

   // Matrix #1 (1st fragment)
	if(Frag_Method == 3 || Frag_Method == 4) // 0 = random fragments; // 2 = weight by gaussian to produce more fragments that have 1/2 of parent
		setrandomseed/BETR 0.25
		Prob1_matrix = abs(enoise(0.25))
		for(m=1;m<=NCarbons;m+=1)
			for(n=1;n<=Noxygens;n+=1)
				ProbIndex = (m-1)*Noxygens + (n-1)
				Prob1_Matrix[][n+1,][ProbIndex] = (n < Noxygens-1 ? 0 : Prob1_Matrix[p][q][ProbIndex])//(n < Noxygens-1 ? 0 : Prob1_Matrix[p][q][ProbIndex])
				Prob1_Matrix[,m-1][][ProbIndex] = 0
				Prob1_Matrix[][][ProbIndex] = (numtype(RefMatrix[p][q]) == 2 ? 0 : Prob1_Matrix[p][q][ProbIndex])
				Prob1_Matrix[][0][ProbIndex] = 0 // 051312 --> cannot have fragments with zero oxygens
				Prob1_Subset = Prob1_Matrix[p][q][ProbIndex]
				WaveStats/Q Prob1_Subset
				Prob1_Matrix[][][ProbIndex] = (V_sum == 0 ? 0 : Prob1_Matrix[p][q][ProbIndex]/V_sum)	
			endfor
		endfor
		
		if(Frag_Method == 4) // use this to weight by gaussian to produce more fragments that have 1/2 of parent
			make/o/d/n=(Ncarbons) Weighting_wv
			variable A = 1
			variable sigma = 2*(12/Ncarbons)
			for(n=1;n<=Ncarbons;n+=1)
				sigma = 3*((Ncarbons-(n-1))/Ncarbons) // allow width to vary with Nc
				weighting_wv = A*exp(-(((Ncarbons+(n-1))/2-x)/sigma)^2) // gaussian weighting
				prob1_matrix[][][(n-1)*Noxygens,n*Noxygens-1] *= weighting_wv[p] // apply gaussian weighting
			endfor
			// renormalize
			for(m=1;m<=NCarbons;m+=1)
			for(n=1;n<=Noxygens;n+=1)
				ProbIndex = (m-1)*Noxygens + (n-1)
				Prob1_Subset = Prob1_Matrix[p][q][ProbIndex]
				WaveStats/Q Prob1_Subset
				Prob1_Matrix[][][ProbIndex] = (V_sum == 0 ? 0 : Prob1_Matrix[p][q][ProbIndex]/V_sum)	
			endfor
			endfor
		endif
	
	elseif(Frag_Method == 2) // Use this to only produce fragments with 1 carbon atom
		Prob1_Matrix[][][] = 0
		Prob1_Matrix[Ncarbons-1][0,2][] = 1/3
		Prob1_Matrix[Ncarbons-1][0,1][;Noxygens] = 0.5
		Prob1_Matrix[Ncarbons-1][2][;Noxygens] = 0
	elseif(Frag_Method == 1) // Equal probabilities
		for(m=1;m<=Ncarbons;m+=1)
			for(n=1;n<=Noxygens;n+=1)
				ProbIndex = (m-1)*Noxygens + (n-1)
				Prob1_Subset = nan
				for(i=1;i<=Ncarbons;i+=1)
					for(k=1;k<=Noxygens;k+=1)
						Prob1_Subset[i-1][k-1] = (i>m) && (k<=n+1 && k>1) && ((n-k+1)/(i-m))<2? 1 : nan
					endfor
				endfor
				Prob1_Subset = numtype(RefMatrix[p][q])==2 ? nan : Prob1_Subset
				wavestats/q Prob1_Subset
				Prob1_Subset *= 1/V_npnts
				Prob1_Subset = numtype(Prob1_Subset)==2 ? 0 : Prob1_Subset
				Prob1_Matrix[][][ProbIndex] = Prob1_Subset[p][q]				
			endfor
		endfor
	
	elseif(Frag_Method == 0) // Better Random Probabilities
		setrandomseed/BETR 0.25
		Prob1_matrix = abs(enoise(0.25))
		for(m=1;m<=Ncarbons;m+=1)
			for(n=1;n<=Noxygens;n+=1)
				ProbIndex = (m-1)*Noxygens + (n-1)
				Prob1_Subset = nan
				for(i=1;i<=Ncarbons;i+=1)
					for(k=1;k<=Noxygens;k+=1)
						Prob1_Subset[i-1][k-1] = (i>m) && (k<=n+1 && k>1) && ((n-k+1)/(i-m))<2 ? Prob1_Matrix[p][q][ProbIndex] : 0
					endfor
				endfor
				Prob1_Subset = numtype(RefMatrix[p][q])==2 ? 0 : Prob1_Subset
				wavestats/q Prob1_Subset
				Prob1_Subset = V_sum == 0 ? 0 : Prob1_Subset/V_sum
				//Prob1_Subset = numtype(Prob1_Subset)==2 ? 0 : Prob1_Subset
				Prob1_Matrix[][][ProbIndex] = Prob1_Subset[p][q]				
			endfor
		endfor
	else	// some stuff
		for(m=1;m<=Ncarbons;m+=1)
			for(n=1;n<=Noxygens;n+=1)
				ProbIndex = (m-1)*Noxygens + (n-1)
				Prob1_Subset = nan
				for(i=1;i<=Ncarbons;i+=1)
					for(k=1;k<=Noxygens;k+=1)
						Prob1_Subset[i-1][k-1] = (i>m) && (k<=n+1 && k>1) && ((n-k+1)/(i-m))<2? 1 : nan
					endfor
				endfor
				Prob1_Subset = numtype(RefMatrix[p][q])==2 ? nan : Prob1_Subset
				Prob1_Subset[][3,] = nan
				wavestats/q Prob1_Subset
				Prob1_Subset *= 1/V_npnts
				Prob1_Subset = numtype(Prob1_Subset)==2 ? 0 : Prob1_Subset
				Prob1_Matrix[][][ProbIndex] = Prob1_Subset[p][q]				
			endfor
		endfor
	endif
   // Matrix #2 (complementary fragment to have mass balance)	
	for(m=1;m<=Ncarbons;m+=1)
		for(n=1;n<=Noxygens;n+=1)
			ProbIndex = (m-1)*Noxygens + (n-1)
			for(i=1;i<=Ncarbons;i+=1)
				for(k=1;k<=Noxygens;k+=1)
					Prob2_Matrix[i-1][k][ProbIndex] = Prob1_Matrix[Ncarbons - ((i-1) - (m-1))][(n-1) - (k-2)][ProbIndex] // 051312 - changed "k-1" to "k" to effectively have 1 oxygen added to these molecules
				endfor
			endfor
			Prob2_Matrix[][n+1,][ProbIndex] = (n < Noxygens-1 ? 0 : Prob2_Matrix[p][q][ProbIndex])
			Prob2_Matrix[,m-1][][ProbIndex] = 0
			Prob2_Matrix[][][ProbIndex] = (numtype(RefMatrix[p][q]) == 2 ? 0 : Prob2_Matrix[p][q][ProbIndex])
		endfor
	endfor
  // Determine fragmentation operator
	if(Pfrag_type==0) // Pfrag = Noxygens*cfrag
		Frag1_Matrix = FragSlope *O_matrix * (C_matrix/12)^(0/3)
		Frag1_Matrix = Frag1_Matrix < 0 ? 0 : Frag1_Matrix
	elseif(Pfrag_type==1) // Pfrag = (O:C)^(mfrag)
		Frag1_Matrix = (OC_matrix)^(FragSlope)
	elseif(Pfrag_type==2) // Pfrag = (Noxygens+1)*cfrag
		Frag1_Matrix = FragSlope * (O_matrix+1)
	elseif(Pfrag_type==3) // constant
		Frag1_Matrix = FragSlope
	elseif(Pfrag_type==4) // Pfrag = (O:C)*qfrag
		Frag1_Matrix = FragSlope*(1-log(C_matrix/(O_matrix+1)))//log(O_matrix+1)//(OC_matrix)*FragSlope
		Frag1_Matrix = Frag1_Matrix > 1 ? 1 : Frag1_Matrix
		Frag1_matrix = Frag1_Matrix < 0 ? 0 : Frag1_matrix
	endif
	Frag1_Matrix = (Frag1_Matrix > 1 || numtype(Frag1_Matrix) == 2 ? 1 : Frag1_Matrix)
//	OxygensPerRxn_Matrix = 1

	// Rename these to deal with varying number of carbon atoms between precursors
	duplicate/o Prob1_matrix Prob1_matrix_small
	duplicate/o Prob2_matrix Prob2_matrix_small
	duplicate/o Frag1_matrix Frag1_matrix_small
//	duplicate/o Frag1_array Frag1_array_small
//	duplicate/o Prob_matrix Prob_Matrix_small
	
	make/o/d/n=(nCmax,Noxygens,nCmax*Noxygens) Prob1_Matrix = 0
	make/o/d/n=(nCmax,Noxygens,nCmax*Noxygens) Prob2_Matrix = 0
	make/o/d/n=(nCmax,Noxygens) Frag1_Matrix = 0
	make/o/d/n=(nCmax*Noxygens) Frag1_Array = 0
	make/o/d/n=(nCmax,Noxygens,nCmax*Noxygens) Prob_Matrix = 0
	Prob1_matrix[nCmax-Ncarbons,][][nCmax*Noxygens-Ncarbons*Noxygens,] = Prob1_matrix_small[p-(nCmax-Ncarbons)][q][r-(nCmax*Noxygens-Ncarbons*Noxygens)]
	Prob2_matrix[nCmax-Ncarbons,][][nCmax*Noxygens-Ncarbons*Noxygens,] = Prob2_matrix_small[p-(nCmax-Ncarbons)][q][r-(nCmax*Noxygens-Ncarbons*Noxygens)]
	Frag1_matrix[nCmax-Ncarbons,][] = Frag1_matrix_small[p-(nCmax-Ncarbons)][q]
	
//	for(i=0;i<Ncarbons;i+=1)
//		Frag1_array[i*nOxygens,(i+1)*nOxygens-1] = Frag1_Matrix[i][p-i*nOxygens]
//	endfor
	for(i=0;i<nCmax;i+=1)
		Frag1_array[i*nOxygens,(i+1)*nOxygens-1] = Frag1_Matrix[i][p-i*nOxygens]
	endfor

End

//***********************************************************************************************
Function SOM_CreateAtomMatrices(Ncarbons,Noxygens,H_per_O,Hadjustment)
// 5. Populate O:C, H:C, Oxidation State and MW Matrices
// updated on 10/24/2015 to allow for creation in specific folders, if desired
	variable Ncarbons
	variable Noxygens
	variable H_per_O
	variable Hadjustment
		
	make/o/d/n=(Ncarbons,Noxygens) C_Matrix = NaN
	make/o/d/n=(Ncarbons,Noxygens) O_Matrix = NaN
	make/o/d/n=(Ncarbons,Noxygens) H_Matrix = NaN
	make/o/d/n=(Ncarbons,Noxygens) OC_Matrix = NaN
	make/o/d/n=(Ncarbons,Noxygens) HC_Matrix = NaN
	make/o/d/n=(Ncarbons,Noxygens) Ox_Matrix = NaN
	make/o/d/n=(Ncarbons,Noxygens) MW_Matrix = NaN
	
	C_Matrix[][] = Ncarbons-x
	O_Matrix[][] = y
	   //5a. account for oxidation state difference. may end up with fractional numbers of hydrogens (added 01/16/12)
	H_Matrix[][] = (Ncarbons-x-2)*2 + 2*3 - y*H_per_O - Hadjustment
	H_Matrix[][] = (H_matrix <0 ? 0 : H_matrix)
	O_Matrix[][] = (O_Matrix > 2*C_Matrix+0 ? NaN : O_Matrix)
	C_Matrix[][] = numtype(O_Matrix)==2 ? NaN : C_Matrix
	O_Matrix[Ncarbons-1][0] = 0; O_Matrix[Ncarbons-1][1] = 1;O_Matrix[Ncarbons-1][2] = 2;O_Matrix[Ncarbons-1][3,] = nan
		
	OC_Matrix[][] = O_Matrix/C_Matrix
	HC_Matrix[][] = H_Matrix/C_Matrix
	Ox_Matrix[][] = 2*OC_Matrix - HC_Matrix
	MW_Matrix[][] = 12*C_Matrix + 16*O_Matrix + H_Matrix
	
End

//***************************************************************************************
Function SOM_Heterogeneous(index,PM_matrix,Np,Tvar,gammaOH,timestep,FragSlope,ProbOx1,ProbOx2,ProbOx3,ProbOx4,seedConcMolec)
	variable index
	wave PM_matrix	// particle molecules matrix, molecules/cm^3
	variable Np,Tvar, gammaOH, timestep, FragSlope
	variable ProbOx1,ProbOx2,ProbOx3,ProbOx4	// oxygen addition probabilities
	variable seedConcMolec	// seed concentration in molecules/cm^3
	
	variable Ncarbons = dimsize(PM_Matrix,0)
	variable nOxygens = dimsize(PM_Matrix,1)
	
	NVAR AddMultipleOxygens
	
	wave Dp_time
	wave hetchem_matrix_minus
	wave hetchem_matrix_plus
	make/o/d/n=(nCarbons,nOxygens)/FREE HetChem_Matrix_Plus_Ox = 0
	wave RxnStep_matrix_frag
	wave Prob1_matrix
	wave Prob2_matrix
	wave Frag1_matrix
	wave GasMolecules_matrix
	wave TotalMolecules_matrix
	wave OH_wave
	
	wave Prob_Matrix
	wave RxnStep_array_Minus
	wave Frag1_Array
	wave ParticleMoleFraction_Matrix
	
	variable OH_Flux
	variable OH_rate
	variable MW_OH = 17/6.022e23	// molecular weight in g/molecule
	variable Diff_OH = 0.3e-5			// m^2/s
	variable MFP_OH = 3*Diff_OH/sqrt((8*1.381e-23*Tvar)/(pi*MW_OH/1000))	// mean free path in meters, for Fuchs-Sutugin correction
	variable Kn_OH = 2*MFP_OH/(Dp_Time[index-1]*1e-9)	// knudsen number
	variable Beta_OH = (1+Kn_OH)/(1+Kn_OH+0.75/Kn_OH)	// Fuchs transition regime correction (kinetic --> transition)
	variable n, i
	variable probindex

	OH_Flux = Beta_OH*(OH_wave[index-1]*1e6)*1.381e-23*Tvar/sqrt(2*pi* (17/1000/6.022e23)*1.381e-23*Tvar) // OH flux (molecules m^-2 s^-1)
	OH_Rate = gammaOH * (pi*(dp_time[index-1]*1e-9)^2) * OH_Flux // OH collision rate (molecules/s per particle)
	wavestats/Q PM_matrix
//	ParticleMoleFraction_Matrix = PM_matrix/V_sum//ParticleMolecules_Fraction = PM_matrix/V_sum
	ParticleMoleFraction_Matrix = PM_Matrix/(V_sum+SeedConcMolec)	// changed on 8/23/13 to include seed
	// loss term
	HetChem_Matrix_minus[][] = OH_Rate*ParticleMoleFraction_Matrix[p][q] * Np * timestep * (1/1e6) // (molecules/s.particle)*(particles/m^3)*(s)*(m^3/cm^3) = (molecules/cm^3)
	HetChem_Matrix_minus[][Noxygens-1] = 0 // needed as limit on rxn.
	HetChem_Matrix_minus[Ncarbons-1][2] = 0 // ..CO2 is unreactive
	HetChem_Matrix_minus = (numtype(HetChem_Matrix_minus)==2 ? 0 : HetChem_Matrix_minus) // Turn NaN's to zero's' for later math.
	HetChem_Matrix_minus = HetChem_matrix_minus > PM_matrix ? PM_matrix : HetChem_Matrix_Minus
	// gain term; assumes reactions add 1 oxygen
	HetChem_Matrix_Plus = 0

	if(AddMultipleOxygens==0)
		// add only 1 oxygen
			HetChem_Matrix_Plus[][1,] = HetChem_Matrix_Minus[p][q-1]*(1-Frag1_Matrix[p][q-1])	// added on 07/11/13
	else
		// add oxygens based on Pfunc, updated on 07/11/13
			HetChem_Matrix_Plus_Ox = 0
			HetChem_Matrix_Plus_Ox[][1,] = HetChem_Matrix_Minus[p][q-1]*ProbOx1*(1-Frag1_Matrix[p][q-1])
			HetChem_Matrix_Plus += HetChem_Matrix_Plus_Ox
			HetChem_Matrix_Plus_Ox[][2,] = HetChem_Matrix_Minus[p][q-2]*ProbOx1*(1-Frag1_Matrix[p][q-2])
			HetChem_Matrix_Plus_Ox[][Noxygens-1] += (HetChem_Matrix_Minus[p][Noxygens-2]*(1-Frag1_Matrix[p][Noxygens-2]))*ProbOx2
			HetChem_Matrix_Plus += HetChem_Matrix_Plus_Ox
			HetChem_Matrix_Plus_Ox[][3,] = HetChem_Matrix_Minus[p][q-3]*ProbOx1*(1-Frag1_Matrix[p][q-3])
			HetChem_Matrix_Plus_Ox[][Noxygens-1] += (HetChem_Matrix_Minus[p][Noxygens-2]*(1-Frag1_Matrix[p][Noxygens-2])+HetChem_Matrix_Minus[p][Noxygens-3]*(1-Frag1_Matrix[p][Noxygens-3]))*ProbOx3
			HetChem_Matrix_Plus += HetChem_Matrix_Plus_Ox
			HetChem_Matrix_Plus_Ox[][4,] = HetChem_Matrix_Minus[p][q-4]*ProbOx1*(1-Frag1_Matrix[p][q-4])
			HetChem_Matrix_Plus_Ox[][Noxygens-1] += (HetChem_Matrix_Minus[p][Noxygens-2]*(1-Frag1_Matrix[p][Noxygens-2])+HetChem_Matrix_Minus[p][Noxygens-3]*(1-Frag1_Matrix[p][Noxygens-3])+HetChem_Matrix_Minus[p][Noxygens-4]*(1-Frag1_Matrix[p][Noxygens-4]))*ProbOx4
			HetChem_Matrix_Plus += HetChem_Matrix_Plus_Ox
	endif
	HetChem_Matrix_plus = (numtype(HetChem_Matrix_plus)==2 ? 0 : HetChem_Matrix_plus)
	if(index==1)
		HetChem_Matrix_minus = 0
		HetChem_Matrix_plus = 0
	endif
	// Account for Loss 
	PM_matrix[][] -= HetChem_Matrix_minus[p][q] // account for loss
	
	// Account for formation & fragmentation: add all molecules directly to the gas-phase
	IF(FragSlope != 0) // only run this if fragmentation indeed occurs.
		RxnStep_Matrix_Frag = 0
		for(i=0;i<Ncarbons;i+=1)
			RxnStep_array_Minus[i*nOxygens,(i+1)*nOxygens-1] = HetChem_Matrix_Minus[i][p-i*nOxygens]
		endfor
		multithread prob_matrix = (Prob1_Matrix+Prob2_Matrix)*Frag1_Array[r]*RxnStep_Array_Minus[r]
		MatrixOp/O RxnStep_Matrix_Frag = sumBeams(prob_matrix)
			PM_Matrix[][] += RxnStep_Matrix_Frag + HetChem_Matrix_Plus	// modified on 07/11/13
	ELSE // No fragmentation
		PM_matrix[][] += HetChem_Matrix_plus // molecules formed per step, accounting for different numbers of oxygens added per step
	ENDIF
	
	TotalMolecules_Matrix = GasMolecules_Matrix + PM_Matrix
END

//*******************************************************************************************
Function SOM_Oligomerization(PM,krxn,timestep,Np,VolumeOrgPerParticle,[FirstTime,method])
	// this version does not treat oxygens correctly
	wave PM	// particleMolecules_Matrix [molecules/cm^3]
	wave krxn	// krxn matrix [cm^3/molecules.s] (2D matrix, currently, for simplicity)
	variable timestep	// timestep [seconds]
	variable Np // particle number concentration [particles/cm^3]
	variable VolumeOrgPerParticle // per particle volume [m^3/p]
	variable FirstTime	// 0 = no, 1 = yes, is this the first time the calculatino is being performed?
	variable method	// which method to use?
	
	if(ParamIsDefault(FirstTime))
		FirstTime = 1
	endif
	if(ParamisDefault(method))
		method = 0
	endif
	
	variable nC = dimsize(PM,0)
	variable nO = dimsize(PM,1)
	variable nC_parent = nC/2
	variable nO_parent = nO/2
	variable TotalMoleculesPerParticle
	variable CurrentConc
	variable Gain
	variable i, j, k, m
	
	// Create waves, if first time. Otherwise, they should exist.
	if(FirstTime==1)
		make/o/d/n=(nC_parent,nO_parent) PM_subset = 0
		make/o/d/n=(nC_parent,nO_parent) PM_Loss = 0
		make/o/d/n=(nC_parent,nO_parent) PM_Gain = 0
		make/o/d/n=(nC,nO) PM_Loss_Full = 0
		make/o/d/n=(nC,nO) PM_Gain_Full = 0
		make/o/d/n=(nO_parent) ConcByOx1, ConcByOx2
		make/o/d/n=(nC,nO) PM_BigMatrix = 0
	endif

	variable PreCalcFactor	
	variable sumO
	variable sumC
	
	// Loss from reaction
	variable krxn_temp = krxn[nC_parent][0]
	PM_Loss = 0
	PM_Loss_Full = 0
	PM_subset[][] = PM[p+nC_parent][q]/(Np*VolumeOrgPerParticle*1e6)	// molecules/cm^3 --> molecules/cm^3 of particle
	PM_subset = numtype(PM_Subset)==2 ? 0 : PM_Subset
	wavestats/q PM_subset
	TotalMoleculesPerParticle = V_sum

	if(method==0)
		PM_Loss = krxn_temp*PM_Subset*(TotalMoleculesPerParticle)*timestep*(Np*VolumeOrgPerParticle*1e6)	// molecules lost per cm^3
		PM_Loss = PM[p+nC_parent][q]<PM_Loss ? PM[p+nC_parent][q] : PM_Loss
		PM_Gain = 0
		matrixop/o ConcSum = sumRows(PM_Subset)
		for(i=0;i<nC_parent;i+=1)
			for(j=0;j<nC_parent;j+=1)
				PM_Gain[j][3] += 0.5*krxn_temp*ConcSum[i]*ConcSum[j]*timestep*(Np*VolumeOrgPerParticle*1e6)
			endfor
		endfor
		PM[0,nC_parent-1][] += PM_Gain[p][q]
		PM[nC_parent,nC][] -= PM_Loss[p-nC_parent][q]
	elseif(method==1)
		PM_Loss_Full = krxn*PM*(TotalMoleculesPerParticle)*timestep
		PM_Loss_Full = PM<PM_Loss_Full ? PM : PM_Loss_Full
		variable current
		PreCalcFactor = 0.5^2*timestep/(Np*VolumeOrgPerParticle*1e6)
		PM_Gain_Full = 0
		for(i=nC_parent;i<nC;i+=1)
			for(j=0;j<nO;j+=1)
				for(k=nC_parent;k<nC;k+=1)
					for(m=0;m<nO;m+=1)
						sumO = min(j+m,nO)
						sumC = nC-((nC-i)+(nC-k))
						current = (krxn[i][j]+krxn[k][m])*PM[i][j]*PM[k][m]*PreCalcFactor
						if(numtype(current)!=2)
							PM_Gain_Full[sumC][sumO] += current
						endif
					endfor
				endfor
			endfor
		endfor
		PM += PM_Gain_Full - PM_Loss_Full
	elseif(method==2)
		PM_Loss_Full = krxn*PM*(TotalMoleculesPerParticle)*timestep
		PM_Loss_Full = PM<PM_Loss_Full ? PM : PM_Loss_Full
		PM_BigMatrix = 0
		PreCalcFactor = 0.5*timestep*(Np*VolumeOrgPerParticle*1e6)
		for(i=0;i<nC_parent;i+=1)
			for(j=0;j<nO_parent;j+=1)
				for(k=0;k<nC_parent;k+=1)
					for(m=0;m<nO_parent;m+=1)
						//sumC = (i+k)//nC_oligs - ((nC_noOligs-i) + (nC_noOligs-k)) // INDEX of carbon, not actual number
						sumC = (i+k)//nC_oligs - ((nC_noOligs-i) + (nC_noOligs-k)) // INDEX of carbon, not actual number
						sumO = j+m	// INDEX and number of oxygens
						if(sumC<nC_parent)
							PM_BigMatrix[sumC][sumO] += PreCalcFactor*((krxn[i+nC_parent][j]+krxn[k+nC_parent][m])/2)*PM_Subset[i][j]*PM_Subset[k][m]// PreCalcFactor*1e-20*PM_noOligs_noNaNs[i][j]*PM_noOligs_noNaNs[k][m]
						endif
					endfor
				endfor
			endfor
		endfor
		PM[][] += PM_BigMatrix[p][q] - PM_Loss_Full[p][q]
	endif
End

//*******************************************************************************************
Function SOM_OligomerRateCoef(nC,nO,[krxn_base,OligsAreStable])
	variable nC	// # of carbon atoms total (usually 2 x parent carbon number)
	variable nO	// max # of oxygen atoms
	variable krxn_base	// base reaction coefficient (cm^3/molecules.s)
	variable OligsAreStable	// condition that will either allow oligomers to react or not. default is no reaction
	
	if(ParamIsDefault(krxn_base))
		krxn_base = 1e-26// cm^3/molecules.s
	endif
	if(ParamIsDefault(OligsAreStable))
		OligsAreStable = 1
	endif
		
	make/o/d/n=(nC,nO) krxn_matrix_olig
	note krxn_matrix_olig "Rate coefficients for oligomerization in particle phase [cm^3/molecules.s]"
	
	krxn_matrix_olig = krxn_base//*(y+1)
	if(OligsAreStable==1)
		krxn_matrix_olig[0,Nc/2-1][] = 0
	endif
	// don't allow CO2 to react
//	krxn_matrix_olig[dimsize(krxn_matrix_olig,0)-2,dimsize(krxn_matrix_olig,0)-1][] = 0
End

//**************************************************************************************
Function SOM_GasPhaseWallLoss(logCstar_matrix,O_matrix,kwg_on,[Cw_base,kwg_increase_perO,Krechmer_Cw])
	// Based on Matsunaga and Ziemann, or with update from Krechmer (EST, 2016)
	// assume that effective wall volume concentration is related to No
	wave logCstar_matrix	// log of C*, with C* in ug/m^3
	wave O_matrix	// matrix of oxygen atoms
	variable kwg_on	// gas-phase wall-loss rate for parent hydrocarbon in 1/s
	variable Cw_base // wall "concentration" in ug/m^3...value for alkanes from M&Z
	variable kwg_increase_perO	// 0 if kwg_on is constant, some other number if kwg_on varies with oxygen content
	variable Krechmer_Cw // allow Cw to vary with C* according to Krechmer et al., 2016
	
	if(ParamIsDefault(Cw_base))
		Cw_base = 17e3
	else
		Cw_base *= 1e3
	endif
	if(ParamIsDefault(kwg_increase_perO))
		kwg_increase_perO = 0
	endif
	if(ParamIsDefault(Krechmer_Cw))
		Krechmer_Cw = 0
	endif
	
	variable deltaCw_perO = 0e3	// by what factor does Cw increase upon addition of an oxygen
	variable Cw_limit = 20000e3	// max wall "concentration" in ug/m^3
		
	variable nC = dimsize(logCstar_matrix,0)
	variable nO = dimsize(logCstar_matrix,1)
	
	make/o/d/n=(nC,nO) Cw_wave = nan
	note/k Cw_wave "Effective wall concentration, in ug/m3"
	make/o/d/n=(nC,nO)/FREE kwg_on_matrix = kwg_on
	kwg_on_matrix = kwg_on + kwg_on * kwg_increase_perO * (O_matrix[p][q])
	kwg_on_matrix = kwg_on_matrix < 0 ? 0 : kwg_on_matrix
	make/o/d/n=(nC,nO) kwg_off = nan
	Cw_wave = Cw_base+deltaCw_perO*(O_matrix)
	Cw_wave = Cw_wave > Cw_limit ? Cw_limit : Cw_wave
	
	if(Krechmer_Cw != 0)
		// Krechmer relationship tied to SIMPOL estimated C* values; they note that these may be off by an order of magnitude
		// Allow for adjustment of initial relationship by multiplying by Cw_base as a scaling factor
		Cw_wave = 16*(Cw_base/1000)*(10^logCstar_matrix)^0.6
		Cw_wave = Cw_wave > 10000 ? 10000 : Cw_wave
		Cw_wave = Cw_wave < 16*(Cw_base/1000) ? 16*(Cw_base/1000) : Cw_wave
//		Cw_wave = 160*(10^logCstar_matrix)^0.6
//		Cw_wave = Cw_wave > 10000 ? 10000 : Cw_wave
//		Cw_wave = Cw_wave < 160 ? 160 : Cw_wave
	endif
	
//	kwg_off = kwg_on*10^(logCstar_matrix)/Cw_wave
	kwg_off = kwg_on_matrix*10^(logCstar_matrix)/Cw_wave

//	killwaves/z Cw_wave
End

//*********************************************************************
Function SOM_SetOHconc(timewv,OHconc_0,OH_scale,idex)
	wave timewv
	variable OHconc_0, OH_scale, idex
	
	variable OHconc_t
	
		OHconc_t = OHconc_0*exp(-timewv[idex-1]*oh_scale)
		return OHconc_t
end

//****************************************************************
Function SOM_setO3conc(wv,timewv,Temperature,idex,O3_conc)
	wave wv,timewv
	variable Temperature,idex,O3_conc
	
	NVAR usescalingforoxidant
	NVAR OH_scale
	
	if(usescalingforoxidant==0)
		wv[idex-1] = O3_conc*1e-9 * (1e-6*760*133.322/(1.381e-23*Temperature))			
	else // set function here
		if(timewv[idex-1]<10)
			wv[idex-1] += (O3_conc*1e-9 * (1e-6*760*133.322/(1.381e-23*Temperature)))  * exp(-timewv[idex-1]*oh_scale)//(1-exp(-3*timewv[idex-1]/OH_scale))
			wv[idex-1] += (O3_conc*1e-9 * (1e-6*760*133.322/(1.381e-23*Temperature)))  * (1-exp(-3*timewv[idex-1]/OH_scale))
		else
			wv[idex-1] = wv[idex-2]*0.985
		endif
	endif
	wv[idex-1] = wv[idex-1] < 0 ? 0 : wv[idex-1]
End

//**************************************************************************************
Function SOM_deltaHvap([logCstar_Matrix,ConstantDHvap])
	// from 09/03/13, to calculate matrix of vaporization enthalpies for vapor pressure adjustment
	wave logCstar_Matrix	// matrix of logCstar values (in log(ug/m^3))
	variable ConstantDHvap	// 0 if variable, >0 if constant value
	
	if(paramisdefault(logCstar_matrix))
		wave logCstar_matrix
	endif
	if(paramisdefault(constantDHvap))
		constantDHvap = 0
	endif
	
	variable nRows = dimsize(logCstar_matrix,0)
	variable nColumns = dimsize(logCstar_matrix,1)
	if(numtype(nRows)==2)
		abort "your logCstar_matrix wave does not exist. Please try again."
	endif
	make/o/d/n=(nRows,nColumns) deltaHvap_matrix
	
	if(constantDHvap==0)
		// use relationship from Epstein et al. (2009), modified by Cappa et al. (2010)
		deltaHvap_matrix = 131 - 11*logCstar_matrix
		deltaHvap_matrix = deltaHvap_matrix > 200 ? 200 : deltaHvap_matrix
		deltaHvap_matrix = deltaHvap_matrix < 30 ? 30 : deltaHvap_matrix
	endif
End

//**************************************************************************************
// This kills a lot of waves generated when running SOM
// Comment things out (do not delete) if you want particular waves to remain for diagnostic purposes
// should be alphabetical
Function SOM_KillWaves()

	killwaves/z beta_matrix
	KillWaves/z C_matrix
	killwaves/z ChiSqTest
	killwaves/z Chi2wave
	killwaves/z Coa
	KillWaves/z concwave
	killwaves/z ConcByOx1
	killwaves/z ConcByOx2
	killwaves/z Cshort
	killwaves/z Cstar_matrix
	KillWaves/z Ctot_init_wave
	Killwaves/z Cw_wave // 7.4.3
//	Killwaves/z deltaHC_time
	killwaves/z deltadynamicpartitioning // 7.4.3
	killwaves/z deltadynamicpartitioningpre // 7.4.3
	killwaves/z deltaHC_time
	killwaves/z deltaHC_time_interp // 7.4.3
	Killwaves/z deltaHCfrac_time
	KillWaves/z deltaHvap_matrix
	KillWaves/z deltaOx_matrix
	killwaves/z deltaVOC_ppb
	killwaves/z deltaVOC_ug
	killwaves/z dDpdt_off
	killwaves/z dDpdt_on
	killwaves/z diameter
	killwaves/z diameter_time
	killwaves/z Diffusivity_matrix
	killwaves/z dlogDp_time
//	killwaves/z dNdlogDp
//	killwaves/z dNdlogDp_time
	killwaves/z dNdlogDp_INT
	killwaves/z dSdlogDp
	killwaves/z dSdlogDp_INT
	killwaves/z dVdlogDp
	killwaves/z dVdlogDp_time
	Killwaves/z FitError_matrix // 7.4.3
	killwaves/z FitQuitReason_matrix // 7.4.3
	KillWaves/z FracAero
	KillWaves/z Frag1_Array
	KillWaves/z Frag1_Matrix
	killwaves/z Fuchs_3D // 7.4.3
	killwaves/z Gas // 7.4.3
	KillWaves/z GasFrac_time_log
	killwaves/z GasMass_matrix
//	killwaves/z GasMass_time
//	killwaves/z GasMolecules_parent
//	killwaves/z GasMolecules_time
	killwaves/z GasMolecules_matrix
	killwaves/z GasMolecules_matrix_current // 7.4.3
	killwaves/z GasMolecules_matrix_current2 // 7.4.3
	killwaves/z GasMolecules_matrix_init // 7.4.3
	killwaves/z GasMolecules_matrix_log
	killwaves/z GasMolecules_Parent // 7.4.3
	KillWaves/z H_Matrix
	KillWaves/z H_matrix_Temp
	killwaves/z HC_matrix
//	killwaves/z HC_ppb_time
//	killwaves/z H2C_time
	KillWaves/Z HetChem_Matrix_Minus
	KillWaves/Z HetChem_Matrix_Plus
	Killwaves/z HOM_time
	killwaves/z InitialGuess_matrix
	KillWaves/Z INT_result
	KillWaves/Z INT_results
	killwaves/z Kelvin // 7.4.3.
	killwaves/z Knudsen_3D // 7.4.3
	killwaves/z Knudsen_Matrix
	killwaves/Z krxn_matrix
	killwaves/z krxn_matrix_olig
	killwaves/z kwall_wv
	killwaves/z kwall_wv_im
	killwaves/z kwg_off
//	killwaves/z lifetimes
	KillWaves/Z logCstar_Matrix
	killwaves/z logCstar_matrix_ref // 7.4.3
	killwaves/z logCstar_mean_time
	killwaves/z logDp
	killwaves/z logDp_time
	killwaves/z MassSeedPerBin // 7.4.3
	killwaves/z MassSeedPerBin_noseed // 7.4.3
	killwaves/z meanfreepath_matrix
	killwaves/z MoleculesSeedPerBin
	killwaves/z Molecules_time
	killwaves/z MoleFraction2D // 7.4.3
	killwaves/z MoleFraction3D // 7.4.3
	killwaves/z MW_matrix
	killwaves/z nC // 7.4.3
	killwaves/z Ncarbons_wave
	killwaves/z Nc_ave_time
	killwaves/z Noxygens_wave
	killwaves/z No_ave_time
//	killwaves/z Np_time
	killwaves/z NpPerBin_time
	Killwaves/z NucleationRate_Time
	killwaves/z NumParticlesByBin // 7.4.3
	killwaves/z O2C_time_interp // 7.4.3
	KillWaves/z O_Matrix
	killwaves/z OC_matrix
	killwaves/z OCseed_time
	killwaves/z OH_exp_interp
	killwaves/z OH_short
	killwaves/z OH_time_experiment
	killwaves/z OH_wave_fixed
	killwaves/z Ox_matrix
	Killwaves/z OxygensPerRxn_matrix
	KillWaves/Z Particle_div_Total
	killwaves/z Particle_Fraction
	killwaves/z Particle_Fraction_log
	killWaves/z particle_temp
	killwaves/z ParticleFrac_matrix
	killwaves/z ParticleMass_matrix
//	killwaves/z ParticleMass_time
//	killwaves/z ParticleMolecules_time
	killwaves/z ParticleMolecules_matrix
	killwaves/z ParticleMolecules_matrix_cur // 7.4.3
	KillWaves/Z ParticleMoleFraction_matrix
	killwaves/z ParticleMoleFraction_time
	KillWaves/Z ParticleVolume_matrix
	killwaves/z ParticleWLR_wave
	killwaves/z Phi_GPP
	killwaves/z PM_gain
	killwaves/z PM_gain_full
	killwaves/z PM_loss
	killwaves/z PM_loss_full
	killwaves/z PM_subset
	killwaves/z PM3D_matrix
	killwaves/z PM3D_Matrix_current // 7.4.3
	killwaves/z PM3D_matrix_init // 7.4.3
	killwaves/z PMM3D_matrix
	killwaves/z PMF3D_matrix
	killwaves/z PM_BigMatrix // 7.4.3
	killwaves/z PM_matrix_init // 7.4.3
	killwaves/z Prob_matrix
	KillWaves/z Prob1_array
	KillWaves/z Prob1_Matrix
	killwaves/z Prob1_Matrix_tot
	KillWaves/z Prob2_array
	KillWaves/z Prob2_matrix
	Killwaves/z Prob2_matrix_tot
	KillWaves/z Prob2_temp
	KillWaves/z ProbOx
	Killwaves/z ProbTot_matrix
	KillWaves/Z RxnStep_Matrix_frag
	KillWaves/z RxnStep_array_minus
	KillWaves/z RxnStep_Matrix_Kinetic_Time
	KillWaves/Z RxnStep_Matrix_minus
	KillWaves/Z RxnStep_Matrix_plus
	KillWaves/Z RxnStep_Matrix_plus1
	KillWaves/Z RxnStep_Matrix_plus2
	KillWaves/Z RxnStep_Matrix_plus3
	KillWaves/Z RxnStep_Matrix_plus4
	KillWaves/Z RxnStepKinetic3D_cond
	KillWaves/Z RxnStepKinetic3D_Evap
	KillWaves/Z RxnStepKinetic3D_Pre1
	KillWaves/Z RxnStepKinetic3D_Pre2
	KillWaves/Z SatConc_matrix
	killwaves/z SatConc_Matrix_3D // 7.4.3
	killwaves/z SatConc_matrix_adj // 7.4.3
	killwaves/z SatConc_Mass
	killwaves/z SatConc_temp
	killwaves/z SatRatio_2D // 7.4.3
	killwaves/z SatRatio_3D // 7.4.3
	KillWaves/Z SeedConc_time
	KillWaves/z SummedMatrix
	killwaves/z SumDeltaDP
	killwaves/z SumGas
	killwaves/z SumKinetics // 7.4.3
	killwaves/z SumParticles // 7.4.3
	killwaves/z SumTotal // 7.4.3
	killwaves/z tarray // 7.4.3
	killwaves/z temperature_time // 7.4.3
	killwaves/z TempWave
	KillWaves/z Testing
	KillWaves/z TestW
	killwaves/z thetestwave
	killwaves/z TimeSecs // 7.4.3
	killwaves/z TimeStep_wave // 7.4.3
	killwaves/z TimeW_secs
	killwaves/z TimeW_minutes
	killwaves/z tmatrix_2D // 7.4.3
	killwaves/z tmatrix_3D // 7.4.3
	killwaves/z tot_temp
	killwaves/z TotalCarbonAtoms
	killwaves/z TotalHydrogenAtoms
	killwaves/z TotalOxygenAtoms
	KillWaves/Z TotalMass_Matrix
	killwaves/z TotalMass_time
	killwaves/z TotalMolecules_matrix
	killwaves/z totalmolecules_matrix_current // 7.4.3
	killwaves/z totalmolecules_matrix_init // 7.4.3
	KillWaves/Z TotalMolecules_Matrix_Log
	Killwaves/Z TotalMolecules_Time
	killwaves/z TotalWallMolecules_current
	killwaves/z VolumeOrgPerBin
	killwaves/z VOlumeSeedPerBin
	killwaves/z WallMass_matrix // 7.4.3
	killwaves/z WallMolecules_matrix
	killwaves/z WallMolecules_time
	killwaves/z WallParticleMolec_Matrix
	killwaves/z WPM3D_matrix
	killwaves/z Yield2_time
End

/////////////////////////////////////////////////////////////////////////////////////
Function SetOH_from_Experiment(krxn)
	variable krxn

	NVAR Nsteps
	NVAR timestep
	wave VOC_experiment_ppb
	wave Time_experiment_VOC
	variable npnts_VOC = numpnts(VOC_experiment_ppb)
	
	make/o/d/n=(Nsteps) OH_time_experiment = 0
	make/o/d/n=(npnts_VOC) OH_short
	make/o/d/n=(npnts_VOC) idexwave
	
	variable VOCfrac_exp
	variable VOCfrac_calc
	variable Ratio
	variable OH_var = 0
	variable OH_step = 0.5e5
	variable deltat
	variable m, n
	variable quitvar
	variable t_counter, counter, last
	
	for(m=1;m<=npnts_VOC-1;m+=1)
		OH_var = 0
		quitvar = 0
		VOCfrac_exp = VOC_experiment_ppb[m]/VOC_experiment_ppb[m-1]
		deltat = (Time_experiment_VOC[m] - Time_experiment_voc[m-1])*60*60 // convert hours to seconds		
		do		
			VOCfrac_calc = exp(-1*krxn*OH_var*deltat) 
			OH_var += OH_step		
		if(VOCfrac_exp > 1 || OH_var > 1e8  || VOCfrac_calc < VOCfrac_exp)
			quitvar = 1
		endif	
		while(quitvar == 0)	
		OH_short[m-1] = OH_var		
	endfor
	
		wavestats/Q OH_short
		OH_short[npnts_VOC-1] = OH_short[npnts_VOC-2]
	
	killwaves/Z 

		t_counter = 0
		counter = 0	
		quitvar = 0
		last = 0
	for(m=1;m<=npnts_VOC;m+=1)
			quitvar = 0
		do
			t_counter += timestep/60/60
			counter +=1
			if(t_counter > Time_experiment_VOC[m])
				quitvar = 1
			endif
		while(quitvar == 0)
		OH_time_experiment[last,counter] = OH_short[m-1]
		last = counter	
	endfor
	
	wavestats/q OH_time_experiment
	OH_time_Experiment[last,V_npnts] = OH_short[npnts_VOC]	
	killwaves OH_short, idexwave
End

//************************************************************
function CopySOMresults2folder(SaveFolderStr,[Ncarbons])
	string SaveFolderStr
	variable Ncarbons
	
	string namestr1, namestr2
	//wave Ncarbons
	NameStr1 = "root:"+SaveFolderStr
	print "Data saved in " + namestr1
	if(datafolderexists(NameStr1) == 0)
		NewDataFolder $NameStr1
	endif
	NameStr2 = NameStr1 + ":Yield_time" //+ num2str(Ncarbons) //+ "_" + num2str(FragSlope*100)
	wave Yield_time
	duplicate/O Yield_time $NameStr2
	NameStr2 = NameStr1 + ":Coa_time" //+ num2str(Ncarbons) //+ "_" + num2str(FragSlope*100)
	wave Coa_Time
	duplicate/O Coa_time $NameStr2
	NameStr2 = NameStr1 + ":O2C_time"// + num2str(Ncarbons) //+ "_" + num2str(FragSlope*100)
	wave O2C_time
	duplicate/O O2C_time $NameStr2
	NameStr2 = NameStr1 + ":deltaHC_time"// + num2str(Ncarbons) //+ "_" + num2str(FragSlope*100)
	wave deltaHC_time
	duplicate/O deltaHC_time $NameStr2
	NameStr2 = NameStr1 + ":deltaHCfrac_time"// + num2str(Ncarbons) //+ "_" + num2str(FragSlope*100)
	wave HC_ppb_time
	duplicate/O HC_ppb_time $NameStr2
	NameStr2 = NameStr1 + ":RunInfo" //+ num2str(Ncarbons)
	wave RunInfo
	duplicate/O RunInfo $NameStr2
	NameStr2 = NameStr1 + ":RunValues"// + num2str(Ncarbons)
	wave RuNValues
	duplicate/O RunValues $NameStr2
	NameStr2 = NameStr1 + ":ParticleMass_time"// + num2str(Ncarbons)
	wave ParticleMass_time
	duplicate/O ParticleMass_time $NameStr2
	NameStr2 = NameStr1 + ":GasMass_time"//
	wave GasMass_time
	duplicate/O GasMass_time $NameStr2
	NameStr2 = NameStr1 + ":TimeW"
	wave TimeW
	duplicate/o TimeW $NameStr2
End

//************************************************************************************
Function SaveStuff()
	newpath/O/Q pathname, "C:\Users\Chris\Documents\Research\SOA\SOM\Reactions\savefolder"
	wave Coa_time
	wave GasMolecules_Time
	wave gasmolecules_time_norm
	wave OC_time
	wave H2C_time
	wave HC_ppb_time
	wave Ncarbons_Wave
	wave Noxygens_Wave
	wave ParticleMolecules_Time
	wave particlemolecules_time_norm
	wave TimeW
	wave TimeW_minutes
	wave Yield_time
	Save/C/P=pathname/O Coa_time,GasMolecules_Time,gasmolecules_time_norm,OC_time,H2C_time,HC_ppb_time,Ncarbons_Wave,Noxygens_Wave,ParticleMolecules_Time,particlemolecules_time_norm,TimeW,TimeW_minutes,Yield_time
End

//*******************************************************
Function Run_SOM(ctrlName) : ButtonControl
	String ctrlName
	
	som_v1()

End

//*******************************************************
Function Run_SOM_Exp(ctrlName) : ButtonControl
	String ctrlName
	
	wave Time_Experiment = root:Time_Experiment
	som_v1(fittime=Time_Experiment)

End

//*******************************************************
Function ButtonProc_1(ctrlName) : ButtonControl
	String ctrlName
	
	SaveDataToFolder()
End



//***********************************************************
Function SaveDataToFolder([dfSave])
	string dfSave	// default = use panel value
	
	setdatafolder root:
	
	if(ParamIsDefault(dfSave))
		SVAR SaveFolderStr	// use folder from panel
		dfSave = SaveFolderStr
	endif

	string dffpBase = "root:SavedModelRuns"
	string dffpSave = dffpBase +":"+dfSave
	
	print "Data saved in " + dffpSave
	if(datafolderexists(dffpBase) == 0)
		NewDataFolder $dffpBase
	endif
	if(datafolderexists(dffpSave)==0)
		NewDataFolder $dffpSave
	endif
	
	string Waves2Save = "RunInfo;RunValues;TimeW;Coa_time;O2C_time;H2C_time;HC_ppb_time;OH_wave;Dp_time;"
	Waves2Save += "ParticleMolecules_time;GasMolecules_time;Yield_time;deltaHC_time;deltaHCfrac_time;"
	Waves2Save += "OH_exposure;Lifetimes"
	
	variable i
	string OldWaveStr, NewWaveStr
	variable nWaves = itemsinlist(waves2save,";")
	Make/o/t/n=1 $(dffpSave+":RunTime")
	wave/T RunTime = $(dffpSave+":RunTime")
	RunTime[0] = "Run performed at " + time() + " on " + date()
	for(i=0;i<nWaves;i+=1)
		OldWaveStr = "root:"+stringfromlist(i,waves2save,";")
		NewWaveStr = dffpSave + ":"+ stringfromlist(i,waves2save,";")
		duplicate/O $OldWaveStr, $NewWaveStr
	endfor	
end			

//***********************************************************************************************************
Function Graph_SpeciesByCarbon(wvstr,nC,[NewGraph])
	string wvstr	// name of 3D wave
	variable nC	// carbon number to consider
	variable NewGraph	// would you like to make a new graph?
	
	wave wv = $wvstr
	
	if(ParamIsDefault(NewGraph))
		NewGraph = 1
	endif
	
	variable nCarbons = dimsize(wv,0)
	variable nOxygens = dimsize(wv,1)
	variable nCindex = nCarbons-nC-1

	Make/O/N=26/FREE pc_red = {  65280, 50000, 39168, 52224, 65280, 0, 16384, 0, 0, 0, 0, 0, 0, 32768, 39168, 30464, 21760, 0, 34816, 43520, 52224, 39168, 29440, 65280, 26112, 0}
	Make/O/N=26/FREE pc_green={   0, 17480, 0, 0, 32512, 0, 48896, 17408, 0, 65280, 52224, 26112, 65280, 65280, 39168, 30464, 21760, 0, 34816, 43520, 34816, 0, 0, 43520, 17408, 0}
	Make/O/N=26/FREE pc_blue={    17408, 0, 15616, 20736, 16384, 65280, 65280, 26112, 52224, 65280, 0, 0, 0, 65280, 0, 30464, 21760, 0, 34816, 43520, 0, 31232, 58880, 0, 0, 0}

	pc_red = {0, 25535, 2, 0, 39321, 48059, 65535, 0, 16385, 65535, 0, 65535, 2, 0, 39321, 48059, 65535, 0, 16385, 65535}
	pc_green = {0, 16385, 39321, 0, 1, 48059, 32768, 65535, 65535, 32768, 0, 16385, 39321, 0, 1, 48059, 32768, 65535, 65535, 32768}
	pc_blue = {0, 16385, 1, 65535, 31457, 48059, 32768, 0, 65535, 58981, 0, 16385, 1, 65535, 31457, 48059, 32768, 0, 65535, 58981}
	
	if(NewGraph==1)
		display
	endif
	variable i
	string this_trace
	variable bottomVal, topVal
	for(i=0;i<nOxygens;i+=1)
		appendtograph wv[nCindex][i][]
		if(i==0 && NewGraph==1)
			this_trace = wvstr
		else
			this_trace = wvstr+"#"+num2str(i)
			ModifyGraph rgb($this_trace)=(pc_red[i],pc_green[i],pc_blue[i])
		endif
	endfor
	ModifyGraph lsize=2	
End

//***********************************************************************************************************
Function Graph_AllSpeciesByCarbon(wvstr,[NewGraph,minCarbons,maxOxygens])
	string wvstr
	variable NewGraph
	variable minCarbons
	variable maxOxygens
	
	wave wv = $wvstr
	
	if(ParamIsDefault(NewGraph))
		NewGraph = 1
	endif
	
	variable nCarbons = dimsize(wv,0)
	variable nOxygens = dimsize(wv,1)
	
	if(ParamIsDefault(minCarbons))
		nCarbons = nCarbons
	else
		nCarbons = nCarbons-minCarbons
	endif
	if(ParamIsDefault(MaxOxygens))
		nOxygens = nOxygens
	else
		nOxygens = MaxOxygens
	endif

	Make/O/N=26/FREE pc_red = {  65280, 52224, 39168, 52224, 65280, 0, 16384, 0, 0, 0, 0, 0, 0, 32768, 39168, 30464, 21760, 0, 34816, 43520, 52224, 39168, 29440, 65280, 26112, 0}
	Make/O/N=26/FREE pc_green={   0, 17480, 0, 0, 32512, 0, 48896, 17408, 0, 65280, 52224, 26112, 65280, 65280, 39168, 30464, 21760, 0, 34816, 43520, 34816, 0, 0, 43520, 17408, 0}
	Make/O/N=26/FREE pc_blue={    0, 0, 15616, 20736, 16384, 65280, 65280, 26112, 52224, 65280, 0, 0, 0, 65280, 0, 30464, 21760, 0, 34816, 43520, 0, 31232, 58880, 0, 0, 0}
	pc_red = {5000, 35535, 2, 0, 39321, 48059, 65535, 0, 16385, 65535, 0, 65535, 2, 0, 39321, 48059, 65535, 0, 16385, 65535}
	pc_green = {0, 16385, 39321, 0, 1, 48059, 32768, 65535, 65535, 32768, 0, 16385, 39321, 0, 1, 48059, 32768, 65535, 65535, 32768}
	pc_blue = {10000, 16385, 1, 65535, 31457, 48059, 32768, 0, 65535, 58981, 0, 16385, 1, 65535, 31457, 48059, 32768, 0, 65535, 58981}

	if(NewGraph==1)
		display /W=(35.25,42.5,503.25,736.25)
	endif
	variable i, j
	variable counter = 0
	string this_trace
	string this_axis
	variable bottomVal, topVal
	
	for(j=0;j<nCarbons;j+=1)
		bottomVal = (j/nCarbons)
		topVal = ((j+1)/nCarbons)-0.01
		for(i=0;i<nOxygens;i+=1)
			if(j==0)
				this_axis = "left"
				appendtograph wv[j][i][]
			else
				this_axis = "L_"+num2str(j)
				appendtograph/L=$this_axis wv[j][i][]
			endif
		
			if(counter==0)
				this_trace = wvstr
			else
				this_trace = wvstr+"#"+num2str(counter)
				ModifyGraph rgb($this_trace)=(pc_red[i],pc_green[i],pc_blue[i])
			endif
			counter +=1
		endfor
		ModifyGraph axisEnab($this_axis)={bottomVal,TopVal}
		ModifyGraph freePos($this_axis)={0,bottom}
		ModifyGraph tick($this_axis)=2,zero($this_axis)=1,standoff($this_axis)=0,mirror($this_axis)=1
	endfor
	ModifyGraph lsize=2
	
End

//*******************************************************************************
Function PopulateRunLog()

	setdatafolder root:
	NVAR Nparticles // number of particles (for heterogeneous chemistry)
	NVAR Pfrag_type // 0 = cfrag*Nox; 1 = (O:C)^mfrag
	NVAR ProbOx1, ProbOx2, ProbOx3, ProbOx4 // probability of adding X oxygens--> is renormalized later
	NVAR OHconc // [OH] in molecules/cm^3
	NVAR Ncarbons // number of carbon atoms in parent molecule
	NVAR FragSlope // either mfrag or cfrag, dependnig on fragmentation method
	NVAR delta_logCstar_perO // change in logCstar per oxygen added
	NVAR ctot_ppm // concentration of parent hydrocarbon (gas + particle) in ppm
	NVAR H_per_O // number of H atoms lost per O added (van Krevelen slope)
	NVAR krxn_parent // set specific value for parent compound
	NVAR SeedVolConc
	
	wave RunLog = root:savedmodelruns:RunLog
	
	string RunLogParams = "Time;Date;RunNumber; Nc; Pfunc1; Pfunc2; Pfunc3; Pfunc4; cfrag;DLVP;"
	RunLogParams += "HC_ppb;OHconc;Np;SeedConc; PfragMethod;kOHparent;HetChem;Oligomer;"
	
	variable nitems = itemsinlist(RunLogParams,";")
	variable nRows = dimsize(RunLog,0)
	redimension/N=(nRows+1,-1) RunLog
	
	string thing = "Nparticles"
	variable i
	for(i=0;i<nItems;i+=1)
		RunLog[nRows][i] = 1//$("Nparticles")
	endfor
	
End
//***************
Function CreateNewRunLog()

	string RunLogParams = "Time;Date;RunNumber; Nc; Pfunc1; Pfunc2; Pfunc3; Pfunc4; cfrag;DLVP;"
	RunLogParams += "HC_ppb;OHconc;Np;SeedConc; PfragMethod;kOHparent;HetChem;Oligomer;"
	
	variable nitems = itemsinlist(RunLogParams,";")
	make/o/d/n=(1,nItems) $("root:SavedModelRuns:RunLog")
	wave RunLog = $("root:SavedModelRuns:RunLog")
	
	edit runlog
	variable i
	string name
	for(i=0;i<nItems;i+=1)
		name = stringfromlist(i,RunLogParams,";")
		setdimlabel 1, i, $name, RunLog
	endfor
end

Function GraphSavedData(folderstr)
	string folderstr
	
	string cdf = getdatafolder(1)
	string dffp = "root:SavedModelRuns:"+folderstr
	if(datafolderexists(dffp)==0)
		abort "data folder does not exist"
	endif
	setdatafolder $dffp
	
	wave Coa_time
	wave O2C_time
	wave Dp_Time
	
	Display /W=(465,66.5,1005.75,509) Coa_time
	AppendToGraph/L=L_O2C O2C_time
	AppendToGraph/R Dp_time
	ModifyGraph lSize=2
	ModifyGraph lStyle(Dp_time)=3
	ModifyGraph rgb(Coa_time)=(0,0,0),rgb(O2C_time)=(65280,0,0),rgb(Dp_time)=(0,15872,65280)
	ModifyGraph tick=2
	ModifyGraph mirror(bottom)=1,mirror(L_O2C)=1
	ModifyGraph lblMargin(right)=9
	ModifyGraph standoff(left)=0,standoff(bottom)=0
	ModifyGraph axOffset(left)=-1.88889
	ModifyGraph axRGB(right)=(0,15872,65280)
	ModifyGraph tlblRGB(right)=(0,15872,65280)
	ModifyGraph alblRGB(right)=(0,15872,65280)
	ModifyGraph lblPos(left)=43,lblPos(L_O2C)=40
	ModifyGraph lblLatPos(left)=-2,lblLatPos(L_O2C)=-1,lblLatPos(right)=-1
	ModifyGraph freePos(L_O2C)={0,bottom}
	ModifyGraph axisEnab(left)={0,0.64}
	ModifyGraph axisEnab(L_O2C)={0.66,1}
	ModifyGraph axisEnab(right)={0,0.64}
	Label left "C\\BOA\\M (\\U)"
	Label bottom "Reaction Time (\\U)"
	Label right "Particle Diameter w/seed (\\U)"
	SetAxis L_O2C 0.2,*
	TextBox/C/N=text0/F=0/A=RB folderstr
	
	setdatafolder $cdf
End

//***************************************************************	
Function GetFractionByCarbonNumber()

	wave PM = root:ParticleMolecules_Time
	
	variable nC = dimsize(PM,0)
	variable nO = dimsize(PM,1)
	variable nTime = dimsize(PM,2)
	variable step = dimdelta(PM,2)
	variable offset = dimoffset(PM,2)
	NVAR H_per_O
	NVAR Hadjustment
	SOM_CreateAtomMatrices(nC,nO,H_per_O,Hadjustment)
	wave MW_matrix
	
	make/o/d/n=(nTime,nC) root:MassFractionByCarbon
	wave MassFractionByCarbon = root:MassFractionByCarbon
	setscale/P x, offset, step, "hours", MassFractionByCarbon
	make/o/d/n=(nO)/FREE ConcByCarbonAndTime
	make/o/d/n=(nC,nO)/FREE ConcByTime
	
	variable i, j
	variable TotalMass, MassByCarbon
	
	for(i=0;i<nTime;i+=1)
		ConcByTime = PM[p][q][i]*MW_matrix[p][q]
		wavestats/q ConcByTime
		TotalMass = V_Sum
		for(j=0;j<nC;j+=1)
			ConcByCarbonAndTime = PM[j][p][i]*MW_Matrix[j][p]
			wavestats/q ConcByCarbonAndTime
			MassByCarbon = V_sum
			MassFractionByCarbon[i][j] = MassByCarbon/TotalMass
		endfor
	endfor
	SOM_KillWaves()
End
//*****************************
Function GraphFractionByCarbonNumber()

	wave MF = root:massFractionByCarbon
	variable nC = dimsize(MF,1)
	
	display
	variable i
	for(i=0;i<nC;i+=1)
		appendtograph MF[][i]
	endfor

End
			
//*************************************************************************************************************************
// Calculate gas-particle equilbirum based on "Pankow" theory (i.e. Raoult's Law)
// Original from 02/07/09
// Modified on 09/03/13
Function EqmCalc(ConcMatrix,ConcSat25CForTD,MWwave,[Coa_guess,SeedConc,SeedMW,units,Tvar,deltaHvap_matrix])
	wave ConcMatrix // in molecules/cm^3 
	wave ConcSat25CForTD // saturation vapor pressure in molecules/m^3
	wave MWwave // molecular weight matrix in g/mol
	variable Coa_guess	// initial guess, in molecules/cm^3 or ug/m^3
	variable SeedConc // seed concentration in molecules/cm^3 or ug/m^3
	variable SeedMW // seed MW, if absorbing
	string units	// either "molecules" or "mass". "molecules" is default
	variable Tvar	// temperature, in K --> added 09/03/13
	wave deltaHvap_matrix	// matrix of deltaHvap values, in kJ/mol --> added 09/03/13

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
		wavestats/Q ConcMatrix
		Coa_Guess = V_sum/4
	endif
	variable Tadjust
	if(ParamIsDefault(Tvar) || ParamIsDefault(deltaHvap_matrix))	// added 09/03/13
		Tadjust = 0
	else
		Tadjust = 1
	endif
	
	Duplicate/O/D ConcMatrix, ConcWave, FracAero, Coa
	Duplicate/O/D ConcSat25CforTD SatConc_temp
	variable Coa_tol	// tolerance on result
	variable Seed
	
	if(stringmatch(units,"molecules"))
		ConcWave = ConcMatrix * 1e6		// convert molecules/cm^3 to molecules/m^3
		ConcWave = numtype(ConcWave)==2 ? 0 : ConcWave
		SatConc_temp = ConcSat25CforTD	// no units change
		Coa_tol = 1	// tolerance is 1 molecule
		Seed = SeedConc * 1e6 // convert molecules/cm^3 to molecules/m^3
	elseif(stringmatch(units,"mass"))
		ConcWave = ConcMatrix * 1e6 * 1e6 * MWwave / 6.022e23	// convert molecules/cm^3 to ug/m^3
		ConcWave = numtype(ConcWave)==2 ? 0 : ConcWave
		SatConc_temp = ConcSat25CforTD * 1e6 * MWwave / 6.022e23	// convert molecules/m^3 to ug/m^3
		Coa_tol = 1e-6	// tolerance is 1e-6 ug/m^3
		Seed = SeedConc * 1e6 * 1e6 * SeedMW / 6.022e23
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
	Coa_guess += Seed
	Variable Coa_calc
	Variable Coa_dif
	variable counter=0	// added 08/31/13
	
	do
		FracAero = 1/(1+(SatConc_temp)/Coa_guess)
		Coa = FracAero * ConcWave
		wavestats/Q Coa
		Coa_calc = V_Sum + Seed//sum(Coa)
		Coa_dif = abs(Coa_calc - Coa_guess)
		if(Coa_dif < Coa_tol)
			break
		endif
		Coa_guess = Coa_calc
		counter += 1
		if(counter==5)	// reset tolerance based on current guess; this speeds things up since the tolerance desired depends on the amount of OA there is.
			Coa_tol = max(Coa_tol,0.01*(Coa_calc-Seed))
		endif
	while(1)
	
	if(stringmatch(units,"molecules"))
		Coa = Coa * 1e-6 // convert molecules/m^3 to molecules/cm^3...does not include seed mass
	elseif(stringmatch(units,"mass"))
		Coa = Coa * 6.022e23 / (1e6 * 1e6 * MWwave)	// convert ug/m^3 to molecules/cm^3
	endif
End	
			
//********************************************************************************************
Function SetFromUCR(Expt,WLCyn)
	// Use this to grab data from the particular experiments in the folder UCR_Data and to move them to root as Coa_experiment, etc.
	string Expt	// string with name of UCR experiment...note that these names are not particularly telling
	variable WLCyn // wall loss correction? 0 = no, 1 = yes
	
	string dfdatafp = "root:UCR_data:"+Expt
	if(datafolderexists(dfdatafp)!=1)
		abort "no data set"
	endif
	if(WLCyn==0)
		wave PMvol = $(dfdatafp+":PMvol")	// PM volume concentration
	else
		wave PMvol = $(dfdatafp+":PMvol_corr")	// PM volume concentration
	endif
	wave PMtime	 = $(dfdatafp+":PMtime") // experiment time (in hours)
	wave OH_obs = $(dfdatafp+":OH_exp")
	wave OH_obs_time = $(dfdatafp+":OH_exp_time")
	wave RunConstants = $(dfdatafp+":RunConstants")
	wave PMnum = $(dfdatafp+":PMnum_corr")
	wavestats/q PMnum
	variable maxnum = V_max
//	
	setdatafolder root:
	make/o/d/n=(numpnts(PMvol)) Coa_experiment = PMvol*1.3	// multiply by assumed density
	if(WLCyn==0)
		note Coa_experiment "No wall loss correction"
	else
		note Coa_experiment "UCR wall loss correction"
	endif
	make/o/d/n=(numpnts(PMtime)) Time_Experiment = PMtime
	make/o/d/n=(numpnts(OH_obs)) OH_exp = OH_obs
	make/o/d/n=(numpnts(OH_obs_time)) OH_exp_time = OH_obs_time
	
	NVAR ctot_ppm
	ctot_ppm = RunConstants[0]
	NVAR krxn_parent 
	krxn_parent = RunConstants[1]
	NVAR Ncarbons
	Ncarbons = RunConstants[2]
	NVAR Hadjustment
	Hadjustment = RunConstants[3]
	NVAR UseScalingForOxidant
	UseScalingForOxidant = 2 // use the observed OH time series
	NVAR ParticleWallLoss
	if(WLCyn==0)
		ParticleWallLoss = 3.5
	else
		ParticleWallLoss = 0
	endif
	NVAR MaxTime_Hours
	wavestats/q Time_Experiment
	MaxTime_Hours = V_max + 0.1
	NVAR Nparticles
	Nparticles = MaxNum
	NVAR SeedVolConc
	SeedVolConc = 1e-4
	
	// only keep data after lights are turned on (i.e. after t = 0)
	variable crossat
	FindLevel/Q Time_Experiment, 0
	if(V_flag==0)
		CrossAt = ceil(V_LevelX)
		deletepoints 0, CrossAt, Coa_experiment, Time_Experiment
	endif
	// make error and masking waves
	make/o/d/n=(numpnts(Coa_experiment)) Coa_experiment_err = max(0.3,0.08*Coa_experiment)
	make/o/d/n=(numpnts(Coa_experiment)) Coa_experiment_mask = 1

	
	FindLevel/Q OH_exp_time, 0
	if(V_flag==0)
		CrossAt = ceil(V_LevelX)
		deletepoints 0, CrossAt, OH_exp, OH_exp_time
	endif	
End

//********************************************************************************************
Function/S GetDataFolderName_Caltech(Compound,NOx)
	string Compound // relates to name of compound
	string NOx // low or high

	string df_base = "root:Caltech_data"
	string df_summary = "root:DataSummary"
	string df_data
	variable i
	
	strswitch(compound)
		Case "dodecane":
			if(stringmatch(NOx,"low"))
				df_data = df_base + ":alkanes:C12:low_NOx:dodecane"
			else
				df_data = df_base + ":alkanes:C12:high_NOx:dodecane"
			endif
			print df_data
			break
		Case "HexylCycloHexane":
			if(stringmatch(NOx,"low"))
				df_data = df_base + ":alkanes:C12:low_NOx:HexylCycloHexane"
			else
				df_data = df_base + ":alkanes:C12:high_NOx:HexylCycloHexane"
			endif
			print df_data
			break
		Case "Cyclododecane":
			if(stringmatch(NOx,"low"))
				df_data = df_base + ":alkanes:C12:low_NOx:Cyclododecane"
			else
				df_data = df_base + ":alkanes:C12:high_NOx:Cyclododecane"
			endif
			print df_data
			break
		Case "Methylundecane":
			if(stringmatch(NOx,"low"))
				df_data = df_base + ":alkanes:C12:low_NOx:Methylundecane"
			else
				df_data = df_base + ":alkanes:C12:high_NOx:Methylundecane"
			endif
			print df_data
			break
		Case "aPinene":
			if(stringmatch(NOx,"low"))
				df_data = df_base + ":biogenics:OH:aPinene_lowNOx"
			elseif(stringmatch(NOx,"high"))
				df_data = df_base + ":biogenics:OH:aPinene_highNOx"
			elseif(stringmatch(NOx,"O3"))
				df_data = df_base + ":biogenics:O3:aPinene"
			endif
			print df_data
			break	
		Case "mXylene":
			if(stringmatch(NOx,"low"))
				df_data = df_base + ":aromatics:mXylene_101011_lowNOx"
			elseif(stringmatch(NOx,"low2"))
				df_data = df_base + ":aromatics:mXylene_061027_lowNOx"
			elseif(stringmatch(NOx,"high"))
				df_data = df_base + ":aromatics:mXylene_highNOx"
				abort "High NOx data does not exist for m-Xylene"
			elseif(stringmatch(NOx,"high2"))
				df_data = df_base + ":aromatics:mXylene_061005_highNOx"
			else
				print "you entered something wrong. try again."
				break 
			endif
			print df_data
			break	
		Case "Toluene":
			if(stringmatch(NOx,"low"))
				df_data = df_base + ":aromatics:Toluene_lowNOx"
				abort "low NOx data does not exist for toluene. Try low2."
			elseif(stringmatch(NOx,"low2"))
				df_data = df_base + ":aromatics:Toluene_061024_lowNOx"
			elseif(stringmatch(NOx,"high"))
				df_data = df_base + ":aromatics:Toluene_071208_highNOx"
			elseif(stringmatch(NOx,"high2"))
				df_data = df_base + ":aromatics:Toluene_061014_highNOx"
			else
				print "you entered something wrong. try again."
				break 
			endif
			print df_data
			break	
		Case "Naphthalene":
			if(stringmatch(NOx,"low"))
				df_data = df_base + ":aromatics:Naphthalene_lowNOx"
			elseif(stringmatch(NOx,"high"))
				df_data = df_base + ":aromatics:Naphthalene_highNOx"
			endif
			print df_data
			break					
		Case "Benzene":
			if(stringmatch(NOx,"low"))
				df_data = df_base + ":aromatics:Benzene_061104_lowNOx"
			elseif(stringmatch(NOx,"high"))
				df_data = df_base + ":aromatics:Benzene_070115_highNOx"
			else
				print "you entered something wrong. try again."
				break 
			endif
			print df_data
			break	
		Case "Isoprene":
			if(stringmatch(NOx,"low"))
				df_data = df_base + ":biogenics:OH:isoprene_lowNOx"			
			elseif(stringmatch(NOx,"high"))
				df_data = df_base + ":biogenics:OH:isoprene_highNOx"
			else
				print "please try again. compound or NOx doesn't exist"
			endif
			print df_data
			break
		Case "Methacrolein":
			if(stringmatch(NOx,"high"))
				df_data = df_base + ":biogenics:OH:Methacrolein_highNOx"
			else
				print "please try again. there is no lowNOx methacrolein data"
			endif
			print df_data
			break
		Case "Acrolein":
			if(stringmatch(NOx,"high"))
				df_data = df_base + ":biogenics:OH:acrolein_highNOx"
			else
				print "please try again. there is no lowNOx acrolein data"
			endif
			print df_data
			break
		Case "Loza2013":
			wave/t ExpLabel = $(df_base + ":alkanes:Loza2013:CompoundIdentifier")
			wave/t ExpDate = $(df_base + ":alkanes:Loza2013:ExpDate")
			variable indexmatch = -1
			for(i=0;i<numpnts(ExpLabel);i+=1)
				if(stringmatch(ExpLabel[i],NOx)==1)
					indexmatch = i
					break
				endif
			endfor
			if(indexmatch!=-1)
				df_data = df_base + ":alkanes:loza2013:expdata_"+ExpDate[indexmatch]
				print df_data
			else
				print "you fool! try again"
			endif
			break				
		Case "Toluene2013":
			if(stringmatch(NOx,"high*"))
				df_data = df_base + ":aromatics:Toluene2013:high_NOx"
				wave/t ExpDate = $(df_base + ":aromatics:Toluene2013:highNOxNames")
			elseif(stringmatch(NOx,"low*"))
				df_data = df_base + ":aromatics:Toluene2013:low_NOx"
				wave/t ExpDate = $(df_base + ":aromatics:Toluene2013:LowNOxNames")
			endif
			variable counter
			string countstr
			for(i=0;i<numpnts(ExpDate);i+=1)
				countstr = "*"+num2istr(i)
				if(stringmatch(NOx,countstr))
					df_data += ":"+ExpDate[i]
				endif
			endfor
			print df_data
			break	
		default:
			print "you have entered a compound name that does not exist"
			abort
		endswitch
		
		return df_data
End

//**************************************************************************************
Function SetFromCaltech(compound,NOx)
	string compound
	string NOx // "low" or "high"
	
	setdatafolder root:
	
	string df_data = GetDataFolderName_Caltech(compound,NOx)
	string df_summary = "root:DataSummary"
	string df_SOM = "root"
	
	string wavenames_Caltech = "Time_experiment;Coa_experiment;Coa_experiment_err;Coa_experiment_mask;O2C_experiment;O2C_experiment_err;O2C_experiment_mask;"
	wavenames_Caltech += "Time_VOC;VOC_ppm"
	string wavenames_SOM = wavenames_Caltech
	variable nWaves = itemsinlist(wavenames_Caltech)
	variable nRows
	variable i
	
	for(i=0;i<nWaves;i+=1)
		string CaltechStr = (df_data+":"+stringfromlist(i,wavenames_Caltech,";"))
		string SOMstr = (df_SOM+":"+stringfromlist(i,wavenames_SOM,";"))
		wave Caltech = $CaltechStr
		nRows = numpnts(Caltech)
		make/o/d/n=(nRows) $SOMstr
		wave SOM = $SOMstr
		SOM = Caltech
	endfor
	
	wave RunParameters = $(df_data+":RunParameters")
	variable nParams = numpnts(RunParameters)
	string RunInfoString
	if(stringmatch(NOx,"O3")==1)
		RunInfoString = "Dodecane;Ctot_ppm;krxn_parent;NCarbons;Hadjustment;UseScalingForOxidant;OHconc;OH_scale;H_per_O;Nparticles;SeedVolConc;O3_yn;O3_conc"
	else
		RunInfoString = "Dodecane;Ctot_ppm;krxn_parent;NCarbons;Hadjustment;UseScalingForOxidant;OHconc;OH_scale;H_per_O;Nparticles;SeedVolConc"
		NVAR O3_yn
		O3_yn = 0
	endif
	
	string current
	for(i=1;i<itemsinlist(RunInfoString);i+=1)
		current = stringfromlist(i,RunInfostring,";")
		NVAR var = $(current)
		var = RunParameters[i]
		if(stringmatch(current,"UseScalingForOxidant") && var==2)
			if(stringmatch(NOx,"O3"))
				wave Ox = $(df_data+":O3_experiment")
				wave Time_Ox = $(df_data+":Time_Experiment")
			else
				wave Ox = $(df_data+":OH_exp")
				wave Time_Ox = $(df_data+":Time_VOC")
			endif
			if(waveexists(Ox)==0)
				print "I don't think your oxidant wave exists. You should check this."
			endif
			duplicate/o/d Ox, $(df_SOM+":OH_exp")
			duplicate/o/d Time_Ox, $(df_SOM+":OH_exp_time")
		endif
	endfor

	NVAR MaxTime_Hours
	wave Time_Experiment
	wavestats/q Time_Experiment
	MaxTime_Hours = V_max + 0.1
	NVAR Nparticles
		
	setdatafolder root:
End

//**************************************************************************************
Function SetFromCaltech_Seed(compound,NOx)
	string compound
	string NOx // "low" or "high"
	
	setdatafolder root:
	
	string df_data = GetDataFolderName_Caltech(compound,NOx)
	string df_summary = "root:DataSummary"
	string df_SOM = "root"
	
	string wavenames_Caltech = "Time_experiment;Coa_experiment;Coa_experiment_err;Coa_experiment_mask;O2C_experiment;O2C_experiment_err;O2C_experiment_mask;"
	wavenames_Caltech += "Time_VOC;VOC_ppm"
	string wavenames_SOM = wavenames_Caltech
	variable nWaves = itemsinlist(wavenames_Caltech)
	variable nRows
	variable i
	
	for(i=0;i<nWaves;i+=1)
		string CaltechStr = (df_data+":"+stringfromlist(i,wavenames_Caltech,";"))
		string SOMstr = (df_SOM+":"+stringfromlist(i,wavenames_SOM,";"))
		wave Caltech = $CaltechStr
		nRows = numpnts(Caltech)
		make/o/d/n=(nRows) $SOMstr
		wave SOM = $SOMstr
		SOM = Caltech
	endfor
	
	wave RunParameters = $(df_data+":RunParameters")
	variable nParams = numpnts(RunParameters)
	string RunInfoString
	if(stringmatch(NOx,"O3")==1)
		RunInfoString = "Dodecane;Ctot_ppm;krxn_parent;NCarbons;Hadjustment;UseScalingForOxidant;OHconc;OH_scale;H_per_O;Nparticles;SeedVolConc;O3_yn;O3_conc"
	else
		RunInfoString = "Dodecane;Ctot_ppm;krxn_parent;NCarbons;Hadjustment;UseScalingForOxidant;OHconc;OH_scale;H_per_O;Nparticles;SeedVolConc;SeedDiameter;SizeSpread"
		NVAR O3_yn
		O3_yn = 0
	endif
	
	string current
	for(i=1;i<itemsinlist(RunInfoString);i+=1)
		current = stringfromlist(i,RunInfostring,";")
		NVAR var = $(current)
		var = RunParameters[i]
		if(stringmatch(current,"UseScalingForOxidant") && var==2)
			if(stringmatch(NOx,"O3"))
				wave Ox = $(df_data+":O3_experiment")
				wave Time_Ox = $(df_data+":Time_Experiment")
			else
				wave Ox = $(df_data+":OH_exp")
				wave Time_Ox = $(df_data+":Time_VOC")
			endif
			if(waveexists(Ox)==0)
				print "I don't think your oxidant wave exists. You should check this."
			endif
			duplicate/o/d Ox, $(df_SOM+":OH_exp")
			duplicate/o/d Time_Ox, $(df_SOM+":OH_exp_time")
		endif
	endfor

	NVAR MaxTime_Hours
	wave Time_Experiment
	wavestats/q Time_Experiment
	MaxTime_Hours = V_max + 0.1
	NVAR Nparticles
		
	setdatafolder root:
End

//**************************************************************************************
Function SetFromCaltech2(compound,UV,SA)
// For dealing with low/high UV data, specific to alpha-pinene experiments
	string compound
	string UV // "low" or "high"
	string SA // low, med or high
	
	setdatafolder root:
	
	string df_data = GetDataFolderName_Caltech2(compound,UV,SA)
	string df_summary = "root:DataSummary"
	string df_SOM = "root"
	
	string wavenames_Caltech = "Time_experiment;Coa_experiment;Coa_experiment_mask;Coa_experiment_err;"
	wavenames_Caltech += "Time_VOC;VOC_ppm;deltaHC_ugm3"
	string wavenames_SOM = wavenames_Caltech
	variable nWaves = itemsinlist(wavenames_Caltech)
	variable nRows
	variable i
	
	for(i=0;i<nWaves;i+=1)
		string CaltechStr = (df_data+":"+stringfromlist(i,wavenames_Caltech,";"))
		string SOMstr = (df_SOM+":"+stringfromlist(i,wavenames_SOM,";"))
		wave Caltech = $CaltechStr
		nRows = numpnts(Caltech)
		make/o/d/n=(nRows) $SOMstr
		wave SOM = $SOMstr
		SOM = Caltech
	endfor
	
	wave RunParameters = $(df_data+":RunParameters")
	variable nParams = numpnts(RunParameters)
	string RunInfoString

	RunInfoString = "Dodecane;Ctot_ppm;krxn_parent;NCarbons;Hadjustment;UseScalingForOxidant;OHconc;OH_scale;H_per_O;Nparticles;SeedVolConc;SeedDiameter"
	
	string current
	for(i=1;i<itemsinlist(RunInfoString);i+=1)
		current = stringfromlist(i,RunInfostring,";")
		NVAR var = $(current)
		var = RunParameters[i]
		if(stringmatch(current,"UseScalingForOxidant") && var==2)
				wave Ox = $(df_data+":OH_exp")
				wave Time_Ox = $(df_data+":Time_VOC")
			if(waveexists(Ox)==0)
				print "I don't think your oxidant wave exists. You should check this."
			endif
			duplicate/o/d Ox, $(df_SOM+":OH_exp")
			duplicate/o/d Time_Ox, $(df_SOM+":OH_exp_time")
		endif
	endfor

	NVAR MaxTime_Hours
	wave Time_Experiment
	wavestats/q Time_Experiment
	MaxTime_Hours = V_max + 0.1
	NVAR Nparticles
		
	setdatafolder root:
End

//********************************************************************************************
Function/S GetDataFolderName_Caltech2(Compound,UV,SA)
// For dealing with low/high UV data, specific to alpha-pinene experiments from 2016
	string Compound // relates to name of compound
	string UV // low or high
	string SA // low, med, or high

	string df_base = "root:Caltech_data"
	string df_summary = "root:DataSummary"
	string df_data
	variable i
	
	strswitch(compound)
		Case "aPinene":
			if(stringmatch(UV,"low")&&stringmatch(SA,"low"))
				df_data = df_base + ":biogenics:OH:aPinene_2016:lowUV_lowSA"
			elseif(stringmatch(UV,"low")&&stringmatch(SA,"high"))
				df_data = df_base + ":biogenics:OH:aPinene_2016:lowUV_highSA"
			elseif(stringmatch(UV,"high")&&stringmatch(SA,"low"))
				df_data = df_base + ":biogenics:OH:aPinene_2016:highUV_lowSA"
			elseif(stringmatch(UV,"high")&&stringmatch(SA,"med"))
				df_data = df_base + ":biogenics:OH:aPinene_2016:highUV_medSA"
			elseif(stringmatch(UV,"high")&&stringmatch(SA,"high"))
				df_data = df_base + ":biogenics:OH:aPinene_2016:highUV_highSA"
			endif
//			print df_data
			break		
		default:
			print "you have entered a compound name that does not exist"
			abort
		endswitch
		
		return df_data
End

//***************************************************************************************************
Function SOM_BatchFitCaltechData(startpos,stoppos)
	variable startpos, stoppos

	string df_FitInfo = "root:FitInfo"
	string df_Global = df_FitInfo+":NewGlobalFitSetup"
	string df_data = "root"
	
	// waves that contain information for individual fitting
	setdatafolder $df_FitInfo
		wave/T CompoundName = $(df_FitInfo+":compoundname")	// name of compound to load and then fit
		wave/T NOx = $(df_FitInfo+":NOx")			// NOx condition of fit. "low", "high" or ""
		wave/T RxnWave = $(df_FitInfo+":RxnWave")
		wave/T SaveWave = $(df_FitInfo+":SaveWave")
		wave/T FitTo = $(df_FitInfo+":FitTo")
		wave/T FragType = $(df_FitInfo+":FragType")
		wave Het_wave = $(df_FitInfo+":Het_Wave")
		wave Np_wave = $(df_FitInfo+":Np_wave")
		wave DpSeed_wave = $(df_FitInfo+":DpSeed_wave")
		wave GammaOH_wave = $(df_FitInfo+":GammaOH_wave")
		wave InitialGuesses = $(df_FitInfo+":InitialGuesses")
		wave FitResults = $(df_FitInfo+":FitResults")
		wave O2Cerr_scale = $(df_FitInfo+":O2Cerr_scale")
		wave ChiSq = $(df_FitInfo+":ChiSq")
		wave ChiSq_O2C = $(df_FitInfo+":ChiSq_O2C")
		wave RunTime = $(df_FitInfo+":RunTime")
		wave alphawave = $(df_FitInfo+":alpha_wave")
		wave DynamicPartitioning = $(df_FitInfo+":dynamicPartitioning_wave")
		wave WLRgas_wave = $(df_FitInfo+":WLRgas_wave")	// gas-phase wall loss rate (1/s)
		wave Cwall_wave = $(df_FitInfo+":Cwall_wave") // Cwall in mg/m^3
		wave kOH_wave = $(df_FitInfo+":kOH_wave") // method number
		
	// Fit Waves (coefficients and constraints) for "normal" fit
	setdatafolder root:
		wave cwave = root:coefwave
		wave uncertainty = root:w_sigma
		Make/O/T/N=12 root:T_Constraints
		wave/T T_constraints = root:T_constraints
		T_Constraints[0] = {"K0 > 0.01","K0 < 10","K1 > .7","K1 < 2.5","K2 > 10e-4","K2 < 3","K3 > 10e-4","K3 < 3","K4 > 10e-4","K4 < 3","K5 > 10e-4","K5 < 3"}
//		T_Constraints[10] = {"K5 > 1e-4","K5 < 1"}
		
	// For Global Fitting to Coa and O2C
//		Make/O/T/N=12 T_constraints_Global
	setdatafolder $df_Global// root:FitInfo:NewGlobalFitSetup
		//make/o/t/n=2 $(df_Global+":NewGF_FitFuncNames") = {"allatonce_FitCoa","allatonce_FitO2C"}
		wave/T FitFuncNames = $(df_Global+":NewGF_FitFuncNames")
		//make/o/t/n=(2,4) $(df_Global+":NewGF_DataSetsList")
		wave/T DataSets = $(df_Global+":NewGF_DataSetsList") // assumes data are in Coa_experiment, O2C_experiment and that there are error waves
		//DataSets[][] = {{"root:Coa_experiment","root:O2C_experiment"},{"root:Time_experiment","root:Time_experiment"},{"root:Coa_experiment_err","root:O2C_experiment_err"},{"root:Coa_experiment_mask","root:O2C_experiment_mask"}}
		wave CoefDataSetLinkage = $(df_Global+":CoefDataSetLinkage")
		wave GlobalCoefWave = $(df_Global+":NewGF_CoefWave")
		GlobalCoefWave[][0] = cwave[p]
		GlobalCoefWave[][1] = {0, 0, 0, 0, 0, 0}
		GlobalCoefWave[][2] = {1e-3, 1e-3, 1e-4, 1e-4, 1e-4, 1e-4}
		wave/T CoefNames = NewGF_coefficientnames
		wave/T ConstraintWave = T_Constraints//SimpleConstraintsListWave // same constraints as T_Constraints above, just listed differently
		variable options_var = 32+64 // NewGFOptionWTISSTD=use std, not 1/std; NewGFOptionMAKE_FIT_WAVES
		variable FitCurvePoints_var // set below
		variable DoAlertsOnError_var = 0
		variable maxIters_var = 12
		wave W_sigma_global = W_sigma
	
	setdatafolder root:
		
	// Global Variables to set
		NVAR Pfrag_type
		NVAR hetchem
		NVAR Nparticles
		NVAR GammaOH
		SVAR SaveFolderStr
		NVAR gasWLR	// gas-wall loss rate (1/s)
		NVAR Cwall		// effective wall concentration (mg/m^3)
		NVAR timestep
		NVAR alpha
		NVAR kineticmasstransfer
		NVAR SeedDiameter
		NVAR krxn_method
	
	variable i, j
	string current
	string cdf
	
	for(i=startpos;i<=stoppos;i+=1)
	
	// Set things up
		// set reaction conditions
			SetFromCaltech(CompoundName[i],NOx[i])
			setdatafolder root:
			if(stringmatch(NOx[i],"low"))
				timestep = 150
			else
				timestep = 75
			endif
		// scale O2C error, if desired
			wave O2Cerr = root:O2C_experiment_err
			O2Cerr *= O2Cerr_scale[i]
		// set initial guesses for frag, dLVP, Pfunc
			for(j=0;j<6;j+=1)
				cwave[j] = InitialGuesses[i][j]
			endfor
			GlobalCoefWave[][0] = cwave[p]
		// select fragmentation parameterization
			if(stringmatch(FragType[i],"cfrag")==1)
				Pfrag_type = 0 // cfrag
			elseif(stringmatch(FragType[i],"mfrag")==1)
				Pfrag_type = 1 // mfrag
			elseif(stringmatch(FragType[i],"cfrag+")==1)
				Pfrag_type = 2 // cfrag + 1
			elseif(stringmatch(Fragtype[i],"constant")==1)
				Pfrag_type = 3 // constant
			endif
		// select gas-phase kinetics kOH method
			krxn_method = kOH_wave[i]
		// select whether to include heterogenous chemistry
			if(het_wave[i]==1)
				hetchem = 1
				GammaOH = 1//gammaoh_wave[i]
			else
				hetchem = 0
			endif	
		// Set particle diameter and number concentration
			SeedDiameter = DpSeed_wave[i]
			Nparticles = Np_wave[i]
		// Set gas-phase wall loss rate and Cwall
			gasWLR = WLRgas_wave[i]
			Cwall = Cwall_wave[i]
		// Set mass transfer method (eqm. or dynamic) and alpha (if dynamic)
			kineticmasstransfer = DynamicPartitioning[i]
			if(kineticmasstransfer==0)
				alpha = 1
			else
				alpha = alphawave[i]
			endif
		// Set folder to save results into
			SaveFolderStr = "Run"+num2istr(i)//SaveWave[m-1]
		
		// Do Fit
			DoWindow/F sompanel
			DoWindow/F FitGraph
			current = ""
			current = CompoundName[i] + "_" + FitTo[i] + "_" + FragType[i] + "_" + num2str(Het_Wave[i]) + "_#" + num2str(i)
			textbox/C/N=text0/A=MC current
			print current
			DoUpdate
		
			// select what type of fit to do
			if(stringmatch(FitTo[i],"Coa")==1)
				//FuncFit/NTHR=0/H="000000" allatonce_FitCoa CoefWave  Coa_experiment /X=Time_Experiment /W=Coa_experiment_err /I=1 /E=epsilon /D /C=T_Constraints 
				//FuncFit/NTHR=0/H="000000" allatonce_FitCoa CoefWave  Coa_experiment /X=Time_Experiment /W=Coa_experiment_err /I=1 /D /C=T_Constraints 
				wave Coa_experiment = $(df_data+":Coa_experiment")
				wave Time_Experiment = $(df_data+":Time_experiment")
				wave Coa_experiment_err = $(df_data+":Coa_experiment_err") 
				FuncFit/NTHR=0/H="000000" allatonce_FitCoa cwave  Coa_experiment /X=Time_Experiment /W=Coa_experiment_err /I=1 /D /C=T_Constraints 
			elseif(stringmatch(FitTo[i],"Coa_OC")==1) // for fitting simultaneously to Coa and O:C
				wave Coa_wv = root:Coa_experiment
				FitCurvePoints_var = numpnts(Coa_experiment)
				Duplicate/O W_sigma_global root:Packages:NewGlobalFit:W_sigma
				DoNewGlobalFit(FitFuncNames, DataSets, CoefDataSetLinkage, GlobalCoefWave, CoefNames, ConstraintWave, Options_var, FitCurvePoints_var, DoAlertsOnError_var, maxIters=maxIters_var)
			else
				print "you are obviously stupid. You can only fit to Coa or Coa_OC."
			endif
			
		// refresh these waves, just in case
			wave cwave = root:coefwave
			wave uncertainty = root:w_sigma

			if(stringmatch(FitTo[i],"Coa_OC")==1)
				wave cwave_Global = GlobalFitCoefficients
				cwave = cwave_Global
			endif
			
		// Run one more time with best-fit results
			NVAR FragSlope // either mfrag or cfrag, dependnig on fragmentation method
			NVAR delta_logCstar_perO // change in logCstar per oxygen added
			NVAR ProbOx1, ProbOx2, ProbOx3, ProbOx4
		
		make/o/d/n=4 cshort
		cshort = cwave[x+2]
		wavestats/q cshort
		cshort /= V_sum
		cwave[2,] = cshort[x-2]
		FitResults[i][][0] = cwave[q]
		FitResults[i][][1] = uncertainty[q]
		NVAR CS = root:ChiSq_Coa	
		ChiSq[i] = CS
		NVAR CS2 = root:ChiSq_O2C
		ChiSq_O2C[i] = CS2
		runtime[i] = datetime
		DoUpdate
	endfor
End

//***************************************************************************************************
Function SOM_BatchFit(startpos,stoppos,dataset)
// For batch fitting of data (e.g. different data sets, or the same data sets with different assumptions)
// Currently setup to work only with data from Caltech or UCR, stored in particular folders; can easily be extended
// The inputs for this are setup 
// 08/27/13: updated to include Cwall_wave and WLRg_scaling_wave
	variable startpos, stoppos
	string dataset // name of dataset to use, "caltech", "UCR", ""

	string df_FitInfo = "root:FitInfo"
	string df_Global = df_FitInfo+":NewGlobalFitSetup"
	string df_data = "root"
	
	// waves that contain information for individual fitting, all stored in df_FitInfo
	setdatafolder $df_FitInfo
		wave/T CompoundName = $(df_FitInfo+":compoundname")	// name of compound to load and then fit
		wave/T NOx = $(df_FitInfo+":NOx")			// NOx condition of fit. "low", "high" or ""
		wave WLRgas_wave = $(df_FitInfo+":WLRgas_wave")	// gas-phase wall loss rate (1/s)
	
		wave/T RxnWave = $(df_FitInfo+":RxnWave")
		wave/T SaveWave = $(df_FitInfo+":SaveWave")
		wave/T FitTo = $(df_FitInfo+":FitTo")
		wave/T FragType = $(df_FitInfo+":FragType")
		wave Het_wave = $(df_FitInfo+":Het_Wave")
		wave Np_wave = $(df_FitInfo+":Np_wave")
		wave GammaOH_wave = $(df_FitInfo+":GammaOH_wave")
		wave InitialGuesses = $(df_FitInfo+":InitialGuesses")
		wave FitResults = $(df_FitInfo+":FitResults")
		wave O2Cerr_scale = $(df_FitInfo+":O2Cerr_scale")
		wave ChiSq = $(df_FitInfo+":ChiSq")
		wave ChiSq_O2C = $(df_FitInfo+":ChiSq_O2C")
		wave RunTime = $(df_FitInfo+":RunTime")
		wave Cwall_wave = $(df_FitInfo+":Cwall_wave")
		wave WLRg_scaling_wave = $(df_FitInfo+":WLRg_scaling_wave")
		wave SPMwave = $(df_FitInfo+":SPMwave")
		wave kOHwave = $(df_FitInfo+":kOH_wave")
		wave ErrorWave = $(df_FitInfo+":ErrorWave")
		wave MaskWave = $(df_FitInfo+":MaskWave")
		wave/T holdStrWave = $(df_FitInfo+":holdStrWave")
		wave DpSeed_wave = $(df_FitInfo+":DpSeed_wave")						// added 11/12/13
		wave DynamicPartitioning_wave = $(df_FitInfo+":DynamicPartitioning_wave")	// added 11/12/13
		wave alpha_wave = $(df_FitInfo+":alpha_wave")							// added 11/12/13
		wave WLRp_wave = $(df_FitInfo+":WLRp_wave")							// added 08/27/14
	
	// Fit Waves (coefficients and constraints) for "normal" fit
	setdatafolder root:
		wave cwave = $(df_data+":coefwave")
		wave uncertainty = $(df_data+":w_sigma")
		Make/O/T/N=12 $(df_Data+":T_Constraints")
		wave/T T_constraints = $(df_Data+":T_Constraints")
		T_Constraints[0] = {"K0 > 0.01","K0 < 1","K1 > .7","K1 < 2.5","K2 > 1e-4","K2 < 1","K3 > 1e-4","K3 < 1","K4 > 1e-4","K4 < 1","K5 > 1e-4","K5 < 1"}
		
	// For Global Fitting to Coa and O2C (uncommon now)
	setdatafolder $df_Global
		wave/T FitFuncNames = $(df_Global+":NewGF_FitFuncNames")
		wave/T DataSets = $(df_Global+":NewGF_DataSetsList") // assumes data are in Coa_experiment, O2C_experiment and that there are error waves
		wave CoefDataSetLinkage = $(df_Global+":CoefDataSetLinkage")
		wave GlobalCoefWave = $(df_Global+":NewGF_CoefWave")
		GlobalCoefWave[][0] = cwave[p]
		GlobalCoefWave[][1] = {0, 0, 0, 0, 0, 0}
		GlobalCoefWave[][2] = {1e-3, 1e-3, 1e-4, 1e-4, 1e-4, 1e-4}
		wave/T CoefNames = $(df_Global+":NewGF_coefficientnames")
		wave/T ConstraintWave = $(df_Global+":T_Constraints//SimpleConstraintsListWave") // same constraints as T_Constraints above, just listed differently
		variable options_var = 32+64 // NewGFOptionWTISSTD=use std, not 1/std; NewGFOptionMAKE_FIT_WAVES
		variable FitCurvePoints_var // set below
		variable DoAlertsOnError_var = 0
		variable maxIters_var = 7
		wave W_sigma_global = W_sigma
	
	setdatafolder root:
		
	// Global Variables to set; these all live in root
		NVAR Pfrag_type
		NVAR hetchem
		NVAR Nparticles
		NVAR GammaOH
		SVAR SaveFolderStr
		NVAR gasWLR	// gas-wall loss rate (1/s)
		NVAR Cwall		// effective wall OA concentration (mg/m^3); added 8/27/13
		NVAR kwall_gas_scaling	// composition dependent wall loss rate??; added 08/27/13
		NVAR timestep
		NVAR SPM
		NVAR krxn_method
		NVAR alpha					// added 11/12/13
		NVAR SeedDiameter			// added 11/12/13
		NVAR KineticMassTransfer	// added 11/12/13
		NVAR ParticleWallLoss		// added 08/27/14 // This is a scaling factor applied to the size dependent wall-loss
		
	variable i, j
	string current
	variable V_fiterror = 0	// suppress abort on error
	variable V_fitMaxIters = 20
	variable reset = 0
	
	for(i=startpos;i<=stoppos;i+=1)
		// set reaction conditions
			if(stringmatch(dataset,"caltech"))
				SetFromCaltech(CompoundName[i],NOx[i])
			elseif(stringmatch(dataset,"UCR"))
				SetFromUCR(CompoundName[i],1)	// note that "compoundname" here is actually an experiment name; the "1" indicates use of wall-loss corrected data
			else
				abort "no dataset selected"
			endif
			if(reset==1)
				initialguesses[i][] = fitresults[i-1][q]
			endif
			setdatafolder root:
			if(stringmatch(NOx[i],"low"))
				timestep = 150
			elseif(stringmatch(NOx[i],"high"))
				timestep = 75
			elseif(stringmatch(dataset,"UCR"))
				timestep = 100
			else
				timestep = 100
			endif
		// scale O2C error, if desired
			wave O2Cerr = root:O2C_experiment_err
			O2Cerr *= O2Cerr_scale[i]
		// set initial guesses for frag, dLVP, Pfunc
			for(j=0;j<6;j+=1)
				cwave[j] = InitialGuesses[i][j]
			endfor
			GlobalCoefWave[][0] = cwave[p]
		// select fragmentation parameterization
			if(stringmatch(FragType[i],"cfrag")==1)
				Pfrag_type = 0 // cfrag
			elseif(stringmatch(FragType[i],"mfrag")==1)
				Pfrag_type = 1 // mfrag
			elseif(stringmatch(FragType[i],"cfrag+")==1)
				Pfrag_type = 2 // cfrag + 1
			elseif(stringmatch(Fragtype[i],"constant")==1)
				Pfrag_type = 3 // constant
			elseif(stringmatch(Fragtype[i],"qfrag")==1)
				Pfrag_type = 4 // constant
			endif
		// select whether to include heterogenous chemistry
			if(het_wave[i]==1)
				hetchem = 1
				//Nparticles = Np_wave[i]
				GammaOH = 1//gammaoh_wave[i]
			else
				hetchem = 0
				//Nparticles = Np_wave[i]
			endif	
		// Set gas-phase wall loss rate
			gasWLR = WLRgas_wave[i]
			Cwall = Cwall_wave[i]
			kwall_gas_scaling = WLRg_scaling_wave[i]
		// sequential partitioning
			SPM = SPMwave[i]
		// kOH method
			krxn_method = kOHwave[i]
		// Set folder to save results into
			SaveFolderStr = "Run"+num2istr(i)
		// Set some other things, added 11/12/13
			SeedDiameter = DpSeed_wave[i]
			Nparticles = Np_wave[i]
			KineticMassTransfer = DynamicPartitioning_Wave[i]
			alpha = alpha_wave[i]
		// Particle Wall Loss (monodisperse) // added 08/27/14
			ParticleWallLoss = WLRp_wave[i]
		// Do Fit
			DoWindow/F sompanel
			DoWindow/F FitWindow
			current = ""
			current = CompoundName[i] + "_" + FitTo[i] + "_" + FragType[i] + "_" + num2str(Het_Wave[i]) + "_#" + num2str(i)
			textbox/C/N=text0/A=MC current
			print current
			DoUpdate
		
			if(stringmatch(FragType[i],"cfrag"))
				Make/O/T/N=12 root:T_Constraints
				wave/T T_constraints = root:T_constraints
				T_Constraints[0] = {"K0 > 0.001","K0 < 0.999","K1 > .7","K1 < 2.5","K2 > 1e-3","K2 < 3","K3 > 1e-3","K3 < 3","K4 > 1e-3","K4 < 3","K5 > 1e-3","K5 < 3"}
			elseif(stringmatch(FragType[i],"mfrag"))
				Make/O/T/N=12 root:T_Constraints
				wave/T T_constraints = root:T_constraints
				T_Constraints[0] = {"K0 > 0.01","K0 < 5","K1 > .7","K1 < 2.5","K2 > 1e-3","K2 < 3","K3 > 1e-3","K3 < 3","K4 > 2e-3","K4 < 3","K5 > 1e-3","K5 < 3"}
			elseif(stringmatch(FragType[i],"qfrag"))
				Make/O/T/N=12 root:T_Constraints
				wave/T T_constraints = root:T_constraints
				T_Constraints[0] = {"K0 > 0.001","K0 < 100","K1 > .7","K1 < 2.5","K2 > 1e-3","K2 < 3","K3 > 1e-3","K3 < 3","K4 > 1e-3","K4 < 3","K5 > 1e-3","K5 < 3"}
			endif
			if(stringmatch(NOx[i],"O3"))
				Make/O/T/N=10 root:T_Constraints
				wave/T T_constraints = root:T_constraints
				T_Constraints[0] = {"K1 > .7","K1 < 2.5","K2 > 1e-3","K2 < 3","K3 > 1e-3","K3 < 3","K4 > 1e-3","K4 < 3","K5 > 1e-3","K5 < 3"}
			endif
//			if(holdDLVPwave[i]==1)
			if(stringmatch(holdStrWave[i],"DLVP"))
				if(stringmatch(NOx[i],"O3"))
					deletepoints 0, 2, T_constraints
				else
					deletepoints 2, 2, T_constraints
				endif
			elseif(stringmatch(holdStrWave[i],"Pfrag"))
				if(stringmatch(NOx[i],"O3")!=1)
					deletepoints 0,2, T_constraints
				endif
			endif
			// select what type of fit to do
			if(stringmatch(FitTo[i],"Coa")==1)
				string HoldStr
				if(stringmatch(NOx[i],"O3")==1 && stringmatch(holdStrWave[i],"none"))
					HoldStr = "100000"
				elseif(stringmatch(NOx[i],"O3")==1 && stringmatch(holdStrWave[i],"DLVP"))
					HoldStr = "110000"
				elseif(stringmatch(NOx[i],"O3")!=1 && stringmatch(holdStrWave[i],"DLVP"))
					HoldStr = "010000"
				elseif(stringmatch(NOx[i],"O3")!=1 && stringmatch(holdStrWave[i],"Pfrag"))
					HoldStr = "100000"
				else
					HoldStr = "000000"
				endif
//				if(stringmatch(NOx[i],"O3")==1 && holdDLVPwave[i]==0)
//					HoldStr = "100000"
//				elseif(stringmatch(NOx[i],"O3")==1 && holdDLVPwave[i]==1)
//					HoldStr = "110000"
//				elseif(stringmatch(NOx[i],"O3")!=1 && holdDLVPwave[i]==1)
//					HoldStr = "010000"
//				else
//					HoldStr = "000000"
//				endif//				if(stringmatch(NOx[i],"O3")==1)
//					if(ErrorWave[i]==1 && MaskWave[i]==1)
//						FuncFit/NTHR=0/H="100000" allatonce_FitCoa CoefWave  Coa_experiment /X=Time_Experiment /W=Coa_experiment_err /M=Coa_experiment_mask /I=1 /D /C=T_Constraints 
//					elseif(ErrorWave[i]==1 && MaskWave[i] == 0)
//						FuncFit/NTHR=0/H="100000" allatonce_FitCoa CoefWave  Coa_experiment /X=Time_Experiment /W=Coa_experiment_err /I=1 /D /C=T_Constraints 
//					elseif(ErrorWave[i]==0 && MaskWave[i]==1)
//						FuncFit/NTHR=0/H="100000" allatonce_FitCoa CoefWave  Coa_experiment /X=Time_Experiment /M=Coa_experiment_mask /I=1 /D /C=T_Constraints 
//					else
//						FuncFit/NTHR=0/H="100000" allatonce_FitCoa CoefWave  Coa_experiment /X=Time_Experiment /I=1 /D /C=T_Constraints 
//					endif
//				else
//					if(ErrorWave[i]==1 && MaskWave[i]==1)
//						FuncFit/NTHR=0/H="000000" allatonce_FitCoa CoefWave  Coa_experiment /X=Time_Experiment /W=Coa_experiment_err /M=Coa_experiment_mask /I=1 /D /C=T_Constraints 
//					elseif(ErrorWave[i]==1 && MaskWave[i] == 0)
//						FuncFit/NTHR=0/H="000000" allatonce_FitCoa CoefWave  Coa_experiment /X=Time_Experiment /W=Coa_experiment_err /I=1 /D /C=T_Constraints 					
//					elseif(ErrorWave[i]==0 && MaskWave[i]==1)
//						FuncFit/NTHR=0/H="000000" allatonce_FitCoa CoefWave  Coa_experiment /X=Time_Experiment /M=Coa_experiment_mask /I=1 /D /C=T_Constraints 					
//					else
//						FuncFit/NTHR=0/H="000000" allatonce_FitCoa CoefWave  Coa_experiment /X=Time_Experiment /I=1 /D /C=T_Constraints 
//					endif
//				endif
				maxIters_var = 20 // may want to change this
				wave Coa_experiment = $(df_Data+":Coa_experiment")
				wave Time_Experiment = $(df_Data+":Time_experiment")
				wave Coa_experiment_err = $(df_Data+":Coa_experiment_err")
				wave Coa_experiment_mask = $(df_Data+":Coa_experiment_mask")
				if(ErrorWave[i]==1 && MaskWave[i]==1)
					FuncFit/NTHR=0/H=holdStr allatonce_FitCoa cwave  Coa_experiment /X=Time_Experiment /W=Coa_experiment_err /M=Coa_experiment_mask /I=1 /D /C=T_Constraints 
				elseif(ErrorWave[i]==1 && MaskWave[i] == 0)
					FuncFit/NTHR=0/H=holdStr allatonce_FitCoa cwave  Coa_experiment /X=Time_Experiment /W=Coa_experiment_err /I=1 /D /C=T_Constraints 					
				elseif(ErrorWave[i]==0 && MaskWave[i]==1)
					FuncFit/NTHR=0/H=holdStr allatonce_FitCoa cwave  Coa_experiment /X=Time_Experiment /M=Coa_experiment_mask /I=1 /D /C=T_Constraints 					
				else
					FuncFit/NTHR=0/H=holdStr allatonce_FitCoa cwave  Coa_experiment /X=Time_Experiment /I=1 /D /C=T_Constraints 
				endif
			elseif(stringmatch(FitTo[i],"Coa_OC")==1) // for fitting simultaneously to Coa and O:C
				wave Coa_wv = root:Coa_experiment
				FitCurvePoints_var = numpnts(Coa_experiment)
				Duplicate/O W_sigma_global root:Packages:NewGlobalFit:W_sigma
				DoNewGlobalFit(FitFuncNames, DataSets, CoefDataSetLinkage, GlobalCoefWave, CoefNames, ConstraintWave, Options_var, FitCurvePoints_var, DoAlertsOnError_var, maxIters=maxIters_var)
			else
				print "oops. You can only fit to Coa or Coa_OC."
			endif
			
		// refresh these waves, just in case
			wave cwave = root:coefwave
			wave uncertainty = root:w_sigma

			if(stringmatch(FitTo[i],"Coa_OC")==1)
				wave cwave_Global = GlobalFitCoefficients
				cwave = cwave_Global
			endif
			
		if(V_fiterror!=0)
			cwave = 0
			V_fiterror = 0
		endif
			
	// Populate FitResults, etc. waves with results	
		make/o/d/n=4 cshort
		cshort = cwave[x+2]
		wavestats/q cshort
		cshort /= V_sum
		cwave[2,] = cshort[x-2]
		FitResults[i][][0] = cwave[q]
		FitResults[i][][1] = uncertainty[q]
		NVAR CS = root:ChiSq_Coa	
		ChiSq[i] = CS
		NVAR CS2 = root:ChiSq_O2C
		ChiSq_O2C[i] = CS2
		runtime[i] = datetime
		RecalculateChi2(i)
//		GetStaticParams("2P",0,i,i)
//		GetStaticParams("VBS",0,i,i)
		DoUpdate
	endfor
End

//**************************************************************************************************
Function SOM_BatchFit_FitTwice(startpos,stoppos,dataset)
	// 4/18/2014 - written by CDC
	// there seems to be some issues with the "final" fit from a single fit not necessarily being the "best" fit.
	// this function runs a fit once based on some initial guesses, but then takes the results from the fit and 
	// performs a second fit using the results from the first fit as initial guesses
	variable startpos,stoppos
	string dataset

	string df_FitInfo = "root:FitInfo"
	string df_data = "root"

	variable i
	for(i=startpos;i<=stoppos;i+=1)
		// First run
		SOM_BatchFit_CoaOnly(i,i,dataset)
		// Second run, only if first run was successful
		wave FitResults = $(df_FitInfo+":FitResults")	// note, 3D matrix
		wave InitialGuesses = $(df_FitInfo+":InitialGuesses") // 2D matrix
		make/o/d/n=(6)/FREE FitResultsTemp
		FitResultsTemp = FitResults[i][p][0]
		wavestats/q FitResultsTemp
		if(V_numnans==0 && V_min > 0)
			InitialGuesses[i][] = FitResults[i][q][0]
			SOM_BatchFit_CoaOnly(i,i,dataset)
		endif
	endfor	
End


//***************************************************************************************************
Function SOM_BatchFit_CoaOnly(startpos,stoppos,dataset)
	// 04/18/2014 - written by CDC
	// Simplified version of SOM_BatchFit that assumes that one only wants to fit Coa, not Coa and O2C
	// thus, no need to access global fit parameters
	variable startpos, stoppos
	string dataset // name of dataset to use, "caltech", "UCR", ""

	string df_FitInfo = "root:FitInfo"
	string df_Global = df_FitInfo+":NewGlobalFitSetup"
	string df_data = "root"
	
	// waves that contain information for individual fitting
	setdatafolder $df_FitInfo
		wave/T CompoundName = $(df_FitInfo+":compoundname")	// name of compound to load and then fit
		wave/T NOx = $(df_FitInfo+":NOx")			// NOx condition of fit. "low", "high" or ""
		wave WLRgas_wave = $(df_FitInfo+":WLRgas_wave")	// gas-phase wall loss rate (1/s)
	
		wave/T RxnWave = $(df_FitInfo+":RxnWave")
		wave/T SaveWave = $(df_FitInfo+":SaveWave")
		wave/T FitTo = $(df_FitInfo+":FitTo")
		wave/T FragType = $(df_FitInfo+":FragType")
		wave Het_wave = $(df_FitInfo+":Het_Wave")
		wave Np_wave = $(df_FitInfo+":Np_wave")
		wave GammaOH_wave = $(df_FitInfo+":GammaOH_wave")
		wave InitialGuesses = $(df_FitInfo+":InitialGuesses")
		wave FitResults = $(df_FitInfo+":FitResults")
		wave O2Cerr_scale = $(df_FitInfo+":O2Cerr_scale")
		wave ChiSq = $(df_FitInfo+":ChiSq")
		wave ChiSq_O2C = $(df_FitInfo+":ChiSq_O2C")
		wave RunTime = $(df_FitInfo+":RunTime")
		wave Cwall_wave = $(df_FitInfo+":Cwall_wave")
		wave WLRg_scaling_wave = $(df_FitInfo+":WLRg_scaling_wave")
		wave SPMwave = $(df_FitInfo+":SPMwave")
		wave kOHwave = $(df_FitInfo+":kOH_wave")
		wave ErrorWave = $(df_FitInfo+":ErrorWave")
		wave MaskWave = $(df_FitInfo+":MaskWave")
		wave/T holdStrWave = $(df_FitInfo+":holdStrWave")
		wave DpSeed_wave = $(df_FitInfo+":DpSeed_wave")						// added 11/12/13
		wave DynamicPartitioning_wave = $(df_FitInfo+":DynamicPartitioning_wave")	// added 11/12/13
		wave alpha_wave = $(df_FitInfo+":alpha_wave")							// added 11/12/13
	
	// Fit Waves (coefficients and constraints) for "normal" fit
	setdatafolder root:
		make/o/d/n=(6) $(df_data+":coefwave")
		wave cwave = $(df_data+":coefwave")
		wave uncertainty = $(df_data+":w_sigma")
		Make/O/T/N=12 $(df_Data+":T_Constraints")
		wave/T T_constraints = $(df_Data+":T_Constraints")
		T_Constraints[0] = {"K0 > 0.01","K0 < 1","K1 > .7","K1 < 2.5","K2 > 1e-4","K2 < 1","K3 > 1e-4","K3 < 1","K4 > 1e-4","K4 < 1","K5 > 1e-4","K5 < 1"}
		
	// Global Variables to set
		NVAR Pfrag_type
		NVAR hetchem
		NVAR Nparticles
		NVAR GammaOH
		SVAR SaveFolderStr
		NVAR gasWLR	// gas-wall loss rate (1/s)
		NVAR Cwall		// effective wall OA concentration (mg/m^3); added 8/27/13
		NVAR kwall_gas_scaling	// composition dependent wall loss rate??; added 08/27/13
		NVAR timestep
		NVAR SPM
		NVAR krxn_method
		NVAR alpha					// added 11/12/13
		NVAR SeedDiameter			// added 11/12/13
		NVAR KineticMassTransfer	// added 11/12/13
		
	
	variable i, j
	string current
	variable V_fiterror = 0	// suppress abort on error
	variable V_fitMaxIters = 20
	variable reset = 0
	
	for(i=startpos;i<=stoppos;i+=1)
		// set reaction conditions
			if(stringmatch(dataset,"caltech"))
				SetFromCaltech(CompoundName[i],NOx[i])
			elseif(stringmatch(dataset,"UCR"))
				SetFromUCR(CompoundName[i],1)	// note that "compoundname" here is actually an experiment name; the "1" indicates use of wall-loss corrected data
			else
				abort "no dataset selected"
			endif
			if(reset==1)
				initialguesses[i][] = fitresults[i-1][q]
			endif
			setdatafolder root:
			if(stringmatch(NOx[i],"low"))
				timestep = 100
			elseif(stringmatch(NOx[i],"high"))
				timestep = 75
			elseif(stringmatch(dataset,"UCR"))
				timestep = 100
			else
				timestep = 100
			endif
		// set initial guesses for frag, dLVP, Pfunc
			for(j=0;j<6;j+=1)
				cwave[j] = InitialGuesses[i][j]
			endfor
		// select fragmentation parameterization
			if(stringmatch(FragType[i],"cfrag")==1)
				Pfrag_type = 0 // cfrag
			elseif(stringmatch(FragType[i],"mfrag")==1)
				Pfrag_type = 1 // mfrag
			elseif(stringmatch(FragType[i],"cfrag+")==1)
				Pfrag_type = 2 // cfrag + 1
			elseif(stringmatch(Fragtype[i],"constant")==1)
				Pfrag_type = 3 // constant
			elseif(stringmatch(Fragtype[i],"qfrag")==1)
				Pfrag_type = 4 // constant
			endif
		// select whether to include heterogenous chemistry
			if(het_wave[i]==1)
				hetchem = 1
				GammaOH = 1//gammaoh_wave[i]
			else
				hetchem = 0
			endif	
		// Set gas-phase wall loss rate
			gasWLR = WLRgas_wave[i]
			Cwall = Cwall_wave[i]
			kwall_gas_scaling = WLRg_scaling_wave[i]
		// sequential partitioning
			SPM = SPMwave[i]
		// kOH method
			krxn_method = kOHwave[i]
		// Set folder to save results into
			SaveFolderStr = "Run"+num2istr(i)
		// Set some other things, added 11/12/13
			SeedDiameter = DpSeed_wave[i]
			Nparticles = Np_wave[i]
			KineticMassTransfer = DynamicPartitioning_Wave[i]
			alpha = alpha_wave[i]
		// Do Fit
			DoWindow/F sompanel
			DoWindow/F FitGraph
			current = ""
			current = CompoundName[i] + "_" + FitTo[i] + "_" + FragType[i] + "_" + num2str(Het_Wave[i]) + "_#" + num2str(i)
			textbox/C/N=text0/A=MC current
			print current
			DoUpdate
		
			if(stringmatch(FragType[i],"cfrag"))
				Make/O/T/N=12 root:T_Constraints
				wave/T T_constraints = root:T_constraints
				T_Constraints[0] = {"K0 > 0.001","K0 < 0.999","K1 > .7","K1 < 2.5","K2 > 1e-3","K2 < 3","K3 > 1e-3","K3 < 3","K4 > 1e-3","K4 < 3","K5 > 1e-3","K5 < 3"}
			elseif(stringmatch(FragType[i],"mfrag"))
				Make/O/T/N=12 root:T_Constraints
				wave/T T_constraints = root:T_constraints
				T_Constraints[0] = {"K0 > 0.01","K0 < 5","K1 > .7","K1 < 2.5","K2 > 1e-3","K2 < 3","K3 > 1e-3","K3 < 3","K4 > 2e-3","K4 < 3","K5 > 1e-3","K5 < 3"}
			elseif(stringmatch(FragType[i],"qfrag"))
				Make/O/T/N=12 root:T_Constraints
				wave/T T_constraints = root:T_constraints
				T_Constraints[0] = {"K0 > 0.001","K0 < 100","K1 > .7","K1 < 2.5","K2 > 1e-3","K2 < 3","K3 > 1e-3","K3 < 3","K4 > 1e-3","K4 < 3","K5 > 1e-3","K5 < 3"}
			endif
			if(stringmatch(NOx[i],"O3"))
				Make/O/T/N=10 root:T_Constraints
				wave/T T_constraints = root:T_constraints
				T_Constraints[0] = {"K1 > .7","K1 < 2.5","K2 > 1e-3","K2 < 3","K3 > 1e-3","K3 < 3","K4 > 1e-3","K4 < 3","K5 > 1e-3","K5 < 3"}
			endif
			if(stringmatch(holdStrWave[i],"DLVP"))
				if(stringmatch(NOx[i],"O3"))
					deletepoints 0, 2, T_constraints
				else
					deletepoints 2, 2, T_constraints
				endif
			elseif(stringmatch(holdStrWave[i],"Pfrag"))
				if(stringmatch(NOx[i],"O3")!=1)
					deletepoints 0,2, T_constraints
				endif
			endif
			// select what type of fit to do
			if(stringmatch(FitTo[i],"Coa")==1)
				string HoldStr
				if(stringmatch(NOx[i],"O3")==1 && stringmatch(holdStrWave[i],"none"))
					HoldStr = "100000"
				elseif(stringmatch(NOx[i],"O3")==1 && stringmatch(holdStrWave[i],"DLVP"))
					HoldStr = "110000"
				elseif(stringmatch(NOx[i],"O3")!=1 && stringmatch(holdStrWave[i],"DLVP"))
					HoldStr = "010000"
				elseif(stringmatch(NOx[i],"O3")!=1 && stringmatch(holdStrWave[i],"Pfrag"))
					HoldStr = "100000"
				else
					HoldStr = "000000"
				endif

				wave Coa_experiment = $(df_Data+":Coa_experiment")
				wave Time_Experiment = $(df_Data+":Time_experiment")
				wave Coa_experiment_err = $(df_Data+":Coa_experiment_err")
				wave Coa_experiment_mask = $(df_Data+":Coa_experiment_mask")
				if(ErrorWave[i]==1 && MaskWave[i]==1)
					FuncFit/NTHR=0/H=holdStr allatonce_FitCoa cwave  Coa_experiment /X=Time_Experiment /W=Coa_experiment_err /M=Coa_experiment_mask /I=1 /D /C=T_Constraints 
				elseif(ErrorWave[i]==1 && MaskWave[i] == 0)
					FuncFit/NTHR=0/H=holdStr allatonce_FitCoa cwave  Coa_experiment /X=Time_Experiment /W=Coa_experiment_err /I=1 /D /C=T_Constraints 					
				elseif(ErrorWave[i]==0 && MaskWave[i]==1)
					FuncFit/NTHR=0/H=holdStr allatonce_FitCoa cwave  Coa_experiment /X=Time_Experiment /M=Coa_experiment_mask /I=1 /D /C=T_Constraints 					
				else
					FuncFit/NTHR=0/H=holdStr allatonce_FitCoa cwave  Coa_experiment /X=Time_Experiment /I=1 /D /C=T_Constraints 
				endif
			else
				print "you are obviously stupid. You can only fit to Coa with this function."
			endif
			
		// refresh these waves, just in case
			wave cwave = $(df_data+":coefwave")
			wave uncertainty = root:w_sigma
		
		if(V_fiterror!=0)
			cwave = 0
			V_fiterror = 0
		endif
			
	// Populate FitResults, etc. waves with results	
		wavestats/q/r=[2,5] cwave
		cwave[2,5]/=V_sum
		FitResults[i][][0] = cwave[q]
		FitResults[i][][1] = uncertainty[q]
		NVAR CS = root:ChiSq_Coa	
		ChiSq[i] = CS
		NVAR CS2 = root:ChiSq_O2C
		ChiSq_O2C[i] = nan
		runtime[i] = datetime
		RecalculateChi2_Batch(i,i)
//		GetStaticParams("2P",0,i,i)
//		GetStaticParams("VBS",0,i,i)
		DoUpdate
	endfor
	KillWaves/z T_constraints, coefwave, Coa_time_interp, deltaHC_time_interp, O2C_time_interp, w_sigma
End

//***************************************************************************************************
Function SOM_BatchFit_logCoaOnly(startpos,stoppos,dataset)
	// 04/18/2014 - written by CDC
	// Simplified version of SOM_BatchFit that assumes that one only wants to fit Coa, not Coa and O2C
	// thus, no need to access global fit parameters
	variable startpos, stoppos
	string dataset // name of dataset to use, "caltech", "UCR", ""

	string df_FitInfo = "root:FitInfo"
	string df_Global = df_FitInfo+":NewGlobalFitSetup"
	string df_data = "root"
	
	// waves that contain information for individual fitting
	setdatafolder $df_FitInfo
		wave/T CompoundName = $(df_FitInfo+":compoundname")	// name of compound to load and then fit
		wave/T NOx = $(df_FitInfo+":NOx")			// NOx condition of fit. "low", "high" or ""
		wave WLRgas_wave = $(df_FitInfo+":WLRgas_wave")	// gas-phase wall loss rate (1/s)
	
		wave/T RxnWave = $(df_FitInfo+":RxnWave")
		wave/T SaveWave = $(df_FitInfo+":SaveWave")
		wave/T FitTo = $(df_FitInfo+":FitTo")
		wave/T FragType = $(df_FitInfo+":FragType")
		wave Het_wave = $(df_FitInfo+":Het_Wave")
		wave Np_wave = $(df_FitInfo+":Np_wave")
		wave GammaOH_wave = $(df_FitInfo+":GammaOH_wave")
		wave InitialGuesses = $(df_FitInfo+":InitialGuesses")
		wave FitResults = $(df_FitInfo+":FitResults")
		wave O2Cerr_scale = $(df_FitInfo+":O2Cerr_scale")
		wave ChiSq = $(df_FitInfo+":ChiSq")
		wave ChiSq_O2C = $(df_FitInfo+":ChiSq_O2C")
		wave RunTime = $(df_FitInfo+":RunTime")
		wave Cwall_wave = $(df_FitInfo+":Cwall_wave")
		wave WLRg_scaling_wave = $(df_FitInfo+":WLRg_scaling_wave")
		wave SPMwave = $(df_FitInfo+":SPMwave")
		wave kOHwave = $(df_FitInfo+":kOH_wave")
		wave ErrorWave = $(df_FitInfo+":ErrorWave")
		wave MaskWave = $(df_FitInfo+":MaskWave")
		wave/T holdStrWave = $(df_FitInfo+":holdStrWave")
		wave DpSeed_wave = $(df_FitInfo+":DpSeed_wave")						// added 11/12/13
		wave DynamicPartitioning_wave = $(df_FitInfo+":DynamicPartitioning_wave")	// added 11/12/13
		wave alpha_wave = $(df_FitInfo+":alpha_wave")							// added 11/12/13
	
	// Fit Waves (coefficients and constraints) for "normal" fit
	setdatafolder root:
		make/o/d/n=(6) $(df_data+":coefwave")
		wave cwave = $(df_data+":coefwave")
		wave uncertainty = $(df_data+":w_sigma")
		Make/O/T/N=12 $(df_Data+":T_Constraints")
		wave/T T_constraints = $(df_Data+":T_Constraints")
		T_Constraints[0] = {"K0 > 0.01","K0 < 1","K1 > .7","K1 < 2.5","K2 > 1e-4","K2 < 1","K3 > 1e-4","K3 < 1","K4 > 1e-4","K4 < 1","K5 > 1e-4","K5 < 1"}
		
	// Global Variables to set
		NVAR Pfrag_type
		NVAR hetchem
		NVAR Nparticles
		NVAR GammaOH
		SVAR SaveFolderStr
		NVAR gasWLR	// gas-wall loss rate (1/s)
		NVAR Cwall		// effective wall OA concentration (mg/m^3); added 8/27/13
		NVAR kwall_gas_scaling	// composition dependent wall loss rate??; added 08/27/13
		NVAR timestep
		NVAR SPM
		NVAR krxn_method
		NVAR alpha					// added 11/12/13
		NVAR SeedDiameter			// added 11/12/13
		NVAR KineticMassTransfer	// added 11/12/13
		
	
	variable i, j
	string current
	variable V_fiterror = 0	// suppress abort on error
	variable V_fitMaxIters = 20
	variable reset = 0
	
	for(i=startpos;i<=stoppos;i+=1)
		// set reaction conditions
			if(stringmatch(dataset,"caltech"))
				SetFromCaltech(CompoundName[i],NOx[i])
			elseif(stringmatch(dataset,"UCR"))
				SetFromUCR(CompoundName[i],1)	// note that "compoundname" here is actually an experiment name; the "1" indicates use of wall-loss corrected data
			else
				abort "no dataset selected"
			endif
			if(reset==1)
				initialguesses[i][] = fitresults[i-1][q]
			endif
			setdatafolder root:
			if(stringmatch(NOx[i],"low"))
				timestep = 100
			elseif(stringmatch(NOx[i],"high"))
				timestep = 75
			elseif(stringmatch(dataset,"UCR"))
				timestep = 100
			else
				timestep = 100
			endif
		// set initial guesses for frag, dLVP, Pfunc
			for(j=0;j<6;j+=1)
				cwave[j] = InitialGuesses[i][j]
			endfor
		// select fragmentation parameterization
			if(stringmatch(FragType[i],"cfrag")==1)
				Pfrag_type = 0 // cfrag
			elseif(stringmatch(FragType[i],"mfrag")==1)
				Pfrag_type = 1 // mfrag
			elseif(stringmatch(FragType[i],"cfrag+")==1)
				Pfrag_type = 2 // cfrag + 1
			elseif(stringmatch(Fragtype[i],"constant")==1)
				Pfrag_type = 3 // constant
			elseif(stringmatch(Fragtype[i],"qfrag")==1)
				Pfrag_type = 4 // constant
			endif
		// select whether to include heterogenous chemistry
			if(het_wave[i]==1)
				hetchem = 1
				GammaOH = 1//gammaoh_wave[i]
			else
				hetchem = 0
			endif	
		// Set gas-phase wall loss rate
			gasWLR = WLRgas_wave[i]
			Cwall = Cwall_wave[i]
			kwall_gas_scaling = WLRg_scaling_wave[i]
		// sequential partitioning
			SPM = SPMwave[i]
		// kOH method
			krxn_method = kOHwave[i]
		// Set folder to save results into
			SaveFolderStr = "Run"+num2istr(i)
		// Set some other things, added 11/12/13
			SeedDiameter = DpSeed_wave[i]
			Nparticles = Np_wave[i]
			KineticMassTransfer = DynamicPartitioning_Wave[i]
			alpha = alpha_wave[i]
		// Do Fit
			DoWindow/F sompanel
			DoWindow/F FitGraph
			current = ""
			current = CompoundName[i] + "_" + FitTo[i] + "_" + FragType[i] + "_" + num2str(Het_Wave[i]) + "_#" + num2str(i)
			textbox/C/N=text0/A=MC current
			print current
			DoUpdate
		
			if(stringmatch(FragType[i],"cfrag"))
				Make/O/T/N=12 root:T_Constraints
				wave/T T_constraints = root:T_constraints
				T_Constraints[0] = {"K0 > 0.001","K0 < 0.999","K1 > .7","K1 < 2.5","K2 > 1e-3","K2 < 3","K3 > 1e-3","K3 < 3","K4 > 1e-3","K4 < 3","K5 > 1e-3","K5 < 3"}
			elseif(stringmatch(FragType[i],"mfrag"))
				Make/O/T/N=12 root:T_Constraints
				wave/T T_constraints = root:T_constraints
				T_Constraints[0] = {"K0 > 0.01","K0 < 5","K1 > .7","K1 < 2.5","K2 > 1e-3","K2 < 3","K3 > 1e-3","K3 < 3","K4 > 2e-3","K4 < 3","K5 > 1e-3","K5 < 3"}
			elseif(stringmatch(FragType[i],"qfrag"))
				Make/O/T/N=12 root:T_Constraints
				wave/T T_constraints = root:T_constraints
				T_Constraints[0] = {"K0 > 0.001","K0 < 100","K1 > .7","K1 < 2.5","K2 > 1e-3","K2 < 3","K3 > 1e-3","K3 < 3","K4 > 1e-3","K4 < 3","K5 > 1e-3","K5 < 3"}
			endif
			if(stringmatch(NOx[i],"O3"))
				Make/O/T/N=10 root:T_Constraints
				wave/T T_constraints = root:T_constraints
				T_Constraints[0] = {"K1 > .7","K1 < 2.5","K2 > 1e-3","K2 < 3","K3 > 1e-3","K3 < 3","K4 > 1e-3","K4 < 3","K5 > 1e-3","K5 < 3"}
			endif
			if(stringmatch(holdStrWave[i],"DLVP"))
				if(stringmatch(NOx[i],"O3"))
					deletepoints 0, 2, T_constraints
				else
					deletepoints 2, 2, T_constraints
				endif
			elseif(stringmatch(holdStrWave[i],"Pfrag"))
				if(stringmatch(NOx[i],"O3")!=1)
					deletepoints 0,2, T_constraints
				endif
			endif
			// select what type of fit to do
			if(stringmatch(FitTo[i],"Coa")==1)
				string HoldStr
				if(stringmatch(NOx[i],"O3")==1 && stringmatch(holdStrWave[i],"none"))
					HoldStr = "100000"
				elseif(stringmatch(NOx[i],"O3")==1 && stringmatch(holdStrWave[i],"DLVP"))
					HoldStr = "110000"
				elseif(stringmatch(NOx[i],"O3")!=1 && stringmatch(holdStrWave[i],"DLVP"))
					HoldStr = "010000"
				elseif(stringmatch(NOx[i],"O3")!=1 && stringmatch(holdStrWave[i],"Pfrag"))
					HoldStr = "100000"
				else
					HoldStr = "000000"
				endif

				wave Coa_experiment = $(df_Data+":Coa_experiment")
				wave Time_Experiment = $(df_Data+":Time_experiment")
				wave Coa_experiment_err = $(df_Data+":Coa_experiment_err")
				wave Coa_experiment_mask = $(df_Data+":Coa_experiment_mask")
				// get log(Coa)
				wave logCoa = $(df_Data+":logOffset_Coa_exp")
				wave logCoa_err = $(df_Data+":logOffset_Coa_unc")
				if(waveexists(logCoa))
					redimension/n=(numpnts(Coa_experiment)) logCoa, logCoa_err
					logCoa = nan
					logCoa_err = nan
				else
					make/o/d/n=(numpnts(Coa_experiment)) $(df_Data+":logOffset_Coa_exp")
					wave logCoa = $(df_Data+":logOffset_Coa_exp")
					make/o/d/n=(numpnts(Coa_experiment)) $(df_Data+":logOffset_Coa_unc")
					wave logCoa_err = $(df_Data+":logOffset_Coa_unc")
				endif
					variable myoffset = 2
					logCoa = log(Coa_experiment+myoffset)
					logCoa_err = Coa_experiment_err/(Coa_experiment+myoffset)
				
				if(ErrorWave[i]==1 && MaskWave[i]==1)
					FuncFit/NTHR=0/H=holdStr allatonce_FitLogOffsetCoa cwave  logCoa /X=Time_Experiment /W=logCoa_err /M=Coa_experiment_mask /I=1 /D /C=T_Constraints 
				elseif(ErrorWave[i]==1 && MaskWave[i] == 0)
					FuncFit/NTHR=0/H=holdStr allatonce_FitLogOffsetCoa cwave  logCoa /X=Time_Experiment /W=logCoa_err /I=1 /D /C=T_Constraints 					
				elseif(ErrorWave[i]==0 && MaskWave[i]==1)
					FuncFit/NTHR=0/H=holdStr allatonce_FitLogOffsetCoa cwave  logCoa /X=Time_Experiment /M=Coa_experiment_mask /I=1 /D /C=T_Constraints 					
				else
					FuncFit/NTHR=0/H=holdStr allatonce_FitLogOffsetCoa cwave  logCoa /X=Time_Experiment /I=1 /D /C=T_Constraints 
				endif
			else
				print "you are obviously stupid. You can only fit to Coa with this function."
			endif
			
		// refresh these waves, just in case
			wave cwave = $(df_data+":coefwave")
			wave uncertainty = root:w_sigma
		
		if(V_fiterror!=0)
			cwave = 0
			V_fiterror = 0
		endif
			
	// Populate FitResults, etc. waves with results	
		wavestats/q/r=[2,5] cwave
		cwave[2,5]/=V_sum
		FitResults[i][][0] = cwave[q]
		FitResults[i][][1] = uncertainty[q]
		NVAR CS = root:ChiSq_Coa	
		ChiSq[i] = CS
		NVAR CS2 = root:ChiSq_O2C
		ChiSq_O2C[i] = nan
		runtime[i] = datetime
		RecalculateChi2_Batch(i,i)
//		GetStaticParams("2P",0,i,i)
//		GetStaticParams("VBS",0,i,i)
		DoUpdate
	endfor
	KillWaves/z T_constraints, coefwave, Coa_time_interp, deltaHC_time_interp, O2C_time_interp, w_sigma
End

//****************************************************************************
function makenewmovie(movienamestr,wvstr)
	string movienamestr
	string wvstr
	
	wave scalewave = Coa_time
	Display;AppendImage $wvstr;DelayUpdate
	wavestats/q $wvstr
	wave wv = $wvstr
	ModifyImage $wvstr ctab= {wv[0][0][0],*,Rainbow,1};DelayUpdate
	ModifyImage $wvstr ctab= {wv[0][0][0],*,Rainbow,1};DelayUpdate
//	ModifyImage $wvstr ctab= {V_max*1e-6,V_max/2,Rainbow,1};DelayUpdate
//	ModifyImage $wvstr ctab= {1e-3,1000,Rainbow,1};DelayUpdate
	ModifyImage $wvstr minRGB=(43520,43520,43520),maxRGB=0
	ModifyImage $wvstr log=1
	setaxis bottom, 1.5,12.5
	setaxis left, -0.5, 7.5
	
	variable nRows = dimsize($wvstr,0)
	variable nCols = dimsize($wvstr,1)
	make/o/d/n=(nRows,nCols)/FREE moviematrix
	variable n_frames 
	variable m
	n_frames = dimsize($wvstr,2)
	
	
	NewPath/O pathname, "C::Users:Chris Cappa:Documents:Research:SAPRC-SOM:GECKO:movies"
//	NewPath/O pathname, "C::Users:Chris Cappa:Documents:Research:SAPRC-SOM:Caltech Data:movies"
		
	NewMovie/P=pathname/O/Z as movienamestr
	string current
	
	for(m=1;m<=n_frames;m+=1)
		moviematrix = wv[p][q][m-1]
		moviematrix[][8,] = 0
		wavestats/q moviematrix
		ModifyImage $wvstr ctab= {V_max*5e-4,V_max*0.5,Rainbow,1};DelayUpdate
//		ModifyImage $wvstr ctab= {wv[0][0][0]*1e-4,wv[0][0][0]*0.05,Rainbow,1};DelayUpdate
//		ModifyImage $wvstr ctab= {scalewave[m-1]*1e-4,scalewave[m-1]*0.2,Rainbow,1};DelayUpdate
		ModifyImage $wvstr plane= m-1; DoUpdate
		AddMovieFrame
	endfor
	
	CloseMovie	
end

//****************************************************************************
function makedualmovie(movienamestr,wvstr1,wvstr2)
	string movienamestr
	string wvstr1
	string wvstr2
	
	wave scalewave = Coa_time
	wave wv1 = $wvstr1
	wave wv2 = $wvstr2

	Display;AppendImage $wvstr1;DelayUpdate
	ModifyImage $wvstr1 ctab= {*,*,Rainbow,1};DelayUpdate
	ModifyImage $wvstr1 plane= 364
	AppendImage/L=L2/B=B2 $wvstr2;DelayUpdate
	ModifyImage $wvstr2 ctab= {*,*,Rainbow,1};DelayUpdate
	ModifyImage $wvstr2 plane= 364
	ModifyGraph lblPos(B2)=37,axisEnab(bottom)={0,0.47},axisEnab(B2)={0.53,1};DelayUpdate
	ModifyGraph freePos(L2)={0,B2},freePos(B2)=0;DelayUpdate
	Label left "Number of Oxygens\\u#2";DelayUpdate
	Label bottom "Number of Carbons\\u#2";DelayUpdate
	Label L2 "\\u#2";DelayUpdate
	Label B2 "Number of Carbons\\u#2"

	ModifyImage $wvstr1 minRGB=(43520,43520,43520),maxRGB=0
	ModifyImage $wvstr1 log=1
	ModifyImage $wvstr2 minRGB=(43520,43520,43520),maxRGB=0
	ModifyImage $wvstr2 log=1
	setaxis bottom, 1.5,12.5
	setaxis left, -0.5, 7.5
	setaxis B2, 1.5,12.5
	setaxis L2, -0.5, 7.5
	ModifyGraph width=576, height = 252
	
	variable nRows = dimsize($wvstr1,0)
	variable nCols = dimsize($wvstr1,1)
	make/o/d/n=(nRows,nCols)/FREE moviematrix
	variable n_frames 
	variable m
	n_frames = dimsize($wvstr1,2)
	
	
	NewPath/O pathname, "C::Users:Chris Cappa:Documents:Research:SAPRC-SOM:GECKO:movies"
////	NewPath/O pathname, "C::Users:Chris Cappa:Documents:Research:SAPRC-SOM:Caltech Data:movies"
//		
	NewMovie/P=pathname/O/Z as movienamestr
	string current
//	
	for(m=1;m<=n_frames;m+=1)
		moviematrix = wv1[p][q][m-1]
		wavestats/q moviematrix
		ModifyImage $wvstr1 ctab= {V_max*1e-4,V_max*0.5,Rainbow,1};DelayUpdate
		moviematrix = wv2[p][q][m-1]
		wavestats/q moviematrix
		ModifyImage $wvstr2 ctab= {V_max*1e-4,V_max*0.5,Rainbow,1};DelayUpdate

		ModifyImage $wvstr1 plane= m-1; DoUpdate
		ModifyImage $wvstr2 plane= m-1; DoUpdate
		AddMovieFrame
	endfor
	
	CloseMovie	
end

//**************************************************************
// Calculate the sum of all rows/columns in a given layer for a 3D wave having n layers
// Return result as a 1D wave, with n rows
function SumMatrixByLayers(mat,newwavestr)
	wave mat
	string newwavestr
	
	variable nRows = dimsize(mat,0)
	variable nColumns = dimsize(mat,1)
	variable nLayers = dimsize(mat,2)
	
	make/o/d/n=(nRows,nColumns)/FREE SingleMat
	make/o/d/n=(nLayers) $newWaveStr
	wave myWave = $newWaveStr
	note/k MyWave "Sum of each layer from wave" 
	variable i
	for(i=0;i<nLayers;i+=1)
		SingleMat = mat[p][q][i]
		wavestats/q SingleMat
		myWave[i] = V_sum
	endfor
end

//*************************************************************
function SOM_SetParamsFromFit(fitnum)
	// function to set parameters for a run based on results of fitting
	variable fitnum // index of best fit from wave "FitResults" in folder "FitInfo"
	
	setdatafolder root:
	NVAR FragSlope // either mfrag or cfrag, depending on fragmentation method
	NVAR delta_logCstar_perO // change in logCstar per oxygen added
	NVAR ProbOx1, ProbOx2, ProbOx3, ProbOx4	// probability of adding some number of oxygen atoms
	NVAR gasWLR			// gas phase wall loss first order rate coefficient (1/s)
	NVAR Cwall				// wall equivalent OA concentration, mg/m^3
	NVAR kwall_gas_scaling	// scaling factor for gasWLR
	NVAR Pfrag_type			// Pfrag = ?
	NVAR SmallFragments	// Random, equal or small fragments
	NVAR SPM				// sequential partitioning model
	NVAR krxn_method		// method to calculate kOH matrix
	NVAR alpha				// mass accommodation coefficient for dynamic partitioning
	NVAR KineticMassTransfer // 0 = instantaneous equilibrium, 1 = dynamic partitioning
	NVAR SeedDiameter		// seed diameter
	NVAR Nparticles			// seed particle number concentration
	
	string fitfolder = "root:fitinfo:"
	wave FitResults = $(fitfolder+"FitResults")
	wave Cwall_wave = $(fitfolder+"Cwall_wave")
	wave WLRg_scaling_wave = $(fitfolder+"WLRg_scaling_wave")
	wave WLRgas_wave = $(fitfolder+"WLRgas_wave")
	wave/T FragType = $(fitfolder+"FragType")
	wave SPMwave = $(fitfolder+"SPMwave")
	wave kOHwave = $(fitfolder+"kOH_wave")
	wave alpha_wave = $(fitfolder+"alpha_wave")
	wave DynamicPartitioning_wave = $(fitfolder+"dynamicpartitioning_wave")
	wave Np_wave = $(fitfolder+"Np_wave")
	wave DpSeed_wave = $(fitfolder+"DpSeed_wave")
	
	FragSlope = FitResults[fitnum][0]
	delta_logCstar_perO = FitResults[fitnum][1]
	ProbOx1 = FitResults[fitnum][2]
	ProbOx2 = FitResults[fitnum][3]
	ProbOx3 = FitResults[fitnum][4]
	ProbOx4 = FitResults[fitnum][5]
	gasWLR = WLRgas_wave[fitnum]
	Cwall = Cwall_wave[fitnum]
	gasWLR = WLRgas_wave[fitnum]
	SPM = SPMwave[fitnum]
	krxn_method = kOHwave[fitnum]
	alpha = alpha_wave[fitnum]
	kineticmasstransfer = dynamicpartitioning_wave[fitnum]
	Nparticles = Np_wave[fitnum]
	SeedDiameter = DpSeed_wave[fitnum]
	
	if(stringmatch(FragType[fitnum],"mfrag"))
		Pfrag_type = 1
	elseif(stringmatch(FragType[fitnum],"cfrag"))
		Pfrag_type = 0
	endif	
End

//******************************************************************
Function CalcGasWallLossBias(kwall,kwall_scaling,Cwall_var)
	// Function to calculate the time-dependent "wall loss bias" based on given inputs and current SOM params
	// Needs to be cleaned up and commented
	variable kwall		// wall loss rate coefficient (1/s)
	variable kwall_scaling	// scaling factor for wall loss rate coefficient
	variable Cwall_var	// effective wall OA concentration (mg/m^3)
	
	setdatafolder root:
	
	NVAR gasWLR
	NVAR Cwall
	NVAR kwall_gas_scaling
	
//	gasWLR = kwall
//	Cwall = Cwall_var
//	kwall_gas_scaling = kwall_scaling
	gasWLR = 0
	kwall_gas_scaling = 0
	
	SOM_v1(quiet=1)
	wave Coa_time
	wave O2C_time
	wave TimeW
	variable nRows = numpnts(Coa_time)
	variable deltatime = TimeW[1]-TimeW[0]
	make/o/d/n=(nRows)/FREE Coa_wall = Coa_time
	make/o/d/n=(nRows)/FREE O2C_wall = O2C_time
	make/o/d/n=(nRows) R_soa
	setscale/P x, 0, deltaTime, "hours", R_soa
	note R_soa "Ratio between Coa with no walls and with walls"
	make/o/d/n=(nRows) R_O2C
	setscale/P x, 0, deltaTime, "unitless", R_O2C
	note R_O2C "Difference between O2C with no walls and with walls"
	
//	gasWLR = 0
//	kwall_gas_scaling = 0
	gasWLR = kwall
	Cwall = Cwall_var
	kwall_gas_scaling = kwall_scaling
	
	SOM_v1(quiet=1)
	wave Coa_time
	wave O2C_time
//	R_soa = Coa_time/Coa_wall
//	R_O2C = O2C_time - O2C_wall
	R_soa = Coa_wall/Coa_time
	R_O2C = O2C_wall - O2C_time

End

//**********************************************************************
// not important
Function BatchRun_OH()

	make/o/d/n=(21) root:FitInfo:OHconc
	wave OHconcWave = root:FitInfo:OHconc
	OHconcWave = 1e6*1.28^x
	make/o/d/n=(21) root:FitInfo:VOCconc
	wave VOCconcWave = root:FitInfo:VOCconc
	VOCconcWave = (1*1.275^x)/1000
	
	make/o/d/n=(20,20) root:FitInfo:Rsoa_OH, root:FitInfo:Csoa_OH
	wave Rsoa_OH = root:FitInfo:Rsoa_OH
	Rsoa_OH = nan
	wave Csoa_OH = root:FitInfo:Csoa_OH
	Csoa_OH = nan
	
	NVAR OHconc
	NVAR UseScalingForOxidant
	UseScalingForOxidant = 0
	NVAR ctot_ppm
	
	variable i, j
	for(i=0;i<numpnts(OHconcWave)-1;i+=1)
		for(j=0;j<numpnts(VOCconcWave)-1;j+=1)
			OHconc = OHconcWave[i]
			Ctot_ppm = VOCconcWave[j]
			CalcGasWallLossBias(1e-4,0,10)
			wave Coa_time
			wave R_soa
			findlevel/P/Q Coa_time, 0.5
			wavestats/q Coa_time
			Csoa_OH[i][j] = V_max
			if(V_max > 0.5)
				wavestats/q/r=[V_levelX,] R_soa
				Rsoa_OH[i][j] = V_avg
			endif
//			Rsoa_OH[i][j] = V_sdev
			//SOM_v1()
		endfor
	endfor
//		rsoa_oh[][][2] = rsoa_oh[p][q][0]+rsoa_oh[p][q][1]
//		rsoa_oh[][][3] = rsoa_oh[p][q][0]-rsoa_oh[p][q][1]
	
End

//***********************************************************************
// not important
Function BatchRun_WallBias(start,stop)
	variable start
	variable stop
	
	setdatafolder root:
	variable i
	variable kwall_current
	wave R_soa_avg = root:fitinfo:R_soa_avg
	wave R_soa_sdev = root:fitinfo:R_soa_sdev
	wave/t compoundname = root:fitinfo:compoundname
	wave/t nox = root:fitinfo:nox
	wave maskwave = root:fitinfo:maskwave
	variable method=0
		
	for(i=start;i<=stop;i+=1)
		setfromcaltech(compoundname[i],nox[i])
		SOM_SetParamsFromFit(i)
		NVAR gasWLR
		NVAR Cwall
		NVAR kwall_gas_scaling
		CalcGasWallLossBias(gasWLR,kwall_gas_scaling,Cwall)
		wave R_soa = root:R_soa
		wave Coa_time
		findlevel/P/Q Coa_time, 0.5
		method = maskwave[i]
		if(method==0)
			wavestats/q/r=[V_levelX,] R_soa
		else
			// trying something new
			variable first = V_levelX
			wave Coa_experiment_mask, Coa_experiment
			wavestats/q Coa_experiment_mask
			findlevel/P/Q Coa_time, Coa_experiment[V_sum]
			variable last = V_levelX
			wavestats/q/r=[first,last] R_soa
			print first,last
			//
//	
		endif
		R_soa_avg[i] = V_avg
		R_soa_sdev[i] = V_sdev
	endfor
End

//***********************************************************************
// not important
Function GasWallLoss_HCconc()

	wave HCconc = root:GasWallLoss:HCconc_ppb
	wave Rsoa = root:GasWallLoss:R_SOA_avg_HCconc
	
	setdatafolder root:
	NVAR Ctot_ppm
	variable i
	variable nConc = numpnts(HCconc)
	for(i=0;i<nConc;i+=1)
		Ctot_ppm = HCconc[i]/1000
		CalcGasWallLossBias(1e-4,0,10)
		wave Coa_time
		wave R_SOA
		findlevel/Q/P Coa_time, 0.5
		wavestats/q/r=[V_levelX,] r_soa
		if(V_flag==0)
			Rsoa[i] = V_avg
		else
			Rsoa[i] = nan
		endif
	endfor
End

//**************************************************************************
Function RecalculateChi2(FitNum)
	variable FitNum
	
	string df_fitinfo = "root:fitinfo"
	string df_data = "root"
	SOM_SetParamsFromFit(FitNum)
	wave/T CompoundName = $(df_fitinfo+":CompoundName")
	wave/T NOx = $(df_fitinfo+":NOx")
	wave ChiSq = $(df_fitinfo+":ChiSq")
	SetFromCaltech(CompoundName[FitNum],NOx[FitNum])
	DoWindow/F sompanel
	DoWindow/F FitGraph
	wave Time_experiment = $(df_data+":Time_Experiment")
	som_v1(quiet=1,fittime=Time_Experiment)
	wave Coa_time = Coa_time_interp
	wave Coa_experiment
	wave Coa_experiment_err
	make/o/d/n=(numpnts(Coa_time))/FREE chi2wave
	chi2wave = (coa_experiment-coa_time)^2/coa_experiment_err^2
	wavestats/q chi2wave
	ChiSq[FitNum] = (1/(numpnts(Coa_time)-6-1))*V_sum
	
End

//**************************************************************************
Function RecalculateChi2_single()
	// assumes that you have already set the experiment that you want and run SOM
	// requires good uncertainty estimates (Coa_experiment_err) for this to work well
	
	setdatafolder root:
	wave Time_experiment // The experimental timewave
	wave Coa_time = Coa_time_interp // the SOM timewave
	wave Coa_experiment // the observed Coa time-series
	wave Coa_experiment_err // the errors associated with the observed Coa time-series
	wave Coa_experiment_mask // mask for the observations; 0 = mask; 1 = okay
	make/o/d/n=(numpnts(Coa_time))/FREE chi2wave
	chi2wave = (coa_experiment-coa_time)^2/coa_experiment_err^2
	chi2wave *= Coa_experiment_mask // set masked values to zero
	wavestats/q chi2wave
	variable ChiSq_var
	ChiSq_var = (1/(numpnts(Coa_time)-6-1))*V_sum
	Return ChiSq_var
	
End

//**************************************************************************
// Recalculate chi-square values after doing fits
// Sometimes necessary b/c this isn't always calculated correctly the first time
// start/stop relate to indices of fits
Function RecalculateChi2_Batch(start,stop)
	variable start // 
	variable stop
	
	variable i
	
	string df_fitinfo = "root:fitinfo"
	string df_data = "root"
	wave/T CompoundName = $(df_fitinfo+":CompoundName")
	wave/T NOx = $(df_fitinfo+":NOx")
	wave ChiSq = $(df_fitinfo+":ChiSq")
	wave MaskWave = $(df_fitInfo+":MaskWave")
	wave Np_wave = $(df_fitinfo+":Np_wave")
	wave alpha_wave = $(df_fitinfo+":alpha_wave")
	wave DynamicPartitioning_wave = $(df_FitInfo+":DynamicPartitioning_wave")
	wave DpSeed_wave = $(df_FitInfo+":DpSeed_wave")
	variable npnts
	
	for(i=start;i<=stop;i+=1)
		SOM_SetParamsFromFit(i)		
		SetFromCaltech(CompoundName[i],NOx[i])
//		DoWindow/F sompanel
		DoWindow/F FitGraph
		//
		NVAR KineticMassTransfer
		KineticMassTransfer = DynamicPartitioning_wave[i]
		NVAR alpha
		alpha = alpha_wave[i]
		NVAR Nparticles
		Nparticles = Np_wave[i]
		NVAR SeedDiameter
		SeedDiameter = DpSeed_wave[i]
		//
		wave Time_Experiment = $(df_data+":Time_Experiment")
		som_v1(quiet=1,fittime=Time_Experiment)
		wave Coa_time = Coa_time_interp
		wave Coa_experiment
		wave Coa_experiment_err
		npnts = numpnts(Coa_time)
		if(MaskWave[i]==1)
			wave Coa_experiment_mask
			wavestats/q Coa_experiment_mask
			npnts = V_Sum
		endif
		make/o/d/n=(numpnts(Coa_time))/FREE chi2wave
		Coa_experiment_err = max(0.5,0.15*coa_experiment)
		chi2wave = (coa_experiment-coa_time)^2/coa_experiment_err^2
		wavestats/q/r=[0,npnts] chi2wave
		ChiSq[i] = (1/(nPnts-6-1))*V_sum
	endfor
End

//********************************************************************************
// Fit yield vs. deltaHC data to determine VBS or 2-product parameters
// Lots of things hardcoded, so not very flexible
// determine separate fits for original data and for "vapor wall loss" corrected data
Function GetStaticParams(VBSor2P,start,stop)
	string VBSor2P	// "2P" or "VBS"
	variable start
	variable stop
	
	variable i
	
	string df_2P = "root:VBSand2P"
	string df_fitinfo = "root:fitinfo"
	
	wave/T CompoundName = $(df_fitinfo+":CompoundName")
	wave/T NOx = $(df_fitinfo+":NOx")
	wave ChiSq = $(df_fitinfo+":ChiSq")
	wave MaskWave = $(df_fitInfo+":MaskWave")
	wave WLRgas_wave = $(df_fitInfo+":WLRgas_wave")
	variable nRows = numpnts(WLRgas_wave)
	make/o/d/n=(4) $(df_2P+":CoefWave")
	wave CoefWave = $(df_2P+":CoefWave")
	make/o/t/n=(8) $(df_2P+":T_constraints")
	wave/T T_constraints = $(df_2P+":T_constraints")
	if(stringmatch(VBSor2P,"2P"))
		// 2 product
		make/o/d/n=(nRows,5,2) $(df_2P+":TwoProductFits")
		wave FitResultsNoWalls = $(df_2P+":TwoProductFits")
		FitResultsNoWalls = nan
		make/o/d/n=(nRows,5,2) $(df_2P+":TwoProductFitsWalls")
		wave FitResultsWalls = $(df_2P+":TwoProductFitsWalls")
		FitResultsWalls = nan
		T_Constraints[0] = {"K0 > 1e-3","K0 < 5","K1 > 1e-3","K1 < 1000","K2 > 1e-3","K2 < 5","K3 > 1e-3","K3 < 10000"}
		CoefWave[0] = {0.1, 5, 1, 100}
	else
		// VBS
		make/o/d/n=(nRows,5,2) $(df_2P+":VBSfits")
		wave FitResultsNoWalls = $(df_2P+":VBSfits")
		FitResultsNoWalls = nan
		make/o/d/n=(nRows,5,2) $(df_2P+":VBSfitsWalls")
		wave FitResultsWalls = $(df_2P+":VBSfitsWalls")
		FitResultsWalls = nan
		T_Constraints[0] = {"K0 > 1e-3","K0 < 5","K1 > 1e-3","K1 < 5","K2 > 1e-3","K2 < 5","K3 > 1e-3","K3 < 5"}
		CoefWave[0] = {0.6,0.4,0.1,0.05}
	endif
	variable npnts
	variable V_FitMaxIters = 80
//	variable V_fitErrror = 0
	
	for(i=start;i<=stop;i+=1)
		setdatafolder root:
		SOM_SetParamsFromFit(i)		
		SetFromCaltech(CompoundName[i],NOx[i])
		if(WLRgas_wave[i]>0)
			// First do calculation with wall loss included
			DoWindow/F FitGraph
			som_v1(quiet=1)
			wave Coa_time	// in ug/m^3
			wave deltaHC_time	// in ug/m^3
			duplicate/o Coa_time $(df_2P+":Coa_time")
			wave deltaM = $(df_2P+":Coa_time")
			duplicate/o deltaHC_time $(df_2P+":deltaHC_time")
			wave deltaHC = $(df_2P+":deltaHC_time")
			setdatafolder $df_2P
			wave/T T_constraints		
			if(stringmatch(VBSor2P,"2P"))
				FuncFit/NTHR=0/N/Q Fit2Product CoefWave  deltaM /X=deltaHC /D /C=T_Constraints
			else
				FuncFit/NTHR=0/N/Q FitVBS CoefWave  deltaM /X=deltaHC /D /C=T_Constraints
			endif
			wave CoefWave, W_sigma
			FitResultsWalls[i][][0] = CoefWave[q]
			FitResultsWalls[i][][1] = W_sigma[q]
			FitResultsWalls[i][4][0] = V_chisq/(numpnts(deltaM)-6-1)
			// next do the calculation where there are no walls
			setdatafolder root:
			NVAR gasWLR
			gasWLR = 0
			som_v1(quiet=1)
			wave Coa_time	// in ug/m^3
			wave deltaHC_time	// in ug/m^3
			duplicate/o Coa_time $(df_2P+":Coa_time")
			wave deltaM = $(df_2P+":Coa_time")
			duplicate/o deltaHC_time $(df_2P+":deltaHC_time")
			wave deltaHC = $(df_2P+":deltaHC_time")
			setdatafolder $df_2P
			wave CoefWave
			wave/T T_constraints
			if(stringmatch(VBSor2P,"2P"))
				FuncFit/NTHR=0/N/Q Fit2Product CoefWave  deltaM /X=deltaHC /D /C=T_Constraints
			else
				FuncFit/NTHR=0/N/Q FitVBS CoefWave  deltaM /X=deltaHC /D /C=T_Constraints
			endif			
			wave CoefWave, W_sigma
			FitResultsNoWalls[i][][0] = CoefWave[q]
			FitResultsNoWalls[i][][1] = W_sigma[q]
			FitResultsNoWalls[i][4][0] = V_chisq/(numpnts(deltaM)-6-1)
		endif
	endfor
End


//************************************************************************************
Function RunAbunchOfThings(start,stop,low_high)
// This will run calculations for data in the Toluene2013 folder from Caltech
// Based on the best fit values of ones favorite fit, Coa and Yields will be calculated for each experiment
// Useful for recalculating things after fits have been been performed, including determining the influence of vapor wall loss
// Too many things are hardwired, and this is not sufficiently flexible
// The results will be saved in a folder.
	variable start
	variable stop
	string low_high	// either "low" or "high"
	
	variable i, j
	variable nexpts = 6
	string myexp
		
	string df_fitinfo = "root:fitinfo"
	wave/T CompoundName = $(df_fitinfo+":CompoundName")
	wave/T NOx = $(df_fitinfo+":NOx")
	wave ChiSq = $(df_fitinfo+":ChiSq")
	wave MaskWave = $(df_fitInfo+":MaskWave")
	variable npnts
	
	string df_results = "root:FitResults"
	string df_results_FitNum
	string df_results_fitnum_expnum
	wave ChiSq_Summary = $(df_results+":ChiSq_Summary")
	wave Yield_Summary = $(df_results+":Yield_Summary")
	wave Rwall_Summary = $(df_results+":Rwall_Summary")
	make/o/d/n=6/FREE DpSeed = {5,107,157,210,275,251}
	NVAR SeedDiameter = root:SeedDiameter
	
	for(i=start;i<=stop;i+=1)
		for(j=0;j<nexpts;j+=1)
			myexp = low_high+num2istr(j)
			df_results_fitnum = df_results+":Run"+num2istr(i)				// a new folder
			df_results_fitnum_expnum = df_results_fitnum + ":" + myexp	// another new folder
			
			SOM_SetParamsFromFit(i)		
			SetFromCaltech(CompoundName[i],myexp)
			SeedDiameter = DpSeed[j]
	//		DoWindow/F sompanel
			DoWindow/F FitGraph
			wave Time_experiment = root:Time_Experiment
			som_v1(quiet=1,fittime=Time_Experiment)
			wave Coa_time = Coa_time_interp
			wave Coa_experiment
			wave Coa_experiment_err
			npnts = numpnts(Coa_time)
			if(MaskWave[i]==1)
				wave Coa_experiment_mask
				wavestats/q Coa_experiment_mask
				npnts = V_Sum
			endif
			make/o/d/n=(numpnts(Coa_time))/FREE chi2wave
			Coa_experiment_err = max(0.5,0.15*coa_experiment)
			chi2wave = 0
			chi2wave = (coa_experiment-coa_time)^2/coa_experiment_err^2
			
			wave Coa_time = Coa_time
			wave Yield_time = Yield_time
			if(datafolderexists(df_results_fitnum)==0)
				newdatafolder $df_results_fitnum
				make/o/d/n=(nexpts) $df_results_fitnum+":ChiSq_Coa"
				make/o/d/n=(nexpts,2) $df_results_fitnum+":Rwall"
				make/o/d/n=(nexpts,2) $df_results_fitnum+":Yield_Final"
			endif
			wave ChiSq_Coa = $df_results_fitnum+":ChiSq_Coa"
			wave Rwall = $df_results_fitnum+":Rwall"
			wave Yield_final = $df_results_fitnum+":Yield_Final"
			wavestats/q/r=[0,npnts] chi2wave
			//ChiSq[i] = (1/(nPnts-6-1))*V_sum
			ChiSq_Summary[i][j] = (1/(nPnts-6-1))*V_sum
			ChiSq_Coa[j] = (1/(nPnts-6-1))*V_sum
			Rwall[j][0] = j
			Yield_Final[j][0] = Yield_time[numpnts(Yield_time)-1]*(100/92)
			Yield_Summary[i][j][0] = Yield_time[numpnts(Yield_time)-1]*(100/92)
			newdatafolder/O $df_results_fitnum_expnum
			duplicate/O Coa_time, $(df_results_fitnum_expnum+":Coa_time")
			duplicate/O Yield_time, $(df_results_fitnum_expnum+":Yield_time")
			wave Coa_wall = $(df_results_fitnum_expnum+":Coa_time")
			DoUpdate
			
			// now calculate for kwall = 0
			NVAR kwall = root:gasWLR
			kwall = 0
			som_v1(quiet=1,fittime=Time_Experiment)
			wave Coa_time = root:Coa_time
			wave Yield_time = root:Yield_Time
			duplicate/O Coa_time, Rwall_time
			Rwall_time = Coa_time/Coa_wall
			duplicate/O Coa_time, $(df_results_fitnum_expnum+":Coa_time_nowall")
			duplicate/O Yield_time, $(df_results_fitnum_expnum+":Yield_time_nowall")
			duplicate/O Rwall_time, $(df_results_fitnum_expnum+":Rwall_time_nowall")
			findlevel/P/Q Coa_wall, 0.5
			wavestats/q/r=[V_levelX,] Rwall_time
			//wavestats/q Rwall_time
			Rwall[j][0] = V_avg
			Rwall[j][1] = V_sdev
			Rwall_summary[i][j][0] = V_avg
			Rwall_summary[i][j][1] = V_sdev
			Yield_Final[j][1] = Yield_time[numpnts(Yield_time)-1]*(100/92)
			Yield_Summary[i][j][1] = Yield_time[numpnts(Yield_time)-1]*(100/92)
			DoUpdate
		endfor
	endfor
End

Window Graph2() : Graph
	PauseUpdate; Silent 1		// building window...
	String fldrSav0= GetDataFolder(1)
	SetDataFolder root:GasWallLoss:dodecane_example:
	Display /W=(420.75,136.25,774,416) :::Caltech_Data:Alkanes:C12:Low_NOx:Dodecane:Coa_experiment vs :::Caltech_Data:Alkanes:C12:Low_NOx:Dodecane:Time_Experiment
	AppendToGraph Coa_time_wl,Coa_time_noWL
	AppendToGraph/R R_soa
	SetDataFolder fldrSav0
	ModifyGraph margin(top)=7
	ModifyGraph mode(Coa_experiment)=3
	ModifyGraph marker(Coa_experiment)=19
	ModifyGraph lSize(Coa_time_wl)=2,lSize(Coa_time_noWL)=2,lSize(R_soa)=2
	ModifyGraph lStyle(Coa_time_wl)=3,lStyle(R_soa)=7
	ModifyGraph rgb(Coa_experiment)=(0,0,0),rgb(Coa_time_wl)=(65280,43520,0),rgb(Coa_time_noWL)=(0,0,63232)
	ModifyGraph zColor(R_soa)={:GasWallLoss:dodecane_example:Coa_time_wl,0.5,0.5,Rainbow}
	ModifyGraph zColorMax(R_soa)=(65280,0,0)
	ModifyGraph zColorMin(R_soa)=NaN
	ModifyGraph tick=2
	ModifyGraph mirror(bottom)=1
	ModifyGraph fSize(left)=12,fSize(bottom)=12
	ModifyGraph standoff=0
	ModifyGraph axOffset(left)=-2.33333,axOffset(bottom)=-0.722222,axOffset(right)=-3.28571
	ModifyGraph axThick(left)=1.5,axThick(bottom)=1.5
	ModifyGraph axRGB(right)=(65280,0,0)
	ModifyGraph tlblRGB(right)=(65280,0,0)
	ModifyGraph alblRGB(right)=(65280,0,0)
	Label left "\\Z13\\[0[SOA] (\\F'Symbol'm\\F'Arial'g m\\S-3\\M)"
	Label bottom "\\Z13\\[0Reaction Time (hours)"
	Label right "\\u#2\\Z13\\[0R\\BSOA"
	SetAxis left 0,*
	SetAxis right 1,7
	Legend/C/N=text0/J/X=4.08/Y=54.57 "\\s(Coa_experiment) Observed\r\\s(Coa_time_WL) with wall loss\r\\s(Coa_time_noWL) without wall loss\r\\s(R_soa) R\\BSOA\\M"
EndMacro

Window Graph23() : Graph
	PauseUpdate; Silent 1		// building window...
	Display /W=(151.5,44.75,708.75,434.75)
	AppendMatrixContour :FitInfo:Rsoa_OH_VOC_matrix
	ModifyContour Rsoa_OH_VOC_matrix manLevels=:FitInfo:wave0, ctabLines={*,*,Rainbow,1}
	ModifyContour Rsoa_OH_VOC_matrix logLines=1, labelSigDigits=2
	AppendImage/B=Bimage/L=Limage :FitInfo:Rsoa_OH_VOC_matrix vs {:FitInfo:OHconc,:FitInfo:VOCconc}
	ModifyImage Rsoa_OH_VOC_matrix ctab= {1,*,Rainbow,1}
	ModifyImage Rsoa_OH_VOC_matrix log= 1
	ModifyGraph lStyle=3
	ModifyGraph rgb=(30464,30464,30464)
	ModifyGraph log(Limage)=1,log(Bimage)=1
	ModifyGraph tick(bottom)=3
	ModifyGraph mirror(left)=0,mirror(bottom)=0,mirror(Limage)=2,mirror(Bimage)=2
	ModifyGraph noLabel(left)=2,noLabel(bottom)=2
	ModifyGraph axRGB(left)=(65535,65535,65535),axRGB(bottom)=(65535,65535,65535)
	ModifyGraph lblPos(left)=42,lblPos(bottom)=37,lblPos(Limage)=51,lblPos(Bimage)=48
	ModifyGraph lblLatPos(Limage)=2
	ModifyGraph freePos(Limage)=0
	ModifyGraph freePos(Bimage)=0
	Label Limage "[toluene] (ppb)"
	Label Bimage "[OH] (molecules cm\\S-3\\M)"
	Tag/C/N=acontour0/I=1/AO=2/F=0/Z=1/B=2/X=0.00/Y=0.00/L=0 'Rsoa_OH_VOC_matri=1.31281', 9, " \\{\"%.2g\",tagVal(3)} "
	ColorScale/C/N=text0/F=0/B=2/A=MC/X=-21.23/Y=-41.06 image=Rsoa_OH_VOC_matrix
	ColorScale/C/N=text0 vert=0, heightPct=4, width=200, tickLen=3, log=1
	ColorScale/C/N=text0 lblMargin=15, minor=1
	AppendText "R\\BSOA\\M"
EndMacro

Window graph17() : Graph
	PauseUpdate; Silent 1		// building window...
	Display /W=(35.25,42.5,556.5,431)
	AppendImage logCstar_Matrix vs {Ncarbons_Wave,Noxygens_Wave}
	ModifyImage logCstar_Matrix ctab= {0,*,Rainbow,1}
	ModifyImage logCstar_Matrix minRGB=(43520,43520,43520)
	ModifyGraph margin(top)=50
	ModifyGraph mirror=2
	ModifyGraph fSize=14
	ModifyGraph axOffset(left)=-0.7
	ModifyGraph axThick=1.5
	Label left "Number of Oxygen Atoms"
	Label bottom "Number of Carbon Atoms"
	ColorScale/C/N=text0/F=0/A=RT/X=11.79/Y=-16.02 image=logCstar_Matrix, vert=0
EndMacro

Window Graph24() : Graph
	PauseUpdate; Silent 1		// building window...
	String fldrSav0= GetDataFolder(1)
	SetDataFolder root:FitInfo:
	Display /W=(608.25,137,1146.75,544.25) alpha_wave[433,*][0] vs WLRgas_wave[433,*]
	SetDataFolder fldrSav0
	ModifyGraph mode=3
	ModifyGraph marker=19
	ModifyGraph msize=8
	ModifyGraph zColor(alpha_wave)={:FitResults:ChiSq_Sum[433,477][1],*,*,Rainbow}
	ModifyGraph logZColor=1
	ModifyGraph log=1
	ModifyGraph tick=2
	ModifyGraph mirror=1
	ModifyGraph fSize=14
	ModifyGraph standoff=0
	ModifyGraph axOffset(left)=-2.4,axOffset(bottom)=-0.333333
	ModifyGraph axThick=1.5
	Label left "Mass Accommodation Coefficient"
	Label bottom "\\f02k\\f00\\Bw\\M (s\\S-1\\M)"
	SetAxis left 0.0008,1
	ShowInfo
EndMacro

//**************************************************************************
Function GetTraces(start,stop)
	variable start
	variable stop
	
	variable i
	
	string df_myresults = "root:myresults"
	string df_data = "root"
	string newwave
	
	string df_fitinfo = "root:fitinfo"
	wave/T CompoundName = $(df_fitinfo+":CompoundName")
	wave/T NOx = $(df_fitinfo+":NOx")
	wave ChiSq = $(df_fitinfo+":ChiSq")
	wave MaskWave = $(df_fitInfo+":MaskWave")
	wave Np_wave = $(df_fitinfo+":Np_wave")
	wave alpha_wave = $(df_fitinfo+":alpha_wave")
	wave DynamicPartitioning_wave = $(df_FitInfo+":DynamicPartitioning_wave")
	wave DpSeed_wave = $(df_FitInfo+":DpSeed_wave")
	variable npnts
	
	for(i=start;i<=stop;i+=1)
		SOM_SetParamsFromFit(i)		
		SetFromCaltech(CompoundName[i],NOx[i])
//		DoWindow/F sompanel
		DoWindow/F FitGraph
		//
		NVAR KineticMassTransfer
		KineticMassTransfer = DynamicPartitioning_wave[i]
		NVAR alpha
		alpha = alpha_wave[i]
		NVAR Nparticles
		Nparticles = Np_wave[i]
		NVAR SeedDiameter
		SeedDiameter = DpSeed_wave[i]
		//
		wave Time_Experiment = $(df_Data+":Time_Experiment")
		som_v1(quiet=1,fittime=Time_Experiment)
		wave Coa_time = Coa_time
		
		newwave = df_myresults+":Coa_"+CompoundName[i]+"_"+NOx[i]+"_"+"wall"
		duplicate/o Coa_time $(newwave)
		
		NVAR gasWLR 
		gasWLR = 0
		som_v1(quiet=1,fittime=Time_Experiment)
		wave Coa_time = Coa_time
		newwave = df_myresults+":Coa_"+CompoundName[i]+"_"+NOx[i]+"_"+"nowall"
		duplicate/o Coa_time $(newwave)
	endfor
End

//****************************************************************************************
// Calculate log normal distributions based on input width, number median diameter, and total particles
// Based on diameter step size, rather than number of steps
Function MakeLogNormalDistn2(sigma,NMD,Np,[start,finish,delta])
	variable sigma	// log-normal spread
	variable NMD	// number median diameter in nm
	variable Np	// number concentration (p/cc is preferred)
	variable start	// first diameter
	variable finish	// last diameter
	variable delta	// diameter step size
	
	if(ParamIsDefault(start))
		start = 10
	endif	
	if(ParamIsDefault(finish))
		finish = 800
	endif
	if(ParamIsDefault(delta))
		delta = 5
	endif
	
	variable length
	length = (finish-start)/delta+1
	
	make/o/d/n=(length) Diameter, dNdlogDp, dSdlogDp, dVdlogDp
	variable m, n
	
	Diameter = start + delta*x
	dNdlogDp = (Np/(log(sigma)*sqrt(2*pi)))*exp(-(log(Diameter)-log(NMD))^2/(2*(log(sigma))^2))		
	dSdlogDp = pi*(Diameter^2)*dNdlogDp
	dVdlogDp = (pi/6)*(Diameter^3)*dNdlogDp	
end

//******************************************************************************
// Calculate log normal distributions based on input width, number median diameter, and total particles
// Written on 09/29/14 in support of polydisperse SOM simulations
Function MakeLogNormalDistn(sigma,NMD,Np,[nbins])
	variable sigma	// log-normal spread
	variable NMD	// number median diameter in nm
	variable Np	// number concentration (p/cc is preferred)
	variable nbins
	
	if(ParamIsDefault(nbins))
		nbins = 7
	endif
	
	variable start	// first diameter
	variable finish	// last diameter
	start = 10^(log(NMD)-3.7*log(sigma))
	finish = 10^(log(NMD)+3.7*log(sigma))
	variable delta
	delta = (log(finish)-log(start))/(nbins-1)
		
	make/o/d/n=(nbins) Diameter, dNdlogDp, dSdlogDp, dVdlogDp
	
	if(nbins<=1)
		Diameter = NMD
		delta = 1
		dNdlogDp = Np // (Np/(log(sigma)*sqrt(2*pi)))*exp(-(log(Diameter)-log(NMD))^2/(2*(log(sigma))^2))		
		dSdlogDp = pi*(Diameter^2)*dNdlogDp
		dVdlogDp = (pi/6)*(Diameter^3)*dNdlogDp
	else
		Diameter = 10^(log(start)+delta*x)
		dNdlogDp = (Np/(log(sigma)*sqrt(2*pi)))*exp(-(log(Diameter)-log(NMD))^2/(2*(log(sigma))^2))		
		dSdlogDp = pi*(Diameter^2)*dNdlogDp
		dVdlogDp = (pi/6)*(Diameter^3)*dNdlogDp // nm3/cm3
	endif
	
	return delta		
end

///////////////////////////////////////////////////////////////
// Graphs the calculated size distributions
Function GraphSizeDists(dXdlogDp,Dp,[stepsize])
	wave dxdlogdp
	wave dp
	variable stepsize
	
	if(paramisdefault(stepsize))
		stepsize=1
	endif
	
	variable nsteps = dimsize(dxdlogdp,0)
	variable i
	display
	for(i=0;i<nsteps;i+=stepsize)
		appendtograph dxdlogdp[i][] vs dp[i][]
	endfor
End

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Will take the SOM particle species and bin them according to the VBS (e.g. log10 bins)
Function BinSOM2VBS(particle_wave, logCstar_wave)
	wave particle_wave, logCstar_wave
	
	variable i
	variable Cstar_low = -6
	variable Cstar_high = 9
	variable npnts = (Cstar_high-Cstar_low)+1
	make/o/d/n=(npnts) VBS = 0
	setscale/P x Cstar_low, 1, VBS
	
	variable nC = dimsize(particle_wave,0)
	variable nO = dimsize(particle_wave,1)
	variable ntot = nC*nO
	
	for(i=0;i<npnts;i+=1)
		extract particle_wave, destWave, (Cstar_low+i-0.5) < logCstar_wave && (Cstar_low+i+0.5) > logCstar_wave
		if(numpnts(destwave)>0)
			wavestats/q destwave
			VBS[i] = V_sum
		endif
	endfor
	wavestats/q VBS
	VBS /= V_Sum
End

//***********************************************************************************
// Will reset the parameters (variables) in the SOM panel based on the results from fitting
// stored in various waves in the "BestFitParams" folder
Function SetParamsFromBestFit(num)
	variable num // index of row to set
	
	string df_BF = "root:BestFitParams"
	
	setdatafolder root:
	NVAR gasWLR // kwall
	NVAR delta_logCstar_perO
	NVAR FragSlope
	NVAR ProbOx1, ProbOx2, ProbOx3, ProbOx4
	NVAR Ncarbons
	
	setdatafolder $df_BF
	wave kwall
	wave mfrag
	wave deltaLVP
	wave pOx1, pOx2, pOx3, pOx4
	wave NcarbonsWave 
	NcarbonsWave = Ncarbons
	wave/T VOC
	wave/T NOx
	
	gasWLR = kwall[num]
	FragSlope = mfrag[num]
	delta_logCstar_perO = deltaLVP[num]
	probox1 = pox1[num]
	probox2 = pox2[num]
	probox3 = pox3[num]
	probox4 = pox4[num]
	Ncarbons = NcarbonsWave[num]
	Set_kOH(VOC[num]) // set krxn_parent
	
	print VOC[num] + "-" + NOx[num]
	setdatafolder root:
End

//***********************************************************************************************************
Function Set_kOH(compound)
	string compound // name of compound
	
	NVAR kOH = root:krxn_parent // in molecules/cm^3
	
	if(stringmatch(compound,"dodecane") || stringmatch(compound,"methylundecane") || stringmatch(compound,"cyclododecane") || stringmatch(compound,"hexylcyclohexane"))
		kOH = 1.34e-11
	elseif(stringmatch(compound,"benzene"))
		kOH = 1.22e-12
	elseif(stringmatch(compound,"toluene*") || stringmatch(compound,"toluene2013"))
		kOH = 5.63e-12
	elseif(stringmatch(compound,"mXylene"))
		kOH = 2.31e-11
	elseif(stringmatch(compound,"naphthalene"))
		kOH = 2.44e-11
	elseif(stringmatch(compound,"aPinene") || stringmatch(compound,"sesquiterpene"))
		kOH = 5.3e-11
	elseif(stringmatch(compound,"isoprene") || stringmatch(compound,"isoprene_alt")) // isoprene_alt added with v7.3.5
		kOH = 1e-10
	else
		kOH = 1e-11
		print "kOH set to default value"
	endif
	
	return kOH
End



//*********************************************************************************
Window sompanel_old() : Panel
	PauseUpdate; Silent 1		// building window...
	NewPanel /W=(161,84,508,882)
	ShowTools/A
	SetDrawLayer UserBack
	SetDrawEnv fillfgc= (65280,54528,48896),fillbgc= (65280,54528,48896)
	DrawRect 16,699,317,742
	DrawRect 152,101,317,176
	DrawText 161,119,"Heterogeneous Chem, etc."
	SetDrawEnv fillfgc= (48896,65280,65280)
	DrawRect 152,3,317,98
	DrawText 196,19,"Fragmentation"
	SetDrawEnv fillfgc= (65280,65280,48896)
	DrawRect 14,3,140,155
	DrawText 28,23,"General Operation"
	DrawRect 15,541,142,594
	SetDrawEnv fillfgc= (65280,48896,48896)
	DrawRect 14,427,142,537
	SetDrawEnv fsize= 11
	DrawText 19,445,"Oxygen Probability Array"
	DrawRect 14,159,141,320
	DrawText 18,176,"Kinetics and Oxidants"
	SetDrawEnv fillfgc= (65280,48896,48896)
	DrawRect 152,397,316,551
	DrawText 182,416,"Model Parameters"
	SetDrawEnv fillfgc= (57344,65280,48896)
	DrawRect 17,749,318,783
	SetDrawEnv fsize= 8
	DrawText 157,30,"0 = No x FragSlope"
	SetDrawEnv fsize= 8
	DrawText 157,40,"1 = (O:C)^FragSlope"
	SetDrawEnv fsize= 8
	DrawText 244,87,"1 = Equal"
	SetDrawEnv fsize= 8
	DrawText 244,77,"0 = Random"
	SetDrawEnv dash= 3
	DrawLine 169,65,298,65
	SetDrawEnv fsize= 8
	DrawText 106,187,"0 = OH"
	SetDrawEnv fsize= 8
	DrawText 106,195,"1 = O3"
	SetDrawEnv fsize= 8
	DrawText 244,97,"2 = Small Frags"
	SetDrawEnv fsize= 7
	DrawText 93,203,"Constant=0"
	SetDrawEnv fsize= 7
	DrawText 93,209,"Scaling=1"
	SetDrawEnv fsize= 7
	DrawText 93,215,"Experiment=2"
	SetDrawEnv fsize= 8
	DrawText 234,30,"2 = (No+1) x FragSlope"
	SetDrawEnv fsize= 8
	DrawText 234,40,"3 = FragSlope"
	SetDrawEnv fsize= 8
	DrawText 248,145,"0=1oxygen"
	SetDrawEnv fsize= 8
	DrawText 248,156,"1=multiple oxygens"
	SetDrawEnv fillfgc= (65280,65280,48896),fillbgc= (65280,65280,48896)
	DrawRect 152,182,317,335
	SetDrawEnv fillfgc= (57344,65280,48896)
	DrawRect 14,324,142,424
	DrawText 30,342,"Wall Loss Terms"
	DrawLine 21,258,135,258
	SetDrawEnv fsize= 7
	DrawText 111,250,"use 4"
	SetDrawEnv fsize= 8
	DrawText 234,51,"4 = new method"
	SetDrawEnv fillfgc= (57344,65280,48896)
	DrawRect 152,339,317,394
	DrawText 177,356,"Dynamic Partitioning"
	SetDrawEnv fillfgc= (65280,65280,48896),fillbgc= (65280,65280,48896)
	DrawRect 151,555,316,595
	DrawText 176,199,"Size Distribution Info"
	SetDrawEnv fillfgc= (65280,54528,48896),fillbgc= (65280,54528,48896)
	DrawRect 15,598,316,641
	SetDrawEnv fillfgc= (65280,48896,48896)
	DrawRect 15,645,316,695
	DrawText 109,662,"Single Precursor"
	DrawText 110,716,"Multiple Precursor"
	SetVariable setvar0,pos={156,199},size={152,16},title="# Particles (p/cm3)"
	SetVariable setvar0,value= Nparticles
	SetVariable setvar1,pos={180,299},size={128,16},title="Density (g/cm3)"
	SetVariable setvar1,value= Density_Base
	SetVariable setvar2,pos={54,25},size={82,16},limits={-inf,inf,0},value= Nsteps
	SetVariable setvar3,pos={32,43},size={106,16},title="Timestep (s)"
	SetVariable setvar3,value= timestep
	SetVariable setvar4,pos={22,62},size={116,16},title="Dilution (%/hr)"
	SetVariable setvar4,value= DilutionVar
	SetVariable setvar5,pos={20,80},size={118,16},title="VP adjustment"
	SetVariable setvar5,value= logCstar_adjustmentfactor
	CheckBox check0,pos={169,122},size={59,14},title="ON/OFF",variable= hetchem
	SetVariable setvar6,pos={157,156},size={110,16},title="Gamma OH",value= gammaOH
	CheckBox check1,pos={17,545},size={123,14},title="\\Z10Sequential Partitioning"
	CheckBox check1,fSize=11,variable= SPM
	SetVariable setvar7,pos={18,452},size={121,16},value= ProbOx1
	SetVariable setvar8,pos={18,472},size={121,16},value= ProbOx2
	SetVariable setvar9,pos={18,494},size={121,16},value= ProbOx3
	SetVariable setvar08,pos={18,515},size={121,16},value= ProbOx4
	SetVariable setvar10,pos={23,262},size={101,16},title="[OH]",value= OHconc
	SetVariable setvar11,pos={16,300},size={123,16},title="ScalingFactor"
	SetVariable setvar11,value= OH_scale
	SetVariable setvar12,pos={167,418},size={122,16},title="# carbon atoms"
	SetVariable setvar12,value= Ncarbons
	SetVariable setvar13,pos={189,437},size={100,16}
	SetVariable setvar13,limits={0,inf,0.01},value= FragSlope
	SetVariable setvar14,pos={199,456},size={90,16},title="dlVP"
	SetVariable setvar14,limits={0,3,0.01},value= delta_logCstar_perO
	SetVariable setvar15,pos={177,475},size={112,16},title="[VOC] (ppm)"
	SetVariable setvar15,value= ctot_ppm
	Button button0,pos={32,666},size={100,26},proc=Run_SOM,title="Run Model"
	Button button0,fSize=16,fStyle=1,fColor=(48896,52992,65280)
	Button button0,valueColor=(65280,0,0)
	SetVariable setvar16,pos={207,494},size={82,16},title="H-per-O",value= H_per_O
	SetVariable setvar17,pos={114,759},size={187,16},title="Folder"
	SetVariable setvar17,value= SaveFolderStr
	Button button1,pos={26,757},size={80,20},proc=ButtonProc_1,title="Save Results"
	SetVariable setvar19,pos={160,44},size={62,16},title="Pfrag",value= Pfrag_type
	SetVariable setvar20,pos={157,70},size={72,16},title="Method"
	SetVariable setvar20,value= small_fragments
	SetVariable setvar21,pos={182,513},size={107,16},title="H adjustment"
	SetVariable setvar21,value= Hadjustment
	SetVariable setvar22,pos={17,218},size={121,16},title="krxn(parent)"
	SetVariable setvar22,value= krxn_parent
	SetVariable setvar23,pos={19,178},size={82,16},title="OH or O3?",value= O3_yn
	SetVariable setvar25,pos={23,281},size={95,16},title="[O3] (ppb)",value= O3_conc
	SetVariable setvar26,pos={20,198},size={71,16},title="Method",fSize=12
	SetVariable setvar26,value= UseScalingForOxidant
	SetVariable setvar06,pos={18,99},size={120,16},title="RunTime (hrs)"
	SetVariable setvar06,value= MaxTime_Hours
	SetVariable setvar07,pos={161,138},size={82,16},title="Oxygens?"
	SetVariable setvar07,value= AddMultipleOxygens
	SetVariable setvar01,pos={155,279},size={160,16},title="Seed Conc (um3/cm3)"
	SetVariable setvar01,limits={-inf,inf,0},value= SeedVolConc
	Button button3,pos={150,660},size={149,31},proc=Run_SOM_Exp,title="Run Model w/ Exp"
	Button button3,fSize=16,fStyle=1,fColor=(48896,52992,65280)
	Button button3,valueColor=(65280,0,0)
	CheckBox check3,pos={190,559},size={89,14},title="Oligomerization"
	CheckBox check3,variable= OligomerizationIsOn
	SetVariable setvar24,pos={173,575},size={115,16},title="krxn,base\\M"
	SetVariable setvar24,value= krxn_base_olig
	SetVariable setvar27,pos={19,343},size={117,16},title="WLR(g) (1/s)"
	SetVariable setvar27,value= gasWLR
	CheckBox check4,pos={17,561},size={81,14},title="Eqm. Method"
	CheckBox check4,variable= EqmMethod
	SetVariable setvar28,pos={20,403},size={117,16},title="WLR(p) (scale)"
	SetVariable setvar28,value= ParticleWallLoss
	CheckBox check5,pos={190,318},size={97,14},title="Absorbing seed?"
	CheckBox check5,variable= AbsorbingSeed,side= 1
	SetVariable setvar29,pos={19,380},size={117,20},title="C\\Bwall\\M (mg/m3)"
	SetVariable setvar29,value= Cwall
	SetVariable setvar30,pos={19,361},size={117,16},title="WLR(g) scaling"
	SetVariable setvar30,value= kwall_gas_scaling
	SetVariable setvar31,pos={17,238},size={91,16},title="krxn Method",fSize=12
	SetVariable setvar31,value= krxn_method
	SetVariable setvar09,pos={18,117},size={120,16},title="Temperature (K)"
	SetVariable setvar09,value= Temp_G
	SetVariable setvar18,pos={20,136},size={118,16},title="Pressure (atm)"
	SetVariable setvar18,value= Pressure_G
	SetVariable setvar32,pos={169,532},size={120,16},title="DHvap (kJ/mol)"
	SetVariable setvar32,limits={0,inf,1},value= DHvap_G
	CheckBox check6,pos={17,576},size={80,14},title="Turn off SOA",variable= NoSOA
	SetVariable setvar02,pos={262,218},size={53,20},title="\\F'Symbol's\\F'Arial'\\Bg"
	SetVariable setvar02,value= SizeSpread
	SetVariable setvar03,pos={155,218},size={106,20},title="\\Z08D\\Bp,seed\\M\\Z08(nm)"
	SetVariable setvar03,value= SeedDiameter
	SetVariable setvar33,pos={168,357},size={127,16},title="Eqm=0; Dyn=1"
	SetVariable setvar33,value= KineticMassTransfer
	SetVariable setvar34,pos={186,375},size={109,16},title="Acc. coeff."
	SetVariable setvar34,value= alpha
	SetVariable setvar04,pos={158,260},size={148,16},title="Seed SA (um2/cm3)"
	SetVariable setvar04,limits={-inf,inf,0},value= SeedSurfaceArea
	CheckBox check7,pos={189,241},size={83,14},title="Polydisperse?"
	CheckBox check7,variable= polydisperse,side= 1
	SetVariable setvar35,pos={39,602},size={132,16},title="Bag Volume (m^3)"
	SetVariable setvar35,limits={0,inf,1},value= CSTR_Volume
	SetVariable setvar36,pos={53,621},size={119,16},title="Flow Rate (lpm)"
	SetVariable setvar36,limits={0,inf,1},value= CSTR_FlowRate
	SetVariable setvar05,pos={179,603},size={93,16},title="Tau (h)"
	SetVariable setvar05,limits={-inf,inf,0},value= CSTR_ResidenceTime_hr
	SetVariable setvar37,pos={53,621},size={119,16},title="Flow Rate (lpm)"
	SetVariable setvar37,limits={0,inf,1},value= CSTR_FlowRate
	SetVariable setvar38,pos={177,622},size={97,16},title="# Times"
	SetVariable setvar38,limits={0,inf,1},value= CSTR_NumTimes
	Button button2,pos={107,717},size={100,20},proc=Run_SOM_MP,title="Run Model"
	Button button2,fSize=16,fStyle=1,fColor=(48896,52992,65280)
	Button button2,valueColor=(65280,0,0)
EndMacro

Function RunSimulations(start,stop)
// currently assumes kwall = 0.0001, 0.00025 or 0
// Use this to rerun simulations based on previously determined best fit parameters
// using Caltech chamber conditions as input
	variable start, stop

	string df_bfp = "root:BestFitParams"
	string df_results = "root:SimulationResults"
	string df_no = df_results+":NoWallLoss"
	string df_low = df_results+":LowWallLoss"
	string df_high = df_results+":HighWallLoss"
	string df_new
	
	if(datafolderexists(df_results)==0)
		newdatafolder $df_results
		newdatafolder $df_no
		newdatafolder $df_low
		newdatafolder $df_high
		setdatafolder $df_no
			newdatafolder lowNOx
			newdatafolder highNOx
		setdatafolder $df_low
			newdatafolder lowNOx
			newdatafolder highNOx
		setdatafolder $df_high
			newdatafolder lowNOx
			newdatafolder highNOx
	endif
	
	wave/T VOCnames = $(df_bfp+":VOC")
	wave/T NOx = $(df_bfp+":NOx")
	wave kwall = $(df_bfp+":kwall")
	
	variable i
	string VOCstr, NOxStr
	variable kwall_current
	
	for(i=start;i<=stop;i+=1)
		VOCstr = VOCnames[i]
		NOxstr = NOx[i]
		kwall_current = kwall[i]
		setfromcaltech(VOCstr,NOxStr) // set the reaction conditions
		setparamsfrombestfit(i) // set the parameters
		// set wall loss to zero
		NVAR kwallVar = root:gasWLR
		kwallVar = 0
		SOM_V1(quiet=1) // run som
		
		wave Coa_time, Yield_time, deltaHC_time, O2C_time, Lifetimes, TimeW
		// need to get mask wave and interpolate to simulation timebase
		wave Coa_experiment_mask, Time_experiment
		Interpolate2/T=1/N=200/I=3/Y=Coa_Mask/X=TimeW Time_experiment, Coa_experiment_mask
		wave Coa_Mask
		Coa_mask = Coa_mask < 0.5 ? 0 : 1
		
		if(stringmatch(NOxstr,"low*")==1)
			NOxstr = "lowNOx"
		else
			NOxstr = "highNOx"
		endif
		if(kwall_current==0)
			df_new = df_no+":"+NOxstr+":"+VOCstr
		elseif(kwall_current==0.0001)
			df_new = df_low+":"+NOxstr+":"+VOCstr
		elseif(kwall_current==0.00025)
			df_new = df_high+":"+NOxstr+":"+VOCstr
		endif
		if(datafolderexists(df_new)==0)
			newdatafolder $df_new
		endif
		duplicate/o Coa_time $(df_new+":Coa_time")
		duplicate/o Yield_time $(df_new+":Yield_time")
		duplicate/o deltaHC_time $(df_new+":deltaHC_time")
		duplicate/o O2C_time $(df_new+":O2C_time")
		duplicate/o Lifetimes $(df_new+":Lifetimes")
		duplicate/o TimeW $(df_new+":TimeW")
		duplicate/o Coa_Mask $(df_new+":Coa_mask")
	endfor
	
End


Function RunAlkanes()

	variable NcStart = 6
	variable NcStop = 13
	variable Ctot_start = 2 // ppm
	variable Ctot_scale = 1.8
	variable Ctot_current = Ctot_start
	string myCoa
	
	NVAR Ncarbons = root:Ncarbons
	NVAR ctot_ppm = root:Ctot_ppm
	
	variable i
	for(i=Ncstart;i<=NcStop;i+=1)
		Ncarbons = i
		Ctot_ppm = Ctot_current
		
		SOM_V1(quiet=1)
		wave Coa_time, deltaHC_time
		duplicate/o Coa_time $("Coa_time_Nc"+num2istr(Ncarbons))
		duplicate/o deltaHC_time $("deltaHC_time_Nc"+num2istr(Ncarbons))
		Ctot_current /= Ctot_scale

	endfor	
	
End


//***********************************************************
Function RunAndSave()

	setdatafolder root:
	string dfSave = "dodecane"
	
	
	string dffpBase = "root:Results"
	string dffpSave = dffpBase +":"+dfSave
	
	print "Data saved in " + dffpSave
	if(datafolderexists(dffpBase) == 0)
		NewDataFolder $dffpBase
	endif
	if(datafolderexists(dffpSave)==0)
		NewDataFolder $dffpSave
	endif
	
	string Waves2Save = "RunInfo;RunValues;TimeW;Coa_time;O2C_time;H2C_time;HC_ppb_time;OH_wave;Dp_time;"
	Waves2Save += "Yield_time;deltaHC_time;deltaHCfrac_time;"
	Waves2Save += "OH_exposure;Lifetimes"
	
	variable i, j
	string OldWaveStr, NewWaveStr
	variable nWaves = itemsinlist(waves2save,";")
	Make/o/t/n=1 $(dffpSave+":RunTime")
	wave/T RunTime = $(dffpSave+":RunTime")

	NVAR Ctot_ppm
	variable nConc = 13
	variable C1 = 0.000144 // ppm
	variable Ctot_current = C1
	
	for(j=0;j<nConc;j+=1)
		Ctot_current = C1*2^j
		Ctot_ppm = Ctot_current
		dffpSave = dffpBase+":"+dfSave+":Run"+num2istr(j)
		if(datafolderexists(dffpSave)==0)
			NewDataFolder $dffpSave
		endif
		som_v1()
	
		RunTime[0] = "Run performed at " + time() + " on " + date()
		for(i=0;i<nWaves;i+=1)
			OldWaveStr = "root:"+stringfromlist(i,waves2save,";")
			NewWaveStr = dffpSave + ":"+ stringfromlist(i,waves2save,";")
			duplicate/O $OldWaveStr, $NewWaveStr
		endfor	
	endfor
end		

Function GraphIt()

	setdatafolder root:
	string dfbase = "root:results:dodecane:Run"
	string w1 = "Coa_time"
	string w2 = "Yield_time"
	variable npnts = 13
	variable i
	display
	for(i=0;i<npnts;i+=1)
		setdatafolder $(dfbase+num2istr(i))
		appendtograph $w2 vs $w1
	endfor
ENd
	

End

Function SOM_CSTR()

	NVAR CSTR_residencetime_hr
	NVAR CSTR_numtimes
	
	variable i_CSTR
	variable deltat_CSTR
	make/o/d/n=(CSTR_numtimes) CSTR_ProbDist
	make/o/d/n=(CSTR_numtimes) CSTR_Times, CSTR_CumProb
	
//	deltat_CSTR = CSTR_residencetime_hr/(CSTR_numtimes*2)
//	CSTR_Times = deltat_CSTR*(x+0.5)^1.8
//	deltat_CSTR = CSTR_residencetime_hr/(CSTR_numtimes)
//	CSTR_Times = deltat_CSTR*(x*1.5)+0.2
	deltat_CSTR = CSTR_residencetime_hr/(CSTR_numtimes*2)
	CSTR_times = 0.2+deltat_CSTR*(x^2)
	CSTR_CumProb = 1-exp(-CSTR_times/CSTR_residencetime_hr)
	for(i_CSTR=0;i_CSTR<CSTR_numtimes;i_CSTR+=1)
		
		if(i_CSTR==0)
			CSTR_ProbDist[i_CSTR] = CSTR_CumProb[0]
		elseif(i_CSTR==CSTR_numtimes-1)
			CSTR_ProbDist[i_CSTR] = 1-CSTR_CumProb[CSTR_numtimes-1]
		else
			CSTR_ProbDist[i_CSTR] = CSTR_CumProb[i_CSTR] - CSTR_CumProb[i_CSTR-1]
		endif
	endfor
	
//	CSTR_ProbDist[CSTR_numtimes-1] = CSTR_ProbDist[CSTR_numtimes-2]*0.9			
	
	wavestats/q CSTR_ProbDist
	CSTR_ProbDist/=V_sum
//	CSTR_ProbDist = (1/CSTR_residencetime_hr)*exp(-CSTR_times/CSTR_residencetime_hr)

// Run SOM multiple times, weight and average results
	NVAR MaxTime_hours
	NVAR nSizeBinsG
	make/o/d/n=(CSTR_numtimes) Coa_CSTR_wave=nan
	make/o/d/n=(CSTR_numtimes) Coa_CSTR_weighted_wave = nan
	make/o/d/n=(CSTR_numtimes) deltaHC_CSTR_wave=nan
	make/o/d/n=(CSTR_numtimes) Yield_CSTR_wave = nan
	make/o/d/n=(CSTR_numtimes,nSizeBinsG) Dp_CSTR_time = nan
	make/o/d/n=(CSTR_numtimes,nSizeBinsG) dNdlogDp_CSTR_time = nan
	
	for(i_CSTR=0;i_CSTR<CSTR_numtimes;i_CSTR+=1)
		MaxTime_hours = CSTR_Times[i_CSTR]
		SOM_v1(quiet=1)
		// get results and store them
		wave Coa_time // ug/m3
		wave deltaHC_time // ug/m3
		wave Diameter_Time // nm; 2D wave
		wave dNdlogDp_time // p/cc; 2D wave
		if(waveexists(Diameter_time)==0)
			abort "You are likely killing your diameter waves...need to update"
		endif
		variable npnts_thisrun = numpnts(Coa_time)-1
		Coa_CSTR_wave[i_CSTR] = Coa_time[npnts_thisrun]
		deltaHC_CSTR_wave[i_CSTR] = deltaHC_time[npnts_thisrun]
		Dp_CSTR_time[i_CSTR][] = Diameter_time[npnts_thisrun][q]
		dNdlogDp_CSTR_time[i_CSTR][] = dNdlogDp_time[npnts_thisrun][q]
		print "CSTR run " + num2str(i_CSTR+1) + " of " + num2str(CSTR_numtimes)
	endfor
	
	Coa_CSTR_weighted_wave = Coa_CSTR_wave*CSTR_Probdist
	wavestats/q Coa_CSTR_weighted_wave
	print "Weighted average SOA = " + num2str(V_sum) + " ug/m3"
	// deal with size distributions
	dNdlogDp_CSTR_time *= CSTR_probdist[p] // weight by probability distribution
	make/o/d/n=(nSizeBinsG) Dp_CSTR_single=nan, dNdlogDp_CSTR_single = nan
	variable nbins = 100
	make/o/d/n=(nbins) Dp_CSTR, dNdlogDp_CSTR=0
	make/o/d/n=(nbins,CSTR_numtimes) dNdlogDp_CSTR_2D = 0
	wavestats/q Dp_CSTR_time
	Dp_CSTR = V_min + x*(V_max-V_min)/(nbins-1)
//	for(i_CSTR=0;i_CSTR<CSTR_numtimes;i_CSTR+=1)
//		dp_cstr_single = dp_cstr_time[i_CSTR][p]
//		dndlogdp_cstr_single = dndlogdp_cstr_time[i_CSTR][p]
//		Interpolate2/T=2/I=3/Y=dNdlogDp_Interp/X=Dp_CSTR Dp_CSTR_single, dNdlogDp_CSTR_single
//		dNdlogDp_Interp = dNdlogDp_Interp < 0 ? 0 : dNdlogDp_Interp
//		dNdlogDp_CSTR += dNdlogDp_Interp
//		dNdlogDp_CSTR_2D[][i_CSTR] = dNdlogDp_interp[p]
//	endfor
	
	killwaves/z dNdlogDP_CSTR_single, dNdlogDp_CSTR_time, dNdlogDp_interp
	killwaves/z dP_CSTR_single, Dp_CSTR_time
	killwaves/z CSTR_CumProb, CSTR_ProbDist, CSTR_times, deltaHC_CSTR_wave
	killwaves/z dDpdt_off, dDpdt_on

End

function tester()
	wave dndlogdp_cstr_time
	wave cstr_probdist
	nvar nsizebinsg
	nvar cstr_numtimes
	variable i_cstr
	wave dp_cstr_time
	wave dndlogdp_cstr_time
	
//	dNdlogDp_CSTR_time *= CSTR_probdist[p] // weight by probability distribution
	make/o/d/n=(nSizeBinsG) Dp_CSTR_single=nan, dNdlogDp_CSTR_single = nan
	variable nbins = 100
	make/o/d/n=(nbins) Dp_CSTR, dNdlogDp_CSTR=0
	make/o/d/n=(nbins,CSTR_numtimes) dNdlogDp_CSTR_2D = 0
	wavestats/q Dp_CSTR_time
	Dp_CSTR = V_min + x*(V_max-V_min)/(nbins-1)
	for(i_CSTR=0;i_CSTR<CSTR_numtimes;i_CSTR+=1)
		dp_cstr_single = dp_cstr_time[i_CSTR][p]
		dndlogdp_cstr_single = dndlogdp_cstr_time[i_CSTR][p]
		Interpolate2/T=2/I=3/Y=dNdlogDp_Interp/X=Dp_CSTR Dp_CSTR_single, dNdlogDp_CSTR_single
		dNdlogDp_Interp = dNdlogDp_Interp < 0 ? 0 : dNdlogDp_Interp
		dNdlogDp_CSTR += dNdlogDp_Interp
		dNdlogDp_CSTR_2D[][i_CSTR] = dNdlogDp_interp[p]
	endfor
	
//	killwaves/z dNdlogDP_CSTR_single, dNdlogDp_CSTR_time, dNdlogDp_interp
//	killwaves/z dP_CSTR_single, Dp_CSTR_time
	
end

Function GraphSizeDist_CSTR()

	wave dndlogdp_cstr_2d
	wave dp_cstr
	variable npnts = dimsize(dndlogdp_cstr_2D,1)
	variable i
	display
	for(i=0;i<npnts;i+=1)
		appendtograph dndlogdp_cstr_2d[][i] vs dp_cstr
	endfor
end

function graphme()

	wave diameter_time
	variable npnts = dimsize(diameter_time,1)
	variable i
	display
	for(i=0;i<npnts;i+=1)
		appendtograph diameter_time[][i]
	endfor
end

Function TestingYields()

	NVAR maxtime_hours
	NVAR timestep
	variable npnts = (maxtime_hours*3600)/timestep + 1
	print npnts
	wave yield_time
	variable ntime = npnts//numpnts(yield_time)
	variable nseed = 20
	make/o/d/n=(ntime,nseed) yield_time_seed, coa_time_seed
	make/o/d/n=(nseed) myseed = 0.1+x
	make/o/d/n=(nseed) myseed_im
	variable i
	NVAR Nparticles
	Nparticles = 140
	for(i=0;i<nseed;i+=1)
		SOM_v1(quiet=1)
		wave yield_time, coa_time
		yield_time_seed[][i] = yield_time[p]
		coa_time_seed[][i] = coa_time[p]
		myseed[i] = Nparticles/1400 // ~1400 p/cc when Dp = 80 nm with sigma = 1.6 is about 1 ug/m3
		myseed_im[i] = Nparticles/1400
		Nparticles *= 1.3
		myseed_im[i+1] = Nparticles/1400
	endfor
End	

Function GraphYields()

	wave yield = yield_time_seed_abs_tol_no
	wave coa = coa_time_seed_abs_tol_no
	variable i
	variable npnts = dimsize(yield,1)
	display
	for(i=0;i<npnts;i+=1)
		appendtograph yield[][i] vs coa[][i]
	endfor
end
Window Graph18() : Graph
	PauseUpdate; Silent 1		// building window...
	Display /W=(811.5,281,1308,675.5) Yield_seed_no,Yield_seed_low,Yield_seed_high vs myseed
	AppendToGraph/R Yield_seed_high_DIF,Yield_seed_low_DIF,Yield_seed_no_DIF vs myseed
	ModifyGraph lSize=3
	ModifyGraph lStyle(Yield_seed_high_DIF)=3,lStyle(Yield_seed_low_DIF)=3,lStyle(Yield_seed_no_DIF)=3
	ModifyGraph rgb(Yield_seed_low)=(0,0,52224),rgb(Yield_seed_high)=(0,0,0),rgb(Yield_seed_high_DIF)=(0,0,0)
	ModifyGraph rgb(Yield_seed_low_DIF)=(0,0,52224)
	ModifyGraph tick=2
	ModifyGraph mirror(bottom)=1
	ModifyGraph fSize=14
	ModifyGraph standoff(left)=0,standoff(bottom)=0
	ModifyGraph axOffset(left)=-2,axOffset(bottom)=-0.7,axOffset(right)=-2
	ModifyGraph axThick=1.5
	Label left "Yield"
	Label bottom "[Initial Seed] (\\F'symbol'm\\F'arial'g m\\S-3\\M)"
	Label right "Yield Susceptibility"
	SetAxis left 0,1.2
	SetAxis bottom 0,15
	SetAxis right 0,*
EndMacro

//setparamsfrombestfit(11)
//setparamsfrombestfit(34)
//setparamsfrombestfit(58)
Function YieldsAndWL()

	NVAR maxtime_hours
	NVAR timestep
	variable npnts = (maxtime_hours*3600)/timestep + 1

	variable nrows = npnts
	variable ncols = 20
	make/o/d/n=(nrows,ncols) yield_time_seed_no, yield_time_seed_low, yield_time_seed_high
	make/o/d/n=(nrows,ncols) Coa_time_seed_no, Coa_time_seed_low, Coa_time_seed_high
	make/o/d/n=(ncols) yield_seed_no, yield_seed_low, yield_seed_high
	NVAR gasWLR
	make/o/d/n=3 myindices = {51,4,27} // {58,11,34} // {54,7,30} // 
	
	// No losses
	setparamsfrombestfit(myindices[0])
	gasWLR = 0
	testingyields()
	wave yield_time_seed
	wave coa_time_seed
	yield_time_seed_no = yield_time_seed[p][q]
	Coa_time_seed_no = Coa_time_seed[p][q]
	yield_seed_no = yield_time_seed_no[npnts][p]
	// low losses
	setparamsfrombestfit(myindices[1])
	gasWLR = 0
	testingyields()
	wave yield_time_seed
	wave coa_time_seed
	yield_time_seed_low = yield_time_seed[p][q]
	Coa_time_seed_low = Coa_time_seed[p][q]
	yield_seed_low = yield_time_seed_low[npnts][p]
	// high losses
	setparamsfrombestfit(myindices[2])
	gasWLR = 0
	testingyields()
	wave yield_time_seed
	wave coa_time_seed
	yield_time_seed_high = yield_time_seed[p][q]
	Coa_time_seed_high = Coa_time_seed[p][q]
	yield_seed_high = yield_time_seed_high[npnts][p]
	// other stuff
	wave myseed
	Differentiate Yield_seed_no/X=myseed/D=Yield_seed_no_DIF
	Differentiate Yield_seed_low/X=myseed/D=Yield_seed_low_DIF
	Differentiate Yield_seed_high/X=myseed/D=Yield_seed_high_DIF
	make/o/d/n=(npnts,ncols) Rwall_high, Rwall_low, Yield_Ratio
	setscale/P x, 0, (60/60/60), "hours", Yield_Ratio
	setscale/P x, 0, (60/60/60), "hours", Rwall_high
	setscale/P x, 0, (60/60/60), "hours", Rwall_low
	Rwall_high = coa_time_seed_high[p][q]/coa_time_seed_no[p][q]
	Rwall_low = coa_time_seed_low[p][q]/coa_time_seed_no[p][q]
	
	wave deltaHC_time
	make/o/d/n=(npnts) DDHC_time
	DDHC_time = deltaHC_time[p]-deltaHC_time[p-1]
	Make/o/d/n=(npnts,ncols) Yield_time_seed_no_instant, Yield_time_seed_low_instant, Yield_time_seed_high_instant
	Yield_time_seed_no_instant = (Coa_time_seed_no[p][q] - Coa_time_seed_no[p-1][q])/DDHC_time[p]
	Yield_time_seed_low_instant = (Coa_time_seed_low[p][q] - Coa_time_seed_low[p-1][q])/DDHC_time[p]
	Yield_time_seed_high_instant = (Coa_time_seed_high[p][q] - Coa_time_seed_high[p-1][q])/DDHC_time[p]
//	Differentiate/DIM=0  Coa_time_seed_no/X=deltaHC_time/D=Coa_time_seed_no_DIF
//	Differentiate/DIM=0  Coa_time_seed_low/X=deltaHC_time/D=Coa_time_seed_low_DIF
//	Differentiate/DIM=0  Coa_time_seed_high/X=deltaHC_time/D=Coa_time_seed_high_DIF
	yield_ratio = Yield_time_seed_high_instant/Yield_time_seed_no_instant
End	

Function GraphRwall()

	wave Rwall = Rwall_low
	variable i
	variable npnts = dimsize(Rwall,1)
	display
	for(i=0;i<npnts;i+=1)
		appendtograph Rwall[][i]
	endfor
end
