## Notes on SOM v7 (through v 7.3.4)

Prof. Chris Cappa ([cdcappa@ucdavis.edu](mailto:cdcappa@ucdavis.edu))

University of California, Davis

Last update: 07/05/17

This &quot;manual&quot; describes how to run simulations and fit data using the SOM model of Cappa and Wilson (ACP, 2012). It discusses the following:

1. Setting up a SOM simulation using the panel
2. Running a single SOM simulation
  1. Single Component
  2. Multi-component
3. Accessing previously determined fit parameters for various VOCs
4. Accessing &quot;historical&quot; Caltech data
5. Fitting a single chamber observation using SOM
  1. Fitting assuming no vapor wall losses
  2. Determining an optimal kwall and alpha

This &quot;manual&quot; and the program may be difficult to use without a short tutorial. But you are welcome to try.

# Running SOM

The IGOR SOM is based on the model described in Cappa and Wilson (ACP, 2012) and Zhang et al. (PNAS, 2014). The various reactions/pathways are simulated using a Forward Euler method to solve the differential equations.

The IGOR SOM is primarily operated from the &quot;sompanel,&quot; where various inputs can be adjusted. Boxes with up/down arrows to the right are adjustable. The various controls are described below.

## SOM Procedures

Need to make sure that you have at least the following procedures loaded. The first set should be in the Igor file by default as they are embedded:

- Procedure Window (the default)
- OtherSOMStuff
- SOM\_Fitting
- aPinene2015

The next set may need to be loaded separately (and may not even be necessary), as they are external procedures.

- Global\_Utils\_CDC.ipf
- 2015GeneralMacros.ipf
- SizeDistributionProcessing\_v1.0.0.ipf

## SOM Functions

There are a variety of functions used when running SOM. The two key ones are:

- SOM\_V1() = the main SOM program for a single-component simulation
- SOM\_MC() = the main SOM program for a multi-component simulation

There are a variety of sub-functions that are accessed. These include:

- SOM\_CreateAtomMatrices(nC,nO,H\_per\_O,Hadjustment) = creates matrices with number of atoms for SOM species (e.g. oxygens, carbons, hydrogens, molecular weight)
- SOM\_RateCoefficients() = determine rate coefficients for SOM species (e.g. CxOy) based on particular rules
- SOM\_OligomerRateCoef() = determine rate coefficients for oligomerization…in development
- SOM\_Fragmentation() = determines the fragmentation matrices
- SOM\_ Fragmentation\_MP() = same as SOM\_ Fragmentation but for multiple components
- SOM\_GasPhaseWallLoss() = determines the species-specific wall desorption rates
- SOM\_ParticleWallLossRate() = determines particle wall loss rates
- SOM\_Heterogeneous() = calculates heterogeneous reaction rates
- EqmCalc() = calculates equilibrium gas-particle distribution based on input conditions
- EqmCalc\_MultipleVOC() = same as EqmCalc() but set up to work with multiple VOCs
- SOM\_Oligomerization() = runs oligomerization reactions…in development
- SOM\_KillWaves() = function to kill unnecessary waves at end of run
- SOM\_KillWaves\_MultiComponent() = kills additional unnecessary waves in compound-specific folders
- Makelognormaldistn() = function to make a log normal distribution based on the specified diameter, number, spread and number of bins
- MakeSOMwaves() = for multicomponent SOM only; makes folders to hold compound-specific information

And there are some functions to work with observational data and/or process data/graph model results:

- Graph\_MC\_yield() = graphs the compound-specific mass yield for multi-component simulations
- Graph\_MC\_GandP() = graphs the compound specific VOC decays and SOA formation as stacked plots for multi-component simulations
- Graph\_MCsoa() = graphs the compound specific SOA formation as stacked plots for multi-component simulations
- SOM\_SetupSingleComponent() = set SOM parameters based on stored fit
- SOM\_SetupMulticomponent() = set SOM parameters based on stored fits for multi-component simulations
- SetFromCaltech() = accesses historical data from the Caltech chamber
- SetFromCaltech2() accesses historical data from the Caltech chamber, specifically data from McVay et al. (2015)
- SOM\_BatchFit() = used to perform fit and store results for multiple experiments
- Allatonce\_FitCoa() = function used to fit a single experiment

# SOM Inputs

## General Operation

- **Timestep** : the timestep for the chemistry calculations.
  - Typically set to 60 seconds
  - When dynamic partitioning is used the mass transfer uses a timestep that is \&lt;= the specified timestep, which is important if one wants a stable simulation. This can be very short, meaning long simulations, when alpha is large (\&gt; 0.1).
- **Dilution** : rate of dilution, not really important for chamber experiments
- **VP adjustment** : just leave this at 1. This allows for adjustment of the parent VOC vapor pressure from the &quot;standard&quot; case (which is based on alkane backbones)
- **Runtime** : total run time. The number of steps (Nsteps) is calculated based on this and the timestep
- Nsteps: this is calculated, not an input. It is reported in the panel for reference.
- **Temperature** (K): self-explanatory
- **Pressure** (atm): self-explanatory

## Kinetics and Oxidants

- **OH or O3** : For simulations using OH or O3 as a reactant. The SOM is really set up for OH reactions (set to zero). O3 can be selected (set to 1), but this is not completely robust.
- **Method** : constant = 0; scaling = 1; experiment = 2
  - Constant = constant OH at the concentration set in the [OH] box
  - Scaling = will use a scaling factor to allow for an exponential decrease in [OH] from the initial [OH] with time. If Method = 1, then set the **scaling factor** input.
  - Experiment = will use a wave with time-varying OH concentrations for the calculations. The wave that is used is named &quot;OH\_exp&quot; and must have an associated time wave with name &quot;OH\_time&quot;. These are interpolated to the appropriate time steps for use in the simulation. This is useful for simulating chamber experiments and you know how [OH] varied with time.
- **krxn(parent)**: units = cm^3 molecules^‑1 s^‑1; The rate coefficient for the parent VOC
  - set this to zero to use a calculated krxn based on SAR relationships (the base case assumes an alkane).
  - Enter a value if you want your parent species have a different rate coefficient than the SAR relationships. Useful for simulation of aromatics, alkenes, etc. (i.e. those things that don&#39;t look like alkanes).
- **[OH]** = OH concentration in molecules/cm^3. Necessary if Method = 0 or 1. If Method = 1 this is the initial [OH].

## Wall Loss Terms

- **WLR(g)** = gas-phase wall loss rate in 1/s. Currently it is assumed that this is the same for all species.
- **WLR scaling** = scaling factor (legacy, set to zero)
- **Cwall** = effective absorbing mass concentration (mg/m^3) of the wall
- **Krechmer C**w**?** = decide whether to assume that the specified Cwall is the same for all species or whether the effective Cwall varies with C\* (as suggested by Krechmer et al. (ES&amp;T, 2016)).
  - Unchecked = constant Cwall
  - Checked = use Krechmer relationship
  - The default Krechmer relationship is scaled according to the **Cwall** variable
    - If Cwall = 1, default relationship
    - If Cwall = 10, leads to larger overall Cwall values
- WLR(p) = scaling factor for particle wall loss rate (usually set to zero). If set to 1, particles deposit irreversibly to the wall following the functional form in Loza et al. (ACP, 2010). This is a multiplicative factor onto the Loza et al. relationship (so setting to e.g. 2 means that particles are lost twice as fast as the observations suggest

## Size Distribution Info

- This sets the properties of the seed particles (if they exist). Only important when doing dynamic partitioning (and actually may not work with eqm. Partitioning).
- **Polydisperse?** : check box to decide whether to use a single size particle or assume a distribution of particle sizes. Note that things slow down considerably for polydisperse. But can be used to examine changes in shape of size distribution. If &quot;yes&quot; (checked) then a log-normal distribution is created that has **number of bins** having a geometric diameter **Dp,seed** with a spread of **sigmag** and a total number concentration of **number of particles** (in p/cm3).
- **No. of Bins:** the number of sizes to use if a polydisperse distribution is assumed. Use of 7 seems reasonable in most cases to get the general idea.
- **No. of particles:** (p/cm3) the number concentration of particles. The total seed surface area and volume depends on this. You can adjust the number of particles and the seed diameter to get the &quot;correct&quot; (i.e. observed) seed surface area when doing monodisperse particles.
- **Dp,seed**: in nm. Geometric median diameter for distribution or monodisperse size
- **sigmag**: the size distribution spread
- **Density:** (g/cm3) the assumed material density of the seed. Not really important unless you are assuming an absorbing seed.
- **Absorbing seed:** no check = seed does not participate in partitioning (i.e. does not mix). Check = seed participates in partitioning (i.e. mixes)
- Seed SA: the calculated seed surface area concentration
- Seed Conc: the calculated seed volume concentration

## Dynamic Partitioning

- **Eqm/Dyn** : select instantaneous equilibrium (0) or dynamic (1) partitioning. Dynamic runs a lot slower than eqm. How slow depends on what the accommodation coefficient is set to. Larger accommodation coefficients correspond to slower calculations. Very slow if you have the accommodation coefficient \&gt; 0.1. But in this case you can probably just assume equilibrium partitioning with minimal problems.
- **Acc. Coeff** : accommodation coefficient (only used when Eqm/Dyn = 1)

## SOM Parameters

- **# carbon atoms** : for a single component simulation, this is the number of carbon atoms of the parent VOC compound
- **DeltaLVP:** the assumed decrease in volatility to occur for each oxygen that is added. This is in log(C\*) units, where C\* is micrograms/m3.
  - typical range = 1-2.5
- **mfrag** the fragmentation parameter, that determines the fragmentation probability. Smaller numbers correspond to greater fragmentation. Set to 10 if you want (essentially) zero fragmentation. Set to 0.01 if you want a lot of fragmentation. Do not set to zero.
- Oxygen probability array: the probabilities of adding some number of oxygen atoms per reaction. These must sum to 1. However, you can enter whatever numbers you want and they will be normalized in the calculation.
  - **ProbOx1** = probability of adding 1 oxygen upon reaction
  - **ProbOx2** = …2 oxygens
  - **ProbOx3** = …3 oxygens
  - **ProbOx4** = …4 oxygens
- **[VOC] (ppm):** the concentration of the parent VOC, for single component simulations
- **H-per-O:** the number of hydrogens assumed to be lost per oxygen added. Default is 1. This affects the MW calculations, and thus the vapor pressures (slightly). It mostly impacts the calculated H:C, basically sets the value.
- **H-adjustment:**  ignore this and set to 0. can be used to adjust the parent MW to account for &quot;missing&quot; hydrogens relative to an alkane, for example as would happen for an aromatic compound. This can be set to the exact value, but can be left at zero too. If it is left at zero, the other fit parameters will simply be a little different than they would be if this is set to some other value. I suggest just leaving at 0.
- **DHvap (kJ/mol):** ignore this and set to 0. Potentially important for T-dependent simulations. Characterizes the variation in vapor pressures with temperature. 0 assumes the relationship given by Epstein et al.

## Heterogeneous Chem, etc.

- ON/OFF button: turn heterogeneous chemistry on or off
- Oxygens? = number of oxygens to add per heterogeneous reaction. 0 = 1 oxygen, 1 = multiple oxygens based on Oxygen Probability Array. The use of multiple oxygens is not yet robust so use with caution
- Gamma OH = OH effective uptake coefficient

## Other

- **Isoprene is Special:** there have been some arguments made (especially by Jose Jimenez) that we know &quot;enough&quot; about isoprene chemistry that we should be treating it more explicitly and with more knowns/constraints. The &quot;isoprene is special&quot; button overrides the default SOM multigeneration oxidation scheme to treat a few species uniquely. Details are given in Hodzic et al. (ACP, 2016). Note that the fit parameters for use with Isoprene Is Special need to be entered by hand.
- **Oligomerization** : This is a work in progress. Leave the check box unchecked
- **Sequential partitioning** : Check to use the &quot;sequential partitioning model&quot; from Cappa and Wilson (2011). It works just fine with eqm. Partitioning, but probably doesn&#39;t work correctly with the dynamic partitioning. A work in progress.
- **Turn Off SOA** : does not allow for SOA to form so that one can examine vapor wall loss only.
- **First Gen Only** : Turns off multi-generational oxidation, only allowing for production of &quot;first generation&quot; products. This is essentially the same as assuming things react with O3 (if there are not multiple double bonds and OH is not generated). Can be used to examine the influence that the multi-generational ageing has on SOA formation.
- **Eqm**. **Method** : (LEGACY; now set to default in code, not panel) When instantaneous equilibrium is used (or really dynamic partitioning too) one can choose to calculate the Raoult&#39;s Law vapor pressure depression based on mass concentrations or molar concentrations. The default is to have this checked and use the mass concentrations (after Donahue et al., ES&amp;T, 2006)

## Single Component Setup

- **Compound:** If you want to use best fit parameters (i.e. the SOM parameters) for some previously fit dataset you can select the compound that you want here. This is used in conjunction with the **Set Params** button and the **VWL** and **NOx** conditions under Fit Conditions.
- **Set Params:** click this to have the krxn(parent), number of carbon atoms, LVP, mfrag and ProbOx values set based on previous fits to observations. Is linked to **Compound** , **VWL** and **NOx**
- **Run SOM** : click this to run the model based on the parameters that have been set. This will run the function SOM\_v1().
- **Run SOM &amp; Interp** : click this to run the model based on the parameters that have been set, but then to interpolate the final results to the same timebase as a given experiment. There must exist a wave with name &quot;Time\_Exp&quot; that has the time (in hours) to which you want to interpolate the model results. This is useful and important when doing data fitting.

## Fit Condition

- You can access SOM parameters based on previously performed fits (c.f. Zhang et al., PNAS, 2014 and Jathar et al., GMD, 2015). The Fit Condition specifies what the assumption was regarding vapor wall losses during fitting and the NOx condition.
- **VWL:** if you want to setup to do calculations based on some previous fit, this specifies the assumption regarding vapor wall losses used in doing the original fit.
  - **Zero:** the fit was done assuming no wall loss (kwall = 0)
  - **Low:** the fit was done assuming &quot;medium&quot; wall loss (kwall = 1e-4)
  - **High:** the fit was done assuming &quot;high&quot; wall loss (kwall = 2.5e-4)
- **NOx:** fits to observations were performed for data collected under &quot;high&quot; and &quot;low&quot; NOx conditions. This specifies which fit to use.

# Multicomponent Setup

Some additional details on the multi-component SOM are available in Notebook0 in the Igor experiment

- It is possible to perform simulations assuming that you have more than one parent compound. This requires some user effort to get things set up correctly
- **# compounds:** the number of parent compounds to consider
- **Max alkane nC:** if alkanes are being considered, this is the maximum number of carbon atoms to be used. This is because if you want to have a multicomponent simulation with alkanes it is actually quite easy to do (although hard to track the evolution on a compound-by-compound basis). If you specify e.g. **max alkane nC** = 10, then you can have a multicomponent system with alkane precursors having carbon numbers ranging from 1 to 10.
- **Setup:** Will pop up a dialog and table to allow you to specify the details of your multicomponent system. More details are given below.
  - **oo**** Dialog:**
    - **Compounds:** select the names of the compounds that you want to consider.
    - **Concentrations:** enter the concentrations of the compounds (in ppb). Enter these as a semi-colon separated string (no spaces)
  - **oo**** Table:**
    - An example for a 4-component system with isoprene, benzene and n-alkanes and branched alkanes (with an imposed upper limit of 10 carbons, for this example) is available in the pdf
    - **ParameterNames** : just identifies the parameter names you are working with
    - **SOMParams** : a matrix that holds the SOM parameters (e.g. # carbons, mfrag) for each compound. The assumed compound name is indicated. These are prepopulated with values based on prior fits and the specified fit conditions (VWL and NOx). You can adjust as desired
    - **InitialConcMatrix** : (units = ppm) This is an _p_ x _q_ matrix, where _p_ is the maximum number of carbon atoms for any of the species being considered and _q_ is the number of species being considered.
      - Non-alkanes: For all compounds except alkanes, the initial concentration is entered in the row corresponding to the carbon number of the species of interest (counted from the bottom). These will automatically be entered for you based on the concentrations you entered in the pop-up dialog. But you can adjust them here without having to run the setup again.
      - Alkanes: For alkanes, the concentration in each row is the concentration of an alkane with a number of carbon atoms equal to |_p_-Max alkane nC+1|. So, if **Max alkane nC** = 10, the first row is for the species with 10 carbon atoms, the second for species with 9, etc. This is populated automatically with a distribution based on Zhao et al. (ES&amp;T, 2015) where the total concentration across all species is equal to the value that you entered in the pop-up dialog. If you only want one alkane compound (e.g. only the C10 alkane) then just set all other rows in that column to zero.
- **Run MC-SOM:** will run a simulation for a multicomponent system. Runs SOM\_MC().

# SOM Outputs

## Single Component SOM

### Key Waves

When the single component SOM is run a number of waves are generated. The key outputs are:

- **Coa\_time** = SOA mass concentration (ug/m^3)
- **O2C\_time** = O:C ratio
- **HC\_ppb\_time** = parent hydrocarbon concentration in ppb.
- **deltaHC\_time** = amount of parent VOC reacted, in ppb
- **TimeW** = time wave for run, in hours

If you chose to interpolate waves to an experimental time base the key output waves are:

- **Coa\_time\_interp** = SOA mass concentration (ug/m^3)
- **O2C\_time\_interp** = O:C ratio
- **deltaHC\_time\_interp** = amount of parent VOC reacted, in ppb

### Other Waves

A number of other waves are generated that store the results from a simulation. Many more are generated and killed at the end. The killing of waves is controlled by the SOM\_KillWaves() function. If you want something to show up at the end comment it out in the SOM\_KillWaves() function

- **Diameter\_time** = time-series of particle diameters for each size bin
- **Diameter\_final** = diameter by bin for last time point in the run
- **dNdlogDp\_time** = time series of dN/dlogDp (p/cc) for polydisperse simulations
- **dNdlogDp\_final** = final size distribution
- **GasMass\_time** = 3D matrix with time-series of each SOM species in the gas phase
- **ParticleMass\_time** = 3D matrix with time-series of each SOM species in the particle phase
- **GasMass\_Matrix** = 2D matrix of mass concentrations for each gas-phase species at last step
- **ParticleMass\_Matrix** = 2D matrix of mass conc. for each particle-phase species at last step

# Comparing to and Fitting Data

## Experimental Data

If one wants to fit data (or compare to data), the data must live in the root folder and have the following names

### Particle Phase

- Coa\_experiment = time-series of SOA concentrations (ug/m3)
- Coa\_experiment\_err = time-series of SOA concentration uncertainties (not necessarily used, but populate nonetheless
- Coa\_experiment\_mask = wave that indicates whether to include the associated Coa\_experiment data point as part of the fitting. 1 = include. 0 = exclude
- Time\_experiment = wave with time elapsed, in hours; must be same length as Coa\_experiment
- O2C\_experiment = time-series of SOA O:C ratios
- O2C\_experiment\_err = time-series of O:C uncertainties
- O2C\_experiment\_mask = same as Coa\_experiment\_mask

### Gas Phase

- VOC\_ppm = time-series of VOC concentrations; does not need to be same length as Coa\_experiment; must match Time\_VOC wave
- Time\_VOC = wave with elapsed time for VOC data; in hours
- OH\_exp = time-series of OH concentrations (probably derived from the VOC decay waves separately); units = molecules/cm3
- OH\_exp\_time = wave with elapsed time for OH concentrations; in hours
  - Generally try and make sure that the longest time here is longer than the run time in the model. For example, if the last point is as 9.3 h then you should only run for 9.29 h max.

## Accessing Historical Data

Historical data from Caltech are stored in the folder root:Caltech\_Data. This folder contains subfolders with the data this includes:

- Alkanes
  - C12 compounds
    - Dodecane
    - Cyclododecane
    - HexylCycloHexane
    - Methylundecane
- Aromatics
  - Benzene
  - mXylene
  - naphthalene
  - Toluene
- Biogenics
  - aPinene
  - isoprene

There exist data for both high/low NOx photooxidation experiments. You can access the data and have the e.g. Coa\_experiment waves in root populated using the function

**Setting &quot;Experiment&quot; values based on a given experiment**

1. You can grab the experimental conditions for a given experiment by using the function &quot;SetFromCaltech&quot; or &quot;SetFromCaltech2&quot;.
2. SetFromCaltech: This works for &quot;historical&quot; data that were collected under either high or low NOx conditions. Entries would look like e.g. SetFromCaltech(&quot;dodecane&quot;,&quot;low&quot;) for a low-NOx experiment using dodecane as the precursor VOC
3. SetFromCaltech2: This works for the data in McVay et al. (2016) performed for alpha-pinene under variable seed conditions and variable UV condtions. Entries would look like SetFromCaltech2(&quot;aPinene&quot;,&quot;low&quot;,&quot;low&quot;) where the first &quot;low&quot; refers to the UV condition and the second &quot;low&quot; to the seed concentration. This can be run from the command line or from a function.
  1. Data are stored in root:Caltech\_data:biogenics:OH:aPinene\_2016 and associated subfolders
    1. Data in subfolders are separated according to the UV and seed surface area during a given experiment
  2. Acceptable combinations of UV and seed conditions are:
    1. Low, low
    2. Low, high
    3. High, low
    4. High, med
    5. High, high
  3. The observed Coa time series will be copied to root:Coa\_experiment and root:Time\_experiment
    1. The &quot;mask&quot; wave that is stored in the data folder will be copied to Coa\_experiment\_mask
    2. The &quot;uncertainty&quot; wave that is stored in the data folder will be copied to Coa\_experiment\_err
  4. The particle number and seed diameter will be set to values that will give the correct seed surface area assuming monodisperse particles
  5. The OH method will be set to &quot;2&quot; (= experimental value) and the observationally constrained OH time series will be copied to root:OH\_exp and root:OH\_exp\_time
  6. The observed VOC decay time series will be copied to root:VOC\_ppm and root:Time\_VOC

## Fitting Data – One compound and one condition

- Igor has a built in fit routine that can be accessed either from dropdown menus or from the command line/functions. Data can be fit to the SOM using either.
- To fit to a function that is not built-in to Igor (such as SOM) requires that one write a function that defines a new fit routine. There is a procedure window titled &quot;SOM\_Fitting&quot;. The fit routines can be found here.
  - To fit just concentration data, use the &quot;allAtOnce\_FitCoa&quot; function.

### Fitting using the dropdown menu

- This is most easily done by putting the data that you want to fit on a graph and going from there. So, start by making a graph of Coa\_experiment vs time\_experiment (or just use the &quot;fitGraph&quot; that is already created).
- Make sure that the data you want to deal with is stored in the root folder in the waves Coa\_experiment and time\_experiment.
- Create a wave in root (if it doesn&#39;t already exist) called CoefWave (or whatever you want to call it) using the command make/o/d/n=(6) CoefWave = {0.5,1.75,0.25,0.25,0.25,0.25}. This wave stores the fit parameters, in this case {mfrag, DLVP, Pfunc1, Pfunc2, Pfunc3, Pfunc4}.
- Bring the graph to the front.
- Go to &quot;Analysis:Curve Fitting&quot;. A pop up dialog will come up.
  - On the &quot;Function and Data&quot; tab:
    - Select &quot;allatonce\_FitCoa&quot; from the &quot;functions&quot; dropdown.
    - Select &quot;From Target&quot;
      - This should set your &quot;Y data&quot; to &quot;Coa\_experiment&quot; and your &quot;X data&quot; to &quot;Time\_experiment&quot;. If this doesn&#39;t happen, you can set them by hand.
  - On the &quot;Data Options&quot; tab you can choose to restrict the range over which the data are fit. I usually have a separate wave (Coa\_experiment\_mask) that indicates whether data should be included or not. In this wave 1 = include and 0 = exclude. If you want to fit all of the data then this doesn&#39;t matter and you can leave it blank.
  - Click on the &quot;Coefficients&quot; tab. The first time you do this you might get a warning saying that it needs you to enter values for a user defined function. This is okay.
    - Under &quot;Coefficient Wave&quot; select the CoefWave that you generated above.
      - Enter your initial guesses here. When in doubt, try the following:
      - Mfrag = 1
      - DLVP = 1.8
      - Pfunc1-4 = 0.25
    - Under &quot;Constraints&quot; select &quot;From Coefficients List&quot;.
      - This puts limits on the values for e.g. DLVP, mfrag, Pfunc. This is important as none of these can be negative.
      - For mfrag (i.e. CoefWave\_0) use the range 0.01 to 10
      - For DLVP (CoefWave\_1) use the range 0.7 to 2.5
      - For Pfunc1-4 (CoefWave\_2-5) use the range 0.01 to 100. Don&#39;t worry that the largest value is \&gt; 1…this gets renormalized in the code.
    - I typically click on &quot;Graph Now&quot; to test whether things seem to be working, or to decide whether I want to change my initial guesses.
  - Once you&#39;re happy with your inputs, click on &quot;Do It&quot; on the bottom left. This will start things running. Depending on what type of calculations you are doing this can be &quot;fast&quot; or very slow. The dynamic partitioning is comparably slow, especially when alpha is large.
- Important: Make sure that the &quot;RunTime&quot; in the panel is longer than max value in the &quot;Time\_Experiment&quot; wave. Otherwise you will run into errors.

## Fitting Data – Multiple compounds, one condition at a time

THIS ONLY WORKS FOR FITTING THE HISTORICAL CALTECH DATA. But the process can be adapted.

This works if you have stored your data (concentrations vs. time, OH vs. time, etc.) in folders. There is a table called &quot;FitSetup&quot; where you will specify the compound to fit

- Create a new line (or over-write an old line) in the &quot;FitSetup&quot; table. Note that there are a lot of different waves in this table, but they must all be filled in with the correct information.
- Waves to populate before fitting:
  - &quot;CompoundName&quot; is the compound name associated with the appropriate folder where the data are stored.
  - &quot;NOx&quot; is the NOx condition (low vs. high). Note that &quot;NOx&quot; might not simply be low or high as I had to accommodate situations where there are multiple data sets for a given compound for a given NOx condition. Each label (e.g. low vs. low2) has an associated data folder, that can be checked by looking at the GetDataFolderName\_Caltech(Compound,NOx) function.
  - HoldWave = ignore this and set to 0
  - MaskWave: 0 = no mask applied during fitting; 1 = mask (Coa\_experiment\_mask) applied during fitting
  - FitTo: Coa = fit to Coa\_experiment only (use this); &quot;Coa and O2C&quot; = fit to both Coa and O2C (don&#39;t use this)
  - O2Cerr\_scale: set to 1 and ignore
  - FragType: decide whether to use mfrag (Pfrag = (O:C)^mfrag) or cfrag (Pfrag = No\*cfrag). Use mfrag
  - Het\_Wave: 0 = no heterogeneous chemistry; 1 = heterogeneous chemistry
  - WLRgas\_wave: vapor wall loss rate in 1/s
  - InitialGuesses: initial guesses for mfrag, DLVP, Pfunc1, Pfunc2, Pfunc3, Pfunc4
  - NpWave = seed particle number concentration (p/cc)
  - GammaOH\_wave: OH uptake coefficient for use when Het\_Wave = 1
  - Cwall\_Wave: effective wall concentration (mg/m^3)
  - WLRg\_scaling: some legacy scaling factor. Set to 0.
  - WLRp\_wave: 0 = no particle loss; 1 = particle loss according to equation given in SOM mechanism (derived based on comparison with real Caltech chamber observations, i.e. Loza et al.)
  - SPMwave: to use the sequential partitioning model, set this to 1. Default is zero. Note that SPM has not been confirmed to work correctly with dynamic partitioning.
  - kOH\_wave: what method to use for specifiying the rate coefficient matrix. Use 4.
  - ErrorWave: legacy. Ignore and set to 0.
  - HoldStr: set to none (unless you really want to hold some of the coefficients constant)
  - DpSeed\_Wave: seed particle diameter
  - Alpha\_wave: mass accommodation coefficient. Only used if DynamicPartitioning\_Wave = 1.
  - DynamicPartitioning\_Wave: 0 = instantaneous eqm; 1 = dynamic partitioning.
- Waves that get populated upon fitting
  - FitResults: results for mfrag, DLVP, Pfunc1, Pfunc2, Pfunc3, Pfunc4
  - ChiSq: chisq associated with best fit…actually, this is always calculated incorrectly for a reason that I still haven&#39;t figured out. You can get the correct value by running &quot;RecalculateChi2\_Batch(start,stop)&quot;
  - RunTime: when did you do the fit
  - ChiSq\_O2C: chisq associated with O2C fit. Ignore this.
- To do the fit:
  - Run the function SOM\_BatchFit(startpos,stoppos,dataset) where startpos and stoppos are the indices associated with the data that you want to fit and dataset = &quot;Caltech&quot;. For example, based on what is currently set up in the FitSetup table, if you run SOM\_BatchFitCaltechData(0,4) this will perform fits for dodecane, methylundecane, cyclododecane and hexylcyclohexane low NOx experiments. This function will automatically grab data from the appropriate folder and copy it to Coa\_experiment and Time\_experiment in the root folder. (It also grabs O2C\_experiment and a few other things.) Once you set this going it should start to fit dodecane. Once it completes that fit it will start on the next compound, etc., etc.
  - When you run this it also grabs the starting concentration from the &quot;RunParameters&quot; wave that exists in each data folder. It also grabs some additional information from this wave, including how to treat the [OH] (e.g. constant, time-dependent according to some function or based on some input wave), the number of carbon atoms in that molecule, the rate coefficient associated with degradation of the parent species (important for aromatics and alkenes).
  - Usually the fit algorithm does pretty well in finding a solution. But it can get lost or stuck, especially if you give a poor initial guess. So watch out for errors, which stop the fitting routine. To kill a fit use &quot;ctrl-shift-break&quot;. The normal &quot;abort&quot; button on the bottom left of the Igor screen does not seem to work all the time.
- After you&#39;ve performed the fits, if you want to reload data for a given compound _and_ you want to show the SOM simulation results for the compound based on a given fit, then you can run RecalculateChi2\_Batch(start,stop), where start and stop serve the same purpose as startpos and stoppos above.

## Fitting Data – Determine optimal kwall and alpha (accommodation coefficient)

The goal is to determine the kwall and alpha value pair that provides for the best fit of SOM to the observations. This can be done for a single experiment (i.e. a given UV/seed pair) or one can aim to find the (kwall,alpha) that gives the best global fit for a given UV condition. Descriptions of the functions that can be used to perform batch fitting are given below. These functions are stored in the &quot;aPinene2016&quot; procedure window.
1. **BatchFit\_SOM\_kwall\_alpha:** A _single_ experiment (UV/seed combination) can be fit using the function BatchFit\_SOM\_kwall\_alpha(compound,[NOx,UV,seed,reset]). The observations will be fit using SOM\_v1, assuming dynamic mass transfer.
  1. The inputs are
    1. Compound always = &quot;aPinene&quot;
    2. UV can be either &quot;low&quot; or &quot;high&quot;
    3. Seed can be &quot;low&quot;, &quot;med&quot; or &quot;high&quot;
    4. Reset is by default 0, meaning no reset. If Reset = 1, the results matrices are re-created and set to default values.
  2. The kwall and alpha values to be considered are hard wired in the code and can be changed. The min and max values considered are specified, along with the total number of points between (and including) the min/max values
  3. The maximum number of iterations per (kwall,alpha) pair is hard wired at 20
  4. Saved information includes
    1. ChiSq\_matrix = the V\_chisq results from fitting
    2. FitResults\_matrix = the SOM fit parameters determined
    3. FitError\_matrix = indicates whether an error occurred
    4. FitQuitReason\_matrix = indicates the reason that the fit ended, which can be used to identify which (kwall,alpha) pairs led to the max iterations being reached (or some other error)
  5. Note: when alpha \&gt; 0.1, the fitting proceeds very slowly because SOM runs very slowly (due to the need for a small time-step to deal with fast mass transfer). However, alpha = 0.1 typically gives very similar results to alpha = 1 because gas-particle mass transfer is sufficiently fast when alpha = 0.1 that the gas and particles (approximately) reach equilibrium on short (\&lt; 1 min) timescales. Thus, the code is currently hardwired to use alpha\_max = 0.1, but then to add one additional run using alpha = 1 but where equilibrium partitioning is assumed.
2. **BatchFit\_SOM\_aP2016()****:** All of the experiments can be individually fit using the function BatchFit\_SOM\_aPinene2016(). This function loops over the different experimental conditions (i.e. UV/seed combinations) and runs BatchFit\_SOM\_kwall\_alpha() for each condition. The resulting results matrices have &quot;\_XY&quot; appended to the end, where X = the UV condition (either L or H) and Y = the seed condition (L, M or H). The user can add an optional additional suffix to the results matrix names by changing the &quot;special\_suffix&quot; string. If you have this as &quot;&quot; then no suffix is appended. The use of special\_suffix can be helpful if you want to run a batchfit for alternative conditions (e.g. using the Jimenez/Krechmer Cw parameterization).
3. **BatchFit\_SOM\_kwall\_alpha\_Refit:** If you have performed a batch fit for a single experiment and find that one (or more) of the (kwall,alpha) pairs gave a particularly bad result you can rerun the fitting and update the results matrices for just these (kwall,alpha) pairs. This function relies on the &quot;Refit\_matrix&quot; wave. Re-fitting will be performed for only those elements that are set to 1. One can use this function to restart a run that was aborted early or to redo some select subset of (kwall,alpha) pairs. This can only be run after BatchFit\_SOM\_kwall\_alpha has been run at least once in &quot;reset&quot; mode, as it relies on waves having already been created. You can use the &quot;special\_suffix&quot; string.
4. **BatchFit\_SOM\_aP2016\_compare**** :** Once you have run BatchFit\_SOM\_aP2016() to determine the best fit for every (kwall,alpha) pair for a given experiment, you can run this function to see how well these best fits reproduce the complementary experiment for a given UV condition. For example, the SOM parameter sets determined from fitting of the &quot;low UV, low SA&quot; experiment are used to simulate the complementary &quot;low UV, high SA&quot; experiment. Chi square values are calculated for each of the experiments and for the combination. You can use the &quot;special\_suffix&quot; string.
  1. For low UV, two matrices are created that have 3 layers. Layer 0 corresponds to the ChiSq values for the experiment to which the SOM fit was originally performed (low or high SA). Layer 1 corresponds to the ChiSq values for the complementary experiment (high or low SA). And Layer 2 corresponds to the sum of the ChiSq&#39;s from the first two layers.
  2. For high UV, it is similar to low UV except there are now 4 layers to account for there being three seed conditions (low, medium, high), with the last layer the sum over the first three.
  3. The resulting chisq matrices are named according to the data set that was originally fit. For example, ChiSq\_matrix\_LL\_compare, ChiSq\_matrix\_LH\_compare, etc., and where the &quot;\_compare&quot; indicates that you are comparing between experiments.

To perform a batch &quot;optimization&quot; for every experimental condition, the following should be done.

1. Run BatchFit\_SOM\_aP2016(). This will run through every condition (UV,Seed pair) and determine a matrix of chi-square values associated with each (kwall,alpha) pair. The best-fit SOM parameters for each (kwall,alpha) pair are also stored. The results waves are appended with e.g. &quot;\_LL&quot; to indicate the UV and seed condition.
  1. This runs the BatchFit\_SOM\_kwall\_alpha() function automatically for each UV,seed pair
  2. This may take a very long time to complete
2. Run BatchFit\_SOM\_aP2016\_compare() for each condition considered. This function is set up to operate only on a single experimental condition, so will need to be run 5 times (or you can write a function that calls this 5 times and loops through the different experimental conditions).
3. Make image plots showing the results (i.e. the chi-square matrices) for each condition.

When starting out, it may be useful to perform the optimization using a limited set of values for (kwall,alpha). You can examine the general behavior using a coarse resolution, and then increase the resolution (i.e. number of points) later. The number of points is set using the nsteps\_kwall and nsteps\_alpha variables in BatchFit\_SOM\_kwall\_alpha(). Starting values of 5 for both seems reasonable, but higher resolution (or a more limited range) should ultimately be considered. Using nsteps\_kwall = 13 and nsteps\_alpha = 14 seem reasonable for &quot;high resolution&quot; optimization.

The optimization process is slow. During an optimization some fits proceed fast while others progress only slowly, due either to a greater number of iterations being required or a smaller timestep being required.  You can run multiple optimizations (e.g. for the different cases above) simultaneously by running multiple instances of Igor. Igor does not operate like e.g. Excel in that you cannot by default have multiple instances open. However, this can be circumvented by holding ctrl while opening a file. This opens up a second instance by running the executable again. You can actually have the same file open twice, although this is not recommended. Instead, it is recommended that you duplicate your file (\*.pxp) of interest and then open the original and the copy. Most computers seem to handle two instances reasonably well, and since most computers now have dual cores they can multi-task so running two optimizations does not slow down one. However, if you start to run too many instances of Igor things can get unstable. Igor does not recover like Microsoft products, so if something crashes and you have not saved recently then your work is lost. (You could add the SaveExperiment command to run at various times during the optimization so that the file is automatically saved. Just don&#39;t have this run too often.) In any case, you could, for example, perform the &quot;default&quot; optimization simultaneous with the &quot;Krechmer alternative&quot; optimization.

When doing a &quot;batch fit&quot; the &quot;abort&quot; button does not work. Normally, the abort button (that pops up in the lower-left of Igor when a procedure is running) will stop the function from executing. But this doesn&#39;t work when doing fitting (for some reason I still don&#39;t understand). However, you can force a quit by using ctrl+break or ctrl+pause. (On my laptop it is ctrl+break. On my keyboard it is ctrl+pause, I think because there is no break button and pause = break on the keyboard.) This can be useful if you want to force stop an optimization to make some change and then restart.
