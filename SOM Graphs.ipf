#pragma rtGlobals=3		// Use modern global access method and strict wave access.

//*******************************************************************
// Functions to recreate/create key graphs/panels/tables
// v7.4.3 (see "SOM procs.ipf" for notes)

//****************************************************************
// Quick recreate of FitWindow graph
Function MakeGraph_FitWindow()

	DoWindow FitWindow
	if(V_flag==1)
		DoWindow/F FitWindow
		KillWindow FitWindow
	endif
	Execute "FitWindow()"

End

//**************************************************************************************************
// Graph of observed OA vs. time, with modeled value. Also includes calculated O:C and yield vs. time
Window FitWindow() : Graph
	PauseUpdate; Silent 1		// building window...
	Display /W=(162.75,84.5,612.75,523.25) Coa_experiment vs Time_experiment as "FitWindow"
	AppendToGraph Coa_time,fit_Coa_experiment
	AppendToGraph/L=L_yield Yield_time
	AppendToGraph/L=L_O2C O2C_time
	ModifyGraph mode(Coa_experiment)=3
	ModifyGraph marker(Coa_experiment)=19
	ModifyGraph lSize(Coa_time)=2,lSize(Yield_time)=3,lSize(O2C_time)=3
	ModifyGraph rgb(Coa_time)=(0,0,0),rgb(fit_Coa_experiment)=(0,0,65535),rgb(Yield_time)=(8704,8704,8704)
	ModifyGraph rgb(O2C_time)=(8704,8704,8704)
	ModifyGraph hideTrace(fit_Coa_experiment)=1
	ModifyGraph zColor(Coa_experiment)={Coa_experiment_mask,0.5,0.5,Grays}
	ModifyGraph zColorMax(Coa_experiment)=(52224,0,0)
	ModifyGraph zColorMin(Coa_experiment)=(52224,52224,52224)
	ModifyGraph tick(left)=2,tick(bottom)=2
	ModifyGraph zero(left)=1
	ModifyGraph mirror(left)=1,mirror(bottom)=1
	ModifyGraph fSize=14
	ModifyGraph standoff(left)=0,standoff(bottom)=0
	ModifyGraph axOffset(left)=-2,axOffset(bottom)=-0.7
	ModifyGraph axThick=1.5
	ModifyGraph lblPos(left)=55,lblPos(L_yield)=55,lblPos(L_O2C)=55
	ModifyGraph freePos(L_yield)={0,bottom}
	ModifyGraph freePos(L_O2C)={0,bottom}
	ModifyGraph axisEnab(left)={0,0.6}
	ModifyGraph axisEnab(L_yield)={0.64,0.8}
	ModifyGraph axisEnab(L_O2C)={0.84,1}
	Label left "SOA Concentration (\\F'symbol'm\\F'arial'g m\\S-3\\M)"
	Label bottom "Reaction Time (hours)"
	Label L_yield "Yield"
	ErrorBars Coa_experiment Y,wave=(Coa_experiment_err,Coa_experiment_err)
EndMacro

//****************************************************************
// Quick recreate of SOM Panel 
Function MakePanel_SOMpanel()

	DoWindow SOMpanel
	if(V_flag==1)
		DoWindow/F SOMpanel
		KillWindow SOMpanel
	endif
	Execute "SOMpanel()"

End

//***********************************************************************************************************************
// Recreate the SOM panel
Window SOMpanel() : Panel
	PauseUpdate; Silent 1		// building window...
	NewPanel /W=(261,63,809,743) as "SOMpanel"
	ShowTools/A
	SetDrawLayer UserBack
	DrawRect 171,577,356,661
	DrawText 187,596,"\\f01Heterogeneous Chem"
	SetDrawEnv fillfgc= (65280,65280,48896)
	DrawRect 14,3,164,170
	DrawText 28,23,"\\f01General Operation"
	DrawRect 15,472,165,522
	DrawRect 14,173,164,318
	DrawText 16,190,"\\f01Kinetics and Oxidants"
	DrawText 192,485,"Model Parameters"
	SetDrawEnv fsize= 8
	DrawText 132,201,"0 = OH"
	SetDrawEnv fsize= 8
	DrawText 132,209,"1 = O3"
	SetDrawEnv fsize= 7
	DrawText 116,220,"Constant=0"
	SetDrawEnv fsize= 7
	DrawText 116,226,"Scaling=1"
	SetDrawEnv fsize= 7
	DrawText 116,232,"Experiment=2"
	SetDrawEnv fsize= 8
	DrawText 275,626,"0=1oxygen"
	SetDrawEnv fsize= 8
	DrawText 275,637,"1=multiple oxygens"
	SetDrawEnv fillfgc= (65280,65280,48896),fillbgc= (65280,65280,48896)
	DrawRect 170,4,355,181
	SetDrawEnv fillfgc= (57344,65280,48896)
	DrawRect 14,322,165,422
	DrawText 30,339,"\\f01Wall Loss Terms"
	DrawLine 21,272,135,272
	SetDrawEnv fillfgc= (57344,65280,48896)
	DrawRect 171,184,356,239
	DrawText 195,201,"\\f01Dynamic Partitioning"
	SetDrawEnv fillfgc= (65280,65280,48896),fillbgc= (65280,65280,48896)
	DrawRect 14,427,165,467
	DrawText 194,21,"\\f01Size Distribution Info"
	SetDrawEnv fillfgc= (65280,54528,48896),fillbgc= (65280,54528,48896)
	DrawRect 15,528,165,628
	DrawLine 176,140,328,140
	SetDrawEnv fillfgc= (65280,48896,48896)
	DrawRect 172,244,355,573
	DrawText 201,262,"\\f01SOM Parameters"
	DrawText 197,277,"(single component)"
	DrawLine 178,427,329,427
	SetDrawEnv fsize= 8
	DrawText 255,523,"0 = Epstein relationship"
	DrawText 58,545,"\\f01 CFSTR"
	SetDrawEnv fillfgc= (65280,48896,48896)
	DrawRect 361,6,530,128
	DrawText 380,24,"\\f01Multicomponent Setup"
	SetDrawEnv fillfgc= (65280,54528,48896),fillbgc= (65280,54528,48896)
	DrawRect 361,183,530,339
	DrawText 375,202,"\\f01Single Component Setup"
	SetDrawEnv fillfgc= (65280,65280,48896),fillbgc= (65280,65280,48896)
	DrawRect 361,132,530,178
	DrawText 410,149,"\\f01Fit Condition"
	SetVariable setvar0,pos={174.00,41.00},size={163.00,19.00},title="# Particles (p/cm3)"
	SetVariable setvar0,value= Nparticles
	SetVariable setvar1,pos={206.00,82.00},size={142.00,19.00},title="Density (g/cm3)"
	SetVariable setvar1,value= Density_Base
	SetVariable setvar2,pos={54.00,61.00},size={82.00,19.00}
	SetVariable setvar2,limits={-inf,inf,0},value= Nsteps
	SetVariable setvar3,pos={27.00,24.00},size={125.00,19.00},title="Timestep (s)"
	SetVariable setvar3,value= timestep
	SetVariable setvar4,pos={22.00,79.00},size={132.00,19.00},title="Dilution (%/hr)"
	SetVariable setvar4,value= DilutionVar
	SetVariable setvar5,pos={20.00,97.00},size={134.00,19.00},title="VP adjustment"
	SetVariable setvar5,value= logCstar_adjustmentfactor
	CheckBox check0,pos={213.00,598.00},size={59.00,16.00},title="ON/OFF"
	CheckBox check0,variable= hetchem
	SetVariable setvar6,pos={176.00,637.00},size={110.00,19.00},title="Gamma OH"
	SetVariable setvar6,value= gammaOH
	CheckBox check1,pos={17.00,475.00},size={116.00,14.00},title="\\Z10Sequential Partitioning"
	CheckBox check1,fSize=11,variable= SPM
	SetVariable setvar7,pos={216.00,342.00},size={134.00,19.00},value= ProbOx1
	SetVariable setvar8,pos={216.00,363.00},size={134.00,19.00},value= ProbOx2
	SetVariable setvar9,pos={216.00,384.00},size={134.00,19.00},value= ProbOx3
	SetVariable setvar08,pos={216.00,406.00},size={134.00,19.00},value= ProbOx4
	SetVariable setvar10,pos={47.00,277.00},size={102.00,19.00},title="[OH]"
	SetVariable setvar10,value= OHconc
	SetVariable setvar11,pos={16.00,231.00},size={123.00,19.00},title="ScalingFactor"
	SetVariable setvar11,value= OH_scale
	SetVariable setvar12,pos={178.00,278.00},size={172.00,19.00},title="# carbon atoms"
	SetVariable setvar12,value= Ncarbons
	SetVariable setvar13,pos={231.00,321.00},size={119.00,21.00},title="m\\Bfrag\\M"
	SetVariable setvar13,limits={0,inf,0.01},value= FragSlope
	SetVariable setvar14,pos={233.00,297.00},size={117.00,23.00},title="\\F'symbol'D\\F'arial'lVP"
	SetVariable setvar14,limits={0,3,0.01},value= delta_logCstar_perO
	SetVariable setvar15,pos={184.00,433.00},size={149.00,19.00},title="[VOC] (ppm)"
	SetVariable setvar15,value= ctot_ppm
	Button button0,pos={393.00,260.00},size={104.00,25.00},proc=Run_SOM,title="Run SOM"
	Button button0,fSize=16,fStyle=1,fColor=(48896,52992,65280)
	Button button0,valueColor=(65280,0,0)
	SetVariable setvar16,pos={259.00,453.00},size={91.00,19.00},title="H-per-O"
	SetVariable setvar16,value= H_per_O
	SetVariable setvar21,pos={231.00,473.00},size={119.00,19.00},title="H adjustment"
	SetVariable setvar21,value= Hadjustment
	SetVariable setvar22,pos={18.00,251.00},size={141.00,19.00},title="krxn(parent)"
	SetVariable setvar22,value= krxn_parent
	SetVariable setvar23,pos={19.00,192.00},size={101.00,19.00},title="OH or O3?"
	SetVariable setvar23,value= O3_yn
	SetVariable setvar25,pos={19.00,296.00},size={129.00,19.00},title="[O3] (ppb)"
	SetVariable setvar25,value= O3_conc
	SetVariable setvar26,pos={20.00,212.00},size={85.00,19.00},title="Method"
	SetVariable setvar26,fSize=12,value= UseScalingForOxidant
	SetVariable setvar06,pos={18.00,42.00},size={134.00,19.00},title="RunTime (hrs)"
	SetVariable setvar06,value= MaxTime_Hours
	SetVariable setvar07,pos={180.00,616.00},size={91.00,19.00},title="Oxygens?"
	SetVariable setvar07,value= AddMultipleOxygens
	SetVariable setvar01,pos={172.00,162.00},size={174.00,19.00},title="Seed Conc (um3/cm3)"
	SetVariable setvar01,limits={-inf,inf,0},value= SeedVolConc
	Button button3,pos={366.00,288.00},size={156.00,25.00},proc=Run_SOM_Exp,title="Run SOM & Interp"
	Button button3,fSize=16,fStyle=1,fColor=(48896,52992,65280)
	Button button3,valueColor=(65280,0,0)
	CheckBox check3,pos={34.00,431.00},size={101.00,16.00},title="Oligomerization"
	CheckBox check3,variable= OligomerizationIsOn
	SetVariable setvar24,pos={21.00,447.00},size={115.00,19.00},title="krxn,base\\M"
	SetVariable setvar24,value= krxn_base_olig
	SetVariable setvar27,pos={19.00,340.00},size={138.00,19.00},title="WLR(g) (1/s)"
	SetVariable setvar27,value= gasWLR
	SetVariable setvar28,pos={20.00,400.00},size={136.00,19.00},title="WLR(p) (scale)"
	SetVariable setvar28,value= ParticleWallLoss
	CheckBox check5,pos={174.00,97.00},size={71.00,32.00},title="Absorbing\r seed?   "
	CheckBox check5,variable= AbsorbingSeed,side= 1
	SetVariable setvar29,pos={19.00,360.00},size={138.00,21.00},title="C\\Bwall\\M (mg/m3)"
	SetVariable setvar29,value= Cwall
	SetVariable setvar09,pos={18.00,115.00},size={135.00,19.00},title="Temperature (K)"
	SetVariable setvar09,value= Temp_G
	SetVariable setvar18,pos={35.00,133.00},size={118.00,19.00},title="Pressure (atm)"
	SetVariable setvar18,value= Pressure_G
	SetVariable setvar32,pos={217.00,493.00},size={133.00,19.00},title="DHvap (kJ/mol)"
	SetVariable setvar32,limits={0,inf,1},value= DHvap_G
	CheckBox check6,pos={17.00,491.00},size={84.00,16.00},title="Turn off SOA"
	CheckBox check6,variable= NoSOA
	SetVariable setvar02,pos={288.00,60.00},size={59.00,23.00},title="\\F'Symbol's\\F'Arial'\\Bg"
	SetVariable setvar02,value= SizeSpread
	SetVariable setvar03,pos={173.00,60.00},size={106.00,19.00},title="\\Z08D\\Bp,seed\\M\\Z08(nm)"
	SetVariable setvar03,value= SeedDiameter
	SetVariable setvar33,pos={187.00,199.00},size={127.00,19.00},title="Eqm=0; Dyn=1"
	SetVariable setvar33,value= KineticMassTransfer
	SetVariable setvar34,pos={205.00,218.00},size={109.00,19.00},title="Acc. coeff."
	SetVariable setvar34,value= alpha
	SetVariable setvar04,pos={184.00,144.00},size={162.00,19.00},title="Seed SA (um2/cm3)"
	SetVariable setvar04,limits={-inf,inf,0},value= SeedSurfaceArea
	CheckBox check7,pos={175.00,24.00},size={87.00,16.00},title="Polydisperse?"
	CheckBox check7,variable= polydisperse,side= 1
	SetVariable setvar35,pos={19.00,547.00},size={130.00,19.00},title="Bag Vol. (m3)"
	SetVariable setvar35,limits={0,inf,1},value= CSTR_Volume
	SetVariable setvar05,pos={45.00,606.00},size={102.00,19.00},title="Tau (h)"
	SetVariable setvar05,limits={-inf,inf,0},value= CSTR_ResidenceTime_hr
	SetVariable setvar37,pos={18.00,566.00},size={131.00,19.00},title="Flow Rate (lpm)"
	SetVariable setvar37,limits={0,inf,1},value= CSTR_FlowRate
	SetVariable setvar38,pos={43.00,586.00},size={106.00,19.00},title="# Times"
	SetVariable setvar38,limits={0,inf,1},value= CSTR_NumTimes
	Button button2,pos={389.00,96.00},size={115.00,25.00},proc=Run_SOM_MP,title="Run MC-SOM"
	Button button2,fSize=16,fStyle=1,fColor=(48896,52992,65280)
	Button button2,valueColor=(65280,0,0)
	SetVariable setvar19,pos={400.00,27.00},size={121.00,19.00},title="# Compounds"
	SetVariable setvar19,limits={1,10,1},value= nCompoundClasses
	PopupMenu popup0,pos={451.00,152.00},size={68.00,20.00},proc=PopMenuProc_NOx,title="NOx"
	PopupMenu popup0,mode=1,popvalue="low",value= #"\"low;high\""
	PopupMenu VWLcondition,pos={369.00,152.00},size={73.00,20.00},proc=PopMenuProc_VWL,title="VWL"
	PopupMenu VWLcondition,mode=3,popvalue="high",value= #"\"zero;low;high\""
	SetVariable setvar20,pos={394.00,46.00},size={127.00,19.00},title="Max alkane nC"
	SetVariable setvar20,limits={1,35,1},value= maxCforAlkanes
	Button button4,pos={413.00,68.00},size={65.00,25.00},proc=SetParams_MPSOM,title="Setup"
	Button button4,fSize=16,fStyle=1,fColor=(48896,52992,65280)
	Button button4,valueColor=(65280,0,0)
	PopupMenu Compound,pos={378.00,205.00},size={124.00,20.00},proc=PopMenuProc_SingleCompound,title="Compound"
	PopupMenu Compound,mode=5,popvalue="toluene",value= #"\"isoprene;aPinene;sesq;benzene;toluene;xylenes;naphthalene;n-alkane;b-alkane;c-alkane;isoprene_Alt\""
	Button button5,pos={393.00,231.00},size={104.00,25.00},proc=SetParams_SOM,title="Set Params"
	Button button5,fSize=16,fStyle=1,fColor=(48896,52992,65280)
	Button button5,valueColor=(65280,0,0)
	CheckBox check2,pos={17.00,506.00},size={90.00,16.00},title="First Gen Only"
	CheckBox check2,variable= FirstGenProductsOnly
	SetVariable setvar17,pos={271.00,22.00},size={73.00,19.00},title="# bins"
	SetVariable setvar17,limits={1,50,1},value= nSizeBinsG
	CheckBox check8,pos={388.00,319.00},size={119.00,16.00},title="Isoprene Is Special?"
	CheckBox check8,variable= IsopreneIsSpecial,side= 1
	CheckBox check9,pos={33.00,381.00},size={88.00,18.00},title="Krechmer C\\Bw\\M?"
	CheckBox check9,variable= Krechmer_Cw,side= 1
	SetVariable setvar30,pos={206.00,526.00},size={144.00,19.00},title="# O atoms in parent"
	SetVariable setvar30,value= Nox_precursor
	CheckBox check4,pos={20.00,152.00},size={117.00,16.00},title="Yield exclude POA?"
	CheckBox check4,variable= YieldCalc,side= 1
	SetVariable setvar31,pos={227.00,545.00},size={123.00,19.00},title="Cut off small nC"
	SetVariable setvar31,value= CutOffSmallStuff
	SetVariable setvar36,pos={188.00,122.00},size={142.00,19.00},title="NucleationTime (h)"
	SetVariable setvar36,value= NucleationTimeEmpirical
	SetVariable setvar39,pos={244.00,103.00},size={105.00,19.00},title=" Min Seed"
	SetVariable setvar39,value= minSeed
EndMacro

// LEGACY
Window FitGraph() : Graph
	PauseUpdate; Silent 1		// building window...
	Display /W=(125.25,71.75,669,497)/L=L_O2C O2C_experiment vs Time_experiment
	AppendToGraph Coa_experiment vs Time_experiment
	AppendToGraph Coa_time
	AppendToGraph/L=L_O2C O2C_time
	AppendToGraph/R Dp_time
	AppendToGraph/L=L_VOC HC_ppb_time
	AppendToGraph/R=R_VOC Lifetimes
	ModifyGraph mode(Coa_experiment)=3
	ModifyGraph marker(Coa_experiment)=8
	ModifyGraph lSize=2
	ModifyGraph lStyle(Dp_time)=3
	ModifyGraph rgb(O2C_experiment)=(65280,43520,0),rgb(Coa_experiment)=(65280,43520,0)
	ModifyGraph rgb(Coa_time)=(0,0,0),rgb(O2C_time)=(0,0,0),rgb(Dp_time)=(0,15872,65280)
	ModifyGraph rgb(HC_ppb_time)=(0,0,0),rgb(Lifetimes)=(65280,0,0)
	ModifyGraph opaque(Coa_experiment)=1
	ModifyGraph zColor(Coa_experiment)={Coa_experiment_mask,0.5,0.5,Grays}
	ModifyGraph zColorMax(Coa_experiment)=(65280,43520,0)
	ModifyGraph zColorMin(Coa_experiment)=(52224,52224,52224)
	ModifyGraph tick=2
	ModifyGraph zero(L_VOC)=1,zero(R_VOC)=1
	ModifyGraph mirror(L_O2C)=1,mirror(bottom)=1
	ModifyGraph lblMargin(right)=9
	ModifyGraph standoff(bottom)=0,standoff(left)=0
	ModifyGraph axOffset(left)=-1.88889
	ModifyGraph axRGB(right)=(0,15872,65280),axRGB(R_VOC)=(65280,0,0)
	ModifyGraph tlblRGB(right)=(0,15872,65280),tlblRGB(R_VOC)=(65280,0,0)
	ModifyGraph alblRGB(right)=(0,15872,65280),alblRGB(R_VOC)=(65280,0,0)
	ModifyGraph lblPos(L_O2C)=55,lblPos(left)=43,lblPos(right)=51,lblPos(L_VOC)=48,lblPos(R_VOC)=42
	ModifyGraph lblLatPos(left)=-2,lblLatPos(right)=-1,lblLatPos(L_VOC)=2,lblLatPos(R_VOC)=1
	ModifyGraph freePos(L_O2C)={0,bottom}
	ModifyGraph freePos(L_VOC)=0
	ModifyGraph freePos(R_VOC)=0
	ModifyGraph axisEnab(L_O2C)={0.51,0.74}
	ModifyGraph axisEnab(left)={0,0.49}
	ModifyGraph axisEnab(right)={0,0.49}
	ModifyGraph axisEnab(L_VOC)={0.76,1}
	ModifyGraph axisEnab(R_VOC)={0.76,1}
	Label L_O2C "O:C"
	Label bottom "Reaction Time (\\U)"
	Label left "C\\BOA\\M (\\U)"
	Label right "Particle Diameter w/seed (\\U)"
	Label L_VOC "VOC"
	Label R_VOC "Lifetimes"
	ErrorBars Coa_experiment Y,wave=(Coa_experiment_err,Coa_experiment_err)
	Cursor/P A Coa_experiment 141;Cursor/P B HC_ppb_time 462
	ShowInfo
	TextBox/C/N=text0/A=MC/X=3.02/Y=4.50 "isoprene_Coa_mfrag_0_#33"
EndMacro

//****************************************************************
// Quick recreate of FitSetup table
Function MakeTable_FitSetup()

	DoWindow FitSetup
	if(V_flag==1)
		DoWindow/F FitSetup
		KillWindow FitSetup
	endif
	Execute "FitSetup()"

End

// Table containing waves that hold inputs and results when performing fits
Window FitSetup() : Table
	PauseUpdate; Silent 1		// building window...
	String fldrSav0= GetDataFolder(1)
	SetDataFolder root:FItInfo:
	Edit/W=(24,137.75,912.75,649.25) CompoundName,NOx,HoldDLVPwave,MaskWave,FitTo,O2Cerr_scale as "FitSetup"
	AppendToTable FragType,Het_Wave,WLRgas_wave,InitialGuesses,FitResults,ChiSq,RunTime
	AppendToTable DpSeed_wave,Np_wave,DynamicPartitioning_wave,alpha_wave,Cwall_wave
	AppendToTable GammaOH_wave,WLRg_scaling_wave,WLRp_wave,kOH_wave,HoldStrWave,ErrorWave
	AppendToTable SPMwave,R_soa_avg,ChiSq_O2C
	ModifyTable format(Point)=1,width(CompoundName)=68,alignment(NOx)=1,width(NOx)=32
	ModifyTable width(HoldDLVPwave)=23,width(MaskWave)=32,alignment(FitTo)=1,width(FitTo)=40
	ModifyTable width(O2Cerr_scale)=32,alignment(FragType)=1,width(FragType)=32,width(Het_Wave)=20
	ModifyTable width(WLRgas_wave)=43,width(InitialGuesses)=35,sigDigits(FitResults)=3
	ModifyTable width(FitResults)=38,width(ChiSq)=42,format(RunTime)=8,width(RunTime)=72
	ModifyTable width(DpSeed_wave)=33,width(Np_wave)=38,width(DynamicPartitioning_wave)=27
	ModifyTable width(alpha_wave)=29,width(Cwall_wave)=26,width(GammaOH_wave)=22,width(WLRg_scaling_wave)=29
	ModifyTable width(WLRp_wave)=30,width(kOH_wave)=25,width(HoldStrWave)=30,width(ErrorWave)=29
	ModifyTable width(SPMwave)=28,width(R_soa_avg)=43,width(ChiSq_O2C)=32
	SetDataFolder fldrSav0
EndMacro
