import os

dirList = os.listdir("./")
dirs = [dir for dir in dirList if dir.startswith("00")]


for dir in dirs:
    file = open(f"{dir}/params.dat", 'r')
    params = []

    for line in file:
        params.append(line.split()[1])

    RHIC_Diject = f"""# -*- ThePEG-repository -*-

read snippets/PPCollider.in

##################################################
# Technical parameters for this run
##################################################
cd /Herwig/Generators
#ls
#cd EventGenerator
#ls

#create ThePEG::RandomEngineGlue /Herwig/RandomGlue 
#set LHCGenerator:RandomNumberGenerator /Herwig/RandomGlue 
#set LHCGenerator:NumberOfEvents 100
#set LHCGenerator:PrintEvent 0
#set LHCGenerator:MaxErrors 10
#set LHCGenerator:DebugLevel 0
#set LHCGenerator:DumpPeriod -1
#set LHCGenerator:DebugEvent 0

#sqrts 
set EventGenerator:EventHandler:LuminosityFunction:Energy 200.0

# Set the PDF 
#cd /Herwig/Particles 
#set p+:PDF /Herwig/Partons/cmsPDFSet 
#set pbar-:PDF /Herwig/Partons/cmsPDFSet 
#set K0:Width 1e300*GeV 
#set Kbar0:Width 1e300*GeV

# Declare stable particles 
cd /Herwig/Particles 
set mu-:Stable Stable 
set mu+:Stable Stable 
set Sigma-:Stable Stable 
set Sigmabar+:Stable Stable 
set Lambda0:Stable Stable 
set Lambdabar0:Stable Stable 
set Sigma+:Stable Stable 
set Sigmabar-:Stable Stable 
set Xi-:Stable Stable 
set Xibar+:Stable Stable 
set Xi0:Stable Stable 
set Xibar0:Stable Stable 
set Omega-:Stable Stable 
set Omegabar+:Stable Stable 
set pi+:Stable Stable 
set pi-:Stable Stable 
set K+:Stable Stable 
set K-:Stable Stable 
set K_S0:Stable Stable 
set K_L0:Stable Stable 

# order Alpha s
cd /Herwig
create Herwig::O2AlphaS O2AlphaS
set Model:QCD/RunningAlphaS O2AlphaS

#cd /Herwig/Partons
#create ThePEG::LHAPDF myPDFset ThePEGLHAPDF.so
#set myPDFset:PDFName cteq6ll.LHpdf
#set myPDFset:RemnantHandler HadronRemnants
#set /Herwig/Particles/p+:PDF myPDFset
#set /Herwig/Particles/pbar-:PDF myPDFset

#cd /Herwig/Partons 
#create ThePEG::LHAPDF CTEQ6L1 ThePEGLHAPDF.so 
#set CTEQ6L1:PDFName cteq6ll.LHpdf 
#set CTEQ6L1:RemnantHandler HadronRemnants 
#cp CTEQ6L1 cmsPDFSet 

# EE4C Tune 
set /Herwig/UnderlyingEvent/MPIHandler:EnergyExtrapolation Power
set /Herwig/UnderlyingEvent/MPIHandler:ReferenceScale 200.*GeV
set /Herwig/UnderlyingEvent/MPIHandler:Power {params[0]}
#set /Herwig/UnderlyingEvent/MPIHandler:pTmin0 3.319255e+00*GeV
set /Herwig/UnderlyingEvent/MPIHandler:pTmin0 {params[1]}*GeV
#Colour reconnection settings
set /Herwig/Hadronization/ColourReconnector:ColourReconnection Yes
set /Herwig/Hadronization/ColourReconnector:ReconnectionProbability {params[2]}pt
#Colour Disrupt settings
set /Herwig/Partons/RemnantDecayer:colourDisrupt {params[3]}
#inverse hadron radius
set /Herwig/UnderlyingEvent/MPIHandler:InvRadius {params[4]}
#MPI model settings
set /Herwig/UnderlyingEvent/MPIHandler:softInt Yes
set /Herwig/UnderlyingEvent/MPIHandler:twoComp Yes
set /Herwig/UnderlyingEvent/MPIHandler:DLmode 2

# EE5C Tune 
#set /Herwig/UnderlyingEvent/MPIHandler:EnergyExtrapolation Power
#set /Herwig/UnderlyingEvent/MPIHandler:ReferenceScale 7000.*GeV
#set /Herwig/UnderlyingEvent/MPIHandler:Power 3.705288e-01
##set /Herwig/UnderlyingEvent/MPIHandler:Power 0.33
#set /Herwig/UnderlyingEvent/MPIHandler:pTmin0 6.660252e+00*GeV
##Colour reconnection settings
#set /Herwig/Hadronization/ColourReconnector:ColourReconnection Yes
##set /Herwig/Hadronization/ColourReconnector:ReconnectionProbability 0.49
#set /Herwig/Hadronization/ColourReconnector:ReconnectionProbability 5.278926e-01
##Colour Disrupt settings
##set /Herwig/Partons/RemnantDecayer:colourDisrupt 0.80 			  
#set /Herwig/Partons/RemnantDecayer:colourDisrupt 6.284222e-01
##inverse hadron radius
##set /Herwig/UnderlyingEvent/MPIHandler:InvRadius 2.30
#set /Herwig/UnderlyingEvent/MPIHandler:InvRadius 2.254998e+00	 		  
##MPI model settings
#set /Herwig/UnderlyingEvent/MPIHandler:softInt Yes
#set /Herwig/UnderlyingEvent/MPIHandler:twoComp Yes
#set /Herwig/UnderlyingEvent/MPIHandler:DLmode 2

# MinBias 
#cd /Herwig/MatrixElements/
#insert SimpleQCD:MatrixElements[0] MEMinBias
##tell the MPI handler that the primary subprocess is QCD jet production
#set /Herwig/UnderlyingEvent/MPIHandler:IdenticalToUE 0
##Cuts for soft-inclusive events
#cd /Herwig/Cuts
#set JetKtCut:MinKT 0.0*GeV 
#set QCDCuts:MHatMin 0.0*GeV 
#set QCDCuts:X1Min 0.01
#set QCDCuts:X2Min 0.01 

# Dijet QCD 
#mkdir /Herwig/Weights
#cd /Herwig/Weights
#create ThePEG::ReweightMinPT reweightMinPT ReweightMinPT.so
cd /Herwig/MatrixElements/
insert SubProcess:MatrixElements[0] MEQCD2to2
#insert SimpleQCD:Preweights[0] /Herwig/Weights/reweightMinPT
#cd /
#set /Herwig/Weights/reweightMinPT:Power 5.5
#set /Herwig/Weights/reweightMinPT:Scale 5*GeV
cd /Herwig/Cuts/
set JetKtCut:MinKT 2.0*GeV
#set JetKtCut:MaxKT 60.0*GeV
# the QCDCuts apparently doesnt exist 
#set QCDCuts:MHatMin 0.0*GeV
#set /Herwig/UnderlyingEvent/MPIHandler:IdenticalToUE 0

# Setting up Hadronization/MPI 
cd /Herwig/EventHandlers
#set EventHandler:CascadeHandler:MPIHandler NULL
#set EventHandler:DecayHandler NULL
#set EventHandler:HadronizationHandler NULL
#set EventHandler:CascadeHandler NULL

# for HepMC Output 
insert /Herwig/Generators/EventGenerator:AnalysisHandlers[0] /Herwig/Analysis/HepMCFile
set /Herwig/Analysis/HepMCFile:PrintEvent 1000
set /Herwig/Analysis/HepMCFile:Format GenEvent
set /Herwig/Analysis/HepMCFile:Units GeV_mm

# plot the event 
#cd /Herwig/Generators
#set /Herwig/Analysis/Plot:EventNumber 5
#insert EventGenerator:AnalysisHandlers 0 /Herwig/Analysis/Plot

# Create the run script 
do /Herwig/MatrixElements/Matchbox/Factory:ProductionMode
cd /Herwig/Generators
saverun RHIC_Dijet_10M EventGenerator
"""

    out = open(f"{dir}/RHIC_Dijet.in", 'w')
    out.writelines(RHIC_Diject)
    out.close()

    file.close()    