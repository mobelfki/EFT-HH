from MadGraphControl.MadGraphUtils import *
from os.path import join as pathjoin 
import os
import sys
from Rivet_i.Rivet_iConf import Rivet_i

fcard = open('proc_card_mg5.dat','w')

safe_factor=2.
nevents=safe_factor*runArgs.maxEvents
runName='run_01'

#==================================================================================

# Merging parameters in MadGraph. removed xqcut

process='pp>hh'

#==================================================================================

if len(str(runArgs.runNumber)) != 5:
    raise RuntimeError("Not a valid run number: must have 5 digits (first digit = production mode; last 2 pairs: XXYY=BSM (19-82), XX00=interference, 0000=SM)")

#==================================================================================
# Set SM or interference or BSM term; SM run number ends with 00
NPmode=''
if runArgs.runNumber % 10000 == 0: 
    print "Simulate SM term"
    NPmode='NP^2==0 NP=2'
elif runArgs.runNumber % 100 == 0:
    print "Simulate interference term"
    NPmode='NP^2==2 NP=2'
else:
    print "Simulate BSM term"
    NPmode='NP^2==4 NP=2'

if len(NPmode)==0:
    raise RuntimeError("Do not recognise which term of the cross section to generate (SM, interference or BSM)")

#==================================================================================
# Write process according to run number
if str(runArgs.runNumber)[:1] == '1': # ggHH
    print "Generate gg-->HH events."
    processCommand="""
generate g g > h h  QED=2 QCD=0 [QCD] """+NPmode+""" @0
    """
elif str(runArgs.runNumber)[:1] == '2': #VBF
    print "Generate pp-->HHjj events."
    processCommand="""
generate p p > h h j j $$ z w+ w- / a j QED=4 """ +NPmode+""" @0
    """ 
else:
    raise RuntimeError("Not a valid production mode (first digit of run number 1-5)")

fcard.write("""
import model /pbs/home/m/mbelfkir/EFT/SMEFTatNLO-NLO
set gauge unitary
define p = p g u c d s u~ c~ d~ s~"""
+processCommand+"""
output -f""")
fcard.close()


#==================================================================================

# Which parameters to modify (default: 0)
#execfile('/pbs/home/m/mbelfkir/EFT/jobOptions/EFTparameterVariations_HH.py')
non_zero_params = dict()
#non_zero_params.update(varyEFTParam(int(str(runArgs.runNumber)[1:3]),1.0))
#non_zero_params.update(varyEFTParam(int(str(runArgs.runNumber)[3:5]),1.0))
#print non_zero_params

#==================================================================================

# General settings
beamEnergy=-999
if hasattr(runArgs,'ecmEnergy'):
    beamEnergy = runArgs.ecmEnergy / 2.
else: 
    raise RuntimeError("No center of mass energy found.")

#==================================================================================

#Fetch default LO run_card.dat and set parameters
extras = { 'lhe_version':'3.0', 
           'pdlabel':'lhapdf',
	   'use_syst'    : "False",
	   'dynamical_scale_choice' : 4.0
}
 
#==================================================================================
process_dir = new_process()

#==================================================================================

new_run_card = build_run_card(run_card_old=get_default_runcard(proc_dir=process_dir),run_card_new='new_run_card.dat',xqcut=0,nevts=nevents,rand_seed=runArgs.randomSeed,beamEnergy=beamEnergy,scalefact=0.5,extras=extras)

#==================================================================================

#==================================================================================
# check if the process is NLO

isNLO = is_NLO_run(process_dir)
print('Is NLO : ', isNLO)

#==================================================================================

param_card=open('param_card.dat','w')
param_card_default=open('/pbs/home/m/mbelfkir/EFT/jobOptions/param_card_SMEFTatNLO-NLO.dat','r')

for line in param_card_default:
	param_card.write(line)

param_card.close()
param_card_default.close()

generate(run_card_loc='new_run_card.dat',param_card_loc='param_card.dat',mode=2,njobs=32,proc_dir=process_dir,run_name=runName,nevents=nevents,random_seed=runArgs.randomSeed)

arrange_output(run_name=runName,proc_dir=process_dir,outputDS=runName+'._00001.events.tar.gz',lhe_version=3,saveProcDir=True)  

#==================================================================================
# Shower 
evgenConfig.description = 'MadGraphSMEFT_DiHiggs'
evgenConfig.keywords+=['Higgs']
evgenConfig.inputfilecheck = runName
runArgs.inputGeneratorFile=runName+'._00001.events.tar.gz'

include("MC15JobOptions/Pythia8_A14_NNPDF23LO_EvtGen_Common.py")
include("MC15JobOptions/Pythia8_MadGraph.py")

"""
if str(runArgs.runNumber)[:1] == '1':
    PYTHIA8_nJetMax=nJetMax
    PYTHIA8_Dparameter=dparameter
    PYTHIA8_Process=process
    PYTHIA8_TMS=ktdurham
    PYTHIA8_nQuarksMerge=maxjetflavor
    include("MC15JobOptions/Pythia8_CKKWL_kTMerge.py")
    genSeq.Pythia8.Commands+=["Merging:mayRemoveDecayProducts=on"]
"""

histoFileName="RivetOutput"
from GaudiSvc.GaudiSvcConf import THistSvc
ServiceMgr += THistSvc()
ServiceMgr.THistSvc.Output = ["Rivet DATAFILE='"+histoFileName+".root' OPT='RECREATE'"]

if str(runArgs.runNumber)[:1] == '1': # ggH+bbH                                                                                                                                      
    os.environ["PRODMODE"]="GGF"
elif str(runArgs.runNumber)[:1] == '2':# VBF                                                                                                                                         
    os.environ["PRODMODE"]="VBF"
else : 
    raise RuntimeError("Prodcution mode not defined ")
                                                                                                                               

from Rivet_i.Rivet_iConf import Rivet_i
genSeq += Rivet_i()
genSeq.Rivet_i.AnalysisPath = os.environ['PWD']
genSeq.Rivet_i.Analyses += [ 'DiHiggsRivet' ]
genSeq.Rivet_i.RunName = ""
genSeq.Rivet_i.HistoFile = "OutputHistos"


