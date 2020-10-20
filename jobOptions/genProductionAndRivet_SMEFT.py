from MadGraphControl.MadGraphUtils import *
from os.path import join as pathjoin 
import os
import sys
from Rivet_i.Rivet_iConf import Rivet_i

fcard = open('proc_card_mg5.dat','w')

safe_factor=2.0
nevents=safe_factor*runArgs.maxEvents
runName='run_01'

#==================================================================================

# Merging parameters in MadGraph. removed xqcut
ickkw=0

ptj=20.0
ptb=20.0

ptl=0.0
etal=10

drll=0.05
drjj=0.05
drbb=0.05
drbj=0.05

pdfwgt='T'

# Additional merging parameters in Pythia
ktdurham=30
dparameter=0.4
nJetMax=2 
process='pp>h'

maxjetflavor=5  # should be 5 if b included;

mbottom=0
ymb=0

#==================================================================================

if len(str(runArgs.runNumber)) != 5:
    raise RuntimeError("Not a valid run number: must have 5 digits (first digit = production mode; last 2 pairs: XXYY=BSM (19-82), XX00=interference, 0000=SM)")

#==================================================================================
# Set SM or interference or BSM term; SM run number ends with 00
NPmode=''
if runArgs.runNumber % 10000 == 0: 
    print "Simulate SM term"
    NPmode='NP^2==0 NP=1'
elif runArgs.runNumber % 100 == 0:
    print "Simulate interference term"
    NPmode='NP^2==1 NP=1'
else:
    print "Simulate BSM term"
    NPmode='NP^2==2 NP=1'

if len(NPmode)==0:
    raise RuntimeError("Do not recognise which term of the cross section to generate (SM, interference or BSM)")

#==================================================================================
# Write process according to run number
if str(runArgs.runNumber)[:1] == '1': # ggH+bbH
    print "Generate ggH+bbH events."
    processCommand="""
generate p p > h QED=1 """+NPmode+""" @0
add process p p > h j QED=1 """+NPmode+""" @1
add process p p > h j j QED=1 """+NPmode+""" @2
add process p p > h b b~ QED=1 """+NPmode+""" @3
    """

elif str(runArgs.runNumber)[:1] == '2': # VBF+VHhad
    print "Generate VBF+VHhad events."
    processCommand="""
generate p p > h j j QCD=0 """+NPmode+""" @0
    """

elif str(runArgs.runNumber)[:1] == '3': # ZHlep
    print "Generate ZHlep events."
    processCommand="""
generate p p > h l+ l- """+NPmode+""" @0
add process p p > h ta+ ta- """+NPmode+""" @1
add process p p > h vl vl~ """+NPmode+""" @2
    """

elif str(runArgs.runNumber)[:1] == '4': # WHlep
    print "Generate WHlep events."
    processCommand="""
generate p p > h l+ vl """+NPmode+""" @0
add process p p > h l- vl~ """+NPmode+""" @1
    """

elif str(runArgs.runNumber)[:1] == '5': # ttH
    print "Generate ttH events."
    processCommand="""
generate p p > h t t~ """+NPmode+""" @0
    """
 
elif str(runArgs.runNumber)[:1] == '6': # tHjb
    print "Generate tHjb events."
    processCommand="""
generate p p > h t b~ j """+NPmode+""" @0
add process p p > h t~ b j """+NPmode+""" @1
    """

elif str(runArgs.runNumber)[:1] == '7': # tHW
    print "Generate tHW events. Using 5FS, mb=0."
    ymb=0
    processCommand="""
generate p p > h t w- """+NPmode+""" @0
add process p p > h t~ w+ """+NPmode+""" @1
    """

elif str(runArgs.runNumber)[:1] == '8': # LO HH
	print "Generate HH events."
	processCommand="""
generate p p > h h QCD=0 """+NPmode+""" @0
	"""

else:
    raise RuntimeError("Not a valid production mode (first digit of run number 1-5)")

fcard.write("""
import model /pbs/home/m/mbelfkir/EFT/SMEFTsim_A_U35_MwScheme_UFO_v2_1-massless
set gauge unitary
define p = p b b~"""
+processCommand+"""
output -f""")
fcard.close()


#==================================================================================

# Which parameters to modify (default: 0)
execfile('/pbs/home/m/mbelfkir/EFT/jobOptions/EFTparameterVariations.py')
non_zero_params = dict()
non_zero_params.update(varyEFTParam(int(str(runArgs.runNumber)[1:3]),1.0))
non_zero_params.update(varyEFTParam(int(str(runArgs.runNumber)[3:5]),1.0))
print non_zero_params

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
           'cut_decays':'F', 
           'pdlabel':'lhapdf',
           'lhaid':'90400',
#           'lhaid':'260000',
           'ickkw'       : ickkw,
           'drll'        : drll,
           'drjj'        : drjj,
           'drbb'        : drbb,
           'drbj'        : drbj,
           'ptj'         : ptj,
           'ptb'         : ptb,
           'ptl'         : ptl,
           'etal'        : etal,
           'pdfwgt'      : pdfwgt,
           'maxjetflavor': maxjetflavor,
           'dparameter'  : dparameter,
           'ktdurham'    : ktdurham,
           'use_syst'    : "False"
}
 
#==================================================================================
process_dir = new_process()

#==================================================================================

build_run_card(run_card_old=get_default_runcard(proc_dir=process_dir),run_card_new='run_card.dat',xqcut=0.0,nevts=nevents,rand_seed=runArgs.randomSeed,beamEnergy=beamEnergy,extras=extras)

print_cards()

#==================================================================================

param_card=open('param_card.dat','w')
param_card_default=open('/pbs/home/m/mbelfkir/EFT/jobOptions/param_card_SMEFTsim_A_U35_MwScheme_UFO_v2_1-massless_default.dat','r')
"""
for line in param_card_default:
    if len(line.strip())==0 or line.strip()[0]=='#':
        param_card.write(line)
        continue
    spl=line.split('#')
    if len(spl)>1 and spl[1].strip().lower() in [p.lower() for p in non_zero_params]:
        param_card.write(spl[0].split()[0]+' '+str(non_zero_params[p])+' #'+spl[1])
    elif len(spl)>1 and 'MB' in spl[1]:
        param_card.write('    5 '+str(mbottom)+' # '+spl[1])
    elif len(spl)>1 and 'ymb' in spl[1]:
        param_card.write('    5 '+str(ymb)+' # '+spl[1])
    else:
        param_card.write(line)

param_card.close()
param_card_default.close()
"""
generate(run_card_loc='run_card.dat',param_card_loc='param_card.dat',proc_dir=process_dir,run_name=runName,nevents=nevents,random_seed=runArgs.randomSeed)
arrange_output(run_name=runName,proc_dir=process_dir,outputDS=runName+'._00001.events.tar.gz',lhe_version=3,saveProcDir=True)  

#==================================================================================
# Shower 
evgenConfig.description = 'MadGraphSMEFT_Higgs'
evgenConfig.keywords+=['Higgs']
evgenConfig.inputfilecheck = runName
runArgs.inputGeneratorFile=runName+'._00001.events.tar.gz'

include("MC15JobOptions/Pythia8_A14_NNPDF23LO_EvtGen_Common.py")
include("MC15JobOptions/Pythia8_MadGraph.py")

if str(runArgs.runNumber)[:1] == '1':
    PYTHIA8_nJetMax=nJetMax
    PYTHIA8_Dparameter=dparameter
    PYTHIA8_Process=process
    PYTHIA8_TMS=ktdurham
    PYTHIA8_nQuarksMerge=maxjetflavor
    include("MC15JobOptions/Pythia8_CKKWL_kTMerge.py")
    genSeq.Pythia8.Commands+=["Merging:mayRemoveDecayProducts=on"]


histoFileName="rivetFile"
from GaudiSvc.GaudiSvcConf import THistSvc
ServiceMgr += THistSvc()
ServiceMgr.THistSvc.Output = ["Rivet DATAFILE='"+histoFileName+".root' OPT='RECREATE'"]

if str(runArgs.runNumber)[:1] == '1': # ggH+bbH                                                                                                                                      
    os.environ["HIGGSPRODMODE"]="GGF"
elif str(runArgs.runNumber)[:1] == '2':# VBF                                                                                                                                         
    os.environ["HIGGSPRODMODE"]="VBF"
elif str(runArgs.runNumber)[:1] == '3':# ZHlep                                                                                                                                       
    os.environ["HIGGSPRODMODE"]="QQ2ZH"
elif str(runArgs.runNumber)[:1] == '4':# WHlep                                                                                                                                       
    os.environ["HIGGSPRODMODE"]="WH"
elif str(runArgs.runNumber)[:1] == '5':# TTH                                                                                                                                         
    os.environ["HIGGSPRODMODE"]="TTH"
#elif str(runArgs.runNumber)[:1] == '8':
#    os.environ["HIGGSPRODMODE"]="HH"
else:
    os.environ["HIGGSPRODMODE"]="TH" #TH and THW                                                                                                                                     

from Rivet_i.Rivet_iConf import Rivet_i
genSeq += Rivet_i()
genSeq.Rivet_i.AnalysisPath = os.environ['PWD']
genSeq.Rivet_i.Analyses += [ 'HiggsTemplateCrossSectionsStage12' ]
genSeq.Rivet_i.RunName = ""
genSeq.Rivet_i.HistoFile = "histoFile"


