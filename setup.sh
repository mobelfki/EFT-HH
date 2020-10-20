export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase
source ${ATLAS_LOCAL_ROOT_BASE}/user/atlasLocalSetup.sh
asetup 21.6.10,AthGeneration,slc6
source ${LCG_RELEASE_BASE}/LCG_88/MCGenerators/rivet/${RIVETVER}/${LCG_PLATFORM}/rivetenv.sh
export PYTHONPATH=/pbs/home/m/mbelfkir/EFT/SMEFTatNLO:$PYTHONPATH

