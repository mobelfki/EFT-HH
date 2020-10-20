export runNumber=10000
export seed=1234
export maxEvent=10000
export jobOption=/pbs/home/m/mbelfkir/EFT/jobOptions/genProductionAndRivet_HH.py

if [ -d "${runNumber}" ]
then 
	echo "${runNumber} exist! "
	rm -r ${runNumber}; mkdir ${runNumber}; echo "$runNumber re-created! "
else
	mkdir ${runNumber}; echo "$runNumber created! "
fi
cd ${runNumber}
pwd
cp /pbs/home/m/mbelfkir/EFT/jobOptions/DiHiggsRivet.* ./
rivet-buildplugin RivetDiHiggsRivet.so DiHiggsRivet.cc

Generate_tf.py --ecmEnergy=13000. --maxEvents=$maxEvent --runNumber=$runNumber --firstEvent=1 --randomSeed=$seed --outputEVNTFile=mc.${runNumber}_${maxEvent}.EVNT.root --jobConfig=$jobOption


