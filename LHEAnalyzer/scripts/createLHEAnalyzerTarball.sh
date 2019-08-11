#!/bin/sh


echo "Running $0"

TARFILE="lheanalyzer.tar"
echo "SCRAM_ARCH: ${SCRAM_ARCH}"

HERE=$(pwd)
echo "The tarball will appear in $HERE"

pushd $CMSSW_BASE

createLHEAnalyzerInputTarball.sh "$@"
EXTRAS=""
if [[ -e lheanalyzer_inputs ]];then
  EXTRAS="lheanalyzer_inputs"
fi


tar Jcvf ${TARFILE} \
lib \
biglib \
cfipython \
config \
external \
bin \
src/ZZMatrixElement \
src/CMSDataTools \
src/MelaAnalytics \
src/LHEAnalyzer \
${EXTRAS} \
--exclude=lib/${SCRAM_ARCH}/* \
--exclude=src/LHEAnalyzer/LHEAnalyzer/test/output \
--exclude=src/LHEAnalyzer/LHEAnalyzer/test/tmp* \
--exclude=src/LHEAnalyzer/LHEAnalyzer/test/temp* \
--exclude=src/LHEAnalyzer/LHEAnalyzer/tmp* \
--exclude=src/LHEAnalyzer/LHEAnalyzer/temp* \
--exclude=src/LHEAnalyzer/LHEAnalyzer/*.sh \
--exclude=src/LHEAnalyzer/LHEAnalyzer/*.root \
--exclude=src/LHEAnalyzer/LHEAnalyzer/*.lhe \
--exclude=src/LHEAnalyzer/LHEAnalyzer/*.txt \
--exclude=src/LHEAnalyzer/LHEAnalyzer/*.lst \
--exclude=src/LHEAnalyzer/LHEAnalyzer/*.log \
--exclude=src/ZZMatrixElement/MELA/test/reference \
--exclude={.git,.gitignore,*.pyc,*.pcm,*.d,*.tar,libmcfm*,libcollier*,*.so,*.o,*.f~,*mod,mstw*.dat,*.err,*.sub,log*,Logs,DONE}

mv $TARFILE $HERE/
rm -f lheanalyzer_inputs # If the inputs tar doesn't exist, it is still ok.

popd
