#!/bin/bash 

echo "Script for testing lordec programs"
echo $PWD

progCor="lordec-correct"
progStat="lordec-stat"
progGraph="lordec-build-SR-graph"

dataDir="./DATA"
resDir="RES"
pDir="../build/tools"

# Short Read files
SRfileL="ill-test-5K-1.fa"
SRfileR="ill-test-5K-2.fa"
twoSRfiles=${SRfileR},${SRfileL}

# SR input as file containing SR filenames
metaSRfile="my2-ill-test-files.txt"

# LR input
PACfile="pacbio-test.fa.gz"

# parameters
k=19
thr=3

#result files
resSingle=../${resDir}/"pacbio-corrected-single"
resDouble=../${resDir}/"pacbio-corrected-double"
resMeta=../${resDir}/"pacbio-corrected-meta"
builtGraph=../${resDir}/"my2-ill-test-graph"

if [ -d ${resDir} ] ; then
    echo " results directory exists ", ${resDir}
    ls ${resDir}

    pushd ${resDir}
    if [ -f *.log ] ; then
	echo "cleaning existing log files"
	rm *.log
    fi
    if [ -f *.fa ] ; then
	echo "cleaning existing result files"
	rm *.fa
    fi
    popd
else # if [ ! -d ${resDir} ] ; then
    echo "creating results directory ", ${resDir}
    mkdir ${resDir}
fi

pushd ${dataDir}

# test with single SR file
${pDir}/${progCor} -2 ${SRfileL} -k ${k} -s ${thr} -i ${PACfile} -o ${resSingle}.fa &>  ${resSingle}.log
if [ ! $? ] ; then
    echo "Pb correction single file: ", $?,  ${SRfileL}, ${resSingle}.log
else
    echo "Correction from single SR file Ok; created ", ${resSingle}.fa, ${resSingle}.log, *.h5
    #rm *.h5
fi

# test with double SR file
${pDir}/${progCor} -2 ${twoSRfiles} -k ${k} -s ${thr} -i ${PACfile} -o ${resDouble}.fa &>  ${resDouble}.log
if [ ! $? ] ; then
    echo "Pb correction double file: ", $?,  ${twoSRfiles}, ${resDouble}.log
else
    echo "Correction from double SR file Ok; created ", ${resDouble}.fa, ${resDouble}.log, *.h5
    #rm *.h5
fi

# test with meta SR file
${pDir}/${progCor} -2 ${metaSRfile} -k ${k} -s ${thr} -i ${PACfile} -o ${resMeta}.fa &>  ${resMeta}.log
if [ ! $? ] ; then
    echo "Pb correction meta file: ", $?,  ${metaSRfile}, ${resMeta}.log
else
    echo "Correction from meta SR file Ok; created ", ${resMeta}.fa, ${resMeta}.log, *.h5
    #rm *.h5
fi

# build graph test
${pDir}/${progGraph} -2 ${metaSRfile} -k ${k} -s ${thr} -g ${builtGraph}.h5 &>  ${builtGraph}.log
if [ ! $? ] ; then
    echo "Pb building graph from meta file: ", $?,  ${metaSRfile}, ${builtGraph}.log, *.h5
else
    echo "Building graph from meta SR file Ok; created ", *.h5
    #rm *.h5
fi

popd
echo "End of tests"
exit 0


