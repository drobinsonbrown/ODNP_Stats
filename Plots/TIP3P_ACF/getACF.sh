#!/bin/bash

INFILE=${1}

echo ${INFILE}
original=$(pwd)
cd ../../${INFILE}/TIP3P/T290/analysis/
cwd=$(pwd)
scriptDIR=/home/dennis/Documents/PhD/Winter20/ODNP_PAPER/scripts
echo $scriptDIR

cd Hbonds
python $scriptDIR/bootstrap_ACF.py HBCorr
cp avg* $original/avgHBCorr_${INFILE}_290.txt
cp tau_HBCorr.txt $original/tau_HBCorr_${INFILE}_290.txt

cd $cwd/Survival
python $scriptDIR/bootstrap_ACF.py numberCorr
cp avg* $original/avgnumberCorr_${INFILE}_290.txt
cp tau_numberCorr.txt $original/tau_numberCorr_${INFILE}_290.txt

cd $cwd/OACF
python $scriptDIR/bootstrap_ACF.py OrientCorr_P1
python $scriptDIR/bootstrap_ACF.py OrientCorr_P2
cp avgOrientCorr_P1.txt $original/avgOrientCorr_P1_${INFILE}_290.txt
cp avgOrientCorr_P2.txt $original/avgOrientCorr_P2_${INFILE}_290.txt
cp tau_OrientCorr_P1.txt $original/tau_OrientCorr_P1_${INFILE}_290.txt
cp tau_OrientCorr_P2.txt $original/tau_OrientCorr_P2_${INFILE}_290.txt

cd $original

