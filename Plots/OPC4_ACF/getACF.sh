#!/bin/bash

INFILE=${1}

for i in {265,275,290,310,340}
do
    #echo ${INFILE}
    original=$(pwd)
    cd ../../${INFILE}/OPC4/T$i/analysis/
    cwd=$(pwd)
    scriptDIR=/home/dennis/Documents/PhD/Winter20/ODNP_PAPER/scripts
    echo $scriptDIR

    cd Hbonds
    python $scriptDIR/bootstrap_ACF.py HBCorr
    cp avg* $original/avgHBCorr_${INFILE}_$i.txt
    cp tau_HBCorr.txt $original/tau_HBCorr_${INFILE}_$i.txt
    
    cd $cwd/Survival
    python $scriptDIR/bootstrap_ACF.py numberCorr
    python $scriptDIR/bootstrap_ACF.py weightedCorr
    cp avgnumber* $original/avgnumberCorr_${INFILE}_$i.txt
    cp avgweighted* $original/avgweightedCorr_${INFILE}_$i.txt
    cp tau_numberCorr.txt $original/tau_numberCorr_${INFILE}_$i.txt
    cp tau_weightedCorr.txt $original/tau_weightedCorr_${INFILE}_$i.txt
    
    cd $cwd/OACF
    python $scriptDIR/bootstrap_ACF.py OrientCorr_P1
    python $scriptDIR/bootstrap_ACF.py OrientCorr_P2
    cp avgOrientCorr_P1.txt $original/avgOrientCorr_P1_${INFILE}_$i.txt
    cp avgOrientCorr_P2.txt $original/avgOrientCorr_P2_${INFILE}_$i.txt
    cp tau_OrientCorr_P1.txt $original/tau_OrientCorr_P1_${INFILE}_$i.txt
    cp tau_OrientCorr_P2.txt $original/tau_OrientCorr_P2_${INFILE}_$i.txt
    
    cd $original
done

