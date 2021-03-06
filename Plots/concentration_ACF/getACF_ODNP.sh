#!/bin/bash

INFILE=${1}

for i in {1,3,10,50,100}
do
    echo ${INFILE}
    original=$(pwd)
    cd ../../${INFILE}/OPC4/T290/concentration/N$i/analysis/
    cwd=$(pwd)
    scriptDIR=/home/dennis/Documents/PhD/Winter20/ODNP_PAPER/scripts
    echo $scriptDIR

    cd ODNP
    
    echo $i
    python $scriptDIR/bootstrap_ACF.py autocorrF0
    python $scriptDIR/bootstrap_ODNP.py autocorrF0
    
    #python $scriptDIR/bootstrap_ACF.py autocorrF1
    #python $scriptDIR/bootstrap_ODNP.py autocorrF1

    #python $scriptDIR/bootstrap_ACF.py autocorrF2
    #python $scriptDIR/bootstrap_ODNP.py autocorrF2

    cp avgautocorrF0.txt $original/avgautocorrF0_${INFILE}_$i.txt
    #cp avgautocorrF1.txt $original/avgautocorrF1_${INFILE}_$i.txt
    #cp avgautocorrF2.txt $original/avgautocorrF2_${INFILE}_$i.txt
    cp tau_autocorrF0.txt $original/tau_autocorrF0_${INFILE}_$i.txt
    cp xi.txt $original/xi_${INFILE}_$i.txt
    cp k_sigma.txt $original/k_sigma_${INFILE}_$i.txt

    cd $original
done

