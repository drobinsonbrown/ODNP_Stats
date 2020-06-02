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

    cd ODNP
    python $scriptDIR/bootstrap_ACF.py autocorrF0
    python $scriptDIR/bootstrap_ODNP.py autocorrF0
    cp avg* $original/avgautocorrF0_${INFILE}_$i.txt
    cp tau_autocorrF0.txt $original/tau_autocorrF0_${INFILE}_$i.txt
    cp xi.txt $original/xi_${INFILE}_$i.txt
    cp k_sigma.txt $original/k_sigma${INFILE}_$i.txt

    cd $original
done

