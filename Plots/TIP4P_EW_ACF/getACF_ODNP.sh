#!/bin/bash

INFILE=${1}

#echo ${INFILE}
original=$(pwd)
cd ../../${INFILE}/TIP4P_EW/T290/analysis/
cwd=$(pwd)
scriptDIR=/home/dennis/Documents/PhD/Winter20/ODNP_PAPER/scripts
echo $scriptDIR

cd ODNP
python $scriptDIR/bootstrap_ACF.py autocorrF0
python $scriptDIR/bootstrap_ODNP.py autocorrF0
cp avg* $original/avgautocorrF0_${INFILE}_290.txt
cp tau_autocorrF0.txt $original/tau_autocorrF0_${INFILE}_290.txt
cp xi.txt $original/xi_${INFILE}_290.txt
cp k_sigma.txt $original/k_sigma${INFILE}_290.txt

cd $original

