#!/bin/bash

for i in {0,0.01,0.03,0.1,0.3,1,3,6,8,10,15,30,60,100}
do
    #echo ${INFILE}
    original=$(pwd)
    cd ../gammaStudy/gamma_$i/
    cwd=$(pwd)
    scriptDIR=/home/dennis/Documents/PhD/Winter20/ODNP_PAPER/scripts
    echo $scriptDIR
    
    #python $scriptDIR/bootstrap_ACF.py autocorrF0
    python $scriptDIR/bootstrap_ODNP.py autocorrF0
    cp avgautocorrF0.txt $original/avgautocorrF0_$i.txt
    cp xi.txt $original/xi_$i.txt
    cp k_sigma.txt $original/k_sigma_$i.txt

    cd $original
done

