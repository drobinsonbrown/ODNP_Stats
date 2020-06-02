for i in {0..9}
do

    cp ../nve$i/numberCorr.txt ./numberCorr_$i.txt
    cp ../nve$i/HBCorr.txt ./HBCorr_$i.txt
    cp ../nve$i/OrientCorr_P1.txt ./OrientCorr_P1_$i.txt
    cp ../nve$i/OrientCorr_P2.txt ./OrientCorr_P2_$i.txt
    cp ../nve$i/autocorrF0.txt ./autocorrF0_$i.txt
    cp ../nve$i/autocorrF1.txt ./autocorrF1_$i.txt
    cp ../nve$i/autocorrF2.txt ./autocorrF2_$i.txt
    mv numberCorr* ./Survival
    mv HBCorr* ./Hbonds
    mv OrientCorr* ./OACF
    mv autocorr* ./ODNP

done
