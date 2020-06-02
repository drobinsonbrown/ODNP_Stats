for i in {0..45}
do
    sshpass -p 'Drob101618776' scp drobins@pod.cnsi.ucsb.edu:/home/drobins/PEO_TEMPO_GPU/ODNP_PAPER/4HTEMPO/OPC4/T290/nCorr/autocorr* ./nve$i

done
