for i in {0..9}
do
    sshpass -p 'Drob101618776' scp drobins@pod.cnsi.ucsb.edu:/home/drobins/PEO_TEMPO_GPU/ODNP_PAPER/4HTEMPO/TIP4P_EW/T290/nve$i/*.txt ./nve$i

        sshpass -p 'Drob101618776' scp drobins@pod.cnsi.ucsb.edu:/home/drobins/PEO_TEMPO_GPU/ODNP_PAPER/4HTEMPO/TIP4P_EW/T290/nve$i/*.png ./nve$i
    
done
