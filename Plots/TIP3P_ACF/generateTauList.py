import numpy as np
import matplotlib.pyplot as plt
import glob
import sys,os
import re

def tryint(s):
    try:
        return int(s)
    except ValueError:
        return s

def alphanum_key(s):
    """ turn a string into a list of string and a number chunks.
        " z23a" -> ["z",23,"a"]
    """
    #pattern='([0-9]+)'
    pattern='(\d{3})\D' 
    return [ tryint(c) for c in re.split(pattern,s) ]

def main(args):    
    try:
        prefix = args[0]
    except IndexError:
        prefix = 'Water'
    
    #if prefix=='Water':
    #    kinds = ['*number*','*P1*','*P2*','*HB*']
    #else:
    kinds = ['*number*','*P1*','*P2*','*HB*','*F0*']
    #kinds = ['*number*','*P1*','*P2*','*HB*']
    
    taus = []
    for kind in kinds:
        names = glob.glob('tau*'+kind+'_'+prefix+'*.txt')
        tau = np.zeros(len(names))

        # define a list to add names to
        namelist = []
        
        for i,name in enumerate(names):
            namelist.append(name)

        namelist.sort(key=alphanum_key)

        for i,name in enumerate(namelist):
            tau[i]=np.loadtxt(name)[1]
        
        taus.append(tau)
    taus=np.hstack(taus)
    np.savetxt('tau'+prefix+'.txt',taus)
    
if __name__ == "__main__":
  main(sys.argv[1:])
