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
        prefix = '4HTEMPO'
    
    ks = []
    names = glob.glob('k_sigma_'+prefix+'*.txt')
    k = np.zeros(len(names))

    # define a list to add names to
    namelist = []
        
    for i,name in enumerate(names):
        namelist.append(name)

    
    for i,name in enumerate(namelist):
        k[i]=np.loadtxt(name)[1]
        
    ks.append(k)
    ks=np.hstack(ks)
    np.savetxt('k_sigma'+prefix+'.txt',ks)
    
if __name__ == "__main__":
  main(sys.argv[1:])
