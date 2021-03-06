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
    
    xis = []
    names = glob.glob('xi_'+prefix+'*.txt')
    xi = np.zeros(len(names))

    # define a list to add names to
    namelist = []
        
    for i,name in enumerate(names):
        namelist.append(name)

    namelist.sort(key=alphanum_key)

    for i,name in enumerate(namelist):
        xi[i]=np.loadtxt(name)[1]
        
    xis.append(xi)
    xis=np.hstack(xis)
    np.savetxt('xi_'+prefix+'.txt',xis)
    
if __name__ == "__main__":
  main(sys.argv[1:])
