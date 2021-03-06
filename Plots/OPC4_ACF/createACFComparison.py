import numpy as np
import matplotlib.pyplot as plt
import glob
import sys,os
import re

def main(args):    
    try:
        prefix = args[0]
    except IndexError:
        prefix = 'avgautocorrF0'
        
    # load ACF data
    if ('number' in prefix or 'weight' in prefix):
        yaxis = '$P_{Survival}(t)$'
    elif 'HB' in prefix:
        yaxis = '$P_{HB}(t)$'
    elif 'P1' in prefix:
        yaxis = '$C_{1}(t)$'
    elif 'P2' in prefix:
        yaxis = '$C_{2}(t)$'
    elif 'F0' in prefix:
        yaxis = '$C_{ODNP}(t)$'
    
    names = glob.glob(prefix+'*290.txt')
    for name in names:
        data = np.loadtxt(name)
        if data.shape[0]>data.shape[1]:
            time = data[:,0]
            corr = data[:,1]
        else:
            time = data[0,:]
            corr = data[1,:]

        if 'Water' in name:
            plt.plot(time,corr,'r',label='Water')
        elif '4HTEMPO' in name:
            plt.plot(time,corr,'b',label='4-OH-TEMPO')
        else:
            plt.plot(time,corr,'g',label='PEO12')

    plt.xlabel('t [ps]')
    plt.ylabel(yaxis)
    plt.xlim(time[0]-0.3,time[len(time)-1]+1)
    plt.ylim(0,1)
    plt.legend(frameon=False)
    plt.savefig(prefix+'Comparison'+'.png',dpi=1000)
    plt.savefig(prefix+'Comparison'+'.eps',format='eps')
    plt.show()
    plt.close()

if __name__ == "__main__":
  main(sys.argv[1:])
