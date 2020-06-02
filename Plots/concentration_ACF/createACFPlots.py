import numpy as np
import matplotlib.pyplot as plt
import glob
import sys,os
import re

def wfrac(N):
    return (800+(N-1)*44*12)/(800+(N-1)*44*12+6000*18)

def main(args):    
    try:
        prefix = args[0]
        system = args[1]
    except IndexError:
        prefix = 'avgautocorrF0'
        system = 'PEO12'
        
    # load ACF data
    if 'number' in prefix:
        yaxis = '$P_{Survival}(t)$'
    elif 'HB' in prefix:
        yaxis = '$P_{HB}(t)$'
    elif 'P1' in prefix:
        yaxis = '$C_{1}(t)$'
    elif 'P2' in prefix:
        yaxis = '$C_{2}(t)$'
    elif 'F0' in prefix:
        yaxis = '$C_{ODNP}(t)$'
    
    names = glob.glob(prefix+'_'+system+'*.txt')
    N = np.zeros(len(names))
    corr = np.zeros((len(names),np.max((np.loadtxt(names[0]).shape))))
    for i,name in enumerate(names):
        #pattern = r'\D(\d{$d})\D'
        pattern = r'-?\d+\.?\d*'
        numset = re.findall(pattern,name)
        N[i] = float(numset[len(numset)-1])
        data = np.loadtxt(name)
        if data.shape[0]>data.shape[1]:
            time = data[:,0]
            corr[i,:] = data[:,1]
        else:
            time = data[0,:]
            corr[i,:] = data[1,:]

        
    indices = np.argsort(N)
    N = N[indices]
    corr = corr[indices,:]
    for j,Num in enumerate(N):
        plt.plot(time,corr[j,:],label='w='+'{:.3f}'.format(wfrac(Num)))

    plt.xlabel('t [ps]',fontsize=14)
    plt.ylabel(yaxis,fontsize=14)
    plt.xlim(time[0]-0.3,time[len(time)-1])
    plt.ylim(0,1)
    plt.legend(frameon=False)
    plt.savefig(prefix+'_'+system+'.png',dpi=1000)
    plt.savefig(prefix+'_'+system+'.eps',format='eps')
    #plt.show()
    plt.close()

if __name__ == "__main__":
  main(sys.argv[1:])
