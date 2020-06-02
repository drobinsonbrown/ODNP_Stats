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
    except IndexError:
        prefix = 'k_sigma'
        
    yaxis = r'$k_{\sigma}$'
         
    names = glob.glob(prefix+'_'+'PEO12'+'*.txt')
    length = len(names)
    N = np.zeros(length)
    k_sigma = np.zeros(length)
    k_sigma_low = np.zeros(length)
    k_sigma_high = np.zeros(length)
    for i,name in enumerate(names):
        pattern = '-?\d+\.?\d*'
        setnums = re.findall(pattern, name)
        N[i] = float(setnums[len(setnums)-1])
        data = np.loadtxt(name)
        k_sigma[i] = data[1]
      
        k_sigma_low[i] = data[0]
        k_sigma_high[i] = data[2]

    newInd = np.argsort(N)
    N = N[newInd]
    k_sigma = k_sigma[newInd]
    print(k_sigma)
    k_sigma_high = k_sigma_low[newInd]
    k_sigma_low = k_sigma_low[newInd]

    fig, ax1 = plt.subplots()

    output = np.transpose(np.vstack((wfrac(N)
                                     ,k_sigma,k_sigma_low,k_sigma_high)))
    print(output)
    np.savetxt('KvW.txt',output)
    
    for i,x in enumerate(wfrac(N)):
        
        ax1.errorbar(x,k_sigma[i],yerr=[[k_sigma_low[i]],[k_sigma_high[i]]],fmt='g.-'
                     ,capsize=2)
    ax1.plot(wfrac(N),k_sigma,'g.--',label=prefix)
    ax1.set_xlabel('weight fraction',fontsize=14)
    ax1.set_ylabel(yaxis,fontsize=14)
    #ax1.hold(True)

    plt.savefig(prefix+'.png',dpi=1000)
    plt.savefig(prefix+'.eps',format='eps')
    plt.show()
        
if __name__ == "__main__":
  main(sys.argv[1:])
