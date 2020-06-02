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
        prefix = 'tau_autocorrF0'
        
    # load ACF data
    if 'number' in prefix:
        yaxis = r'$\tau_{Survival}(t)$'
    elif 'HB' in prefix:
        yaxis = r'$\tau_{HB}(t)$'
    elif 'P1' in prefix:
        yaxis = r'$\tau_{1}(t)$'
    elif 'P2' in prefix:
        yaxis = r'$\tau_{2}(t)$'
    elif 'F0' in prefix:
        yaxis = r'$\tau_{ODNP}(t)$'
        
    names = glob.glob(prefix+'_'+'PEO12'+'*.txt')
    length = len(names)
    N = np.zeros(length)
    tau = np.zeros(length)
    tau_low = np.zeros(length)
    tau_high = np.zeros(length)
    for i,name in enumerate(names):
        pattern = '-?\d+\.?\d*'
        setnums = re.findall(pattern, name)
        N[i] = float(setnums[len(setnums)-1])
        data = np.loadtxt(name)
        tau[i] = data[1]
      
        tau_low[i] = data[0]
        tau_high[i] = data[2]

    newInd = np.argsort(N)
    N = N[newInd]
    tau = tau[newInd]
    tau_high = tau_low[newInd]
    tau_low = tau_low[newInd]

    fig, ax1 = plt.subplots()

    print(tau)
    output = np.transpose(np.vstack((wfrac(N)
                                     ,tau,tau_low,tau_high)))
    print(output)
    np.savetxt('TauvW.txt',output)
    
    for i,x in enumerate(wfrac(N)):
        
        ax1.errorbar(x,tau[i],yerr=[[tau_low[i]],[tau_high[i]]],fmt='g.-'
                     ,capsize=2)
    ax1.plot(wfrac(N),tau,'g.--',label=prefix)
    ax1.set_xlabel('weight fraction',fontsize=14)
    ax1.set_ylabel(yaxis,fontsize=14)
    #ax1.hold(True)

    plt.savefig(prefix+'.png',dpi=1000)
    plt.savefig(prefix+'.eps',format='eps')
    plt.show()
        
if __name__ == "__main__":
  main(sys.argv[1:])
