import numpy as np
import matplotlib.pyplot as plt
import glob
import sys,os
import re

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
        
    names = glob.glob(prefix+'_'+'Water'+'*.txt')
    length = len(names)
    T = np.zeros(length)
    tau = np.zeros(length)
    tau_low = np.zeros(length)
    tau_high = np.zeros(length)
    for i,name in enumerate(names):
        pattern = '(\d{3})\D'
        T[i] = float(re.findall(pattern, name)[0])
        data = np.loadtxt(name)
        tau[i] = data[1]
      
        tau_low[i] = data[0]
        tau_high[i] = data[2]

    newInd = np.argsort(T)
    T = T[newInd]
    tau = tau[newInd]
    tau_high = tau_low[newInd]
    tau_low = tau_low[newInd]

    fig, ax1 = plt.subplots()
 
    for i,x in enumerate(T):
        
        ax1.errorbar(x,tau[i],yerr=[[tau_low[i]],[tau_high[i]]],fmt='r.-',
                     label='Bulk Water',capsize=2)
    ax1.plot(T,tau,'r.--',label=prefix)
    ax1.set_xlabel('T [K]',fontsize=14)
    ax1.set_ylabel(yaxis,fontsize=14)
    #ax1.hold(True)
    
    names = glob.glob(prefix+'_'+'4HTEMPO'+'*.txt')
    length = len(names)
    T = np.zeros(length)
    tau = np.zeros(length)
    tau_low = np.zeros(length)
    tau_high = np.zeros(length)
    for i,name in enumerate(names):
        pattern = '(\d{3})\D'
        T[i] = float(re.findall(pattern, name)[0])
        data = np.loadtxt(name)
        tau[i] = data[1]
      
        tau_low[i] = data[0]
        tau_high[i] = data[2]

    newInd = np.argsort(T)
    T = T[newInd]
    tau = tau[newInd]
    tau_high = tau_low[newInd]
    tau_low = tau_low[newInd]

 
    for i,x in enumerate(T):
        
        ax1.errorbar(x,tau[i],yerr=[[tau_low[i]],[tau_high[i]]],fmt='b.-',
                     label=prefix,capsize=2)

    ax1.plot(T,tau,'b.--',label=prefix)
    ax1.set_xlabel('T [K]',fontsize=14)
    ax1.set_ylabel(yaxis,fontsize=14)
    #plt.hold(True)
    
    names = glob.glob(prefix+'_'+'PEO12'+'*.txt')
    length = len(names)
    T = np.zeros(length)
    tau = np.zeros(length)
    tau_low = np.zeros(length)
    tau_high = np.zeros(length)
    for i,name in enumerate(names):
        pattern = '(\d{3})\D'
        T[i] = float(re.findall(pattern, name)[0])
        data = np.loadtxt(name)
        tau[i] = data[1]
      
        tau_low[i] = data[0]
        tau_high[i] = data[2]

    newInd = np.argsort(T)
    T = T[newInd]
    tau = tau[newInd]
    tau_high = tau_low[newInd]
    tau_low = tau_low[newInd]

    for i,x in enumerate(T):
        
        ax1.errorbar(x,tau[i],yerr=[[tau_low[i]],[tau_high[i]]],fmt='g.-',
                     label=prefix,capsize=2)

    ax1.plot(T,tau,'g.--',label=prefix)
    ax1.set_xlabel('T [K]',fontsize=14)
    ax1.set_ylabel(yaxis,fontsize=14)
    
    plt.savefig(prefix+'.png')
    plt.savefig(prefix+'.eps',format='eps')
    plt.show()
    
if __name__ == "__main__":
  main(sys.argv[1:])
