import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import glob
import sys,os
import re

matplotlib.rcParams.update({'figure.figsize':[8.0,6.0]})
matplotlib.rcParams.update({'font.size':18.0})



def main(args):    
    try:
        prefix = args[0]
    except IndexError:
        prefix = 'tau_autocorrF0'

    yaxis = r'$k_{\sigma}$'
        
    names = glob.glob(prefix+'_'+'Water'+'*.txt')
    length = len(names)
    T = np.zeros(length)
    k_sigma = np.zeros(length)
    k_sigma_low = np.zeros(length)
    k_sigma_high = np.zeros(length)
    for i,name in enumerate(names):
        pattern = '(\d{3})\D'
        T[i] = float(re.findall(pattern, name)[0])
        data = np.loadtxt(name)
        k_sigma[i] = data[1]
      
        k_sigma_low[i] = data[0]
        k_sigma_high[i] = data[2]

    newInd = np.argsort(T)
    T = T[newInd]
    k_sigma = k_sigma[newInd]
    k_sigma_high = k_sigma_low[newInd]
    k_sigma_low = k_sigma_low[newInd]

    fig, ax1 = plt.subplots()
    
    names = glob.glob(prefix+'4HTEMPO'+'*.txt')
    print(names)
    length = len(names)
    T = np.zeros(length)
    k_sigma = np.zeros(length)
    k_sigma_low = np.zeros(length)
    k_sigma_high = np.zeros(length)
    for i,name in enumerate(names):
        pattern = '(\d{3})\D'
        T[i] = float(re.findall(pattern, name)[0])
        data = np.loadtxt(name)
        k_sigma[i] = data[1]
      
        k_sigma_low[i] = data[0]
        k_sigma_high[i] = data[2]

    newInd = np.argsort(T)
    T = T[newInd]
    k_sigma = k_sigma[newInd]
    k_sigma_high = k_sigma_low[newInd]
    k_sigma_low = k_sigma_low[newInd]
   
    output = np.transpose(np.vstack((T,k_sigma,k_sigma_low,k_sigma_high)))     
    np.savetxt('KvT_4HTEMPO.txt',output)
 
    for i,x in enumerate(1/T):
        
        ax1.errorbar(x,k_sigma[i],yerr=[[k_sigma_low[i]],[k_sigma_high[i]]],fmt='b.-',
                     label=prefix,capsize=2)

    ax1.plot(1/T,k_sigma,'b.--',label=prefix)
    ax1.set_xlabel('1/T [$K^{-1}$]')
    ax1.set_ylabel(yaxis)
    #plt.hold(True)
    
    names = glob.glob(prefix+'PEO12'+'*.txt')
    length = len(names)
    T = np.zeros(length)
    k_sigma = np.zeros(length)
    k_sigma_low = np.zeros(length)
    k_sigma_high = np.zeros(length)
    for i,name in enumerate(names):
        pattern = '(\d{3})\D'
        T[i] = float(re.findall(pattern, name)[0])
        data = np.loadtxt(name)
        k_sigma[i] = data[1]
      
        k_sigma_low[i] = data[0]
        k_sigma_high[i] = data[2]

    newInd = np.argsort(T)
    T = T[newInd]
    k_sigma = k_sigma[newInd]
    k_sigma_high = k_sigma_low[newInd]
    k_sigma_low = k_sigma_low[newInd]

    output = np.transpose(np.vstack((T,k_sigma,k_sigma_low,k_sigma_high)))     
    np.savetxt('KvT_PEO12.txt',output)
    
    for i,x in enumerate(1/T):
        
        ax1.errorbar(x,k_sigma[i],yerr=[[k_sigma_low[i]],[k_sigma_high[i]]],fmt='g.-',
                     label=prefix,capsize=2)

    ax1.plot(1/T,k_sigma,'g.--',label=prefix)
    ax1.set_xlabel('1/T [$K^{-1}$]')
    ax1.set_ylabel(yaxis)
    
    plt.savefig(prefix+'.png',dpi=1000)
    plt.savefig(prefix+'.eps',format='eps')
    plt.show()
    
if __name__ == "__main__":
  main(sys.argv[1:])
