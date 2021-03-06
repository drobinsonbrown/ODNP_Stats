import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import glob
import sys,os
import re

matplotlib.rcParams.update({'figure.figsize':[8.0,6.0]})
matplotlib.rcParams.update({'font.size':18.0})
matplotlib.rcParams.update({'xtick.labelsize':14.0})
matplotlib.rcParams.update({'ytick.labelsize':14.0})

def main(args):    
    try:
        prefix = args[0]
    except IndexError:
        prefix = 'xi'

    yaxis = r'$\xi$'

    names = glob.glob(prefix+'_'+'4HTEMPO'+'*.txt')
    length = len(names)
    T = np.zeros(length)
    xi = np.zeros(length)
    xi_low = np.zeros(length)
    xi_high = np.zeros(length)
    for i,name in enumerate(names):
        pattern = '(\d{3})\D'
        T[i] = float(re.findall(pattern, name)[0])
        data = np.loadtxt(name)
        xi[i] = data[1]
      
        xi_low[i] = data[0]
        xi_high[i] = data[2]

    newInd = np.argsort(T)
    T = T[newInd]
    xi = xi[newInd]
    xi_high = xi_low[newInd]
    xi_low = xi_low[newInd]


    fig, ax1 = plt.subplots()
 
    #for i,x in enumerate(T):
        
        #ax1.errorbar(x,xi[i],yerr=[[xi_low[i]],[xi_high[i]]],fmt='r.-',
        #             label='Bulk Water',capsize=2)
    #ax1.plot(T,xi,'r.--',label=prefix)
    #ax1.set_xlabel('T [K]',fontsize=14)
    #ax1.set_ylabel(yaxis,fontsize=14)
    #ax1.hold(True)
    
    names = glob.glob(prefix+'_'+'4HTEMPO'+'*.txt')
    length = len(names)
    T = np.zeros(length)
    xi = np.zeros(length)
    xi_low = np.zeros(length)
    xi_high = np.zeros(length)
    for i,name in enumerate(names):
        pattern = '(\d{3})\D'
        T[i] = float(re.findall(pattern, name)[0])
        data = np.loadtxt(name)
        xi[i] = data[1]
      
        xi_low[i] = data[0]
        xi_high[i] = data[2]

    newInd = np.argsort(T)
    T = T[newInd]
    xi = xi[newInd]
    xi_high = xi_low[newInd]
    xi_low = xi_low[newInd]

    output = np.transpose(np.vstack((T,xi,xi_low,xi_high)))     
    np.savetxt('XivT_4HTEMPO.txt',output)
    
    for i,x in enumerate(1/T):
        
        ax1.errorbar(x,xi[i],yerr=[[xi_low[i]],[xi_high[i]]],fmt='b.-',
                     label=prefix,capsize=2)

    ax1.plot(1/T,xi,'b.--',label=prefix)
    ax1.set_xlabel('1/T [$K^{-1}$]')
    ax1.set_ylabel(yaxis)
    #plt.hold(True)
    
    names = glob.glob(prefix+'_'+'PEO12'+'*.txt')
    length = len(names)
    T = np.zeros(length)
    xi = np.zeros(length)
    xi_low = np.zeros(length)
    xi_high = np.zeros(length)
    for i,name in enumerate(names):
        pattern = '(\d{3})\D'
        T[i] = float(re.findall(pattern, name)[0])
        data = np.loadtxt(name)
        xi[i] = data[1]
      
        xi_low[i] = data[0]
        xi_high[i] = data[2]

    newInd = np.argsort(T)
    T = T[newInd]
    xi = xi[newInd]
    xi_high = xi_low[newInd]
    xi_low = xi_low[newInd]

    output = np.transpose(np.vstack((T,xi,xi_low,xi_high)))     
    np.savetxt('XivT_PEO12.txt',output)
        
    for i,x in enumerate(1/T):
        
        ax1.errorbar(x,xi[i],yerr=[[xi_low[i]],[xi_high[i]]],fmt='g.-',
                     label=prefix,capsize=2)

    ax1.plot(1/T,xi,'g.--',label=prefix)
    ax1.set_xlabel('1/T [$K^{-1}$]')
    ax1.set_ylabel(yaxis)
    
    plt.savefig(prefix+'.png',dpi=1000)
    plt.savefig(prefix+'.eps',format='eps')
    plt.show()
    
if __name__ == "__main__":
  main(sys.argv[1:])
