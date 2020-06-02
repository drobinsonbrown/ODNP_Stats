import numpy as np
import matplotlib.pyplot as plt
import glob
import sys,os
import re
import matplotlib

matplotlib.rcParams.update({'figure.figsize':[8.0,6.0]})
matplotlib.rcParams.update({'font.size':18.0})
matplotlib.rcParams.update({'xtick.labelsize':14.0})
matplotlib.rcParams.update({'ytick.labelsize':14.0})

def main(args):    
    try:
        prefix = args[0]
        system = args[1]
    except IndexError:
        prefix = 'tau_autocorrF0'
        system = 'Water'
        
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
    
    T = []
    tau = []
    tau_low = []
    tau_high = []
    xi = []
    xi_low = []
    xi_high = []

    system = 'Water'
    names = glob.glob(prefix+'_'+system+'*.txt')
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

    Tnew = []
    names = glob.glob('xi_'+system+'*.txt')    
    for name in names:
        pattern = '(\d{3})\D'
        Tnew.append(float(re.findall(pattern, name)[0]))
        data = np.loadtxt(name)
        xi.append(data[1])
        xi_low.append(data[0])
        xi_high.append(data[2])
        
    #T = np.array(T)
    Tnew = np.array(Tnew)
    xi = np.array(xi)
    xi_high = np.array(xi_high)
    xi_low = np.array(xi_low)
    
    newInd = np.argsort(Tnew)
    Tnew = Tnew[newInd]
    xi = xi[newInd]
    xi_high = xi_low[newInd]
    xi_low = xi_low[newInd]
    fig, ax1 = plt.subplots()
    ax2 = ax1.twinx()

    for i,x in enumerate(T):
        
        ax1.errorbar(x,tau[i],yerr=[[tau_low[i]],[tau_high[i]]],fmt='b.-',
                     label=prefix,capsize=2)
        ax2.errorbar(x,xi[i],yerr=[[xi_low[i]],[xi_high[i]]],fmt='r*-',
                     label=prefix,capsize=2)
    ax1.plot(T,tau,'b.--',label=prefix)
    ax2.plot(T,xi,'r.--',label=prefix)
    ax1.set_xlabel('T [K]')#,fontsize=14)
    ax1.set_ylabel(yaxis)#,fontsize=14)
    ax2.set_ylabel(r'$\xi$')#,fontsize=14)

    fig.savefig(prefix+'_'+system+'.png')
    fig.savefig(prefix+'_'+system+'.eps',format='eps')
    plt.show()

    T = []
    tau = []
    tau_low = []
    tau_high = []
    xi = []
    xi_low = []
    xi_high = []

    system = 'PEO12'
    names = glob.glob(prefix+'_'+system+'*.txt')
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

    Tnew = []
    names = glob.glob('xi_'+system+'*.txt')    
    for name in names:
        pattern = '(\d{3})\D'
        Tnew.append(float(re.findall(pattern, name)[0]))
        data = np.loadtxt(name)
        xi.append(data[1])
        xi_low.append(data[0])
        xi_high.append(data[2])
        
    #T = np.array(T)
    Tnew = np.array(Tnew)
    xi = np.array(xi)
    xi_high = np.array(xi_high)
    xi_low = np.array(xi_low)
    
    newInd = np.argsort(Tnew)
    Tnew = Tnew[newInd]
    xi = xi[newInd]
    xi_high = xi_low[newInd]
    xi_low = xi_low[newInd]
    fig, ax1 = plt.subplots()
    ax2 = ax1.twinx()

    for i,x in enumerate(T):
        
        ax1.errorbar(x,tau[i],yerr=[[tau_low[i]],[tau_high[i]]],fmt='b.-',
                     label=prefix,capsize=2)
        ax2.errorbar(x,xi[i],yerr=[[xi_low[i]],[xi_high[i]]],fmt='r*-',
                     label=prefix,capsize=2)
    ax1.plot(T,tau,'b.--',label=prefix)
    ax2.plot(T,xi,'r.--',label=prefix)
    ax1.set_xlabel('T [K]')#,fontsize=14)
    ax1.set_ylabel(yaxis)#,fontsize=14)
    ax2.set_ylabel(r'$\xi$')#,fontsize=14)

    fig.savefig(prefix+'_'+system+'.png')
    fig.savefig(prefix+'_'+system+'.eps',format='eps')
    plt.show()

    

if __name__ == "__main__":
  main(sys.argv[1:])
