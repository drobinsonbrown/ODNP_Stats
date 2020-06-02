import numpy as np
import matplotlib.pyplot as plt
import glob
import sys,os
import re

def main(args):    
    try:
        prefix = args[0]
    except IndexError:
        prefix = 'xi'
        
    yaxis = r'$\xi$'
         
    names = glob.glob(prefix+'*.txt')
    length = len(names)
    gamma = np.zeros(length)
    xi = np.zeros(length)
    xi_low = np.zeros(length)
    xi_high = np.zeros(length)
    for i,name in enumerate(names):
        pattern = '-?\d+\.?\d*'
        setnums = re.findall(pattern, name)
        gamma[i] = float(setnums[len(setnums)-1])
        data = np.loadtxt(name)
        xi[i] = data[1]
      
        xi_low[i] = data[0]
        xi_high[i] = data[2]

    newInd = np.argsort(gamma)
    gamma = gamma[newInd]
    xi = xi[newInd]
    #print(xi)
    xi_high = xi_low[newInd]
    xi_low = xi_low[newInd]

    fig, ax1 = plt.subplots()

    D =	np.array([1.94,1.94,1.86,1.842,1.841,1.7,1.4,
                  1.1,1.0,0.88,0.7,0.43,0.25,0.16])
    invD = 1.0/D    
    output = np.transpose(np.vstack((invD,xi,xi_low,xi_high)))
    print(output)
    
    np.savetxt('XivInvD.txt',output)
    for i,x in enumerate(gamma):
        
        ax1.errorbar(x,xi[i],yerr=[[xi_low[i]],[xi_high[i]]],fmt='b.-'
                     ,capsize=2)
    ax1.plot(gamma,xi,'b.--',label=prefix)
    ax1.set_xlabel('collision frequency [$ps^{-1}$]',fontsize=14)
    ax1.set_ylabel(yaxis,fontsize=14)
    #ax1.hold(True)

    plt.savefig(prefix+'.png',dpi=1000)
    plt.savefig(prefix+'.eps',format='eps')
    plt.show()
        
if __name__ == "__main__":
  main(sys.argv[1:])
