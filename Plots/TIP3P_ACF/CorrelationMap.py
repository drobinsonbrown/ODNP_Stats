import numpy as np
import matplotlib.pyplot as plt
import glob
import sys,os
import re

def main(args):    
    try:
        prefix = args[0]
    except IndexError:
        prefix = 'Water'

    tau = np.transpose(np.vstack((np.loadtxt('tauWater.txt'),
                                  np.loadtxt('tau4HTEMPO.txt'),
                                  np.loadtxt('tauPEO12.txt'))))

    names = ['Survival','P1','P2','HB']

    fig0, ax0 = plt.subplots()
    #fig1, ax1 = plt.subplots()
    #fig2, ax2 = plt.subplots()
    #fig3, ax3 = plt.subplots()
    count=0
    for j in range(tau.shape[1]):
        
        for i in range(tau.shape[0]):
       
            if i==0:
                for k in range(tau.shape[1]):
                    count=count+1
                    ax0.plot(tau[i,k],tau[i,j],'o')
            else:
                break
    plt.show()
    print(count)          
    
if __name__ == "__main__":
  main(sys.argv[1:])
