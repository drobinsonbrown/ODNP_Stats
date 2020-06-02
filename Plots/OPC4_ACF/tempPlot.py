import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import glob
import sys,os
import re

matplotlib.rcParams.update({'figure.figsize':[4.0,3.0]})
matplotlib.rcParams.update({'font.size':9.0})
matplotlib.rcParams.update({'axes.labelsize':9.0})
matplotlib.rcParams.update({'legend.fontsize':9.0})
matplotlib.rcParams.update({'xtick.labelsize':7.0})
matplotlib.rcParams.update({'ytick.labelsize':7.0})


tauData = np.loadtxt('TauvT_PEO12.txt')
kData = np.loadtxt('KvT_PEO12.txt')
xiData = np.loadtxt('XivT_PEO12.txt')

wfrac = tauData[:,0]

tau = tauData[:,1]
tau_low = tauData[:,2]
tau_high = tauData[:,3]

k = kData[:,1]
k_low = kData[:,2]
k_high = kData[:,3]

xi = xiData[:,1]*100
xi_low = xiData[:,2]*100
xi_high = xiData[:,3]*100

print(xi)

fig, axs = plt.subplots(3,1,figsize=(4.0,9.0))

for i,x in enumerate(wfrac):
    print(i)
    print(x)
    axs[0].errorbar(x,tau[i],yerr=[[tau_low[i]],[tau_high[i]]],fmt='g.-'
                 ,capsize=2)
    axs[1].errorbar(x,k[i],yerr=[[k_low[i]],[k_high[i]]],fmt='g.-'
                 ,capsize=2)
    axs[2].errorbar(x,xi[i],yerr=[[xi_low[i]],[xi_high[i]]],fmt='g.-'
                 ,capsize=2)

axs[0].plot(wfrac,tau,'g.--')
axs[1].plot(wfrac,k,'g.--')
axs[2].plot(wfrac,xi,'g.--')

axs[0].set_ylabel(r'$\tau_{Survival}$')
axs[1].set_ylabel(r'$k_{\sigma}$')
axs[2].set_ylabel(r'$\xi$')

#axs[0].set_xlabel('weight fraction of PEO',fontsize=9)
#axs[1].set_xlabel('weight fraction of PEO',fontsize=9)
axs[2].set_xlabel('weight fraction of PEO')

plt.savefig('concPlot.png',dpi=1000)
#    plt.savefig(prefix+'.eps',format='eps')
#    plt.show()


matplotlib.rcParams.update({'figure.figsize':[4.0,3.0]})
matplotlib.rcParams.update({'font.size':9.0})
matplotlib.rcParams.update({'axes.labelsize':9.0})
matplotlib.rcParams.update({'legend.fontsize':9.0})
matplotlib.rcParams.update({'xtick.labelsize':7.0})
matplotlib.rcParams.update({'ytick.labelsize':7.0})

fig, axs2 = plt.subplots(2,1,figsize=(4.0,6.0))

for i,x in enumerate(tau):
    print(i)
    print(x)
    axs2[0].errorbar(x,xi[i],
                     xerr=[[tau_low[i]],[tau_high[i]]],
                     yerr=[[xi_low[i]],[xi_high[i]]],
                     fmt='g.-'
                     ,capsize=2)
    axs2[1].errorbar(x,k[i],
                     xerr=[[tau_low[i]],[tau_high[i]]],
                     yerr=[[k_low[i]],[k_high[i]]],
                     fmt='g.-'
                     ,capsize=2)
    
zxi = np.polyfit(tau,xi,1)
pxi = np.poly1d(zxi)
zk = np.polyfit(tau,k,1)
pk = np.poly1d(zk)

axs2[0].plot(tau,pxi(tau),'g--')
axs2[1].plot(tau,pk(tau),'g--')

axs2[1].set_ylabel(r'$k_{\sigma} [s^{-1}]$')
axs2[0].set_ylabel(r'$\xi$')

#axs2[0].set_xlabel('weight fraction of PEO',fontsize=9)
#axs2[1].set_xlabel('weight fraction of PEO',fontsize=9)
axs2[1].set_xlabel(r'$\tau_{survival} [ps]$')

R_xi = 1.0 - np.sum(np.array([ (xi[i]-pxi(tau[i]))**2.0 for i in range(len(xi)) ]))/np.sum(np.array([ (xi[i]-np.mean(xi))**2.0 for i in range(len(xi)) ]))

R_k = 1.0 - np.sum(np.array([ (k[i]-pk(tau[i]))**2.0 for i in range(len(k)) ]))/np.sum(np.array([ (k[i]-np.mean(k))**2.0 for i in range(len(k)) ]))

print(R_xi,R_k)

plt.savefig('tauPlot.png',dpi=1000)
