import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import glob
import sys,os
import re

matplotlib.rcParams.update({'figure.figsize':[8.0,6.0]})
matplotlib.rcParams.update({'font.size':18.0})
matplotlib.rcParams.update({'legend.fontsize':18.0})
matplotlib.rcParams.update({'xtick.labelsize':14.0})
matplotlib.rcParams.update({'ytick.labelsize':14.0})


kModelData = np.loadtxt('ksigma_model.txt')
xiModelData = np.loadtxt('ksi_model.txt')
invDModel = xiModelData[0,:]


xiModel = xiModelData[1,:]
kModel = kModelData[1,:]

kData = np.loadtxt('KvInvD.txt')
xiData = np.loadtxt('XivInvD.txt')

gamma = np.array([0,0.01,0.03,0.1,0.3,1.0,3.0,6.0,8.0,10.0,15.0,30.0,60.0,100.0])
invD = kData[:,0]*2.3

xi = xiData[:,1]*100
xi_low = xiData[:,2]*100
xi_high = xiData[:,3]*100

k = kData[:,1]
k_low = kData[:,2]
k_high = kData[:,3]

D = 1.0/invD*2.3

fig2,ax2 = plt.subplots()#figsize=(4.0,3.0))
ax2.plot(gamma,D)
ax2.set_xlabel(r'collision frequency [$ps^{-1}$]')
ax2.set_ylabel(r'$Dx10^{-5}[cm^{2}/s]$')
plt.savefig('DvG.png',dpi=1000)

plt.tight_layout()
plt.close()

fig3,ax3 = plt.subplots()#figsize=(4.0,3.0))

for i,x in enumerate(invD):
    print(i)
    print(x)
    ax3.errorbar(x,k[i]/k[0],
                 yerr=[[k_low[i]/k[0]],[k_high[i]/k[0]]],
                 color='tab:orange',
                 linestyle='dashed'
                 ,capsize=2)
    ax3.errorbar(x,xi[i]/xi[0],
                 yerr=[[xi_low[i]/xi[0]],[xi_high[i]/xi[0]]],fmt='b.-'
                 ,capsize=2)

ax3.set_xlabel(r'$D(298.15K)/D_{local}$')
ax3.plot(invD,k/k[0],
         color='tab:orange',linestyle='dashed',
         label=r'$k_{\sigma}/k_{\sigma,bulk}$')
ax3.plot(invD,xi/xi[0],'b.--',label=r'$\xi/\xi_{bulk}$')
ax3.plot(invDModel,kModel,color='tab:orange',
         label=r'$k^{FFHS}_{\sigma}/k^{FFHS}_{\sigma,bulk}$')
ax3.plot(invDModel,xiModel,'b',label=r'$\xi^{FFHS}/\xi^{FFHS}_{bulk}$')
ax3.set_xlim([min(invD),max(invD)])

box = ax3.get_position()
ax3.set_position([box.x0,box.y0 + box.height*0.2, box.width, box.height*0.9])

ax3.legend(loc='upper center',fancybox=True,
           bbox_to_anchor=(0.5,-0.12),
           shadow=True,ncol=2)
#ax3.set_ylabel(r'$Dx10^{-5}[cm^{2}/s]$')
#plt.legend(loc='upper center')
#plt.tight_layout()
plt.savefig('modelCompare.png',dpi=1000)


matplotlib.rcParams.update({'figure.figsize':[4.0,3.0]})
matplotlib.rcParams.update({'font.size':9.0})

fig, axs = plt.subplots(2,1,figsize=(4.0,6.0))

for i,x in enumerate(invD):
    print(i)
    print(x)
    axs[0].errorbar(x,k[i],yerr=[[k_low[i]],[k_high[i]]],fmt='b.-'
                 ,capsize=2)
    axs[1].errorbar(x,xi[i],yerr=[[xi_low[i]],[xi_high[i]]],fmt='b.-'
                 ,capsize=2)

axs[0].plot(invD,k,'b.--')
axs[1].plot(invD,xi,'b.--')

axs[0].set_ylabel(r'$k_{\sigma}$')
axs[1].set_ylabel(r'$\xi$')

#axs[0].set_xlabel('weight fraction of PEO',fontsize=9)
#axs[1].set_xlabel('weight fraction of PEO',fontsize=9)
axs[1].set_xlabel(r'$D(298.15K)/D$')

plt.tight_layout()

plt.savefig('gammaPlot.png',dpi=1000)
#    plt.savefig(prefix+'.eps',format='eps')
#    plt.show()

