import numpy as np
import matplotlib.pyplot as plt
import glob
import sys,os
import re

# normalize the time axes
t_HB = np.loadtxt('avgHBCorr_Water_290.txt')[0,:]
t_HB = t_HB/t_HB[len(t_HB)-1]
PEO12_HB_265 = np.loadtxt('avgHBCorr_PEO12_265.txt')[1,:]
PEO12_HB_275 = np.loadtxt('avgHBCorr_PEO12_275.txt')[1,:]
PEO12_HB_290 = np.loadtxt('avgHBCorr_PEO12_290.txt')[1,:]
PEO12_HB_310 = np.loadtxt('avgHBCorr_PEO12_310.txt')[1,:]
PEO12_HB_340 = np.loadtxt('avgHBCorr_PEO12_340.txt')[1,:]

t_number = np.loadtxt('avgnumberCorr_Water_290.txt')[0,:]
t_number = t_number/t_number[len(t_number)-1]
PEO12_number_265 = np.loadtxt('avgnumberCorr_PEO12_265.txt')[1,:]
PEO12_number_275 = np.loadtxt('avgnumberCorr_PEO12_275.txt')[1,:]
PEO12_number_290 = np.loadtxt('avgnumberCorr_PEO12_290.txt')[1,:]
PEO12_number_310 = np.loadtxt('avgnumberCorr_PEO12_310.txt')[1,:]
PEO12_number_340 = np.loadtxt('avgnumberCorr_PEO12_340.txt')[1,:]

t_OACF = np.loadtxt('avgOrientCorr_P1_Water_290.txt')[0,:]
#print(t_OACF[len(t_OACF)-1])
t_OACF = t_OACF/t_OACF[len(t_OACF)-1]
#print(t_OACF.shape)
PEO12_OACF_265 = np.loadtxt('avgOrientCorr_P1_PEO12_265.txt')[1,:]
PEO12_OACF_275 = np.loadtxt('avgOrientCorr_P1_PEO12_275.txt')[1,:]
PEO12_OACF_290 = np.loadtxt('avgOrientCorr_P1_PEO12_290.txt')[1,:]
PEO12_OACF_310 = np.loadtxt('avgOrientCorr_P1_PEO12_310.txt')[1,:]
PEO12_OACF_340 = np.loadtxt('avgOrientCorr_P1_PEO12_340.txt')[1,:]

t_ODNP = np.loadtxt('avgautocorrF0_Water_290.txt')[0,:]
#print(t_ODNP.shape)
t_ODNP = t_ODNP/t_ODNP[len(t_ODNP)-1]
PEO12_ODNP_265 = np.loadtxt('avgautocorrF0_PEO12_265.txt')[1,:]
PEO12_ODNP_275 = np.loadtxt('avgautocorrF0_PEO12_275.txt')[1,:]
PEO12_ODNP_290 = np.loadtxt('avgautocorrF0_PEO12_290.txt')[1,:]
PEO12_ODNP_310 = np.loadtxt('avgautocorrF0_PEO12_310.txt')[1,:]
PEO12_ODNP_340 = np.loadtxt('avgautocorrF0_PEO12_340.txt')[1,:]
#print(t_ODNP.shape,PEO12_ODNP.shape,TEMPO_ODNP.shape)

# four-panel
fig, axs = plt.subplots(2, 2, sharex='col', sharey='row',
                        gridspec_kw={'hspace': 0, 'wspace': 0})
(ax1, ax2), (ax3, ax4) = axs

line_labels = ['PEO12','4-OH-TEMPO','Bulk Water']

ax1.plot(t_number, PEO12_number_265, color=plt.cm.Blues(0.2))
ax1.plot(t_number, PEO12_number_275, color=plt.cm.Blues(0.4))
ax1.plot(t_number, PEO12_number_290, color=plt.cm.Blues(0.6))
ax1.plot(t_number, PEO12_number_310, color=plt.cm.Blues(0.8))
ax1.plot(t_number, PEO12_number_340, color=plt.cm.Blues(1.0))

ax1.set_ylabel('$P_{survival}(t)$',fontsize=12)
#ax1.legend(frameon=False)

ax2.plot(t_OACF, PEO12_OACF_265, color=plt.cm.Blues(0.2))
ax2.plot(t_OACF, PEO12_OACF_275, color=plt.cm.Blues(0.4))
ax2.plot(t_OACF, PEO12_OACF_290, color=plt.cm.Blues(0.6))
ax2.plot(t_OACF, PEO12_OACF_310, color=plt.cm.Blues(0.8))
ax2.plot(t_OACF, PEO12_OACF_340, color=plt.cm.Blues(1.0))

ax2.yaxis.tick_right()
ax2.yaxis.set_label_position("right")
ax2.set_ylabel('$C_{OACF}(t)$',fontsize=12)
#ax2.legend(frameon=False)

ax3.plot(t_HB, PEO12_HB_265, color=plt.cm.Blues(0.2))
ax3.plot(t_HB, PEO12_HB_275, color=plt.cm.Blues(0.4))
ax3.plot(t_HB, PEO12_HB_290, color=plt.cm.Blues(0.6))
ax3.plot(t_HB, PEO12_HB_310, color=plt.cm.Blues(0.8))
ax3.plot(t_HB, PEO12_HB_340, color=plt.cm.Blues(1.0))

ax3.set_ylabel('$P_{HB}(t)$',fontsize=12)
#ax3.legend(frameon=False)

ax4.plot(t_ODNP, PEO12_ODNP_265, color=plt.cm.Blues(0.2))
ax4.plot(t_ODNP, PEO12_ODNP_275, color=plt.cm.Blues(0.4))
ax4.plot(t_ODNP, PEO12_ODNP_290, color=plt.cm.Blues(0.6))
ax4.plot(t_ODNP, PEO12_ODNP_310, color=plt.cm.Blues(0.8))
ax4.plot(t_ODNP, PEO12_ODNP_340, color=plt.cm.Blues(1.0))

ax4.yaxis.tick_right()
ax4.yaxis.set_label_position("right")
ax4.set_ylabel('$C_{ODNP}(t)$',fontsize=12)
#ax4.legend(frameon=False)

label = ['(a)','(b)','(c)','(d)']
for i,ax in enumerate(axs.flat):
    ax.text(0.05,0.95,label[i],fontsize=12)
    ax.set_xlabel('$t$',fontsize=12)
    ax.set_xticklabels([])
    
#fig.legend([ax4],labels=line_labels,frameon=False,loc='upper center')
#plt.legend(frameon=False)
plt.savefig('corrsTemp.png',dpi=1000)
plt.savefig('corrsTemp.eps',format='eps',dpi=1000)
plt.show()
