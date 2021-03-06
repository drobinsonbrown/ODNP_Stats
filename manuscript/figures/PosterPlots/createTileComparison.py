import numpy as np
import matplotlib.pyplot as plt
import glob
import sys,os
import re

# normalize the time axes
t_HB = np.loadtxt('avgHBCorr_Water_290.txt')[0,:]
t_HB = t_HB/t_HB[len(t_HB)-1]
water_HB = np.loadtxt('avgHBCorr_Water_290.txt')[1,:]
TEMPO_HB = np.loadtxt('avgHBCorr_4HTEMPO_290.txt')[1,:]
PEO12_HB = np.loadtxt('avgHBCorr_PEO12_290.txt')[1,:]

t_number = np.loadtxt('avgnumberCorr_Water_290.txt')[0,:]
t_number = t_number/t_number[len(t_number)-1]
water_number = np.loadtxt('avgnumberCorr_Water_290.txt')[1,:]
TEMPO_number = np.loadtxt('avgnumberCorr_4HTEMPO_290.txt')[1,:]
PEO12_number = np.loadtxt('avgnumberCorr_PEO12_290.txt')[1,:]

t_OACF = np.loadtxt('avgOrientCorr_P1_Water_290.txt')[0,:]
#print(t_OACF[len(t_OACF)-1])
t_OACF = t_OACF/t_OACF[len(t_OACF)-1]
#print(t_OACF.shape)
water_OACF = np.loadtxt('avgOrientCorr_P1_Water_290.txt')[1,:]
TEMPO_OACF = np.loadtxt('avgOrientCorr_P1_4HTEMPO_290.txt')[1,:]
PEO12_OACF = np.loadtxt('avgOrientCorr_P1_PEO12_290.txt')[1,:]


t_ODNP = np.loadtxt('avgautocorrF0_Water_290.txt')[0,:]
#print(t_ODNP.shape)
t_ODNP = t_ODNP/t_ODNP[len(t_ODNP)-1]
water_ODNP = np.loadtxt('avgautocorrF0_Water_290.txt')[1,:]
TEMPO_ODNP = np.loadtxt('avgautocorrF0_4HTEMPO_290.txt')[1,:]
PEO12_ODNP = np.loadtxt('avgautocorrF0_PEO12_290.txt')[1,:]
#print(t_ODNP.shape,PEO12_ODNP.shape,TEMPO_ODNP.shape)

# four-panel
fig, axs = plt.subplots(2, 2, sharex='col', sharey='row',
                        gridspec_kw={'hspace': 0, 'wspace': 0})
(ax1, ax2), (ax3, ax4) = axs

line_labels = ['PEO12','4-OH-TEMPO','Bulk Water']

ax1.plot(t_number, PEO12_number, 'tab:green', linestyle='dashdot', label='PEO12')
ax1.plot(t_number, TEMPO_number, 'tab:blue', linestyle='solid', label='4-OH-TEMPO')
ax1.plot(t_number, water_number, 'tab:red', linestyle='dashed', label='Bulk Water')
ax1.set_ylabel('$P_{survival}(t)$',fontsize=12)
#ax1.legend(frameon=False)

ax2.plot(t_OACF, PEO12_OACF, 'tab:green', linestyle='dashdot', label='PEO12')
ax2.plot(t_OACF, TEMPO_OACF, 'tab:blue', linestyle='solid', label='4-OH-TEMPO')
ax2.plot(t_OACF, water_OACF, 'tab:red', linestyle='dashed', label='Bulk Water')
ax2.yaxis.tick_right()
ax2.yaxis.set_label_position("right")
ax2.set_ylabel('$C_{OACF}(t)$',fontsize=12)
#ax2.legend(frameon=False)

ax3.plot(t_HB, PEO12_HB, 'tab:green', linestyle='dashdot', label='PEO12')
ax3.plot(t_HB, TEMPO_HB, 'tab:blue', linestyle='solid', label='4-OH-TEMPO')
ax3.plot(t_HB, water_HB, 'tab:red', linestyle='dashdot', label='Bulk Water')
ax3.set_ylabel('$P_{HB}(t)$',fontsize=12)
#ax3.legend(frameon=False)

ax4.plot(t_ODNP, PEO12_ODNP, 'tab:green', linestyle='dashdot', label='PEO12')
ax4.plot(t_ODNP, TEMPO_ODNP, 'tab:blue', linestyle='solid', label='4-OH-TEMPO')
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
plt.savefig('corrs.png',dpi=1000)
plt.savefig('corrs.eps',format='eps',dpi=1000)
plt.show()
