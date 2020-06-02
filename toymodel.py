import matplotlib.pyplot as plt
import matplotlib.pylab as pylab
import numpy as np

D_local = 2.3e-9
D_sl = 4.1e-10
d = 4.7e-10
freq = 94.0e9
taurot = 5e-9
d_free = 3e-10
rot_weighting = 0.5 #percentage/100 of bound waters

t = np.linspace(100e-12,20e-9,100) #t_rot axis
taxis = [1,5,10,15,20,25] #set tick locations
tylim = [0.0,1.05] # lim for t_rot plot Y axis

D = np.r_[D_local*3:D_local/28:10000j] # d_local axis
daxis = [1, 5, 10, 15, 20, 25, 28] #set tick locations
dylim1 = [0.0,1.65] # ylim for k_sig and ksi Y axis
dylim2 = [0.8,4.1] # ylim for k_low Y axis

f = np.logspace(6,12,100) # frequency axis for spectral density plot
normalization_ratio = 2e26 # jf normalization ratio

#plot settings
plt_params = {'legend.fontsize': 'small',
          'figure.figsize': (12, 6),
         'axes.labelsize': 'small',
         'axes.titlesize':'small',
         'xtick.labelsize':'small',
         'ytick.labelsize':'small',
         'font.weight':'heavy'}
pylab.rcParams.update(plt_params)


def tau(Dw = D_local, Dsl = D_sl, d = d,**extrakwargs):
    return d**2/(Dw+Dsl)
    
def J_ffhs(f,**kwargs):
    z = np.sqrt(1j*2*np.pi*f*tau(**kwargs))
    return c_prime(**kwargs)*np.real((1.+z/4.)/(1.+z+4./9.*z**2+1./9.*z**3))
    
def J_rot(f,tau_rot = taurot,**extrakwargs):
    return np.real(tau_rot/(1+2j*np.pi*f*tau_rot))
    
def J(f,rot_weight = rot_weighting,**kwargs):
    J_ffhs_at_zero = J_ffhs(0,**kwargs)
    J_rot_at_zero = J_rot(0,**kwargs)
    return (1.0-rot_weight)*J_ffhs(f,**kwargs)+rot_weight*J_rot(f,**kwargs)*normalization_ratio
    
def ksigma(f,**kwargs):
    fH = 1.51671e-3*f
    return 6*J(f-fH,**kwargs)-J(f+fH,**kwargs)
    
def krho(f,**kwargs):
    fH = 1.51671e-3*f
    return 6*J(f-fH,**kwargs)+3*J(fH,**kwargs)+J(f+fH,**kwargs)
    
def klow(f,**kwargs):
    return 5./3.*krho(f,**kwargs)-7./3.*ksigma(f,**kwargs)
    
def c_prime(Dw = D_local, Dsl = D_sl, d = d,**extrakwargs):
    return 1.0/d/(Dw+Dsl)
    

xibulk = ksigma(freq,d=d_free)/klow(freq,d=d_free)

ksig_t = ksigma(freq,tau_rot = t)/ksigma(freq,d = d_free)
#krho_t = krho(freq,tau_rot = t)/krho(freq,d = d_free)
ksi_t = ksigma(freq,tau_rot = t)/krho(freq,tau_rot = t)/xibulk
klow_t = klow(freq,tau_rot = t)/klow(freq,d = d_free)

ksig_d = ksigma(freq,Dw=D,d=d_free)/ksigma(freq,d=d_free)
#krho_d = krho(freq,Dw=D,d=d_free)/krho(freq,d=d_free)
ksi_d = ksigma(freq,Dw=D,d=d_free)/krho(freq,Dw=D,d=d_free)/xibulk
klow_d = klow(freq,Dw=D,d=d_free)/klow(freq,d=d_free)


J_normalization = J(0)

Jffhs = (1.0-rot_weighting)*J_ffhs(f)/J_normalization
Jrot = rot_weighting*normalization_ratio*J_rot(f)/J_normalization
Jtot = J(f)/J_normalization

ksig_pt = J(freq)/J_normalization
klow_pt = J(freq*1.51671e-3)/J_normalization

## plotting

fig = plt.figure()
#t_plot = plt.subplot(131)
t_plot = plt.Axes(fig, [0.06, .1, .25, .80])
fig.add_axes(t_plot)
t_plot.plot(t/1e-9, ksig_t, 'r-', label = r'$k_\sigma / k_{\sigma,bulk}$')
t_plot.plot(t/1e-9, ksi_t, 'g-', label = r'$\xi / \xi_{bulk}$')
t_plot.plot(t/1e-9, klow_t, 'b-', label = r'$k_{low}/k_{low,bulk}$')
t_plot.set_xlabel(r'$\tau_{rot}$ $/$ $ns$')
t_plot.set_xticks(taxis)
t_plot.set_ylim(tylim)
t_plot.legend()

#d_plot = plt.subplot(132)
d_plot = plt.Axes(fig, [.39, .1, .25, .80])
fig.add_axes(d_plot)
d_plot.plot(D_local/D, ksig_d, 'r-', label = r'$k_\sigma / k_{\sigma,bulk}$')
d_plot.plot(D_local/D, ksi_d, 'g-', label = r'$\xi / \xi_{bulk}$')
d_plot2 = d_plot.twinx()
d_plot2.plot(D_local/D, klow_d, 'b-', label = r'$k_{low}/k_{low,bulk}$')
d_plot2.plot(D_local/D,J(0,Dw=D,d=d_free)/J(0,d = d_free),'b--', label = r'$J(0)/\left( J(0)|_{d,D\rightarrow bulk} \right)$')
d_plot.set_xlabel(r'$%g$$e^{-9}$ $ m^{2}/s $ $/$ $D_{local}$'%(D_local*1e9))
d_plot.set_xticks(daxis)
d_plot.set_ylim(dylim1)
d_plot2.set_ylim(dylim2)
d_plot2.tick_params(axis='y', colors='blue')
d_plot.legend(loc='upper right')
d_plot2.legend(loc='upper left')

#j_plot = plt.subplot(133)
j_plot = plt.Axes(fig, [.73, .1, .25, .80])
fig.add_axes(j_plot)
j_plot.plot(np.log(f), Jffhs, 'r:',linewidth=2,label = r'$J_{FFHS}$')
j_plot.plot(np.log(f), Jrot, 'b:',linewidth=2,label = r'$J_{rot}$')
j_plot.plot(np.log(f), Jtot, 'k--',linewidth=2,label = r'$J_{total}$')

j_plot.plot(np.log(freq),J(freq)/J_normalization,'ro')
j_plot.annotate(r"$k_\sigma$",xy=(np.log(freq),J(freq)/J_normalization),
    xycoords = 'data',
    textcoords = 'offset points',
    xytext = (15,20),
    color = 'r',
    arrowprops=dict(arrowstyle='fancy',
        fc=np.r_[1,0,0],
        ec="none",
        alpha = 0.5,
        connectionstyle="arc3,rad=0.2"))
j_plot.axvline(x=np.log(freq),color='r',alpha = 0.25,linewidth=2)

j_plot.plot(np.log(freq*1.51671e-3),J(freq*1.51671e-3)/J_normalization,'bo')
j_plot.annotate(r"$k_{low}$",xy=(np.log(freq*1.51671e-3),J(freq*1.51671e-3)/J_normalization),
    xycoords = 'data',
    textcoords = 'offset points',
    xytext = (15,-20),
    color = 'b',
    arrowprops=dict(arrowstyle='fancy',
        fc=np.r_[0,0,1],
        ec="none",
        alpha = 0.5,
        connectionstyle="arc3,rad=0.2"))
j_plot.axvline(x=np.log(freq*1.51671e-3),color='b',alpha = 0.25,linewidth=2)

j_plot.set_xlabel(r'$Frequency$ $(Hz)$')
j_plot.set_xticks([13.8155,16.1181,18.4207,20.7233,23.0259,25.3284,27.6310])
j_plot.set_xticklabels(['10$^{6}$','10$^{7}$','10$^{8}$','10$^{9}$','10$^{10}$','10$^{11}$','10$^{12}$'])
j_plot.legend()

plt.show()
