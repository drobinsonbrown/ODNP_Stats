from pyspecdata import *
from matplotlib.offsetbox import TextArea,AnnotationBbox
Dw_default = 2.3e-9 # Brandon's
Dsl_default = 4.1e-10 # Brandon's
d_default = 0.47e-9 # hodges + bryant, surface label for bilayer at 298K
normalization_ratio = 2e26 # the hand-picked weighting between the rotational and the translational
tau_rot_default = 5e-9

obs(r"""To explain some comments above more concretely,
    let's imagine a somewhat unrealistically simplistic
    (but easy to picture) system where half of the water molecules
    obey the FFHS model, and half of them obey the rotational
    spectral density
    \begin{equation}
        J_{rot}(f,\tau_{rot})
        =
        \Re
        \left\{
        \frac{\tau}{1+i2\pi f \tau_{rot}}
        \right\}
    \end{equation}
    (which comes from assuming a correlation function of
    \linebreak$\exp\left( \frac{-t}{\tau_{rot}} \right)$, and matches with the first term in BPP theory, as it should)
    and then
    investigate the dependence of \klow and \ksigma on 
    both $D_{local}$ and $\tau_{rot}$.
I hand-pick a relative normalization constant to give a reasonable weighting between the two spectral density functions, so that for the vesicle $D_{local}$ and $d$ above.
I get a total spectral density function that looks like this (for $d=%0.2f$\Ang, $\Dw=%g \mbox{m}^2/s$, $\tau_{rot}=%0.2f\ns$:
    """%(d_default/1e-10,Dw_default,tau_rot_default))
d_free = 0.3e-9 # just what I need to give me what I want for tau
def tau(Dw = Dw_default, Dsl = Dsl_default, d = d_default,**extrakwargs):
    return d**2/(Dw+Dsl)
def J_ffhs(f,**kwargs):
    z = sqrt(1j*2*pi*f*tau(**kwargs))
    return c_prime(**kwargs)*real((1.+z/4.)/(1.+z+4./9.*z**2+1./9.*z**3))
def J_rot(f,tau_rot = tau_rot_default,**extrakwargs):
    return real(tau_rot/(1+2j*pi*f*tau_rot))
def J(f,rot_weighting = 0.5,**kwargs):
    J_ffhs_at_zero = J_ffhs(0,**kwargs)
    J_rot_at_zero = J_rot(0,**kwargs)
    return (1.0-rot_weighting)*J_ffhs(f,**kwargs)+rot_weighting*J_rot(f,**kwargs)*normalization_ratio
def ksigma(f,**kwargs):
    fH = 1.51671e-3*f
    return 6*J(f-fH,**kwargs)-J(f+fH,**kwargs)
def krho(f,**kwargs):
    fH = 1.51671e-3*f
    return 6*J(f-fH,**kwargs)+3*J(fH,**kwargs)+J(f+fH,**kwargs)
def klow(f,**kwargs):
    return 5./3.*krho(f,**kwargs)-7./3.*ksigma(f,**kwargs)
def c_prime(Dw = Dw_default, Dsl = Dsl_default, d = d_default,**extrakwargs):
    return 1.0/d/(Dw+Dsl) # don't include comments
fl = figlist_var(width = 0.99,showbox = False)
fl.next('J_vs_omega')
f = logspace(6,12,100)
J_normalization = J(0)
plot(f,J(f)/J_normalization,'k',
        label = r'$J_{total}$ d = %0.2f$\AA$'%(d_default*1e10),
        linewidth=4)
plot(f,0.5*normalization_ratio*J_rot(f)/J_normalization,'b:',
        linewidth=2,label = r'$J_{rot}$ d = %0.2f$\AA$'%(d_default*1e10))
plot(f,0.5*J_ffhs(f)/J_normalization,'r--',
        linewidth=2,label = r'$J_{FFHS}$ d = %0.2f$\AA$'%(d_default*1e10))
#{{{ explicitly show the points
plot(9.8e9,J(9.8e9)/J_normalization,'go')
annotate(r"$k_\sigma$",xy=(9.8e9,J(9.8e9)/J_normalization),
    xycoords = 'data',
    textcoords = 'offset points', # I think this means offset in points
    xytext = (60,20),
    arrowprops=dict(arrowstyle='fancy',
        fc=r_[0,1,0],
        ec="none",
        alpha = 0.5,
        connectionstyle="arc3,rad=0.3"))
axvline(x=9.8e9,color='g',alpha = 0.25,linewidth=2)#or set vspan with two parameters
plot(9.8e9*1.51671e-3,J(9.8e9*1.51671e-3)/J_normalization,'go')
annotate(r"$k_{low}$",xy=(9.8e9*1.51671e-3,J(9.8e9*1.51671e-3)/J_normalization),
    xycoords = 'data',
    textcoords = 'offset points', # I think this means offset in points
    xytext = (60,-20),
    arrowprops=dict(arrowstyle='fancy',
        fc=r_[0,1.0,0],
        ec="none",
        alpha = 0.5,
        connectionstyle="arc3,rad=0.3"))
axvline(x=9.8e9*1.51671e-3,color='g',alpha = 0.25,linewidth=2)#or set vspan with two parameters
#}}}
expand_y()
autolegend()
fl.text(r"""\o{Thus, in my ultra-simplified model,
    the \Dw of the FFHS $J_{FFHS}$ controls the speed of the fast-timescale waters,
    while the $\tau_{rot}$ of the rotational $J_{rot}$
    controls the speed of the slow-timescale waters.

So, let's start by varying the speed of the fast water
    by changing \Dw (keeping $d = %0.2f$\Ang, as appropriate for a surface label, and keeping $\tau_{rot}=%0.2f$\ns).}"""%(1e10*d_default,tau_rot_default/1e-9))
fl.setprops(autopad = False)
fl.next('vs_D')
#D = linspace(0.05*Dw_default,3*Dw_default,100)
#D = logspace(log10(0.05*Dw_default),log10(3*Dw_default),100)
D = r_[Dw_default:Dw_default/28.0:600j] # the values I use in the retardation factor plot, except that I start at a retardation factor of 2
ax1 = gca()
plot(Dw_default/D,ksigma(9.8e9,Dw=D,d=d_free)/ksigma(9.8e9,d=d_free),'r-',
        label = r'$k_\sigma / k_{\sigma,bulk}$',ax = ax1)
autopad_figure()
ax2 = twinx()
plot(Dw_default/D,klow(9.8e9,Dw=D,d=d_free)/klow(9.8e9,d=d_free),'b-',
        label = r'$k_{low}/k_{low,bulk}$',ax = ax2)
plot(Dw_default/D,J(0,Dw=D,d=d_free)/J(0,d = d_free),'k--',
        label = r'$J(0 MHz)/\left( J(0 MHz)|_{d,D\rightarrow bulk} \right)$',
        ax = ax2)
xibulk = ksigma(9.8e9,d=d_free)/klow(9.8e9,d=d_free)
plot(Dw_default/D,ksigma(9.8e9,Dw=D,d=d_free)/krho(9.8e9,Dw=D,d=d_free)/xibulk,'g-',
    label = r'$\xi / \xi_{free}$',ax = ax1)
ax1.set_xlabel(r'%g $m^2/s$ / $D_{local}$'%(Dw_default))
autolegend(ax = ax1)
autolegend(ax = ax2)
fl.setprops(autopad = True)
fl.text(r"""\o{We see that both \klow and \ksigma are sensitive to $D_{local}$,
    which determines the FFHS translation.
The value of \klow is changing almost entirely because the value of $c'$ is
    changing.
I show this by showing that \klow tracks $J(0)$
    fairly closely
    (as a side note, $J(0)$ should have a straightforward relationship
    with $T_2^{-1}$).
As a sidenote,
    I can also compare \klow and \ksigma to the spectral density values that they are supposed to approximate (i.e. due to $J(\omega_e\pm\omega_H)\approx J(\omega_e)$ represent),
    and I see that they match nearly identically.}""")
fl.text(r"\o{On the other hand, only \klow varies as we change $\tau_{rot}$, and as we expect, it peaks where $\tau_{rot}$ matches the $1/2 \pi 15 \MHz$=%0.2f ns. This plot emphasizes the \textbf{disadvantage} of $\xi$: \textit{it changes as we change the slow-motional timescale}. By contrast, the \textbf{advantage} of \ksigma is that it does not.}"%(1./2/pi/15e6/1e-9))
fl.next('vs_taurot')
#tau_varied = logspace(log10(100e-12),log10(100e-9),100)
tau_varied = linspace(100e-12,20e-9,100)
#print 'ksigma is ',lsafen(ksigma(9.8e9,tau_rot = tau_varied)/ksigma(9.8e9,d = d_free))
#print 'klow is ',lsafen(klow(9.8e9,tau_rot = tau_varied)/klow(9.8e9,d = d_free))
plot(tau_varied/1e-9,klow(9.8e9,tau_rot = tau_varied)/klow(9.8e9,d = d_free),label = r'$k_{low} / k_{low,bulk}$ d = %0.2f$\AA$'%(d_default*1e10))
plot(tau_varied/1e-9,ksigma(9.8e9,tau_rot = tau_varied)/krho(9.8e9,tau_rot = tau_varied)/xibulk,label = r'$\xi/\xi_{free}$ d = %0.2f$\AA$'%(d_default*1e10))
plot(tau_varied/1e-9,ksigma(9.8e9,tau_rot = tau_varied)/ksigma(9.8e9,d = d_free),label = r'$k_\sigma / k_{\sigma,bulk}$ d = %0.2f$\AA$'%(d_default*1e10))
axis('tight')
ax = gca()
saveaxes = ax.get_ylim() # so that the T2 doesn't blow out the plot
axis('tight')
expand_y()
xlabel(r'$\tau_{rot}$ / ns')
autolegend()
fl.show('ksigmaklowcomp_121007.pdf')
