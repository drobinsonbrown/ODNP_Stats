
# coding: utf-8

# 

get_ipython().magic(u'load_ext pyspecdata.ipy')
rcParams['image.cmap'] = 'jet'


# 

Dw_default = 2.3e-9 # diffusivity of water
Dsl_default = 4.1e-10 # Brandon's
d_default = 0.47e-9 # hodges + bryant, surface label for bilayer at 298K
d_free = 0.3e-9 # just what I need to give me what I want for tau
d_buried = 0.8e-9 # this is made up
tau_rot_default = 5e-9
default_normalization_ratio = normalization_ratio = 2e26 # the hand-picked weighting between the rotational and the translational


# 

def tau(Dw = Dw_default, Dsl = Dsl_default, d = d_default,tau_rot = 'junk',normalization_ratio = 'junk'):
    return d**2/(Dw+Dsl)
def J_ffhs(f,**kwargs):
    z = sqrt(1j*2*pi*f*tau(**kwargs))
    c = c_prime(**kwargs)
    retval = ((1.+z/4.)/(1.+z+4./9.*z**2+1./9.*z**3)).real * c
    return retval
def ksigma(f,**kwargs):
    fH = 1.51671e-3*f
    return 6*J(f-fH,**kwargs)-J(f+fH,**kwargs)
def krho(f,**kwargs):
    fH = 1.51671e-3*f
    return 6*J(f-fH,**kwargs)+3*J(fH,**kwargs)+J(f+fH,**kwargs)
def klow(f,**kwargs):
    return 5./3.*krho(f,**kwargs)-7./3.*ksigma(f,**kwargs)
def c_prime(Dw = Dw_default, Dsl = Dsl_default, d = d_default,tau_rot = 'junk',normalization_ratio = 'junk'):
    # this is bennati's k', which is multiplied by 7J(omega_s,tau_D)+3J(omega_I,tau_D) to get the diffusion
    #kprime = 32000.*pi/405.(mu_0/4./pi)**2*N_A*C*(gammabar_H*g_e*mu_B)**2*S*(S+1)/(d*(D_M+D_L))
    # remember that hbar*omega/B_0 = g*mu_B = gammabar_e*hbar, so 
    # also, eliminate C, so that this gives a relaxivity, rather than a relaxation rate
    gammabar_e = gammabar_H/1.51671e-3
    S = 0.5
    kprime = 32000.*pi/405
    kprime *= (mu_0/4./pi)**2
    kprime *= N_A
    kprime *= (hbar*2*pi*gammabar_H*2*pi*gammabar_e)**2 # here, I am converting from my usual circular units to angular units
    kprime *= S*(S+1)
    kprime /= d*(Dsl+Dw)
    #print 'for: <d> ={}\Ang <D>={}'.format(d.mean()/1e-10,Dsl+Dw.mean()),
    #print 'kprime is {}--{}'.format(kprime.min(),kprime.max())
    return kprime
def J_rot(f,tau_rot = tau_rot_default,**extrakwargs):
    return (tau_rot/(1+2j*pi*f*tau_rot)).real
def J(f,rot_weighting = 0.5,normalization_ratio = default_normalization_ratio,**kwargs):
    if rot_weighting != 0:
        J_ffhs_at_zero = J_ffhs(0,**kwargs)
        J_rot_at_zero = J_rot(0,**kwargs)
        return (1.0-rot_weighting)*J_ffhs(f,**kwargs)+rot_weighting*J_rot(f,**kwargs)*normalization_ratio
    else:
        #print "yes, rot weighting is 0:"
        return J_ffhs(f,**kwargs)


# 

f = logspace(6,12)


# # Normalization
# 
# The spectral density at zero frequency is the integral of the correlation function -- so it should be order of the product of the correlation time and the square of the interaction strength.
# 
# We should rationalize this later, but for now, let's set these to similar orders of magnitude.
# 
# Here, I show that this is the case over a reasonable range of values

# 

fl=figlist_var()
f = 9.8e9
#f = logspace(5,12,100)
#D = Dw_default*logspace(0,-3,100)
D_array = Dw_default*linspace(1,1e-3,100)
d_name = '$D/D_{bulk}$'
D = nddata(D_array,[d_name]).labels(d_name,D_array/Dw_default)
t_norm = 140.
J_ffhs(0.,Dw=D)/t_norm


# 

tau_name = r'$\tau_{rot}$'
tau_rot = linspace(0.5e-9,12e-9,100)
tau_rot = nddata(tau_rot,[tau_name]).labels(tau_name, tau_rot).set_units(tau_name,'s')


# the following corresponds to interactions of equal strength

# 

r_norm = 1e-9
J_rot(0.,tau_rot=tau_rot)/r_norm


# 

f = logspace(4,10,100)
f = nddata(f,[r'$\nu$']).labels(r'$\nu$',f)
(J_ffhs(f)/t_norm+J_rot(f)/r_norm).name(r'vs. frequency')


# ## Now put them together

# 

esr = 9.8e9
nmr = esr * gammabar_H / gammabar_e


# 

print "nmr is",nmr/1e6


# 

(J_ffhs(esr,Dw=D)/t_norm+J_rot(esr,tau_rot=tau_rot)/r_norm).name(r'magnitude of $k_{\sigma}$')


# 

(J_ffhs(nmr,Dw=D)/t_norm+J_rot(nmr,tau_rot=tau_rot)/r_norm).name(r'magnitude of $k_{low}$')


# 



