import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

### if there are multiple correlation functions
#dummy = np.loadtxt('autocorrF0_0.txt')
#time = dummy[:,0] 

#load ncorr vals
#corrs = np.zeros((0,len(time)))
#for i in np.int64(np.linspace(0,10)):
#    corrs[i,:] = np.loadtxt('autocorrF0_'+str(i)+'.txt')[:,1]

#corr = np.mean(corrs,axis=0)

### if there is a single correlation fucntion
time = np.loadtxt('autocorrF0.txt')[:,0]
corr = np.loadtxt('autocorrF0.txt')[:,1]
corr_mean = np.mean(corr)
sst = np.sum((corr-corr_mean)**2)

def fitfunc(t,tau):
    return np.exp(-t/tau)

popt,pcov = curve_fit(fitfunc,time,corr,p0=[0.3])
corrfit1 = fitfunc(time,popt[0])
ssr1 = np.sum((corr-corrfit1)**2)


def fitfunc(t,tau,beta):
    return np.exp(-(t/tau)**beta)

popt,pcov = curve_fit(fitfunc,time,corr,p0=[0.3,1.1])
corrfit1_stretch = fitfunc(time,popt[0],popt[1])
ssr1_stretch = np.sum((corr-corrfit1_stretch)**2)

def fitfunc(t,a1,tau1,tau2):
    return a1*np.exp(-t/tau1)+(1-a1)*np.exp(-t/tau2)

popt,pcov = curve_fit(fitfunc,time,corr,p0=[0.3,0.1,0.1])
corrfit2 = fitfunc(time,popt[0],popt[1],popt[2])
ssr2 = np.sum((corr-corrfit2)**2)

def fitfunc(t,a1,tau1,tau2,beta1,beta2):
    return a1*np.exp(-(t/tau1)**beta1) + (1-a1)*np.exp(-(t/tau2)**beta2)

popt,pcov = curve_fit(fitfunc,time,corr,p0=[0.3,0.1,0.1,1.1,1.1])
corrfit2_stretch = fitfunc(time,popt[0],popt[1],popt[2],popt[3],popt[4])
ssr2_stretch = np.sum((corr-corrfit2_stretch)**2)

def fitfunc(t,a1,a2,tau1,tau2,tau3):
    return a1*np.exp(-t/tau1)+a2*np.exp(-t/tau2)+(1.0-a1-a2)*np.exp(-t/tau3)

popt,pcov = curve_fit(fitfunc,time,corr,p0=[0.3,0.3,0.1,0.1,0.1])
corrfit3 = fitfunc(time,popt[0],popt[1],popt[2],popt[3],popt[4])
ssr3 = np.sum((corr-corrfit3)**2)

def var_J(w):
    term_a1 = (popt[2]/(1+w**2*popt[2]**2))**2
    term_a2 = (popt[3]/(1+w**2*popt[3]**2))**2
    term_a3 = (popt[4]/(1+w**2*popt[4]**2))**2
    term_t1 = popt[0]**2
    term_t2 = popt[1]**2
    term_t3 = (1-popt[0]-popt[1])**2
    var_a1 = pcov[0,0]
    var_a2 = pcov[1,1]
    var_a3 = pcov[0,0] + pcov[1,1]
    var_t1 = pcov[2,2]
    var_t2 = pcov[3,3]
    var_t3 = pcov[4,4]
    return term_a1*var_a1 + term_a2*var_a2 + term_a3*var_a3 + term_t1*var_t1 + term_t2*var_t2 + term_t3*var_t3


plt.plot(time,corr,label='average')
plt.plot(time,corrfit1,label='fit exp')
plt.plot(time,corrfit2,label='fit biexp')
plt.plot(time,corrfit3,label='fit triexp')
plt.plot(time,corrfit1_stretch,label='fit stretch')
plt.plot(time,corrfit2_stretch,label='fit bi-stretch')
plt.xlabel('Time (ps)')
plt.ylabel('$C(t)/C(0)$')
plt.legend()
plt.savefig('corr.png')
plt.close()

plt.loglog(time,corr,label='average')
plt.loglog(time,corrfit1,label='fit exp')
plt.loglog(time,corrfit2,label='fit biexp')
plt.loglog(time,corrfit3,label='fit triexp')
plt.loglog(time,corrfit1_stretch,label='fit stretch')
plt.loglog(time,corrfit2_stretch,label='fit bi-stretch')
plt.xlabel('Time (ps)')
plt.ylabel('$C(t)/C(0)$')
plt.legend()
plt.savefig('corr_log.png')
plt.close()

plt.plot(time,corr,label='average')
plt.plot(time,corrfit1,label='fit exp')
plt.xlabel('Time (ps)')
plt.ylabel('$C(t)/C(0)$')
plt.legend()
plt.savefig('corr_exp.png')
plt.close()

plt.plot(time,corr,label='average')
plt.plot(time,corrfit2,label='fit biexp')
plt.xlabel('Time (ps)')
plt.ylabel('$C(t)/C(0)$')
plt.legend()
plt.savefig('corr_biexp.png')
plt.close()

plt.plot(time,corr,label='average')
plt.plot(time,corrfit3,label='fit triexp')
plt.xlabel('Time (ps)')
plt.ylabel('$C(t)/C(0)$')
plt.legend()
plt.savefig('corr_triexp.png')
plt.close()

plt.plot(time,corr,label='average')
plt.plot(time,corrfit1_stretch,label='fit stretch')
plt.xlabel('Time (ps)')
plt.ylabel('$C(t)/C(0)$')
plt.legend()
plt.savefig('corr_stretch.png')
plt.close()

plt.plot(time,corr,label='average')
plt.plot(time,corrfit2_stretch,label='fit bi-stretch')
plt.xlabel('Time (ps)')
plt.ylabel('$C(t)/C(0)$')
plt.legend()
plt.savefig('corr_bistretch.png')
plt.close()

plt.loglog(time,corr,label='average')
plt.loglog(time,corrfit1,label='fit exp')
plt.xlabel('Time (ps)')
plt.ylabel('$C(t)/C(0)$')
plt.legend()
plt.savefig('corr_log_exp.png')
plt.close()

plt.loglog(time,corr,label='average')
plt.loglog(time,corrfit2,label='fit biexp')
plt.xlabel('Time (ps)')
plt.ylabel('$C(t)/C(0)$')
plt.legend()
plt.savefig('corr_log_biexp.png')
plt.close()

plt.loglog(time,corr,label='average')
plt.loglog(time,corrfit3,label='fit triexp')
plt.xlabel('Time (ps)')
plt.ylabel('$C(t)/C(0)$')
plt.legend()
plt.savefig('corr_log_triexp.png')
plt.close()

plt.loglog(time,corr,label='average')
plt.loglog(time,corrfit1_stretch,label='fit stretch')
plt.xlabel('Time (ps)')
plt.ylabel('$C(t)/C(0)$')
plt.legend()
plt.savefig('corr_log_stretch.png')
plt.close()

plt.loglog(time,corr,label='average')
plt.loglog(time,corrfit2_stretch,label='fit bi-stretch')
plt.xlabel('Time (ps)')
plt.ylabel('$C(t)/C(0)$')
plt.legend()
plt.savefig('corr_log_bistretch.png')
plt.close()

a1 = popt[0]
a2 = popt[1]
a3 = 1.0-a1-a2
tau1 = popt[2]
tau2 = popt[3]
tau3 = popt[4]

def spectral(freq):
    return (a1*tau1)/(1+freq**2.0*tau1**2) + (a2*tau2)/(1+freq**2.0*tau2**2) + (a3*tau3)/(1+freq**2.0*tau3**2)

freq = np.linspace(-4,4,100)

omegaS = np.array([1.0/17.0,1.0/1.7,1.0/0.61])
omegaI = np.array([1.0/11000.0,1.0/1100.0,1.0/400.0])
B = np.array([0.342,3.35,9.2])

xi = 5*spectral(omegaS)/(3*spectral(omegaI)+7*spectral(omegaS))
k_sig = 5.0*spectral(omegaS)
print(xi)
print(k_sig)
#print(var_J(omegaI)**(0.5)/spectral(omegaI))
#print(var_J(omegaS)**(0.5)/spectral(omegaS))
var_xi = 125.0/9.0/spectral(omegaI)**2*var_J(omegaS) + 125.0/2401.0/spectral(omegaS)**2*var_J(omegaI)
xi_err = np.sqrt(var_xi)
#print(xi_err)

xi = [33.3,3.38,0.87]
xi_err = np.array([[33.3-32.7,3.38-3.22,0.87-0.82],[34.2-33.3,3.67-3.38,0.95-0.87]])
#print(xi_err.shape)

plt.plot(freq,spectral(freq))
plt.xlabel('$\omega/ps^{-1}$')
plt.ylabel('$J(\omega)/ps$')
plt.savefig('spectral.png')
plt.close()

xi_exp = [36.0,6.0,2.2]
xi_exp_err = [2.0,2.0,0.6]

plt.plot(B,xi,label='MD Simulations',color='black')
plt.plot(B,xi_exp,label='NMRD experiment',color='blue')
plt.plot(np.array([0,0.342]),np.array([50.0,33.3]),'k-.')
plt.plot(np.array([0,0.342]),np.array([50.0,36.0]),'b-.')
plt.errorbar(B,xi,xi_err,color='black',capsize=6.0)
plt.errorbar(B,xi_exp,xi_exp_err,color='blue',capsize=6.0)
plt.xlabel('$B/T$')
plt.ylabel('$\\xi$')
plt.ylim([0, 55.0])
plt.xlim([0, 9.3])
plt.legend()
plt.savefig('cfactor_vs_B.png')
plt.close()
