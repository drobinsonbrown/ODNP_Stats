import numpy as np
import sys, os
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import glob

### model to fit the ACF to

def fitfunc(t,a1,a2,tau1,tau2,tau3):
    return a1*np.exp(-t/tau1)+a2*np.exp(-t/tau2)+(1.0-a1-a2)*np.exp(-t/tau3)

### spectral density function resulting from tri-exponential model
def spectral(freq,a1,a2,tau1,tau2,tau3):
    a3 = 1.0-a1-a2
    return (a1*tau1)/(1+freq**2.0*tau1**2) + (a2*tau2)/(1+freq**2.0*tau2**2) + (a3*tau3)/(1+freq**2.0*tau3**2)

### integral of fitting function
def fitint(a1,a2,tau1,tau2,tau3):
    a3 = 1.0 - a1 - a2
    return a1*tau1 + a2*tau2 + a3*tau3

### bootstrapping function
def bootsRelax(nsamp,nresamp,time,Corrs):
    index = 0
    meanCorr = np.zeros((nresamp,len(time)))
    meanCorrFit = np.zeros((nresamp,len(time)))
    tau = np.zeros(nresamp)
    while index < nresamp:
        # select nresamp number of average ACFs
        idx = np.random.choice(Corrs.shape[0],nsamp)
        meanCorr[index,:] = np.mean(Corrs[idx,:],axis=0)
        plt.plot(time,meanCorr[index,:])
        
        # perform a fit to the ACF
        pguess = [0.1,0.5,20.0,5.0,0.2]
        popt,pcov = curve_fit(fitfunc,time,meanCorr[index],p0 = pguess)
       
        meanCorrFit[index,:] = fitfunc(time,popt[0],popt[1],popt[2],popt[3],popt[4])
        # analytical fourier transform considering the fit parameters
        tau[index] = fitint(popt[0],popt[1],popt[2],popt[3],popt[4])
        index+=1
    plt.show()
    plt.close()
    return meanCorr, tau

### relaxation time bootstrapping function
def bootsODNP(nsamp,nresamp,time,Corrs):
    index = 0
    meanCorr = np.zeros((nresamp,len(time)))
    meanCorrFit = np.zeros((nresamp,len(time)))
    k_sigma = np.zeros(nresamp)
    xi = np.zeros(nresamp)
    while index < nresamp:
        # select nresamp number of average ACFs
        idx = np.random.choice(Corrs.shape[0],nsamp)
        meanCorr[index,:] = np.mean(Corrs[idx,:],axis=0)
        plt.plot(time,meanCorr[index,:])

        # perform a fit to the ACF
        pguess = [0.1,0.5,20.0,5.0,0.2]
        popt,pcov = curve_fit(fitfunc,time,meanCorr[index],p0 = pguess)
       
        meanCorrFit[index,:] = fitfunc(time,popt[0],popt[1],popt[2],popt[3],popt[4])
        # analytical integral considering the fit parameters
        larmorS = 1/17.0
        larmorI = 1/11000.0
        JS = spectral(larmorS,popt[0],popt[1],popt[2],popt[3],popt[4])
        JI = spectral(larmorI,popt[0],popt[1],popt[2],popt[3],popt[4])
        k_sigma[index] = 5*JS
        xi[index] = 5.0*JS/(3.0*JI+7.0*JS)
        index+=1
    plt.show()
    plt.close()
    return meanCorr, tau

def CI(nsamp,nresamp,y):
    y = np.sort(y)

    # define the CI indices for 95% CI
    upperind = int(np.ceil(nresamp*.975))-1
    lowerind = int(np.floor(nresamp*.125))-1

    # Find means of y
    yMean = (y[int(nresamp/2)-1]+y[int(nresamp/2)])/2

    # find CI of y
    yCI = np.array([y[lowerind],y[upperind]])

    print('The lower CI value, mean and upper CI value for a 95% confidence interval:\n')
    print(yCI[0],yMean,yCI[1])
    return y

def main(args):
    try:
        prefix = args[0]
    except IndexError:
        prefix = 'autocorrF0'

    # define bootstrapping procedure parameters
    nsamp = len(glob.glob(prefix+'_*.txt')) # sample size
    nresamp = 500 # number of resampling steps

    # BootStrap!!!
    np.random.seed(123456)
    index = 0

    # load ACF samples
    names = glob.glob(prefix+'*.txt')
    time = np.loadtxt(names[0])[:,0]
    Corrs = np.zeros((nsamp,len(time)))
    for ind,name in enumerate(names):
        Corrs[ind,:] = np.loadtxt(name)[:,1]

    [meanCorr,tau]=bootsRelax(nsamp,nresamp,time,Corrs)

    # save a representation of average ACF
    Corr=np.mean(meanCorr,axis=0)
    np.savetxt('avg'+prefix+'.txt',(time,Corr))

    tau = CI(nsamp,nresamp,tau)

    nbins = 10

    plt.hist(tau,bins=nbins)
    plt.xlabel('Mean $\tau_{trans}$')
    plt.ylabel('frequency')
    plt.title('Number of resamples: {}\n Number of data points per resample: {} \n'.format(nresamp,nsamp)) 
    plt.savefig('tau_stats.png')
    plt.close()

if __name__ == "__main__":
  main(sys.argv[1:])
