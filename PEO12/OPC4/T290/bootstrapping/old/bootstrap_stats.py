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

### bootstrapping function
def boots(nsamp,nresamp,time,Corrs):
    index = 0
    meanCorr = np.zeros((nresamp,len(time)))
    meanCorrFit = np.zeros((nresamp,len(time)))
    while index < nresamp:
        ### select nresamp number of average ACFs
        idx = np.random.choice(Corrs.shape[0],nsamp)
        meanCorr[index,:] = np.mean(Corrs[idx,:],axis=0)
        plt.plot(time,meanCorr[index,:])

        ### perform a fit to the ACF
        pguess = [0.2,0.3,20.0,10.0,1.0]
        popt,pcov = curve_fit(fitfunc,time,meanCorr[index],p0 = pguess)
        meanCorrFit[index,:] = fitfunc(time,popt[0],popt[1],popt[2],popt[3],popt[4])

        ### analytical fourier transform considering the fit parameters
        larmorS = 1/17.0
        larmorI = 1/11000.0
        JS = spectral(larmorS,popt[0],popt[1],popt[2],popt[3],popt[4])
        JI = spectral(larmorI,popt[0],popt[1],popt[2],popt[3],popt[4])
        k_sigma[index] = 5*JS
        xi[index] = 5.0*JS/(3.0*JI+7.0*JS)
        index+=1

def main(args):
    try:
        prefix = args[0]
    except IndexError:
        prefix = 'autocorrF0'

    # define bootstrapping procedure parameters
    nsamp = len(glob.glob(prefix+'_*.txt')) # sample size
    nresamp = 1000 # number of resampling steps

    # BootStrap!!!

    np.random.seed(123456)
    index = 0

    # load ACF samples
    names = glob.glob(prefix+'*.txt') 
    time = np.loadtxt(names[0])[:,0]
    Corrs = np.zeros((nsamp,len(time)))
    for ind,name in enumerate(names):
        Corrs[ind,:] = np.loadtxt(name)[:,1]
    
    meanCorr = np.zeros((nresamp,len(time)))
    meanCorrFit = np.zeros((nresamp,len(time)))
    k_sigma = np.zeros(nresamp)
    xi = np.zeros(nresamp)
    while index < nresamp:
        ### select nresamp number of average ACFs
        idx = np.random.choice(Corrs.shape[0],nsamp)
        meanCorr[index,:] = np.mean(Corrs[idx,:],axis=0)
        plt.plot(time,meanCorr[index,:])

        ### perform a fit to the ACF
        pguess = [0.2,0.3,20.0,10.0,1.0]
        popt,pcov = curve_fit(fitfunc,time,meanCorr[index],p0 = pguess)
        meanCorrFit[index,:] = fitfunc(time,popt[0],popt[1],popt[2],popt[3],popt[4])

        ### analytical fourier transform considering the fit parameters
        larmorS = 1/17.0
        larmorI = 1/11000.0
        JS = spectral(larmorS,popt[0],popt[1],popt[2],popt[3],popt[4])
        JI = spectral(larmorI,popt[0],popt[1],popt[2],popt[3],popt[4])
        k_sigma[index] = 5*JS
        xi[index] = 5.0*JS/(3.0*JI+7.0*JS)
        index+=1
    plt.show()
    plt.close()
    Corr=np.mean(meanCorr,axis=0)
    np.savetxt('avg'+prefix+'.txt',(time,Corr))
    print(np.mean(k_sigma))

    k_sigma = np.sort(k_sigma)
    xi = np.sort(xi)

    # define the CI indices for 95% CI
    upperind = int(np.ceil(nresamp*.975))-1
    lowerind = int(np.floor(nresamp*.125))-1
    print(upperind,lowerind)

    # Find means for each
    k_sigmaMean = (k_sigma[int(nresamp/2)-1]+k_sigma[int(nresamp/2)])/2
    xiMean = (xi[int(nresamp/2)-1]+xi[int(nresamp/2)])/2*100

    # find CI for each
    k_sigmaCI = np.array([k_sigma[lowerind],k_sigma[upperind]])
    xiCI = 100*np.array([xi[lowerind],xi[upperind]])

    print('The lower CI value, mean and upper CI value for a 95% confidence interval:\n')
    print('k_sigma')
    print(k_sigmaCI[0],k_sigmaMean,k_sigmaCI[1])
    print('xi')
    print(xiCI[0],xiMean,xiCI[1])


    nbins = 10

    plt.hist(k_sigma,bins=nbins)
    plt.xlabel('Mean $k_{\sigma}$')
    plt.ylabel('frequency')
    plt.title('Number of resamples: {}\n Number of data points per resample: {} \n'.format(nresamp,nsamp)) 
    plt.savefig('k_sigma_stats.png')
    plt.close()

    plt.hist(xi,bins=nbins)
    plt.xlabel('Mean xi')
    plt.ylabel('frequency')
    plt.title('Number of resamples: {}\n Number of data points per resample: {} \n'.format(nresamp,nsamp)) 
    plt.savefig('xi_stats.png')
    plt.close()

if __name__ == "__main__":
  main(sys.argv[1:])
