import numpy as np
import sys, os
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from lmfit import Minimizer, Parameters, report_fit
import glob

### model to fit the ACF to

def fitfunc(t,a1,a2,tau1,tau2,tau3):
    return a1*np.exp(-t/tau1)+a2*np.exp(-t/tau2)+(1.0-a1-a2)*np.exp(-t/tau3)

### define objective function: returns array to be minimized
def fcn2min(params,t,data):
    """Model a tri-exponential and subtract data."""
    a1 = params['a1']
    a2 = params['a2']
    a3 = 1.0 - a1 - a2
    tau1 = params['tau1']
    tau2 = params['tau2']
    tau3 = params['tau3']
    model = a1*np.exp(-t/tau1) + a2*np.exp(-t/tau2) + a3*np.exp(-t/tau3)
    return model-data

### spectral density function resulting from tri-exponential model
def spectral(freq,a1,a2,tau1,tau2,tau3):
    a3 = 1.0-a1-a2
    return (a1*tau1)/(1+freq**2.0*tau1**2) + (a2*tau2)/(1+freq**2.0*tau2**2) + (a3*tau3)/(1+freq**2.0*tau3**2)

### bootstrapping function
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
        pguess = [0.3,0.4,20.0,10.0,1.0]
        popt,pcov = curve_fit(fitfunc,time,meanCorr[index],p0 = pguess,bounds=((0,0,0,0,0),(1.0,1.0,100.0,100.0,100)))
        meanCorrFit[index,:] = fitfunc(time,popt[0],popt[1],popt[2],popt[3],popt[4])
        # create a set of parameters
        params = Parameters()
        params.add('a1',value=0.3,min=0,max=1.0)
        params.add('a2',value=0.3,min=0,max=1.0)
        params.add('tau1',value=10.0,min=0)
        params.add('tau2',value=5.0,min=0)
        params.add('tau3',value=1.0,min=0)

        # do fit, here with the default leastsq algorithm
        minner = Minimizer(fcn2min,params,fcn_args=(time,meanCorr[index]))
        result = minner.minimize()

        # calculate final result
        final = meanCorr[index] + result.residual
        #report_fit(result)
        #print(popt)
        
        # analytical fourier transform considering the fit parameters
        larmorS = 1/1.70
        larmorI = 1/1100.0
        JS = spectral(larmorS,popt[0],popt[1],popt[2],popt[3],popt[4])
        JI = spectral(larmorI,popt[0],popt[1],popt[2],popt[3],popt[4])
        k_sigma[index] = 5*JS
        xi[index] = 5.0*JS/(3.0*JI+7.0*JS)
        index+=1
    plt.close()
    return meanCorr, k_sigma, xi
def CI(nsamp,nresamp,y):
    y = np.sort(y)

    # define the CI indices for 95% CI
    upperind = int(np.ceil(nresamp*.975))-1
    lowerind = int(np.floor(nresamp*.125))-1

    # Find means of y
    yMean = (y[int(nresamp/2)-1]+y[int(nresamp/2)])/2

    # find CI of y
    yCI = np.array([y[lowerind],y[upperind]])
    error = np.array([yMean-y[lowerind],yMean,y[upperind]-yMean])
    
    print('The lower CI value, mean and upper CI value for a 95% confidence interval:\n')
    print(yCI[0],yMean,yCI[1])
    return y,error

def main(args):
    try:
        prefix = args[0]
    except IndexError:
        prefix = 'autocorrF0'

    # define bootstrapping procedure parameters
    nsamp = len(glob.glob(prefix+'*.txt')) # sample size
    nresamp = 1000 # number of resampling steps

    # BootStrap!!!
    np.random.seed(151321213)
    index = 0

    # load ACF samples
    names = glob.glob(prefix+'*.txt') 
    time = np.loadtxt(names[0])[:,0]
    print(nsamp)
    Corrs = np.zeros((nsamp,len(time)))
    for ind,name in enumerate(names):
        Corrs[ind,:] = np.loadtxt(name)[:,1]

    [meanCorr,k_sigma,xi]=bootsODNP(nsamp,nresamp,time,Corrs)

    # save a representation of average ACF
    Corr=np.mean(meanCorr,axis=0)
    np.savetxt('avg'+prefix+'.txt',(time,Corr))

    [k_sigma,error_k_sigma] = CI(nsamp,nresamp,k_sigma)

    np.savetxt('k_sigma.txt',error_k_sigma)

    [xi,error_xi] = CI(nsamp,nresamp,xi)

    np.savetxt('xi.txt',error_xi)
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
