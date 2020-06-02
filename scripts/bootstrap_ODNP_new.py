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
def bootsODNP(nsamp,nresamp,time,Corrs0,Corrs1,Corrs2):
    # define larmor frequencies necessary for computation
    larmorS = 1/17.1
    larmorI = 1/11000.0
    larmorDiff = larmorS-larmorI
    larmorSum = larmorS+larmorI
    
    index = 0
    meanCorr0 = np.zeros((nresamp,len(time)))
    meanCorrFit0 = np.zeros((nresamp,len(time)))
    meanCorr1 = np.zeros((nresamp,len(time)))
    meanCorrFit1 = np.zeros((nresamp,len(time)))
    meanCorr2 = np.zeros((nresamp,len(time)))
    meanCorrFit2 = np.zeros((nresamp,len(time)))
    k_sigma = np.zeros(nresamp)
    xi = np.zeros(nresamp)
    while index < nresamp:
        # select nresamp number of average ACFs
        idx = np.random.choice(Corrs0.shape[0],nsamp)
        meanCorr0[index,:] = np.mean(Corrs0[idx,:],axis=0)
        plt.plot(time,meanCorr0[index,:])

        # perform a fit to the ACF
        pguess = [0.3,0.4,20.0,10.0,1.0]
        popt,pcov = curve_fit(fitfunc,time,meanCorr0[index],p0 = pguess,
                              bounds=((0,0,0,0,0),(1.0,1.0,200.0,200.0,200.0)))
        meanCorrFit0[index,:] = fitfunc(time,popt[0],popt[1],popt[2],popt[3],popt[4])
    
        J0_diff = spectral(larmorDiff,popt[0],popt[1],popt[2],popt[3],popt[4])

        meanCorr1[index,:] = np.mean(Corrs1[idx,:],axis=0)
        plt.plot(time,meanCorr1[index,:])

        # perform a fit to the ACF
        pguess = [0.3,0.4,20.0,10.0,1.0]
        popt,pcov = curve_fit(fitfunc,time,meanCorr1[index],p0 = pguess,
                              bounds=((0,0,0,0,0),(1.0,1.0,200.0,200.0,200.0)))
        meanCorrFit1[index,:] = fitfunc(time,popt[0],popt[1],popt[2],popt[3],popt[4])
        J1_I = spectral(larmorI,popt[0],popt[1],popt[2],popt[3],popt[4])
        
        meanCorr2[index,:] = np.mean(Corrs2[idx,:],axis=0)
        plt.plot(time,meanCorr2[index,:])

        # perform a fit to the ACF
        pguess = [0.3,0.4,20.0,10.0,1.0]
        popt,pcov = curve_fit(fitfunc,time,meanCorr2[index],p0 = pguess,
                              bounds=((0,0,0,0,0),(1.0,1.0,200.0,200.0,200.0)))
        meanCorrFit2[index,:] = fitfunc(time,popt[0],popt[1],popt[2],popt[3],popt[4])

        J2_sum = spectral(larmorSum,popt[0],popt[1],popt[2],popt[3],popt[4])
        
        k_sigma[index] = 6.0*J2_sum - J0_diff
        xi[index] = k_sigma[index]/(J0_diff + 3.0*J1_I + 6.0*J2_sum)
        index+=1
    plt.close()
    return meanCorr0, meanCorr1, meanCorr2, k_sigma, xi
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
    nsamp = len(glob.glob(prefix+'*F0*'+'*.txt')) # sample size
    nresamp = 10000 # number of resampling steps

    # BootStrap!!!
    np.random.seed(151321213)
    index = 0

    # load ACF samples

    names0 = glob.glob(prefix+'*F0*'+'*.txt')
    names0.sort()
    names1 = glob.glob(prefix+'*F1*'+'*.txt')
    names1.sort()
    names2 = glob.glob(prefix+'*F2*'+'*.txt')
    names2.sort()
    
    time = np.loadtxt(names0[0])[:,0]
    Corrs0 = np.zeros((nsamp,len(time)))
    Corrs1 = np.zeros((nsamp,len(time)))
    Corrs2 = np.zeros((nsamp,len(time)))
   
    for ind,name in enumerate(names0):
        Corrs0[ind,:] = np.loadtxt(name)[:,1]
    for ind,name in enumerate(names1):
        Corrs1[ind,:] = np.loadtxt(name)[:,1]
    for ind,name in enumerate(names2):
        Corrs2[ind,:] = np.loadtxt(name)[:,1]
        
    [meanCorr0,meanCorr1,meanCorr2,k_sigma,xi] = bootsODNP(nsamp,nresamp,time,Corrs0,Corrs1,Corrs2)

    # save a representation of average ACF
    Corr0=np.mean(meanCorr0,axis=0)
    Corr1=np.mean(meanCorr1,axis=0)
    Corr2=np.mean(meanCorr2,axis=0)
    np.savetxt('avg'+prefix+'F0'+'.txt',(time,Corr0))
    np.savetxt('avg'+prefix+'F1'+'.txt',(time,Corr1))
    np.savetxt('avg'+prefix+'F2'+'.txt',(time,Corr2))

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
