import numpy as np
import sys, os
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import glob

### model to fit the ACF to

def fitfunc(t,m,b):
    return m*t+b

### bootstrapping function
def bootsTransport(nsamp,nresamp,time,rmsds):
    index = 0
    meanRMSD = np.zeros((nresamp,len(time)))
    meanRMSDFit = np.zeros((nresamp,len(time)))
    D = np.zeros(nresamp)
    while index < nresamp:
        # select nresamp number of average ACFs
        idx = np.random.choice(rmsds.shape[0],nsamp)
        meanRMSD[index,:] = np.mean(rmsds[idx,:]/6000.0,axis=0)
        plt.plot(time,meanRMSD[index,:])
        
        # perform a fit to the ACF
        pguess = [0.003,0.1]
        popt,pcov = curve_fit(fitfunc,time,meanRMSD[index],p0 = pguess)
       
        meanRMSDFit[index,:] = fitfunc(time,popt[0],popt[1])
        D[index]=popt[0]/6.0/4.0
        index+=1
    plt.show()
    plt.close()
    return meanRMSD, D

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
    print(yMean-yCI[0],yMean,yCI[1]-yMean)
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
    print(names)
    time = np.loadtxt(names[0])[0,:]#,100:1000]
    print(time)
    RMSDs = np.zeros((nsamp,len(time)))
    for ind,name in enumerate(names):
        RMSDs[ind,:] = np.loadtxt(name)[1,:]#100:1000]

    [meanRMSD,D]=bootsTransport(nsamp,nresamp,time,RMSDs)

    # save a representation of average ACF
    RMSD=np.mean(meanRMSD,axis=0)
    np.savetxt('avg'+prefix+'.txt',(time,RMSD))

    tau = CI(nsamp,nresamp,D)

    nbins = 10

    plt.hist(tau,bins=nbins)
    plt.xlabel('Mean $D$')
    plt.ylabel('frequency')
    plt.title('Number of resamples: {}\n Number of data points per resample: {} \n'.format(nresamp,nsamp)) 
    plt.savefig('D_stats.png')
    plt.close()

if __name__ == "__main__":
  main(sys.argv[1:])
