import sys, os
import copy
import numpy as np
import waterlib as wl
import parmed as pmd
import pytraj as pt
from scipy.interpolate import UnivariateSpline as UVS
from scipy.optimize import curve_fit
import matplotlib
try:
  os.environ["DISPLAY"]
except KeyError:
  showPlots = False
  matplotlib.use('Agg')
import matplotlib.pyplot as plt
#This defines a functions for calculating dipolar coupling factors from MD simulation similar to ODNP
#An example script that loops through a trajectory and implements the functions here is also included
#This is all based on the work of Sezer and coworkers, mainly the 2009 paper


def Ffuncs(distvecs):
  """Given vectors of distances (x,y,z) from the radical to nuclei, 
     calculates all Fi functions (see Sezer, 2009) that contributes to dipolar
     coupling. Computing all of them takes about the same time as computing 
     just one, so just return all of them.
     Also note that how you compute the distance is now up to you...
     Could be minimum image, or other convention.
     Inputs:
       distvecs - Nx3 array of x, y, and z distances from radical for all nuclei
     Outputs:
       fzerofunc - the value of the F0 function (as complex objects)
                   for every position in pos relative to its minimum 
                   image distance to refpos
       fonefunc - the values of the F1 function
       ftwofunc - the values of the F2 function
  """
  #Need squared distances in each dimension
  sqpos = distvecs**2

  #Compute radial squared distance from the radical
  rsq = np.sum(sqpos, axis=1)

  #And distance to the fifth power as well
  rfive = rsq*rsq*np.sqrt(rsq)

  #Compute the F0 function
  fzerofunc = np.zeros(distvecs.shape[0], dtype=complex)
  fzerofunc[:].real = np.sqrt(3.0/2.0) * (rsq - 3.0*sqpos[:,2]) / rfive

  #Compute the F1 function
  fonefunc = np.zeros(distvecs.shape[0], dtype=complex)
  fonefunc[:].real = 3.0 * distvecs[:,2] * distvecs[:,0] / rfive
  fonefunc[:].imag = 3.0 * distvecs[:,2] * distvecs[:,1] / rfive

  #Compute the F2 function
  ftwofunc = np.zeros(distvecs.shape[0], dtype=complex)
  ftwofunc[:].real = -(3.0/2.0) * (sqpos[:,0] - sqpos[:,1]) / rfive
  ftwofunc[:].imag = 3.0 * distvecs[:,0] * distvecs[:,1] / rfive
  
  #Return the function values, noting that the real part is in the first column, 
  #and the imaginary part is in the second column
  return fzerofunc, fonefunc, ftwofunc


def timeCorrFunc(trajchunk, radNames='@Or,Nr', nucNames='@HW'):
  """Given a pytraj trajectory object and atom names to treat as the O and N in the
     radical, as well as the names of water hydrogens, computes the time
     autocorrelation functions for all of the F functions from Sezer, 2009.
    
     Inputs:
       trajchunk - pytraj trajectory object
       radNames - (default '@Or,Nr') string specifying atom names of radical atoms
       nucNames - (default '@HW') string specifying atom names of water oxygens
     Outputs:
       autocorrFzero - time autocorrelation for F0 function (will have real and 
                       imaginary components in general, but F0 is only real)
       autocorrFone - time autocorrelation for F1 function
       autocorrFtwo - time autocorrelation for F2 function
  """
  #Get the indices of atoms with the desired names
  radInds = trajchunk.top.select(radNames)
  nucInds = trajchunk.top.select(nucNames)

  #We will use all of the trajectory frames
  #So just take the first frame and compute the starting positions
  #and F functions to treat as time zero
  pos = np.array(trajchunk[0].xyz)
  
  #Define single radical position as average of nitrogen and oxygen
  radPos = np.average(pos[radInds], axis=0)

  #And get the nuclei (hydrogen) positions
  nucPos = pos[nucInds]
 
  #Define box dimensions
  boxDims = np.array(trajchunk[0].box.values[:3])

  #Get minimum image distances in first frame, which we will add to
  nucPos = wl.reimage(nucPos, radPos, boxDims) - radPos

  #Get the F functions based on distances of nuclei relative to radical
  initFzero, initFone, initFtwo = Ffuncs(nucPos)

  #Set up arrays to store average values of the time correlation functions
  autocorrFzero = np.zeros(len(trajchunk), dtype=complex)
  autocorrFone = np.zeros(len(trajchunk), dtype=complex)
  autocorrFtwo = np.zeros(len(trajchunk), dtype=complex)

  #Also need to store a reference to the previous timestep nuclear-radical distances
  #This allows us to have distances larger than half the box length, preventing
  #periodic effects due to use of the minimum image convention
  nucPosOld = copy.deepcopy(nucPos)

  #And loop over the trajectory frames
  for t, frame in enumerate(trajchunk):

    thispos = np.array(frame.xyz)
    thisbox = np.array(frame.box.values[:3])
    thisrad = np.average(thispos[radInds], axis=0)
    thisnuc = thispos[nucInds]
    thisnuc = wl.reimage(thisnuc, thisrad, thisbox) - thisrad
    
    #Before computing F functions, need to ADD to nuclear distances based on shift since last frame 
    #Can accomplish by reimaging each nuclear position around its old position (relative to the radical)
    #Then taking difference gives amount (and direction) to add
    #Assumes no nuclei could possibly (physically) move half the box length (relative to the radical) between frames
    thisdists = np.zeros(thisnuc.shape)
    for k, nuc in enumerate(thisnuc):
      thisdists[k] = wl.reimage(np.array([nuc]), nucPosOld[k], thisbox)[0] - nucPosOld[k]
    nucPos += thisdists
    nucPosOld = thisnuc

    thisFzero, thisFone, thisFtwo = Ffuncs(nucPos)

    #Will average over all hydrogens to get the autocorrelation function average
    autocorrFzero[t] = np.average(initFzero * np.conj(thisFzero))
    autocorrFone[t] = np.average(initFone * np.conj(thisFone))
    autocorrFtwo[t] = np.average(initFtwo * np.conj(thisFtwo))
    
  #Return the autocorrelation functions
  return autocorrFzero, autocorrFone, autocorrFtwo


def compute(topFile,trajFile):
  """Runs the dynamics code for a trajectory in chunks.
     Inputs:
       args - (optional) list of arguments, first is topology file, next is trajectory
     Outputs:
       Only files holding the autocorrelation functions.
  """

  #Load in the files of interest
  #I generally find parmed does a better job loading different topologies, then pytraj is ok with parmed objects
  top = pmd.load_file(topFile)
  top = pt.load_parmed(top, traj=False)
  traj = pt.iterload(trajFile, top)
  #traj = traj[10000:,:,:] # discard the first nanosecond of NVE production
  #Define how we will loop over the trajectory in chunks
  timestep = 0.1
  chunkSize = 1000 #100 ps if saving frames every 0.1 ps
  chunkSpace = 100 #10 ps should be spaced far enough apart for uncorrelated configurations
  nChunks = 0 #Keep track of number of chunks to do averaging right
  timeVals = np.arange(0.0, (chunkSize+1)*timestep, timestep)
  #Set up arrays to store correlation functions
  corrFzero = np.zeros(chunkSize, dtype=complex)
  corrFone = np.zeros(chunkSize, dtype=complex)
  corrFtwo = np.zeros(chunkSize, dtype=complex)

  #Loop over the trajectory in chunks to look at dynamics
  for i in range(0, len(traj), chunkSpace):

    if i+chunkSize > len(traj):
      break

    thisZero, thisOne, thisTwo = timeCorrFunc(traj[i:i+chunkSize], radNames='@Or,Nr', nucNames='@HW')

    corrFzero += thisZero
    corrFone += thisOne
    corrFtwo += thisTwo

    nChunks += 1

  #Finish averaging and save raw correlation function data
  corrFzero /= float(nChunks)
  corrFone /= float(nChunks)
  corrFtwo /= float(nChunks)
  corrFzero = corrFzero/max(corrFzero)
  corrFone = corrFone/max(corrFone)
  corrFtwo = corrFtwo/max(corrFtwo) 
  
  # Define the frequency arguments used to determine k_sigma and
  # the coupling factor
  larmorS = 1.0/17.0 #From Sezer, 2009 at field strength of 0.342 Tesla (change as needed) 
  larmorI = 1.0/11000.0 #In ps - looks like may need to really collect correlations for a long time...
                        #It will depend on the field strength
                        #Could maybe get more bang for buck by also going backwards in time
                        #To accurately fit and evaluate FFT of correlation function, will need to go longer
                        #than the inverse of the smallest larmour frequency
  freqDiff = larmorI - larmorS
  freqAdd = larmorI + larmorS
  #magGyroI = 1.0 #Need to look up what these are, or just reporting rates divided by these
  #magGyroS = 1.0
  #magPermeability = 1.2566370614E-06 #in N/A^2
  #planckConstant = 6.62607004E-34  #in m^2kg/s
  #deltaFac = magPermeability * planckConstant * magGyroI * magGyroS / (4.0*np.pi)
  deltaFac = 1.0
  freqRange = np.fft.fftfreq(chunkSize, d=timestep)
  
  # Carry out correlation function fit a la Sezer's paper where tau1_inv, tau2_inv and tau3_inv are the inverses of the characteristic times in Sezer's tri-exponential fit 
  def fitfunc(t,a1,a2,tau1,tau2,tau3):
    return a1*np.exp(-t/tau1) + a2*np.exp(-t/tau2) + (1-a1-a2)*np.exp(-t/tau3)
  times = timeVals[:chunkSize]

# Plot correlation functions data
  plt.plot(times,corrFzero.real,'k-',label='C0')
  plt.plot(times,corrFone.real,'r-',label='C1')
  plt.plot(times,corrFtwo.real,'b-',label='C2')
  plt.xlabel('time [ps]')
  plt.ylabel('C(t)/C(0)')
  plt.legend()
  plt.savefig('acf.png')
  plt.close()



  np.savetxt('autocorrF0.txt', np.stack([timeVals[:chunkSize],corrFzero.real,corrFzero.imag],axis=1), header='Time (ps)      real(F_0)      im(F_0)',fmt="%.3e")
  np.savetxt('autocorrF1.txt', np.stack([timeVals[:chunkSize],corrFone.real,corrFone.imag],axis=1), header='Time (ps)      real(F_1)      im(F_1)',fmt="%.3e")
  np.savetxt('autocorrF2.txt', np.stack([timeVals[:chunkSize],corrFtwo.real,corrFtwo.imag],axis=1), header='Time (ps)      real(F_2)      im(F_2)',fmt="%.3e")




  #Next we want to actually compute the relaxation rates of interest
  #To do this, need to Fourier Transform the correlation functions (and only take real part)
  freqVals = np.fft.fftfreq(chunkSize, d=timestep)
  fftFzero = np.fft.fft(corrFzero).real
  fftFone = np.fft.fft(corrFone).real
  fftFtwo = np.fft.fft(corrFtwo).real
  
  #Save the FFTs
  np.savetxt('fft_F0.txt', np.transpose(np.vstack((freqVals,fftFzero))), header='Frequency (1.0/ps)      real(fft(F_0))',fmt="%.3e")
  np.savetxt('fft_F1.txt', np.transpose(np.vstack((freqVals,fftFone))), header='Frequency (1.0/ps)      real(fft(F_1))',fmt="%.3e")
  np.savetxt('fft_F2.txt', np.transpose(np.vstack((freqVals,fftFtwo))), header='Frequency (1.0/ps)      real(fft(F_2))',fmt="%.3e")
 
  #Interpolate fft's at frequencies defined above
  sortind = np.argsort(freqVals)
  freqVals = freqVals[sortind]
  fftFzero = fftFzero[sortind]
  fftFone = fftFone[sortind]
  fftFtwo = fftFtwo[sortind]
  fftFzeroInterp = UVS(freqVals,fftFzero)
  fftFoneInterp = UVS(freqVals,fftFone)
  fftFtwoInterp = UVS(freqVals,fftFtwo)

  #Now find frequencies from FFT closest to Larmor and FFT values there
  diffInd = np.argmin(abs(freqRange - freqDiff))
  addInd = np.argmin(abs(freqRange - freqAdd))
  larmorInd = np.argmin(abs(freqRange - larmorI))
  larmorSInd = np.argmin(abs(freqRange - larmorS))
  
  # Plot raw, fitted, and splined fft data
  plt.plot(freqVals,fftFzero,'ko',label='fft data')
  plt.plot( np.linspace(min(freqVals),max(freqVals),100),fftFzeroInterp(np.linspace(min(freqVals),max(freqVals),100)),'b-',label='splined fft data')
  


# compute rates and coupling factor using our method
  krho = ((deltaFac**2)/12.0) * (fftFzeroInterp(freqDiff) + 3.0*fftFone[larmorInd] + 6.0*fftFtwoInterp(freqAdd))
  ksig = ((deltaFac**2)/12.0) * (6.0*fftFtwoInterp(freqAdd)-fftFzeroInterp(freqDiff))
  ksi = ksig/krho

# compute coupling factor using Sezer method (note: added a factor of 2.0 because I think it's more mathematically rigorous, must check by looking at sigma and rho values)
  larmorS = np.array([1.0/17.0,1.0/1.7,1.0/0.61])
  larmorI = np.array([1.0/11000.0,1.0/1100.0,1.0/400.0])
