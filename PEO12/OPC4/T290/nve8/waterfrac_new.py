import sys, os
import copy
import numpy as np
import waterlib as wl
import parmed as pmd
import mdtraj as md
from scipy.interpolate import UnivariateSpline as UVS
from scipy.optimize import curve_fit
from scipy.special import iv
from scipy.special import kv
from scipy.special import spherical_in
from scipy.special import spherical_kn
from scipy import integrate
import matplotlib
try:
  os.environ["DISPLAY"]
except KeyError:
  showPlots = False
  matplotlib.use('Agg')
import matplotlib.pyplot as plt


def survival(topFile,trajFile):
  """Runs the dynamics code for a trajectory in chunks.
     Inputs:
       args - (optional) list of arguments, first is topology file, next is trajectory
     Outputs:
       Only files holding the autocorrelation functions.
  """
  
  #Load in the files of interest
  #I generally find parmed does a better job loading different topologies, then pytraj is ok with parmed objects (the above statement is not true a.t.m.)
  pi = 3.14159265358979323846
  
  #identify centroid and water nucleus locations
  top = md.load(topFile).topology
  OrInd = top.select('name Or') #np.append(top.select("name O4"),[top.select("name O1")])
  OInds = top.select("name O")
  nucInds = top.select("name HW") 

  cutoff = 0.8
  query = OrInd

  #traj = traj[np.int64(len(Pos)/4.0):]
  dt = 0.01 # time interval in picosecons
  chunkSpace = 1000
  chunkSize = 1000
  nChunks = 0
  number = np.zeros(chunkSize)
  
  # define the time origin and total length of the trajectory                  
  timeOrigin=0
  maxSize=0
  for i,chunk in enumerate(md.iterload(trajFile,top=topFile,chunk=chunkSize)):
    maxSize = maxSize+len(chunk)
  
  while timeOrigin<=int(maxSize-chunkSize):

    for i,chunk in enumerate(md.iterload(trajFile,top=topFile,skip=timeOrigin,
                                       chunk=chunkSize)):
      print('the index of the for loop is {}'.format(i))
      if i!=0:
        break
      else:
        print(chunk)
        timeOrigin = timeOrigin+chunkSpace
        boxL = chunk.unitcell_lengths
        boxV = boxL[0,0]**3

        #get density of the box
        n_waters = chunk.n_residues-1

        BulkDens = n_waters/boxV

        init_traj = chunk[0]
        trajLength = len(chunk.xyz)
        init_neighbors = md.compute_neighbors(init_traj,cutoff,query,
                                              haystack_indices=OInds)[0]
        # indices of < cutoff water O's and TEMPO atoms
        restrictInds = np.concatenate((OrInd,init_neighbors)) 

        trajchunk = chunk.restrict_atoms(restrictInds)
        init_traj = chunk[0]
        init_neighbors = md.compute_neighbors(init_traj,cutoff,[0])[0]
        haystack = init_neighbors[len(OrInd):]

        init_neighbors = md.compute_neighbors(init_traj,cutoff,[0],
                                              haystack_indices=haystack)[0]
    
        init_number = len(init_neighbors)
        thisnumber = np.zeros(trajLength)
        thisnumber[0] = init_number

        for t,frame in enumerate(chunk[1:]):
          t=t+1
          neighbors = md.compute_neighbors(frame,cutoff,[0],
                                           haystack_indices=haystack)[0]
          thisnumber[t] = np.sum(np.array([ 1 if neighbors[i] in 
                                            init_neighbors else 0 
                                            for i in range(len(neighbors))]))
          init_neighbors = [init_neighbors[i] if init_neighbors[i] in 
                            neighbors else 0 for i in 
                            range(len(init_neighbors))]
        number += thisnumber
        nChunks += 1
  print(nChunks)
  number /= float(nChunks)
  number /= number[0]
  print(len(number))
  timevals =np.linspace(0,(len(number)-1))*dt
  print(len(timevals))
  print(timevals)
  #time = chunk.time
  #timevals = (time-time[0])*dt

  def fitfunc(t,a1,tau1,tau2):
    return a1*np.exp(-t/tau1)+(1-a1)*np.exp(-t/tau2)

  popt,pcov = curve_fit(fitfunc,timevals,number,p0=[0.1,0.1,0.1])
  number_fit = fitfunc(timevals,popt[0],popt[1],popt[2])

  np.savetxt('numberCorr.txt',(timevals,number))
  np.savetxt('numberCorr_fit.txt',(timevals,number_fit))
  plt.plot(timevals,number)
  plt.xlabel('time (ps)')
  plt.ylabel('fraction of waters')
  plt.savefig('numberCorr.png')
  plt.close()

  plt.plot(timevals,number,label='data')
  plt.plot(timevals,number_fit,label='fit')
  plt.xlabel('time (ps)')
  plt.ylabel('fraction of waters')
  plt.legend()
  plt.savefig('numberCorr_fit.png')
