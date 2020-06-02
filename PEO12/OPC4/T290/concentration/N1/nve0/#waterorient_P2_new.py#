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
from scipy.special import legendre
from scipy import integrate
import matplotlib
try:
  os.environ["DISPLAY"]
except KeyError:
  showPlots = False
  matplotlib.use('Agg')
import matplotlib.pyplot as plt


def P2(topFile,trajFile):
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
  OrInd = top.select('name Or')
  OInds = top.select("name O")
  nucInds = top.select("name HW") #hydrogen atoms

  dt = 0.1 # time interval in picosecons
  chunkSpace = 100
  chunkSize = 1000
  nChunks = 0
  Corr = np.zeros(chunkSize)

  cutoff = 0.8
  query = OrInd

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
        chunkTime = chunk.time
        timevals = (chunkTime-chunkTime[0])
        BoxL = chunk.unitcell_lengths

        init_traj = chunk[0]
        init_Pos = init_traj.xyz
        trajLength = len(chunk.xyz)

        OInds = md.compute_neighbors(init_traj,cutoff,query,
                                   haystack_indices=OInds)[0]
        nucInds1 = [OInd+1 for OInd in OInds]
        nucInds2 = [OInd+2 for OInd in OInds]
        nucInds = nucInds1+nucInds2
        nucInds.sort()

        init_Hvec = init_traj.xyz[0,nucInds,:]
        init_Ovec = init_traj.xyz[0,OInds,:]
    
        init_u = wl.watohvec(init_Ovec,init_Hvec,BoxL[0,:])
        dot_u = np.array([np.dot(init_u[i,:],init_u[i,:]) 
                      for i in range(len(init_u))])
        init_P2 = 0.5*(3.0*dot_u**2.0-1.0)
 
        thisCorr = np.zeros(trajLength)
        thisCorr[0] = np.sum(init_P2)
    
        for t,frame in enumerate(chunk[1:]):
          t=t+1

          Hvec = frame.xyz[0,nucInds,:]
          Ovec = frame.xyz[0,OInds,:]
          u = wl.watohvec(Ovec,Hvec,BoxL[0,:])
          #u = np.concatenate((u[0],u[1]))
          dot_u = np.array([np.dot(u[i,:],init_u[i,:]) 
                            for i in range(len(init_u))])
          P2 = 0.5*(3.0*dot_u**2.0-1.0)
          thisCorr[t] = np.sum(P2)
      
        Corr += thisCorr
        nChunks += 1
        tempCorr = Corr/float(nChunks)
        np.savetxt('OrientCorr_P2.txt',(timevals,tempCorr))
        plt.plot(timevals,tempCorr)
        plt.xlabel('time (ps)')
        plt.ylabel('Orientational ACF')
        plt.savefig('OrientCorr_P2.png')
        plt.close()

  print(nChunks)
  Corr /= float(nChunks)
  Corr /= np.abs(Corr[0])
  timevals = np.linspace(0,len(Corr)-1,len(Corr))*dt

  def fitfunc(t,tau):
    return np.exp(-t/tau)

  popt,pcov = curve_fit(fitfunc,timevals,Corr,p0=[0.1])
  Corr_fit = fitfunc(timevals,popt[0])

  np.savetxt('OrientCorr_P2.txt',(timevals,Corr))
  np.savetxt('OrientCorrfit_P2.txt',(timevals,Corr_fit))
  plt.plot(timevals,Corr)
  plt.xlabel('time (ps)')
  plt.ylabel('Orientational ACF')
  plt.savefig('OrientCorr_P2.png')
  plt.close()

  plt.plot(timevals,Corr,label='data')
  plt.plot(timevals,Corr_fit,label='fit')
  plt.xlabel('time (ps)')
  plt.ylabel('Orientational ACF')
  plt.legend()
  plt.savefig('OrientCorrfit_P2.png')
