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


def HBsurvival(topFile,trajFile):
  """Runs the dynamics code for a trajectory in chunks.
     Inputs:
       args - (optional) list of arguments, first is topology file, next is trajectory
     Outputs:
       Only files holding the autocorrelation functions.
  """
  
  #try:
  #  topFile = args[0]
  #except IndexError:
  #  topFile = 'EO59wTEMPO_24296.pdb'

  #try:
  #  trajFile = args[1]
  #except IndexError:
  #  trajFile = 'nve_production_output_wrapped0.dcd'

  #Load in the files of interest
  #I generally find parmed does a better job loading different topologies, then pytraj is ok with parmed objects (the above statement is not true a.t.m.)
  pi = 3.14159265358979323846
  
  #identify centroid and water nucleus locations
  top = md.load(topFile).topology
  OrInd = top.select('name Or') #np.append(top.select("name O4"),[top.select("name O1")])
  OInds = top.select("name O")
  nucInds = top.select("name HW") 

  don = []
  for i,ind in enumerate(OInds):
    don.append(ind)
    don.append(ind)
  don = np.array([don])
  donh = nucInds
  acc = OInds

  cutoff = 0.8
  query = OrInd

  dt = 0.01 # time interval in picosecons
  chunkSpace = 10
  chunkSize = 100
  nChunks = 1000
  number = np.zeros(chunkSize)

  timeOrigin = 0
  maxSize = 0

  for i,chunk in enumerate(md.iterload(trajFile,top=topFile,chunk=chunkSize)):
    maxSize= maxSize+len(chunk)

  j = 0
  
  while timeOrigin<=int(maxSize-chunkSize):
    for i,chunk in enumerate(md.iterload(trajFile,top=topFile,
                                         skip=timeOrigin,
                                         chunk=chunkSize)):
      print('the index of the for loop is {}'.format(i))
      if i!=0:
        break
      else:
        print(chunk)
        time = chunk.time
        timevals = (time-time[0])*dt
        timeOrigin = timeOrigin+chunkSpace
        print(timeOrigin)
        don = []
        for i,ind in enumerate(OInds):
          don.append(ind)
          don.append(ind)
        don = np.array([don])
        donh = nucInds
        acc = OInds
        boxL = chunk.unitcell_lengths
        BoxL = boxL[0,:]
        boxV = boxL[0,0]**3

        #get density of the box
        n_waters = chunk.n_residues-1

        BulkDens = n_waters/boxV

        init_traj = chunk[0]
        trajLength = len(chunk.xyz)

        thisnumber = np.zeros(trajLength)

        acc = md.compute_neighbors(init_traj,cutoff,query,
                                   haystack_indices=acc)[0]

        donh = []
        for i,ind in enumerate(acc):
          donh.append(ind+1)
          donh.append(ind+2)
        donh = np.array([donh])

        don = []
        for i,ind in enumerate(acc):
          don.append(ind)
          don.append(ind)
        don = np.array([don])
    
        accpos = init_traj.xyz[0,acc,:]
        donpos = init_traj.xyz[0,don,:]
        donhpos = init_traj.xyz[0,donh,:]
        init_hb = wl.generalhbonds(accpos,donpos,donhpos,BoxL,0.35,120.0)
        init_hb = np.array([init_hb])[0]
       
        thisnumber[0] = np.sum(init_hb)
        print('printing thisnumber variable {}'.format(thisnumber[0]))
    
        for t,frame in enumerate(chunk[1:]):
          t=t+1

          # calculate position of acceptor O, donor O and donor H
          accpos = frame.xyz[0,acc,:]
          donpos = frame.xyz[0,don,:]
          donhpos = frame.xyz[0,donh,:]

          # find the hydrogen bond matrix accposxdonhpos dimensional
          hb = wl.generalhbonds(accpos,donpos,donhpos,BoxL,0.35,120.0)
          hb = np.array([hb])[0]
          hb = [ [init_hb[i,j] if hb[i,j]==init_hb[i,j] 
                  else 0 for j in range(hb.shape[1])] 
                 for i in range(hb.shape[0])]
          hb = np.array([hb])[0,:,:]
      
          thisnumber[t] = np.sum(hb)
    
        number += thisnumber
        nChunks += 1
        np.savetxt('HBCorr.txt',timevals,number/nChunks)
        print('printing the number variable {}'.format(number))
    j = j+1
  number /= float(nChunks)
  number /= number[0]

  print("the number of chunks is {}".format(nChunks))
  print("The number of origins is {}".format(j))

  def fitfunc(t,tau):
    return np.exp(-t/tau)

  popt,pcov = curve_fit(fitfunc,timevals,number,p0=[0.05])
  number_fit = fitfunc(timevals,popt[0])

  print("{}".format(popt[0]))

  np.savetxt('HBCorr.txt',(timevals,number))
  np.savetxt('HBCorr_fit.txt',(timevals,number_fit))
  plt.plot(timevals,number)
  plt.xlabel('time (ps)')
  plt.ylabel('fraction of HBs')
  plt.savefig('HBCorr.png')
  plt.close()

  plt.plot(timevals,number,label='data')
  plt.plot(timevals,number_fit,label='fit')
  plt.xlabel('time (ps)')
  plt.ylabel('fraction of HBs')
  plt.legend()
  plt.savefig('HBCorr_fit.png')
