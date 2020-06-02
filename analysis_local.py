from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from sys import stdout
from parmed import gromacs
from parmed import unit as u
#from dcdreporter import DCDReporter
import parmed as pmd
import time
from openmmtools import *
import mdtraj as md
import mdtraj.reporters
import numpy as np
import pytraj as pt
#import odnp_dynamics as odnp
import odnp_dynamics_local as odnp


### (insert other dynamics)
trajName = 'nve.nc'
topName = 'membraneSolvated.pdb'
#hbonds.HBsurvival(topName,trajName)
#waterorient_P1.P1(topName,trajName)
#waterorient_P2.P2(topName,trajName)
#waterfrac.survival(topName,trajName)
odnp.compute(topName,trajName)
