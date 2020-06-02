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
import waterorient_P1_new as waterorient_P1
import waterorient_P2_new as waterorient_P2
import hbonds_int_new as hbonds
import waterfrac_new as waterfrac
import odnp_dynamics as odnp


### (insert other dynamics)
trajName = 'nve_production_output_wrapped_short.dcd'
topName = 'PEO12_SL_6000.pdb'
hbonds.HBsurvival(topName,trajName)
#waterorient_P1.P1(topName,trajName)
#waterorient_P2.P2(topName,trajName)
#waterfrac.survival(topName,trajName)
#odnp.compute(topName,trajName)
