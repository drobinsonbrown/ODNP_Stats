import polymerlib as pl
import pytraj as pt
import numpy as np
import matplotlib.pyplot as plt

topName = 'PEO12_SL_6000.pdb'
trajName = 'nve_production_output_wrapped18_short.dcd'
traj = pt.load(trajName,top=topName)
traj
BoxL = np.array(traj[0].box.values[:3])/10.0

index_c1 = traj.top.select('@C1')
index_c2 = traj.top.select('@C5')
#print(index_c1,index_c2)

#Pos1_c = np.array([])
#Pos2_c = np.array([])

#for i in range(len(traj)):
#    Pos1_c = np.append(Pos1_c,[ traj.xyz[i,index_c1,0], traj.xyz[i,index_c1,1], traj.xyz[i,index_c1,2] ],axis=0)
#    Pos2_c = np.append(Pos2_c,[ traj.xyz[i,index_c2,0], traj.xyz[i,index_c2,1] , traj.xyz[i,index_c2,2] ], axis=0)

Pos1_c = traj.xyz[:,index_c1,:]
Pos2_c = traj.xyz[:,index_c2,:]

Ree = pl.ree(Pos1_c,Pos2_c,BoxL)
print(Ree)
print(np.mean(Ree))

# plot P(r)
plt.hist(Ree,bins=25)
plt.show()
