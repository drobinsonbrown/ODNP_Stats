import pytraj as pt
traj = pt.load('nve_production_output_wrapped17_short.dcd',top='PEO12_SL_6000.pdb')
traj = pt.autoimage(traj)
pt.save('nve_production_output_unwrapped17_short.dcd',traj)
