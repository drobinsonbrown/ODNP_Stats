import pytraj as pt
traj = pt.load('nve_production_output_wrapped11_short.dcd',top='PEO12_SL_6000.pdb')
traj = pt.autoimage(traj)
pt.save('nve_production_output_unwrapped11_short.dcd',traj)
