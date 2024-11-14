#calc_velocities.py
"""
    Assumes coordinates in angstroms
    Use minimum image convention to correct for pbc
    Return normalized vectors

    return dictionary with {atomic index: (vel_x, vel_y, vel_z)} 

"""

from MDAnalysis.analysis.distances import dist
import numpy as np
import sys

def calc_velocities(obj, u, atoms, dt):
    avg_vel = list()
    dims = u.dimensions
    Lx = dims[0] * 10**-10
    Ly = dims[1] * 10**-10
    Lz = dims[2] * 10**-10
    result_array = np.zeros((len(atoms), len(atoms)))
    
    if obj.run_length == 'n':
        start = obj.start_frame
        end = obj.end_frame+1
    else:
        start = 0
        end = len(u.trajectory)
   
    avg_vels = []
    len_atoms_25 = int(len(atoms) / 10)
    count = 0

    print(f'Calculating Velocities for {len(atoms)} atoms')
    u.trajectory[start]
    cords0 = atoms.positions * 10**-10
    for ts in range(start+1, end):
        if ts % 1000 == 0: print(f'Calculating step: {ts}')
        u.trajectory[ts]
        cords1 = atoms.positions * 10**-10
       
        diff = cords1 - cords0

        x = diff[:,0]
        y = diff[:,1]
        z = diff[:,2]

        img_x = np.abs((x/Lx).round())
        img_y = np.abs((y/Ly).round())
        img_z = np.abs((z/Lz).round())
       
        vel_x = (x - (Lx * img_x)) / dt
        vel_y = (y - (Ly * img_y)) / dt
        vel_z = (z - (Lz * img_z)) / dt
      
        #vels = np.zeros((len(vel_x), 3), dtype=float)
        #vels[:, 0] = vel_x  
        #vels[:, 1] = vel_y
        #vels[:, 2] = vel_z

        #row_norms = np.linalg.norm(vels, axis=1, keepdims=True)
        #velocities = vels / row_norms

        #avg_vel_x = np.average(velocities[:, 0])
        #avg_vel_y = np.average(velocities[:, 1])
        #avg_vel_z = np.average(velocities[:, 2])

        avg_vel_x = np.average(vel_x)
        avg_vel_y = np.average(vel_y)
        avg_vel_z = np.average(vel_z)

        avg_vels.append((avg_vel_x, avg_vel_y, avg_vel_z))
        
        cords0 = cords1

    print('Done')

    return avg_vels

