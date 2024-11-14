#calcDipole.py
"""
    Calculate dipole for graph
    edge weight
"""

from MDAnalysis.analysis.distances import distance_array
from MDAnalysis.lib.distances import calc_angles
from math import cos
import numpy as np

def calcDipole(h_atom, a_atom, d_atom, box_dims):
    #Fix dipoles to hydrogen origin
    i = distance_array(h_atom.position, d_atom.position, box=box_dims)[0][0]
    j = distance_array(h_atom.position, a_atom.position, box=box_dims)[0][0]
    ui = (h_atom.charge + d_atom.charge) / i
    uj = (h_atom.charge + a_atom.charge) / j
    r = distance_array(d_atom.position, a_atom.position, box=box_dims)[0][0]
    ti = calc_angles(h_atom.position, d_atom.position, a_atom.position,
            box=box_dims) * 57.2957795
    tj = calc_angles(h_atom.position, a_atom.position, d_atom.position,
            box=box_dims) * 57.2957795
    tij = calc_angles(d_atom.position, h_atom.position, a_atom.position,
            box=box_dims) * 57.2957795
    #print(i, j, ui, uj, r)
    #print(ti, tj, tij)

    dist_depend = 4*np.pi*(.008854)*(r**3)    #angstroms
    angular_depend = (cos(tij) - 3*cos(ti)*cos(tj))
    
    #Output energy in joules - 1e-12 scale,picoJoules
    Vr = -1 * ((ui*uj)/dist_depend) * angular_depend * 10**-10
    Vr = round(float('{:.25f}'.format(Vr)) * 1e12, 4)
    #print(Vr, ui*uj, dist_depend, angular_depend)
    
    return Vr

