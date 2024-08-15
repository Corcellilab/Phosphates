#makeResids.py
"""

    Set molecule resids if not alreay

"""

import MDAnalysis as mda
import numpy as np

from aggPy.network2 import makeGraphDict
from aggPy.network2 import Dijkstra

def makeResids(u):      #universe object
    atom_ej = np.zeros([len(u.atoms)+1, len(u.atoms)+1])
    for atom in u.atoms:
        if len(atom.bonds) == 0:
            atom_ej[atom.id, atom.id] = 1
        else:
            for bond in atom.bonds:
                atom_ej[bond[0].id, bond[1].id] = 1
                atom_ej[bond[1].id, bond[0].id] = 1
    
    top = makeGraphDict(atom_ej)
    nodes = top.keys()
    path = []
    for source_node in nodes:
        path.append(Dijkstra(top, nodes, source_node))

    mols = []
    t = [list(p.keys()) for p in path]
    t = sorted(t, key=lambda x: len(x))
    t = [set(x) for x in t]
    for i in range(len(t)):
        a = t.pop(0)
        if not any([a.issubset(x) for x in t]):
            mols.append(a)
    print(f'Molecules: {len(mols)}') 
    u.atoms[0].residue.resid = 0

    resid_num = 0
    for mol in mols:
        resid_num += 1
        newres = u.add_Residue(segment=u.segments, resid=resid_num, resname=resid_num, resnum=resid_num)
        for atom in mol:
            u.atoms[atom-1].residue = newres

######
