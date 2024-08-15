#adjacency.py
"""

   Setup adjacency matrices

"""

import numpy as np

def makeAdj(selected_atoms):
  atom_adj = np.zeros((len(selected_atoms), len(selected_atoms)), dtype=float)
  selected_resids = [i.resid for i in selected_atoms]
  selected_resids = list(set(selected_resids))
  resid_adj = np.zeros((len(selected_resids), len(selected_resids)), dtype=float)

  atom_adj_index = []
  resid_adj_index = []
 
  for atom in selected_atoms:
      if len(atom.bonds) == 0:
        atom_adj_index.append(atom.id)
        x1 = atom_adj_index.index(atom.id)

        resid_adj_index.append(atom.resid)
        x1 = resid_adj_index.index(atom.resid)
      
      for bond in atom.bonds:
      #for bond in selected_atoms.bonds:
        if bond[0].id not in atom_adj_index:
          atom_adj_index.append(bond[0].id)
        if bond[1].id not in atom_adj_index:
          atom_adj_index.append(bond[1].id)
        
        x1 = atom_adj_index.index(bond[0].id)  #
        x2 = atom_adj_index.index(bond[1].id)
        atom_adj[x1, x2] = 1
        atom_adj[x2, x1] = 1

        if bond[0].resid not in resid_adj_index:
          resid_adj_index.append(bond[0].resid)
        if bond[1].resid not in resid_adj_index:
          resid_adj_index.append(bond[1].resid)

  universe_adj = (atom_adj_index, atom_adj, resid_adj_index, resid_adj)
  return universe_adj

