#CalcEjs.py
"""

    Class definition for ej calculation

"""

import MDAnalysis as mda
from MDAnalysis.analysis.base import AnalysisBase
import pandas as pd
import networkx as nx
import numpy as np
import sys

from .calcDipole import calcDipole

class CalcEjs(AnalysisBase):
    def __init__(self, g1, g2, name, **kwargs):
        self.name = name
        self.g1 = g1
        self.g2 = g2
        self.min_dist = kwargs[name]['distance_cutoff'][0]
        self.max_dist = kwargs[name]['distance_cutoff'][1]
        self.pre_allo = np.zeros([len(g1),len(g2)])
        atomgroup = g1.union(g2)
        self.universe = kwargs['universe']

        if 'XYZ' in kwargs[name].keys(): 
            self.XYZ = True
            self.xyz_name = kwargs[name]['XYZ']
        else: self.XYZ = False
        
        if 'angle_cutoff' in kwargs[name].keys(): 
            self.hbond = True
            self.angle_min = kwargs[name]['angle_cutoff'][0]
            self.angle_max = kwargs[name]['angle_cutoff'][1]
        else: 
            self.hbond = False
            self.angle_min = 0
            self.angle_max = 0
        
        if 'calc_dipole' in kwargs[name].keys(): self.calc_dipole = True
        else: self.calc_dipole = False

        super(CalcEjs, self).__init__(atomgroup.universe.trajectory)

    @staticmethod
    def remove_empty_nodes(dictionary, keys):
        for k in list(keys):
            if dictionary[k] == {}:
                del dictionary[k]
        return dictionary

    @staticmethod
    def write_xyz(self, xyz_lst, outfile):
        atoms_to_write = self.universe.select_atoms(f'resid {" ".join(map(str, xyz_lst))}')
        with open(outfile, mode='a') as f:
            f.write(f'{len(atoms_to_write)}\n')
            f.write(f'\tStep: {self.count}\n')
            for atom in atoms_to_write:
                f.write(f'{atom.element}\t')
                f.write(f'{atom.position[0]} {atom.position[1]} {atom.position[2]}')
                f.write('\n')

    @staticmethod
    def prep_adj(ind1, ind2):
        nodes = list(ind1) + list(ind2)
        edges = [{} for _ in range(len(nodes))]
        return dict(zip(nodes, edges))

    def _prepare(self):
        self.count = 0
        self.results.graphs = {}
        self.results.geometric = {'Distances':[], 'Angles':[]}
        if self.XYZ == True:
            with open(self.xyz_name, mode='w') as f:
                f.close()

    def _single_frame(self):
        if self.count % 100 == 0: print(self.count)
        xyz_lst = []

        distances = mda.lib.distances.distance_array(self.g1.positions, self.g2.positions,
                box=self.universe.dimensions, result=self.pre_allo)
        #Return the indices where distances meets network requirement
        ej_lst = np.column_stack(np.where((distances >= self.min_dist) & (distances <= self.max_dist)))
       
        #Define empty atomic adj matrix and empty molecular adj matrix
        vr_adj_dict = self.prep_adj(self.g1.indices, self.g2.indices)
        vr_mol_adj_dict = self.prep_adj(self.g1.resindices, self.g2.resindices)
        atom_ejs = self.prep_adj(self.g1.indices, self.g2.indices)
        filtered_ejs = self.prep_adj(self.g1.resindices, self.g2.resindices)
        atom_donor_adj = self.prep_adj(self.g1.indices, self.g2.indices)
        mol_donor_adj = self.prep_adj(self.g1.resindices, self.g2.resindices)

        #Convert the indices to resids
        #Can make adj matrix to add weights based on Hbond strength
        for ej in ej_lst:
            hydrogen = self.g1[ej[0]]
            acceptor = self.g2[ej[1]]
            resid1 = hydrogen.resindex
            resid2 = acceptor.resindex
            a1 = hydrogen.index
            a2 = acceptor.index
            if resid1 != resid2:
                if self.hbond == False:
                    #Add connection to adj matricies
                    #filtered = undirected weighted connections
                    try: filtered_ejs[resid1][resid2]["weight"] += 1
                    except KeyError: filtered_ejs[resid1][resid2] = {"weight": 1}
                    try: filtered_ejs[resid2][resid1]["weight"] += 1
                    except KeyError: filtered_ejs[resid2][resid1] = {"weight": 1}

                    try: atom_ejs[a1][a2]["weight"] = 1
                    except KeyError: atom_ejs[a1][a2] = {"weight": 1}
                    try: atom_ejs[a2][a1]["weight"] = 1
                    except KeyError: atom_ejs[a2][a1] = {"weight": 1}
                   
                    #atom_donor = directed weighted between H to acceptor
                    try: atom_donor_adj[a1][a2]["weight"] += 1
                    except KeyError: atom_donor_adj[a1][a2] = {"weight": 1}

                    #mol_donor = directed weighted between molecules
                    try: mol_donor_adj[resid1][resid2]["weight"] += 1
                    except KeyError: mol_donor_adj[resid1][resid2] = {"weight": 1}

                    #Add distance to histogram in returned result 
                    self.results.geometric['Distances'].append(distances[ej[0],ej[1]])
                    xyz_lst.append(resid1)
                    xyz_lst.append(resid2)
                
                elif self.hbond == True:   #Calculate three point V(r) and connectivity with angles
                    donor = hydrogen.bonds[0].partner(hydrogen)
                    angle = mda.lib.distances.calc_angles(donor.position,hydrogen.position,acceptor.position,
                        box=self.universe.dimensions) * 57.2957795
                    
                    if (angle <= self.angle_max) and (angle >= self.angle_min):
                        #Add connection to adj matricies
                        #filtered = undirected weighted connections
                        try: filtered_ejs[resid1][resid2]["weight"] += 1
                        except KeyError: filtered_ejs[resid1][resid2] = {"weight": 1}
                        try: filtered_ejs[resid2][resid1]["weight"] += 1
                        except KeyError: filtered_ejs[resid2][resid1] = {"weight": 1}

                        try: atom_ejs[a1][a2]["weight"] = 1
                        except KeyError: atom_ejs[a1][a2] = {"weight": 1}
                        try: atom_ejs[a2][a1]["weight"] = 1
                        except KeyError: atom_ejs[a2][a1] = {"weight": 1}

                        #atom_donor = directed weighted between H to acceptor
                        try: atom_donor_adj[a1][a2]["weight"] += 1
                        except KeyError: atom_donor_adj[a1][a2] = {"weight": 1}

                        #mol_donor = directed weighted between molecules
                        try: mol_donor_adj[resid1][resid2]["weight"] += 1
                        except KeyError: mol_donor_adj[resid1][resid2] = {"weight": 1} 
                        
                        #Add distance+angle to histogram in returned result 
                        self.results.geometric['Distances'].append(distances[ej[0],ej[1]])
                        self.results.geometric['Angles'].append(angle)
                        xyz_lst.append(hydrogen.resindex)
                        xyz_lst.append(acceptor.resindex)

                        #Calculate Dipole-Dipole potential
                        if self.calc_dipole == True:
                            Vr = calcDipole(hydrogen, acceptor, donor,
                                    self.universe.dimensions)
                            
                            #Add dip-dip potential as weight in atomic adj matrix
                            vr_adj_dict[hydrogen.index][acceptor.index] = {"weight": Vr}
                            vr_adj_dict[acceptor.index][hydrogen.index] = {"weight": -1*Vr}
                            
                            #Add dip-dip potential as weight in molecular adj matrix
                            try: vr_mol_adj_dict[hydrogen.resindex][acceptor.resindex]["weight"] += Vr
                            except KeyError: vr_mol_adj_dict[hydrogen.resindex][acceptor.resindex] = {"weight": Vr}
                            try: vr_mol_adj_dict[acceptor.resindex][hydrogen.resindex]["weight"] += -1*Vr
                            except KeyError: vr_mol_adj_dict[acceptor.resindex][hydrogen.resindex] = {"weight": -1*Vr}

        if self.XYZ == True:
            if xyz_lst != []: self.write_xyz(self, xyz_lst, self.xyz_name)

        #Keep nodes without edges and make nx graph
        graphs = {}
        #filtered_ejs = self.remove_empty_nodes(filtered_ejs, filtered_ejs.keys())
        graphs['simple_resind'] = filtered_ejs
        graphs['simple_atomic'] = atom_ejs
        graphs['dir_atomic'] = atom_donor_adj
        graphs['dir_resind'] = mol_donor_adj
        if self.hbond == True and self.calc_dipole == True:
            vr_adj_dict = self.remove_empty_nodes(vr_adj_dict, vr_adj_dict.keys())
            graphs['atomic_potential'] = vr_adj_dict
            vr_mol_adj_dict = self.remove_empty_nodes(vr_mol_adj_dict, vr_mol_adj_dict.keys())
            graphs['molecular_potential'] = vr_mol_adj_dict
        
        self.results.graphs[self.count] = graphs
        self.count += 1
    
    def _conclude(self):
        df = pd.DataFrame.from_dict(self.results.graphs,
                                    orient='index',
                                    columns=list(self.results.graphs[0].keys()),
                                   )
        df.to_csv(f'{self.name}.csv', index=False)

