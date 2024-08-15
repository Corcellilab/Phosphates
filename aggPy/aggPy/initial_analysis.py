import MDAnalysis as mda
from MDAnalysis.transformations.wrap import (unwrap, wrap)
from aggPy.makeResids import makeResids
import numpy as np
import sys

def init_analysis(self):
    print(f'Data file: {self.data_file} \n trj_file: {self.trj_file} \nmd_form: {self.md_form}\n')#atom_style={self.atom_style}')
    
    u = mda.Universe(self.data_file, self.trj_file, format=self.md_form)#, atom_style=self.atom_style)
    
    workflow = []
    box = u.dimensions
    try:
        if box == None:
            box = [float(x) for x in list(dims.split(','))]
            box = np.array(box)
            transform = mda.transformations.boxdimensions.set_dimensions(box)
            workflow.append(transform)
    except ValueError:
        pass

    print(f'Box Dimensions: {box}')
    ag = u.atoms
    workflow.append(mda.transformations.unwrap(ag))

    u.trajectory.add_transformations(*workflow)

    print(f'Hydrogen Bonding Parameters: {self.HbondParm}')
    print(f'Hydrogen atom type: {self.hydrogen_type}')
    print(f'Donors type: {self.donors_type}')
    print(f'Acceptors type: {self.acceptors_type}')

    hydrogens = u.select_atoms(self.hydrogen_type)
    donors = u.select_atoms(self.donors_type)
    acceptors = u.select_atoms(self.acceptors_type)

    if self.merge_a_d == 'y':
        acceptors = acceptors.union(donors)

    pre_allo = np.zeros([len(hydrogens),len(acceptors)])

    print(f'Total atoms: {len(u.atoms)}')
    print(f'{len(hydrogens)} hydrogens of type {set(hydrogens.types)} found')
    print(f'{len(acceptors)} acceptors type {set(acceptors.types)} found')
    print(f'{len(donors)} donors type {set(donors.types)} found')

    name_data = False
    if self.make_resids == 'y':
        makeResids(u)
        topology = u.atoms.groupby('resids')
        mols = {}
        for k,v in topology.items():
            try:
                names = [atom.name for atom in v]
                name = ''.join(names)
                try: mols[name] += 1
                except KeyError: mols[name] = 1
                name_data = True
            except mda.exceptions.NoDataError:
                break
    else:
        topology = u.atoms.groupby(f'{self.mol_delim}s')
        mols = {}
        for k,v in topology.items():
            try:
                names = [atom.name for atom in v]
                name = ''.join(names)
                try: mols[name] += 1
                except KeyError: mols[name] = 1
                name_data = True
            except mda.exceptions.NoDataError:
                break
    if name_data == True:
        print(f'Molecules found:\n{mols}')


    self.universe = u
    self.box = box
    self.hydrogens = hydrogens
    self.donors = donors
    self.acceptors = acceptors
    self.pre_allo = pre_allo
    self.topology = topology
    self.selected_atoms = self.universe.select_atoms(
            f'(same residue as {self.hydrogen_type}) and (same residue as {self.donors_type}) and (same residue as {self.acceptors_type})'
                ,sorted=True) 

#####
