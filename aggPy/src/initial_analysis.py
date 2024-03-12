import MDAnalysis as mda
from MDAnalysis.transformations.wrap import (unwrap, wrap)
from MDAnalysis.analysis.base import (AnalysisBase, AnalysisFromFunction, analysis_class)
import numpy as np

def init_analysis(self):
    print(f'Data file: {self.data_file} \n trj_file: {self.trj_file} \n md_form: {self.md_form}')
    
    u = mda.Universe(self.data_file, self.trj_file, format=self.md_form)
    
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
    print(f'{len(u.trajectory)} steps')

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

    topology = u.atoms.groupby(f'{self.mol_deliminator}s') 

    self.universe = u
    self.box = box
    self.hydrogens = hydrogens
    self.donors = donors
    self.acceptors = acceptors
    self.pre_allo = pre_allo
    self.topology = topology

#####
