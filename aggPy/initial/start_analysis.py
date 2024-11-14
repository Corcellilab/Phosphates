#start_analysis.py
"""

    Start analysis objects

"""

import MDAnalysis as mda
from MDAnalysis.transformations.wrap import (unwrap, wrap)

from .make_resids import make_resids

def start_analysis(self):
    self.universe = mda.Universe(self.data_file, self.trj_file, format=self.md_form)#, atom_style=self.atom_style)
    self.dt = self.universe.trajectory.dt

    print(f'Timestep read: {self.universe.trajectory.dt} ps')
    print(f'Length of Traj: {len(self.universe.trajectory)} steps')
    
    workflow = []
    box = self.universe.dimensions
    try:
        if box == None:
            box = [float(x) for x in list(dims.split(','))]
            box = np.array(box)
            transform = mda.transformations.boxdimensions.set_dimensions(box)
            workflow.append(transform)
    except ValueError:
        pass

    print(f'Box Dimensions: {box}')
    #workflow.append(mda.transformations.unwrap(self.universe.select_atoms('all')))
    
    self.universe.trajectory.add_transformations(*workflow)

    if self.make_resids == "y": 
        make_resids(self.universe)

