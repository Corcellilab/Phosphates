#get_velocities.py
"""

    Get velocities of specified atoms

"""

import numpy as np
from .calc_velocities import calc_velocities

def get_velocities(self, dt):
    #iterate through vacf_X and {g1: , g2: }
    result = {}
    for k,v in vars(self).items():
        name = k.split('_')
        if name[0] == 'vacf':
            atoms = {}
            for group, selection in v.items():
                s = getattr(self, selection)
                atoms[group] = self.universe.select_atoms(s)
            g = atoms['group1']
            for group in atoms.values():
                g += group
            ####
            result = calc_velocities(self, self.universe, g, dt)
            setattr(self, k, result)

