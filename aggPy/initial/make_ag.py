#make_ag.py
"""

    Make atom groups from MDAnalysis

"""

import MDAnalysis as mda

def make_ag(self):
    for k,v in vars(self).items():
        name = k.split('_')
        if name[0] == 'ag':
            setattr(self, k, self.universe.select_atoms(v))
            try:
                print(f'{name}: {set(getattr(self, k).names)}')
            except mda.exceptions.NoDataError:
                print(f'{name}: {set(getattr(self, k).types)}')
    
