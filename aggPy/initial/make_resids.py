#make_resids.py
"""

    Define new resids to deliminate molecules
    uses the MDA fragments to deliminate

    Input: MDAnalysis universe obj
    
    updates universe resids

"""

import MDAnalysis as mda

def make_resids(u):
    fragments = u.atoms.fragments
    molecule = 1

    for fragment in fragments:
        newres = u.add_Residue(segment=u.segments, resid=molecule, resname=molecule, resnum=molecule)
        for atom in fragment:
            atom.residue = newres
        molecule += 1

