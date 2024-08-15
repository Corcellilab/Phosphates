def atom_mapping(atoms):
    per_table = {1.0078: 'H', 4.003: 'He', 6.94: 'Li', 9.012: 'Be', 10.81: 'B', 
            12.01: 'C', 14.01: 'N', 16.0: 'O', 19.0: 'F', 20.18: 'Ne', 22.93: 'Na', 
            24.31: 'Mg', 26.98: 'Al', 28.09: 'Si', 30.97: 'P', 32.06: 'S', 35.45: 'Cl', 
            39.95: 'Ar'}
    atom_map = {}
    atom_types = atoms.groupby('types')
    for k,v in atom_types.items():
        masses = list(set(v.masses))
        if len(masses) > 1: break
        try:
            atom_map[k] = per_table[masses[0]]
        except KeyError:
            closest = min(per_table, key=lambda x:abs(x-masses[0]))
            atom_map[k] = per_table[closest]
    return atom_map


