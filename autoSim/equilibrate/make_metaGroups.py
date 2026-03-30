
#####
def make_metaGroups(types,):
    min_lines = []
    with open('system.in.groups','w') as f:
        t_atoms = 0

        solvent_groups = []
        molecule_groups = []
        solvents = f' group solvent union'
        for k in types.keys():
            name = k.split('_')
            if name[0] == 'solvent':
                solvent_groups.append(types[k])
                t_atoms += len(types[types[k]])
                min_lines.append(f' group {types[k]} type {t_atoms}')
                min_lines.append(f' group {types[k]} include molecule')
                if name[1] == 0:
                    min_lines.append(f' group solvent type {t_atoms}')
                else:
                    solvents += f' {types[k]}'

        min_lines.append(solvents)
        min_lines.append('')

        solutes = f'# group solute union'
        for k in types.keys():
            name = k.split('_')
            if name[0] == 'molecule':
                molecule_groups.append(types[k])
                t_atoms += len(types[types[k]])
                min_lines.append(f'# group {types[k]} type {t_atoms}')
                min_lines.append(f'# group {types[k]} include molecule')
                if name[1] == 0:
                    min_lines.append(f' group solute type {t_atoms}')
                else:
                    solutes += f' {types[k]}'
                
        min_lines.append(solutes)
        min_lines.append('')

        min_lines.append(' group solute subtract all solvent')
        min_lines.append('')

        for line in min_lines:
            f.write(''.join(line))
            f.write('\n')
#####
if __name__ == '__main__':
    import json
    
    with open('meta_atomTypes.json') as f:
        types = json.load(f)

    make_metaGroups(types)

