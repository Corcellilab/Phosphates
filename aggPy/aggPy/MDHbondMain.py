import MDAnalysis as mda
import numpy as np
import aggPy.network2   

def hbond(topology=[],donors=[],hydrogens=[],acceptors=[],HbondParm=[],box=[],pre_allo=[],mol_delim='resid',XYZ=[False]):
    if XYZ[0] == True:
        xyz_string = ''
        atom_map = XYZ[1]
        XYZfilename = XYZ[2]
        num_atoms = 0

    dist = mda.lib.distances.distance_array(hydrogens.positions,acceptors.positions,
            box=box,result=pre_allo)
    index = np.column_stack(np.where((dist >= HbondParm[0]) & (dist <= HbondParm[1])))

    output = {'Distance':[],'Angle':[],'Network':[]}
    for i in index:
        if hydrogens[i[0]].resid != acceptors[i[1]].resid:
            hydrogen = hydrogens[i[0]]
            acceptor = acceptors[i[1]]
            distance = dist[i[0]][i[1]]
            bonded = hydrogen.bonds[0]
            h_bond_atom = bonded.partner(hydrogen)
            angle = (mda.lib.distances.calc_angles(h_bond_atom.position,hydrogen.position,acceptor.position,
                box=box))         #radians
            if angle > HbondParm[2]:
                output['Distance'].append(distance)
                output['Angle'].append(angle)
                output['Network'].append((hydrogen.resid,acceptor.resid))

            if XYZ[0] == True:
                rez = [h_bond_atom.resid, acceptor.resid]
                for i in range(0,len(rez)):
                    rez_h = rez[i]
                    num_atoms += len(topology[rez_h])
                    for atom in topology[rez_h]:
                        symbol = atom_map[atom.type]
                        pos = atom.position
                        xyz_string += f'{symbol} {pos[0]} {pos[1]} {pos[2]}\n'
    if XYZ[0] == True:
        with open(XYZfilename,'a') as f:
            f.write(f'{num_atoms}\n')
            f.write(f'\t comment\n')
            f.write(xyz_string)

    if output['Network'] != []:
        net = aggPy.network2(output['Network'])
        output['Coordinations'] = net['Coordinations']
        output['Aggregate Size'] = net['Aggregate Size']
        output['Aggregate resids'] = net['Aggregate resids']
        output['Node Degree'] = net['Degree']

    class HbondCalc:
        def __init__(self, output):
            self.__dict__.update(**output)

    return HbondCalc(output)



