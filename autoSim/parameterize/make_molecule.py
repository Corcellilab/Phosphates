from rdkit import Chem
from rdkit.Chem import AllChem 

def make_molecule(smiles=None, output='myMol'):
    if smiles == None:
        print('No smiles string found...')
        return 1

    mol = Chem.MolFromSmiles(smiles)

    mol = Chem.AddHs(mol)

    AllChem.EmbedMolecule(mol, AllChem.ETKDG())

    AllChem.UFFOptimizeMolecule(mol)

    xyz = Chem.MolToXYZBlock(mol)

    pdb = Chem.MolToPDBBlock(mol)

    with open(f'{output}.xyz','w') as f:
        f.write(xyz)

    with open(f'{output}.pdb','w') as f:
        f.write(pdb)

    return f'{output}.xyz'

#####
if __name__ == '__main__':
    make_molecule(
            smiles = "CN(C)C1=CC=C(/C=[NH+]/C2=CC=C(N(C)C)C=C2)C=C1",
            output='myTest',
        )


