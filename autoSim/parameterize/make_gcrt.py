#"""
# Sample master_gcrt file:
#    [['%chk=sym.chk'], ['%NProcs=4'], ['%mem=16GB'], ['#HF/6-31G*', 'SCF=tight', 'Test', 'Pop=MK', 'iop(6/33=2)', 'iop(6/42=6)', 'opt'], ['#', 'iop(6/50=1)'], [], ['remark'], [    ], ['charge', 'spin'], ['coords'], [], ['sym.gesp'], [], ['sym.gesp'], []]
#
#    %chk=sym.chk
#    %NProcs=4
#    %mem=16GB
#    #HF/6-31G* SCF=tight Test Pop=MK iop(6/33=2) iop(6/42=6) opt
#    # iop(6/50=1)
#
#    remark
#
#    charge   spin
#    coords
#
#    sym.gesp
#
#    sym.gesp 
#
#
#"""


from .read_file import read_file

#####

def make_gcrt(xyz_file='', master_gcrt='', output='output',charge=0, spin=1):
    with open(xyz_file,'r') as f:
        coords = [line.split() for line in f]
    coords = coords[2:]

    #with open(master_gcrt,'r') as f:
    #    gcrt = [line.split() for line in f]

    gcrt = read_file(master_gcrt)

    ind = gcrt.index(['charge','spin'])
    gcrt[ind][0] = str(charge)
    gcrt[ind][1] = str(spin)

    ind = gcrt.index(['%chk=sym.chk'])
    gcrt[ind][0] = f'%chk={output}.chk'

    ind = gcrt.index(['sym.gesp'])
    gcrt[ind][0] = f'{output}.gesp'
    gcrt[ind+2][0] = f'{output}.gesp'

    ind = gcrt.index(['coords'])
    del gcrt[ind]
    [gcrt.insert(ind,c) for c in coords]

    with open(f'{output}.gcrt','w') as f:
        for item in gcrt:
            f.write(' '.join(item))
            f.write('\n')

    return f'{output}.gcrt'

######
if __name__ == "__main__":
    make_gcrt(
            xyz_file='myMol.xyz',
            master_gcrt='master.gcrt',
            output='output.gcrt',
            charge=0,
            spin=1,
        )
