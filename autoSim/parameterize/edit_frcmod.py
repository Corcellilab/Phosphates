
def edit_frcmod(master_frcmod='master.frcmod',frcmod=''):
    add_mod = {
            'MASS':['__','0','0'],
            'BOND':['__-__','0','0'],
            'ANGLE':['__-__-__','0','0','0'],
            'DIHE':['__-__-__-__','1','1','1','1'],
            'IMPROPER':['__-__-__-__','0','0','0'],
            'NONBON':['__','0','0'],
        }
    with open(frcmod,'r') as f:
        lines = [line.strip().split() for line in f]

    for i, line in enumerate(lines):
        try:
            if line[0] in add_mod.keys() and lines[i+1] == []:
                lines[i+1] = add_mod[line[0]]
                lines.insert(i+2, [])
        except IndexError:
            pass

    with open(frcmod,'w') as f:
        for line in lines:
            f.write(' ' .join(line))
            f.write('\n')

#####

if __name__ == '__main__':
    edit_frcmod(frcmod='myTest.frcmod',)

