

def edit_submit(master_submit='', output='', pe='smp 4', q='l', N='myJob', modules=[], submit_lines=''):
    with open(master_submit, 'r') as f:
        lines = [line.strip().split() for line in f]

    ind = lines.index([])
    lines = lines[:ind]

    for line in lines[1:]:
        if line[1] == '-pe':
            pe = pe.split()
            line[2] = pe[0]
            line[3] = pe[1]
        elif line[1] == '-q':
            line[2] = q
        elif line[1] == '-N':
            line[2] = N

    lines.append([])
    for module in modules:
        lines.append([f'module load {module}'])

    lines.append([])
    for submit_line in submit_lines:
        lines.append([f'{submit_line}'])

    with open(output, 'w') as f:
        for line in lines:
            f.write(' ' .join(line))
            f.write('\n')

#####
if __name__ == '__main__':
    modules = [
            'lammps',
            'gaussian',
            'python',
        ]

    submit_lines = [
            'lmp -in in.prod',
            'g16 myTest.gcrt',
            'python my_script.py',
        ]

    edit_submit(
            master_submit='master_submit.sh',
            output='mySub.sh',
            modules=modules,
            submit_lines=submit_lines,
            pe='smp 2',
            q='myQ',
            N='myN',
        )

