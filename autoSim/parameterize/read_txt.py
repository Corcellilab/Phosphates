

with open('sym_amine.txt','r') as f:
    lines = [line.strip().split() for line in f]

with open('sym_amine.xyz','w') as f:
    f.write(f'{len(lines)-2}\n')
    f.write(f'\tcomment\n')
    for line in lines:
        try:
            name = line[0].split(':')[1][0].capitalize()
            x = line[-3]
            y = line[-2]
            z = line[-1]
            f.write(f'{name} {x} {y} {z}\n')
        except IndexError:
            pass



