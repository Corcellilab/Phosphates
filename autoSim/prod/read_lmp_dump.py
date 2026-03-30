import numpy as np

def read_lmp_dump(filename):
    data = {}
    with open(filename,'r') as f:
        for line in f:
            l = line.strip().split()
            if l[0] == '#': continue
            ind = l.index('=')
            labels = l[:ind]
            values = l[ind+1:]
            for i,label in enumerate(labels):
                try: data[label].append(float(values[i]))
                except KeyError: data[label] = [float(values[i])]

    return data

#####
if __name__ == '__main__':
    filename = 'defrost.txt'
    read_lmp_dump(filename)

