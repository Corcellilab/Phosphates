#gen_sys.py

import random
import scipy
import numpy as np
import sys

num_wat = 2500
box_side = 50
min_dist = 1.5


with open("system.lt","w") as f:
    f.write('# -- System --\n\n')
    f.write('import "./spce.lt"\n')
    f.write('import "sym_guest.lt"\n')
    f.write('import "cb7.lt"\n')
    f.write('import "cl.lt"\n\n')
    f.write('write_once("Data Boundary") {\n')
    f.write(f'\t 0 {box_side} xlo xhi\n')
    f.write(f'\t 0 {box_side} ylo yhi\n')
    f.write(f'\t 0 {box_side} zlo zhi\n')
    f.write('}\n\n')

    f.write(f'cb7 = new CB7.rot(45, 0,1,0).move(30,20,25)\n')
    f.write(f'guest = new sym_guest.move(15,25,25)\n')
    f.write(f'cl = new cl.move(25,5,0)\n')

    cords = [[0,0,0]]
    for i in range(num_wat):
        while True:
            x, y, z = (random.uniform(0, box_side) for _ in range(3))
            b = [x,y,z]
            dist = scipy.spatial.distance.cdist(cords, np.array([b]))
            if np.any(dist < min_dist):
                pass
            else:
                cords.append(b)  
                theta = random.uniform(0, 360)
                f.write(f"wat_{i} = new SPCE.rot({theta}, 1,0,0).move({x}, {y}, {z})\n")
                break

