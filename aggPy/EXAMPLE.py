import matplotlib.pyplot as plt
import sys

from aggPy import *

hbonds = Analysis('../MDin.json', 'hbond')
hbonds.hCoordination()

x, y = Analysis.timeCorr(hbonds)

plt.plot(x, y)
plt.show()

i = Analysis.aggregate(hbonds, 'Distance')
plt.hist(i, bins=100)
plt.show()

j = Analysis.std_dev(hbonds, 'Coordinations')
print(f'Std Dev Coords: {j}')

k = Analysis.average(hbonds, 'Coordinations')
print(f'Avg Coords: {k}')
