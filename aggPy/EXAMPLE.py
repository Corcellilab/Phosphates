import matplotlib.pyplot as plt
from aggPy import *

hbonds = Analysis('example.json', 'hbond')
hbonds.hCoordination()

dists = hbonds.aggregate('Distance')
plt.hist(dists, bins=100)
plt.show()

hbonds.std_dev('Coordinations',bin_width=25)
print(hbonds.CoordinationsStdev)

hbonds.average('Coordinations')
print(hbonds.CoordinationsAvg)

mol_ct, sys_ct = hbonds.timeCorr()
plt.plot(list(sys_ct.keys()), list(sys_ct.values()))
plt.show()
