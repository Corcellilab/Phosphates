import matplotlib.pyplot as plt
from aggPy import *

hbonds = Analysis('MDin.json', 'hbond')
hbonds.hCoordination()

i = Analysis.aggregate(hbonds, 'Distance')
plt.hist(i, bins=100)
plt.show()

j = Analysis.std_dev(hbonds, 'Coordinations')
print(j)

x = Analysis.average(hbonds, 'Coordinations')
print(x)

mol_ct, sys_ct = hbonds.timeCorr(hbonds)
plt.plot(list(sys_ct.keys()), list(sys_ct.values()))
plt.show()

