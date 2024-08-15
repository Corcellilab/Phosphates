  hbonds = Analysis('./MDin.json', 'hbond')
  hbonds.hCoordination()

  n = Analysis.aggregate(hbonds, 'Cluster Coeff')

  e = Analysis.aggregate(hbonds, 'Entropy')

  i = Analysis.aggregate(hbonds, 'Distance')
  plt.hist(i, bins=100)
  plt.show()

  a = Analysis.aggregate(hbonds, 'Eigs')
  plt.hist(a, bins=100)
  plt.show()

  b = Analysis.average(hbonds, 'Aggregate Size')
  #print(f'Std Dev Length: {b}')

  k = Analysis.average(hbonds, 'Aggregate Size')
  #print(f'Avg Length: {k}')

  j = Analysis.std_dev(hbonds, 'Coordinations')
  #print(f'Std Dev Coords: {j}')

  x = Analysis.average(hbonds, 'Coordinations')
  plt.scatter(x.keys(), x.values())
  plt.show()
  #print(f'Avg Coords: {x}')
