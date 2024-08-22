import numpy as np
from numpy import inf
from math import log2
import sys

###
def avg_min_path(path):
    agg_length = 0
    for p in path:
        length = sum(p.values())
        factor = 1/(len(p.keys())*(len(p.keys())-1))
        agg_length += length

    return agg_length / (len(p)*(len(p)-1))

###
def cluster_size_label(path):
    cluster_size = {}
    clusters = []
    for agg in path:
        ag = set(list(agg.keys()))
        if ag not in clusters:
            try: cluster_size[len(ag)] += 1
            except KeyError: cluster_size[len(ag)] = 1
            clusters.append(ag)
    
    return clusters, cluster_size

###
def cluster_coeff(graphDict):
    Cj = 0
    for node, edges in graphDict.items():
        edge = 0
        degree = len(edges.keys())
        if degree == 1: continue
        neighborhood = list(edges.keys())
        for neighbor in neighborhood:
            neigh_edges = set(list(graphDict[neighbor].keys()))
            edge += len(neigh_edges.intersection(set(neighborhood)))

        Cj += (2*edge) / (degree*(degree-1))

    return (Cj / len(graphDict.keys()))

###
def entropyCalc(Coordinations):
    Sg = 0
    norm = sum(Coordinations.values())
    for degree, count in Coordinations.items():
        pj = count / norm
        Sg += -1*(pj * log2(pj))
    
    return Sg

###
def Dijkstra(graph, nodes, start):
  unvisited = {node: None for node in nodes}
  visited = {}
  currentDistance = 0
  unvisited[start] = currentDistance

  while True:
    for neighbour, distance in graph[start].items():
      if neighbour not in unvisited: continue
      newDistance = currentDistance + distance
      if unvisited[neighbour] is None or unvisited[neighbour] > newDistance:
        unvisited[neighbour] = newDistance
    visited[start] = currentDistance
    del unvisited[start]
    if not unvisited: break
    candidates = []
    candidates = [node for node in unvisited.items() if node[1]]
    if candidates == []:
      return visited
    else:
      start, currentDistance = sorted(candidates, key = lambda x: x[1])[0]

###
def makeGraphDict(resid_adj):
  ej = np.column_stack(np.nonzero(resid_adj))
  resid_adj = 1/resid_adj                       #2 conncections = 0.5 ej weight
 
 #print(np.allclose(resid_adj,resid_adj.T, 1e-05, 1e-8))
  
  G = {}
  for i, row in enumerate(resid_adj):
    connections = list(np.where((row < inf) & (row > 0)))
    for connect in connections:
      for c in connect:
        try: G[i][c] = row[c]
        except KeyError: G[i] = {c: row[c]}
        try: G[c][i] = resid_adj[c][i]
        except KeyError: G[c] = {i: resid_adj[c][i]}
  
  return dict(sorted(G.items()))

####

def network2(atom_ej, universe_adj):
    atom_adj_index = universe_adj[0].copy()
    atom_adj = universe_adj[1].copy()
    resid_adj_index = universe_adj[2].copy()
    resid_adj = universe_adj[3].copy()
   
    for edge in atom_ej:
      node1 = atom_adj_index.index(edge[0].id)
      node2 = atom_adj_index.index(edge[1].id)
      atom_adj[node1, node2] += 1
      atom_adj[node2, node1] += 1

      node1 = resid_adj_index.index(edge[0].resid)
      node2 = resid_adj_index.index(edge[1].resid)
      resid_adj[node1, node2] += 1                  #Allows values > 1
      resid_adj[node2, node1] += 1

    #Degree for resids - Can also be obtained from graphDict
    Coordinations = {}
    for i in range(0, len(resid_adj_index)):
        count = np.count_nonzero(resid_adj[i,:])
        try: Coordinations[count] += 1
        except KeyError: Coordinations[count] = 1

    #Shortest path for each node to other nodes ONLY if node is in aggregate
    graphDict = makeGraphDict(resid_adj)
    nodes = graphDict.keys()
    path = []
    for source_node in nodes:
        path.append(Dijkstra(graphDict, nodes, source_node))
    
    #Get aggregate labels and sizes
    cluster_label, cluster_size = cluster_size_label(path)

    #Get aggregate labels and sizes
    output = {}
    if path[0] == None:
        cluster_label = 0
        cluster_size = 0
        output['Coordinations'] = Coordinations
        output['Aggregate Size'] = [0]
        output['Aggregate resids'] = [0]
        output['Cluster Coeff'] = 0.0
        output['Path'] = [0]
        output['Eigs'] = np.linalg.eigvalsh(resid_adj)
        output['Entropy'] = 0.0
    else:
        cluster_label, cluster_size = cluster_size_label(path)
        output['Coordinations'] = Coordinations
        output['Aggregate Size'] = cluster_size
        output['Aggregate resids'] = cluster_label
        output['Cluster Coeff'] = cluster_coeff(graphDict)
        output['Path'] = avg_min_path(path)
        output['Eigs'] = np.linalg.eigvalsh(resid_adj)
        output['Entropy'] = entropyCalc(Coordinations)

    return output

####

