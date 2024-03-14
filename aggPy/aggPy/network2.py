import numpy as np

###
def dij(adj_list, visited, node, result, key):
    visited.add(node)
    result[key].append(node)
    for neighbor in adj_list[node]:
        if neighbor not in visited:
            dij(adj_list, visited, neighbor, result, key)

def connectivity(data, nodes):
    data = np.tril(data)
    z = np.nonzero(data)
    agg_tuple = [tuple((row,col)) for row, col in zip(z[1],z[0])]

    adj_list = {k: [] for k in range(len(data))}
    for x, y in agg_tuple:
        try:
            adj_list[x].append(y)
            adj_list[y].append(x)
        except KeyError:
            adj_list[x] = [y]
            adj_list[y] = [x]

    result = {k: [] for k in range(len(data))}
    visited = set()
    for node in adj_list:
        if node not in visited:
            dij(adj_list, visited, node, result, node)

    out = result.values()
    out = [i for i in out if i != []]
    
    output = []
    for i in out:
        lst = [nodes[index] for index in i]
        output.append(lst)

    return output

###

def network2(lst):
    edge_lst = list({*map(tuple, map(sorted, lst))})        #Removes duplicates and bifurcated - Unique hbonds to resid
    nodes_total = [item for sublist in edge_lst for item in sublist]
    nodes = sorted(set(nodes_total))
    
    arra = np.zeros((len(nodes),len(nodes)),dtype=int)
    for x,y in edge_lst:
       row = nodes.index(x)
       col = nodes.index(y)
       arra[row,col] = 1
       arra[col,row] = 1    
    
    Coordinations = {}
    Chain = {}
    degree = {}
    #print(arra)
    for i in range(0,len(nodes)):
        count = np.count_nonzero(arra[i,:])
        try:
            Coordinations[count] += 1
        except KeyError:
            Coordinations[count] = 1
        degree[nodes[i]] = count

    chains = connectivity(arra,nodes)
    
    for i in chains:
        try:
            Chain[len(i)] += 1
        except KeyError:
            Chain[len(i)] = 1
    
    output = {}
    output['Coordinations'] = Coordinations
    output['Aggregate Size'] = Chain
    output['Aggregate resids'] = chains
    output['Degree'] = degree

    return output

#lst = [('11','12'),('12','13'),('12','14')
#        ,('15','16'),('16','15'),
#        ]       #will be resid's

#lst = [('a','b'),('a','c'),('d','e')]

#c = network2(lst)
#print(c)

