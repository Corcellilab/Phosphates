#time_corr_Ct.py
"""

    Time correlation for Hbond
    1 or 0 values for property
    Degree time correlation

"""

import matplotlib.pyplot as plt
import numpy as np
import networkx as nx
import pandas as pd
import sys

def get_time_series(network, matrix):
    # Return the time series values for a given graph
    # edges: dict with {(n1, n2): np.array[1 or 0 if Hbond between two nodes exists]}
    # degrees: dict with {n1: np.array[degree at given step],}
  
    print('Calculating Time Series...')

    df = pd.read_csv(f'{network}.csv')[matrix]
    run_length = len(df)

    edges = dict()
    degrees = dict()
    for i in range(0, run_length):
        if i % 500 == 0: print(i)
        G = nx.Graph(eval(df[i]))
        connections = dict(sorted(nx.to_dict_of_lists(G).items()))

        for n0, node_lst in connections.items():
            if node_lst == []: continue
            try:
                degrees[n0][i] = G.degree[n0]
            except KeyError:
                degrees[n0] = np.zeros(run_length, dtype=int)
                degrees[n0][i] = G.degree[n0]
            
            for node in node_lst:
                if int(n0) < int(node):
                    try:
                        edges[(n0, node)][i] = 1
                    except KeyError:
                        edges[(n0, node)] = np.zeros(run_length, dtype=int)
                        edges[(n0, node)][i] = 1
                else:
                    try:
                        edges[(node, n0)][i] = 1
                    except KeyError:
                        edges[(node, n0)] = np.zeros(run_length, dtype=int)
                        edges[(node, n0)][i] = 1

    return edges, degrees

######

def time_corr_zero_and_one(edges, max_dt):
    values = [[] for i in range(max_dt)]
    keys = [i for i in range(max_dt)]
    data = dict(zip(keys, values))

    for bond, series in edges.items():
        t_length = 1
        nonzeros = list(np.nonzero(series))[0]
        for i in range(len(nonzeros)):
            c_value = nonzeros[i]
            try:  n_value = nonzeros[i+1]
            except IndexError: break
            if c_value + 1 == n_value:
                t_length += 1
            else:
                for k in data.keys():
                    if k <= t_length:
                        data[k].append(1)
                    else:
                        data[k].append(0)
    return data

#####

def time_corr_degrees(degrees, max_dt, max_degree):
    from copy import deepcopy

    values = [[] for i in range(max_dt)]
    keys = [i for i in range(max_dt)]
    indiv_data = dict(zip(keys, values))

    ds = [i for i in range(1, max_degree)]
    group_data = {}
    for i in ds:
        group_data[i] = deepcopy(indiv_data)

    for bond, series in degrees.items():
        max_degree = max(series)
        while max_degree >= 1:
            t_length = 1
            degree_matches = list(np.where(series==max_degree)[0])
            for i in range(len(degree_matches)):
                c_value = degree_matches[i]
                try:  n_value = degree_matches[i+1]
                except IndexError: break
                if c_value + 1 == n_value:
                    t_length += 1
                else:
                    for k in group_data[max_degree].keys():
                        if k <= t_length:
                            group_data[max_degree][k].append(1)
                        else:
                            group_data[max_degree][k].append(0)
            
            max_degree -= 1

    return group_data

######

def get_xy(diction, dt):
    x = list()
    y = list()
    for k, v in diction.items():
        y.append(np.average(v))
        x.append(k*dt)
    return x, y

######

def binary(network, matrix, min_dt=0, max_dt=100, skip_dt=1, t_prop=1,
        dt=1e-15):

    #edges: dict with Hbond 1 or 0 series
    #degrees: dict with node: degree series
    edges, degrees = get_time_series(network, matrix) 

    print('Calculating Correlation...')
    
    #data = time_corr_degrees(degrees, max_dt, 11)
    #for degree, ct in data.items():
    #    if ct[0] != []:
    #        x, y = get_xy(ct, dt)
    #        plt.plot(x, y, label=degree)
    
    data = time_corr_zero_and_one(edges, max_dt)
    
    x, y = get_xy(data, dt)
    
    csv_info = dict(zip(x,y))
    df = pd.DataFrame([csv_info])
    df.to_csv(f'{matrix}Ct.csv')

    return df

######

if __name__ == "__main__":
    #binary('network_1', 'simple_atomic')
    binary(input('Csv prefix: '), input('Adj mat: '))
