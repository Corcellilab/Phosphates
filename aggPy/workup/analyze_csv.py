#analyze_csvs.py

import pandas as pd
import numpy as np
import networkx as nx

from .calc_graph_props import *

####

def analyze_csvs(network_lst, node_lst, matrix,
                diGraph=False, rm_empty_nodes=False, weight=False,
                calc_path_length=False, calc_eigs=False,
                ):
    matrices = {}
    for i in range(len(network_lst)):
        df = pd.read_csv(f'{network_lst[i]}.csv')
        matrices[i] = df[:][matrix]

    results = {}
    length = int(0.25*len(matrices[0]))
    for i in range(len(matrices[0])):       #For every step
        if i%length == 0: print(f'Step {i}')
        if diGraph == False: G = nx.Graph()
        elif diGraph == True: G = nx.DiGraph()

        nodes_of_interest = list()
        for k in matrices.keys():
            if diGraph == False: g = nx.Graph(eval(matrices[k][i]))
            elif diGraph == True: g = nx.DiGraph(eval(matrices[k][i]))
            try:
                k == network_lst.index(node_lst)
                nodes_of_interest += g.nodes
            except ValueError:
                if node_lst == 'all':
                    nodes_of_interest += g.nodes
            G = nx.compose(G, g)
        
        results[i] = calc_graph_props(G, nodes_of_interest,
                    diGraph=diGraph, rm_empty_nodes=rm_empty_nodes, weight=weight,
                    calc_path_length=calc_path_length, calc_eigs=calc_eigs,
                   )

    return pd.DataFrame.from_dict(results, orient='index')

####
