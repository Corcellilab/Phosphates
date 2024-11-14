#get_dir_degrees.py
""""

    Get degree of nodes in directed graph
    depending on specified row or column

"""

import networkx as nx
import numpy as np
import sys
import matplotlib.pyplot as plt

def get_dir_degrees(G, nodes_of_interest, weight=False):
    adj_matrix = nx.to_numpy_array(G)
    
    in_degrees = []
    out_degrees = []
    tot_degrees = []
    for node in nodes_of_interest:
        in_degree = G.in_degree(node, weight=weight)
        out_degree = G.out_degree(node, weight=weight)
        tot_degree = in_degree + out_degree 
        if in_degree != 0: in_degrees.append(in_degree)
        if out_degree != 0: out_degrees.append(out_degree)
        if tot_degree != 0: tot_degrees.append(tot_degree)
    
    return in_degrees, out_degrees, tot_degrees

