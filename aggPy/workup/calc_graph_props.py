#calc_graph_props
"""
    Calculate network props at step
"""

import networkx as nx
from numpy import log2

from .get_degree_dist import *
from .get_dir_degrees import *
from .get_graph_entropy import *
from .get_avg_length import *
from .get_eigs import *

def calc_graph_props(G, nodes_of_interest, diGraph=False, rm_empty_nodes=False,
        weight=False, calc_path_length=False, calc_eigs=False,):
    
    path_lengths = list()
    entropy = list()
    eigs = list()
    square_clusters = list()
    cluster_coeffs = list()
    avg_clustering = list()
    ej_weights = list()
    hydrogen_degrees = list()
    acceptor_degrees = list()
    tot_degrees = list()
    
    #Average Edge Weight
    for edge in G.edges(data=True):
        ej_weights.append(edge[2]['weight'])

    degrees = get_degree_dist(G, nodes_of_interest)

    #Remove empty nodes
    if rm_empty_nodes == True:
        nodes_with_no_edges = [node for node in G.nodes() if G.degree(node) == 0]
        G.remove_nodes_from(nodes_with_no_edges) 
        nodes_of_interest = list(set(nodes_of_interest).intersection(set(G.nodes)))
        if diGraph == False:
            if calc_path_length == True: 
                path_lengths = get_avg_length(G, nodes_of_interest, weight=weight)
            if calc_eigs == True: 
                eigs = get_eigs(G, nodes_of_interest)

    if len(G.nodes) > 0:
        square_clusters = list(nx.square_clustering(G, nodes=nodes_of_interest).values())
        cluster_coeffs = list(nx.clustering(G, nodes=nodes_of_interest).values())
        avg_clustering = nx.average_clustering(G, nodes=nodes_of_interest, count_zeros=True)
        entropy = get_graph_entropy(G, nodes_of_interest, weight=weight)    

    if diGraph == True:
        in_d, out_d, tot_d = get_dir_degrees(G, nodes_of_interest, weight=weight)
        hydrogen_degrees += out_d
        acceptor_degrees += in_d
        tot_degrees += tot_d

    return {
            'degrees': degrees,
            'edge_weights': ej_weights, #'avg_edge_weight': np.mean(ej_weights),
            'path_lengths': path_lengths, #'avg_path_length': np.mean(path_lengths),
            'eigs': eigs, 
            'square_clusters': square_clusters, #'avg_square_cluster': np.mean(square_clusters),
            'triangle_clusters': cluster_coeffs, #'avg_tri_cluster': np.mean(cluster_coeffs),
            #'cluster_coeff': avg_clustering,
            'entropy': entropy,
            'hydrogen_degrees': hydrogen_degrees,
            'acceptor_degrees': acceptor_degrees,
            'tot_degrees': tot_degrees,
           }

