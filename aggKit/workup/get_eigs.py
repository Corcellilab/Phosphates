#get_eigs.py
"""

    Get the eigenvalue of a cluster
    from a networkx graph

"""

from scipy import linalg as LA
import networkx as nx
import numpy as np
import sys

def get_eigs(graph, nodes_of_interest):
    clusters = nx.connected_components(graph)
    edges = graph.edges
    e_values = list()
    for cluster in clusters:
        cluster_adj = list()
        for edge in edges:
            if (edge[0] in cluster) or (edge[1] in cluster):
                if (edge[0] in nodes_of_interest) or (edge[1] in nodes_of_interest):
                    cluster_adj.append(edge)
        G = nx.Graph(cluster_adj)
        eigs = nx.normalized_laplacian_matrix(G)
        e = np.linalg.eigvals(eigs.toarray())
        e_values.append(list(e))

    return [x for xs in e_values for x in xs]

