#get_avg_length.py
"""

    Return the average minimum path length
    from a networkx graph object

"""

import networkx as nx

def get_avg_length(graph, nodes_of_interest, weight=False):
    c = nx.connected_components(graph)
    clusters = [list(cluster) for cluster in c]
    avg_cluster_path = []
    for cluster in clusters:
        path_length = 0
        for i in range(len(cluster)):
            source_node = cluster[i]
            if source_node in nodes_of_interest:    #Nodes_of_interst term
                for j in range(i+1, len(cluster)):
                    target_node = cluster[j]
                    if weight == False:
                        path_length += len(nx.shortest_path(graph, source_node, target_node))
                    else:
                        path_length += len(nx.shortest_path(graph, source_node, target_node, weight='weight'))
        N = len(cluster)
        factor = 2 / (N*(N-1))
        avg_cluster_path.append(factor * path_length)
    
    return avg_cluster_path

