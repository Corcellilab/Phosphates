#get_graph_entropy.py
"""

    Returns scalar of entropy
    from networkx graph obj

"""

from numpy import log2

def get_graph_entropy(graph, nodes_of_interest, weight):
    degrees = {}
    for node in nodes_of_interest:
        try: degrees[graph.degree(node, weight=weight)] += 1
        except KeyError: degrees[graph.degree[node]] = 1

    S = 0
    norm = sum(degrees.values())
    for v in degrees.values():
        p = v / norm
        S += p*log2(p)

    return -1*S

