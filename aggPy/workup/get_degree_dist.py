#get_degree_dist.py
"""

    Returns the normalized frequency of node degree 
    that occurs 

    returns: {degree: frequency,}

"""

def get_degree_dist(graph, nodes_of_interest):
    degrees = {}
    for node in nodes_of_interest:
        try: degrees[graph.degree[node]] += 1
        except KeyError: degrees[graph.degree[node]] = 1

    num_mols = len(nodes_of_interest)
    for k,v in degrees.items():
        degrees[k] = v/num_mols
    
    return degrees
