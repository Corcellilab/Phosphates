
import networkx as nx
import pandas as pd
import matplotlib.pyplot as plt
import sys

import argparse

from analyze_csv import *

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-v' ,'--verbose', action='store_true')
    parser.add_argument('-n' ,'--networks', action='append', help='Network prefix of csv file produced')
    parser.add_argument('-l' ,'--node_lst', type=str, default='all',
                        help='The nodes in a network that are desired. Use "all" to look at all nodes.')
    parser.add_argument('-o' ,'--outfile', type=str, default='results', help='Outfile name for csvs')
    parser.add_argument('-adj' ,'--adj_mat', action='store', 
                        help='Type of adjacency matrix to analyze from each network. Types include "simple_resind", "dir_resind", etc.')
    parser.add_argument('-di','--diGraph', action='store_true', help='Make a graph directed in analysis by providing flag')
    parser.add_argument('-rm','--rm_empty_nodes', action='store_true',
                        help='Provide flag to remove empty nodes from analysis. Needed for eigenvalue and path analysis.')
    parser.add_argument('-w','--weight', action='store_true', help='Use the edge weights in analysis by providing flag.')
    parser.add_argument('-p','--path', action='store_true', help='Calculate the average path length between two nodes in network.')
    parser.add_argument('-e','--eigs', action='store_true', help='Calculate eigenvalues of connected clusters.')
    args = parser.parse_args()

    network_lst = args.networks
    node_lst = args.node_lst
    outfile = args.outfile

    kwargs = {"diGraph": args.diGraph, 
              "rm_empty_nodes": args.rm_empty_nodes, 
              "weight": args.weight,
              "calc_path_length": args.path, 
              "calc_eigs": args.eigs,
              }
    
    results = analyze_csvs(network_lst, node_lst, args.adj_mat, **kwargs)
    dump_analyze_csvs(results, args.outfile)

