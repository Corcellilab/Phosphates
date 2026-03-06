from aggKit import *
from aggKit.workup import analyze_csvs
from aggKit.workup import dump_analyze_csvs
import sys

"""

  This is a sample script that can be run with
  $python sample_run.py

"""

json_file = 'MDin.json'
main(json_file)  #Generates csvs

kwargs = {
          "diGraph": False,
          "rm_empty_nodes": True,
          "weight": False,
          "calc_path_length": False,
          "calc_eigs": False,
         }

network_lst = [
        'network_1',     #csv prefixes
        ]

node_lst = 'all'            #or 'network_1'
adj_mat = 'simple_resind'   #adj_mat to analyze
outfile = 'myResults'

#Calc graph properties - results = pandas dataframe object
results = analyze_csvs(network_lst, node_lst, adj_mat, **kwargs)

#Calc avg of graph props - dumps dataframe 
print(results)
results = dump_analyze_csvs(results, outfile)
print(results)

#Calc C(t) of any hbond
#time_corr_Ct.binary(network_lst[0], adj_mat,
#        min_dt=0, max_dt=100, skip_dt=1,
#        t_prop=1, dt=1e-15)

#Run command line interactable plotting tool for produced csvs
#plot_csv()

