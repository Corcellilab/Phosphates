import aggPy
import aggPy.workup


json_file = 'MDin.json'
analysis = 'spectra'
aggPy.main(json_file, analysis)  #Generates csvs

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
results = aggPy.workup.analyze_csvs(network_lst, node_lst, adj_mat, **kwargs)

#Calc avg of graph props - dumps dataframe 
aggPy.workup.dump_analyze_csvs(results, outfile)

#Calc C(t) of any hbond
aggPy.time_corr_Ct.binary(network_lst[0], adj_mat,
        min_dt=0, max_dt=100, skip_dt=1,
        t_prop=1, dt=1e-15)

aggPy.plot_csv()
