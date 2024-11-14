#main.py
"""

    Run the given analysis
    from the given json file

"""

from MDAnalysis.coordinates.XYZ import XYZWriter
import matplotlib.pyplot as plt
import numpy as np
import networkx as nx
import json
import pandas as pd
import sys

from .initial.initial import Analysis

from .network_analysis.make_ejs import make_ej

from .vacf.time_corr import time_corr_vector
from .vacf.fourier import ftransform

def main(json_file, analysis):
    data = Analysis(json_file, analysis)
    if analysis == 'network':
        data.make_ags()
        data.make_ejs()
        return data

    elif analysis == 'spectra':
        dt = data.dt
        if data.run_length == 'n':
            run_length = (data.end_frame - data.start_frame) / data.frame_dt
        elif data.run_length == 'y':
            run_length = len(data.universe.trajectory) / data.frame_dt

        data.get_velocities(dt)
        
        x, y = time_corr_vector(
                data.vacf_1,
                min_dt=0, max_dt=run_length, skip_dt=1,
                t_prop=100
                )

        corr = dict(zip(x,y))
        with open(f'vacfCt.txt','w') as f:
            f.write(str(corr))

        x, y = ftransform(
                x, y,
                dt=dt
                )

        spec = dict(zip(x,y))
        with open('spectra.txt','w') as f:
            f.write(str(spec))

        return spec
    
######

#if __name__ == "__main__":
#    json_file = sys.argv[1]
#    analysis = sys.argv[2]
#    data = main(json_file, analysis)

    #graphs = [
    #        (data.network_1.results.graphs, 'simple_resind'),
    #        (data.network_2.results.graphs, 'simple_resind'),
    #        (data.network_1.results.graphs, 'dir_resind'),
    #        (data.network_2.results.graphs, 'simple_resind'),
    #        (data.network_1.results.graphs, 'molecular_potential'),    #set calc_dipole: 0 in MDin.json
    #        (data.network_1.results.graphs, 'atomic_potential'),
    #        (data.network_1.results.graphs, 'simple_atomic'),
    #        ]

    #time_corr_bin(
    #        (data.network_1.results.graphs, 'simple_resind'), 
    #        run_length,
    #        min_dt = 0,
    #        max_dt = 20000,
    #        skip_dt = 1,
    #        t_prop = 1
    #        )

########

