#main.py
"""

    Run the given analysis
    from the given json file

"""

#from MDAnalysis.coordinates.XYZ import XYZWriter
#import matplotlib.pyplot as plt
#import numpy as np
#import networkx as nx
#import json
#import pandas as pd
#import sys

from .initial.initial import Analysis

from .network_analysis.make_ejs import make_ej

######

def main(json_file):
    data = Analysis(json_file, 'network')
    data.make_ags()
    data.make_ejs()
    return data

######

