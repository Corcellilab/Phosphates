#initial.py
"""
    
    Make Analysis object
    Read in json file
    Define self attributes

"""

import json

from .start_analysis import start_analysis
from ..network_analysis.make_ejs import make_ej
from .make_ag import make_ag 

class Analysis:
    def __init__(self, json_file, key):
        with open(json_file) as f:
            args = json.load(f)
            initial = args['init'][0]
            args = args[key][0]
        self.__dict__.update(**initial)
        self.__dict__.update(**args)
        start_analysis(self)

    def make_ags(self):
        make_ag(self)

    def make_ejs(self):
        make_ej(self)

