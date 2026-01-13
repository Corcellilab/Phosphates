#make_ejs.py

"""

    Main runner to calculate edges
    Can have multiple selections of edge list criteria

"""

from MDAnalysis.analysis.base import (AnalysisBase, AnalysisFromFunction, analysis_class)
import numpy as np
import networkx as nx

from .CalcEjs import CalcEjs 

def make_ej(self):
    for k,v in vars(self).items():
        name = k.split('_')
        if name[0] == 'network':
            pre_allo = np.zeros([len(v['group1']),len(v['group2'])])
            g1 = getattr(self, v['group1'])
            g2 = getattr(self, v['group2'])
           
            edges = CalcEjs(g1, g2, k, **vars(self))
            
            if self.run_length == 'n':
                length = int((self.end_frame - self.start_frame) / self.frame_dt)
                edges.run(start=self.start_frame, stop=self.end_frame,
                        step=self.frame_dt, verbose=self.verbose)
            else:
                length = int(len(self.universe.trajectory) / self.frame_dt)
                edges.run(step=self.frame_dt, verbose=self.verbose)

            setattr(self, k, edges)

######

