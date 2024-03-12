from MDAnalysis.analysis.base import (AnalysisBase, AnalysisFromFunction, analysis_class)

import numpy as np
import json

from MDHbondMain import hbond
from atom_map import atom_mapping
from initial_analysis import init_analysis 
from timeCorr import timeCorr
from workup import aggregate, average, std_dev

class Analysis:
    def __init__(self, json_file, key):
        with open(json_file) as f:
            args = json.load(f)
            initial = args['init'][0]
            args = args[key][0]
        self.__dict__.update(**initial)
        self.__dict__.update(**args)

    def hCoordination(self):
        init_analysis(self)
        self.atom_map = atom_mapping(self.universe.atoms)
        self.XYZ = [x for x in self.XYZ.split(',')]
        self.HbondParm = [float(x) for x in self.HbondParm.split(',')]
        
        hbonds =  AnalysisFromFunction(hbond,self.universe.trajectory,
            topology = self.topology,
            donors = self.donors,
            hydrogens = self.hydrogens,
            acceptors = self.acceptors,
            HbondParm = self.HbondParm,
            box = self.box,
            pre_allo = self.pre_allo,
            mol_delim = self.mol_deliminator,
            XYZ = [eval(self.XYZ[0]), self.atom_map, f'{self.XYZ[1].strip()}.xyz']
        )
        
        if self.run_length == 'n':
            self.length = int((self.end_frame - self.start_frame / self.frame_dt))
            hbonds.run(start=self.start_frame, stop=self.end_frame, step=self.frame_dt)
        else:
            self.length = int(len(self.universe.trajectory))
            hbonds.run()

        results = hbonds.results.timeseries
        output = {}
        for i in range(0,len(results)):
            output[f'result{i}'] = results[i]
        self.__dict__.update(**output)

    def vibSpectra(self):
        return self.bond

    def aggregate(self, variable):
        return aggregate(self, variable)

    def average(self, variable):
        average(self, variable)

    def std_dev(self, variable, bin_width=1):
        std_dev(self, variable, bin_width)
    
    def timeCorr(self, min_dt=1, max_dt=10, skip_dt=1, t_prop=1, num_mols=1):
        return timeCorr(self, min_dt, max_dt, skip_dt, t_prop, num_mols)

    def graphing(self, data):
        return None

######

hbonds = Analysis('MDin.json', 'hbond')
freq = Analysis('MDin.json', 'freq') 

hbonds.hCoordination()

#mol_ct, sys_ct = hbonds.timeCorr()
#print(sys_ct)
#plt.plot(list(sys_ct.keys()), list(sys_ct.values()))
#plt.show()

#hbonds.std_dev('Coordinations', bin_width=25)
#print(hbonds.CoordinationsStdev)
#hbonds.average('Coordinations')
#print(hbonds.CoordinationsAvg)

#dists = hbonds.aggregate('Distance')
#angles = hbonds.aggregate('Angle')
#angles = np.array(angles)*(180/np.pi)
#plt.hist(dists, bins=100)
#plt.gca().invert_xaxis()
#plt.show()

#c = hbonds.aggregate('Coordinations')
#print('#####')
#print(c)

#hbonds.averaging('Coordinations')
#print(hbonds.CoordinationsAvg)
#hbonds.averaging('Angle')
#print(hbonds.AngleAvg)


