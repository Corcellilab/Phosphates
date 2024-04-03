import numpy as np

#########
def aggregate(self, variable):
    T = type(getattr(self.result0, variable))
    if T == list:
        variables = []
        for i in range(0, self.length):
            try:
                var = getattr(getattr(self, f'result{i}'), variable)
                variables.append(var)
            except AttributeError: print('Attr Error')
        return [v for va in variables for v in va]
    elif T == dict:
        v0 = getattr(self.result0, variable)
        for i in range(1, self.length):
            try:
                dN = getattr(getattr(self, f'result{i}'), variable)
                for k,v in dN.items():
                    try: v0[k] += v
                    except KeyError: v0[k] = v
            except AttributeError: print('Attr Error')
        return v0
    else:
        print('Variable format not dict or list')

#########
def average(self, variable):
    agg = self.aggregate(self, variable)
    T = type(getattr(self.result0, variable))
    if T == dict:
        for k,v in agg.items():
            agg[k] = v/self.length
    elif T == list:
        agg = sum(agg) / self.length

    setattr(self, f'{variable}Avg', agg)
    return agg

#########
def std_dev(self, variable, bin_width=1):
    cn_group = {}
    for i in range(0, self.length):
        for k in getattr(getattr(self, f'result{i}'), variable).keys():
            try: cn_group[k].append(getattr(getattr(self, f'result{i}'), variable)[k])
            except KeyError: cn_group[k] = [getattr(getattr(self, f'result{i}'), variable)[k]]

    trj_stdev = {}
    for k in cn_group.keys():
        bin_values = []
        for t in range(0, self.length, bin_width):
            x = np.mean(cn_group[k][t:t+bin_width])
            bin_values.append(x)
        bin_values = np.array(bin_values)
        bin_values[np.isnan(bin_values)] = 0
        trj_stdev[k] = np.std(bin_values)

    setattr(self, f'{variable}Stdev', trj_stdev)
    return trj_stdev

##########


