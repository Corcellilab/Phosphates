import numpy as np

def timeCorr(self, min_dt=1, max_dt=20, skip_dt=1, t_prop=1, num_mols=1):
    counter = 0
    time_corr = {}
    for i in range(0, self.length):
        total = getattr(getattr(self, f'result{i}'), 'Aggregate resids')
        IDs = [item for sublist in total for item in sublist]
        for ID in IDs:
            try:
                time_corr[ID][counter] = 1
            except KeyError:
                time_corr[ID] = np.zeros(self.length, dtype=int)
                time_corr[ID][counter] = 1
        counter += 1

    mol_counter = 0
    Time_Corr = {}
    for k,v in time_corr.items():
        mol_counter += 1
        for dt in range(min_dt,max_dt,skip_dt):       #range(1,int(tot_steps/2)+1,10)
            t = 0
            l = np.array([],dtype=int)
            while t < self.length:
                a = v[t:t+dt]
                A = np.prod(a)
                l = np.append(l,A)
                t += t_prop                 #Spacer propagation

            l = np.sum(l) / len(l)
            try: Time_Corr[k][dt] = round(l/Time_Corr[k][1], 5)
            except KeyError: Time_Corr[k] = {dt:l}

        Time_Corr[k][1] = 1            #Normalization

    data = {}
    for v in Time_Corr.values():        #Averages over all molecules
        for key,value in v.items():
            try:
                data[key] += value/num_mols
            except KeyError:
                data[key] = value/num_mols
    print(str(mol_counter) + ' molecules time correlated')
    datum = {}
    for k,v in data.items():
        datum[k] = v/data[1]
    #return {'Molecule C(t)':Time_Corr,'System Averaged C(t)':datum}
    return Time_Corr, datum

