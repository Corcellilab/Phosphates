from itertools import islice
import itertools
import math
import sys
import time
from scipy.spatial import distance_matrix
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
from collections import Counter
from collections import defaultdict

import json

########
def getMass(atom_type):
    if atom_type == 1:
        return 30.973762
    elif atom_type == 2 or atom_type == 3 or atom_type == 6:
        return 15.999
    elif atom_type == 4 or atom_type == 7:
        return 1.00784
    elif atom_type == 5:
        return 22.989769

########
def avgNested(dictionary):
    for key,value in dictionary.items():
        total = 0
        for k,v in value.items():
            total += v
        for k,v in value.items():
            dictionary[key][k] = v/total
    return dictionary
########
def angleCalc(x,y,z,q,w,e,i,j,k):
    A = [x-q,y-w,z-e]
    B = [i-q,j-w,k-e]
    a = np.array(A)
    b = np.array(B)
    
    bm = np.sqrt(b.dot(b))
    am = np.sqrt(a.dot(a))
    d = np.dot(A,B) / (bm*am)
    D = math.degrees(math.acos(d))
    
    return D
#####

#####
def distanceCalc(x, y, z, nxt_x, nxt_y, nxt_z):
    return ((float(x)-float(nxt_x))**2 + (float(y)-float(nxt_y))**2 + ((float(z)-float(nxt_z))**2))**(1/2)

#####
def datreader(filename):     #Returns dictionary of OH bonds
    #print('Start datreader')
    with open('periodicTable.txt') as f:
        table = eval(str(f.readline()))

    with open(filename) as dat:
        lines = (line.rstrip('r\n') for line in dat)
        lines = [line for line in lines]

        atomTypes = lines[3].split()
        atomTypes = int(atomTypes[0])

        num_bonds = lines[4].split()
        num_bonds = int(num_bonds[0])

        for i in range(0,len(lines)):
            if lines[i] == 'Masses':
                atom_lst = lines[i+2:(i+atomTypes+2)]
                atomMap = []
                for value in atom_lst:
                    pair = value.split()
                    atomMap.append((int(float(pair[1])),int(pair[0])))
                print(atomMap)
                typeMap = ['0']
                for value in atomMap:
                    sym = table[value[0]]
                    typeMap.append((sym,value[0]))
            if lines[i] == 'Bonds':
                diction = {}
                bond = lines[i+2:i+num_bonds+2]
                for b in bond:
                    data = b.split()
                    if len(data) == 4:
                        diction[str(data[-2])] = str(data[-1])
                        diction[str(data[-1])] = str(data[-2])

    return [diction,typeMap]
#####
#def dimer_typing(lst, contacts):
#    p1 = np.zeros(7,dtype=list)
#    p2 = np.zeros(7,dtype=list)
    
#    contacts_lst
#    for i in lst:
#        if i[0] ==  
#        in contacts_lst or i[1] in conct:
#            np.insert = []

#    dimer = np.array(distance_matrix(np.array(donors),np.array(acceptors))) 

#    p = {1:{1:0,2:0},2:{1:0,2:0}}
#    for phos in range(1,3):
#        row

#    sys.exit(1)

##########
def XYZwriter(cordpath, lst, i, filename):
    step_xyz = str()
    #step_xyz = full string to input
    counter = 0
    for cord in lst:
        if cord[5] in cordpath: #or cord[1] == '5':
            counter += 1
            atom_type = cord[1]                
            if atom_type == '1':
                atom_type = 'P' 
            elif atom_type == '2':
                atom_type = 'O'
            elif atom_type == '3':
                atom_type = 'O'
            elif atom_type == '4':
                atom_type = 'H'
            elif atom_type == '5':
                atom_type = 'Na'
            x = cord[2]
            y = cord[3]
            z = cord[4]
            step_xyz += atom_type + '  ' + x + '  ' + y + '  ' + z + '\n'
    with open(filename,'a') as new_file:
        new_file.write(str(counter) + '\n')
        new_file.write('     ' + str(i) + '\n')
        new_file.write(str(step_xyz))

#########
def radiiGyra(network, points, box_length):
    gyra_distances = []
    for i in network:                              #i = one aggregate
        if len(i) >= 2:
            M = 0
            indy = []
            theta = {}
            E = {0:[],1:[],2:[]}
            C = {0:[],1:[],2:[]}
            COM = []
            for particle_id in points.keys():
                if points[particle_id].mol_id in i:
                    atom = points[particle_id]
                    M += points[particle_id].mass
                    indy.append(particle_id)
                
                    for _ in range(0,3):
                        theta[_] = 2*math.pi*atom.cord[_] / box_length
                        E[_].append(atom.mass*math.cos(theta[_]))
                        C[_].append(atom.mass*math.sin(theta[_]))

            for _ in range(0,3):
                E[_] = (1/M)*sum(E[_])          #E_bar C_bar theta_bar
                C[_] = (1/M)*sum(C[_])
                theta = math.atan2(-1*C[_],-1*E[_]) + math.pi
                COM.append((theta*box_length) / (2*math.pi))
            
            for cord in indy:
                og = points[cord]
                dist = distanceCalc(*COM,*og.cord)**(1/2)
                if dist > box_length/2:
                    new_cords = {0: og.cord[0], 1: og.cord[1], 2: og.cord[2]}
                    for _ in range(0,3):
                        og_dist = ((og.cord[_] - COM[_])**2)**(1/2)
                        dif1 = ((og.unwrap[_][0] - COM[_])**2)**(1/2)      #cord + b
                        dif2 = ((og.unwrap[_][1] - COM[_])**2)**(1/2)      #cord - b
                        if dif1 <= og_dist:
                            new_cords[_] += box_length
                        elif dif2 <= og_dist:
                            new_cords[_] -= box_length
                    
                    gyra_distances.append(distanceCalc(new_cords[0],new_cords[1],new_cords[2], *COM)**(1/2))
    
    return gyra_distances

#########
def network(lst, p, branches=False):
    chains = []
    checked = []
    neighs = []
    cordination = {}
    net = {}
    G = nx.Graph(lst)
    
    for value in p:
        neighbor = list(nx.neighbors(G, value))
        neighbor.append(value)
        neighs.append(neighbor)
        net[value] = len(neighbor)

        path = nx.node_connected_component(G, value)
        chains.append(list(path))
    
    [x.sort() for x in neighs]
    neighs.sort()
    neighs = list(neighs for neighs,_ in itertools.groupby(neighs))
    for i in neighs:
        n = str(len(i))
        try:
            cordination[n] += 1
        except KeyError:
            cordination[n] = 1

    [checked.append(sorted(c, key=lambda x: float(x))) for c in chains]
    checked.sort()
    checked = list(checked for checked,_ in itertools.groupby(checked))
   
    p = list(set([item for sublist in checked for item in sublist]))
   
    if branches == True:
        edges = {}
        terminii = []
        for i in p:
            edges[i] = G.degree(i)
        for item in checked:
            if len(item) > 0:
                edge = []
                for it in item:
                    edge.append(edges.get(it))
                terminii.append([len(item),edge.count(1)])
    else:
        terminii = []
    
    return [checked,cordination,terminii,net]

########
def stepHbondCalc(lst, box_length, origin, bonds, donor_id, acceptor_id, Angles=[False,0], Gyra=False, branches=False, HbondParm=2.5):
    class Cords:
            def __init__(self, particle_id, atom_type, x, y, z, mol_id, ix=0, iy=0, iz=0):
                self.particle_id = particle_id
                self.atom_type = atom_type
                self.cord = [float(x)-origin, float(y)-origin, float(z)-origin]
                self.mol_id = mol_id
                self.unwrap = {0: [self.cord[0]+box_length,self.cord[0]-box_length],
                        1: [self.cord[0]+box_length,self.cord[0]-box_length],
                        2: [self.cord[0]+box_length,self.cord[0]-box_length]}
                self.iflags = [int(ix), int(iy), int(iz)]
                self.mass = getMass(int(atom_type))

    angleParm = Angles[1]
    Angles = Angles[0]

    molecules = {}
    for i in lst:
        molecules[i[0]] = Cords(*i)                       #Dict of [particle_ID]: TrjInfo
    
    np.set_printoptions(threshold = np.inf)
    donors = []
    acceptors = []
    d_lst = []
    a_lst = []
    angles = []
    distances = []
    for i in molecules.values():
        if i.atom_type in donor_id:
            donors.append(i.cord)
            d_lst.append(i.particle_id)
        elif i.atom_type in acceptor_id:
            acceptors.append(i.cord)
            a_lst.append(i.particle_id)

    arra = np.array(distance_matrix(np.array(donors),np.array(acceptors)))
    indices = np.column_stack(np.where((arra >= HbondParm[0]) & (arra <= HbondParm[1])))        #[[x,y],[x,y],...] indicies where arra<HbondParm
    
    output = []
    p = []
    for i in indices:                       #Go through each Cords/particle object in array, 
        og = molecules[d_lst[i[0]]]         #Donor
        nxt = molecules[a_lst[i[1]]]        #Acceptor
        if og.mol_id != nxt.mol_id:
            if Angles == True:
                alcohol_O = molecules[bonds[nxt.particle_id]]
                angle = angleCalc(*og.cord,*nxt.cord,*alcohol_O.cord)
                if angle >= angleParm:
                    distances.append(arra[i[0],i[1]])
                    angles.append(angle)
                    output.append((og.mol_id, nxt.mol_id))
                    p.extend((og.mol_id, nxt.mol_id))
            else:
                distances.append(arra[i[0],i[1]])
                output.append((og.mol_id, nxt.mol_id))
                p.extend((og.mol_id, nxt.mol_id))

    indices = np.column_stack(np.where((arra >= box_length-HbondParm[1]) & (arra <= box_length-HbondParm[0])))         #Possible PBC Hbond

    for i in indices:
        og = molecules[d_lst[i[0]]]
        nxt = molecules[a_lst[i[1]]]
        if og.mol_id != nxt.mol_id:
            new_cords = {0: og.cord[0], 1: og.cord[1], 2: og.cord[2]}
            for _ in range(0,3):
                og_dist = ((og.cord[_] - nxt.cord[_])**2)**(1/2)
                dif1 = ((og.unwrap[_][0] - nxt.cord[_])**2)**(1/2)      #cord + b
                dif2 = ((og.unwrap[_][1] - nxt.cord[_])**2)**(1/2)      #cord - b
                if dif1 <= og_dist:
                    new_cords[_] += box_length
                elif dif2 <= og_dist:
                    new_cords[_] -= box_length
            dist = distanceCalc(new_cords[0],new_cords[1],new_cords[2], *nxt.cord)
            if Angles == True and dist <= HbondParm[1] and dist >= HbondParm[0]:
                alcohol_O = molecules[bonds[nxt.particle_id]]
                angle = angleCalc(new_cords[0],new_cords[1],new_cords[2],*nxt.cord,*alcohol_O.cord)
                #print('unwrap',angle)
                if angle >= angleParm:
                    distances.append(dist)
                    angles.append(angle)
                    output.append((og.mol_id,nxt.mol_id))
                    p.extend((og.mol_id,nxt.mol_id))
            elif Angles == False and dist <= HbondParm[1] and dist >= HbondParm[0]:
                distances.append(dist)
                output.append((og.mol_id,nxt.mol_id))
                p.extend((og.mol_id,nxt.mol_id))

    out = Counter(tuple(sorted(t)) for t in output)     #bifurcated still counted
    contacts = [i for i in out if out[i] >= 2]
    #contacts = [item for sublist in contacts for item in sublist]
    #dimer_type = dimer_typing(lst, contacts)

    n = network(output,list(set(p)),branches)       #list(set(p)) takes out bifurcated

    if branches == True:
        branches = n[2]
    else:
        branches = []

    if Gyra == True:
        gyra_distances = radiiGyra(n[0],molecules,box_length)
    else:
        gyra_distances = []

    return [output,p,n[0],n[1],gyra_distances,branches,n[3],angles,distances]       #Return class data type rather than lists

#########
def trjreader(filename, steps, bonds, donor_id, acceptor_id, XYZ=[False], binning=[False], Gyra=False, time_corr=[False,0], branches=False, ion_pairs=[False], Angles=[False,0], HbondParm=2.5):
    with open(filename) as f:
        #print('hbondParm - min:' + str(HbondParm[0]) + ' max: ' + str(HbondParm[1]))
        n = 9
        st = []         #Initial params
        dic = {}
        chain_length = {}

        Na_cordination = defaultdict(dict)
        Na_len = defaultdict(dict) 
        
        distances = []

        if binning[0] == True:
            num_bins = binning[1]
            chain_bins_avg = {}
            cord_bins_avg = {}
            bin_width = steps/num_bins
            #print('bin_width ' + str(bin_width))
            num_bin_actual = 1
            binning = binning[0]
        if time_corr[0] == True:
            mol_cordination = {}
            num_mols = time_corr[1]
            min_dt = time_corr[2]
            max_dt = time_corr[3]
            skip_dt = time_corr[4]
            t_prop = time_corr[5]
            time_corr = time_corr[0]
        if branches == True:
            bran = {}
        if Angles[0] == True:
            angle = []
        if XYZ[0] == True:
            ag_len = XYZ[1]
            xyz_filename = XYZ[2]
            XYZ = XYZ[0]
        if Gyra == True:
            gyra = []

        init = list(islice(f, n))
        lines = (line.rstrip('r\n') for line in init)
        for data in lines:
            st.append(data.split())
        num_atoms = int(st[3][0])
        #print('Number atoms: ' + str(num_atoms))
        box_length = ((float(st[5][0]) - float(st[5][1]))**2)**(1/2)
        origin = float(st[5][0])
        #print('box length ' + str(box_length)) 
       
        tot_steps = 1
        lst = []
        print('Steps = ' + str(steps))
        while tot_steps < steps:     #All calcs under 
            if tot_steps % 10000 == 0:
                print(tot_steps)
            next_n_lines = list(islice(f, num_atoms+9))
            lines = (line.rstrip('r\n') for line in next_n_lines)
            for data in lines:
                lst.append(data.split())
            del lst[-9::]                                               #Reads in every step
            
            calc = stepHbondCalc(lst, box_length, origin, bonds, donor_id, acceptor_id, Angles, Gyra, branches, HbondParm)  #Donor, acceptor
            
            #StepHbondCalc returns [output,p,n[0],n[1],gyra_distances,branches,n[3],angles,distances]
            distances.append(calc[8])

            if Angles[0] == True:
                angle.append(calc[7])

            ##### Calc {P_interactionNumber:{num_Na's cordinated:freq}} #######
            if ion_pairs[0] == True:
                calcP = stepHbondCalc(lst, box_length, origin, bonds, ion_pairs[1], ion_pairs[2], [False,0], Gyra, branches, ion_pairs[3])   #obj for many input and return
                count = Counter([item for tup in calcP[0] for item in tup])
                for k,v in calc[6].items():
                    try:
                        Na_cordination[v][count[k]] += 1
                    except KeyError:
                        Na_cordination[v][count[k]] = 1

                ##### Calc number of unique Na for chain length #####
                x = []
                Calc = [tuple(calc[2][i]) for i in range(0,len(calc[2]))]                       #calc2=AggMol IDs
                for agg in Calc:
                    ion = False
                    for cord in calcP[0]:
                        if cord[0] in agg:
                            ion = True
                            x.append({agg:cord[1]})
                    if ion == False:
                        x.append({agg:[]})
                res = defaultdict(list)
                for d in x:
                    for k,v in d.items():
                        res[k].append(v)
                for k,v in res.items():
                    res[k] = list(np.unique((np.asarray(v,dtype=object))))
                    try:
                        Na_len[len(k)][len(res[k])] += 1
                    except KeyError:
                        Na_len[len(k)][len(res[k])] = 1
            else:
                calcP = [[],[],[]]
             
            if branches == True:                #Ratio of ends to chain length
                b = calc[5]
                for i in b:
                    d = {i[0]:i[1]}
                    try:
                        bran[i[0]][i[1]] += 1  
                    except KeyError:
                        try:
                            bran[i[0]][i[1]] = 1
                        except KeyError:
                            bran[i[0]] = {i[1]:1}
            
            if Gyra == True:
                gyra.append(calc[4])
            
            if XYZ == True and (calcP[2] != [] or calc[2] != []):
                calc_master = [calcP[2] + calc[2]]
                for i in calc_master:
                    for ag in i:
                        if len(ag) >= ag_len:
                        #flat = [item for sublist in i for item in sublist]
                            XYZwriter(ag, lst, tot_steps, xyz_filename) 
            
            if time_corr == True:
                for key,value in calc[6].items():       #time correlation array
                    try:
                        mol_cordination[key] = np.insert(mol_cordination[key],tot_steps,1)
                    except KeyError:
                        steps = int(steps)
                        mol_cordination[key] = np.zeros((steps,), dtype=int)
                        mol_cordination[key] = np.insert(mol_cordination[key],tot_steps,1)
            
            #Bin average dictionary values
            if binning == True and tot_steps % bin_width == 0 and tot_steps != 0:
                num_bin_actual += 1
                chain_bin_avg = chain_length.copy()
                cord_bin_avg = dic.copy()
                for key, value in chain_bin_avg.items():    #Avg cordinations in one step
                    chain_bin_avg[key] = value/bin_width
                for key, value in cord_bin_avg.items():
                    cord_bin_avg[key] = value/bin_width
                
                for key, value in chain_bin_avg.items():    #Add to running counter/total  average
                    try:
                        chain_bins_avg[key] += value
                    except KeyError:
                        chain_bins_avg[key] = value
                for key, value in cord_bin_avg.items():
                    try:
                        cord_bins_avg[key] += value
                    except KeyError:
                        cord_bins_avg[key] = value

            for item in calc[2]:
                length = str(len(item))
                try:
                    chain_length[length] += 1
                except KeyError:
                    chain_length[length] = 1
            for item in calc[3].items():
                try:
                    dic[item[0]] += int(item[1])
                except KeyError:
                    dic[item[0]] = int(item[1])
           
           ###### 
            lst = []
            tot_steps += 1
           ######Next Snapshot

        output = {}

        dists = [item for sublist in distances for item in sublist]
        plt.hist(dists, bins=100)
        fontsize=16
        plt.title('Hydrogen bond Distances',fontsize=fontsize)
        plt.xlabel('Distance (Angstroms)',fontsize=fontsize)
        plt.ylabel('Counts',fontsize=fontsize)
        plt.savefig('Dist.png',dpi=300)
        plt.clf()

        if ion_pairs[0] == True:
            Na_cordination = avgNested(Na_cordination)
            Na_len = avgNested(Na_len)
            #print('len of Agg: #Na: Freq')
            #print(Na_len)
            #print('Cordination: #Na: Freq')
            #print(Na_cordination)
            output['Ion_len'] = Na_len
            output['Ion_cordination'] = Na_cordination

        if Angles[0] == True:
            angles = [item for sublist in angle for item in sublist]
            plt.hist(angles, bins=100)
            plt.title('D-H-A Hydrogen bond Angle',fontsize=fontsize)
            plt.xlabel('Angle (degrees)',fontsize=fontsize)
            plt.ylabel('Counts',fontsize=fontsize)
            plt.savefig('Angles.png',dpi=300)
            plt.clf()

        if branches == True:
            for key,value in bran.items():
                tot = 0
                for v in value.values():
                    tot += v
                for k,v in value.items():
                    value[k] = (v/tot)*100 
                    #print('Chain Length: {Edges: Ratio}')
                #print(bran)
            output['branches'] = bran

        if time_corr == True:
            mol_counter = 0
            Time_Corr = defaultdict(dict)
            for k,v in mol_cordination.items():
                #print('Next Mol')
                mol_counter += 1
                for dt in range(min_dt,max_dt,skip_dt):       #range(1,int(tot_steps/2)+1,10)
                    t = 0
                    l = np.array([],dtype=int)
                    while t < tot_steps:
                        a = v[t:t+dt]
                        A = np.prod(a)
                        l = np.append(l,A)
                        t += t_prop                 #Spacer propagation
                    if dt != 1:               #Normalization
                        l = (np.sum(l) / len(l)) / Time_Corr[k][1]
                        Time_Corr[k][dt] = round(l, 3)
                    else:
                        l = np.sum(l) / len(l)
                        Time_Corr[k][dt] = l

                    Time_Corr[k][1] = 1            #Normalization

            data = {}
            for v in Time_Corr.values():
                #print(v)
                for key,value in v.items():
                    try:
                        data[key] += value/num_mols
                    except KeyError:
                        data[key] = value/num_mols
            output['time_corr'] = data
            #print(str(mol_counter) + ' molecules time correlated')

        if Gyra == True:
            radii = [item for sublist in gyra for item in sublist]
            radii = np.array(radii)
            plt.hist(radii, bins=20)
            plt.title('Distribution of Radius of Gyration')
            plt.xlabel('Radius   (angstroms)')
            plt.ylabel('Counts')
            plt.savefig("radii.png")
            plt.clf()

        if binning == True:
            #print('Number of bins ' + str(num_bin_actual))
            chain_std = {}
            cord_std = {}
            for key, value in chain_bins_avg.items():
                try:
                    chain_std[key] = value/(num_bin_actual**(1/2))
                except KeyError:
                    chain_std[key] = value/(num_bin_actual**(1/2))
            for key, value in cord_bins_avg.items():
                try:
                    cord_std[key] = value/(num_bin_actual**(1/2))
                except KeyError:
                    cord_std[key] = value/(num_bin_actual**(1/2))

            #print('Chain Length')
            #print(chain_length)
            #print('Cordinations')
            #print(dic)
        output['Chain Length'] = chain_length
        output['Cordinations'] = dic

        if binning == True:
            #print('Chain Length std')
            #print(chain_std)
            #print('Cordination Std')
            #print(cord_std)
            output['Chain Length Std'] = chain_std
            output['Cordination Std'] = cord_std
        
        json_f = open('output.json','a')
        #json.dump(output,json_f,indent=len(output.keys()))
        json_string = json.dumps(output)
        json_f.write('\n')
        json_f.write(json_string)
        json_f.close()
        
        #try:
        #    current_data = json.load(json_f)
        #    for k,v in output.items():
        #        x = current_data[k]
        #        y = output[k]
        #        x = Counter(x)
        #        y = Counter(y)
        #        z = dict(Counter(x)+Counter(y))
        #        current_data[k] = z
        #    json.dump(current_data,json_f,indent=len(current_data.keys()))
                        
        #except json.decoder.JSONDecodeError:    #Except empty file
        #    json_f.close()
        #    json_f = open('output.json','w')
        #    json.dump(output,json_f,indent=len(output.keys()))
        #    json_f.close()

        return [dic,chain_length]   # Simple list of IN number in each step

#########


