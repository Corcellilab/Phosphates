import ast
from multiprocessing import Process
from multiprocessing import Queue
from filesplit.split import Split
from filesplit.merge import Merge
from HbondCore import *

def cleanup():
    with open('output.json') as json_f:
        count = 0
        for line in json_f.readlines():
            if line != '\n':
                data = ast.literal_eval(line)
                if count == 0:
                    output = data
                    count += 1
                elif count != 0:
                    for key,value in data.items():
                        for k,v in value.items():
                            try:
                                output[key][k] += data[key][k]
                            except KeyError:
                                output[key][k] = data[key][k]
#        print(output)
    with open('output.json','w') as json_f:
        json.dump(output,json_f,indent=len(output.keys()))

def multiProcs(filename,tot_step,bonds,num_procs,trjargs):
    f = open('output.json','w')
    f.close()

    with open(filename) as input_file:
        data = [input_file.readline().split() for line in range(0,4)]
        num_atoms = int(data[-1][0])
    
    split = Split(filename,'./')

    trunc_steps = tot_step % num_procs
    print('Removing last ' + str(trunc_steps) + ' steps')
    tot_step -= trunc_steps

    steps_per_proc = tot_step/num_procs

    num_lines_perstep = num_atoms + 9 
    num_lines_perproc = num_lines_perstep * (steps_per_proc)

    print('File splitting')
    split.splitdelimiter='.'
    split.bylinecount(num_lines_perproc)
    print('Done splitting')

    print('Rename')
    file_lst = []
    for i in range(1,num_procs+1):
        f = filename.split('.')
        f.insert(2, str(i))
        f = '.'.join(f)
        file_lst.append(f)
    print('Done rename')

    #print('sending to multiple procs')
    if __name__ == "__main__":
        print('sending to multiple procs')
        proc_lst = []
        print('Steps per proc = ' + str(steps_per_proc))
        for i in range(0,num_procs):
            x = Process(target=trjreader, args=(file_lst[i], steps_per_proc, bonds, *trjargs))
            proc_lst.append(x)
        print('Starting jobs')
        for i in proc_lst:
            i.start()
        print('Joining jobs')
        for i in proc_lst:
            i.join()
    
    print('Merging')
    merge = Merge('./','./','MergedTRJ.lammpstrj')
    merge.merge(cleanup=True)
    print('Done Merging')

    return None
#######

XYZ = [False, 2, 'Hbond1.xyz']      #Min aggregate size, filename
binning = [False, 100]               #Number of bins
Gyra = False
time_corr = [False, 192, 1, 15000, 10, 100]      # -, number of molecules, min_spacer_size, max_space_size, spacer_size_interval, spacer_jump
branches = False
ion_pairs = [False,'1','5',[3.15,4.10]]        # -, atom_type1, atom_type2, Distance Check
Angles=[True,150]                            # -, max Hbond Angle
HbondParm=[0,2.5]
acceptor_id = ['2','3']        #Lone pair atom
hyd_id = '4'            #Hydrogen atom

tot_step = 150000

data_file = 'final.Conc1M.data'
trj_file = 'phos.Conc1M.lammpstrj'
num_procs = 4

trjargs = [acceptor_id,hyd_id,XYZ,binning,Gyra,time_corr,branches,ion_pairs,Angles,HbondParm]

bonds = datreader(data_file)              #[{particleID:particleID,pID:pID,...}, typeMap]
bonds = bonds[0]
multiProcs(trj_file,tot_step,bonds,num_procs,trjargs)

cleanup()


