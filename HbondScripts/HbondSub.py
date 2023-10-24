from HbondMulti import *

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

