Example.json file format w/ keywords
------------------------------------
{
  "init": [
    {
        "data_file":            "../postNPT.data",
        "trj_file":             "../modFF.dcd",		
        "md_form":              "LAMMPS",		- Md analysis format
        "dims":                 "",			- Box dimensions
        "mol_deliminator":      "resid",		- Md analysis selection deliminator
        "run_length":           "n",			- y/n , y = run for entire trj, n = windowed analysis
        "start_frame":          0,			
        "end_frame":            50,
        "frame_dt":             1			- jump between frames
        }
  ],
  "hbond": [
    {
        "HbondParm":            "0, 2.13, 2.61" ,	- minimum distance, maximum distance, angle cutoff in radians
        "hydrogen_type":        "type 4",		- Md analysis hydrogen selection
        "donors_type":          "type 2",		- Md analysis donors selection
        "acceptors_type":       "type 3",		- Md analysis acceptors selection
        "merge_a_d":            "y",			- Combine the donor list into the acceptor list
        "XYZ":                  "[False]",		- Output an xyz file of aggregates at each step [True/False, atom_map, XYZfilename]
        "mol_bool":             "y",			- Guess molecule delimination, y/n
        "mol_deliminator":	"resid",		- Md analysis molecule selector
	"dist_bool":            "y",			- Output distances to out file
        "angle_bool":           "y",			- Output angles to out fle
        "output_filename":      "testing"		
        }
  ]
}
------------------------------------

1) make analysis object
	json file has the parameters for the desired analysis
	x = Analysis('file.json', 'hbond')

2) Run analysis - x variable now has more attributes
	x.hCoordination()

3) Workup - Possible key values = 'Distance', 'Angle', 'Coordinations' ,'Aggregate Size', 'Aggregate resids', 'Network', 'Node Degree'
	
	x.aggregate('key')	- Totals values of a key property - return: list of total key values 
				  	
	x.average('key')	- Average value of key each ts - adds x.{key}Avg attribute

	x.std_dev('key', bin_width=1)	- std_dev from binning avg - add x.{key}Stdev attribute

	x.timeCorr()		- Time Correlation - returns: mol_ct, sys_ct - mol_ct=per molecule Ct , sys_ct=total Ct	
	    mol_ct, sys_ct = x.timeCorr()


