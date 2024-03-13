# **Install**

   	pip install aggPy

# **Make Analysis object**

	json file has the parameters for the desired analysis
	x = Analysis('file.json', 'hbond')

# **Run analysis**

	x.hCoordination()
 	x variable has desired analysis attributes
	Possible key values = 'Distance', 'Angle', 'Coordinations' ,'Aggregate Size', 'Aggregate resids', 'Network', 'Node Degree'

# **Workup** -
   ### Totals values of a key property - return: list of total key values
	i = Analysis.aggregate(x, 'key')	- Totals values of a key property - return: list of total key values 

   ### Average value of key each ts - returns same dtype as 'key' dtype
	j = Analysis.average(x, 'key')	- Average value of key each ts - returns same dtype as 'key' dtype

   ### std_dev from binning avg - returns same dtype as 'key' dtype
	k = Analysis.std_dev(x, 'key', bin_width=1)	- std_dev from binning avg - returns same dtype as 'key' dtype

   ### Time Correlation - returns: mol_ct, sys_ct - _mol_ct=per molecule Ct , sys_ct=total Ct
	y = Analysis.timeCorr(x) , x.timeCorr()	 - Time Correlation - returns: mol_ct, sys_ct - _mol_ct=per molecule Ct , sys_ct=total Ct	
	   mol_ct, sys_ct = x.timeCorr()
   
