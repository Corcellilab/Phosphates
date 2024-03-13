# **Install**
   	pip install aggPy

# **Make Analysis object**
	_json file has the parameters for the desired analysis_
	x = Analysis('file.json', 'hbond')

# **Run analysis**
   	_x variable now has more attributes_
	x.hCoordination()

# **Workup** -
   	_Possible key values = 'Distance', 'Angle', 'Coordinations' ,'Aggregate Size', 'Aggregate resids', 'Network', 'Node Degree'_
	
	i = Analysis.aggregate(x, 'key')	- Totals values of a key property - return: list of total key values 
				  	
	j = Analysis.average(x, 'key')	- Average value of key each ts - returns same dtype as 'key' dtype

	k = Analysis.std_dev(x, 'key', bin_width=1)	- std_dev from binning avg - returns same dtype as 'key' dtype

	y = Analysis.timeCorr(x) , x.timeCorr()	 - Time Correlation - returns: mol_ct, sys_ct - _mol_ct=per molecule Ct , sys_ct=total Ct	
	   mol_ct, sys_ct = x.timeCorr()
   
