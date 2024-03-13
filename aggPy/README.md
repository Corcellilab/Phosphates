0) #**Install**
   	pip install aggPy

2) #**Make Analysis object**
	_json file has the parameters for the desired analysis_
	x = Analysis('file.json', 'hbond')

3) #**Run analysis**
   	_x variable now has more attributes_
	x.hCoordination()

5) #**Workup** -
   	_Possible key values = 'Distance', 'Angle', 'Coordinations' ,'Aggregate Size', 'Aggregate resids', 'Network', 'Node Degree'_
	
	x.aggregate('key')	- Totals values of a key property - return: list of total key values 
				  	
	x.average('key')	- Average value of key each ts - adds x.{key}Avg attribute

	x.std_dev('key', bin_width=1)	- std_dev from binning avg - add x.{key}Stdev attribute

	x.timeCorr()		- Time Correlation - returns: mol_ct, sys_ct - mol_ct=per molecule Ct , sys_ct=total Ct	
	    mol_ct, sys_ct = x.timeCorr()
   
