0) install
   	pip install aggPy

2) make analysis object
	json file has the parameters for the desired analysis
	x = Analysis('file.json', 'hbond')

3) Run analysis - x variable now has more attributes
	x.hCoordination()

4) Workup - Possible key values = 'Distance', 'Angle', 'Coordinations' ,'Aggregate Size', 'Aggregate resids', 'Network', 'Node Degree'
	
	x.aggregate('key')	- Totals values of a key property - return: list of total key values 
				  	
	x.average('key')	- Average value of key each ts - adds x.{key}Avg attribute

	x.std_dev('key', bin_width=1)	- std_dev from binning avg - add x.{key}Stdev attribute

	x.timeCorr()		- Time Correlation - returns: mol_ct, sys_ct - mol_ct=per molecule Ct , sys_ct=total Ct	
	    mol_ct, sys_ct = x.timeCorr()
   
