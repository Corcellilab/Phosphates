# Use of aggPy

## **Install**

   	pip install aggPy

## **Make Analysis object**

	x = Analysis('file.json', 'hbond')
 	
  	json file has the parameters for the desired analysis

## **Run analysis**

	x.hCoordination()
 
 	x variable has desired analysis attributes
	Possible key values = 'Distance', 'Angle', 'Coordinations' ,'Aggregate Size', 'Aggregate resids', 'Network', 'Node Degree'

## **Workup**
   ### Totals values of a key property - return: list of total key values
	i = Analysis.aggregate(x, 'key')	 

   ### Average value of key each ts - returns same dtype as 'key' dtype
	j = Analysis.average(x, 'key')	

   ### std_dev from binning avg - returns same dtype as 'key' dtype
	k = Analysis.std_dev(x, 'key', bin_width=1)	

   ### Time Correlation - returns: mol_ct, sys_ct - mol_ct=per molecule Ct , sys_ct=total Ct
	y = Analysis.timeCorr(x) 	 	
	mol_ct, sys_ct = x.timeCorr()
   
