# Use of aggPy

## **Install**

   	pip install aggPy

## **Make Analysis object**

	import aggPy
 	hbonds = Analysis('file.json', 'hbond')
 	
  	json file has the parameters for the desired analysis

## **Run analysis**

	hbonds.hCoordination()
 
 	hbonds variable has desired analysis attributes
	Possible key values = 'Distance', 'Angle', 'Coordinations' ,'Aggregate Size', 'Aggregate resids', 'Network', 'Node Degree'

## **Workup**
   ### Totals values of a key property - return: list of total key values
	i = Analysis.aggregate(hbonds, 'key')	 

   ### Average value of key each ts - returns same dtype as 'key' dtype
	j = Analysis.average(hbonds, 'key')	

   ### std_dev from binning avg - returns same dtype as 'key' dtype
	k = Analysis.std_dev(hbonds, 'key', bin_width=1)	

   ### Time Correlation - returns: x, y
	x, y = Analysis.timeCorr(hbonds) 
 		x = number of timesteps until hbond between molecules no longer exists 
   		y = Normalized counts
