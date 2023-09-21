				Organoid Patterning aLgorithm (OPaL)	

	Author: Abdel Rahman Abdel Fattah
	email: dita.roushdy@gmail.com/aabdelfattah@cemm.at

       
       	Readme (OPaL_One): 
      	       (OPaL_One) returns average source cell (SC) expression outputs from 
       	        a fixed set of condition inputs by using the OPaL_Single_Run function. 
	        Please keep both files (OPaL_One and OPaL_Single_Run) in the same directory.
		Input parameters in OPaL_One, Output parameters can be extracted from the MATLAB workspace

       	Inputs (OPaL_One):
           size:       The number of elements (cells) in the domain (currently set at 46 cells)
           nuc_wid:    The width of a cell's nucleus, this is used to
                       calculate domain lenght (currently set at 6 um)
           time:       The number of time steps (currently set at 10 steps)
           iterations: The number of simulated organoids, this is the
                       number that will be used for averaging (analogous to number of
                       observed organoids in a culture)
           threshold:  An activation signal threshold that defines whether
                       cell remains, becomes or loses the SC identity. The higher the
                       value the more difficult it is for a cell to become an SC
                       (currnetly set at 2.36)
           induction:  Defines the fraction of initial SCs in a domain
                       (currently set at 0.5)
           ca:         Activation signal decay constant (currently set at 0.05)
           a:          Activation signal amplitude (currently set at 1)

`      
	Outputs(OPaL_One):
           Exp_scattering_organoids: Average scattering expression fraction
           Exp_patterning_organoids: Average patterning expression fraction
           Exp_negative_organoids: Average negative expression fraction
           Table_of_outputs: expression, SC output fraction, number of poles 
                             and position of poles for each in-silico organoid 
                             in the iteration
               Column 1: 0 = Scattered, 1 = Patterned, -1 = Negative
               Column 2: SC output fraction (between 0 and 1)
               Column 3: Number of poles detected (scattered organoids
                         still detect a pole)
              Columns 4 to 10: pole position (degs) for every pole
                               detected in the organoid

_______________________________________________________________________________________________________________________
	
	Readme (OPaL_Multiple): 
      	       (OPaL_Multiple) returns 3 matrices for each of the expression phenotype
		Scattered, Patterned and Negative. The matrices map the parameter space as 
		defined by the user input of threshold variation and exponential decay 
		constant variation and used the function defined in OPaL_Multiple_Run.
		Please keep both files (OPaL_Multiple and OPaL_Multiple_Run) in the same directory.
		Input parameters in OPaL_Multiple, Output parameters can be extracted from the MATLAB workspace

	Inputs (OPaL_Multiple):
           size:          The number of elements (cells) in the domain (currently set at 46 cells)
           nuc_wid:       The width of a cell's nucleus, this is used to
                          calculate domain lenght (currently set at 6 um)
           time:          The number of time steps (currently set at 10 steps)
           iterations: The number of simulated organoids, this is the
                       number that will be used for averaging (analogous to number of
                       observed organoids in a culture)
           induction:     Defines the fraction of initial SCs in a domain
                          (currently set at 0.5))
           a:             Activation signal amplitude (currently set at 1)	
	   threshold_ini: Initial threshold value
	   threshold_end: Final threshold value
	   threshold_inc: Threshold increment value
	   ca_ini: Initial exponential decay value
	   ca_end: Final exponential decay value
	   ca_inc: exponential decay increment value

	Outputs(OPaL_Multiple):
           result_matrix_scattering: Average scattering expression fraction for every tested condition
           result_matrix_patterning: Average patterning expression fraction for every tested condition
           result_matrix_negative: Average negative expression fraction for every tested condition
	   NOTE: The rows mark the threshold range while the columns mark the exponential decay constant range
