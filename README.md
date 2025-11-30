# LP_for_PL
C code and Windows binary for the LP method described in the Ms "An approach based on Linear Programming to build experimentally-driven Pump-Leak models" by Luigi Catacuzzeno, Maurizio Gildo Cavaliere and Antonio Michelucci  "


		Linear Programming Method
Software to find unknown parameters in a system of differential equation, and then numerically solve the system.

0 - download all the files and unzip in a folder

1 - in the folder you will find 
	a. a text file named "model.txt"  that is used to input information on the L-P model
	b. an executable named "LPmethod.exe" that process "model.txt" and gives back the result
	c. all the other files represent the source code used to compile "LPmethod.exe". 
NB   	Compilation of "LPmethod.exe" uses the glpk library that needs to be included in your compiler

2 - Modify the "model.txt" file to fit your problem and save it. Inside this file you will find 
		comment lines preceeded by "#" giving you detailed instructions

3 - double click on the "LPmethod.exe" file

4 - a text file "results.txt" should appear reporting the results of the LP solver

5 - If your system has no unknown parameters, a text file "output.txt" should additionally appear reporting the solution of the system of differential equations (the first column is time and in the following columns are the variables in your system)

Good luck!  
