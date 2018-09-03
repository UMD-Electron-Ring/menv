# menv

Envelope integration and matching code written by Hui Li. Based off SPOT by Chris Allen 

# Quick-Start Guide 
Download MENV directories, add to Matlab path (run setpath from Matlab command line)
To launch GUI: run "menv" from command line. Can assemble lattice elements from within GUI, or load .spt file.
To run from command line: clmenv() creates a class object with all MENV methods. See example file run_clmenv.m 

More information is in documentation/menvdoc.pdf



# UPDATES

8/29/2018
Removed I/O operations for temporary files during solves (runtmp.mat, params.mat, etc) resulting in SIGNIFICANT speed-up. 

6/1/2018
An updated version has now been merged with the master branch. Updated features include: 
* scriptable command-line version
* now includes dispersion calculations based on bending radius in dipoles
* A reference trajectory can be defined and used as a target for optimization routines
* Global search and multi-objective optimization algorithms are now included in addition to the nonlinear-least-squares fitting routine previously used
* Tune is now included as a target condition during optimization

