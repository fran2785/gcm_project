This is the code used to structure the Monte Carlo simulation of a set of N particles in a box with periodic boundary conditions in the framework of the Gaussian Core model.
The simulation has its core in Fortran 90 where
- 'gcm_module.f90' sets the useful functions for the model
- 'montecarlo_module.f90' is a general code to perform the basic calculations, like pbc setting, energy calculations and the MC step
- 'gr_module.f90' is a readapted code from Allen and Tildesley and performs the rdf calculation

The simulation is performed in Python by the Fortran functions with the f2py interface.
In Python there are three communicating classes: one for the initialization of the system, one for the actual simulation and the last one for energy and rdf calculation.
The initial configuration in this case is taken from a file where there are the number of particles, the lenght of the simulation box based on the input density and the (N,3) array of initial positions. 
