2D Parallelised Burgers' Equation Solver using C++ and MPI (with time integration done by a Explicit Forward Euler integration scheme).

A Makefile is provided with a few different configurations - the parameters of the problem and number of processors to use can be easily changed here.
The initial velocity field can be changed relatively easily in SetVelField in Burgers.cpp (note that the data is stored and processed in column-major format)
Boundary conditions are set as u, v = 0 at the edges of the domain and cannot be easily changed.

Parameters are stored in the Model class (Model.cpp and Model.h), with data processing done in the Burgers class (Burgers.cpp and Burgers.h), and the functions are called in burgers.cpp.

 




