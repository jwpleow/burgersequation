# 2D Burgers' Equation Solver using C++ and MPI 
- time integration done by a Explicit Forward Euler integration scheme
- can use however many processes you want




## Using the program
A Makefile is provided with a few different configurations - the parameters of the problem and number of processors to use can be easily changed here.
To use the Makefile to execute the program:
1. Place the Makefile and BurgersProg in the same directory
2. In a terminal in that directory, type:
```
make burgpn
```
Windows users would require a terminal that can read Makefiles e.g. NMake.

A few other cases are provided in the Makefile (e.g. Diffusion can be run with make diffpn)


The final velocity field is printed in VelocityFields.txt


The initial velocity field can be changed relatively easily in SetVelField in Burgers.cpp (note that the data is stored and processed in column-major format).
Boundary conditions are set as u, v = 0 at the edges of the domain and cannot be easily changed.

Parameters are stored in the Model class (Model.cpp and Model.h), with data processing done in the Burgers class (Burgers.cpp and Burgers.h), and the functions are called in burgers.cpp.

 




