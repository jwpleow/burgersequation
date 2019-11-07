# Parallel 2D Burgers' Equation Solver using C++ and MPI 
- time integration done by an FTCS scheme
- infinitely scalable (I think...)
- Serial (not-parallel) version is included in the other branch!



## Using the program
A Makefile is provided with a few different configurations - the parameters of the problem and number of processors to use can be changed here.
To use the Makefile to execute the program (e.g. default is set to compile and execute the Burgers' Equation solver for the parameters shown in the makefile):
1. Place the Makefile and BurgersProg in the same directory
2. In a terminal in that directory, type:
```
make
```

A few other cases are provided in the Makefile (e.g. Diffusion can be run with make diff)

The final velocity field is printed in VelocityFields.txt


The initial velocity field can be changed relatively easily in SetVelField in Burgers.cpp (note that the data is stored and processed in column-major format).
Boundary conditions are set as u, v = 0 at the edges of the domain.

Parameters of the problem are stored in the Model class (Model.cpp and Model.h), with data processing done in the Burgers class (Burgers.cpp and Burgers.h), and the functions called in burgers.cpp.

 
