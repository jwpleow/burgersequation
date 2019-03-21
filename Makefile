CXX = mpicxx
CXXFLAGS = -std=c++14 -Wall -O2 -ffast-math -march=native
# warning: -ffast-math does approximations for better computation speed!
default: burg

burgers.o: burgers.cpp
	$(CXX) $(CXXFLAGS) -o burgers.o -c burgers.cpp

Burgers.o: Burgers.cpp Burgers.h
	$(CXX) $(CXXFLAGS) -o Burgers.o -c Burgers.cpp

Model.o: Model.cpp Model.h
	$(CXX) $(CXXFLAGS) -o Model.o -c Model.cpp

compile: burgers.o Burgers.o Model.o
	$(CXX) $(CXXFLAGS) -o BurgersProg burgers.o Burgers.o Model.o -lblas

# Test cases: Arguments are given in the order (after ./BurgersProg): Lx Ly T Nx Ny Nt ax ay b c Px Py 
# Note - must match number of processors (the number after -np) with Px * Py
diff: compile
	mpiexec -np 16 ./BurgersProg 10 10 1 1001 1001 40000 0 0 0 1 4 4
#note diffusion has convergence issues if spatial discretisation is too fine!
advx: compile
	mpiexec -np 16 ./BurgersProg 10 10 1 2001 2001 4000 1 0 0 0 4 4

advy: compile
	mpiexec -np 16 ./BurgersProg 10 10 1 2001 2001 4000 0 1 0 0 4 4

burg: compile
	mpiexec -np 16 ./BurgersProg 10 10 1 2001 2001 4000 1 0.5 1 0.02 4 4

# Rule to clean the source directory
.PHONY: clean
clean:
	-rm -f *.o BurgersProg
