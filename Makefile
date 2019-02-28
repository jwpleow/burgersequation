CXX = mpicxx
CXXFLAGS = -std=c++14 -Wall -O3

default: burg

all: diff advx advx burg

burgers.o: burgers.cpp
	$(CXX) $(CXXFLAGS) -o burgers.o -c burgers.cpp

Burgers.o: Burgers.cpp Burgers.h
	$(CXX) $(CXXFLAGS) -o Burgers.o -c Burgers.cpp

Model.o: Model.cpp Model.h
	$(CXX) $(CXXFLAGS) -o Model.o -c Model.cpp

compile: burgers.o Burgers.o Model.o
	$(CXX) $(CXXFLAGS) -o myProg burgers.o Burgers.o Model.o -lblas

# Test cases: Arguments are given in the order (after ./myProg): Lx Ly T Nx Ny Nt ax ay b c Px Py
diff: compile
	mpiexec -np 2 ./myProg 10 10 1 101 201 4000 0 0 0 1

advx: compile
	mpiexec -np 2 ./myProg 10 10 1 101 201 4000 1 0 0 0

advy: compile
	mpiexec -np 2 ./myProg 10 10 1 101 201 4000 0 1 0 0

burg: compile
	mpiexec -np 12 ./myProg 10 10 1 501 501 4000 1 0.5 1 0.02 4 3

# Rule to clean the source directory
.PHONY: clean
	
clean:
	-rm -f *.o myProg





#test each section for O2, O3 and Os speeds
#burgers.o: burgers.cpp
#    g++ -std=c++14 -Wall -O2 -o burgers.o -c burgers.cpp

#Burgers.o: Burgers.cpp Burgers.h
#    g++ -std=c++14 -Wall -O2 -o Burgers.o -c Burgers.cpp

#Model.o: Model.cpp Model.h
#    g++ -std=c++14 -Wall -O2 -o
