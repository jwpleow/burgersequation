CXX = g++
CXXFLAGS = -std=c++14 -Wall -O3
HDRS = Model.h Burgers.h
OBJS = burgers.o Burgers.o Model.o
LDLIBS = -lblas

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

# Test cases: Arguments are given in the order: Lx Ly T Nx Ny Nt ax ay b c
diff: compile
	./myProg 10 10 1 2001 2001 4000 0 0 0 1

advx: compile
	./myProg 10 10 1 501 501 4000 1 0 0 0

burg2131: compile
	./myProg 10 10 1 21 31 4001 1 0.5 1 0.02

burg501: compile
	./myProg 10 10 1 21 31 4000 1 0.5 1 0.02

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
