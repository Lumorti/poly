CXX=mpic++
CXXFLAGS=-Wfatal-errors -O3
LIBS=-I/usr/include/eigen3

all:
	$(CXX) $(CXXFLAGS) -o solve solver.cpp $(LIBS)
