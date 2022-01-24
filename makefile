CXX=g++
CXXFLAGS=-Wfatal-errors -O3
LIBS=-I/usr/include/eigen3

all: solver gen

solver:
	$(CXX) $(CXXFLAGS) -o solve solver.cpp $(LIBS)

gen:
	$(CXX) $(CXXFLAGS) -o gen gen.cpp $(LIBS)
