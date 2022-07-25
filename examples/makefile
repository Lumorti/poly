CXX=g++
CXXFLAGS=-Wfatal-errors -O3
LIBS=-I/usr/include/eigen3 -fopenmp

all: mub simple matrices

mub:
	$(CXX) $(CXXFLAGS) -o mub mub.cpp $(LIBS)

simple:
	$(CXX) $(CXXFLAGS) -o simple simple.cpp $(LIBS)

matrices:
	$(CXX) $(CXXFLAGS) -o matrices matrices.cpp $(LIBS)

