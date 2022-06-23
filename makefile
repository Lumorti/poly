CXX=g++
CXXFLAGS=-Wfatal-errors -O3
LIBS=-I/usr/include/eigen3 -fopenmp

all: mub

mub:
	$(CXX) $(CXXFLAGS) -o mub mub.cpp $(LIBS)


