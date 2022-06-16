CXX=g++
CXXFLAGS=-Wfatal-errors -O3
LIBS=-I/usr/include/eigen3

all: integral
	
integral:
	$(CXX) $(CXXFLAGS) -o integral integral.cpp $(LIBS)
