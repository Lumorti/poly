CXX=g++
CXXFLAGS=-Wfatal-errors -O3
LIBS=-I/usr/include/eigen3

all: solver gen integral

solver:
	$(CXX) $(CXXFLAGS) -o solver solver.cpp $(LIBS)

gen:
	$(CXX) $(CXXFLAGS) -o gen gen.cpp $(LIBS)
	
integral:
	$(CXX) $(CXXFLAGS) -o integral integral.cpp $(LIBS)
