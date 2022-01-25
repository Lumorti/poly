CXX=g++
CXXFLAGS=-Wfatal-errors -O3
LIBS=-I/usr/include/eigen3

all:
	$(CXX) $(CXXFLAGS) -o solver solver.cpp $(LIBS)
	$(CXX) $(CXXFLAGS) -o gen gen.cpp $(LIBS)

solver:
	$(CXX) $(CXXFLAGS) -o solver solver.cpp $(LIBS)

gen:
	$(CXX) $(CXXFLAGS) -o gen gen.cpp $(LIBS)
