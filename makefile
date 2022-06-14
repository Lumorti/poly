CXX=g++
CXXFLAGS=-Wfatal-errors -O3
LIBS=-I/usr/include/eigen3 -I/usr/include/lbfgspp

all: integral
	
integral:
	$(CXX) $(CXXFLAGS) -o integral integral.cpp $(LIBS)

integral2:
	$(CXX) $(CXXFLAGS) -o integral integral2.cpp $(LIBS)
