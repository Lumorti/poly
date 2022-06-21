CXX=g++
CXXFLAGS=-Wfatal-errors -O3
LIBS=-I/usr/include/eigen3

all: mub

mub:
	$(CXX) $(CXXFLAGS) -o mub mub.cpp $(LIBS)

factor:
	$(CXX) $(CXXFLAGS) -o factor factor.cpp $(LIBS)

