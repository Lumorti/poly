CXX=g++
CXXFLAGS=-Wfatal-errors -O3
LIBS=-I/usr/include/eigen3

all: mub

mub:
	$(CXX) $(CXXFLAGS) -o mub mub.cpp poly.h $(LIBS)

factor:
	$(CXX) $(CXXFLAGS) -o factor factor.cpp poly.h $(LIBS)

