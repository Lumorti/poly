CXX=g++
CXXFLAGS=-Wfatal-errors -O3 -march=native
DEBUGFLAGS=-Wall -g -pg -Wfatal-errors
LIBSEIGEN= -I/usr/include/eigen3
LIBSOMP= -fopenmp
LIBSSCS= -L${SCSHOME}/lib/ -I${SCSHOME}/include/scs/ -lscsdir
LIBSMOSEK= -I${MSKHOME}/h -L${MSKHOME}/bin -Wl,-rpath-link,${MSKHOME}/bin -Wl,-rpath=${MSKHOME}/bin -lmosek64 -lfusion64
LIBSOPTIM= -I${OPTIMHOME}/include/ -L${OPTIMHOME} -loptim
LIBS=$(LIBSEIGEN) $(LIBSOPTIM) $(LIBSMOSEK) $(LIBSSCS) $(LIBSOMP)

all: mub simple sos ising

mubdebug:
	$(CXX) $(DEBUGFLAGS) -o mub mub.cpp $(LIBS) 

mub:
	$(CXX) $(CXXFLAGS) -o mub mub.cpp $(LIBS) 

simple:
	$(CXX) $(CXXFLAGS) -o simple simple.cpp $(LIBS) 

matrices:
	$(CXX) $(CXXFLAGS) -o matrices matrices.cpp $(LIBS) 

sos:
	$(CXX) $(CXXFLAGS) -o sos sos.cpp $(LIBS) 

ising:
	$(CXX) $(CXXFLAGS) -o ising ising.cpp $(LIBS) 


