CXX=g++
CXXFLAGS=-Wfatal-errors -O3 -march=native
DEBUGFLAGS=-Wall -g -pg -Wfatal-errors
LIBSEIGEN= -I/usr/include/eigen3
LIBSOMP= -fopenmp
LIBSSCS= -L${SCSHOME}/lib/ -I${SCSHOME}/include/scs/ -lscsdir
LIBSMOSEK= -I${MSKHOME}/h -L${MSKHOME}/bin -Wl,-rpath-link,${MSKHOME}/bin -Wl,-rpath=${MSKHOME}/bin -lmosek64 -lfusion64
LIBSOPTIM= -I/usr/local/include/optim/ -loptim
LIBS=$(LIBSEIGEN) $(LIBSOPTIM) $(LIBSMOSEK) $(LIBSSCS) $(LIBSOMP)

all: mub simple matrices sos ising

mubdebug:
	$(CXX) $(DEBUGFLAGS) -o mub mub.cpp $(LIBS) 

mub:
	$(CXX) $(CXXFLAGS) -o mub mub.cpp $(LIBS) 

mubComplex:
	$(CXX) $(CXXFLAGS) -o mubComplex mubComplex.cpp $(LIBS) 

testing:
	$(CXX) $(CXXFLAGS) -o testing testing.cpp $(LIBS) 

testingdebug:
	$(CXX) $(DEBUGFLAGS) -o testing testing.cpp $(LIBS) 

simple:
	$(CXX) $(CXXFLAGS) -o simple simple.cpp $(LIBS) 

matrices:
	$(CXX) $(CXXFLAGS) -o matrices matrices.cpp $(LIBS) 

sos:
	$(CXX) $(CXXFLAGS) -o sos sos.cpp $(LIBS) 

ising:
	$(CXX) $(CXXFLAGS) -o ising ising.cpp $(LIBS) 


