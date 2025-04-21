CXX=g++
CXXFLAGS=-Wfatal-errors -O3 -march=native
#CXXFLAGS=-Wfatal-errors -O0 -g 
DEBUGFLAGS=-Wall -g -pg -Wfatal-errors
LIBSEIGEN= -I/usr/include/eigen3
LIBSOMP= -fopenmp
LIBSSCS= -L${SCSHOME}/lib/ -I${SCSHOME}/include/scs/ -lscsdir
LIBSMOSEK= -I${MSKHOME}/h -L${MSKHOME}/bin -Wl,-rpath-link,${MSKHOME}/bin -Wl,-rpath=${MSKHOME}/bin -lmosek64 -lfusion64
LIBSOPTIM= -I${OPTIMHOME}/include/ -L${OPTIMHOME} -loptim
LIBS=$(LIBSEIGEN) $(LIBSOPTIM) $(LIBSMOSEK) $(LIBSSCS) $(LIBSOMP)

all: mub simple sos ising causal

mubdebug: mub.cpp ../poly.h
	$(CXX) $(DEBUGFLAGS) -o mub mub.cpp $(LIBS) 

mub: mub.cpp ../poly.h
	$(CXX) $(CXXFLAGS) -o mub mub.cpp $(LIBS) 

simple: simple.cpp ../poly.h
	$(CXX) $(CXXFLAGS) -o simple simple.cpp $(LIBS) 

matrices: matrices.cpp ../poly.h
	$(CXX) $(CXXFLAGS) -o matrices matrices.cpp $(LIBS) 

sos: sos.cpp ../poly.h
	$(CXX) $(CXXFLAGS) -o sos sos.cpp $(LIBS) 

ising: ising.cpp ../poly.h
	$(CXX) $(CXXFLAGS) -o ising ising.cpp $(LIBS) 

causal: causal.cpp ../poly.h
	$(CXX) $(CXXFLAGS) -o causal causal.cpp $(LIBS)
