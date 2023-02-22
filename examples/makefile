CXX=g++
CXXFLAGS=-Wfatal-errors -O3
DEBUGFLAGS=-Wall -g -pg -Wfatal-errors
LIBS=-I/usr/include/eigen3 -fopenmp
LIBSCS=-L${SCSHOME}/lib/ -I${SCSHOME}/include/scs/ -lscsdir
HDIR=${MSKHOME}/h
LIBDIR=${MSKHOME}/bin
LDLIBSMOSEK= -I$(HDIR) -L$(LIBDIR) -Wl,-rpath-link,$(LIBDIR) -Wl,-rpath=$(LIBDIR) -lmosek64 -lfusion64

all: mub simple matrices sos ising

mubdebug:
	$(CXX) $(DEBUGFLAGS) -o mub mub.cpp $(LIBS) $(LDLIBSMOSEK)

mub:
	$(CXX) $(CXXFLAGS) -o mub mub.cpp $(LIBS) $(LDLIBSMOSEK) $(LIBSCS)

mubComplex:
	$(CXX) $(CXXFLAGS) -o mubComplex mubComplex.cpp $(LIBS) $(LDLIBSMOSEK)

testing:
	$(CXX) $(CXXFLAGS) -o testing testing.cpp $(LIBS) $(LDLIBSMOSEK)

testingdebug:
	$(CXX) $(DEBUGFLAGS) -o testing testing.cpp $(LIBS) $(LDLIBSMOSEK)

simple:
	$(CXX) $(CXXFLAGS) -o simple simple.cpp $(LIBS) $(LDLIBSMOSEK)

matrices:
	$(CXX) $(CXXFLAGS) -o matrices matrices.cpp $(LIBS) $(LDLIBSMOSEK)

sos:
	$(CXX) $(CXXFLAGS) -o sos sos.cpp $(LIBS) $(LDLIBSMOSEK)

ising:
	$(CXX) $(CXXFLAGS) -o ising ising.cpp $(LIBS) $(LDLIBSMOSEK)


