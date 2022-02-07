CXX=g++
CXXFLAGS=-Wfatal-errors -O3 -Wall

all:
	$(CXX) $(CXXFLAGS) -o lmd main.cpp
