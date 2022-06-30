# poly.h

This is a small library for dealing with polynomials in c++.

### Features

- create polynomials or systems of polynomials of abritrary degree
- real and complex polynomials supported
- easy stream output
- perform addition, subtraction and multiplication
- root finding algorithm
- symbolic integration and differentiation
- Hilbert's Nullstellensatz for systems of polynomials

### Compiling

This depends on Eigen (for matrix manipulations) and openmp (for parallel root finding).
```bash
sudo apt install libeigen3-dev libomp-dev
```

Then you need to simply mention "poly.h" in your cpp includes. So you can put "poly.h" in a location with other c++ headers or simply refer directly to the header file. 

### Examples

A simple example showing basic polynomial construction and output is given in the examples folder as "simple.cpp". It can be compiled using make:
```bash
make simple
./simple
```

A more useful example shows how this library can be used to find mutually unbiased bases using the polynomial version of this problem:
```bash
make mub
./mub
```


