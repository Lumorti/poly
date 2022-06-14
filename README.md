
# LMD - Low Memory Nullstellensatz

This is a small executable for finding an effective Hilbert's Nullstellensatz for a system of polynomial equations. Rather than creating an exponentially growing linear system (of which only a few rows might be needed), it iteratively performs a low-memory Guassian elimation on the system of equations.

### Compiling

You'll need basic compiling stuff, although you probably already have it:
```bash
sudo apt-get install make g++
```

Then simply clone the repo, enter the repo, then run:
```bash
make
```

### Usage

To run the example:
```bash
./lmd example.eqn
```

To view the help:
```bash
./lmd -h
```
