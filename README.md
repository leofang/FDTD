# FDTD: solving 1+1D delay PDE
This code is intended to solve 1+1D delay PDE which emerges in waveguide QED problems.

## Requirements
None. This code is written in pure C and conforming the gnu99 standard, so it can be compiled on any modern OS with standard C libraries.

## Features
* Fast and efficient
* Valgrind-clean (no memory leak)
* Proof of concept for numerically solving a PDE which has delay in both dimensions

## Installation
A makefile is provided. After cloning the git repo or downloading the source code, simply type `make` in the same folder to compile, and an executable named `FDTD` will be generated.

## Usage
`./FDTD input_filename`, where `input_filename` is the name of the input file.

## Mandatory input parameters
Only 7 are required: `nx`, `Nx`, `Ny`, `Delta`, `k`, `w0`, and `Gamma`. The first four are simulation-related, and the rest are physics-related. The format of the input file should be one key-value pair per line, with a equal sign separating the key and the value (see the sample file `k0a_0.5pi_on`):
```bash
nx=200
Nx=4000
Ny=40000
Delta=0.0100000000
k=1.5707963268
w0=1.5707963268
Gamma=0.0785398163
```
For futher details (e.g., the layout of the grid, decriptions for various parameters, etc) see the comments in `grid.h`. A Python script is provided in the `utilities` folder for the ease of preparing input.

Currently only one kind of initial condition (two-photon plane wave) is built in. Other forms are planned to be included in the future.

## Output
Currently two files will be generated: `input_filename.re.out` and `input_filename.im.out` (real and imaginary parts of the wavefunction, respectively). A Mathematica notebook is provided in the `utilities` folder for simple plotting purposes.

## License
This code is released under the [WTFPL](http://www.wtfpl.net) license. For the academic uses, citation to *(ref)* is strongly encouraged but not required. Copyright (c) 2016 Leo Fang.
