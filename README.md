# FDTD: solving 1+1D delayed PDE
This code is intended to solve 1+1D complex-valued, delayed PDE which emerges in waveguide-QED problems: ![](https://latex.codecogs.com/gif.latex?%5Cinline%20%5Cbegin%7Balign%7D%20%5Cfrac%7Bd%7D%7Bdt%7D%5Cpsi%28x%2Ct%29%26%3D-%5Cfrac%7Bd%7D%7Bdx%7D%5Cpsi%28x%2Ct%29-%5Cleft%28i%5Comega_0&plus;%5Cfrac%7B%5CGamma%7D%7B2%7D%5Cright%29%5Cpsi%28x%2Ct%29&plus;%5Cfrac%7B%5CGamma%7D%7B2%7D%5Cpsi%28x-2a%2C%20t-2a%29%5Ctheta%28t-2a%29%5Cnonumber%5C%5C%20%26%5Cquad-%5Cfrac%7B%5CGamma%7D%7B2%7D%5Cbiggl%5B%5Cbigl%28%5Cpsi%28-x-2a%2C%20t-x-a%29-%5Cpsi%28-x%2C%20t-x-a%29%5Cbigr%29%5Ctheta%28x&plus;a%29%5Ctheta%28t-x-a%29%5Cnonumber%5C%5C%20%26%5Cquad%5Cquad&plus;%5Cbigl%28%5Cpsi%282a-x%2C%20t-x&plus;a%29-%5Cpsi%28-x%2C%20t-x&plus;a%29%5Cbigr%29%5Ctheta%28x-a%29%5Ctheta%28t-x&plus;a%29%5Cbiggr%5D%5Cnonumber%5C%5C%20%26%5Cquad&plus;%5Csqrt%7B%5Cfrac%7B%5CGamma%7D%7B4%7D%7D%5Cbiggl%5B%5Cchi%28x-t%2C-a-t%2C0%29&plus;%5Cchi%28-a-t%2Cx-t%2C0%29-%5Cchi%28x-t%2Ca-t%2C0%29-%5Cchi%28a-t%2Cx-t%2C0%29%5Cbiggr%5D%5Cnonumber%20%5Cend%7Balign%7D)

## Requirements
None. This code is written in pure C and conforming the gnu99 standard, so it can be compiled on any modern OS with C compilers that conform C99. Tested with gcc on Linux and clang on Mac. Depending on the grid size, however, the code can be highly memory- and storage-demanding; see below.

## Features
* Fast and efficient
* Valgrind-clean (no memory leak)
* Proof of concept for numerically solving a PDE with delay in both dimensions using FDTD

## Installation
A makefile is provided. After cloning the git repo or downloading the source code, simply type `make` in the same folder to compile, and an executable named `FDTD` will be generated.

## Usage
`./FDTD input_filename`, where `input_filename` is the name of the input file that specifies the input parameters, each in one line (see below).

## Input parameters
At least 8 are required: `nx`, `Nx`, `Ny`, `Delta`, `init_cond`, `k`, `w0`, and `gamma`. The first four are simulation-related, and the rest are physics-related. The format of the input file should be one key-value pair per line, with a equal sign separating the key and the value (see the sample file `k0a_0.5pi_on_new`):
```bash
nx=200
Nx=10000
Ny=20000
Delta=0.0100000000
init_cond=1
k=1.5707963268
w0=1.5707963268
gamma=0.0785398163
save_psi=0
save_chi=1
```
For futher details (e.g., the layout of the grid, decriptions for various parameters, etc) see the comments in `grid.h` as well as the documentation. A Python script is provided in the `utilities` folder for the ease of preparing input.

Currently two kinds of initial conditions are built in: two-photon plane wave (set `init_cond=1`) and single-photon exponential wavepacket (`init_cond=2`). For the latter, the (dimensionless) wavepacket width `alpha` needs to be specified. Other kinds of initial conditions can be incorperated into the code easily.

For the ease of post-processing data, several functions for constructing various non-Markovian measures can be calculated on the fly if `measure_NM=1` is set. See the documentation for detail. Note that currently in this situation only the single-photon exponential wavepacket is supported (so remember to set `init_cond=2` and `alpha`).

Other options controlling the behavior of the program can also be given; if not given, the program assumes a default value. Currently all available options are `save_psi` (default=0), `save_psi_binary` (default=0), `save_chi` (default=0), `init_cond` (default=0: invalid), `Tstep` (default=0), and `measure_NM` (default=0).

**WARNING**: depending on the grid size, the memory usage and the output files can be excessively huge. For the former, a quick estimation is 2\*16\*Nx\*Ny/1024^3 (in GB); for the latter, setting `Tstep=30` or larger (write the wavefunction for every Tstep+1 temporal steps) can help reduce significantly the file size.

## Output
Depending on the options, the following files will be generated: 
* `save_psi`: `input_filename.re.out` and `input_filename.im.out` (real and imaginary parts, respectively, of the wavefunction described by the delay PDE). 
* `save_psi_binary`: `input_filename.bin` (the entire wavefunction, complex numbers, written in a binary file).
* `save_chi`: `input_filename.abs_chi.out` (absolute value of the two-photon wavefunction).
* `measure_NM`: `input_filename.re_e0.out`, `input_filename.re_e1.out`, `input_filename.re_mu.out`, their imaginary counterparts, and `input_filename.lambda.out`; see the documentation for their meanings.

Note that (i) these options cannot be simultaneously turned off, or the program would generate nothing; (ii) for the wavefunctions, each row in the output file gives the wavefunction along the x-direction starting from x=-a, and rows are written in the order t=0, t=Tstep+1, t=2(Tstep+1), ...

A Mathematica notebook is provided in the `utilities` folder for simple plotting purposes.

## License
This code is released under the [WTFPL v2](http://www.wtfpl.net). For the academic uses, citation to **(ref)** is *strongly encouraged and acknowledged.* Copyright (C) 2016 Leo Fang.

This program is free software. It comes without any warranty, to the extent permitted by applicable law. You can redistribute it and/or modify it under the terms of the WTFPL, Version 2, as published by Sam Hocevar. See the accompanying LICENSE file or http://www.wtfpl.net/ for more details.
