# FDTD: solving 1+1D delay PDE
The code is intended to solve 1+1D delay PDE which emerges in waveguide QED problems.

## Requirements
None. Any modern OS has provided standard C libraries, and the code just depends on it, as it is written in pure C and conforming the gnu99 standard.

## Installation
A makefile is provided. Type `make` in the same folder to compile, and an executable named `FDTD` will be generated.

## Usage
`./FDTD input_filename`

## Mandatory input parameters
nx
Nx
Ny
Delta
k
w0
Gamma

See the comments in `grid.h` for the layout of the grid and detail of various parameters. A Python script is provided in the `utilities` folder for users' convenience.

## Author
Leo Fang

## License
Not decided yet.
