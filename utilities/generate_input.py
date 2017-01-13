# Copyright (C) 2016 Leo Fang <leofang@phy.duke.edu>
#
# This program is free software. It comes without any warranty,
# to the extent permitted by applicable law. You can redistribute
# it and/or modify it under the terms of the WTFPL, Version 2, as
# published by Sam Hocevar. See the accompanying LICENSE file or
# http://www.wtfpl.net/ for more details.

import sys
from numpy import pi
from math import ceil


########################################
#      Put input parameters below      #
########################################

######### Simulation Parameters ########
# grid size
Delta = 0.01
# number of grids in-between the qubit and its mirror image (nx = 2a/Delta)
nx = 2100
# half of number of (qubit) wavelength in x direction
Nx = 50
# number of (qubit) wavelength in t direction
Ny = 400
# save the result of psi (solution of delay PDE)?
save_psi = 1
# save the result of chi (two-photon wavefunction)?
save_chi = 0
# save psi for every Tstep+1 temporal steps
Tstep = 49

########## Physics Paramters ###########
# initial condition (1: two-photon plane wave; 2: one-photon exponential wavepacket)
init_cond = 2
# exponential tail (used if init_cond=2; dimensionless)
alpha = 0.1
# qubit frequency (in Gamma)
k0 = 20.0
# incident frequency (in Gamma)
k = 20.0
# k0 a = n pi
n = 10.5


########################################
#     Do Not Touch the Code Below!     #
########################################

if len(sys.argv) != 2:
   if len(sys.argv) == 1:
      sys.exit("No filename is given! Abort!")
   else:
      sys.exit("Too many arguments. Abort!")

# convert to number of grid points
Nx = ceil(Nx * nx/n)
Ny = ceil(Ny * nx/n)

# use Delta^-1 as the unit of frequency
if (2.*n*pi/nx) * (k/k0) < pi:
   k = (2.*n*pi/nx) * (k/k0) / Delta
else:
   sys.exit("k is beyond the Nyquist limit. Choose a larger nx for the given n. Abort!")

if 2.*n*pi/nx < pi:
   w0 = 2.*n*pi/nx / Delta
else:
   sys.exit("w0 is beyond the Nyquist limit. Choose a larger nx for the given n. Abort!")

if (2.*n*pi/nx)/k0 < pi:
   Gamma = (2.*n*pi/nx)/k0 / Delta
else:
   sys.exit("Gamma is beyond the Nyquist limit. Choose a larger nx or k0 for the given n. Abort!")

# create the input file
f = open(sys.argv[1], "w")
f.write("nx=%i\n"%nx)
f.write("Nx=%i\n"%Nx)
f.write("Ny=%i\n"%Ny)
f.write("Delta=%.10f\n"%Delta)
f.write("k=%.10f\n"%k)
f.write("w0=%.10f\n"%w0)
f.write("gamma=%.10f\n"%Gamma)
f.write("save_chi=%i\n"%save_chi)
f.write("save_psi=%i\n"%save_psi)
f.write("init_cond=%i\n"%init_cond)
f.write("alpha=%f\n"%alpha)
f.write("Tstep=%i\n"%Tstep)
f.close()
