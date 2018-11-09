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
# grid step size (just give an arbitrary positive number here; see the doc)
Delta = pi/12000.
# desired accuracy (=Delta/lambda, the ratio of step size to wavelength)
dx = 1.0/240.0
# half of number of (qubit) wavelength in x direction
Nx = 270
# number of (qubit) wavelength in t direction
Ny = 270
# save the result of psi (solution of delay PDE)?
save_psi = 1
# save the result of \int dx |psi(x,t)|^2 ?
save_psi_square_integral = 0
# save the result of chi (two-photon wavefunction)?
save_chi = 0
# save the result of chi as a 2D map?
save_chi_map = 0
# whether or not save e0(t) and e1(t)
measure_NM = 0
# save psi for every Tstep+1 temporal steps
Tstep = 49
# number of solvers
Nth = 1

########## Physics Paramters ###########
# initial condition 
# 1: two-photon plane wave;
# 2: one-photon exponential wavepacket;
# 3. two-photon exponential wavepacket
init_cond = 1
# qubit frequency (in Gamma)
k0 = 20.0
# k0 a = n pi
n = 10.0
# exponential tail (used if init_cond=2; dimensionless)
alpha = 0.1
# incident frequency (in Gamma)
k = 20.0
# two identical photons? Used when init_cond=3
identical_photons = 1
# exponential tail for photon #1
alpha1 = 0.1
# exponential tail for photon #2
alpha2 = 0.1
# incident frequency for photon #1
k1 = 20.0
# incident frequency for photon #2
k2 = 20.0

########################################
#     Do Not Touch the Code Below!     #
########################################

if len(sys.argv) != 2:
   if len(sys.argv) == 1:
      sys.exit("Usage: python generate_input.py input_filename")
   else:
      sys.exit("Too many arguments. Abort!")

# compute nx based on accuracy goal
nx = int(round(n / dx))

# convert to number of grid points
Nx = ceil(Nx * nx/n)
Ny = ceil(Ny * nx/n)

# use Delta^-1 as the unit of frequency
if (2.*n*pi/nx) * (k/k0) < pi:
   k_in = (2.*n*pi/nx) * (k/k0) / Delta
   k1_in = (2.*n*pi/nx) * (k1/k0) / Delta
   k2_in = (2.*n*pi/nx) * (k2/k0) / Delta
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
f.write("Delta=%.15E\n"%Delta)
f.write("w0=%.15E\n"%w0)
f.write("gamma=%.15E\n"%Gamma)
f.write("init_cond=%i\n"%init_cond)
f.write("identical_photons=%i\n"%identical_photons)
if init_cond != 3 or (init_cond==3 and identical_photons==1):
   f.write("k=%.15E\n"%k_in)
else: 
   f.write("k1=%.15E\n"%k1_in)
   f.write("k2=%.15E\n"%k2_in)
if init_cond == 2 or (init_cond==3 and identical_photons==1):
   f.write("alpha=%.15E\n"%alpha)
if (init_cond==3 and identical_photons==0):
   f.write("alpha1=%.15E\n"%alpha1)
   f.write("alpha2=%.15E\n"%alpha2)
f.write("save_chi=%i\n"%save_chi)
f.write("save_chi_map=%i\n"%save_chi_map)
f.write("save_psi=%i\n"%save_psi)
f.write("save_psi_square_integral=%i\n"%save_psi_square_integral)
f.write("measure_NM=%i\n"%measure_NM)
f.write("Tstep=%i\n"%Tstep)
f.write("Nth=%i\n"%Nth)
f.close()
