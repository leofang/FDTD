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
nx = 20
# half of number of (qubit) wavelength in x direction
Nx = 10
# number of (qubit) wavelength in t direction
Ny = 200

########## Physics Paramters ###########
# qubit frequency (in Gamma)
k0 = 20.0
# incident frequency (in Gamma)
k = 20.0
# k0 a = n pi
n = 0.5



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
f.write("Gamma=%.10f\n"%Gamma)
f.close()
