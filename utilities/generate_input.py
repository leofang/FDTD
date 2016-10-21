import sys
from numpy import pi
from math import ceil


########################################
#      Put input parameters below      #
########################################

######### Simulation Parameters ########
# nx = 2a/Delta
nx = 20
# horizontal (half) box size (in qubit wavelength)
Nx = 10
# vertical box size (in qubit wavelength)
Ny = 200
# grid size
Delta = 0.01

########## Physics Paramters ###########
# qubit frequency (in Gamma)
k0_n = 20
# incident frequency (in Gamma)
kn = 20
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

# in units of Delta
Nx = ceil(Nx * nx/n)
Ny = ceil(Ny * nx/n)

# in units of Delta^-1
k = (2.*n*pi/nx) * (kn/k0_n)
w0 = 2.*n*pi/nx
Gamma = (2.*n*pi/nx)/k0_n

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
