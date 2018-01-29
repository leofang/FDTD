# Copyright (C) 2016 Leo Fang <leofang@phy.duke.edu>
#
# This program is free software. It comes without any warranty,
# to the extent permitted by applicable law. You can redistribute
# it and/or modify it under the terms of the WTFPL, Version 2, as
# published by Sam Hocevar. See the accompanying LICENSE file or
# http://www.wtfpl.net/ for more details.

import sys, os
import numpy
import datetime
from numpy import pi
from math import ceil


########################################
#      Put input parameters below      #
########################################

######### Simulation Parameters ########
# grid size
Delta = pi/12000.
#Delta=1E-6
# the ratio of grid size to wavelength (Delta/lambda)
dx = 1.0/240.0
# half of number of (qubit) wavelength in x direction
Nx = 260
# number of (qubit) wavelength in t direction
Ny = 260
# save the result of psi (solution of delay PDE)?
save_psi = 1
# save the result of chi (two-photon wavefunction)?
save_chi = 0
# save psi for every Tstep+1 temporal steps
Tstep = 29
# whether or not save e0(t) and e1(t)
measure_NM = 1

########## Physics Paramters ###########
# initial condition (1: two-photon plane wave; 2: one-photon exponential wavepacket)
init_cond = 2
# qubit frequency (in Gamma)
k0 = 20.0
# k0 a = n pi
n = 1.0
# two identical wavapackets?
identical_photons = 1
# incident frequency for identical photons (in Gamma)
k = 20.0
# incident frequency for photon #1 (in Gamma)
k1 = 16.5
# incident frequency for photon #2 (in Gamma)
k2 = 17.5

########## condor configuration ###########
OnlyUseNanoMachines = True
RunOnOSG = False
NumCPU = 32

########################################
#     Do Not Touch the Code Below!     #
########################################
   
########## generate condor submit file ###########
f=open("condor_submit_file","w")

#f.write("Universe = parallel\n") # run MPI jobs ---- seems not working!
#f.write("machine_count = 1\n") # number of cores for a MPI job
f.write("Universe = vanilla\n")
f.write("Executable = /var/phy/project/baranger/yf30/fdtd_test/FDTD/FDTD \n")
f.write("notification=Error\n")
f.write("notify_user = leofang@phy.duke.edu\n")
if (RunOnOSG if "RunOnOSG" in locals() else False): # run on Duke Ci-Connect
   f.write("+ProjectName=\"duke-CMT\"\n")
   f.write("should_transfer_files=YES\n") # output files are transferred
   f.write("when_to_transfer_output = ON_EXIT\n")
   # Periodically retry the jobs every 60 seconds, up to a maximum of 5 retries.
   f.write("periodic_release =  (NumJobStarts < 5) && ((CurrentTime - EnteredCurrentStatus) > 60)\n")
   # Send the job to Held state on failure. 
   f.write("on_exit_hold = (ExitBySignal == True) || (ExitCode != 0)\n")
   # Reconnect jobs faster if disconnected (default is 2400 sec)
   f.write("JobLeaseDuration = 90\n")
else:  # run on physics condor
   f.write("+Department = \"Physics\"\n")
   f.write("should_transfer_files=NO\n") # use Physics shared files system
   requirement = "requirements = ( OpSys == \"LINUX\" && Arch ==\"X86_64\" && FileSystemDomain != \"\" )"
#   # do not use atl or phy-compute machines because they do not have OpenMPI installed
#   for i in [11, 12, 13, 14, 15, 16]:
#      requirement += "(TARGET.Machine != \"atl0%i.phy.duke.edu\")" %i
#      requirement = (requirement+")" if i==16 else requirement+" && ")
#   requirement += " && ("
#   for i in [1, 4]:
#      requirement += "(TARGET.Machine != \"phy-compute-0%i.phy.duke.edu\")" %i
#      requirement = (requirement+")" if i==4 else requirement+" && ")
   # only use nano machines 
   if OnlyUseNanoMachines is True:
      requirement += " && ("
      for i in [6, 7, 8, 9]:
          requirement += "(TARGET.Machine == \"nano0%i.internal.phy.duke.edu\")" %i
          requirement = (requirement+")" if i==9 else requirement+" || ")
   # my code can only run on SL6 machines since the condor was back online on Sep 21, 2016 
#   requirement += "&& OPSYSMAJORVER == 6" 
   f.write(requirement+"\n")
   #f.write("request_memory = 50\n")
#   f.write("getenv=True\n") # important for finding shared libraries!
#f.write("periodic_release = (NumGlobusSubmits < 5) && ((CurrentTime - EnteredCurrentStatus) > (60*60))\n")
#f.write("periodic_hold =  (JobStatus==2) && ((CurrentTime - EnteredCurrentStatus) > (%i*60*60))" % MAX_WALL_TIME)
f.write("\n")

# create a folder for condor logs
if not os.path.exists("condor_log"): # create one if the directory does not exist yet
   os.makedirs("condor_log")

values = range(-11, 12)
#values.remove(-10); values.remove(0); values.remove(10);

# exponential tail (used if init_cond=2; dimensionless)
for alpha in 10.**(numpy.array(values)/10.): # 10.**numpy.array([-1, 0, 1]): # :

   #if len(sys.argv) != 2:
   #   if len(sys.argv) == 1:
   #      sys.exit("No filename is given! Abort!")
   #   else:
   #      sys.exit("Too many arguments. Abort!")
   if abs(alpha-1.0)<1E-10:
      #continue
      alpha = 0.99999 # no effect after continue
   
   #nx = int(n / dx)
   nx = int(round(n / dx))

   time_object = datetime.date.today()
   if identical_photons == 1:
      filename = time_object.strftime("%Y%b%d") + "_n_%.2f_alpha_%.4f_k_%.2f" %(n, alpha, k)
   else:
      #same width: alpha1=alpha2
      filename = time_object.strftime("%Y%b%d") + "_n_%.2f_alpha1_%.4f_alpha2_%.4f_k1_%.2f_k2_%.2f" %(n, alpha, alpha, k1, k2) 
   
   # convert to number of grid points
   Nx_in = ceil(Nx * nx/n)
   Ny_in = ceil(Ny * nx/n)
   
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
   f2 = open(filename, "w")
   f2.write("nx=%i\n"%nx)
   f2.write("Nx=%i\n"%Nx_in)
   f2.write("Ny=%i\n"%Ny_in)
   f2.write("Delta=%.15E\n"%Delta)
   if identical_photons == 1:
      f2.write("k=%.15E\n"%k_in)
   f2.write("w0=%.15E\n"%w0)
   f2.write("gamma=%.15E\n"%Gamma)
   f2.write("save_chi=%i\n"%save_chi)
   f2.write("save_psi=%i\n"%save_psi)
   f2.write("measure_NM=%i\n"%measure_NM)
   f2.write("init_cond=%i\n"%init_cond)
   if identical_photons == 1:
      f2.write("alpha=%f\n"%alpha)
   f2.write("Tstep=%i\n"%Tstep)
   f2.write("identical_photons=%i\n"%identical_photons)
   if identical_photons == 0:
      f2.write("k1=%.15E\n"%k1_in)
      f2.write("k2=%.15E\n"%k2_in)
      f2.write("alpha1=%f\n"%alpha) # same width: alpha1=alpha2
      f2.write("alpha2=%f\n"%alpha) # same width: alpha1=alpha2
   f2.close()

   # condor submit file
   f.write("Initialdir = ./\n")
   f.write("output = condor_log/" + "$(Cluster).$(Process).out\n")
   f.write("error  = condor_log/" + "$(Cluster).$(Process).err\n")
   f.write("Log    = condor_log/" + "$(Cluster).$(Process).log\n")
   f.write("request_cpus = %i\n" % NumCPU)
   #f.write("request_memory = %.fGB\n" %(2*Nx_in*Ny_in*16.*1.1/1E+9) )
   #if (RunOnOSG if "RunOnOSG" in locals() else False): # run on OSG nodes
   #    f.write("transfer_input_files = " + parms['BASENAME'] + ".in.h5")
   #    for i in range(N_ENV if "N_ENV" in locals() else 1):
   #        f.write(', DELTA' + str(i) + ".h5")
   #    f.write('\n')
   #f.write("transfer_input_files = " + parms['BASENAME'] + ".in.h5," + parms['BASENAME'] + ".out.h5\n")
   f.write("Arguments = " + filename + "\n")
   f.write("Queue\n\n")

f.close()
