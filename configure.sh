#!/bin/bash

# Output status if successful
output_status() {
   if [ $? == 0 ]; then
      echo ""
      echo "****************************************"
      echo "Code configured for:"
      echo "$1"
      echo "****************************************"
      echo ""
   fi 
}

autoreconf
automake --add-missing

if [ "$#" -ne 2 ]
then
   echo "ERROR: Two input arguments are required."
   echo "   Arg 1 : MPI library to use during compilation (openmpi or mpich)"
   echo "   Arg 2 : Mode of execution (sim or plot)"
   exit 1
fi

if [ ${2} = "sim" ]
then
   ./configure CXXFLAGS="-Ofast" --with-mpi=${1} --with-execution=PARALLEL --with-trajectory=PARKER --with-time_flow=BACKWARD --with-rkmethod=0 --with-server=BATL --with-server_interp_order=1 --with-server_num_gcs=1
   output_status "    GCR modulation simulations."
elif [ ${2} = "plot" ]
then
   ./configure CXXFLAGS="-Ofast" --with-mpi=${1} --with-execution=PARALLEL --with-trajectory=FIELDLINE --with-time_flow=FORWARD --with-rkmethod=25 --with-server=BATL --with-server_interp_order=1 --with-server_num_gcs=1
   output_status "    MHD results visualization."
else
   echo "ERROR: Second input argument not recognized."
   echo "Valid inputs are:"
   echo "   sim : use this input when configuring the code for test-particle simulations"
   echo "   plot : use this input when configuring the code for plotting routines"
fi