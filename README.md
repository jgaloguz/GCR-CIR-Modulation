# GCR CIR Modulation

This is a specialization of the SPECTRUM software applied to modeling galactic cosmic rays (GCR) in the heliosphere. Specifically, the codes in this repository simulate GCR modulation in corotating interaction regions (CIR), i.e. persistent inner heliosphere plasma structures that form when fast solar wind stream overtakes a slow stream, analyze the results, and generate useful figures for a future scientific publication.

For this application, the magnetohydrodynamic (MHD) model of the inner heliosphere is static in corotating coordinates and computed with the BATS-R-US module of the Space Weather Modeling Framework (https://github.com/SWMFsoftware/SWMF). Instructions on how to generate the MHD results are provided in the folder `generate_CIR`, along with the necessary auxiliary files to do so. 

## Using the code

Before first use, the code must be configured with autotools. After cloning this repository, execute the configure script in the working directory
```
git clone https://github.com/jgaloguz/GCR-CIR-Modulation
cd GCR-CIR-Modulation
./configure.sh <mpi-option>
```
where `<mpi-option>` is either `openmpi` or `mpich`, whichever is installed in your system. You may have to change the permissions of `configure.sh` before you can execute it. You will know the configuration stage ran successfully if a `config.h` file was generated in the working directory.

If the configuration stage runs successfully, you can compile and run codes in the `runs` folder. To compile the code, simply use
```
cd runs
make <name-of-code>
```
Once the code is compiled, you can run it by executing
```
./<name-of-code> <inputs>
```
or
```
mpirun -np <N> <name-of-code> <inputs>
```
where `<name-of-code>` is the name of the C++ file containing the program you wish to compile and execute (without the `.cc` extension), `<N>` is the number of processors to use when running the code in parallel (with MPI), and `<inputs>` are any arguments fed to the programs from the terminal (separated by a space).

The SWMF Block Adaptive Tree Library (BATL) must be cloned, configured, and installed on the same level as this repository. To do so, execute the following commands in the working directory (outside of `GCR-CIR-Modulation`)
```
git clone https://github.com/SWMFsoftware/BATL.git
cd BATL
git clone https://github.com/SWMFsoftware/share.git
git clone https://github.com/SWMFsoftware/util.git
git clone https://github.com/SWMFsoftware/srcBATL.git src
./Config.pl -install -compiler=gfortran
./Config.pl -g=8,8,8 -r=2,2,2 -ng=1 -single
make READAMRLIB
```
Ensure that mpi is enabled when running the last line.

## Important note

**This is NOT the official SPECTRUM repository.** For information about SPECTRUM, go to https://github.com/vflorins/SPECTRUM.
