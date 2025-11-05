# Generate CIRs with SWMF

To generate the magnetohydrodynamic (MHD) solar wind background used in these project, you must install, configure, and run the Space Weather Modeling Framework (SWMF). To do this, execute the following

```
git clone https://github.com/SWMFsoftware/SWMF.git
cd SWMF
Config.pl -install -compiler=gfortran
```
Prior to the installation command, you might have to edit the SWMFsoftware url in `share/Scripts/gitclone` from each instance of `git\@` to `https://`. Several sub-dirs/repos will ask for username/password and fail but this is ok.

You may also need/want the additional data repositories.
```
git clone https://github.com/SWMFsoftware/SWMFSOLAR.git
git clone https://github.com/SWMFsoftware/SWMF_data.git
```

Once SWMF is installed, you can generate the MHD results through the following procedure. In the top-level SWMF folder, run
```
./Config.pl -v=Empty,SC/BATSRUS,IH/BATSRUS
./Config.pl -o=SC:u=Awsom,e=AwsomAnisoPiSA,g=6,8,8,ng=2
./Config.pl -o=IH:u=Awsom,e=AwsomAnisoPiSA,g=8,8,8,ng=2
make SWMF
make rundir
```
This will compile the code for execution and generate a run directory. Next, choose one of the dates for which to run the solar wind MHD model within the `run_params`. The subfolders are named in a `<YYMMDD>` format representing the starting date of the 27-day period that is being simulated up to a steady-state solution in corotating coordinates. Each of these subfolders contains a PFSS harmonics coefficients file `map_<YYMMDD>.dat`, computed from GONG magnetograms (https://gong.nso.edu/data/magmap/archive.html) as well as a parameters file with specific instructions for the simulation `PARAM.in`.

To proceed, rename the `rundir` folder `run_dir_<YYMMDD>` where `<YYMMDD>` is your choice of date from the available options.
```
mv rundir run_cir_<YYMMDD>
```
This step is not critical for SWMF to run, but it is important so that the results are correctly read by the routines set up in the particle tracing components of this repository. Next, place the harmonic coefficients and parameter files in the run directory. Assuming `SWMF` and `GCR-CIR-Modulation` are sibling directories, this can be done with
```
cp ../GCR-CIR-Modulation/generate_CIR/run_params/<YYMMDD>/* run_cir_<YYMMDD>/.
```

The MHD simulation is now ready to run and can be executed in parallel with MPI like so
```
cd run_cir_<YYMMDD>
mpirun –np <N> ./SWMF.exe > run.log
```
which generates the log file `run.log` containing runtime information of the program. We recommend using a large number of cores (hundreds), which will likely require running the simulation in a managed cluster of some sort. Sample scripts to submit jobs on clusters with a variety of work managers can be found in `SWMF/share/JobScripts`.

Once the simulation finishes, the results must be post-processed. This is easily done by executing the following commands from the top-level SWMF directory
```
make PIDL
cd run_cir_<YYMMDD>
./PostProc.pl
```
After this step, the simulation results, found in `SWMF/run_cir_<YYMMDD>/IH/IO2/`, have been converted to binary format, which the particle tracking software can read.