# GCR CIR Modulation

This is a specialization of the SPECTRUM software applied to modeling galactic cosmic rays (GCR) in the heliosphere. Specifically, the codes in this repository simulate GCR modulation in corotating interaction regions (CIR), i.e. persistent inner heliosphere plasma structures that form when fast solar wind stream overtakes a slow stream, analyze the results, and generate useful figures for a future scientific publication.

For this application, the magnetohydrodynamic (MHD) model of the inner heliosphere is static in corotating coordinates and computed with the BATS-R-US module of the Space Weather Modeling Framework (https://github.com/SWMFsoftware/SWMF). Instructions on how to generate the MHD results are provided in the folder `generate_CIR`, along with the necessary auxiliary files to do so.

## Configuring the code

Before first use, the code must be configured with autotools. After cloning this repository, execute the configure script in the working directory
```
git clone https://github.com/jgaloguz/GCR-CIR-Modulation
cd GCR-CIR-Modulation
./configure.sh <mpi-option> <run-type>
```
where `<mpi-option>` is either `openmpi` or `mpich`, whichever is installed in your system, and `<run-type>` is either `sim` for running modulation simulations or `plot` for running visualization scripts on the MHD results. You may have to change the permissions of `configure.sh` before you can execute it. You will know the configuration stage ran successfully if a `config.h` file was generated in the working directory.

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

## Compiling and running the code

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
mpirun -np <NP> <name-of-code> <inputs>
```
where `<name-of-code>` is the name of the C++ file containing the program you wish to compile and execute (without the `.cc` extension), `<NP>` is the number of processors to use when running the code in parallel (with MPI), and `<inputs>` are any arguments fed to the programs from the terminal (separated by a space).

**Plotting Programs**

There are 3 C++ programs meant to be run with the `plot` configuration as explained in the previous section.
They can be compiled and executed as follows
```
make main_plot_background_1D_batl
mpirun -np 4 main_plot_background_1D_batl <YYYYMMDD> <NUM>
```
plots the main MHD quantities, like flow speed, magnetic field, proton density, etc, as well as derived particle quantities like drifts and diffusion coefficients, along radial cuts and the Earth's trajectory for CIR centered around date `<YYYYMMDD>` and a resolution of `<NUM>` points.
Data files are saved in `output_YYYYMMDD/CIR/` and can be visualized with
```
python plot_background_1D_batl.py <YYYYMMDD> --zero_epoch <zero_epoch>
```
where `<zero_epoch>` is a shift related to the difference between the real and simulated CIR arrival time. These can be found in `data/zero_epochs_SWMF.txt` for the MHD simulation that use the files in the `generate_CIR` folder.
Similarly,
```
make main_plot_background_2D_batl
mpirun -np 4 main_plot_background_2D_batl <YYYYMMDD> <NUM>
```
plots the main MHD quantities, like flow speed, magnetic field, proton density, etc, as well as derived particle quantities like drifts and diffusion coefficients, on 2D meridional and equatorial slices for CIR centered around date `<YYYYMMDD>` and a resolution of `<NUM>` points per dimension.
Data files are saved in `output_YYYYMMDD/CIR/` and can be visualized with
```
python plot_background_2D_batl.py <YYYYMMDD>
```
Finally,
```
make main_plot_fieldline_batl
mpirun -np 4 main_plot_fieldline_batl <YYYYMMDD> <NUM>
```
plots up to `<NUM>` fieldlines, traced outward from Earth's position, for CIR centered around date `<YYYYMMDD>`. Fieldlines are saved in `output_YYYYMMDD/fieldlines/`. The program will look for a list of Earth positions in `data/earth_position_YYYYMMDD.dat`. Two are provided as examples for CIRs with arrival dates 20080228 and 20210615.
Earth positions for other events can be generated using
```
python compute_earth_position.py <YYYYMMDD> [--num <NUM>] [--time <TIME>]
```
where `<YYYYMMDD>` has its usual meaning, and the optional parameters `<NUM>` and `<TIME>` indicate either how many points to evaluate along a full rotation or what particular time for which to calculate Earth's location, respectively. In either case, the Earth position(s) in simulation (HGC) coordinates is saved in `output_YYYYMMDD/earth_position_YYYYMMDD.dat`.

Note that all plotting C++ programs executed in parallel required 4 processes.
This is because one is the Master process, which does not actually participate here, two are Server processes, with only one of them being active, and the remaining process is the Worker that computes the plotted quantities with help from the active Server process.

The plotting script
```
python plot_CIR_ACE.py <YYYYMMDD> <zero_epoch>
```
where both inputs mean the same as before, will produce plots comparing simulated and observed solar wind quantities for a particular event.
Also, the plotting script
```
python plot_CIR_diffusion.py <YYYYMMDD> <zero_epoch>
```
will produce plots of the parallel and perpendicular (physics-based) diffusion coefficients along an 8-day interval centered around `<zero_epoch>`.
In both cases, plots are saved in the `output_YYYYMMDD/` folder.

**Simulation Programs**

There is one main simulation program, compiled and executed with
```
make main_run_modulation_simulation
mpirun -np <NP> main_run_modulation_simulation <inputs>
```
where up to 6 `<inputs>` are given by the following list (in order):

- (1) total number of pseudo-trajectories in simulation
- (2) number of pseudo-trajectories per batch
- (3) CIR date in YYYYMMDD format
- (4) days from observed CIR arrival at which to initialize pseudo-trajectories (zero epoch)
- (5) average perpendicular diffusion coefficient
- (6) average parallel diffusion coefficient

A variety of simulation types can be compiled using certain macros throughout the code.
Any simulation will require inputs (1)-(4), and by default the code is prepared to compile a "full" simulation type, which includes all transport terms in the Parker Transport Equation.
Advection can be omitted by commenting out `TRAJ_PARKER_USE_U_ADVECT` in `../src/trajectory_parker.hh` ("nadv").
Drifts can be omitted by commenting out `TRAJ_PARKER_USE_B_DRIFTS` in `../src/trajectory_parker.hh` ("ndft").
The diffusion model can be set through the `CONST_DIFF` macro in `main_run_modulation_simulation.cc`
A value of `0` uses physics-based expressions.
A value of `1` or `2` substitute perpendicular ("ndif0") or parallel ("ndif1") diffusion, respectively, with empirical formulas that do not depend on the CIR structure, while a value of `3` substitutes both ("ndif").
In such cases, all inputs (1)-(6) are required.

Inputs (5) and (6) can be calculated using
```
python compute_avg_diff.py <YYYYMMDD>
```
where `<YYYYMMDD>` has its usual meaning. This command can only be run *after* the 1D plotting C++ program, since it directly averages the values along Earth's trajectory.
Additionally, note that all simulations a bent power-law fit of to the force-field approximation proton spectrum at 5 au hard-coded.
To see a plot of the spectrum fit for year `<YYYY>` in range [2007,2010] or [2018,2021] run
```
python plot_fit_spectrum_5au.py <YYYY>
```

To post-process results for a single simulation, compile and run the following programs
```
make main_post_process_spectrum
./main_post_process_spectrum <YYYYMMDD> <t0>
make main_post_process_time
./main_post_process_time <YYYYMMDD> <t0>
make main_post_process_position
./main_post_process_position <YYYYMMDD> <t0>
```
where `<YYYYMMDD>` and `<t0>` are inputs (3) and (4) for the `main_run_modulation_simulation.cc` code.
After a series of `<NS>` many simulations have been run for equally spaced start times in a range spanning 8 days, these can all be integrated to produce a GCR response over time with
```
make main_post_integrate_spectra_over_time
./main_post_integrate_spectra_over_time <YYYYMMDD> <NS> <low_epoch>
```
where `<low_epoch>` is the lower bound of said 8-day range.
The output files of the particle simulation code are stored in `output_YYYYMMDD/GCR/`, where the post-processing files will look for them.
Prior to plotting them, this folder must be relabeled `output_YYYYMMDD/GCR_<type>`, where `<type>` can be `full`, `nadv`, `ndft`, `ndif`, `ndif0`, or `ndif1`, as described previously.

```
python plot_CIR_modulation.py <YYYYMMDD> <zero_epoch>
```
where the inputs have their usual meaning, will compare the results from the "full", "nadv", "ndft", and "ndif" runs.
Similarly,
```
python plot_CIR_modulation_ndifs.py <YYYYMMDD> <zero_epoch>
```
where the inputs have their usual meaning, will compare the results from the "full", "ndif0", and "ndif1" runs.
For completion,
```
python compare_CIR_residence_times.py data/zero_epochs_SWMF.txt
```
will plot the change in residence times between ndif and full simulations vs the radial diffusion coefficient along Earth's trajectory for all events in `data/zero_epochs_SWMF.txt`, and
```
python compare_CIR_modulation.py data/zero_epochs_SWMF.txt [--phase <PHASE>]
```
will fit the slopes of the observed and simulated GCR percentage change for all events listed in `data/zero_epochs_SWMF.txt`, producing plots of such fits for each event in its respective output folder as well as summary plot of all fits for the chosen phases.
The phases, indexed as 0, 1, or 2 with the optional parameter `<PHASE>` represent different epoch ranges, [-4,-2] (onset), [-1,1] (drop), and [2,4] (recovery). 

## Important note

**This is NOT the official SPECTRUM repository.** For information about SPECTRUM, go to https://github.com/vflorins/SPECTRUM.
