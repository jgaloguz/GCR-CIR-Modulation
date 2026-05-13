#include "src/simulation.hh"
#include "src/distribution_other.hh"
#include "src/background_server_batl.hh"
#include "src/diffusion_other.hh"
#include "src/boundary_time.hh"
#include "src/boundary_space.hh"
#include "src/initial_time.hh"
#include "src/initial_space.hh"
#include "src/initial_momentum.hh"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <filesystem>
#include <cstdlib>

using namespace Spectrum;

#define CONST_DIFF 0

int main(int argc, char** argv)
{
   DataContainer container;

//--------------------------------------------------------------------------------
// Create a simulation object
//--------------------------------------------------------------------------------

   std::unique_ptr<SimulationWorker> simulation;
   simulation = CreateSimulation(argc, argv);

//--------------------------------------------------------------------------------
// Particle type
//--------------------------------------------------------------------------------

   int specie = SPECIES_PROTON_BEAM;
   simulation->SetSpecie(specie);

//--------------------------------------------------------------------------------
// Background
//--------------------------------------------------------------------------------

   container.Clear();

// Initial time
   double t0 = 0.0;
   container.Insert(t0);

// Origin
   container.Insert(gv_zeros);

// Velocity
   container.Insert(gv_zeros);

// Magnetic field
   container.Insert(gv_zeros);

// Effective "mesh" resolution
   double one_au =  GSL_CONST_CGSM_ASTRONOMICAL_UNIT / unit_length_fluid;
   double dmax = 0.1 * one_au;
   container.Insert(dmax);

// Name of file with MHD results
   std::string cir_date = "";
   if (argc > 3) {
      cir_date = argv[3];
      if (MPI_Config::is_master) {
         std::cout << "CIR date: " << cir_date << std::endl;
         std::cout << "Number of servers: " << MPI_Config::boss_comm_size-1 << std::endl;
         std::cout << "Number of workers: " << MPI_Config::n_workers << std::endl;
         std::filesystem::create_directory("output_" + cir_date);
         std::filesystem::create_directory("output_" + cir_date + "/GCR");
      };
   } else {
      if (MPI_Config::is_master) std::cout << "ERROR: No CIR date provided." << std::endl;
      exit(1);
   };
   std::string fname_pattern = "../../SWMF/run_cir_" + cir_date
                             + "/IH/IO2/3d__var_1_n00005000";

   simulation->AddBackground(BackgroundServerBATL(), container, fname_pattern);

//--------------------------------------------------------------------------------
// Time initial condition
//--------------------------------------------------------------------------------

   container.Clear();

// Initial time - set implicitly by the position of Earth each day according to SWMF
   double init_t = 0.0;
   container.Insert(init_t);

   simulation->AddInitial(InitialTimeFixed(), container);

//--------------------------------------------------------------------------------
// Spatial initial condition
//--------------------------------------------------------------------------------

   container.Clear();

// Initial position - set by the position of Earth
   std::string init_time = "0.0";
   if (argc > 4) init_time = argv[4];

   GeoVector earth_pos;
   if (MPI_Config::is_master) {
      std::string command_str = "python compute_earth_position.py " + cir_date + " --time " + init_time;
      const char* command_char = command_str.c_str();
      std::system(command_char);
      std::ifstream EarthPos_file;
      EarthPos_file.open("output_" + cir_date + "/earth_position_" + cir_date + ".dat");
      EarthPos_file >> earth_pos[0];
      EarthPos_file >> earth_pos[1];
      EarthPos_file >> earth_pos[2];
      EarthPos_file.close();
      std::cerr << "earth_pos = "
                << earth_pos
                << " au" << std::endl;
      for (int cpu = 1; cpu < MPI_Config::glob_comm_size; cpu++) {
         MPI_Send(earth_pos.Data(), 3, MPI_DOUBLE, cpu, 1001, MPI_Config::glob_comm);
      };
   }
   else {
      MPI_Recv(earth_pos.Data(), 3, MPI_DOUBLE, 0, 1001, MPI_Config::glob_comm, MPI_STATUS_IGNORE);
   };
   GeoVector init_pos = earth_pos * one_au;
   container.Insert(init_pos);

   simulation->AddInitial(InitialSpaceFixed(), container);

//--------------------------------------------------------------------------------
// Momentum initial condition
//--------------------------------------------------------------------------------

   container.Clear();

// Lower bound for momentum
   double momentum1 = Mom(100.0 * SPC_CONST_CGSM_MEGA_ELECTRON_VOLT / unit_energy_particle, specie);
   container.Insert(momentum1);

// Upper bound for momentum
   double momentum2 = Mom(1000.0 * SPC_CONST_CGSM_MEGA_ELECTRON_VOLT / unit_energy_particle, specie);
   container.Insert(momentum2);

// Log bias
   bool log_bias = true;
   container.Insert(log_bias);

   simulation->AddInitial(InitialMomentumThickShell(), container);

//--------------------------------------------------------------------------------
// Inner boundary
//--------------------------------------------------------------------------------

   container.Clear();

// Max crossings
   int max_crossings_Sun = 1;
   container.Insert(max_crossings_Sun);

// Action
   std::vector<int> actions_Sun;
   actions_Sun.push_back(-1);
   actions_Sun.push_back(-1);
   actions_Sun.push_back(-1);
   container.Insert(actions_Sun);

// Origin
   container.Insert(gv_zeros);

// Radius
   double inner_boundary = 0.1 * one_au;
   container.Insert(inner_boundary);

   simulation->AddBoundary(BoundarySphereReflect(), container);

//--------------------------------------------------------------------------------
// Outer boundary
//--------------------------------------------------------------------------------

   container.Clear();

// Max crossings
   int max_crossings_outer = 1;
   container.Insert(max_crossings_outer);

// Action
   std::vector<int> actions_outer;
   actions_outer.push_back(0);
   actions_outer.push_back(0);
   actions_outer.push_back(0);
   container.Insert(actions_outer);

// Origin
   container.Insert(gv_zeros);

// Radius
   double outer_boundary = 5.0 * one_au;
   container.Insert(outer_boundary);

   simulation->AddBoundary(BoundarySphereAbsorb(), container);

//--------------------------------------------------------------------------------
// Time limit
//--------------------------------------------------------------------------------

   container.Clear();

// Max crossings
   int max_crossings_time = 1;
   container.Insert(max_crossings_time);

// Action
   std::vector<int> actions_time;
   actions_time.push_back(-1);
   actions_time.push_back(-1);
   actions_time.push_back(-1);
   container.Insert(actions_time);
   
// Max duration of the trajectory
   double one_day = 24.0 * 60.0 * 60.0 / unit_time_fluid;
   double maxtime = -one_day * 365.0;
   container.Insert(maxtime);

   simulation->AddBoundary(BoundaryTimeExpire(), container);

//--------------------------------------------------------------------------------
// Diffusion model
//--------------------------------------------------------------------------------

   container.Clear();

#if CONST_DIFF < 3
// Index of region with forward propagating Alfven wave density
   int W_pls_idx = 1;
   container.Insert(W_pls_idx);

// Index of region with backward propagating Alfven wave density
   int W_mns_idx = 2;
   container.Insert(W_mns_idx);

// Constant = correlation_length * sqrt(B)
   double L_perp_times_sqrtB = 150.0 * (1.0e5 / unit_length_fluid) * sqrt(1.0e4 / unit_magnetic_fluid);
   container.Insert(L_perp_times_sqrtB);
#endif

#if CONST_DIFF > 0
   double kap0 = 1.0;
   double kap1 = 1.0;
   if (argc > 6) {
      kap0 = atof(argv[5]) / unit_diffusion_fluid;
      kap1 = atof(argv[6]) / unit_diffusion_fluid;
   };

// Base diffusion coefficient
#if CONST_DIFF == 1
   container.Insert(kap0);
#else
   container.Insert(kap1);
#endif

// Kinetic energy normalization factor
   double T0d = 0.2 * SPC_CONST_CGSM_GIGA_ELECTRON_VOLT / unit_energy_particle;
   container.Insert(T0d);

// Radius normalization factor
   container.Insert(one_au);

// Power law slope for kinetic energy
   double pow_law_Td = 1.0;
   container.Insert(pow_law_Td);

// Power law slope for radius
   double pow_law_r = 1.0;
   container.Insert(pow_law_r);

#if CONST_DIFF == 1
// Which coefficient to make empirical
   int which_kap = 0;
   container.Insert(which_kap);
#elif CONST_DIFF == 2
// Which coefficient to make empirical
   int which_kap = 1;
   container.Insert(which_kap);
#elif CONST_DIFF >= 3
// Ratio of perpendicular to parallel diffusion
   double kap_rat = kap0 / kap1;
   container.Insert(kap_rat);
#endif

#endif

// Pass ownership of "diffusion" to simulation
#if CONST_DIFF == 0
   simulation->AddDiffusion(DiffusionQLT_NLGC_AWSoM(), container);
#elif CONST_DIFF == 1 || CONST_DIFF == 2
   simulation->AddDiffusion(DiffusionQLT_or_NLGC_AWSoM(), container);
#elif CONST_DIFF == 3
   simulation->AddDiffusion(DiffusionKineticEnergyRadialDistancePowerLaw(), container);
#endif

   if (MPI_Config::is_master) {
#if CONST_DIFF > 0
      std::cerr << "CONSTANT DIFFUSION" << std::endl;
#endif
#if CONST_DIFF == 1 || CONST_DIFF == 3
      std::cerr << "k_perp = "
                << kap0 * unit_diffusion_fluid
                << " cm^2 / s" << std::endl;
#endif
#if CONST_DIFF == 2 || CONST_DIFF == 3
      std::cerr << "k_para = "
                << kap1 * unit_diffusion_fluid
                << " cm^2 / s" << std::endl;
#endif
   };

//--------------------------------------------------------------------------------
// Distribution 1 (spectrum)
//--------------------------------------------------------------------------------

   container.Clear();

// Number of bins
   MultiIndex n_bins1(100, 0, 0);
   container.Insert(n_bins1);
   
// Smallest value
   GeoVector minval1(EnrKin(momentum1), 0.0, 0.0);
   container.Insert(minval1);

// Largest value
   GeoVector maxval1(EnrKin(momentum2), 0.0, 0.0);
   container.Insert(maxval1);

// Linear or logarithmic bins
   MultiIndex log_bins1(1, 0, 0);
   container.Insert(log_bins1);

// Add outlying events to the end bins
   MultiIndex bin_outside1(0, 0, 0);
   container.Insert(bin_outside1);

// Physical units of the distro variable
   double unit_distro1 = 1.0;
   container.Insert(unit_distro1);

// Physical units of the bin variable
   GeoVector unit_val1 = {unit_energy_particle, 1.0, 1.0};
   container.Insert(unit_val1);

// Don't keep records
   bool keep_records1 = false;
   container.Insert(keep_records1);

// Normalization for the "hot" boundary
   double J0 = 18.6;
   container.Insert(J0);

//! Characteristic energy
   double T0 = SPC_CONST_CGSM_GIGA_ELECTRON_VOLT / unit_energy_particle;
   container.Insert(T0);

//! Spectral power law
   double pow_law_T = 0.987;
   container.Insert(pow_law_T);

//! Constant value for the "cold" condition
   double val_cold1 = 0.0;
   container.Insert(val_cold1);

//! Bendover energy
   double Tb = 1.04 * T0;
   container.Insert(Tb);

//! Spectral power law after bend
   double pow_law_Tb = -2.81;
   container.Insert(pow_law_Tb);

//! Smoothness of bend
   double bend_smooth = 3.97;
   container.Insert(bend_smooth);

   simulation->AddDistribution(DistributionSpectrumKineticEnergyBentPowerLaw(), container);

//--------------------------------------------------------------------------------
// Distribution 2 (exit time)
//--------------------------------------------------------------------------------

// Parameters for distribution
   container.Clear();

// Number of bins
   MultiIndex n_bins2(50, 0, 0);
   container.Insert(n_bins2);

// Smallest value
   GeoVector minval2(init_t - 10.0 * one_day, 0.0, 0.0);
   container.Insert(minval2);

// Largest value
   GeoVector maxval2(init_t, 0.0, 0.0);
   container.Insert(maxval2);

// Linear or logarithmic bins
   MultiIndex log_bins2(0, 0, 0);
   container.Insert(log_bins2);

// Add outlying events to the end bins
   MultiIndex bin_outside2(1, 0, 0);
   container.Insert(bin_outside2);

// Physical units of the distro variable
   double unit_distro2 = 1.0;
   container.Insert(unit_distro2);

// Physical units of the bin variable
   GeoVector unit_val2 = {unit_time_fluid, 1.0, 1.0};
   container.Insert(unit_val2);

// Keep records
   bool keep_records2 = false;
   container.Insert(keep_records2);

// Value for the "hot" condition
   double val_hot2 = 1.0;
   container.Insert(val_hot2);

// Value for the "cold" condition
   double val_cold2 = 0.0;
   container.Insert(val_cold2);

// Value for which time to bin (final)
   int val_time2 = 1;
   container.Insert(val_time2);

   simulation->AddDistribution(DistributionTimeUniform(), container);

//--------------------------------------------------------------------------------
// Distribution 3 (exit position)
//--------------------------------------------------------------------------------

// Parameters for distribution
   container.Clear();

// Number of bins
   MultiIndex n_bins3(1, 90, 180);
   container.Insert(n_bins3);

// Smallest value
   GeoVector minval3(outer_boundary - 0.1 * one_au, 0.0, 0.0);
   container.Insert(minval3);

// Largest value
   GeoVector maxval3(outer_boundary + 0.1 * one_au, M_PI, M_2PI);
   container.Insert(maxval3);

// Linear or logarithmic bins
   MultiIndex log_bins3(0, 0, 0);
   container.Insert(log_bins3);

// Add outlying events to the end bins
   MultiIndex bin_outside3(1, 0, 0);
   container.Insert(bin_outside3);

// Physical units of the distro variable
   double unit_distro3 = 1.0;
   container.Insert(unit_distro3);

// Physical units of the bin variable
   GeoVector unit_val3 = {unit_length_fluid, 1.0, 1.0};
   container.Insert(unit_val3);

// Keep records
   bool keep_records3 = false;
   container.Insert(keep_records3);

// Value for the "hot" condition
   double val_hot3 = 1.0;
   container.Insert(val_hot3);

// Value for the "cold" condition
   double val_cold3 = 0.0;
   container.Insert(val_cold3);

// Value for which time to bin (final)
   int val_time3 = 1;
   container.Insert(val_time3);

// Value for which coordinate system to use (spherical)
   int val_coord3 = 1;
   container.Insert(val_coord3);

   simulation->AddDistribution(DistributionPositionUniform(), container);

//--------------------------------------------------------------------------------
// Run the simulation
//--------------------------------------------------------------------------------

   int n_traj;
   int batch_size;

   batch_size = n_traj = 1;
   if(argc > 1) n_traj = atoi(argv[1]);
   if(argc > 2) batch_size = atoi(argv[2]);

   std::string simulation_files_prefix = "output_" + cir_date
                                       + "/GCR/cir_gcr_mod_"
                                       + init_time + "_distro_";
   simulation->DistroFileName(simulation_files_prefix);
   simulation->SetTasks(n_traj, batch_size);
   simulation->MainLoop();
   
   return 0;
};
