#include "src/server_config.hh"
#include "src/traj_config.hh"
#include "src/background_server_batl.hh"
#include "src/boundary_time.hh"
#include "src/boundary_space.hh"
#include "src/boundary_momentum.hh"
#include "src/initial_time.hh"
#include "src/initial_space.hh"
#include "src/initial_momentum.hh"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <filesystem>

using namespace Spectrum;

int main(int argc, char** argv)
{
   int active_local_workers, workers_stopped;
   std::string cir_date;

   std::shared_ptr<MPI_Config> mpi_config = std::make_shared<MPI_Config>(argc, argv);
   if (argc > 1) {
      cir_date = argv[1];
      if (mpi_config->is_master) {
         std::cout << "CIR date: " << cir_date << std::endl;
         std::filesystem::create_directory("output_" + cir_date);
         std::filesystem::create_directory("output_" + cir_date + "/fieldlines");
      };
   } else {
      if (mpi_config->is_master) std::cout << "ERROR: No CIR date provided." << std::endl;
      return 1;
   };
   MPI_Barrier(mpi_config->glob_comm);

//--------------------------------------------------------------------------------
// Server
//--------------------------------------------------------------------------------

   std::shared_ptr<ServerBaseBack> server_back = nullptr;

   std::string fname_pattern = "../../SWMF/run_cir_" + cir_date
                             + "/IH/IO2/3d__var_1_n00005000";

   if(mpi_config->is_boss) {
      server_back = std::make_unique<ServerBackType>(fname_pattern);
      active_local_workers = mpi_config->workers_in_node;
      server_back->ServerStart();
   };

   DataContainer container;

//--------------------------------------------------------------------------------
// Trajectory
//--------------------------------------------------------------------------------

   std::unique_ptr<TrajectoryBase> trajectory = std::make_unique<TrajectoryType>();
   
   std::shared_ptr<RNG> rng = std::make_shared<RNG>(time(NULL));
   trajectory->ConnectRNG(rng);
   
   int specie = SPECIES_PROTON_BEAM;
   trajectory->SetSpecie(specie);

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
   double dmax = GSL_CONST_CGSM_ASTRONOMICAL_UNIT / unit_length_fluid;
   container.Insert(dmax);

   trajectory->AddBackground(BackgroundServerCartesian(), container);

//--------------------------------------------------------------------------------
// Time initial condition
//--------------------------------------------------------------------------------

   container.Clear();

// Initial time
   double init_t = 0.0;
   container.Insert(init_t);

   trajectory->AddInitial(InitialTimeFixed(), container);

//--------------------------------------------------------------------------------
// Spatial initial condition
//--------------------------------------------------------------------------------

   container.Clear();

   fname_pattern = "data/earth_position_" + cir_date + ".dat";
   container.Insert(fname_pattern);

   double scale = GSL_CONST_CGSM_ASTRONOMICAL_UNIT / unit_length_fluid;
   container.Insert(scale);

   bool random = false;
   container.Insert(random);

   trajectory->AddInitial(InitialSpaceTable(), container);

//--------------------------------------------------------------------------------
// Momentum initial condition
//--------------------------------------------------------------------------------

   container.Clear();

// Initial momentum
   double GeV_kinetic_energy = 1000.0;
   container.Insert(Mom(GeV_kinetic_energy * SPC_CONST_CGSM_GIGA_ELECTRON_VOLT / unit_energy_particle, specie));

   trajectory->AddInitial(InitialMomentumShell(), container);

//--------------------------------------------------------------------------------
// Time limit
//--------------------------------------------------------------------------------

   container.Clear();

// Max crossings
   int max_crossings_time = 1;
   container.Insert(max_crossings_time);

// Action
   std::vector<int> actions_time;
   container.Insert(actions_time);
   
// Max duration of the trajectory
   double one_day = 60.0 * 60.0 * 24.0 / unit_time_fluid;
   double maxtime = 1.0 * one_day;
   container.Insert(maxtime);

   trajectory->AddBoundary(BoundaryTimeExpire(), container);

//--------------------------------------------------------------------------------
// Space boundary condition 1 (inner boundary)
//--------------------------------------------------------------------------------

   container.Clear();

// Max crossings
   int max_crossings = 1;
   container.Insert(max_crossings);

// Action
   std::vector<int> actions;
   container.Insert(actions);

// Origin
   container.Insert(gv_zeros);

// Radius
   double r_in = 0.1 * GSL_CONST_CGSM_ASTRONOMICAL_UNIT / unit_length_fluid;
   container.Insert(r_in);

   trajectory->AddBoundary(BoundarySphereAbsorb(), container);

//--------------------------------------------------------------------------------
// Space boundary condition 2 (outer boundary)
//--------------------------------------------------------------------------------

   container.Clear();

// Max crossings
   container.Insert(max_crossings);

// Action
   container.Insert(actions);

// Origin
   container.Insert(gv_zeros);

// Radius
   double r_out = 5.0 * GSL_CONST_CGSM_ASTRONOMICAL_UNIT / unit_length_fluid;
   container.Insert(r_out);

   trajectory->AddBoundary(BoundarySphereAbsorb(), container);

//--------------------------------------------------------------------------------

   if (mpi_config->is_boss) {
      while(active_local_workers) {
         workers_stopped = server_back->ServerFunctions();
         active_local_workers -= workers_stopped;
      };
      server_back->ServerFinish();
   }
   else if (mpi_config->is_worker) {
      std::string trajectory_file = "output_" + cir_date + "/fieldlines/fieldline_";
      int reps = 1;
      if (argc > 2) reps = atoi(argv[2]);
      for (int i = 0; i < reps; i++) {
         trajectory->SetStart();
         trajectory->Integrate();
         trajectory->InterpretStatus();
         trajectory->PrintCSV(trajectory_file + std::to_string(i) + ".lines", false);
         std::cout << "Field line outputted to " << trajectory_file + std::to_string(i) + ".lines" << std::endl;
      };
      trajectory->StopBackground();
   };

   std::cout << "CPU " << mpi_config->glob_comm_rank << " exited." << std::endl;

   return 0;
};


