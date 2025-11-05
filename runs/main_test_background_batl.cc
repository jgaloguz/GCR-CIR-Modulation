#include "src/server_config.hh"
#include "src/background_server_batl.hh"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <filesystem>

using namespace Spectrum;

inline GeoVector drift_numer(double r_L, double vel, SpatialData spdata, int specie)
{
   GeoVector drift = (r_L * vel / 3.0) * (spdata.curlB() - 2.0 * (spdata.gradBmag ^ spdata.bhat)) / spdata.Bmag;
// Correct magnitude if necessary
   if (drift.Norm() > 0.5 * vel) {
      drift.Normalize();
      drift *= 0.5 * vel;
   };
   return drift;
};

int main(int argc, char** argv)
{
   int active_local_workers, workers_stopped;
   BackgroundServerBATL background;

   SpatialData spdata;
   double t = 0.0;
   int i,j,k;
   GeoVector pos = gv_zeros, vel = gv_zeros, mom = gv_zeros;
   int specie = SPECIES_PROTON_BEAM;
   std::ofstream Den_file;
   std::ofstream AbsVel_file;
   std::ofstream DivVel_file;
   std::ofstream AbsMag_file;
   std::ofstream PolMag_file;
   std::ofstream dmax_file;
   std::ofstream drift_file;
   std::ofstream HetFlx_file;
   std::ofstream TurEnr_file;
   std::string cir_date;

   std::shared_ptr<MPI_Config> mpi_config = std::make_shared<MPI_Config>(argc, argv);
   if (argc > 1) {
      cir_date = argv[1];
      if (mpi_config->is_master) {
         std::cout << "CIR date: " << cir_date << std::endl;
         std::filesystem::create_directory("output_" + cir_date);
      };
   } else {
      if (mpi_config->is_master) std::cout << "ERROR: No CIR date provided." << std::endl;
      return 1;
   };
   MPI_Barrier(mpi_config->glob_comm);

//--------------------------------------------------------------------------------------------------
// Server
//--------------------------------------------------------------------------------------------------

   std::shared_ptr<ServerBaseBack> server_back = nullptr;

   std::string fname_pattern = "/data001/cosmicrays_vf/Juan/SWMF/run_cir_"
                             + cir_date + "/IH/IO2/3d__var_1_n00005000";

   if (mpi_config->is_boss) {
      server_back = std::make_unique<ServerBackType>(fname_pattern);
      active_local_workers = mpi_config->workers_in_node;
      server_back->ServerStart();
   };

   DataContainer container;

//--------------------------------------------------------------------------------------------------
// Background
//--------------------------------------------------------------------------------------------------

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

   background.SetupObject(container);
   background.SetSpecie(specie);

//--------------------------------------------------------------------------------------------------

   if (mpi_config->is_boss) {
      while(active_local_workers) {
         workers_stopped = server_back->ServerFunctions();
         active_local_workers -= workers_stopped;
      };
      server_back->ServerFinish();
   }
   else if (mpi_config->is_worker) {

      int i, j, k, N = 1000;
      double Rs = 6.957e+10 / unit_length_fluid;
      double x_min = -1075.0 * Rs;
      double y_min = -1075.0 * Rs;
      double z_min = -1075.0 * Rs;
      double dx = 2150.0 * Rs / (N-1);
      double dy = 2150.0 * Rs / (N-1);
      double dz = 2150.0 * Rs / (N-1);
      double polarity, r_L;
      GeoVector drift_vel;
      spdata._mask = BACKGROUND_ALL | BACKGROUND_gradU | BACKGROUND_gradB;

      double one_au = GSL_CONST_CGSM_ASTRONOMICAL_UNIT / unit_length_fluid;
      double one_day = 24.0 * 60.0 * 60.0 / unit_time_fluid;
      pos[0] = 1.0 * one_au;
      mom[0] = Mom(1000.0 * SPC_CONST_CGSM_MEGA_ELECTRON_VOLT / unit_energy_particle, specie);
      vel[0] = Vel(mom[0], specie);
      background.GetFields(t, pos, mom, spdata);
      r_L = LarmorRadius(mom[0], spdata.Bmag, specie);
      std::cout << "|B| @ 1au = "
                << std::setw(18) << spdata.Bmag * unit_magnetic_fluid * 1.0e5 << " nT"
                << std::endl;
      std::cout << "rL (1 GeV) = "
                << std::setw(18) << r_L << " au"
                << std::endl;
      std::cout << "dmax @ 1 au = "
                << std::setw(18) << spdata.dmax << " au"
                << std::endl;

//--------------------------------------------------------------------------------------------------

      std::cout << "2D plots..." << std::endl;

      Den_file.open("output_" + cir_date + "/den_equ_" + cir_date + ".dat");
      AbsVel_file.open("output_" + cir_date + "/vel_equ_" + cir_date + ".dat");
      DivVel_file.open("output_" + cir_date + "/div_vel_equ_" + cir_date + ".dat");
      AbsMag_file.open("output_" + cir_date + "/mag_equ_" + cir_date + ".dat");
      PolMag_file.open("output_" + cir_date + "/pol_equ_" + cir_date + ".dat");
      dmax_file.open("output_" + cir_date + "/dmax_equ_" + cir_date + ".dat");
      drift_file.open("output_" + cir_date + "/drift_equ_" + cir_date + ".dat");
      HetFlx_file.open("output_" + cir_date + "/het_flx_equ_" + cir_date + ".dat");
      TurEnr_file.open("output_" + cir_date + "/tur_enr_equ_" + cir_date + ".dat");
      
      pos = gv_zeros;
      for (i = 0; i < N; i++) {
         pos[0] = x_min + i * dx;
         for (j = 0; j < N; j++) {
            pos[1] = y_min + j * dy;
            if (pos.Norm() < 20.0 * Rs) {
               spdata.n_dens = 0.0;
               spdata.Uvec = gv_zeros;
               spdata.gradUvec = gm_zeros;
               spdata.Bmag = 0.0;
               polarity = 0.0;
               spdata.dmax = 0.0;
               drift_vel = gv_zeros;
               spdata.region = gv_zeros;
            }
            else {
               background.GetFields(t, pos, mom, spdata);
               polarity = (spdata.Bvec * pos >= 0.0 ? 1.0 : -1.0);
               r_L = LarmorRadius(mom[0], spdata.Bmag, specie);
               drift_vel = drift_numer(r_L, vel[0], spdata, specie);
            };
            Den_file << std::setw(18) << spdata.n_dens * unit_density_fluid;
            AbsVel_file << std::setw(18) << spdata.Uvec.Norm() * unit_velocity_fluid;
            DivVel_file << std::setw(18) << spdata.divU() * unit_velocity_fluid / unit_length_fluid;
            AbsMag_file << std::setw(18) << spdata.Bmag * unit_magnetic_fluid;
            PolMag_file << std::setw(18) << polarity;
            dmax_file << std::setw(18) << spdata.dmax;
            drift_file << std::setw(18) << drift_vel.Norm() / vel[0];
            HetFlx_file << std::setw(18) << spdata.region[0];
            TurEnr_file << std::setw(18) << spdata.region[1]+spdata.region[2];
         };
         Den_file << std::endl;
         AbsVel_file << std::endl;
         DivVel_file << std::endl;
         AbsMag_file << std::endl;
         PolMag_file << std::endl;
         dmax_file << std::endl;
         drift_file << std::endl;
         HetFlx_file << std::endl;
         TurEnr_file << std::endl;
      };
      Den_file.close();
      AbsVel_file.close();
      DivVel_file.close();
      AbsMag_file.close();
      PolMag_file.close();
      dmax_file.close();
      drift_file.close();
      HetFlx_file.close();
      TurEnr_file.close();

      Den_file.open("output_" + cir_date + "/den_mer_" + cir_date + ".dat");
      AbsVel_file.open("output_" + cir_date + "/vel_mer_" + cir_date + ".dat");
      DivVel_file.open("output_" + cir_date + "/div_vel_mer_" + cir_date + ".dat");
      AbsMag_file.open("output_" + cir_date + "/mag_mer_" + cir_date + ".dat");
      PolMag_file.open("output_" + cir_date + "/pol_mer_" + cir_date + ".dat");
      dmax_file.open("output_" + cir_date + "/dmax_mer_" + cir_date + ".dat");
      drift_file.open("output_" + cir_date + "/drift_mer_" + cir_date + ".dat");
      HetFlx_file.open("output_" + cir_date + "/het_flx_mer_" + cir_date + ".dat");
      TurEnr_file.open("output_" + cir_date + "/tur_enr_mer_" + cir_date + ".dat");
      
      pos = gv_zeros;
      for (i = 0; i < N; i++) {
         pos[0] = x_min + i * dx;
         for (k = 0; k < N; k++) {
            pos[2] = z_min + k * dz;
            if (pos.Norm() < 20.0 * Rs) {
               spdata.n_dens = 0.0;
               spdata.Uvec = gv_zeros;
               spdata.gradUvec = gm_zeros;
               spdata.Bmag = 0.0;
               polarity = 0.0;
               spdata.dmax = 0.0;
               drift_vel = gv_zeros;
               spdata.region = gv_zeros;
            }
            else {
               background.GetFields(t, pos, mom, spdata);
               polarity = (spdata.Bvec * pos >= 0.0 ? 1.0 : -1.0);
               r_L = LarmorRadius(mom[0], spdata.Bmag, specie);
               drift_vel = drift_numer(r_L, vel[0], spdata, specie);
            };
            Den_file << std::setw(18) << spdata.n_dens * unit_density_fluid;
            AbsVel_file << std::setw(18) << spdata.Uvec.Norm() * unit_velocity_fluid;
            DivVel_file << std::setw(18) << spdata.divU() * unit_velocity_fluid / unit_length_fluid;
            AbsMag_file << std::setw(18) << spdata.Bmag * unit_magnetic_fluid;
            PolMag_file << std::setw(18) << polarity;
            dmax_file << std::setw(18) << spdata.dmax;
            drift_file << std::setw(18) << drift_vel.Norm() / vel[0];
            HetFlx_file << std::setw(18) << spdata.region[0];
            TurEnr_file << std::setw(18) << spdata.region[1]+spdata.region[2];
         };
         Den_file << std::endl;
         AbsVel_file << std::endl;
         DivVel_file << std::endl;
         AbsMag_file << std::endl;
         PolMag_file << std::endl;
         dmax_file << std::endl;
         drift_file << std::endl;
         HetFlx_file << std::endl;
         TurEnr_file << std::endl;
      };
      Den_file.close();
      AbsVel_file.close();
      DivVel_file.close();
      AbsMag_file.close();
      PolMag_file.close();
      dmax_file.close();
      drift_file.close();
      HetFlx_file.close();
      TurEnr_file.close();

//--------------------------------------------------------------------------------------------------

      std::cout << "1D plots..." << std::endl;

      Den_file.open("output_" + cir_date + "/den_1au_" + cir_date + ".dat");
      AbsVel_file.open("output_" + cir_date + "/vel_1au_" + cir_date + ".dat");
      DivVel_file.open("output_" + cir_date + "/div_vel_1au_" + cir_date + ".dat");
      AbsMag_file.open("output_" + cir_date + "/mag_1au_" + cir_date + ".dat");
      PolMag_file.open("output_" + cir_date + "/pol_1au_" + cir_date + ".dat");
      dmax_file.open("output_" + cir_date + "/dmax_1au_" + cir_date + ".dat");
      drift_file.open("output_" + cir_date + "/drift_1au_" + cir_date + ".dat");
      HetFlx_file.open("output_" + cir_date + "/het_flx_1au_" + cir_date + ".dat");
      TurEnr_file.open("output_" + cir_date + "/tur_enr_1au_" + cir_date + ".dat");

      pos[0] = one_au * cos(M_PI_4);
      pos[1] = one_au * sin(M_PI_4);
      pos[2] = 0.0;
      double dt = 27.0 * one_day / N;
      for (i = 0; i < N; i++) {
// Frame rotates CCW in the xy-plane, so steady-state data should be sampled CW
         t = i * dt;

         background.GetFields(t, pos, mom, spdata);
         polarity = (spdata.Bvec * pos >= 0.0 ? 1.0 : -1.0);
         r_L = LarmorRadius(mom[0], spdata.Bmag, specie);
         drift_vel = drift_numer(r_L, vel[0], spdata, specie);

         Den_file << std::setw(18) << t / one_day
                  << std::setw(18) << spdata.n_dens * unit_density_fluid
                  << std::endl;
         AbsVel_file << std::setw(18) << t / one_day 
                     << std::setw(18) << spdata.Uvec.Norm() * unit_velocity_fluid
                     << std::endl;
         DivVel_file << std::setw(18) << t / one_day 
                     << std::setw(18) << spdata.divU() * unit_velocity_fluid / unit_length_fluid
                     << std::endl;
         AbsMag_file << std::setw(18) << t / one_day 
                     << std::setw(18) << spdata.Bmag * unit_magnetic_fluid
                     << std::endl;
         PolMag_file << std::setw(18) << t / one_day 
                     << std::setw(18) << polarity
                     << std::endl;
         dmax_file << std::setw(18) << t / one_day 
                   << std::setw(18) << spdata.dmax
                   << std::endl;
         drift_file << std::setw(18) << t / one_day 
                    << std::setw(18) << drift_vel.Norm() / vel[0]
                    << std::endl;
         HetFlx_file << std::setw(18) << t / one_day 
                     << std::setw(18) << spdata.region[0]
                     << std::endl;
         TurEnr_file << std::setw(18) << t / one_day 
                     << std::setw(18) << spdata.region[1]+spdata.region[2]
                     << std::endl;
      };

      Den_file.close();
      AbsVel_file.close();
      DivVel_file.close();
      AbsMag_file.close();
      PolMag_file.close();
      dmax_file.close();
      drift_file.close();
      HetFlx_file.close();
      TurEnr_file.close();

      background.StopServerFront();
      std::cout << "Background samples outputted to file." << std::endl;
   };

   std::cout << "Node " << mpi_config->glob_comm_rank << " exited." << std::endl;
   return 0;
};


