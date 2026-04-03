#include "src/server_config.hh"
#include "src/background_server_batl.hh"
#include "src/diffusion_other.hh"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <filesystem>

using namespace Spectrum;

const int specie = SPECIES_PROTON_BEAM;
const double one_au = GSL_CONST_CGSM_ASTRONOMICAL_UNIT / unit_length_fluid;
const double one_day = 24.0 * 60.0 * 60.0 / unit_time_fluid;
const double Rs = 6.957e+10 / unit_length_fluid;

inline GeoVector drift_numer(double r_L, double vel, SpatialData spdata)
{
   GeoVector drift = (r_L * vel / 3.0) * (spdata.curlB() - 2.0 * (spdata.gradBmag ^ spdata.bhat)) / spdata.Bmag;
// Correct magnitude if necessary
   if (drift.Norm() > 0.5 * vel) {
      drift.Normalize();
      drift *= 0.5 * vel;
   };
   return drift;
};

inline double correct_Z2(SpatialData spdata, double r)
{
   double Z2 = spdata.region[1] + spdata.region[2];
#if defined(DIFF_CORRECT_Z2)
   Z2 += Z2_diff * (one_au / r);
#endif
   return Z2 / SPC_CONST_CGSM_MASS_PROTON;
};

void print_1D_cuts(MultiIndex direction, std::string cir_date, BackgroundServerBATL *background, DiffusionQLT_NLGC_AWSoM *diffusion, SpatialData spdata, double mom0)
{
   int i, N = 1000;
   double t = 0.0, l0 = 20.0 * Rs, l1 = 5.0 * one_au, dd;
   double rad_vel, polarity, r_L, diff1, diff2;
   GeoVector pos, dir, drift_vel;
   GeoVector mom = gv_zeros, vel = gv_zeros;
   mom[0] = mom0;
   vel[0] = Vel(mom[0], specie);

   std::ofstream Den_file;
   std::ofstream Pth_file;
   std::ofstream AbsVel_file;
   std::ofstream RadVel_file;
   std::ofstream DivVel_file;
   std::ofstream AbsMag_file;
   std::ofstream PolMag_file;
   std::ofstream dmax_file;
   std::ofstream HetFlx_file;
   std::ofstream TurEnr_file;
   std::ofstream drift_file;
   std::ofstream diff1_file;
   std::ofstream diff2_file;

   std::string dir_string = std::to_string(direction.i) + std::to_string(direction.j) + std::to_string(direction.k);
   std::cout << "Direction: " << dir_string << std::endl;

   Den_file.open("output_" + cir_date + "/CIR/den_" + dir_string + "_" + cir_date + ".dat");
      Pth_file.open("output_" + cir_date + "/CIR/pth_" + dir_string + "_" + cir_date + ".dat");
      AbsVel_file.open("output_" + cir_date + "/CIR/vel_" + dir_string + "_" + cir_date + ".dat");
      RadVel_file.open("output_" + cir_date + "/CIR/rad_vel_" + dir_string + "_" + cir_date + ".dat");
      DivVel_file.open("output_" + cir_date + "/CIR/div_vel_" + dir_string + "_" + cir_date + ".dat");
      AbsMag_file.open("output_" + cir_date + "/CIR/mag_" + dir_string + "_" + cir_date + ".dat");
      PolMag_file.open("output_" + cir_date + "/CIR/pol_" + dir_string + "_" + cir_date + ".dat");
      dmax_file.open("output_" + cir_date + "/CIR/dmax_" + dir_string + "_" + cir_date + ".dat");
      HetFlx_file.open("output_" + cir_date + "/CIR/het_flx_" + dir_string + "_" + cir_date + ".dat");
      TurEnr_file.open("output_" + cir_date + "/CIR/tur_enr_" + dir_string + "_" + cir_date + ".dat");
      drift_file.open("output_" + cir_date + "/CIR/drift_" + dir_string + "_" + cir_date + ".dat");
      diff1_file.open("output_" + cir_date + "/CIR/diff1_" + dir_string + "_" + cir_date + ".dat");
      diff2_file.open("output_" + cir_date + "/CIR/diff2_" + dir_string + "_" + cir_date + ".dat");

      dir.x = direction.i;
      dir.y = direction.j;
      dir.z = direction.k;
      dir.Normalize();
      dd = (log(l1) - log(l0)) / (N - 1);
      double dt = 0.0;
      for (i = 0; i < N; i++) {
         pos = dir * exp(log(l0) + dd * i);

         background->GetFields(t, pos, mom, spdata);
         rad_vel = UnitVec(pos) * UnitVec(spdata.Uvec);
         polarity = (spdata.Bvec * pos >= 0.0 ? 1.0 : -1.0);
         r_L = LarmorRadius(mom[0], spdata.Bmag, specie);
         drift_vel = drift_numer(r_L, vel[0], spdata);
         diff1 = diffusion->GetComponent(0, t, pos, mom, spdata);
         diff2 = diffusion->GetComponent(1, t, pos, mom, spdata);   

         Den_file << std::setw(18) << pos.Norm() / one_au
                  << std::setw(18) << spdata.n_dens * unit_number_density_fluid
                  << std::endl;
         Pth_file << std::setw(18) << pos.Norm() / one_au
                  << std::setw(18) << spdata.p_ther * unit_pressure_fluid
                  << std::endl;
         AbsVel_file << std::setw(18) << pos.Norm() / one_au 
                     << std::setw(18) << spdata.Uvec.Norm() * unit_velocity_fluid
                     << std::endl;
         RadVel_file << std::setw(18) << pos.Norm() / one_au 
                     << std::setw(18) << rad_vel
                     << std::endl;
         DivVel_file << std::setw(18) << pos.Norm() / one_au 
                     << std::setw(18) << spdata.divU() * unit_velocity_fluid / unit_length_fluid
                     << std::endl;
         AbsMag_file << std::setw(18) << pos.Norm() / one_au 
                     << std::setw(18) << spdata.Bmag * unit_magnetic_fluid
                     << std::endl;
         PolMag_file << std::setw(18) << pos.Norm() / one_au 
                     << std::setw(18) << polarity
                     << std::endl;
         dmax_file << std::setw(18) << pos.Norm() / one_au 
                   << std::setw(18) << spdata.dmax
                   << std::endl;
         HetFlx_file << std::setw(18) << pos.Norm() / one_au 
                     << std::setw(18) << spdata.region[0]
                     << std::endl;
         TurEnr_file << std::setw(18) << pos.Norm() / one_au 
                     << std::setw(18) << correct_Z2(spdata, pos.Norm())
                     << std::endl;
         drift_file << std::setw(18) << pos.Norm() / one_au 
                    << std::setw(18) << drift_vel.Norm() / vel[0]
                    << std::endl;
         diff1_file << std::setw(18) << pos.Norm() / one_au 
                    << std::setw(18) << diff1 * unit_diffusion_fluid
                    << std::endl;
         diff2_file << std::setw(18) << pos.Norm() / one_au 
                    << std::setw(18) << diff2 * unit_diffusion_fluid
                    << std::endl;
      };

      Den_file.close();
      Pth_file.close();
      AbsVel_file.close();
      RadVel_file.close();
      DivVel_file.close();
      AbsMag_file.close();
      PolMag_file.close();
      dmax_file.close();
      HetFlx_file.close();
      TurEnr_file.close();
      drift_file.close();
      diff1_file.close();
      diff2_file.close();
};

int main(int argc, char** argv)
{
   int active_local_workers, workers_stopped;
   BackgroundServerBATL background;
   DiffusionQLT_NLGC_AWSoM diffusion;

   SpatialData spdata;
   double t = 0.0;
   int i,j,k;
   GeoVector pos = gv_zeros, vel = gv_zeros, mom = gv_zeros;
   std::ofstream Den_file;
   std::ofstream Pth_file;
   std::ofstream AbsVel_file;
   std::ofstream RadVel_file;
   std::ofstream DivVel_file;
   std::ofstream AbsMag_file;
   std::ofstream PolMag_file;
   std::ofstream dmax_file;
   std::ofstream HetFlx_file;
   std::ofstream TurEnr_file;
   std::ofstream drift_file;
   std::ofstream diff1_file;
   std::ofstream diff2_file;
   std::string cir_date;

   std::shared_ptr<MPI_Config> mpi_config = std::make_shared<MPI_Config>(argc, argv);
   if (argc > 1) {
      cir_date = argv[1];
      if (mpi_config->is_master) {
         std::cout << "CIR date: " << cir_date << std::endl;
         std::filesystem::create_directory("output_" + cir_date);
         std::filesystem::create_directory("output_" + cir_date + "/CIR");
         std::cout << "Number of servers: " << mpi_config->boss_comm_size-1 << std::endl;
         std::cout << "Number of workers: " << mpi_config->n_workers << std::endl;
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

   std::string fname_pattern = "/data001/cosmicrays_vf/Juan/SWMF/run_cir_"
                             + cir_date + "/IH/IO2/3d__var_1_n00005000";

   if (mpi_config->is_boss) {
      server_back = std::make_unique<ServerBackType>(fname_pattern);
      active_local_workers = mpi_config->workers_in_node;
      server_back->ServerStart();
   };

   DataContainer container;

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

   background.SetupObject(container);
   background.SetSpecie(specie);

//--------------------------------------------------------------------------------
// Diffusion model
//--------------------------------------------------------------------------------

   container.Clear();

// Index of region with forward propagating Alfven wave density
   int W_pls_idx = 1;
   container.Insert(W_pls_idx);

// Index of region with backward propagating Alfven wave density
   int W_mns_idx = 2;
   container.Insert(W_mns_idx);

// Constant = correlation_length * sqrt(B)
   double L_perp_times_sqrtB = 150.0 * (1.0e5 / unit_length_fluid) * sqrt(1.0e4 / unit_magnetic_fluid);
   container.Insert(L_perp_times_sqrtB);

// Set up diffusion object
   diffusion.SetupObject(container);
   diffusion.SetSpecie(specie);

//--------------------------------------------------------------------------------

   if (mpi_config->is_boss) {
      while(active_local_workers) {
         workers_stopped = server_back->ServerFunctions();
         active_local_workers -= workers_stopped;
      };
      server_back->ServerFinish();
   }
   else if (mpi_config->is_worker) {

      int i, j, k, N = 1000;
      double x_min = -1075.0 * Rs;
      double y_min = -1075.0 * Rs;
      double z_min = -1075.0 * Rs;
      double dx = 2150.0 * Rs / (N-1);
      double dy = 2150.0 * Rs / (N-1);
      double dz = 2150.0 * Rs / (N-1);
      double rad_vel;
      double polarity, r_L;
      double diff1, diff2;
      GeoVector drift_vel;
      spdata._mask = BACKGROUND_ALL | BACKGROUND_gradU | BACKGROUND_gradB;

      pos[0] = 1.0 * one_au;
      mom[0] = Mom(200.0 * SPC_CONST_CGSM_MEGA_ELECTRON_VOLT / unit_energy_particle, specie);
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

//--------------------------------------------------------------------------------

      std::cout << "1D plots..." << std::endl;

      Den_file.open("output_" + cir_date + "/CIR/den_1au_" + cir_date + ".dat");
      Pth_file.open("output_" + cir_date + "/CIR/pth_1au_" + cir_date + ".dat");
      AbsVel_file.open("output_" + cir_date + "/CIR/vel_1au_" + cir_date + ".dat");
      RadVel_file.open("output_" + cir_date + "/CIR/rad_vel_1au_" + cir_date + ".dat");
      DivVel_file.open("output_" + cir_date + "/CIR/div_vel_1au_" + cir_date + ".dat");
      AbsMag_file.open("output_" + cir_date + "/CIR/mag_1au_" + cir_date + ".dat");
      PolMag_file.open("output_" + cir_date + "/CIR/pol_1au_" + cir_date + ".dat");
      dmax_file.open("output_" + cir_date + "/CIR/dmax_1au_" + cir_date + ".dat");
      HetFlx_file.open("output_" + cir_date + "/CIR/het_flx_1au_" + cir_date + ".dat");
      TurEnr_file.open("output_" + cir_date + "/CIR/tur_enr_1au_" + cir_date + ".dat");
      drift_file.open("output_" + cir_date + "/CIR/drift_1au_" + cir_date + ".dat");
      diff1_file.open("output_" + cir_date + "/CIR/diff1_1au_" + cir_date + ".dat");
      diff2_file.open("output_" + cir_date + "/CIR/diff2_1au_" + cir_date + ".dat");

      pos[0] = one_au;
      pos[1] = 0.0;
      pos[2] = 0.0;
      double dt = 27.0 * one_day / N;
      for (i = 0; i < N; i++) {
// Frame rotates CCW in the xy-plane, so steady-state data should be sampled CW
         t = i * dt;

         background.GetFields(t, pos, mom, spdata);
         rad_vel = UnitVec(pos) * UnitVec(spdata.Uvec);
         polarity = (spdata.Bvec * pos >= 0.0 ? 1.0 : -1.0);
         r_L = LarmorRadius(mom[0], spdata.Bmag, specie);
         drift_vel = drift_numer(r_L, vel[0], spdata);
         diff1 = diffusion.GetComponent(0, t, pos, mom, spdata);
         diff2 = diffusion.GetComponent(1, t, pos, mom, spdata);   

         Den_file << std::setw(18) << t / one_day
                  << std::setw(18) << spdata.n_dens * unit_number_density_fluid
                  << std::endl;
         Pth_file << std::setw(18) << t / one_day
                  << std::setw(18) << spdata.p_ther * unit_pressure_fluid
                  << std::endl;
         AbsVel_file << std::setw(18) << t / one_day 
                     << std::setw(18) << spdata.Uvec.Norm() * unit_velocity_fluid
                     << std::endl;
         RadVel_file << std::setw(18) << t / one_day 
                     << std::setw(18) << rad_vel
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
         HetFlx_file << std::setw(18) << t / one_day 
                     << std::setw(18) << spdata.region[0]
                     << std::endl;
         TurEnr_file << std::setw(18) << t / one_day 
                     << std::setw(18) << correct_Z2(spdata, one_au)
                     << std::endl;
         drift_file << std::setw(18) << t / one_day 
                    << std::setw(18) << drift_vel.Norm() / vel[0]
                    << std::endl;
         diff1_file << std::setw(18) << t / one_day 
                    << std::setw(18) << diff1 * unit_diffusion_fluid
                    << std::endl;
         diff2_file << std::setw(18) << t / one_day 
                    << std::setw(18) << diff2 * unit_diffusion_fluid
                    << std::endl;
      };

      Den_file.close();
      Pth_file.close();
      AbsVel_file.close();
      RadVel_file.close();
      DivVel_file.close();
      AbsMag_file.close();
      PolMag_file.close();
      dmax_file.close();
      HetFlx_file.close();
      TurEnr_file.close();
      drift_file.close();
      diff1_file.close();
      diff2_file.close();

      MultiIndex direction;
      direction[0] = 1; direction[1] = 0; direction[2] = 0;
      print_1D_cuts(direction, cir_date, &background, &diffusion, spdata, mom[0]);
      direction[0] = -1; direction[1] = 0; direction[2] = 0;
      print_1D_cuts(direction, cir_date, &background, &diffusion, spdata, mom[0]);
      direction[0] = 0; direction[1] = 1; direction[2] = 0;
      print_1D_cuts(direction, cir_date, &background, &diffusion, spdata, mom[0]);
      direction[0] = 0; direction[1] = -1; direction[2] = 0;
      print_1D_cuts(direction, cir_date, &background, &diffusion, spdata, mom[0]);
      direction[0] = 0; direction[1] = 0; direction[2] = 1;
      print_1D_cuts(direction, cir_date, &background, &diffusion, spdata, mom[0]);
      direction[0] = 0; direction[1] = 0; direction[2] = -1;
      print_1D_cuts(direction, cir_date, &background, &diffusion, spdata, mom[0]);

      background.StopServerFront();
      std::cout << "Background samples outputted to file." << std::endl;
   };

   std::cout << "CPU " << mpi_config->glob_comm_rank << " exited." << std::endl;
   return 0;
};


