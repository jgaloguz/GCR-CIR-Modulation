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

      double one_au = GSL_CONST_CGSM_ASTRONOMICAL_UNIT / unit_length_fluid;
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

      std::cout << "2D plots..." << std::endl;

      Den_file.open("output_" + cir_date + "/CIR/den_equ_" + cir_date + ".dat");
      Pth_file.open("output_" + cir_date + "/CIR/pth_equ_" + cir_date + ".dat");
      AbsVel_file.open("output_" + cir_date + "/CIR/vel_equ_" + cir_date + ".dat");
      RadVel_file.open("output_" + cir_date + "/CIR/rad_vel_equ_" + cir_date + ".dat");
      DivVel_file.open("output_" + cir_date + "/CIR/div_vel_equ_" + cir_date + ".dat");
      AbsMag_file.open("output_" + cir_date + "/CIR/mag_equ_" + cir_date + ".dat");
      PolMag_file.open("output_" + cir_date + "/CIR/pol_equ_" + cir_date + ".dat");
      dmax_file.open("output_" + cir_date + "/CIR/dmax_equ_" + cir_date + ".dat");
      HetFlx_file.open("output_" + cir_date + "/CIR/het_flx_equ_" + cir_date + ".dat");
      TurEnr_file.open("output_" + cir_date + "/CIR/tur_enr_equ_" + cir_date + ".dat");
      drift_file.open("output_" + cir_date + "/CIR/drift_equ_" + cir_date + ".dat");
      diff1_file.open("output_" + cir_date + "/CIR/diff1_equ_" + cir_date + ".dat");
      diff2_file.open("output_" + cir_date + "/CIR/diff2_equ_" + cir_date + ".dat");
      
      pos = gv_zeros;
      for (i = 0; i < N; i++) {
         pos[0] = x_min + i * dx;
         for (j = 0; j < N; j++) {
            pos[1] = y_min + j * dy;
            if (pos.Norm() < 20.0 * Rs) {
               spdata.n_dens = 0.0;
               spdata.p_ther = 0.0;
               spdata.Uvec = gv_zeros;
               rad_vel = 0.0;
               spdata.gradUvec = gm_zeros;
               spdata.Bmag = 0.0;
               polarity = 0.0;
               spdata.dmax = 0.0;
               spdata.region = gv_zeros;
               drift_vel = gv_zeros;
               diff1 = 0.0;
               diff2 = 0.0;
            }
            else {
               background.GetFields(t, pos, mom, spdata);
               rad_vel = UnitVec(pos) * UnitVec(spdata.Uvec);
               polarity = (spdata.Bvec * pos >= 0.0 ? 1.0 : -1.0);
               r_L = LarmorRadius(mom[0], spdata.Bmag, specie);
               drift_vel = drift_numer(r_L, vel[0], spdata);
               diff1 = diffusion.GetComponent(0, t, pos, mom, spdata);
               diff2 = diffusion.GetComponent(1, t, pos, mom, spdata);
            };
            Den_file << std::setw(18) << spdata.n_dens * unit_number_density_fluid;
            Pth_file << std::setw(18) << spdata.p_ther * unit_pressure_fluid;
            AbsVel_file << std::setw(18) << spdata.Uvec.Norm() * unit_velocity_fluid;
            RadVel_file << std::setw(18) << rad_vel;
            DivVel_file << std::setw(18) << spdata.divU() * unit_velocity_fluid / unit_length_fluid;
            AbsMag_file << std::setw(18) << spdata.Bmag * unit_magnetic_fluid;
            PolMag_file << std::setw(18) << polarity;
            dmax_file << std::setw(18) << spdata.dmax;
            HetFlx_file << std::setw(18) << spdata.region[0];
            TurEnr_file << std::setw(18) << correct_Z2(spdata, pos.Norm());
            drift_file << std::setw(18) << drift_vel.Norm() / vel[0];
            diff1_file << std::setw(18) << diff1 * unit_diffusion_fluid;
            diff2_file << std::setw(18) << diff2 * unit_diffusion_fluid;
         };
         Den_file << std::endl;
         Pth_file << std::endl;
         AbsVel_file << std::endl;
         RadVel_file << std::endl;
         DivVel_file << std::endl;
         AbsMag_file << std::endl;
         PolMag_file << std::endl;
         dmax_file << std::endl;
         HetFlx_file << std::endl;
         TurEnr_file << std::endl;
         drift_file << std::endl;
         diff1_file << std::endl;
         diff2_file << std::endl;
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

      Den_file.open("output_" + cir_date + "/CIR/den_mer_" + cir_date + ".dat");
      Pth_file.open("output_" + cir_date + "/CIR/pth_mer_" + cir_date + ".dat");
      AbsVel_file.open("output_" + cir_date + "/CIR/vel_mer_" + cir_date + ".dat");
      RadVel_file.open("output_" + cir_date + "/CIR/rad_vel_mer_" + cir_date + ".dat");
      DivVel_file.open("output_" + cir_date + "/CIR/div_vel_mer_" + cir_date + ".dat");
      AbsMag_file.open("output_" + cir_date + "/CIR/mag_mer_" + cir_date + ".dat");
      PolMag_file.open("output_" + cir_date + "/CIR/pol_mer_" + cir_date + ".dat");
      dmax_file.open("output_" + cir_date + "/CIR/dmax_mer_" + cir_date + ".dat");
      HetFlx_file.open("output_" + cir_date + "/CIR/het_flx_mer_" + cir_date + ".dat");
      TurEnr_file.open("output_" + cir_date + "/CIR/tur_enr_mer_" + cir_date + ".dat");
      drift_file.open("output_" + cir_date + "/CIR/drift_mer_" + cir_date + ".dat");
      diff1_file.open("output_" + cir_date + "/CIR/diff1_mer_" + cir_date + ".dat");
      diff2_file.open("output_" + cir_date + "/CIR/diff2_mer_" + cir_date + ".dat");
      
      pos = gv_zeros;
      for (i = 0; i < N; i++) {
         pos[0] = x_min + i * dx;
         for (k = 0; k < N; k++) {
            pos[2] = z_min + k * dz;
            if (pos.Norm() < 20.0 * Rs) {
               spdata.n_dens = 0.0;
               spdata.p_ther = 0.0;
               spdata.Uvec = gv_zeros;
               rad_vel = 0.0;
               spdata.gradUvec = gm_zeros;
               spdata.Bmag = 0.0;
               polarity = 0.0;
               spdata.dmax = 0.0;
               spdata.region = gv_zeros;
               drift_vel = gv_zeros;
               diff1 = 0.0;
               diff2 = 0.0;
            }
            else {
               background.GetFields(t, pos, mom, spdata);
               rad_vel = UnitVec(pos) * UnitVec(spdata.Uvec);
               polarity = (spdata.Bvec * pos >= 0.0 ? 1.0 : -1.0);
               r_L = LarmorRadius(mom[0], spdata.Bmag, specie);
               drift_vel = drift_numer(r_L, vel[0], spdata);
               diff1 = diffusion.GetComponent(0, t, pos, mom, spdata);
               diff2 = diffusion.GetComponent(1, t, pos, mom, spdata);
            };
            Den_file << std::setw(18) << spdata.n_dens * unit_number_density_fluid;
            Pth_file << std::setw(18) << spdata.p_ther * unit_pressure_fluid;
            AbsVel_file << std::setw(18) << spdata.Uvec.Norm() * unit_velocity_fluid;
            RadVel_file << std::setw(18) << rad_vel;
            DivVel_file << std::setw(18) << spdata.divU() * unit_velocity_fluid / unit_length_fluid;
            AbsMag_file << std::setw(18) << spdata.Bmag * unit_magnetic_fluid;
            PolMag_file << std::setw(18) << polarity;
            dmax_file << std::setw(18) << spdata.dmax;
            HetFlx_file << std::setw(18) << spdata.region[0];
            TurEnr_file << std::setw(18) << correct_Z2(spdata, pos.Norm());
            drift_file << std::setw(18) << drift_vel.Norm() / vel[0];
            diff1_file << std::setw(18) << diff1 * unit_diffusion_fluid;
            diff2_file << std::setw(18) << diff2 * unit_diffusion_fluid;
         };
         Den_file << std::endl;
         Pth_file << std::endl;
         AbsVel_file << std::endl;
         RadVel_file << std::endl;
         DivVel_file << std::endl;
         AbsMag_file << std::endl;
         PolMag_file << std::endl;
         dmax_file << std::endl;
         HetFlx_file << std::endl;
         TurEnr_file << std::endl;
         drift_file << std::endl;
         diff1_file << std::endl;
         diff2_file << std::endl;
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

      background.StopServerFront();
      std::cout << "Background samples outputted to file." << std::endl;
   };

   std::cout << "CPU " << mpi_config->glob_comm_rank << " exited." << std::endl;
   return 0;
};


