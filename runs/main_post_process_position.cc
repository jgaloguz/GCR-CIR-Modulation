#include "src/distribution_other.hh"
#include "common/physics.hh"
#include <fstream>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <filesystem>

using namespace Spectrum;

int main(int argc, char** argv)
{
   std::string cir_date;
   if(argc > 1) cir_date = argv[1];
   std::string init_time;
   if(argc > 2) init_time = argv[2];
   std::ifstream input_spectrum_file;
   std::ofstream output_spectrum_file;
   std::string distroname = "output_" + cir_date + "/GCR/cir_gcr_mod_" + init_time + "_distro_2.out";
   std::string infilename1 = "output_" + cir_date + "/GCR/cir_gcr_mod_" + init_time + "_lat.dat";
   std::string infilename2 = "output_" + cir_date + "/GCR/cir_gcr_mod_" + init_time + "_lon.dat";
   std::string outfilename1 = "output_" + cir_date + "/GCR/cir_gcr_mod_" + init_time + "_lat_pp.dat";
      std::string outfilename2 = "output_" + cir_date + "/GCR/cir_gcr_mod_" + init_time + "_lon_pp.dat";
   std::string line;
   int i, N_lat = 90, N_lon = 180;
   int sum_c_lat[N_lat], total_c_lat = 0;
   int sum_c_lon[N_lon], total_c_lon = 0;
   double lat[N_lat], distro_lat[N_lat], sum_w_lat[N_lat];
   double lon[N_lon], distro_lon[N_lon], sum_w_lon[N_lon];

// Restore distribution
   DistributionPositionUniform distro_obj;
   distro_obj.Restore(distroname);
   distro_obj.Print1D(1, infilename1, true);
   distro_obj.Print1D(2, infilename2, true);

// Open input distro file -- LATITUDE
   input_spectrum_file.open(infilename1);

// Read first two lines of distro file
   std::getline(input_spectrum_file, line);
   std::getline(input_spectrum_file, line);

// Read data
   for(i = 0; i < N_lat; i++) {
      input_spectrum_file >> lat[i];
      input_spectrum_file >> distro_lat[i];
      input_spectrum_file >> sum_w_lat[i];
      input_spectrum_file >> sum_c_lat[i];
      total_c_lat += sum_c_lat[i];
   };

// Close input distro file
   input_spectrum_file.close();

// Open output distro file
   output_spectrum_file.open(outfilename1);

// Output data
   output_spectrum_file << std::setprecision(8);
   for(i = 0; i < N_lat; i++) {
      output_spectrum_file << std::setw(20) << 90.0 - RadToDeg(lat[i]);
      output_spectrum_file << std::setw(20) << (double)sum_c_lat[i] / (double)total_c_lat;
      output_spectrum_file << std::endl;
   };

// Close output distro file
   output_spectrum_file.close();

// Open input distro file -- LONGITUDE
   input_spectrum_file.open(infilename2);

// Read first two lines of distro file
   std::getline(input_spectrum_file, line);
   std::getline(input_spectrum_file, line);

// Read data
   for(i = 0; i < N_lon; i++) {
      input_spectrum_file >> lon[i];
      input_spectrum_file >> distro_lon[i];
      input_spectrum_file >> sum_w_lon[i];
      input_spectrum_file >> sum_c_lon[i];
      total_c_lon += sum_c_lon[i];
   };

// Close input distro file
   input_spectrum_file.close();

// Open output distro file
   output_spectrum_file.open(outfilename2);

// Output data
   output_spectrum_file << std::setprecision(8);
   for(i = 0; i < N_lon; i++) {
      output_spectrum_file << std::setw(20) << RadToDeg(lon[i]);
      output_spectrum_file << std::setw(20) << (double)sum_c_lon[i] / (double)total_c_lon;
      output_spectrum_file << std::endl;
   };

// Close output distro file
   output_spectrum_file.close();

// Delete distro and input files
   std::filesystem::remove(distroname);
   std::filesystem::remove(infilename1);
   std::filesystem::remove(infilename2);

   return 0;
};
