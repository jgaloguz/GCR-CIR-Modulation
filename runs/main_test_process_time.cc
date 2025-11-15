#include "src/distribution_other.hh"
#include "common/physics.hh"
#include <fstream>
#include <iostream>
#include <iomanip>
#include <cmath>

using namespace Spectrum;

int main(int argc, char** argv)
{
   std::string cir_date;
   if(argc > 1) cir_date = argv[1];
   std::string init_time;
   if(argc > 2) init_time = argv[2];
   std::string distroname = "output_" + cir_date + "/cir_gcr_mod_" + init_time + "_distro_1.out";
   std::string infilename = "output_" + cir_date + "/cir_gcr_mod_" + init_time + "_time.dat";
   std::string outfilename = "output_" + cir_date + "/cir_gcr_mod_" + init_time + "_time_pp.dat";
   std::string line;
   int i, N = 100;
   int sum_c[N], total_c = 0;
   double time[N], distro[N], sum_w[N];

// Restore distribution
   DistributionTimeUniform distro_obj;
   distro_obj.Restore(distroname);
   distro_obj.Print1D(0, infilename, true);

// Open input analytic distro file
   std::ifstream input_spectrum_file(infilename);

// Read first two lines of distro file
   std::getline(input_spectrum_file, line);
   std::getline(input_spectrum_file, line);

// Read data
   for(i = 0; i < N; i++) {
      input_spectrum_file >> time[i];
      input_spectrum_file >> distro[i];
      input_spectrum_file >> sum_w[i];
      input_spectrum_file >> sum_c[i];
      total_c += sum_c[i];
   };

// Close input cartesian distro file
   input_spectrum_file.close();

// Open output distro file
   std::ofstream output_spectrum_file(outfilename);

// Output data
   output_spectrum_file << std::setprecision(8);
   for(i = 0; i < N; i++) {
      output_spectrum_file << std::setw(20) << DAYS(time[i]);
      output_spectrum_file << std::setw(20) << (double)sum_c[i] / (double)total_c;
      output_spectrum_file << std::endl;
   };

// Close output distro file
   output_spectrum_file.close();

   return 0;
};
