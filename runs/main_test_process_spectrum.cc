#include "src/distribution_other.hh"
#include "common/physics.hh"
#include <fstream>
#include <iostream>
#include <iomanip>
#include <cmath>

using namespace Spectrum;

const double J0 = 20.4;
const double mu0 = 0.986;
const double T0 = SPC_CONST_CGSM_GIGA_ELECTRON_VOLT;
const double Tb = 1.01 * T0;
const double mu1 = -2.81;
const double ds = 3.98;
const int specie = SPECIES_PROTON_BEAM;

// Unmodulated spectrum
inline double unmod_spectrum(double T) {return J0 * pow(T / T0, mu0) / pow(1.0 + pow(T / Tb, (mu0-mu1)/ds), ds);};

int main(int argc, char** argv)
{
   std::string cir_date;
   if(argc > 1) cir_date = argv[1];
   std::string init_time;
   if(argc > 2) init_time = argv[2];
   std::string distroname = "output_" + cir_date + "/cir_gcr_mod_" + init_time + "_distro_0.out";
   std::string infilename = "output_" + cir_date + "/cir_gcr_mod_" + init_time + "_spec.dat";
   std::string outfilename = "output_" + cir_date + "/cir_gcr_mod_" + init_time + "_spec_pp.dat";
   std::string line;
   int i, N = 100;
   int sum_c[N];
   double energy[N], distro[N], sum_w[N];
   double p2, energy0 = SPC_CONST_CGSM_GIGA_ELECTRON_VOLT;

// Restore distribution
   DistributionSpectrumKineticEnergyBentPowerLaw distro_obj;
   distro_obj.Restore(distroname);
   distro_obj.Print1D(0, infilename, true);

// Open input analytic distro file
   std::ifstream input_spectrum_file(infilename);

// Read first two lines of distro file
   std::getline(input_spectrum_file, line);
   std::getline(input_spectrum_file, line);

// Read data
   for(i = 0; i < N; i++) {
      input_spectrum_file >> energy[i];
      input_spectrum_file >> distro[i];
      input_spectrum_file >> sum_w[i];
      input_spectrum_file >> sum_c[i];
   };

// Close input cartesian distro file
   input_spectrum_file.close();

// Open output distro file
   std::ofstream output_spectrum_file(outfilename);

// Output data
   output_spectrum_file << std::setprecision(8);
   for(i = 0; i < N; i++) {
      p2 = Sqr(Mom(energy[i] / unit_energy_particle, specie));
      output_spectrum_file << std::setw(20) << energy[i] / energy0;
      output_spectrum_file << std::setw(20) << unmod_spectrum(energy[i]);
      output_spectrum_file << std::setw(20) << p2 * distro[i];
      output_spectrum_file << std::endl;
   };

// Close output distro file
   output_spectrum_file.close();

   return 0;
};
