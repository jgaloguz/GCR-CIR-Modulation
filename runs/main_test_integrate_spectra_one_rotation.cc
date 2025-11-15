#include "common/physics.hh"
#include <fstream>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <sstream>

using namespace Spectrum;

const double J0 = 20.4;
const double mu0 = 0.986;
const double T0 = 1000.0;
const double Tb = 1.01 * T0;
const double mu1 = -2.81;
const double ds = 3.98;
const double response = 1.68e-4;

// Unmodulated spectrum
inline double unmod_spectrum(double T) {return J0 * pow(T / T0, mu0) / pow(1.0 + pow(T / Tb, (mu0-mu1)/ds), ds);};

// Integrate differential intensity (trapezoid rule)
double integrate_diff_int(double *J, double *T, int N)
{
   double S = 0.0;
   for(int i = 0; i < N - 1; i++) S += 0.5 * (J[i] + J[i+1]) * (T[i+1] - T[i]);
   return response * S;
};

int main(int argc, char** argv)
{
   int eng, n_eng = 100;
   int day, n_day = 100;
   int t_idx, n_time = 1;
   double time_dbl;
   std::string time_str;
   std::string cir_date;
   if(argc > 1) cir_date = argv[1];
   if(argc > 2) n_time = atoi(argv[2]);
   std::ifstream input_spectrum_file;
   std::ofstream output_spectrum_file;
   std::string infilename;
   std::string outfilename = "output_" + cir_date + "/cir_gcr_mod_sim_rate.dat";
   double energy[n_eng], diff_int[n_eng], sum_w1[n_eng];
   double days, frac, avg_days;
   double skip;

// Open output distro file
   output_spectrum_file.open(outfilename);

// Iterate over segments endpoints
   for(t_idx = 0; t_idx < n_time; t_idx++) {

// Find initial time
      time_dbl = t_idx * 27.0 / n_time;
      std::stringstream ss;
      ss << std::fixed << std::setprecision(1) << time_dbl;
      time_str = ss.str();

// Status message
      std::cerr << "\tInitial time " << time_str << std::endl;

// Open input differential intensity file
      infilename = "output_" + cir_date + "/cir_gcr_mod_" + time_str + "_spec_pp.dat";
      input_spectrum_file.open(infilename);

// Read differential intensity
      for(eng = 0; eng < n_eng; eng++) {
         input_spectrum_file >> energy[eng];
         energy[eng] *= T0;
         input_spectrum_file >> skip;
         input_spectrum_file >> diff_int[eng];
      };

// Close input differential intensity file
      input_spectrum_file.close();

// Open input time distribution file
      infilename = "output_" + cir_date + "/cir_gcr_mod_" + time_str + "_time_pp.dat";
      input_spectrum_file.open(infilename);

// Read fraction of particles for each residence time
      avg_days = 0.0;
      for(day = 0; day < n_day; day++) {
         input_spectrum_file >> days;
         input_spectrum_file >> frac;
         avg_days += -(days - time_dbl) * frac;
      };

// Output integrated differential intensity and average residence time in days
      output_spectrum_file << std::setw(20) << time_dbl
                           << std::setw(20) << integrate_diff_int(diff_int, energy, n_eng)
                           << std::setw(20) << avg_days
                           << std::endl;

// Close input differential intensity file
      input_spectrum_file.close();
   };

// Close output distro file
      output_spectrum_file.close();

// Output 5 au integrated differential intensity for reference
   for(eng = 0; eng < n_eng; eng++) diff_int[eng] = unmod_spectrum(energy[eng]);
   std::cout << std::setw(20) << integrate_diff_int(diff_int, energy, n_eng) << std::endl;

   return 0;
};
