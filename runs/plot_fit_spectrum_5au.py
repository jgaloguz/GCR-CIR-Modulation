# Code to plot and fit LISM spectra

# Import modules
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.interpolate import interp1d
from scipy.integrate import quad
import sys

# Broken power law
def log_broken_pow_law(T, J, T_b, a1, a2, d):
   return np.log(J * (T / 1000)**a1 / (1.0 + (T / T_b)**(a2/d))**d)

def broken_pow_law(T, J, T_b, a1, a2, d):
   return np.exp(log_broken_pow_law(T, J, T_b, a1, a2, d))

# Check specie
if len(sys.argv) < 2:
   print("Error: data set must be specified.")
   print("Please specify a year.")
   exit(1)
   
print("Plotting results for year {:s}.".format(sys.argv[1]))

# Import data
spectrum = np.loadtxt("data/spectrum-{:s}.txt".format(sys.argv[1]), skiprows=3)
energy = spectrum[:,0]
flux = spectrum[:,1]
guess_params = [50.0, 60.0, 0.5, 3.0, 5.0]
low_params = [10.0, 1.0, 0.0, 2.0, 1.0]
high_params = [100.0, 2000.0, 2.0, 4.0, 10.0]

# Fit data
opt_params, cov = curve_fit(log_broken_pow_law, energy, np.log(flux),
                            p0=guess_params, maxfev=10000,
                            bounds=(low_params, high_params))

# Plot
fig = plt.figure(figsize=(8, 10), layout='tight')
ax = fig.add_subplot(111, projection='rectilinear')
ax.loglog(energy, flux, linestyle="", marker="o", color="tab:orange", label="FFA", markersize=6)
ax.loglog(energy, broken_pow_law(energy, *opt_params), color="tab:blue", label="BPL Fit", linewidth=3)
ax.set_xlabel("Kinetic Energy (MeV)", fontsize = 20)
ax.set_ylabel("H$^+$ Flux (m$^2$ s sr MeV)$^{-1}$", fontsize = 20)
ax.tick_params(labelsize=20)
ax.legend(fontsize = 20)

plt.show()
plt.savefig("spectrum-{:s}-fit.png".format(sys.argv[1]))
plt.close(fig)
print("Optimal Fit Parameters:")
print(opt_params)