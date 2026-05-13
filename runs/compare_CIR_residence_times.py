# Code to plot the change in residence times between ndif and full simulations vs the radial diffusion coefficient along Earth's trajectory

# Import libraries
import matplotlib.pyplot as plt
import numpy as np
import argparse
from scipy import stats
from scipy.optimize import curve_fit

# Parse arguments
parser = argparse.ArgumentParser(description="Plot the change in residence times between ndif and full simulations vs the radial diffusion coefficient along Earth's trajectory")
parser.add_argument("input_file", type=str, help="File with dates of CIRs in YYYYMMDD format and time in days relative to CIR observed date to set as the zero epoch in plots.")
args = parser.parse_args()

# Fit function
def line(x, a, b):
   return a * x + b

# Read dates/epochs and allocate data
dates_epochs = np.loadtxt(args.input_file)
n_dates = np.size(dates_epochs, 0)
avg_res_tim_dif = np.zeros(n_dates)
avg_kap = np.zeros(n_dates)

# Iterate over dates
for i in range(n_dates):
# Import data
   date = np.rint(dates_epochs[i,0]).astype(int)
   print("Processing date {:d}".format(date))

# Residence time
   gcr_full_data = np.loadtxt("output_{:d}/GCR_full/cir_gcr_mod_sim_rate_{:d}.dat".format(date, date))
   gcr_ndif_dat = np.loadtxt("output_{:d}/GCR_ndif/cir_gcr_mod_sim_rate_{:d}.dat".format(date, date))
   avg_res_tim_dif[i] = np.mean(gcr_ndif_dat[:,2] - gcr_full_data[:,2])

# Average diffusion
   diff_data = np.loadtxt("output_{:d}/CIR/diff3_1au_{:d}.dat".format(date, date))
   avg_kap[i] = np.mean(diff_data[:,1]) / 1.0e22

# Find linear fit
popt, pcov = curve_fit(line, avg_kap, avg_res_tim_dif)
kap = np.linspace(np.min(avg_kap), np.max(avg_kap))
res_tim_dif = line(kap, *popt)

# Make plot
SC_23_24 = slice(0,8,1)
SC_24_25 = slice(8,14,1)
fig = plt.figure(figsize=(14, 8), layout='tight')
ax = fig.add_subplot(111, projection='rectilinear')

ax.scatter(avg_kap[SC_23_24], avg_res_tim_dif[SC_23_24],
           s=100, c='tab:blue', marker='o', label="SC 23/24")
ax.scatter(avg_kap[SC_24_25], avg_res_tim_dif[SC_24_25],
           s=100, c='tab:orange', marker='s', label="SC 24/25")
ax.plot(kap, res_tim_dif, 'k--', linewidth=3, label='linear fit')

ax.set_xlabel("$\\kappa_{rr}$ ($10^{22}$ cm$^2$ s$^{-1}$)", fontsize = 20)
ax.set_ylabel("$\\Delta$ residence time (days)", fontsize = 20)
ax.tick_params(labelsize=20)
ax.legend(fontsize=20)

plt.savefig("gcr_res_time_diff_vs_perp_diff.png")
plt.close(fig)

corr = stats.pearsonr(avg_kap, avg_res_tim_dif)
print("Correlation coefficient: {:.2f}".format(corr.statistic))