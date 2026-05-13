# Code to compare the slopes of CIR modulation results with observations

# Import libraries
import matplotlib.pyplot as plt
import numpy as np
import argparse
from scipy.optimize import curve_fit

# Parse arguments
parser = argparse.ArgumentParser(description="Compare the slopes of CIR modulation results to observations")
parser.add_argument("input_file", type=str, help="File with dates of CIRs in YYYYMMDD format and time in days relative to CIR observed date to set as the zero epoch in plots.")
parser.add_argument("--phase", type=int, choices=[0, 1, 2], default=1, help="Index for which phase to plot the collective fits. 0: onset, 1: drop, 2: recovery. Default is 1.")
args = parser.parse_args()

# Compute percentage change of array
def percentage_change(array):
   avg = np.mean(array)
   return (array - avg) / avg * 100.0

# Fit function
def line(x, a, b):
   return a * x + b

# Plot configuration for fits
sim_fit_plot_config = {
   'color': 'tab:blue',
   'linestyle': '--',
   'linewidth': 3
}

ACE_fit_plot_config = {
   'color': 'k',
   'linestyle': '--',
   'linewidth': 3
}

sim_fit_plot_config2 = {
   'capsize': 4,
   'linewidth': 3,
   'color':'tab:blue',
   'marker': 'o',
   'markersize': 10
}

ACE_fit_plot_config2 = {
   'capsize': 4,
   'linewidth': 3,
   'color':'k',
   'marker': 's',
   'markersize': 10
}

# Read dates/epochs and allocate data
dates_epochs = np.loadtxt(args.input_file)
n_dates = np.size(dates_epochs, 0)
fit_params = np.zeros([n_dates, 6])
err_params = np.zeros([n_dates, 6])

# Iterate over dates
pct = 0.5
pct_pass = []
for i in range(n_dates):
# Import data
   date = np.rint(dates_epochs[i,0]).astype(int)
   zero_epoch = dates_epochs[i,1]
   print("Processing date {:d}".format(date))

# ACE data
   n_time = 193
   onset_ACE = slice(0, n_time // 4 + 1)
   drop_ACE = slice(3 * n_time // 8, 5 * n_time // 8 + 1)
   recov_ACE = slice(n_time - n_time // 4 - 1, -1)
   ACE_CR = np.loadtxt("data/SI_{:d}_CR.txt".format(date))
   ACE_CR[:,1] = percentage_change(ACE_CR[:,1])

# Simulation data
   n_time = 65
   onset_sim = slice(0, n_time // 4 + 1)
   drop_sim = slice(3 * n_time // 8, 5 * n_time // 8 + 1)
   recov_sim = slice(n_time - n_time // 4 - 1, -1)
   gcr_dat = np.loadtxt("output_{:d}/GCR_full/cir_gcr_mod_sim_rate_{:d}.dat".format(date, date))
   gcr_tim = gcr_dat[:,0] - zero_epoch
   gcr_rat = percentage_change(gcr_dat[:,1])

# Determine whether GCR simulation is "reasonable"
   gcr_max = np.max(gcr_rat[0:n_time // 2 + 1])
   gcr_min = -np.min(gcr_rat[n_time // 2: -1])
   if (gcr_max > pct) and (gcr_min > pct):
      print("\tGCR response passed {:.1f}% thresholds.".format(pct))
      pct_pass.append(i)
   else:
      print("\tGCR response did NOT pass {:.1f}% thresholds.".format(pct))

# Fit data and simulation in onset, drop, and recovery phases
   popt_onset_sim, pcov_onset_sim = curve_fit(line, gcr_tim[onset_sim], gcr_rat[onset_sim])
   popt_drop_sim, pcov_drop_sim = curve_fit(line, gcr_tim[drop_sim], gcr_rat[drop_sim])
   popt_recov_sim, pcov_recov_sim = curve_fit(line, gcr_tim[recov_sim], gcr_rat[recov_sim])
   popt_onset_ACE, pcov_onset_ACE = curve_fit(line, ACE_CR[onset_ACE,0], ACE_CR[onset_ACE,1])
   popt_drop_ACE, pcov_drop_ACE = curve_fit(line, ACE_CR[drop_ACE,0], ACE_CR[drop_ACE,1])
   popt_recov_ACE, pcov_recov_ACE = curve_fit(line, ACE_CR[recov_ACE,0], ACE_CR[recov_ACE,1])

# Save fit parameters
   fit_params[i,0] = popt_onset_sim[0]
   fit_params[i,1] = popt_drop_sim[0]
   fit_params[i,2] = popt_recov_sim[0]
   fit_params[i,3] = popt_onset_ACE[0]
   fit_params[i,4] = popt_drop_ACE[0]
   fit_params[i,5] = popt_recov_ACE[0]
   err_params[i,0] = np.sqrt(pcov_onset_sim[0, 0])
   err_params[i,1] = np.sqrt(pcov_drop_sim[0, 0])
   err_params[i,2] = np.sqrt(pcov_recov_sim[0, 0])
   err_params[i,3] = np.sqrt(pcov_onset_ACE[0, 0])
   err_params[i,4] = np.sqrt(pcov_drop_ACE[0, 0])
   err_params[i,5] = np.sqrt(pcov_recov_ACE[0, 0])

# Make plots of rates and fits
   fig = plt.figure(figsize=(14, 8), layout='tight')
   ax = fig.add_subplot(111, projection='rectilinear')

   ax.scatter(gcr_tim, gcr_rat, s=50, c="tab:blue",
              marker="o", label="full")
   ax.plot(gcr_tim[onset_sim], line(gcr_tim[onset_sim], *popt_onset_sim),
           **sim_fit_plot_config)
   ax.plot(gcr_tim[drop_sim], line(gcr_tim[drop_sim], *popt_drop_sim),
           **sim_fit_plot_config, label="full fit")
   ax.plot(gcr_tim[recov_sim], line(gcr_tim[recov_sim], *popt_recov_sim),
           **sim_fit_plot_config)
   ax.plot(ACE_CR[:,0], ACE_CR[:,1],
           linewidth=3, color='k',
           linestyle='-', label="ACE")
   ax.plot(ACE_CR[onset_ACE,0], line(ACE_CR[onset_ACE,0], *popt_onset_ACE),
           **ACE_fit_plot_config)
   ax.plot(ACE_CR[drop_ACE,0], line(ACE_CR[drop_ACE,0], *popt_drop_ACE),
           **ACE_fit_plot_config, label="ACE fit")
   ax.plot(ACE_CR[recov_ACE,0], line(ACE_CR[recov_ACE,0], *popt_recov_ACE),
           **ACE_fit_plot_config)
   ax.set_xlabel("time (days)", fontsize = 20)
   ax.set_ylabel("rate change (%)", fontsize = 20)
   ax.set_xlim(gcr_tim[0],gcr_tim[-1])
   ax.tick_params(labelsize=20)
   ax.legend(fontsize=20)

   plt.savefig("output_{:d}/gcr_rate_comparison_{:d}.png".format(date, date))
   plt.close(fig)

# Make plots of fit parameters
idx = args.phase
if idx == 0:
   phase = "onset"
elif idx == 1:
   phase = "drop"
elif idx == 2:
   phase = "recovery"

fig = plt.figure(figsize=(14, 8), layout='tight')
ax = fig.add_subplot(111, projection='rectilinear')

k = 0
ax.errorbar([k], [fit_params[pct_pass[0], idx]], yerr=err_params[pct_pass[0], idx],
            **sim_fit_plot_config2, label="sim")
ax.errorbar([k], [fit_params[pct_pass[0], idx+3]], yerr=err_params[pct_pass[0], idx+3],
            **ACE_fit_plot_config2, label="ACE")
for i in range(1, len(pct_pass)):
   j = pct_pass[i]
   k = k + 1
   ax.errorbar([k], [fit_params[j, idx]], yerr=err_params[j, idx],
               **sim_fit_plot_config2)
   ax.errorbar([k], [fit_params[j, idx+3]], yerr=err_params[j, idx+3],
               **ACE_fit_plot_config2)

ax.set_xticks(range(len(pct_pass)))
ax.set_xticklabels(["{:d}".format(pct+1) for pct in pct_pass])
ax.set_xlabel("CIR index", fontsize = 20)
ax.set_ylabel("{:s} phase slope".format(phase), fontsize = 20)
ax.tick_params(labelsize=20)
ax.legend(fontsize=20)

plt.savefig("gcr_rate_fits_comparison_{:s}.png".format(phase))
plt.close(fig)

error = 0
print(" CIR   SIM s  SIM e  ACE s  ACE e")
print("-----  -----  -----  -----  -----")
for i in range(len(pct_pass)):
   j = pct_pass[i]
   print("{:5d}  {:5.2f}  {:5.2f}  {:5.2f}  {:5.2f}".format(j+1, fit_params[j, idx], err_params[j, idx], fit_params[j, idx+3], err_params[j, idx+3]))
   error = error + np.abs(fit_params[j, idx] - fit_params[j, idx+3])
error = error / len(pct_pass)
print("Average error: {:.2f}".format(error))
