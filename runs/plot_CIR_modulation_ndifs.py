# Code to plot CIR modulation results

# Import libraries
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
import argparse

# Parse arguments
parser = argparse.ArgumentParser(description="Plot CIR modulation results")
parser.add_argument("date", type=str, help="Date of CIR in YYYYMMDD format")
parser.add_argument("zero_epoch", type=float, help="Time in days relative to CIR observed date to set as the zero epoch in plots. Possible values are any number in range [-10.0,10.0].")
args = parser.parse_args()

# Compute percentage change of array
def percentage_change(array):
   avg = np.mean(array)
   return (array - avg) / avg * 100.0

# Import data
ACE_CR = np.loadtxt("data/SI_{:s}_CR.txt".format(args.date))
ACE_CR[:,1] = percentage_change(ACE_CR[:,1])

GCR_labels = ["full", "ndif0", "ndif1"]
n_time = 65
n_GCR = len(GCR_labels)
gcr_dat = [np.zeros(n_time) for _ in range(n_GCR)]
gcr_tim = [np.zeros(n_time) for _ in range(n_GCR)]
gcr_rat = [np.zeros(n_time) for _ in range(n_GCR)]
for i in range(n_GCR):
   gcr_dat[i] = np.loadtxt("output_{:s}/GCR_{:s}/cir_gcr_mod_sim_rate_{:s}.dat".format(args.date, GCR_labels[i], args.date))
   gcr_tim[i] = gcr_dat[i][:,0] - args.zero_epoch
   gcr_rat[i] = percentage_change(gcr_dat[i][:,1])
vel_data = np.loadtxt("output_{:s}/CIR/vel_1au_{:s}.dat".format(args.date, args.date))
vel_time = vel_data[:,0] - args.zero_epoch
vel_magn = vel_data[:,1]
vel_magn /= 1.0e5 # cm -> km
pol_data = np.loadtxt("output_{:s}/CIR/pol_1au_{:s}.dat".format(args.date, args.date))
pol_time = pol_data[:,0] - args.zero_epoch
pol_sign = pol_data[:,1]

# Find polarity changes
n_pol = np.size(pol_sign)
pol_chgs = np.zeros(n_pol)
for i in range(n_pol-1):
   pol_chgs[i] = pol_sign[i+1] - pol_sign[i]
pol_chgs[n_pol-1] = pol_sign[0] - pol_sign[n_pol-1]
hcs_pls = np.where(pol_chgs > 1.0)[0]
hcs_mns = np.where(pol_chgs < -1.0)[0]

# Plot GCR rate
fig = plt.figure(figsize=(14, 8), layout='tight')
ax = fig.add_subplot(111, projection='rectilinear')

marker_styles=["o","^","s","X"]
marker_colors=["tab:blue", "deepskyblue", "aqua", "navy"]
ax.scatter(gcr_tim[0], gcr_rat[0], s=50, c=marker_colors[0],
              marker=marker_styles[0], label="full")
ax.scatter(gcr_tim[1], gcr_rat[1], s=50, c=marker_colors[1],
              marker=marker_styles[1], label="ndif$\\perp$")
ax.scatter(gcr_tim[2], gcr_rat[2], s=50, c=marker_colors[2],
              marker=marker_styles[2], label="ndif$\\parallel$")
ax.plot(ACE_CR[:,0], ACE_CR[:,1],
        linewidth=3, color='k',
        linestyle='-', label="ACE")
ax.set_xlabel("time (days)", fontsize = 20)
ax.set_ylabel("rate change (%)", fontsize = 20)
ax.set_xlim(gcr_tim[0][0],gcr_tim[0][-1])
ax.tick_params(labelsize=20)
leg = ax.legend(fontsize=20)
leg.remove()

# Velocity magnitude
axx = ax.twinx()
axx.plot(vel_time, vel_magn,
         linewidth=3, color="tab:orange")
axx.set_ylabel("sim. flow speed (km/s)", fontsize=20,
               color="tab:orange")
axx.tick_params(labelsize=20)
axx.tick_params(axis='y')
axx.tick_params(axis='y', labelcolor="tab:orange")
axx.add_artist(leg)

# Current sheet crossings
for idx in hcs_pls:
   axx.axvline(pol_time[idx], color="g",
               linestyle="--", linewidth=2)
for idx in hcs_mns:
   axx.axvline(pol_time[idx], color="r",
               linestyle="--", linewidth=2)

ax.set_title("CIR {:s}".format(args.date), fontsize=24)
plt.savefig("output_{:s}/gcr_rate_ndifs_{:s}.png".format(args.date, args.date))
plt.close(fig)

# Get correlations
pcc_ndif0 = stats.pearsonr(gcr_rat[0], gcr_rat[1])
pcc_ndif1 = stats.pearsonr(gcr_rat[0], gcr_rat[2])
print("Pearson correlation coefficients relative to full simulation:")
print("   ndif0      ndif1")
print("  ------     ------")
print("   {:4.2f}      {:4.2f}".format(pcc_ndif0.statistic, pcc_ndif1.statistic))

# Get mean absolute error
mae_ndif0 = np.mean(np.abs(gcr_rat[0] - gcr_rat[1]))
mae_ndif1 = np.mean(np.abs(gcr_rat[0] - gcr_rat[2]))
print("Mean absolute error relative to full simulation:")
print("   ndif0      ndif1")
print("  ------     ------")
print("   {:4.2f}      {:4.2f}".format(mae_ndif0, mae_ndif1))