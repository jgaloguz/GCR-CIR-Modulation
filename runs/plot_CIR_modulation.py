# Code to plot CIR modulation results

# Import libraries
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm
import scipy.signal as spsig
import sys

# Check CIR date
if len(sys.argv) > 1:
   cir_date = sys.argv[1]
   print("CIR date: {:s}".format(cir_date))
else:
   print("ERROR: No CIR date provided.")
   sys.exit(1)

# Function to shift center 1D array to a particular index
def periodic_shift(array, index):
   N = np.size(array)
   shifted_array = np.zeros(N)
   N_2 = N//2
   if index < N_2:
      chop = N-index-N_2
   elif index > N_2:
      chop = N-index+N_2
   else:
      chop = 0
   shifted_array[:chop] = array[-chop:]
   shifted_array[chop:] = array[:-chop]
   return shifted_array

# Compute percentage change of array
def percentage_change(array):
   avg = np.mean(array)
   return (array - avg) / avg * 100.0

# Import data
ACE_CR = np.loadtxt("data/SI_{:s}_CR.txt".format(cir_date))
ACE_CR[:,1] = percentage_change(ACE_CR[:,1])

GCR_labels = ["full", "ndft", "nadv", "ndif"]
n_time = 65
t0 = 11.8
n_GCR = len(GCR_labels)
gcr_dat = [np.zeros(n_time) for _ in range(n_GCR)]
gcr_tim = [np.zeros(n_time) for _ in range(n_GCR)]
gcr_rat = [np.zeros(n_time) for _ in range(n_GCR)]
gcr_res = [np.zeros(n_time) for _ in range(n_GCR)]
for i in range(n_GCR):
   gcr_dat[i] = np.loadtxt("output_{:s}/GCR_{:s}/cir_gcr_mod_sim_rate_{:s}.dat".format(cir_date, GCR_labels[i], cir_date))
   gcr_tim[i] = gcr_dat[i][:,0] - t0
   gcr_rat[i] = percentage_change(gcr_dat[i][:,1])
   gcr_res[i] = gcr_dat[i][:,2]
vel_data = np.loadtxt("output_{:s}/CIR/vel_1au_{:s}.dat".format(cir_date, cir_date))
vel_time = vel_data[:,0] - t0
vel_magn = vel_data[:,1]
vel_magn /= 1.0e5 # cm -> km
pol_data = np.loadtxt("output_{:s}/CIR/pol_1au_{:s}.dat".format(cir_date, cir_date))
pol_time = pol_data[:,0] - t0
pol_sign = pol_data[:,1]

# Find polarity changes
n_pol = np.size(pol_sign)
pol_chgs = np.zeros(n_pol)
for i in range(n_pol-1):
   pol_chgs[i] = pol_sign[i+1] - pol_sign[i]
pol_chgs[n_pol-1] = pol_sign[0] - pol_sign[n_pol-1]
hcs_pls = np.where(pol_chgs > 1.0)[0]
hcs_mns = np.where(pol_chgs < -1.0)[0]

# Make plots
fig = plt.figure(figsize=(14, 8), layout='tight')
ax = fig.add_subplot(111, projection='rectilinear')

# GCR rate
marker_styles=["o","^","s","X"]
marker_colors=["tab:blue", "deepskyblue", "aqua", "navy"]
for i in range(n_GCR):
   ax.scatter(gcr_tim[i], gcr_rat[i], s=50, c=marker_colors[i],
              marker=marker_styles[i], label=GCR_labels[i])
ax.plot(ACE_CR[:,0], ACE_CR[:,1],
        linewidth=3, color='k',
        linestyle='-', label="ACE")
ax.set_xlabel("time (days)", fontsize = 20)
ax.set_ylabel("rate change (%)", fontsize = 20)
ax.set_xlim(gcr_tim[i][0],gcr_tim[i][-1])
ax.tick_params(labelsize=20)
ax.legend(fontsize=20)

# Velocity magnitude
axx = ax.twinx()
axx.plot(vel_time, vel_magn,
         linewidth=3, color="tab:orange")
axx.set_ylabel("sim. flow speed (km/s)", fontsize=20,
               color="tab:orange")
axx.tick_params(labelsize=20)
axx.tick_params(axis='y')
axx.tick_params(axis='y', labelcolor="tab:orange")

# Current sheet crossings
for idx in hcs_pls:
   axx.axvline(pol_time[idx], color="g",
               linestyle="--", linewidth=2)
for idx in hcs_mns:
   axx.axvline(pol_time[idx], color="r",
               linestyle="--", linewidth=2)

plt.savefig("output_{:s}/gcr_rate_{:s}.png".format(cir_date, cir_date))
plt.close(fig)

fig = plt.figure(figsize=(14, 8), layout='tight')
ax = fig.add_subplot(111, projection='rectilinear')

# GCR residence time
for i in range(n_GCR):
   ax.scatter(gcr_tim[i], gcr_res[i], s=50, c=marker_colors[i],
              marker=marker_styles[i], label=GCR_labels[i])
ax.set_xlabel("time (days)", fontsize = 20)
ax.set_ylabel("residence time (days)", fontsize = 20)
ax.set_xlim(gcr_tim[i][0],gcr_tim[i][-1])
ax.tick_params(labelsize=20)
ax.legend(fontsize=20)

# Velocity magnitude
axx = ax.twinx()
axx.plot(vel_time, vel_magn,
         linewidth=3, color="tab:orange")
axx.set_ylabel("sim. flow speed (km/s)", fontsize=20,
               color="tab:orange")
axx.tick_params(labelsize=20)
axx.tick_params(axis='y')
axx.tick_params(axis='y', labelcolor="tab:orange")

# Current sheet crossings
for idx in hcs_pls:
   axx.axvline(pol_time[idx], color="g",
               linestyle="--", linewidth=2)
for idx in hcs_mns:
   axx.axvline(pol_time[idx], color="r",
               linestyle="--", linewidth=2)

plt.savefig("output_{:s}/gcr_res_time_{:s}.png".format(cir_date, cir_date))
plt.close(fig)