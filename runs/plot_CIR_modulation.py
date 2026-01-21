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
gcr_dat1 = np.loadtxt("output_{:s}/GCR_full/cir_gcr_mod_sim_rate_{:s}.dat".format(cir_date, cir_date))
gcr_tim1 = gcr_dat1[:,0]
gcr_rat1 = percentage_change(gcr_dat1[:,1])
gcr_res1 = gcr_dat1[:,2]
gcr_dat2 = np.loadtxt("output_{:s}/GCR_ndft/cir_gcr_mod_sim_rate_{:s}.dat".format(cir_date, cir_date))
gcr_tim2 = gcr_dat2[:,0]
gcr_rat2 = percentage_change(gcr_dat2[:,1])
gcr_res2 = gcr_dat2[:,2]
vel_data = np.loadtxt("output_{:s}/CIR/vel_1au_{:s}.dat".format(cir_date, cir_date))
vel_time = vel_data[:,0]
vel_magn = vel_data[:,1]
vel_magn /= 1.0e5 # cm -> km
pol_data = np.loadtxt("output_{:s}/CIR/pol_1au_{:s}.dat".format(cir_date, cir_date))
pol_time = pol_data[:,0]
pol_sign = pol_data[:,1]

# Average velocity
n_vel = np.size(vel_magn)
n_days = round(n_vel * 1.0 / 27.0)
forw_avg_vel = np.zeros(n_vel)
for i in range(n_vel):
   for j in range(i, i+n_days):
      forw_avg_vel[i] += vel_magn[j % n_vel]
back_avg_vel = np.zeros(n_vel)
for i in range(n_vel):
   for j in range(i-n_days, i):
      back_avg_vel[i] += vel_magn[j % n_vel]

# Compute gradient of average velocity
grad_vel = np.zeros(n_vel)
for i in range(n_vel):
   grad_vel[i] = forw_avg_vel[i] - back_avg_vel[i]

# Find max of gradient of velocity
vel_idx = np.argmax(grad_vel)
# Find index in GCR data that is closest to the gradient peak
n_gcr = np.size(gcr_rat1)
gcr_idx = round(n_gcr * vel_idx / n_vel)
# Adjust MHD indices (higher resolution) to match the GCRs
vel_idx = round(n_vel * gcr_idx / n_gcr)
n_pol = np.size(pol_sign)
pol_idx = round(n_pol * gcr_idx / n_gcr)

# Shift arrays
vel_magn = periodic_shift(vel_magn, vel_idx)
pol_sign = periodic_shift(pol_sign, pol_idx)
gcr_rat1 = periodic_shift(gcr_rat1, gcr_idx)
gcr_res1 = periodic_shift(gcr_res1, gcr_idx)
gcr_rat2 = periodic_shift(gcr_rat2, gcr_idx)
gcr_res2 = periodic_shift(gcr_res2, gcr_idx)

# Find polarity changes
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
ax.plot(gcr_tim1, gcr_rat1,
        linewidth=3, color="tab:blue",
        linestyle="-", label="full")
ax.plot(gcr_tim2, gcr_rat2,
        linewidth=3, color="tab:blue",
        linestyle="--", label="no drift")
ax.set_xlabel("time (days)", fontsize = 20)
ax.set_ylabel("rate change (%)", fontsize = 20,
               color="tab:blue")
ax.set_xlim(0, 27.0)
ax.tick_params(labelsize=20)
ax.tick_params(axis='y', labelcolor="tab:blue")
ax.legend(fontsize=20)

# Velocity magnitude
axx = ax.twinx()
axx.plot(vel_time, vel_magn,
         linewidth=3, color="tab:orange")
axx.set_ylabel("Flow speed (km/s)", fontsize=20,
               color="tab:orange")
axx.tick_params(labelsize=20)
axx.tick_params(axis='y')
axx.tick_params(axis='y', labelcolor="tab:orange")

# Current sheet crossings
for idx in hcs_pls:
   axx.axvline(pol_time[idx], color="k",
               linestyle=":", linewidth=2)
for idx in hcs_mns:
   axx.axvline(pol_time[idx], color="r",
               linestyle=":", linewidth=2)

plt.savefig("output_{:s}/gcr_rate_{:s}.png".format(cir_date, cir_date))
plt.close(fig)

fig = plt.figure(figsize=(14, 8), layout='tight')
ax = fig.add_subplot(111, projection='rectilinear')

# GCR residence time
ax.plot(gcr_tim1, gcr_res1,
        linewidth=3, color="tab:blue",
        linestyle="-", label="full")
ax.plot(gcr_tim2, gcr_res2,
        linewidth=3, color="tab:blue",
        linestyle="--", label="no drift")
ax.set_xlabel("time (days)", fontsize = 20)
ax.set_ylabel("residence time (days)", fontsize = 20,
               color="tab:blue")
ax.set_xlim(0, 27.0)
ax.tick_params(labelsize=20)
ax.tick_params(axis='y', labelcolor="tab:blue")
ax.legend(fontsize=20)

# Velocity magnitude
axx = ax.twinx()
axx.plot(vel_time, vel_magn,
         linewidth=3, color="tab:orange")
axx.set_ylabel("Flow speed (km/s)", fontsize=20,
               color="tab:orange")
axx.tick_params(labelsize=20)
axx.tick_params(axis='y')
axx.tick_params(axis='y', labelcolor="tab:orange")

# Current sheet crossings
for idx in hcs_pls:
   axx.axvline(pol_time[idx], color="k",
               linestyle=":", linewidth=2)
for idx in hcs_mns:
   axx.axvline(pol_time[idx], color="r",
               linestyle=":", linewidth=2)

plt.savefig("output_{:s}/gcr_res_time_{:s}.png".format(cir_date, cir_date))
plt.close(fig)