# Code to plot CIR modulation results

# Import libraries
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
import argparse
import os

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

GCR_labels = ["full", "nadv", "ndft", "ndif"]
n_time = 65
n_GCR = len(GCR_labels)
gcr_tim = [np.zeros(n_time) for _ in range(n_GCR)]
gcr_rat = [np.zeros(n_time) for _ in range(n_GCR)]
gcr_res = [np.zeros(n_time) for _ in range(n_GCR)]
gcr_lat = [np.zeros([90,2]) for _ in range(n_GCR)]
gcr_lon = [np.zeros([180,2]) for _ in range(n_GCR)]
epochs = np.linspace(args.zero_epoch-4.0, args.zero_epoch+4.0, num=n_time)
for i in range(n_GCR):
# Spectrum and residence time
   gcr_dat = np.loadtxt("output_{:s}/GCR_{:s}/cir_gcr_mod_sim_rate_{:s}.dat".format(args.date, GCR_labels[i], args.date))
   gcr_tim[i] = gcr_dat[:,0] - args.zero_epoch
   gcr_rat[i] = percentage_change(gcr_dat[:,1])
   gcr_res[i] = gcr_dat[:,2]

# Exit positions
   for j in range(n_time):
      gcr_lat_t = np.loadtxt("output_{:s}/GCR_{:s}/cir_gcr_mod_{:.3f}_lat_pp.dat".format(args.date, GCR_labels[i], epochs[j]))
      gcr_lat[i][:,1] += gcr_lat_t[:,1]
      gcr_lon_t = np.loadtxt("output_{:s}/GCR_{:s}/cir_gcr_mod_{:.3f}_lon_pp.dat".format(args.date, GCR_labels[i], epochs[j]))
      gcr_lon[i][:,1] += gcr_lon_t[:,1]
   gcr_lat[i][:,0] = gcr_lat_t[:,0]
   gcr_lat[i][:,1] /= n_time
   gcr_lon[i][:,0] = gcr_lon_t[:,0]
   gcr_lon[i][:,1] /= n_time

# Solar wind data
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

# Fieldlines at 5au (if present)
fls_5au = False
N_fl = 400
if os.path.isfile("output_{:s}/fieldlines/fieldline_{:d}.lines".format(args.date, N_fl-1)):
   fls_5au = True
   lat_5au = np.zeros(N_fl)
   lon_5au = np.zeros(N_fl)
   for i in range(N_fl):
      file = open("output_{:s}/fieldlines/fieldline_{:d}.lines".format(args.date, i))
      lines = file.readlines()
      file.close()
      pos = lines[-1].split(",")
      x = float(pos[0])
      y = float(pos[1])
      z = float(pos[2])
      r = x**2 + y**2 + z**2
      if (r < 5.0):
         print("ERROR: r < 5 au")
      lat_5au[i] = np.rad2deg(0.5 * np.pi - np.acos(z / r))
      lon_5au[i] = np.rad2deg(np.atan2(y, x))

# Plot GCR rate
fig = plt.figure(figsize=(14, 8), layout='tight')
ax = fig.add_subplot(111, projection='rectilinear')

marker_styles=["o","^","s","X"]
plot_colors=["tab:blue", "deepskyblue", "aqua", "navy"]
line_styles=["-","--","-.",":"]
for i in range(n_GCR):
   ax.scatter(gcr_tim[i], gcr_rat[i], s=50, c=plot_colors[i],
              marker=marker_styles[i], label=GCR_labels[i])
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
plt.savefig("output_{:s}/gcr_rate_{:s}.png".format(args.date, args.date))
plt.close(fig)

# Get extreme GCR responses
gcr_max = np.max(gcr_rat[0][0:n_time // 2 + 1])
gcr_min = np.min(gcr_rat[0][n_time // 2: -1])
print("Extreme values of percentage change:")
print(" max[-4,0]      min[0,4] ")
print(" ---------      -------- ")
print("    {:4.2f}         {:4.2f} ".format(gcr_max, gcr_min))

# Get correlations
pcc_nadv = stats.pearsonr(gcr_rat[0], gcr_rat[1])
pcc_ndft = stats.pearsonr(gcr_rat[0], gcr_rat[2])
pcc_ndif = stats.pearsonr(gcr_rat[0], gcr_rat[3])
print("Pearson correlation coefficients relative to full simulation:")
print("   nadv      ndft      ndif")
print("  ------    ------    ------")
print("   {:4.2f}      {:4.2f}     {:4.2f}".format(pcc_nadv.statistic, pcc_ndft.statistic, pcc_ndif.statistic))

# Get mean absolute error
mae_nadv = np.mean(np.abs(gcr_rat[0] - gcr_rat[1]))
mae_ndft = np.mean(np.abs(gcr_rat[0] - gcr_rat[2]))
mae_ndif = np.mean(np.abs(gcr_rat[0] - gcr_rat[3]))
print("Mean absolute error relative to full simulation:")
print("   nadv      ndft      ndif")
print("  ------    ------    ------")
print("   {:4.2f}      {:4.2f}     {:4.2f}".format(mae_nadv, mae_ndft, mae_ndif))

# Plot GCR residence time
fig = plt.figure(figsize=(14, 8), layout='tight')
ax = fig.add_subplot(111, projection='rectilinear')

for i in range(n_GCR):
   ax.scatter(gcr_tim[i], gcr_res[i], s=50, c=plot_colors[i],
              marker=marker_styles[i], label=GCR_labels[i])
ax.set_xlabel("time (days)", fontsize = 20)
ax.set_ylabel("residence time (days)", fontsize = 20)
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
plt.savefig("output_{:s}/gcr_res_time_{:s}.png".format(args.date, args.date))
plt.close(fig)

# Plot GCR exit latitude and longitude
fig = plt.figure(figsize=(14, 8), layout='tight')
ax = fig.add_subplot(211, projection='rectilinear')

for i in range(n_GCR):
   ax.plot(gcr_lat[i][:,0], gcr_lat[i][:,1],
           linewidth=3, color=plot_colors[i],
           linestyle=line_styles[i], label=GCR_labels[i])
if fls_5au:
   ax.axvspan(np.min(lat_5au), np.max(lat_5au), color='tab:orange', alpha=0.5)
ax.set_xlabel("latitude ($^\\circ$)", fontsize = 20)
ax.set_ylabel("density", fontsize = 20)
ax.set_xlim(-90.0,90.0)
ax.tick_params(labelsize=20)
ax.legend(fontsize=20)
ax.set_title("CIR {:s}".format(args.date), fontsize=24)

ax = fig.add_subplot(212, projection='rectilinear')

for i in range(n_GCR):
   ax.plot(np.concatenate((gcr_lon[i][90:,0]-360, gcr_lon[i][:90,0])), np.concatenate((gcr_lon[i][90:,1], gcr_lon[i][:90,1])),
           linewidth=3, color=plot_colors[i],
           linestyle=line_styles[i], label=GCR_labels[i])
if fls_5au:
   ax.axvspan(np.min(lon_5au), np.max(lon_5au), color='tab:orange', alpha=0.5)
ax.set_xlabel("longitude ($^\\circ$)", fontsize = 20)
ax.set_ylabel("density", fontsize = 20)
ax.set_xlim(-180.0,180.0)
ax.tick_params(labelsize=20)
ax.legend(fontsize=20)

plt.savefig("output_{:s}/exit_lat_lon_{:s}.png".format(args.date, args.date))
plt.close(fig)