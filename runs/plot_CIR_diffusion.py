# Code to plot CIR diffusion coefficients

# Import libraries
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm
import scipy.signal as spsig
import argparse

# Parse arguments
parser = argparse.ArgumentParser(description="Plot 1D time cuts of simulated diffusion coefficients")
parser.add_argument("date", type=str, help="Date of CIR in YYYYMMDD format")
parser.add_argument("zero_epoch", type=float, help="Time in days relative to CIR observed date to set as the zero epoch in plot. Possible values are any number in range [-10.0,10.0].")
args = parser.parse_args()


# Import sim data
Z = np.loadtxt("output_{:s}/CIR/diff1_1au_{:s}.dat".format(args.date, args.date))
N = np.size(Z, 0)
idx_l = int(N * (14.0 + args.zero_epoch - 4.0) / 28.0)
idx_h = int(N * (14.0 + args.zero_epoch + 4.0) / 28.0)
diff = np.zeros((idx_h-idx_l,3))

diff[:,0] = Z[idx_l:idx_h,0] - args.zero_epoch
diff[:,1] = Z[idx_l:idx_h,1] / 1.0e21
Z = np.loadtxt("output_{:s}/CIR/diff2_1au_{:s}.dat".format(args.date, args.date))
diff[:,2] = Z[idx_l:idx_h,1] / 1.0e23

# Make plots
fig = plt.figure(figsize=(14, 8), layout='tight')
ax = fig.add_subplot(111, projection='rectilinear')

# Perpendicular diffusion
ax.plot(diff[:,0], diff[:,1],
        linewidth=3, color='tab:blue')
ax.set_xlabel("time (days)", fontsize = 20)
ax.set_ylabel("$\\kappa_{\\perp}$ (10$^{21}$ cm$^{2}$ s$^{-1}$)", fontsize = 20, color='tab:blue')
ax.set_xlim(diff[0,0],diff[-1,0])
ax.set_ylim(bottom=0.0)
ax.tick_params(labelsize=20)
ax.tick_params(axis='y', labelcolor="tab:blue")

# Parallel diffusion
axx = ax.twinx()
axx.plot(diff[:,0], diff[:,2],
        linewidth=3, color='tab:orange')
axx.set_xlabel("time (days)", fontsize = 20)
axx.set_ylabel("$\\kappa_{\\parallel}$ (10$^{23}$ cm$^{2}$ s$^{-1}$)", fontsize = 20, color='tab:orange')
axx.set_xlim(diff[0,0],diff[-1,0])
axx.set_ylim(bottom=0.0)
axx.tick_params(labelsize=20)
axx.tick_params(axis='y', labelcolor="tab:orange")

ax.set_title("CIR {:s}".format(args.date), fontsize=24)
plt.savefig("output_{:s}/sim_diff_{:s}.png".format(args.date, args.date), dpi=200)
# plt.show()
