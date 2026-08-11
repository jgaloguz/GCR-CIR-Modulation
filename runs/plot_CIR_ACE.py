# Code to plot simulated CIR along with measured ACE profile

# Import libraries
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm
import scipy.signal as spsig
import argparse

# Parse arguments
parser = argparse.ArgumentParser(description="Plot 1D time cuts of simulated and observed solar wind quantities")
parser.add_argument("date", type=str, help="Date of CIR in YYYYMMDD format")
parser.add_argument("zero_epoch", type=float, help="Time in days relative to CIR observed date to set as the zero epoch in plot. Possible values are any number in range [-10.0,10.0].")
args = parser.parse_args()

# Constants
m_p = 1.67262192e-24

# Import ACE data
ACE_SW = np.loadtxt("data/SI_{:s}_SW.txt".format(args.date))
ACE_SW[:,4] = ACE_SW[:,4] / 1000
ACE_SW[:,5] = m_p * ACE_SW[:,3] * ACE_SW[:,5] / 1.0e-10 # G^2 --> nT

labels = ["$|B|$ (nT)", "$|u|$ (km/s)", "$n$ (amu/cc)", "$T$ (kK)", "$\\langle \\delta B^2\\rangle$ (nT$^2$)"]

# Import sim data
Z = np.loadtxt("output_{:s}/CIR/mag_1au_{:s}.dat".format(args.date, args.date))
N = np.size(Z, 0)
idx_l = int(N * (14.0 + args.zero_epoch - 4.0) / 28.0)
idx_h = int(N * (14.0 + args.zero_epoch + 4.0) / 28.0)
SIM_SW = np.zeros((idx_h-idx_l,6))

SIM_SW[:,0] = Z[idx_l:idx_h,0] - args.zero_epoch
SIM_SW[:,1] = Z[idx_l:idx_h,1] / 1.0e-5 # G --> nT
Z = np.loadtxt("output_{:s}/CIR/vel_1au_{:s}.dat".format(args.date, args.date))
SIM_SW[:,2] = Z[idx_l:idx_h,1] / 1.0e5 # cm/s --> km/s
Z = np.loadtxt("output_{:s}/CIR/den_1au_{:s}.dat".format(args.date, args.date))
SIM_SW[:,3] = Z[idx_l:idx_h,1] 
Z = np.loadtxt("output_{:s}/CIR/pth_1au_{:s}.dat".format(args.date, args.date))
SIM_SW[:,4] = Z[idx_l:idx_h,1] / 10.0 / SIM_SW[:,3] / 1.0e6 / 1.380649e-23 / 1.0e3 # T = P / (n*k_B) and dyn/cm^2 --> Pa and cm^-3 --> m^-3 and K --> kK
Z = np.loadtxt("output_{:s}/CIR/tur_enr_1au_{:s}.dat".format(args.date, args.date))
SIM_SW[:,5] = m_p * SIM_SW[:,3] * Z[idx_l:idx_h,1] / 1.0e-10 # G^2 --> nT

# Make plots
fig = plt.figure(figsize=(14, 8), layout='tight')
ax = fig.add_subplot(111, projection='rectilinear')

# GCR rate
fig, axs = plt.subplots(nrows=5, ncols=1, figsize=(10, 18))

for i in range(5):
   if i == 4:
      plot_func = axs[i].semilogy
   else:
      plot_func = axs[i].plot
   plot_func(ACE_SW[:,0], ACE_SW[:,i+1], label="ACE")
   plot_func(SIM_SW[:,0], SIM_SW[:,i+1], label="AWSoM")
   axs[i].set_ylabel(labels[i], fontsize=20)
   axs[i].tick_params(labelsize=16)
   if i != 1 and i != 4:
      axs[i].set_ylim(bottom=0.0)
   axs[i].legend(fontsize=20)
   axs[i].set_xlim(-4.0,4.0)

axs[4].set_xlabel("time (days)", fontsize=20)
axs[0].set_title("CIR {:s}".format(args.date), fontsize=24)
plt.tight_layout()
plt.savefig("output_{:s}/sim_vs_ACE_{:s}.png".format(args.date, args.date), dpi=200)
# plt.show()
