# Code to plot CIR modulation results

# Import libraries
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm
import scipy.signal as spsig
import sys

# Constants
m_p = 1.67262192e-24

# Check CIR date
if len(sys.argv) > 1:
   cir_date = sys.argv[1]
   print("CIR date: {:s}".format(cir_date))
else:
   print("ERROR: No CIR date provided.")
   sys.exit(1)

# Import ACE data
ACE_SW = np.loadtxt("data/SI_{:s}_SW.txt".format(cir_date))
ACE_SW[:,4] = ACE_SW[:,4] / 1000
ACE_SW[:,5] = m_p * ACE_SW[:,3] * ACE_SW[:,5] * 1.0e9

labels = ["$|B|$ (nT)", "$|u|$ (km/s)", "$n$ (amu/cc)", "$T$ (kK)", "$Z^2$ (nerg cm$^{-3}$)"]

# Import sim data
CIR_epoch = 11.8
idx_l = 297
idx_h = 593

SIM_SW = np.zeros((idx_h-idx_l,6))
Z = np.loadtxt("output_{:s}/CIR/mag_1au_{:s}.dat".format(cir_date, cir_date))
SIM_SW[:,0] = Z[idx_l:idx_h,0] - CIR_epoch
SIM_SW[:,1] = Z[idx_l:idx_h,1] / 1.0e-5 # G --> nT
Z = np.loadtxt("output_{:s}/CIR/vel_1au_{:s}.dat".format(cir_date, cir_date))
SIM_SW[:,2] = Z[idx_l:idx_h,1] / 1.0e5 # cm/s --> km/s
Z = np.loadtxt("output_{:s}/CIR/den_1au_{:s}.dat".format(cir_date, cir_date))
SIM_SW[:,3] = Z[idx_l:idx_h,1] 
Z = np.loadtxt("output_{:s}/CIR/pth_1au_{:s}.dat".format(cir_date, cir_date))
SIM_SW[:,4] = Z[idx_l:idx_h,1] / 10.0 / SIM_SW[:,3] / 1.0e6 / 1.380649e-23 / 1.0e3 # T = P / (n*k_B) and dyn/cm^2 --> Pa and cm^-3 --> m^-3 and K --> kK
Z = np.loadtxt("output_{:s}/CIR/tur_enr_1au_{:s}.dat".format(cir_date, cir_date))
SIM_SW[:,5] = m_p * SIM_SW[:,3] * Z[idx_l:idx_h,1] * 1.0e9

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
   plot_func(SIM_SW[:,0], SIM_SW[:,i+1], label="sim")
   axs[i].set_ylabel(labels[i], fontsize=20)
   axs[i].tick_params(labelsize=16)
   if i != 1 and i != 4:
      axs[i].set_ylim(bottom=0.0)
   axs[i].legend(fontsize=20)
   axs[i].set_xlim(-4.0,4.0)

plt.tight_layout()
plt.savefig("output_{:s}/sim_vs_ACE{:s}.png".format(cir_date, cir_date), dpi=200)
# plt.show()
