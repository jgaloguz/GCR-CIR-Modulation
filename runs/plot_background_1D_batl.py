# Code to plot 1D time cuts of solar wind and particle quantities

# Import libraries
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm
import argparse

# Parse arguments
parser = argparse.ArgumentParser(description="Plot 1D time cuts of solar wind and particle quantities")
parser.add_argument("date", type=str, default="YYYYMMDD", help="Date of CIR in YYYYMMDD format")
parser.add_argument("--zero_epoch", type=float, default=-1000.0, help="Time in days relative to CIR observed date to set as the zero epoch in plot. Possible values are any number in range [-10.0,10.0]. If omitted, plots are shown between -14 and +14 days of simulation. If present, plots are shown between epoch-4 and epoch+4 days of simulation.")
args = parser.parse_args()

plot_den = True
plot_vel = True
plot_mag = True
plot_tmp = True
plot_pol = True
plot_dmax = True
plot_tur_enr = True
plot_het_flx = True
plot_rad_vel = True
plot_div_vel = True
plot_drift = True
plot_diff1 = True
plot_diff2 = True
plot_diff3 = True

quants = 0
data1DX = []
data1DY = []
data1DL = []

directions = ["100","-100","010","0-10","001","00-1"]

if -10.0 <= args.zero_epoch and args.zero_epoch <= 10.0:
   epoch1 = args.zero_epoch - 4.0
   epoch2 = args.zero_epoch + 4.0
else:
   epoch1 = -14.0
   epoch2 = 14.0

print("====================")
print("Plotting 1D cuts:")
print("====================")

# Create the pcolormesh plots
# ==================================================
if plot_den:
   quants += 1
   print("Plotting plasma density")

   Z = np.loadtxt("output_{:s}/CIR/den_1au_{:s}.dat".format(args.date, args.date))
   plt.plot(Z[:,0], Z[:,1])
   plt.title('Plasma Density (1 au)')
   plt.xlabel('days')
   plt.ylabel('n (amu/cc)')
   plt.xlim(epoch1, epoch2)

   data1DX.append(Z[:,0])
   data1DY.append(Z[:,1])
   data1DL.append('n (amu/cc)')

   plt.savefig("output_{:s}/CIR/den_1au_{:s}.png".format(args.date, args.date), dpi=200)
   plt.close()

   for direction in directions:
      Z = np.loadtxt("output_{:s}/CIR/den_{:s}_{:s}.dat".format(args.date, direction, args.date))
      plt.loglog(Z[:,0], Z[:,1], label=direction)
   plt.title('Plasma Density')
   plt.xlabel('r (au)')
   plt.ylabel('n (amu/cc)')
   plt.xlim(Z[0,0], Z[-1,0])
   plt.legend()

   plt.savefig("output_{:s}/CIR/den_cuts_{:s}.png".format(args.date, args.date), dpi=200)
   plt.close()

# ==================================================
if plot_vel:
   quants += 1
   print("Plotting plasma flow")

   Z = np.loadtxt("output_{:s}/CIR/vel_1au_{:s}.dat".format(args.date, args.date))
   Z[:,1] = Z[:,1] / 1.0e5 # cm/s --> km/s
   plt.plot(Z[:,0], Z[:,1])
   plt.title('Plasma Flow Magnitude (1 au)')
   plt.xlabel('days')
   plt.ylabel('|u| (km/s)')
   plt.xlim(epoch1, epoch2)

   data1DX.append(Z[:,0])
   data1DY.append(Z[:,1])
   data1DL.append('|u| (km/s)')

   plt.savefig("output_{:s}/CIR/vel_1au_{:s}.png".format(args.date, args.date), dpi=200)
   plt.close()

   for direction in directions:
      Z = np.loadtxt("output_{:s}/CIR/vel_{:s}_{:s}.dat".format(args.date, direction, args.date))
      Z[:,1] = Z[:,1] / 1.0e5 # cm/s --> km/s
      plt.semilogx(Z[:,0], Z[:,1], label=direction)
   plt.title('Plasma Flow Magnitude')
   plt.xlabel('r (au)')
   plt.ylabel('|u| (km/s)')
   plt.xlim(Z[0,0], Z[-1,0])
   plt.legend()

   plt.savefig("output_{:s}/CIR/vel_cuts_{:s}.png".format(args.date, args.date), dpi=200)
   plt.close()

# ==================================================
if plot_mag:
   quants += 1
   print("Plotting magnetic field")

   Z = np.loadtxt("output_{:s}/CIR/mag_1au_{:s}.dat".format(args.date, args.date))
   Z[:,1] = Z[:,1] / 1.0e-5 # G --> nT
   plt.plot(Z[:,0], Z[:,1])
   plt.title('Magnetic Field Magnitude (1 au)')
   plt.xlabel('days')
   plt.ylabel('|B| (nT)')
   plt.xlim(epoch1, epoch2)

   data1DX.append(Z[:,0])
   data1DY.append(Z[:,1])
   data1DL.append('|B| (nT)')

   plt.savefig("output_{:s}/CIR/mag_1au_{:s}.png".format(args.date, args.date), dpi=200)
   plt.close()

   for direction in directions:
      Z = np.loadtxt("output_{:s}/CIR/mag_{:s}_{:s}.dat".format(args.date, direction, args.date))
      Z[:,1] = Z[:,1] / 1.0e-5 # G --> nT
      plt.loglog(Z[:,0], Z[:,1], label=direction)
   plt.title('Magnetic Field Magnitude')
   plt.xlabel('r (au)')
   plt.ylabel('|B| (nT)')
   plt.xlim(Z[0,0], Z[-1,0])
   plt.legend()

   plt.savefig("output_{:s}/CIR/mag_cuts_{:s}.png".format(args.date, args.date), dpi=200)
   plt.close()

# ==================================================
if plot_tmp:
   quants += 1
   print("Plotting plasma temperature")

   Z = np.loadtxt("output_{:s}/CIR/pth_1au_{:s}.dat".format(args.date, args.date))
   P = np.copy(Z)
   P[:,1] = P[:,1] * 0.1 # dyn/cm^2 --> Pa
   n = np.loadtxt("output_{:s}/CIR/den_1au_{:s}.dat".format(args.date, args.date))
   n[:,1] = n[:,1] * 1.0e6 # cm^-3 --> m^-3
   Z[:,1] = P[:,1] / n[:,1] / 1.380649e-23 / 1.0e6 # T = P / (n*k_B) and K --> MK
   Z[:,0] = P[:,0]
   plt.plot(Z[:,0], Z[:,1])
   plt.title('Plasma Temperature (1 au)')
   plt.xlabel('days')
   plt.ylabel('T (MK)')
   plt.xlim(epoch1, epoch2)

   data1DX.append(Z[:,0])
   data1DY.append(Z[:,1])
   data1DL.append('T (MK)')

   plt.savefig("output_{:s}/CIR/tmp_1au_{:s}.png".format(args.date, args.date), dpi=200)
   plt.close()

   for direction in directions:
      Z = np.loadtxt("output_{:s}/CIR/pth_{:s}_{:s}.dat".format(args.date, direction, args.date))
      P = np.copy(Z)
      P[:,1] = P[:,1] * 0.1 # dyn/cm^2 --> Pa
      n = np.loadtxt("output_{:s}/CIR/den_{:s}_{:s}.dat".format(args.date, direction, args.date))
      n[:,1] = n[:,1] * 1.0e6 # cm^-3 --> m^-3
      Z[:,1] = P[:,1] / n[:,1] / 1.380649e-23 / 1.0e6 # T = P / (n*k_B) and K --> MK
      Z[:,0] = P[:,0]
      plt.loglog(Z[:,0], Z[:,1], label=direction)
   plt.title('Plasma Temperature')
   plt.xlabel('r (au)')
   plt.ylabel('T (MK)')
   plt.xlim(Z[0,0], Z[-1,0])
   plt.legend()

   plt.savefig("output_{:s}/CIR/tmp_cuts_{:s}.png".format(args.date, args.date), dpi=200)
   plt.close()

# ==================================================
if plot_pol:
   quants += 1
   print("Plotting magnetic polarity")

   Z = np.loadtxt("output_{:s}/CIR/pol_1au_{:s}.dat".format(args.date, args.date))
   plt.plot(Z[:,0], Z[:,1])
   plt.title('Magnetic Field Polarity (1 au)')
   plt.xlabel('days')
   plt.ylabel('sgn(B) (+out/-in)')
   plt.xlim(epoch1, epoch2)

   data1DX.append(Z[:,0])
   data1DY.append(Z[:,1])
   data1DL.append('sgn(B) (+out/-in)')

   plt.savefig("output_{:s}/CIR/pol_1au_{:s}.png".format(args.date, args.date), dpi=200)
   plt.close()

   for direction in directions:
      Z = np.loadtxt("output_{:s}/CIR/pol_{:s}_{:s}.dat".format(args.date, direction, args.date))
      plt.semilogx(Z[:,0], Z[:,1], label=direction)
   plt.title('Magnetic Field Polarity')
   plt.xlabel('r (au)')
   plt.ylabel('sgn(B) (+out/-in)')
   plt.xlim(Z[0,0], Z[-1,0])
   plt.legend()

   plt.savefig("output_{:s}/CIR/pol_cuts_{:s}.png".format(args.date, args.date), dpi=200)
   plt.close()

# ==================================================
if plot_dmax:
   quants += 1
   print("Plotting cell size")

   Z = np.loadtxt("output_{:s}/CIR/dmax_1au_{:s}.dat".format(args.date, args.date))
   plt.plot(Z[:,0], Z[:,1])
   plt.title('Cell Size (1 au)')
   plt.xlabel('days')
   plt.ylabel('$\\Delta x$ (au)')
   plt.xlim(epoch1, epoch2)

   data1DX.append(Z[:,0])
   data1DY.append(Z[:,1])
   data1DL.append('$\\Delta x$ (au)')

   plt.savefig("output_{:s}/CIR/dmax_1au_{:s}.png".format(args.date, args.date), dpi=200)
   plt.close()

   for direction in directions:
      Z = np.loadtxt("output_{:s}/CIR/dmax_{:s}_{:s}.dat".format(args.date, direction, args.date))
      plt.semilogx(Z[:,0], Z[:,1], label=direction)
   plt.title('Cell Size')
   plt.xlabel('r (au)')
   plt.ylabel('$\\Delta x$ (au)')
   plt.xlim(Z[0,0], Z[-1,0])
   plt.legend()

   plt.savefig("output_{:s}/CIR/dmax_cuts_{:s}.png".format(args.date, args.date), dpi=200)
   plt.close()

# ==================================================
if plot_tur_enr:
   quants += 1
   print("Plotting turbulent energy")

   Z = np.loadtxt("output_{:s}/CIR/tur_enr_1au_{:s}.dat".format(args.date, args.date))
   plt.plot(Z[:,0], Z[:,1])
   plt.title('Alfv\'enic Fluctuation Energy (1 au)')
   plt.xlabel('days')
   plt.ylabel('$Z_A^2$ (cm$^2$ s$^{-2}$)')
   plt.xlim(epoch1, epoch2)

   data1DX.append(Z[:,0])
   data1DY.append(Z[:,1])
   data1DL.append('$Z^2$ (cm$^2$ s$^{-2}$)')

   plt.savefig("output_{:s}/CIR/tur_enr_1au_{:s}.png".format(args.date, args.date), dpi=200)
   plt.close()

   for direction in directions:
      Z = np.loadtxt("output_{:s}/CIR/tur_enr_{:s}_{:s}.dat".format(args.date, direction, args.date))
      plt.loglog(Z[:,0], Z[:,1], label=direction)
   plt.title('Alfv\'enic Fluctuation Energy')
   plt.xlabel('r (au)')
   plt.ylabel('$Z_A^2$ (cm$^2$ s$^{-2}$)')
   plt.xlim(Z[0,0], Z[-1,0])
   plt.legend()

   plt.savefig("output_{:s}/CIR/tur_enr_cuts_{:s}.png".format(args.date, args.date), dpi=200)
   plt.close()

# ==================================================
if plot_het_flx:
   quants += 1
   print("Plotting heat flux speed")

   Z = np.loadtxt("output_{:s}/CIR/het_flx_1au_{:s}.dat".format(args.date, args.date))
   plt.plot(Z[:,0], Z[:,1])
   plt.title('Heat Flux (1 au)')
   plt.xlabel('days')
   plt.ylabel('$Q_e$ (cm$^2$ s$^{-2}$)')
   plt.xlim(epoch1, epoch2)

   data1DX.append(Z[:,0])
   data1DY.append(Z[:,1])
   data1DL.append('$Q_e$ (cm$^2$ s$^{-2}$)')

   plt.savefig("output_{:s}/CIR/het_flx_1au_{:s}.png".format(args.date, args.date), dpi=200)
   plt.close()

   for direction in directions:
      Z = np.loadtxt("output_{:s}/CIR/het_flx_{:s}_{:s}.dat".format(args.date, direction, args.date))
      plt.loglog(Z[:,0], Z[:,1], label=direction)
   plt.title('Heat Flux')
   plt.xlabel('r (au)')
   plt.ylabel('$Q_e$ (cm$^2$ s$^{-2}$)')
   plt.xlim(Z[0,0], Z[-1,0])
   plt.legend()

   plt.savefig("output_{:s}/CIR/het_flx_cuts_{:s}.png".format(args.date, args.date), dpi=200)
   plt.close()

# ==================================================
if plot_rad_vel:
   quants += 1
   print("Plotting fraction of radial plasma flow")

   Z = np.loadtxt("output_{:s}/CIR/rad_vel_1au_{:s}.dat".format(args.date, args.date))
   plt.plot(Z[:,0], Z[:,1])
   plt.title('Fraction of Radial Plasma Flow (1 au)')
   plt.xlabel('days')
   plt.ylabel('$\\hat{r} \\cdot \\hat{u}$')
   plt.xlim(epoch1, epoch2)

   data1DX.append(Z[:,0])
   data1DY.append(Z[:,1])
   data1DL.append('$\\hat{r} \\cdot \\hat{u}$')

   plt.savefig("output_{:s}/CIR/rad_vel_1au_{:s}.png".format(args.date, args.date), dpi=200)
   plt.close()

   for direction in directions:
      Z = np.loadtxt("output_{:s}/CIR/rad_vel_{:s}_{:s}.dat".format(args.date, direction, args.date))
      plt.semilogx(Z[:,0], Z[:,1], label=direction)
   plt.title('Fraction of Radial Plasma Flow')
   plt.xlabel('r (au)')
   plt.ylabel('$\\hat{r} \\cdot \\hat{u}$')
   plt.xlim(Z[0,0], Z[-1,0])
   plt.legend()

   plt.savefig("output_{:s}/CIR/rad_vel_cuts_{:s}.png".format(args.date, args.date), dpi=200)
   plt.close()

# ==================================================
if plot_div_vel:
   quants += 1
   print("Plotting divergence of plasma flow")

   Z = np.loadtxt("output_{:s}/CIR/div_vel_1au_{:s}.dat".format(args.date, args.date))
   plt.plot(Z[:,0], Z[:,1])
   plt.title('Divergence of Plasma Flow (1 au)')
   plt.xlabel('days')
   plt.ylabel('$\\nabla \\cdot u$ (km/s)')
   plt.xlim(epoch1, epoch2)

   data1DX.append(Z[:,0])
   data1DY.append(Z[:,1])
   data1DL.append('$\\nabla \\cdot u$ (1/s)')

   plt.savefig("output_{:s}/CIR/div_vel_1au_{:s}.png".format(args.date, args.date), dpi=200)
   plt.close()

   for direction in directions:
      Z = np.loadtxt("output_{:s}/CIR/div_vel_{:s}_{:s}.dat".format(args.date, direction, args.date))
      plt.loglog(Z[:,0], Z[:,1], label=direction)
   plt.title('Divergence of Plasma Flow')
   plt.xlabel('r (au)')
   plt.ylabel('$\\nabla \\cdot u$ (km/s)')
   plt.xlim(Z[0,0], Z[-1,0])
   plt.legend()

   plt.savefig("output_{:s}/CIR/div_vel_cuts_{:s}.png".format(args.date, args.date), dpi=200)
   plt.close()

# ==================================================
if plot_drift:
   quants += 1
   print("Plotting drift speed")

   Z = np.loadtxt("output_{:s}/CIR/drift_1au_{:s}.dat".format(args.date, args.date))
   plt.plot(Z[:,0], Z[:,1])
   plt.title('Magnetic Drift Speed (1 au)')
   plt.xlabel('days')
   plt.ylabel('$v_d$ (/ $v_p$)')
   plt.xlim(epoch1, epoch2)

   data1DX.append(Z[:,0])
   data1DY.append(Z[:,1])
   data1DL.append('$v_d$ (/ $v_p$)')

   plt.savefig("output_{:s}/CIR/drift_1au_{:s}.png".format(args.date, args.date), dpi=200)
   plt.close()

   for direction in directions:
      Z = np.loadtxt("output_{:s}/CIR/drift_{:s}_{:s}.dat".format(args.date, direction, args.date))
      plt.semilogx(Z[:,0], Z[:,1], label=direction)
   plt.title('Magnetic Drift Speed')
   plt.xlabel('r (au)')
   plt.ylabel('$v_d$ (/ $v_p$)')
   plt.xlim(Z[0,0], Z[-1,0])
   plt.legend()

   plt.savefig("output_{:s}/CIR/drift_cuts_{:s}.png".format(args.date, args.date), dpi=200)
   plt.close()

# ==================================================
if plot_diff1:
   quants += 1
   print("Plotting diffusion 1 (perpendicular)")

   Z = np.loadtxt("output_{:s}/CIR/diff1_1au_{:s}.dat".format(args.date, args.date))
   plt.plot(Z[:,0], Z[:,1])
   plt.title('Diffusion Model Perpendicular (1 au)')
   plt.xlabel('days')
   plt.ylabel('$\\kappa_{\\perp}$')
   plt.xlim(epoch1, epoch2)

   data1DX.append(Z[:,0])
   data1DY.append(Z[:,1])
   data1DL.append('$\\kappa_{\\perp}$')

   plt.savefig("output_{:s}/CIR/diff1_1au_{:s}.png".format(args.date, args.date), dpi=200)
   plt.close()

   for direction in directions:
      Z = np.loadtxt("output_{:s}/CIR/diff1_{:s}_{:s}.dat".format(args.date, direction, args.date))
      plt.loglog(Z[:,0], Z[:,1], label=direction)
   plt.title('Diffusion Model Perpendicular')
   plt.xlabel('r (au)')
   plt.ylabel('$\\kappa_{\\perp}$')
   plt.xlim(Z[0,0], Z[-1,0])
   plt.legend()

   plt.savefig("output_{:s}/CIR/diff1_cuts_{:s}.png".format(args.date, args.date), dpi=200)
   plt.close()

   for direction in directions:
      Z = np.loadtxt("output_{:s}/CIR/diff1_E_{:s}_{:s}.dat".format(args.date, direction, args.date))
      Z[:,0] = Z[:,0] * 200.0 # Unitless --> MeV
      plt.loglog(Z[:,0], Z[:,1], label=direction)
   plt.title('Diffusion Model Perpendicular')
   plt.xlabel('E (MeV)')
   plt.ylabel('$\\kappa_{\\perp}$')
   plt.xlim(Z[0,0], Z[-1,0])
   plt.legend()

   plt.savefig("output_{:s}/CIR/diff1_E_cuts_{:s}.png".format(args.date, args.date), dpi=200)
   plt.close()

# ==================================================
if plot_diff2:
   quants += 1
   print("Plotting diffusion 2 (parallel)")

   Z = np.loadtxt("output_{:s}/CIR/diff2_1au_{:s}.dat".format(args.date, args.date))
   plt.plot(Z[:,0], Z[:,1])
   plt.title('Diffusion Model Parallel (1 au)')
   plt.xlabel('days')
   plt.ylabel('$\\kappa_{\\parallel}$')
   plt.xlim(epoch1, epoch2)

   data1DX.append(Z[:,0])
   data1DY.append(Z[:,1])
   data1DL.append('$\\kappa_{\\parallel}$')

   plt.savefig("output_{:s}/CIR/diff2_1au_{:s}.png".format(args.date, args.date), dpi=200)
   plt.close()

   for direction in directions:
      Z = np.loadtxt("output_{:s}/CIR/diff2_{:s}_{:s}.dat".format(args.date, direction, args.date))
      plt.loglog(Z[:,0], Z[:,1], label=direction)
   plt.title('Diffusion Model Parallel')
   plt.xlabel('r (au)')
   plt.ylabel('$\\kappa_{\\parallel}$')
   plt.xlim(Z[0,0], Z[-1,0])
   plt.legend()

   plt.savefig("output_{:s}/CIR/diff2_cuts_{:s}.png".format(args.date, args.date), dpi=200)
   plt.close()

   for direction in directions:
      Z = np.loadtxt("output_{:s}/CIR/diff2_E_{:s}_{:s}.dat".format(args.date, direction, args.date))
      Z[:,0] = Z[:,0] * 200.0 # Unitless --> MeV
      plt.loglog(Z[:,0], Z[:,1], label=direction)
   plt.title('Diffusion Model Parallel')
   plt.xlabel('E (MeV)')
   plt.ylabel('$\\kappa_{\\parallel}$')
   plt.xlim(Z[0,0], Z[-1,0])
   plt.legend()

   plt.savefig("output_{:s}/CIR/diff2_E_cuts_{:s}.png".format(args.date, args.date), dpi=200)
   plt.close()

# ==================================================
if plot_diff3:
   quants += 1
   print("Plotting diffusion 3 (radial)")

   Z = np.loadtxt("output_{:s}/CIR/diff3_1au_{:s}.dat".format(args.date, args.date))
   plt.plot(Z[:,0], Z[:,1])
   plt.title('Diffusion Model Radial (1 au)')
   plt.xlabel('days')
   plt.ylabel('$\\kappa_{rr}$')
   plt.xlim(epoch1, epoch2)

   data1DX.append(Z[:,0])
   data1DY.append(Z[:,1])
   data1DL.append('$\\kappa_{rr}$')

   plt.savefig("output_{:s}/CIR/diff3_1au_{:s}.png".format(args.date, args.date), dpi=200)
   plt.close()

   for direction in directions:
      Z = np.loadtxt("output_{:s}/CIR/diff3_{:s}_{:s}.dat".format(args.date, direction, args.date))
      plt.loglog(Z[:,0], Z[:,1], label=direction)
   plt.title('Diffusion Model Radial')
   plt.xlabel('r (au)')
   plt.ylabel('$\\kappa_{rr}$')
   plt.xlim(Z[0,0], Z[-1,0])
   plt.legend()

   plt.savefig("output_{:s}/CIR/diff3_cuts_{:s}.png".format(args.date, args.date), dpi=200)
   plt.close()

# ==================================================
print("Plotting", quants, "1D quantities in composite figure")
if quants > 0:
   fig, axs = plt.subplots(nrows=quants, ncols=1, figsize=(10, 3*quants))

   for i in range(quants):
      axs[i].plot(data1DX[i], data1DY[i])
      axs[i].set_ylabel(data1DL[i])
      axs[i].set_xlim(epoch1, epoch2)

   axs[quants-1].set_xlabel('days')
   plt.tight_layout()
   plt.savefig("output_{:s}/CIR/composite_1au_{:s}.png".format(args.date, args.date), dpi=200)
   plt.show()