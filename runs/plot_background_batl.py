# Code to plot 2D slices and 1D time curves of solar wind and particle quantities

# Import libraries
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm
import sys

# Check CIR date
if len(sys.argv) > 1:
   cir_date = sys.argv[1]
   print("CIR date: {:s}".format(cir_date))
else:
   print("ERROR: No CIR date provided.")
   sys.exit(1)

plot_den = True
plot_vel = True
plot_rad_vel = True
plot_div_vel = True
plot_mag = True
plot_Blines = False
plot_pol = True
plot_dmax = True
plot_het_flx = True
plot_tur_enr = True
plot_drift = True
plot_diff1 = True
plot_diff2 = True

quants = 0
data1DX = []
data1DY = []
data1DL = []

# Find index for midpoint of curve
def midpoint(curve_xy):
   N = np.size(curve_xy,0)
   s = 0.0
   for seg in range(N-1):
      s += np.sqrt((curve_xy[seg+1,0]-curve_xy[seg,0])**2 \
                 + (curve_xy[seg+1,1]-curve_xy[seg,1])**2)
   hs = 0.5 * s
   s = 0.0
   for seg in range(N-1):
      s += np.sqrt((curve_xy[seg+1,0]-curve_xy[seg,0])**2 \
                 + (curve_xy[seg+1,1]-curve_xy[seg,1])**2)
      if s > hs:
         break
   return seg

# Create sample data
N = 1000
x = np.linspace(-5.0, 5.0, N)
y = np.linspace(-5.0, 5.0, N)
X, Y = np.meshgrid(x, y)

# Create the pcolormesh plots
# ==================================================
if plot_den:
   quants += 1
   print("Plotting plasma density")

   Z = np.loadtxt("output_{:s}/CIR/den_equ_{:s}.dat".format(cir_date, cir_date))
   cmap = plt.colormaps.get_cmap('magma')
   cmap.set_under("white")
   plt.pcolormesh(X, Y, np.transpose(Z), norm=LogNorm(vmin=np.min(Z[Z > 0.0]), vmax=np.max(Z)),
                           cmap=cmap, shading='gouraud')
   plt.colorbar(label='n (amu/cc)')
   plt.title('Plasma Density (z = 0)')
   plt.xlabel('x (au)')
   plt.ylabel('y (au)')

   plt.savefig("output_{:s}/CIR/den_equ_{:s}.png".format(cir_date, cir_date), dpi=200)
   plt.close()

   Z = np.loadtxt("output_{:s}/CIR/den_mer_{:s}.dat".format(cir_date, cir_date))
   cmap = plt.colormaps.get_cmap('magma')
   cmap.set_under("white")
   plt.pcolormesh(X, Y, np.transpose(Z), norm=LogNorm(vmin=np.min(Z[Z > 0.0]), vmax=np.max(Z)),
                           cmap=cmap, shading='gouraud')
   plt.colorbar(label='n (amu/cc)')
   plt.title('Plasma Density (y = 0)')
   plt.xlabel('x (au)')
   plt.ylabel('z (au)')

   plt.savefig("output_{:s}/CIR/den_mer_{:s}.png".format(cir_date, cir_date), dpi=200)
   plt.close()

   Z = np.loadtxt("output_{:s}/CIR/den_1au_{:s}.dat".format(cir_date, cir_date))
   plt.plot(Z[:,0], Z[:,1])
   plt.title('Plasma Density (1 au)')
   plt.xlabel('days')
   plt.ylabel('n (amu/cc)')
   plt.xlim(0, 27.0)

   data1DX.append(Z[:,0])
   data1DY.append(Z[:,1])
   data1DL.append('n (amu/cc)')

   plt.savefig("output_{:s}/CIR/den_1au_{:s}.png".format(cir_date, cir_date), dpi=200)
   plt.close()

# ==================================================
if plot_vel:
   quants += 1
   print("Plotting plasma flow")

   Z = np.loadtxt("output_{:s}/CIR/vel_equ_{:s}.dat".format(cir_date, cir_date)) / 1.0e5 # cm/s --> km/s
   cmap = plt.colormaps.get_cmap('plasma')
   cmap.set_under("white")
   plt.pcolormesh(X, Y, np.transpose(Z), vmin=np.min(Z[Z > 0.0]), vmax=np.max(Z),
                           cmap=cmap, shading='gouraud')
   plt.colorbar(label='|u| (km/s)')
   plt.title('Plasma Flow Magnitude (z = 0)')
   plt.xlabel('x (au)')
   plt.ylabel('y (au)')

   plt.savefig("output_{:s}/CIR/vel_equ_{:s}.png".format(cir_date, cir_date), dpi=200)
   plt.close()

   Z = np.loadtxt("output_{:s}/CIR/vel_mer_{:s}.dat".format(cir_date, cir_date)) / 1.0e5 # cm/s --> km/s
   cmap = plt.colormaps.get_cmap('plasma')
   cmap.set_under("white")
   plt.pcolormesh(X, Y, np.transpose(Z), vmin=np.min(Z[Z > 0.0]), vmax=np.max(Z),
                           cmap=cmap, shading='gouraud')
   plt.colorbar(label='|u| (km/s)')
   plt.title('Plasma Flow Magnitude (y = 0)')
   plt.xlabel('x (au)')
   plt.ylabel('z (au)')

   plt.savefig("output_{:s}/CIR/vel_mer_{:s}.png".format(cir_date, cir_date), dpi=200)
   plt.close()

   Z = np.loadtxt("output_{:s}/CIR/vel_1au_{:s}.dat".format(cir_date, cir_date))
   Z[:,1] = Z[:,1] / 1.0e5 # cm/s --> km/s
   plt.plot(Z[:,0], Z[:,1])
   plt.title('Plasma Flow Magnitude (km/s)')
   plt.xlabel('days')
   plt.ylabel('|u| (km/s)')
   plt.xlim(0, 27.0)

   data1DX.append(Z[:,0])
   data1DY.append(Z[:,1])
   data1DL.append('|u| (km/s)')

   plt.savefig("output_{:s}/CIR/vel_1au_{:s}.png".format(cir_date, cir_date), dpi=200)
   plt.close()

# ==================================================
if plot_rad_vel:
   quants += 1
   print("Plotting fraction of radial plasma flow")

   Z = np.loadtxt("output_{:s}/CIR/rad_vel_equ_{:s}.dat".format(cir_date, cir_date))
   cmap = plt.colormaps.get_cmap('plasma')
   cmap.set_under("white")
   plt.pcolormesh(X, Y, np.transpose(Z), vmin=np.min(Z[Z > 0.0]), vmax=np.max(Z),
                           cmap=cmap, shading='gouraud')
   plt.colorbar(label='$\\hat{r} \\cdot \\hat{u}$')
   plt.title('Fraction of Radial Plasma Flow (z = 0)')
   plt.xlabel('x (au)')
   plt.ylabel('y (au)')

   plt.savefig("output_{:s}/CIR/rad_vel_equ_{:s}.png".format(cir_date, cir_date), dpi=200)
   plt.close()

   Z = np.loadtxt("output_{:s}/CIR/rad_vel_mer_{:s}.dat".format(cir_date, cir_date))
   cmap = plt.colormaps.get_cmap('plasma')
   cmap.set_under("white")
   plt.pcolormesh(X, Y, np.transpose(Z), vmin=np.min(Z[Z > 0.0]), vmax=np.max(Z),
                           cmap=cmap, shading='gouraud')
   plt.colorbar(label='$\\hat{r} \\cdot \\hat{u}$')
   plt.title('Fraction of Radial Plasma Flow (y = 0)')
   plt.xlabel('x (au)')
   plt.ylabel('z (au)')

   plt.savefig("output_{:s}/CIR/rad_vel_mer_{:s}.png".format(cir_date, cir_date), dpi=200)
   plt.close()

   Z = np.loadtxt("output_{:s}/CIR/rad_vel_1au_{:s}.dat".format(cir_date, cir_date))
   Z[:,1] = Z[:,1]
   plt.plot(Z[:,0], Z[:,1])
   plt.title('Fraction of Radial Plasma Flow')
   plt.xlabel('days')
   plt.ylabel('$\\hat{r} \\cdot \\hat{u}$')
   plt.xlim(0, 27.0)

   data1DX.append(Z[:,0])
   data1DY.append(Z[:,1])
   data1DL.append('$\\hat{r} \\cdot \\hat{u}$')

   plt.savefig("output_{:s}/CIR/rad_vel_1au_{:s}.png".format(cir_date, cir_date), dpi=200)
   plt.close()

# ==================================================
if plot_div_vel:
   quants += 1
   print("Plotting divergence of plasma flow")

   Z = np.loadtxt("output_{:s}/CIR/div_vel_equ_{:s}.dat".format(cir_date, cir_date))
   cmap = plt.colormaps.get_cmap('plasma')
   cmap.set_under("white")
   plt.pcolormesh(X, Y, np.transpose(Z), vmin=np.min(Z), vmax=np.max(Z),
                           cmap=cmap, shading='gouraud')
   plt.colorbar(label='$\\nabla \\cdot u$ (1/s)')
   plt.title('Divergence of Plasma Flow (z = 0)')
   plt.xlabel('x (au)')
   plt.ylabel('y (au)')

   plt.savefig("output_{:s}/CIR/div_vel_equ_{:s}.png".format(cir_date, cir_date), dpi=200)
   plt.close()

   Z = np.loadtxt("output_{:s}/CIR/div_vel_mer_{:s}.dat".format(cir_date, cir_date))
   cmap = plt.colormaps.get_cmap('plasma')
   cmap.set_under("white")
   plt.pcolormesh(X, Y, np.transpose(Z), vmin=np.min(Z), vmax=np.max(Z),
                           cmap=cmap, shading='gouraud')
   plt.colorbar(label='$\\nabla \\cdot u$ (1/s)')
   plt.title('Divergence of Plasma Flow (y = 0)')
   plt.xlabel('x (au)')
   plt.ylabel('z (au)')

   plt.savefig("output_{:s}/CIR/div_vel_mer_{:s}.png".format(cir_date, cir_date), dpi=200)
   plt.close()

   Z = np.loadtxt("output_{:s}/CIR/div_vel_1au_{:s}.dat".format(cir_date, cir_date))
   Z[:,1] = Z[:,1]
   plt.plot(Z[:,0], Z[:,1])
   plt.title('Divergence of Plasma Flow (1/s)')
   plt.xlabel('days')
   plt.ylabel('$\\nabla \\cdot u$ (km/s)')
   plt.xlim(0, 27.0)

   data1DX.append(Z[:,0])
   data1DY.append(Z[:,1])
   data1DL.append('$\\nabla \\cdot u$ (1/s)')

   plt.savefig("output_{:s}/CIR/div_vel_1au_{:s}.png".format(cir_date, cir_date), dpi=200)
   plt.close()

# ==================================================
if plot_mag:
   quants += 1
   print("Plotting magnetic field")

   Z = np.loadtxt("output_{:s}/CIR/mag_equ_{:s}.dat".format(cir_date, cir_date)) / 1.0e-5 # G --> nT
   plt.pcolormesh(X, Y, np.transpose(Z), norm=LogNorm(vmin=0.1, vmax=np.max(Z)),
                           cmap='viridis', shading='gouraud')
   plt.colorbar(label='|B| (nT)')
   plt.title('Magnetic Field Magnitude (z = 0)')
   plt.xlabel('x (au)')
   plt.ylabel('y (au)')

# Plot field-lines if enabled
   if plot_Blines:
      for i in range(1,11):
         filename = "output_{:s}/CIR/main_test_fieldline_{:d}_{:s}.lines".format(cir_date, i, cir_date)
         F = np.loadtxt(filename,delimiter=",")
         m = midpoint(F)
         plt.plot(F[:,0], F[:,1], linewidth=1, color="white")
         plt.arrow(F[m,0], F[m,1], F[m+1,0]-F[m,0], F[m+1,1]-F[m,1],
                   head_width=0.1, color="white")
      plt.xlim(np.min(x),np.max(x))
      plt.ylim(np.min(y),np.max(y))

   plt.savefig("output_{:s}/CIR/mag_equ_{:s}.png".format(cir_date, cir_date), dpi=200)
   plt.close()

   Z = np.loadtxt("output_{:s}/CIR/mag_mer_{:s}.dat".format(cir_date, cir_date)) / 1.0e-5 # G --> nT
   plt.pcolormesh(X, Y, np.transpose(Z), norm=LogNorm(vmin=0.1, vmax=np.max(Z)),
                           cmap='viridis', shading='gouraud')
   plt.colorbar(label='|B| (nT)')
   plt.title('Magnetic Field Magnitude (y = 0)')
   plt.xlabel('x (au)')
   plt.ylabel('y (au)')

   plt.savefig("output_{:s}/CIR/mag_mer_{:s}.png".format(cir_date, cir_date), dpi=200)
   plt.close()

   Z = np.loadtxt("output_{:s}/CIR/mag_1au_{:s}.dat".format(cir_date, cir_date))
   Z[:,1] = Z[:,1] / 1.0e-5 # G --> nT
   plt.plot(Z[:,0], Z[:,1])
   plt.title('Magnetic Field Magnitude (1 au)')
   plt.xlabel('days')
   plt.ylabel('|B| (nT)')
   plt.xlim(0, 27.0)

   data1DX.append(Z[:,0])
   data1DY.append(Z[:,1])
   data1DL.append('|B| (nT)')

   plt.savefig("output_{:s}/CIR/mag_1au_{:s}.png".format(cir_date, cir_date), dpi=200)
   plt.close()

# ==================================================
if plot_pol:
   quants += 1
   print("Plotting magnetic polarity")

   Z = np.loadtxt("output_{:s}/CIR/pol_equ_{:s}.dat".format(cir_date, cir_date))
   plt.pcolormesh(X, Y, np.transpose(Z), cmap='bwr', shading='gouraud')
   plt.colorbar(label='sgn(B) (+out/-in)')
   plt.title('Magnetic Field Polarity (z = 0)')
   plt.xlabel('x (au)')
   plt.ylabel('y (au)')

   plt.savefig("output_{:s}/CIR/pol_equ_{:s}.png".format(cir_date, cir_date), dpi=200)
   plt.close()

   Z = np.loadtxt("output_{:s}/CIR/pol_mer_{:s}.dat".format(cir_date, cir_date))
   plt.pcolormesh(X, Y, np.transpose(Z), cmap='bwr', shading='gouraud')
   plt.colorbar(label='sgn(B) (+out/-in)')
   plt.title('Magnetic Field Polarity (y = 0)')
   plt.xlabel('x (au)')
   plt.ylabel('y (au)')

   plt.savefig("output_{:s}/CIR/pol_mer_{:s}.png".format(cir_date, cir_date), dpi=200)
   plt.close()

   Z = np.loadtxt("output_{:s}/CIR/pol_1au_{:s}.dat".format(cir_date, cir_date))
   plt.plot(Z[:,0], Z[:,1])
   plt.title('Magnetic Field Polarity (1 au)')
   plt.xlabel('days')
   plt.ylabel('sgn(B) (+out/-in)')
   plt.xlim(0, 27.0)

   data1DX.append(Z[:,0])
   data1DY.append(Z[:,1])
   data1DL.append('sgn(B) (+out/-in)')

   plt.savefig("output_{:s}/CIR/pol_1au_{:s}.png".format(cir_date, cir_date), dpi=200)
   plt.close()

# ==================================================
if plot_dmax:
   quants += 1
   print("Plotting cell size")

   Z = np.loadtxt("output_{:s}/CIR/dmax_equ_{:s}.dat".format(cir_date, cir_date))
   plt.pcolormesh(X, Y, np.transpose(Z), norm=LogNorm(vmin=np.min(Z[Z > 0.0]), vmax=np.max(Z)),
                           cmap='inferno', shading='gouraud')
   plt.colorbar(label='$\\Delta x$ (au)')
   plt.title('Cell Size (z = 0)')
   plt.xlabel('x (au)')
   plt.ylabel('y (au)')

   plt.savefig("output_{:s}/CIR/dmax_equ_{:s}.png".format(cir_date, cir_date), dpi=200)
   plt.close()

   Z = np.loadtxt("output_{:s}/CIR/dmax_mer_{:s}.dat".format(cir_date, cir_date))
   plt.pcolormesh(X, Y, np.transpose(Z), norm=LogNorm(vmin=np.min(Z[Z > 0.0]), vmax=np.max(Z)),
                           cmap='inferno', shading='gouraud')
   plt.colorbar(label='$\\Delta x$ (au)')
   plt.title('Cell Size (y = 0)')
   plt.xlabel('x (au)')
   plt.ylabel('y (au)')

   plt.savefig("output_{:s}/CIR/dmax_mer_{:s}.png".format(cir_date, cir_date), dpi=200)
   plt.close()

   Z = np.loadtxt("output_{:s}/CIR/dmax_1au_{:s}.dat".format(cir_date, cir_date))
   plt.plot(Z[:,0], Z[:,1])
   plt.title('Cell Size (1 au)')
   plt.xlabel('days')
   plt.ylabel('$\\Delta x$ (au)')
   plt.xlim(0, 27.0)

   data1DX.append(Z[:,0])
   data1DY.append(Z[:,1])
   data1DL.append('$\\Delta x$ (au)')

   plt.savefig("output_{:s}/CIR/dmax_1au_{:s}.png".format(cir_date, cir_date), dpi=200)
   plt.close()

# ==================================================
if plot_het_flx:
   quants += 1
   print("Plotting heat flux speed")

   Z = np.loadtxt("output_{:s}/CIR/het_flx_equ_{:s}.dat".format(cir_date, cir_date))
   cmap = plt.colormaps.get_cmap('magma')
   cmap.set_under("white")
   plt.pcolormesh(X, Y, np.transpose(Z), vmin=np.min(Z[Z > 0.0]), vmax=np.max(Z),
                           cmap=cmap, shading='gouraud')
   plt.colorbar(label='$Q_e$ (erg)')
   plt.title('Heat Flux (z = 0)')
   plt.xlabel('x (au)')
   plt.ylabel('y (au)')

   plt.savefig("output_{:s}/CIR/het_flx_equ_{:s}.png".format(cir_date, cir_date), dpi=200)
   plt.close()

   Z = np.loadtxt("output_{:s}/CIR/het_flx_mer_{:s}.dat".format(cir_date, cir_date))
   cmap = plt.colormaps.get_cmap('magma')
   cmap.set_under("white")
   plt.pcolormesh(X, Y, np.transpose(Z), vmin=np.min(Z[Z > 0.0]), vmax=np.max(Z),
                           cmap=cmap, shading='gouraud')
   plt.colorbar(label='$Q_e$ (erg)')
   plt.title('Heat Flux (y = 0)')
   plt.xlabel('x (au)')
   plt.ylabel('y (au)')

   plt.savefig("output_{:s}/CIR/het_flx_mer_{:s}.png".format(cir_date, cir_date), dpi=200)
   plt.close()

   Z = np.loadtxt("output_{:s}/CIR/het_flx_1au_{:s}.dat".format(cir_date, cir_date))
   plt.plot(Z[:,0], Z[:,1])
   plt.title('Heat Flux (1 au)')
   plt.xlabel('days')
   plt.ylabel('$Q_e$ (erg)')
   plt.xlim(0, 27.0)

   data1DX.append(Z[:,0])
   data1DY.append(Z[:,1])
   data1DL.append('$Q_e (erg)$')

   plt.savefig("output_{:s}/CIR/het_flx_1au_{:s}.png".format(cir_date, cir_date), dpi=200)
   plt.close()

# ==================================================
if plot_tur_enr:
   quants += 1
   print("Plotting turbulent energy")

   Z = np.loadtxt("output_{:s}/CIR/tur_enr_equ_{:s}.dat".format(cir_date, cir_date))
   cmap = plt.colormaps.get_cmap('cividis')
   cmap.set_under("white")
   plt.pcolormesh(X, Y, np.transpose(Z), vmin=np.min(Z[Z > 0.0]), vmax=np.max(Z),
                           cmap=cmap, shading='gouraud')
   plt.colorbar(label='$\\langle \\delta B^2 \\rangle$ (erg)')
   plt.title('Turbulent Energy (z = 0)')
   plt.xlabel('x (au)')
   plt.ylabel('y (au)')

   plt.savefig("output_{:s}/CIR/tur_enr_equ_{:s}.png".format(cir_date, cir_date), dpi=200)
   plt.close()

   Z = np.loadtxt("output_{:s}/CIR/tur_enr_mer_{:s}.dat".format(cir_date, cir_date))
   cmap = plt.colormaps.get_cmap('cividis')
   cmap.set_under("white")
   plt.pcolormesh(X, Y, np.transpose(Z), vmin=np.min(Z[Z > 0.0]), vmax=np.max(Z),
                           cmap=cmap, shading='gouraud')
   plt.colorbar(label='$\\langle \\delta B^2 \\rangle$ (erg)')
   plt.title('Turbulent Energy (y = 0)')
   plt.xlabel('x (au)')
   plt.ylabel('y (au)')

   plt.savefig("output_{:s}/CIR/tur_enr_mer_{:s}.png".format(cir_date, cir_date), dpi=200)
   plt.close()

   Z = np.loadtxt("output_{:s}/CIR/tur_enr_1au_{:s}.dat".format(cir_date, cir_date))
   plt.plot(Z[:,0], Z[:,1])
   plt.title('Turbulent Energy (1 au)')
   plt.xlabel('days')
   plt.ylabel('$\\langle \\delta B^2 \\rangle$ (erg)')
   plt.xlim(0, 27.0)

   data1DX.append(Z[:,0])
   data1DY.append(Z[:,1])
   data1DL.append('$\\langle \\delta B^2 \\rangle$ (erg)')

   plt.savefig("output_{:s}/CIR/tur_enr_1au_{:s}.png".format(cir_date, cir_date), dpi=200)
   plt.close()

# ==================================================
if plot_drift:
   quants += 1
   print("Plotting drift speed")

   Z = np.loadtxt("output_{:s}/CIR/drift_equ_{:s}.dat".format(cir_date, cir_date))
   cmap = plt.colormaps.get_cmap('plasma')
   cmap.set_under("white")
   plt.pcolormesh(X, Y, np.transpose(Z), vmin=np.min(Z[Z > 0.0]), vmax=np.max(Z),
                           cmap=cmap, shading='gouraud')
   plt.colorbar(label='$v_d$ (/ $v_p$)')
   plt.title('Magnetic Drift Speed (z = 0)')
   plt.xlabel('x (au)')
   plt.ylabel('y (au)')

   plt.savefig("output_{:s}/CIR/drift_equ_{:s}.png".format(cir_date, cir_date), dpi=200)
   plt.close()

   Z = np.loadtxt("output_{:s}/CIR/drift_mer_{:s}.dat".format(cir_date, cir_date))
   cmap = plt.colormaps.get_cmap('plasma')
   cmap.set_under("white")
   plt.pcolormesh(X, Y, np.transpose(Z), vmin=np.min(Z[Z > 0.0]), vmax=np.max(Z),
                           cmap=cmap, shading='gouraud')
   plt.colorbar(label='$v_d$ (/ $v_p$)')
   plt.title('Magnetic Drift Speed (y = 0)')
   plt.xlabel('x (au)')
   plt.ylabel('y (au)')

   plt.savefig("output_{:s}/CIR/drift_mer_{:s}.png".format(cir_date, cir_date), dpi=200)
   plt.close()

   Z = np.loadtxt("output_{:s}/CIR/drift_1au_{:s}.dat".format(cir_date, cir_date))
   plt.plot(Z[:,0], Z[:,1])
   plt.title('Magnetic Drift Speed (1 au)')
   plt.xlabel('days')
   plt.ylabel('$v_d$ (/ $v_p$)')
   plt.xlim(0, 27.0)

   data1DX.append(Z[:,0])
   data1DY.append(Z[:,1])
   data1DL.append('$v_d$ (/ $v_p$)')

   plt.savefig("output_{:s}/CIR/drift_1au_{:s}.png".format(cir_date, cir_date), dpi=200)
   plt.close()

# ==================================================
if plot_diff1:
   quants += 1
   print("Plotting diffusion 1 (perpendicular)")

   Z = np.loadtxt("output_{:s}/CIR/diff1_equ_{:s}.dat".format(cir_date, cir_date))
   cmap = plt.colormaps.get_cmap('plasma')
   cmap.set_under("white")
   plt.pcolormesh(X, Y, np.transpose(Z), norm=LogNorm(vmin=np.min(Z[Z > 0.0]), vmax=np.max(Z)),
                           cmap=cmap, shading='gouraud')
   plt.colorbar(label='$\\kappa_{\\perp}$')
   plt.title('Diffusion Model Perpendicular (z = 0)')
   plt.xlabel('x (au)')
   plt.ylabel('y (au)')

   plt.savefig("output_{:s}/CIR/diff1_equ_{:s}.png".format(cir_date, cir_date), dpi=200)
   plt.close()

   Z = np.loadtxt("output_{:s}/CIR/diff1_mer_{:s}.dat".format(cir_date, cir_date))
   cmap = plt.colormaps.get_cmap('plasma')
   cmap.set_under("white")
   plt.pcolormesh(X, Y, np.transpose(Z), norm=LogNorm(vmin=np.min(Z[Z > 0.0]), vmax=np.max(Z)),
                           cmap=cmap, shading='gouraud')
   plt.colorbar(label='$\\kappa_{\\perp}$')
   plt.title('Diffusion Model Perpendicular (y = 0)')
   plt.xlabel('x (au)')
   plt.ylabel('y (au)')

   plt.savefig("output_{:s}/CIR/diff1_mer_{:s}.png".format(cir_date, cir_date), dpi=200)
   plt.close()

   Z = np.loadtxt("output_{:s}/CIR/diff1_1au_{:s}.dat".format(cir_date, cir_date))
   plt.plot(Z[:,0], Z[:,1])
   plt.title('Diffusion Model Perpendicular (1 au)')
   plt.xlabel('days')
   plt.ylabel('$\\kappa_{\\perp}$')
   plt.xlim(0, 27.0)

   data1DX.append(Z[:,0])
   data1DY.append(Z[:,1])
   data1DL.append('$\\kappa_{\\perp}$')

   plt.savefig("output_{:s}/CIR/diff1_1au_{:s}.png".format(cir_date, cir_date), dpi=200)
   plt.close()

# ==================================================
if plot_diff2:
   quants += 1
   print("Plotting diffusion 2 (parallel)")

   Z = np.loadtxt("output_{:s}/CIR/diff2_equ_{:s}.dat".format(cir_date, cir_date))
   cmap = plt.colormaps.get_cmap('plasma')
   cmap.set_under("white")
   plt.pcolormesh(X, Y, np.transpose(Z), norm=LogNorm(vmin=np.min(Z[Z > 0.0]), vmax=np.max(Z)),
                           cmap=cmap, shading='gouraud')
   plt.colorbar(label='$\\kappa_{\\parallel}$')
   plt.title('Diffusion Model Parallel (z = 0)')
   plt.xlabel('x (au)')
   plt.ylabel('y (au)')

   plt.savefig("output_{:s}/CIR/diff2_equ_{:s}.png".format(cir_date, cir_date), dpi=200)
   plt.close()

   Z = np.loadtxt("output_{:s}/CIR/diff2_mer_{:s}.dat".format(cir_date, cir_date))
   cmap = plt.colormaps.get_cmap('plasma')
   cmap.set_under("white")
   plt.pcolormesh(X, Y, np.transpose(Z), norm=LogNorm(vmin=np.min(Z[Z > 0.0]), vmax=np.max(Z)),
                           cmap=cmap, shading='gouraud')
   plt.colorbar(label='$\\kappa_{\\parallel}$')
   plt.title('Diffusion Model Parallel (y = 0)')
   plt.xlabel('x (au)')
   plt.ylabel('y (au)')

   plt.savefig("output_{:s}/CIR/diff2_mer_{:s}.png".format(cir_date, cir_date), dpi=200)
   plt.close()

   Z = np.loadtxt("output_{:s}/CIR/diff2_1au_{:s}.dat".format(cir_date, cir_date))
   plt.plot(Z[:,0], Z[:,1])
   plt.title('Diffusion Model Parallel (1 au)')
   plt.xlabel('days')
   plt.ylabel('$\\kappa_{\\parallel}$')
   plt.xlim(0, 27.0)

   data1DX.append(Z[:,0])
   data1DY.append(Z[:,1])
   data1DL.append('$\\kappa_{\\parallel}$')

   plt.savefig("output_{:s}/CIR/diff2_1au_{:s}.png".format(cir_date, cir_date), dpi=200)
   plt.close()

# ==================================================
print("Plotting", quants, "1D quantities in composite figure")
if quants > 0:
   fig, axs = plt.subplots(nrows=quants, ncols=1, figsize=(10, 3*quants))

   for i in range(quants):
      axs[i].plot(data1DX[i], data1DY[i])
      axs[i].set_ylabel(data1DL[i])

   axs[quants-1].set_xlabel('days')
   plt.xlim(0, 27.0)
   plt.tight_layout()
   plt.savefig("output_{:s}/CIR/composite_1au_{:s}.png".format(cir_date, cir_date), dpi=200)
   plt.show()