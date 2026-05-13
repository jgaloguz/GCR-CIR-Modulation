# Code to plot 2D slices of solar wind and particle quantities

# Import libraries
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm
import argparse

# Parse arguments
parser = argparse.ArgumentParser(description="Plot 2D slices of solar wind and particle quantities")
parser.add_argument("date", type=str, help="Date of CIR in YYYYMMDD format")
parser.add_argument("num", type=int, help="Number of data points along each dimension")
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
x = np.linspace(-5.0, 5.0, args.num)
y = np.linspace(-5.0, 5.0, args.num)
X, Y = np.meshgrid(x, y)

print("====================")
print("Plotting 2D slices:")
print("====================")

# Create the pcolormesh plots
# ==================================================
if plot_den:
   print("Plotting plasma density")

   Z = np.loadtxt("output_{:s}/CIR/den_equ_{:s}.dat".format(args.date, args.date))
   cmap = plt.colormaps.get_cmap('magma')
   cmap.set_under("white")
   plt.pcolormesh(X, Y, np.transpose(Z), norm=LogNorm(vmin=np.min(Z[Z > 0.0]), vmax=np.max(Z)),
                           cmap=cmap, shading='gouraud')
   plt.colorbar(label='n (amu/cc)')
   plt.title('Plasma Density (z = 0)')
   plt.xlabel('x (au)')
   plt.ylabel('y (au)')

   plt.savefig("output_{:s}/CIR/den_equ_{:s}.png".format(args.date, args.date), dpi=200)
   plt.close()

   Z = np.loadtxt("output_{:s}/CIR/den_mer_{:s}.dat".format(args.date, args.date))
   cmap = plt.colormaps.get_cmap('magma')
   cmap.set_under("white")
   plt.pcolormesh(X, Y, np.transpose(Z), norm=LogNorm(vmin=np.min(Z[Z > 0.0]), vmax=np.max(Z)),
                           cmap=cmap, shading='gouraud')
   plt.colorbar(label='n (amu/cc)')
   plt.title('Plasma Density (y = 0)')
   plt.xlabel('x (au)')
   plt.ylabel('z (au)')

   plt.savefig("output_{:s}/CIR/den_mer_{:s}.png".format(args.date, args.date), dpi=200)
   plt.close()

# ==================================================
if plot_vel:
   print("Plotting plasma flow")

   Z = np.loadtxt("output_{:s}/CIR/vel_equ_{:s}.dat".format(args.date, args.date)) / 1.0e5 # cm/s --> km/s
   cmap = plt.colormaps.get_cmap('plasma')
   cmap.set_under("white")
   plt.pcolormesh(X, Y, np.transpose(Z), vmin=np.min(Z[Z > 0.0]), vmax=np.max(Z),
                           cmap=cmap, shading='gouraud')
   plt.colorbar(label='|u| (km/s)')
   plt.title('Plasma Flow Magnitude (z = 0)')
   plt.xlabel('x (au)')
   plt.ylabel('y (au)')

   plt.savefig("output_{:s}/CIR/vel_equ_{:s}.png".format(args.date, args.date), dpi=200)
   plt.close()

   Z = np.loadtxt("output_{:s}/CIR/vel_mer_{:s}.dat".format(args.date, args.date)) / 1.0e5 # cm/s --> km/s
   cmap = plt.colormaps.get_cmap('plasma')
   cmap.set_under("white")
   plt.pcolormesh(X, Y, np.transpose(Z), vmin=np.min(Z[Z > 0.0]), vmax=np.max(Z),
                           cmap=cmap, shading='gouraud')
   plt.colorbar(label='|u| (km/s)')
   plt.title('Plasma Flow Magnitude (y = 0)')
   plt.xlabel('x (au)')
   plt.ylabel('z (au)')

   plt.savefig("output_{:s}/CIR/vel_mer_{:s}.png".format(args.date, args.date), dpi=200)
   plt.close()

# ==================================================
if plot_mag:
   print("Plotting magnetic field")

   Z = np.loadtxt("output_{:s}/CIR/mag_equ_{:s}.dat".format(args.date, args.date)) / 1.0e-5 # G --> nT
   plt.pcolormesh(X, Y, np.transpose(Z), norm=LogNorm(vmin=0.1, vmax=np.max(Z)),
                           cmap='viridis', shading='gouraud')
   plt.colorbar(label='|B| (nT)')
   plt.title('Magnetic Field Magnitude (z = 0)')
   plt.xlabel('x (au)')
   plt.ylabel('y (au)')

   plt.savefig("output_{:s}/CIR/mag_equ_{:s}.png".format(args.date, args.date), dpi=200)
   plt.close()

   Z = np.loadtxt("output_{:s}/CIR/mag_mer_{:s}.dat".format(args.date, args.date)) / 1.0e-5 # G --> nT
   plt.pcolormesh(X, Y, np.transpose(Z), norm=LogNorm(vmin=0.1, vmax=np.max(Z)),
                           cmap='viridis', shading='gouraud')
   plt.colorbar(label='|B| (nT)')
   plt.title('Magnetic Field Magnitude (y = 0)')
   plt.xlabel('x (au)')
   plt.ylabel('z (au)')

   plt.savefig("output_{:s}/CIR/mag_mer_{:s}.png".format(args.date, args.date), dpi=200)
   plt.close()

# ==================================================
if plot_tmp:
   print("Plotting plasma temperature")

   P = np.loadtxt("output_{:s}/CIR/pth_equ_{:s}.dat".format(args.date, args.date)) * 0.1 # dyn/cm^.2 --> Pa
   n = np.loadtxt("output_{:s}/CIR/den_equ_{:s}.dat".format(args.date, args.date)) * 1.0e6 # cm^-3 --> m^-3
   Z = np.zeros_like(P)
   np.divide(P, n * 1.380649e-23, out=Z, where=n != 0) # T = P / (n*k_B)
   Z = Z / 1.0e6 # K --> MK
   cmap = plt.colormaps.get_cmap('magma')
   cmap.set_under("white")
   plt.pcolormesh(X, Y, np.transpose(Z), norm=LogNorm(vmin=np.min(Z[Z > 0.0]), vmax=np.max(Z)),
                           cmap=cmap, shading='gouraud')
   plt.colorbar(label='T (MK)')
   plt.title('Plasma Temperature (z = 0)')
   plt.xlabel('x (au)')
   plt.ylabel('y (au)')

   plt.savefig("output_{:s}/CIR/tmp_equ_{:s}.png".format(args.date, args.date), dpi=200)
   plt.close()

   P = np.loadtxt("output_{:s}/CIR/pth_mer_{:s}.dat".format(args.date, args.date)) * 0.1 # dyn/cm^.2 --> Pa
   n = np.loadtxt("output_{:s}/CIR/den_mer_{:s}.dat".format(args.date, args.date)) * 1.0e6 # cm^-3 --> m^-3
   Z = np.zeros_like(P)
   np.divide(P, n * 1.380649e-23, out=Z, where=n != 0) # T = P / (n*k_B)
   Z = Z / 1.0e6 # K --> MK
   cmap = plt.colormaps.get_cmap('magma')
   cmap.set_under("white")
   plt.pcolormesh(X, Y, np.transpose(Z), norm=LogNorm(vmin=np.min(Z[Z > 0.0]), vmax=np.max(Z)),
                           cmap=cmap, shading='gouraud')
   plt.colorbar(label='T (MK)')
   plt.title('Plasma Temperature (y = 0)')
   plt.xlabel('x (au)')
   plt.ylabel('z (au)')

   plt.savefig("output_{:s}/CIR/tmp_mer_{:s}.png".format(args.date, args.date), dpi=200)
   plt.close()

# ==================================================
if plot_pol:
   print("Plotting magnetic polarity")

   Z = np.loadtxt("output_{:s}/CIR/pol_equ_{:s}.dat".format(args.date, args.date))
   plt.pcolormesh(X, Y, np.transpose(Z), cmap='bwr', shading='gouraud')
   plt.colorbar(label='sgn(B) (+out/-in)')
   plt.title('Magnetic Field Polarity (z = 0)')
   plt.xlabel('x (au)')
   plt.ylabel('y (au)')

   plt.savefig("output_{:s}/CIR/pol_equ_{:s}.png".format(args.date, args.date), dpi=200)
   plt.close()

   Z = np.loadtxt("output_{:s}/CIR/pol_mer_{:s}.dat".format(args.date, args.date))
   plt.pcolormesh(X, Y, np.transpose(Z), cmap='bwr', shading='gouraud')
   plt.colorbar(label='sgn(B) (+out/-in)')
   plt.title('Magnetic Field Polarity (y = 0)')
   plt.xlabel('x (au)')
   plt.ylabel('z (au)')

   plt.savefig("output_{:s}/CIR/pol_mer_{:s}.png".format(args.date, args.date), dpi=200)
   plt.close()

# ==================================================
if plot_dmax:
   print("Plotting cell size")

   Z = np.loadtxt("output_{:s}/CIR/dmax_equ_{:s}.dat".format(args.date, args.date))
   plt.pcolormesh(X, Y, np.transpose(Z), norm=LogNorm(vmin=np.min(Z[Z > 0.0]), vmax=np.max(Z)),
                           cmap='inferno', shading='gouraud')
   plt.colorbar(label='$\\Delta x$ (au)')
   plt.title('Cell Size (z = 0)')
   plt.xlabel('x (au)')
   plt.ylabel('y (au)')

   plt.savefig("output_{:s}/CIR/dmax_equ_{:s}.png".format(args.date, args.date), dpi=200)
   plt.close()

   Z = np.loadtxt("output_{:s}/CIR/dmax_mer_{:s}.dat".format(args.date, args.date))
   plt.pcolormesh(X, Y, np.transpose(Z), norm=LogNorm(vmin=np.min(Z[Z > 0.0]), vmax=np.max(Z)),
                           cmap='inferno', shading='gouraud')
   plt.colorbar(label='$\\Delta x$ (au)')
   plt.title('Cell Size (y = 0)')
   plt.xlabel('x (au)')
   plt.ylabel('z (au)')

   plt.savefig("output_{:s}/CIR/dmax_mer_{:s}.png".format(args.date, args.date), dpi=200)
   plt.close()

# ==================================================
if plot_tur_enr:
   print("Plotting turbulent energy")

   Z = np.loadtxt("output_{:s}/CIR/tur_enr_equ_{:s}.dat".format(args.date, args.date))
   cmap = plt.colormaps.get_cmap('cividis')
   cmap.set_under("white")
   plt.pcolormesh(X, Y, np.transpose(Z), norm=LogNorm(vmin=np.min(Z[Z > 0.0]), vmax=np.max(Z)),
                           cmap=cmap, shading='gouraud')
   plt.colorbar(label='$Z^2$ (cm$^2$ s$^{-2}$)')
   plt.title('Turbulent Energy (z = 0)')
   plt.xlabel('x (au)')
   plt.ylabel('y (au)')

   plt.savefig("output_{:s}/CIR/tur_enr_equ_{:s}.png".format(args.date, args.date), dpi=200)
   plt.close()

   Z = np.loadtxt("output_{:s}/CIR/tur_enr_mer_{:s}.dat".format(args.date, args.date))
   cmap = plt.colormaps.get_cmap('cividis')
   cmap.set_under("white")
   plt.pcolormesh(X, Y, np.transpose(Z), norm=LogNorm(vmin=np.min(Z[Z > 0.0]), vmax=np.max(Z)),
                           cmap=cmap, shading='gouraud')
   plt.colorbar(label='$Z^2$ (cm$^2$ s$^{-2}$)')
   plt.title('Turbulent Energy (y = 0)')
   plt.xlabel('x (au)')
   plt.ylabel('z (au)')

   plt.savefig("output_{:s}/CIR/tur_enr_mer_{:s}.png".format(args.date, args.date), dpi=200)
   plt.close()

# ==================================================
if plot_het_flx:
   print("Plotting heat flux speed")

   Z = np.loadtxt("output_{:s}/CIR/het_flx_equ_{:s}.dat".format(args.date, args.date))
   cmap = plt.colormaps.get_cmap('magma')
   cmap.set_under("white")
   plt.pcolormesh(X, Y, np.transpose(Z), norm=LogNorm(vmin=np.min(Z[Z > 0.0]), vmax=np.max(Z)),
                           cmap=cmap, shading='gouraud')
   plt.colorbar(label='$Q_e$ (cm$^2$ s$^{-2}$)')
   plt.title('Heat Flux (z = 0)')
   plt.xlabel('x (au)')
   plt.ylabel('y (au)')

   plt.savefig("output_{:s}/CIR/het_flx_equ_{:s}.png".format(args.date, args.date), dpi=200)
   plt.close()

   Z = np.loadtxt("output_{:s}/CIR/het_flx_mer_{:s}.dat".format(args.date, args.date))
   cmap = plt.colormaps.get_cmap('magma')
   cmap.set_under("white")
   plt.pcolormesh(X, Y, np.transpose(Z), norm=LogNorm(vmin=np.min(Z[Z > 0.0]), vmax=np.max(Z)),
                           cmap=cmap, shading='gouraud')
   plt.colorbar(label='$Q_e$ (cm$^2$ s$^{-2}$)')
   plt.title('Heat Flux (y = 0)')
   plt.xlabel('x (au)')
   plt.ylabel('z (au)')

   plt.savefig("output_{:s}/CIR/het_flx_mer_{:s}.png".format(args.date, args.date), dpi=200)
   plt.close()

# ==================================================
if plot_rad_vel:
   print("Plotting fraction of radial plasma flow")

   Z = np.loadtxt("output_{:s}/CIR/rad_vel_equ_{:s}.dat".format(args.date, args.date))
   cmap = plt.colormaps.get_cmap('plasma')
   cmap.set_under("white")
   plt.pcolormesh(X, Y, np.transpose(Z), vmin=np.min(Z[Z > 0.0]), vmax=np.max(Z),
                           cmap=cmap, shading='gouraud')
   plt.colorbar(label='$\\hat{r} \\cdot \\hat{u}$')
   plt.title('Fraction of Radial Plasma Flow (z = 0)')
   plt.xlabel('x (au)')
   plt.ylabel('y (au)')

   plt.savefig("output_{:s}/CIR/rad_vel_equ_{:s}.png".format(args.date, args.date), dpi=200)
   plt.close()

   Z = np.loadtxt("output_{:s}/CIR/rad_vel_mer_{:s}.dat".format(args.date, args.date))
   cmap = plt.colormaps.get_cmap('plasma')
   cmap.set_under("white")
   plt.pcolormesh(X, Y, np.transpose(Z), vmin=np.min(Z[Z > 0.0]), vmax=np.max(Z),
                           cmap=cmap, shading='gouraud')
   plt.colorbar(label='$\\hat{r} \\cdot \\hat{u}$')
   plt.title('Fraction of Radial Plasma Flow (y = 0)')
   plt.xlabel('x (au)')
   plt.ylabel('z (au)')

   plt.savefig("output_{:s}/CIR/rad_vel_mer_{:s}.png".format(args.date, args.date), dpi=200)
   plt.close()

# ==================================================
if plot_div_vel:
   print("Plotting divergence of plasma flow")

   Z = np.loadtxt("output_{:s}/CIR/div_vel_equ_{:s}.dat".format(args.date, args.date))
   cmap = plt.colormaps.get_cmap('plasma')
   cmap.set_under("white")
   plt.pcolormesh(X, Y, np.transpose(Z), vmin=np.min(Z), vmax=np.max(Z),
                           cmap=cmap, shading='gouraud')
   plt.colorbar(label='$\\nabla \\cdot u$ (1/s)')
   plt.title('Divergence of Plasma Flow (z = 0)')
   plt.xlabel('x (au)')
   plt.ylabel('y (au)')

   plt.savefig("output_{:s}/CIR/div_vel_equ_{:s}.png".format(args.date, args.date), dpi=200)
   plt.close()

   Z = np.loadtxt("output_{:s}/CIR/div_vel_mer_{:s}.dat".format(args.date, args.date))
   cmap = plt.colormaps.get_cmap('plasma')
   cmap.set_under("white")
   plt.pcolormesh(X, Y, np.transpose(Z), vmin=np.min(Z), vmax=np.max(Z),
                           cmap=cmap, shading='gouraud')
   plt.colorbar(label='$\\nabla \\cdot u$ (1/s)')
   plt.title('Divergence of Plasma Flow (y = 0)')
   plt.xlabel('x (au)')
   plt.ylabel('z (au)')

   plt.savefig("output_{:s}/CIR/div_vel_mer_{:s}.png".format(args.date, args.date), dpi=200)
   plt.close()

# ==================================================
if plot_drift:
   print("Plotting drift speed")

   Z = np.loadtxt("output_{:s}/CIR/drift_equ_{:s}.dat".format(args.date, args.date))
   cmap = plt.colormaps.get_cmap('plasma')
   cmap.set_under("white")
   plt.pcolormesh(X, Y, np.transpose(Z), vmin=np.min(Z[Z > 0.0]), vmax=np.max(Z),
                           cmap=cmap, shading='gouraud')
   plt.colorbar(label='$v_d$ (/ $v_p$)')
   plt.title('Magnetic Drift Speed (z = 0)')
   plt.xlabel('x (au)')
   plt.ylabel('y (au)')

   plt.savefig("output_{:s}/CIR/drift_equ_{:s}.png".format(args.date, args.date), dpi=200)
   plt.close()

   Z = np.loadtxt("output_{:s}/CIR/drift_mer_{:s}.dat".format(args.date, args.date))
   cmap = plt.colormaps.get_cmap('plasma')
   cmap.set_under("white")
   plt.pcolormesh(X, Y, np.transpose(Z), vmin=np.min(Z[Z > 0.0]), vmax=np.max(Z),
                           cmap=cmap, shading='gouraud')
   plt.colorbar(label='$v_d$ (/ $v_p$)')
   plt.title('Magnetic Drift Speed (y = 0)')
   plt.xlabel('x (au)')
   plt.ylabel('z (au)')

   plt.savefig("output_{:s}/CIR/drift_mer_{:s}.png".format(args.date, args.date), dpi=200)
   plt.close()

# ==================================================
if plot_diff1:
   print("Plotting diffusion 1 (perpendicular)")

   Z = np.loadtxt("output_{:s}/CIR/diff1_equ_{:s}.dat".format(args.date, args.date))
   cmap = plt.colormaps.get_cmap('plasma')
   cmap.set_under("white")
   plt.pcolormesh(X, Y, np.transpose(Z), norm=LogNorm(vmin=np.min(Z[Z > 0.0]), vmax=np.max(Z)),
                           cmap=cmap, shading='gouraud')
   plt.colorbar(label='$\\kappa_{\\perp}$')
   plt.title('Diffusion Model Perpendicular (z = 0)')
   plt.xlabel('x (au)')
   plt.ylabel('y (au)')

   plt.savefig("output_{:s}/CIR/diff1_equ_{:s}.png".format(args.date, args.date), dpi=200)
   plt.close()

   Z = np.loadtxt("output_{:s}/CIR/diff1_mer_{:s}.dat".format(args.date, args.date))
   cmap = plt.colormaps.get_cmap('plasma')
   cmap.set_under("white")
   plt.pcolormesh(X, Y, np.transpose(Z), norm=LogNorm(vmin=np.min(Z[Z > 0.0]), vmax=np.max(Z)),
                           cmap=cmap, shading='gouraud')
   plt.colorbar(label='$\\kappa_{\\perp}$')
   plt.title('Diffusion Model Perpendicular (y = 0)')
   plt.xlabel('x (au)')
   plt.ylabel('z (au)')

   plt.savefig("output_{:s}/CIR/diff1_mer_{:s}.png".format(args.date, args.date), dpi=200)
   plt.close()

# ==================================================
if plot_diff2:
   print("Plotting diffusion 2 (parallel)")

   Z = np.loadtxt("output_{:s}/CIR/diff2_equ_{:s}.dat".format(args.date, args.date))
   cmap = plt.colormaps.get_cmap('plasma')
   cmap.set_under("white")
   plt.pcolormesh(X, Y, np.transpose(Z), norm=LogNorm(vmin=np.min(Z[Z > 0.0]), vmax=np.max(Z)),
                           cmap=cmap, shading='gouraud')
   plt.colorbar(label='$\\kappa_{\\parallel}$')
   plt.title('Diffusion Model Parallel (z = 0)')
   plt.xlabel('x (au)')
   plt.ylabel('y (au)')

   plt.savefig("output_{:s}/CIR/diff2_equ_{:s}.png".format(args.date, args.date), dpi=200)
   plt.close()

   Z = np.loadtxt("output_{:s}/CIR/diff2_mer_{:s}.dat".format(args.date, args.date))
   cmap = plt.colormaps.get_cmap('plasma')
   cmap.set_under("white")
   plt.pcolormesh(X, Y, np.transpose(Z), norm=LogNorm(vmin=np.min(Z[Z > 0.0]), vmax=np.max(Z)),
                           cmap=cmap, shading='gouraud')
   plt.colorbar(label='$\\kappa_{\\parallel}$')
   plt.title('Diffusion Model Parallel (y = 0)')
   plt.xlabel('x (au)')
   plt.ylabel('z (au)')

   plt.savefig("output_{:s}/CIR/diff2_mer_{:s}.png".format(args.date, args.date), dpi=200)
   plt.close()

# ==================================================
if plot_diff3:
   print("Plotting diffusion 3 (radial)")

   Z = np.loadtxt("output_{:s}/CIR/diff3_equ_{:s}.dat".format(args.date, args.date))
   cmap = plt.colormaps.get_cmap('plasma')
   cmap.set_under("white")
   plt.pcolormesh(X, Y, np.transpose(Z), norm=LogNorm(vmin=np.min(Z[Z > 0.0]), vmax=np.max(Z)),
                           cmap=cmap, shading='gouraud')
   plt.colorbar(label='$\\kappa_{rr}$')
   plt.title('Diffusion Model Radial (z = 0)')
   plt.xlabel('x (au)')
   plt.ylabel('y (au)')

   plt.savefig("output_{:s}/CIR/diff3_equ_{:s}.png".format(args.date, args.date), dpi=200)
   plt.close()

   Z = np.loadtxt("output_{:s}/CIR/diff3_mer_{:s}.dat".format(args.date, args.date))
   cmap = plt.colormaps.get_cmap('plasma')
   cmap.set_under("white")
   plt.pcolormesh(X, Y, np.transpose(Z), norm=LogNorm(vmin=np.min(Z[Z > 0.0]), vmax=np.max(Z)),
                           cmap=cmap, shading='gouraud')
   plt.colorbar(label='$\\kappa_{rr}$')
   plt.title('Diffusion Model Radial (y = 0)')
   plt.xlabel('x (au)')
   plt.ylabel('z (au)')

   plt.savefig("output_{:s}/CIR/diff3_mer_{:s}.png".format(args.date, args.date), dpi=200)
   plt.close()