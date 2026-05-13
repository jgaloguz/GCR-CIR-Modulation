# Compute Earth position in HGC frame during SWMF simulation

# Import libraries
import numpy as np
from scipy.interpolate import interp1d
import argparse

# Parse arguments
parser = argparse.ArgumentParser(description="Compute Earth position in HGC frame during SWMF simulation")
parser.add_argument("date", type=str, help="Date of CIR in YYYYMMDD format")
parser.add_argument("--num", type=int, default=0, help="Number of data points to compute table of positions (Default: 0)")
parser.add_argument("--time", type=float, default=0.0, help="Time in days relative to CIR observed date at which to compute position in HGC frame (Default 0.0). Possible values are any number in range [-14.0,14.0]")
args = parser.parse_args()

# Grab position data and make interpolator
Rs_in_au = 6.957e+10 / 1.496e+13
sat_file = np.loadtxt("../../SWMF/run_cir_" + args.date + "/IH/IO2/trj_earth_n00005000.sat", skiprows=2)
n_pts = np.size(sat_file, 0)
t = np.linspace(-14.0, 14.0, num=n_pts)
x = sat_file[:,8] * Rs_in_au
y = sat_file[:,9] * Rs_in_au
z = sat_file[:,10] * Rs_in_au
fx = interp1d(t, x, kind='linear')
fy = interp1d(t, y, kind='linear')
fz = interp1d(t, z, kind='linear')

# Check usage mode
if args.num == 0:    # Single point interpolation
   earth_pos = np.array([fx(args.time),fy(args.time),fz(args.time)])
else:                # Table interpolation
   earth_pos = np.zeros((args.num, 3))
   _t = np.linspace(-14.0, 14.0, num=args.num)
   earth_pos[:, 0] = fx(_t)
   earth_pos[:, 1] = fy(_t)
   earth_pos[:, 2] = fz(_t)

np.savetxt("output_" + args.date + "/earth_position_" + args.date + ".dat", earth_pos)
