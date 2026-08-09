# Code to compute the average diffusion coefficients over a rotation

# Import libraries
import numpy as np
import argparse
from scipy.optimize import curve_fit

# Parse arguments
parser = argparse.ArgumentParser(description="Compute the average diffusion coefficients over a rotation")
parser.add_argument("date", type=str, help="Date of CIR in YYYYMMDD format")
args = parser.parse_args()

# Power law function
def line(x, a, b):
	return a * x + b

# Import data
diff1_data = np.loadtxt("output_{:s}/CIR/diff1_1au_{:s}.dat".format(args.date, args.date))
diff1_magn = diff1_data[:,1]
diff2_data = np.loadtxt("output_{:s}/CIR/diff2_1au_{:s}.dat".format(args.date, args.date))
diff2_magn = diff2_data[:,1]
vel_data = np.loadtxt("output_{:s}/CIR/vel_1au_{:s}.dat".format(args.date, args.date))
vel_magn = vel_data[:,1]

# Report averages over entire time series
print("CIR date: {:s}".format(args.date))

print("kappa_perp arithmetic average (1au) = {:.2e} cm^2 s^-1".format(np.mean(diff1_magn)))
print("kappa_para arithmetic average (1au) = {:.2e} cm^2 s^-1".format(np.mean(diff2_magn)))
print("kappa_perp standard deviation (1au) = {:.2e} cm^2 s^-1".format(np.std(diff1_magn)))
print("kappa_para standard deviation (1au) = {:.2e} cm^2 s^-1".format(np.std(diff2_magn)))

# Fit radial and energy cuts with power laws
directions = ["100", "-100", "010", "0-10", "001", "00-1"]
slope_perp = []
slope_para = []
for direction in directions:
	diff1_data = np.loadtxt("output_{:s}/CIR/diff1_{:s}_{:s}.dat".format(args.date, direction, args.date))
	popt, pcov = curve_fit(line, np.log(diff1_data[:,0]), np.log(diff1_data[:,1]))
	slope_perp.append(popt[0])
	diff2_data = np.loadtxt("output_{:s}/CIR/diff2_{:s}_{:s}.dat".format(args.date, direction, args.date))
	popt, pcov = curve_fit(line, np.log(diff2_data[:,0]), np.log(diff2_data[:,1]))
	slope_para.append(popt[0])
print("mean kappa_perp radial power slope: {:.2f}".format(np.mean(slope_perp)))
print("mean kappa_para radial power slope: {:.2f}".format(np.mean(slope_para)))
slope_perp = []
slope_para = []
for direction in directions:
	diff1_data = np.loadtxt("output_{:s}/CIR/diff1_E_{:s}_{:s}.dat".format(args.date, direction, args.date))
	popt, pcov = curve_fit(line, np.log(diff1_data[:,0]), np.log(diff1_data[:,1]))
	slope_perp.append(popt[0])
	diff2_data = np.loadtxt("output_{:s}/CIR/diff2_E_{:s}_{:s}.dat".format(args.date, direction, args.date))
	popt, pcov = curve_fit(line, np.log(diff2_data[:,0]), np.log(diff2_data[:,1]))
	slope_para.append(popt[0])
print("mean kappa_perp energy power slope: {:.2f}".format(np.mean(slope_perp)))
print("mean kappa_para energy power slope: {:.2f}".format(np.mean(slope_para)))
