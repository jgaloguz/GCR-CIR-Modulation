# Code to compute the average diffusion coefficients over a rotation

# Import libraries
import numpy as np
import sys

# Check CIR date
if len(sys.argv) > 1:
   cir_date = sys.argv[1]
   print("CIR date: {:s}".format(cir_date))
else:
   print("ERROR: No CIR date provided.")
   sys.exit(1)

# Import data
diff1_data = np.loadtxt("output_{:s}/CIR/diff1_1au_{:s}.dat".format(cir_date, cir_date))
diff1_magn = diff1_data[:,1]
diff2_data = np.loadtxt("output_{:s}/CIR/diff2_1au_{:s}.dat".format(cir_date, cir_date))
diff2_magn = diff2_data[:,1]

# Report averages
print("kappa_perp arithmetic average (1au) = {:.2e} cm^2 s^-1".format(np.mean(diff1_magn)))
print("kappa_para arithmetic average (1au) = {:.2e} cm^2 s^-1".format(np.mean(diff2_magn)))

print("kappa_perp geometric average (1au) = {:.2e} cm^2 s^-1".format(np.exp(np.mean(np.log(diff1_magn)))))
print("kappa_para geometric average (1au) = {:.2e} cm^2 s^-1".format(np.exp(np.mean(np.log(diff2_magn)))))