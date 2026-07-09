import sys
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

# use python3 plot_lund.py lund01.dat to run this script, with output lund01_plane.pdf
input_file = Path(sys.argv[1])
output_file = input_file.with_suffix("").with_name(input_file.stem + "_plane.pdf")

# plot
x, y = np.loadtxt(input_file, comments="#", usecols=(2, 3), unpack=True)

plt.hist2d(x, y, bins=80, cmap="magma")
plt.xlabel(r"$\ln(1/\Delta)$")
plt.ylabel(r"$\ln(k_t/\mathrm{GeV})$")
plt.colorbar(label="entries")
plt.tight_layout()
plt.savefig(output_file)