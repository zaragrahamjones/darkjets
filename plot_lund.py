import sys
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

# use python3 plot_lund.py lund01.dat to run this script, with output lund01_lundplane.pdf
input_file = Path(sys.argv[1])
output_file = input_file.with_suffix("").with_name(input_file.stem + "_lundplane.pdf")

# plot
Delta, kt = np.loadtxt(input_file, comments="#", usecols=(2, 3), unpack=True)

plt.hist2d(np.log(1. / Delta), np.log(kt), bins=80, cmap="magma")
plt.xlabel(r"$\ln(1/\Delta)$")
plt.ylabel(r"$\ln(k_t/\mathrm{GeV})$")
plt.colorbar(label="entries")
plt.title(f"Lund plane for {input_file.stem}")
plt.tight_layout()
plt.savefig(output_file)