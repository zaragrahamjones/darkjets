import sys
import textwrap
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

# use python3 plot_lund.py hv.lunddat to run this script, with output hv_lundplane.pdf
input_file = Path(sys.argv[1])
output_file = input_file.with_suffix("").with_name(input_file.stem + "_lundplane.pdf")

xbins = np.linspace(1.0, 6.0, 80)
ybins = np.linspace(-4.0,6.0, 80)

# plot
x, y = np.loadtxt(input_file, comments="#", skiprows=1, usecols=(6,7), unpack=True)

plt.hist2d(x, y, bins=(xbins, ybins), cmap="magma")
plt.xlabel(r"$\ln(1/\Delta)$")
plt.ylabel(r"$\ln(k_t/\mathrm{GeV})$")
plt.colorbar(label="entries")
title = f"Lund plane for\n{input_file.stem}"
title_size = min(12, 500 / max(len(input_file.stem), 1))
plt.title(title, fontsize=title_size)
plt.tight_layout(rect=(0, 0, 1, 0.95))
plt.savefig(output_file)
