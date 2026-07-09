import sys
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

if len(sys.argv) != 3:
    print(f"Usage: python {sys.argv[0]} file1.dat file2.dat")
    sys.exit(1)

file1 = sys.argv[1]
file2 = sys.argv[2]

xbins = np.linspace(0.0, 5.0, 80)
ybins = np.linspace(-2.0, 6.0, 80)

def load_hist(fname):
    x, y = np.loadtxt(fname, comments="#", usecols=(2, 3), unpack=True)
    h, _, _ = np.histogram2d(x, y, bins=(xbins, ybins))
    return h / h.sum()

h1 = load_hist(file1)
h2 = load_hist(file2)

diff = h2 - h1
vmax = np.abs(diff).max()

fig, ax = plt.subplots(figsize=(5, 4))

im = ax.pcolormesh(
    xbins, ybins, diff.T,
    cmap="RdBu_r",
    vmin=-vmax,
    vmax=vmax,
)

name1 = Path(file1).stem
name2 = Path(file2).stem

ax.set_title(f"{name2} - {name1}")
ax.set_xlabel(r"$\ln(1/\Delta)$")
ax.set_ylabel(r"$\ln(k_t/\mathrm{GeV})$")

fig.colorbar(im, ax=ax, label="Normalised entry difference")
plt.tight_layout()

outfile = f"{name2}_minus_{name1}.pdf"
plt.savefig(outfile)
print(f"Saved {outfile}")