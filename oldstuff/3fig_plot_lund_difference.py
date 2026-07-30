import sys
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

if len(sys.argv) < 3:
    print(f"Usage: python {sys.argv[0]} reference.dat comparison1.dat [comparison2.dat ...]")
    sys.exit(1)

reference_file = sys.argv[1]
comparison_files = sys.argv[2:]

def plot_name(fname):
    name = Path(fname).stem
    if name.endswith("_lund"):
        name = name[:-5]
    return name

xbins = np.linspace(0.0, 5.0, 80)
ybins = np.linspace(-2.0, 6.0, 80)

comparison_hists = {}

for fname in comparison_files:
    Delta, kt = np.loadtxt(fname, comments="#", usecols=(4, 5), unpack=True)
    h, _, _ = np.histogram2d(np.log(1/Delta), np.log(kt), bins=(xbins, ybins))
    comparison_hists[plot_name(fname)] = h / h.sum()

Delta, kt = np.loadtxt(reference_file, comments="#", usecols=(4, 5), unpack=True)
h, _, _ = np.histogram2d(np.log(1/Delta), np.log(kt), bins=(xbins, ybins))
reference = h / h.sum()
reference_name = plot_name(reference_file)

fig, axs = plt.subplots(1, len(comparison_files), figsize=(4 * len(comparison_files), 4), sharex=True, sharey=True)
if len(comparison_files) == 1:
    axs = [axs]

for ax, (label, h) in zip(axs, comparison_hists.items()):
    diff = h - reference
    im = ax.pcolormesh(xbins, ybins, diff.T, cmap="RdBu_r", vmin=-diff.max(), vmax=diff.max())
    title = f"{label} - {reference_name}"
    title_size = min(12, 500 / max(len(title), 1))
    ax.set_title(title, fontsize=title_size)
    ax.set_xlabel(r"$\ln(1/\Delta)$")

axs[0].set_ylabel(r"$\ln(k_t/\mathrm{GeV})$")
fig.colorbar(im, ax=axs, label="normalised entry difference")
#plt.tight_layout()
plt.savefig("3fig_lund_plane_differences.pdf")
