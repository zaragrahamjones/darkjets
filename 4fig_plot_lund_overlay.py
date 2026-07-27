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

# Bin centres for contour plots
xc = 0.5 * (xbins[:-1] + xbins[1:])
yc = 0.5 * (ybins[:-1] + ybins[1:])

hists = {}

for fname in [reference_file] + comparison_files:
    Delta, kt = np.loadtxt(fname, comments="#", usecols=(4, 5), unpack=True)
    h, _, _ = np.histogram2d(np.log(1/Delta), np.log(kt), bins=(xbins, ybins))
    hists[plot_name(fname)] = h

reference_name = plot_name(reference_file)
reference = hists[reference_name]

fig, axs = plt.subplots(1, len(comparison_files) + 1, figsize=(4 * (len(comparison_files) + 1), 4), sharex=True, sharey=True)
if len(comparison_files) == 0:
    axs = [axs]

# Plot the reference Lund plane
im_sm = axs[0].pcolormesh(xbins, ybins, reference.T, cmap="inferno", shading="auto")
title_size = min(12, 500 / max(len(reference_name), 1))
axs[0].set_title(reference_name, fontsize=title_size)
axs[0].set_xlabel(r"$\ln(1/\Delta)$")
axs[0].set_ylabel(r"$\ln(k_t/\mathrm{GeV})$")

# Plot overlays for the comparison samples
comparison_labels = [plot_name(fname) for fname in comparison_files]

for ax, comparison_label in zip(axs[1:], comparison_labels):
    # ax.pcolormesh(xbins, ybins, sm.T,
    #             cmap="inferno",
    #             shading="auto",
    #             alpha=1,
    #             label="SM")

    # ax.pcolormesh(xbins, ybins, hists[dm_label].T,
    #             cmap="inferno",
    #             shading="auto",
    #             alpha=0.5,
    #             label=dm_label)

    ax.contour(xc, yc, reference.T,
            levels=3,
            cmap="inferno",
            linewidths=2,
            alpha=1,
            label=reference_name)

    ax.contour(xc, yc, hists[comparison_label].T,
            levels=3,
            cmap="inferno",
            linewidths=2,
            alpha=1,
            label=comparison_label)    
    
    title = f"{comparison_label} vs {reference_name}"
    title_size = min(12, 500 / max(len(title), 1))
    ax.set_title(title, fontsize=title_size)
    ax.set_xlabel(r"$\ln(1/\Delta)$")
    ax.set_ylabel(r"$\ln(k_t/\mathrm{GeV})$")

fig.colorbar(im_sm, ax=axs[0], label="Entries")

#plt.tight_layout()
plt.savefig("4fig_lund_plane_overlay.pdf")
