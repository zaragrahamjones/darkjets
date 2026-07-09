import numpy as np
import matplotlib.pyplot as plt

SM_samples = {
    "SM Z(inv)+j": "lund_sm.dat",
}
DM_samples = {
    "dark g=0.1": "lund_dark_g0p1.dat",
    "dark g=0.5": "lund_dark_g0p5.dat",
    "dark g=1.0": "lund_dark_g1p0.dat",
}

xbins = np.linspace(0.0, 5.0, 80)
ybins = np.linspace(-2.0, 6.0, 80)

DM_hists = {}
SM_hists = {}

for label, fname in DM_samples.items():
    x, y = np.loadtxt(fname, comments="#", usecols=(2, 3), unpack=True)
    h, _, _ = np.histogram2d(x, y, bins=(xbins, ybins))
    DM_hists[label] = h / h.sum()
for label, fname in SM_samples.items():
    x, y = np.loadtxt(fname, comments="#", usecols=(2, 3), unpack=True)
    h, _, _ = np.histogram2d(x, y, bins=(xbins, ybins))
    SM_hists[label] = h / h.sum()

sm = SM_hists["SM Z(inv)+j"]

fig, axs = plt.subplots(1, len(DM_samples), figsize=(4 * len(DM_samples), 4), sharex=True, sharey=True)

for ax, (label, h) in zip(axs, DM_hists.items()):
    diff = h - sm
    im = ax.pcolormesh(xbins, ybins, diff.T, cmap="RdBu_r", vmin=-diff.max(), vmax=diff.max())
    ax.set_title(label + " - SM")
    ax.set_xlabel(r"$\ln(1/\Delta)$")

axs[0].set_ylabel(r"$\ln(k_t/\mathrm{GeV})$")
fig.colorbar(im, ax=axs, label="normalised entry difference")
#plt.tight_layout()
plt.savefig("3fig_lund_plane_differences.pdf")