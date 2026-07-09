import numpy as np
import matplotlib.pyplot as plt

samples = {
    "SM Z(inv)+j": "lund_sm.dat",
    "dark g=0.1": "lund_dark_g0p1.dat",
    "dark g=0.5": "lund_dark_g0p5.dat",
    "dark g=1.0": "lund_dark_g1p0.dat",
}

xbins = np.linspace(0.0, 5.0, 80)
ybins = np.linspace(-2.0, 6.0, 80)

# Bin centres for contour plots
xc = 0.5 * (xbins[:-1] + xbins[1:])
yc = 0.5 * (ybins[:-1] + ybins[1:])

hists = {}

for label, fname in samples.items():
    x, y = np.loadtxt(fname, comments="#", usecols=(2, 3), unpack=True)
    h, _, _ = np.histogram2d(x, y, bins=(xbins, ybins))
    hists[label] = h

sm = hists["SM Z(inv)+j"]

fig, axs = plt.subplots(1, 4, figsize=(16, 4), sharex=True, sharey=True)

# Plot the SM Lund plane
im_sm = axs[0].pcolormesh(xbins, ybins, sm.T, cmap="inferno", shading="auto")
axs[0].set_title("SM Z(inv)+j")
axs[0].set_xlabel(r"$\ln(1/\Delta)$")
axs[0].set_ylabel(r"$\ln(k_t/\mathrm{GeV})$")

# Plot differences for the DM samples
dm_labels = [k for k in samples if k != "SM Z(inv)+j"]

for ax, dm_label in zip(axs[1:], dm_labels):
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

    ax.contour(xc, yc, sm.T,
            levels=3,
            cmap="inferno",
            linewidths=2,
            alpha=1,
            label="SM")

    ax.contour(xc, yc, hists[dm_label].T,
            levels=3,
            cmap="inferno",
            linewidths=2,
            alpha=1,
            label=dm_label)    
    
    ax.set_title(f"{dm_label} vs SM")
    ax.set_xlabel(r"$\ln(1/\Delta)$")
    ax.set_ylabel(r"$\ln(k_t/\mathrm{GeV})$")

fig.colorbar(im_sm, ax=axs[0], label="Entries")

#plt.tight_layout()
plt.savefig("4fig_lund_plane_overlay.pdf")