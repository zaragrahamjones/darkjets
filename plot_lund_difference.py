import argparse
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

def format_value(value):
    return f"{value:g}".replace(".", "p")


def parse_args():
    parser = argparse.ArgumentParser(
        description="Plot the difference between two Lund density histograms."
    )
    parser.add_argument("file1", type=Path)
    parser.add_argument("file2", type=Path)
    parser.add_argument(
        "--jets1",
        choices=("all", "hardest", "pt-min"),
        default="all",
        help="Jet mask applied to file1 before histogramming.",
    )
    parser.add_argument(
        "--jet-pt-min1",
        type=float,
        default=None,
        help="Jet pT threshold for --jets1 pt-min. Defaults to file1's matching .eventdat header.",
    )
    parser.add_argument(
        "--jets2",
        choices=("all", "hardest", "pt-min"),
        default="all",
        help="Jet mask applied to file2 before histogramming.",
    )
    parser.add_argument(
        "--jet-pt-min2",
        type=float,
        default=None,
        help="Jet pT threshold for --jets2 pt-min. Defaults to file2's matching .eventdat header.",
    )
    parser.add_argument(
        "--mask1",
        choices=("none", "dm", "nondm"),
        default="none",
        help="DM mask applied to file1 before histogramming.",
    )
    parser.add_argument(
        "--dm-mask1",
        type=float,
        default=0.75,
        help="Dark energy fraction threshold for --mask1 dm/nondm.",
    )
    parser.add_argument(
        "--mask2",
        choices=("none", "dm", "nondm"),
        default="none",
        help="DM mask applied to file2 before histogramming.",
    )
    parser.add_argument(
        "--dm-mask2",
        type=float,
        default=0.75,
        help="Dark energy fraction threshold for --mask2 dm/nondm.",
    )
    return parser.parse_args()


args = parse_args()
file1 = args.file1
file2 = args.file2

def read_event_metadata(fname):
    event_file = fname.with_suffix(".eventdat")
    with event_file.open() as f:
        f.readline()
        values = f.readline().lstrip("#").split()
        f.readline()
        last_event = -1
        for line in f:
            if line.strip() and not line.startswith("#"):
                last_event = int(line.split()[0])
    return float(values[1]), last_event + 1


def load_jet_data(fname):
    jet_file = fname.with_suffix(".jetdat")
    event, jet, pt, eta, dark_efrac, dark_ptfrac, n_lund = np.loadtxt(
        jet_file, comments="#", usecols=(0, 1, 2, 3, 4, 5, 6), unpack=True
    )
    return event.astype(int), jet.astype(int), pt, eta, dark_efrac, dark_ptfrac, n_lund.astype(int)


event_jet_pt_min1, n_events1 = read_event_metadata(file1)
event_jet_pt_min2, n_events2 = read_event_metadata(file2)
if n_events1 != n_events2:
    raise ValueError(f"Number of events in {file1} ({n_events1}) and {file2} ({n_events2}) do not match.")

jet_pt_min1 = args.jet_pt_min1
jet_pt_min2 = args.jet_pt_min2
if jet_pt_min1 is None:
    jet_pt_min1 = event_jet_pt_min1
if jet_pt_min2 is None:
    jet_pt_min2 = event_jet_pt_min2

xbins = np.linspace(0.8, 8.0, 80)
ybins = np.linspace(-4.0,8.0, 80)

def jet_suffix(jet_kind, jet_pt_min):
    if jet_kind == "all":
        return ""
    if jet_kind == "pt-min":
        return f"jets-ptmin-{format_value(jet_pt_min)}"
    return f"jets-{jet_kind}"


def jet_label(jet_kind, jet_pt_min):
    if jet_kind == "hardest":
        return "hardest jet only"
    if jet_kind == "pt-min":
        return rf"jets with postapplied $p_T > {jet_pt_min:g}$ GeV"
    return ""


def mask_suffix(mask_kind, mask_value):
    if mask_kind == "none":
        return ""
    return f"{mask_kind}-mask-{format_value(mask_value)}"


def mask_label(mask_kind, mask_value):
    if mask_kind == "dm":
        return f"DM fraction > {mask_value:g}"
    if mask_kind == "nondm":
        return f"DM fraction <= {mask_value:g}"
    return ""


def load_density(fname, jet_kind, jet_pt_min, mask_kind, mask_value):
    jet_event, jet_index, jet_pt, _jet_eta, jet_dark_efrac, _jet_dark_ptfrac, _n_lund = load_jet_data(fname)
    event, jet, x, y = np.loadtxt(
        fname, comments="#", skiprows=1, usecols=(0, 1, 4, 5), unpack=True
    )
    event = event.astype(int)
    jet = jet.astype(int)

    if jet_kind == "hardest":
        jet_selection = jet_index == 0
    elif jet_kind == "pt-min":
        jet_selection = jet_pt > jet_pt_min
    else:
        jet_selection = np.ones_like(jet_index, dtype=bool)

    if mask_kind == "dm":
        jet_selection &= jet_dark_efrac > mask_value
    elif mask_kind == "nondm":
        jet_selection &= jet_dark_efrac <= mask_value

    selected_jets = set(zip(jet_event[jet_selection], jet_index[jet_selection]))
    n_jets = len(selected_jets)
    lund_selection = np.array([(event_id, jet_id) in selected_jets
                               for event_id, jet_id in zip(event, jet)])
    x = x[lund_selection]
    y = y[lund_selection]

    h, _, _ = np.histogram2d(x, y, bins=(xbins, ybins))
    bin_area = np.diff(xbins)[:, None] * np.diff(ybins)[None, :]
    if n_jets > 0:
        return h / (n_jets * bin_area), n_jets
    return h, n_jets

h1, n_jets1 = load_density(file1, args.jets1, jet_pt_min1,
                           args.mask1, args.dm_mask1)
h2, n_jets2 = load_density(file2, args.jets2, jet_pt_min2,
                           args.mask2, args.dm_mask2)

diff = h2 - h1
vmax = np.abs(diff).max()

fig, ax = plt.subplots(1, 1, figsize=(5, 4))
name1 = Path(file1).stem
name2 = Path(file2).stem

# Difference
im1 = ax.pcolormesh(
    xbins, ybins, diff.T,
    cmap="RdBu_r",
    vmin=-vmax,
    vmax=vmax,
)

plot_name = f"{name2}_minus_{name1}"
if plot_name.endswith("_lund"):
    plot_name = plot_name[:-5]
event_line = f"{n_jets2} jets from {n_events2} events minus {n_jets1} jets from {n_events1} events"
mask_lines = []
jet1_label = jet_label(args.jets1, jet_pt_min1)
jet2_label = jet_label(args.jets2, jet_pt_min2)
mask1_label = mask_label(args.mask1, args.dm_mask1)
mask2_label = mask_label(args.mask2, args.dm_mask2)
if jet1_label:
    mask_lines.append(f"{name1}: {jet1_label}")
if jet2_label:
    mask_lines.append(f"{name2}: {jet2_label}")
if mask1_label:
    mask_lines.append(f"{name1}: {mask1_label}")
if mask2_label:
    mask_lines.append(f"{name2}: {mask2_label}")
title_lines = [plot_name, event_line] + mask_lines
title = "\n".join(title_lines)
title_size = min(12, 500 / max(len(line) for line in title_lines))
ax.set_title(title, fontsize=title_size)
ax.set_xlabel(r"$\ln(1/\Delta)$")
ax.set_ylabel(r"$\ln(k_t/\mathrm{GeV})$")
fig.colorbar(im1, ax=ax, label=r"$\rho(\Delta, k_t)$ difference")

plt.tight_layout(rect=(0, 0, 1, 0.92))

outfile_parts = [f"{name2}_minus_{name1}"]
outfile_parts.append(f"file1-jets-{args.jets1}")
if args.jets1 == "pt-min":
    outfile_parts.append(f"file1-jetptmin-{format_value(jet_pt_min1)}")
outfile_parts.append(f"file2-jets-{args.jets2}")
if args.jets2 == "pt-min":
    outfile_parts.append(f"file2-jetptmin-{format_value(jet_pt_min2)}")
if args.mask1 == "none":
    outfile_parts.append("file1-mask-none")
else:
    outfile_parts.append(f"file1-{args.mask1}-mask-{format_value(args.dm_mask1)}")
if args.mask2 == "none":
    outfile_parts.append("file2-mask-none")
else:
    outfile_parts.append(f"file2-{args.mask2}-mask-{format_value(args.dm_mask2)}")
outfile = "_".join(outfile_parts) + ".pdf"
plt.savefig(outfile)
print(f"Saved {outfile}")
