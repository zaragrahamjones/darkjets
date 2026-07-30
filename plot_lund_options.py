import argparse
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

def format_value(value):
    return f"{value:g}".replace(".", "p")


def read_event_metadata(input_file):
    event_file = input_file.with_suffix(".eventdat")
    with event_file.open() as f:
        f.readline()
        values = f.readline().lstrip("#").split()
        f.readline()
        last_event = 0
        for line in f:
            if line.strip() and not line.startswith("#"):
                last_event = int(line.split()[0])
    return float(values[1]), last_event + 1


def load_jet_data(input_file):
    jet_file = input_file.with_suffix(".jetdat")
    event, jet, pt, eta, dark_efrac, dark_ptfrac, n_lund = np.loadtxt(
        jet_file, comments="#", usecols=(0, 1, 2, 3, 4, 5, 6), unpack=True
    )
    return event.astype(int), jet.astype(int), pt, eta, dark_efrac, dark_ptfrac, n_lund.astype(int)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Plot selected Lund plane entries from a .lunddat file."
    )
    parser.add_argument("input_file", type=Path)
    parser.add_argument(
        "--jets",
        choices=("hardest", "all", "pt-min"),
        default="all",
        help="Choose jet index 0 only, all jets, or jets above --jet-pt-min.",
    )
    parser.add_argument(
        "--jet-pt-min",
        type=float,
        default=None,
        help="Jet pT threshold used when --jets pt-min. Defaults to the matching .eventdat header.",
    )
    parser.add_argument(
        "--dm-selection",
        choices=("dm", "nondm", "all"),
        default="all",
        help="Select DM-classified, nonDM-classified, or all jets.",
    )
    parser.add_argument(
        "--dm-mask",
        type=float,
        default=0.75,
        help="Dark energy fraction threshold for DM classification.",
    )
    return parser.parse_args()


args = parse_args()
input_file = args.input_file
event_jet_pt_min, n_events = read_event_metadata(input_file)
jet_pt_min = args.jet_pt_min
if jet_pt_min is None:
    jet_pt_min = event_jet_pt_min

name_parts = [
    input_file.stem,
    "lundplane",
    f"jets-{args.jets}",
]
if args.jets == "pt-min":
    name_parts.append(f"jetptmin-{format_value(jet_pt_min)}")
if args.dm_selection != "all":
    name_parts.append(f"{args.dm_selection}-mask-{format_value(args.dm_mask)}")
else:
    name_parts.append("dm-all")
output_file = input_file.with_suffix("").with_name("_".join(name_parts) + ".pdf")

xbins = np.linspace(1.0, 6.0, 80)
ybins = np.linspace(-4.0,6.0, 80)

jet_event, jet_index, jet_pt, _jet_eta, jet_dark_efrac, _jet_dark_ptfrac, _n_lund = load_jet_data(input_file)
event, jet, x, y = np.loadtxt(
    input_file, comments="#", skiprows=1, usecols=(0, 1, 4, 5), unpack=True
)
event = event.astype(int)
jet = jet.astype(int)

if args.jets == "hardest":
    jet_mask = jet_index == 0
elif args.jets == "pt-min":
    jet_mask = jet_pt > jet_pt_min
else:
    jet_mask = np.ones_like(jet_index, dtype=bool)

if args.dm_selection == "dm":
    jet_mask &= jet_dark_efrac > args.dm_mask
elif args.dm_selection == "nondm":
    jet_mask &= jet_dark_efrac <= args.dm_mask

selected_jets = set(zip(jet_event[jet_mask], jet_index[jet_mask]))
n_jets = len(selected_jets)
lund_mask = np.array([(event_id, jet_id) in selected_jets
                      for event_id, jet_id in zip(event, jet)])
x = x[lund_mask]
y = y[lund_mask]

hist, xedges, yedges = np.histogram2d(x, y, bins=(xbins, ybins))
bin_area = np.diff(xedges)[:, None] * np.diff(yedges)[None, :]
if n_jets > 0:
    density = hist / (n_jets * bin_area)
else:
    density = hist

mesh = plt.pcolormesh(xedges, yedges, density.T, cmap="magma")
plt.xlabel(r"$\ln(1/\Delta)$")
plt.ylabel(r"$\ln(k_t/\mathrm{GeV})$")
plt.colorbar(mesh, label=r"$\rho(\Delta, k_t)$")
jet_label = {
    "hardest": "hardest jet only",
    "all": rf"all jets ($p_T > {jet_pt_min:g}$ GeV)",
    "pt-min": rf"jets with postapplied $p_T > {jet_pt_min:g}$ GeV",
}[args.jets]
dm_label = {
    "dm": f"dark hadron energy fraction in jet > {args.dm_mask:g}",
    "nondm": f"non-dark hadron energy fraction in jet <= {args.dm_mask:g}",
    "all": "no dark hadron jet classification",
}[args.dm_selection]
norm_label = f"normalised by {n_jets} selected jets, from {n_events} events"
title = f"Lund density plane for\n{input_file.stem}\n{norm_label}\n{jet_label}; {dm_label}"
#title_size = min(12, 500 / max(len(title), 1))
title_size = 10
plt.title(title, fontsize=title_size)
plt.tight_layout(rect=(0, 0, 1, 0.95))
plt.savefig(output_file)
