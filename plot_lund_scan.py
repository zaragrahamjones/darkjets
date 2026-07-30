import argparse
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
import numpy as np


# Match this vector to the variables in the .eventdat header:
# etacut jetPtMin Z'width Lambda m_qv
# Use 1 for the two variables to scan and 0 for variables kept fixed.
SCAN_VARIABLES = [0, 0, 1, 1, 0]

AXIS_LABELS = {
    "etacut": r"$|\eta|$ cut",
    "jetPtMin": r"jet $p_T^{\min}$ [GeV]",
    "Z'width": r"$Z'$ width",
    "Lambda": r"$\Lambda$ [GeV]",
    "m_qv": r"$m_{q_v}$ [GeV]",
}


def format_value(value):
    return f"{value:g}".replace(".", "p")


def read_event_metadata(input_file):
    event_file = input_file.with_suffix(".eventdat")
    with event_file.open() as f:
        names = f.readline().lstrip("#").split()
        values = [float(value) for value in f.readline().lstrip("#").split()]
        f.readline()
        last_event = 0
        for line in f:
            if line.strip() and not line.startswith("#"):
                last_event = int(line.split()[0])
    return dict(zip(names, values)), values[1], last_event + 1


def lund_density(input_file, jet_pt_min, n_events, args):
    event, jet, dmfrac, clust_pt, x, y = np.loadtxt(
        input_file, comments="#", skiprows=1, usecols=(0, 1, 2, 5, 6, 7), unpack=True
    )

    if args.jets == "hardest":
        mask = jet == 0
    elif args.jets == "pt-min":
        mask = np.zeros_like(jet, dtype=bool)
        for event_id, jet_id in set(zip(event.astype(int), jet.astype(int))):
            jet_mask = (event == event_id) & (jet == jet_id)
            if np.max(clust_pt[jet_mask]) > jet_pt_min:
                mask |= jet_mask
    else:
        mask = np.ones_like(jet, dtype=bool)

    if args.dm_selection == "dm":
        mask &= dmfrac > args.dm_mask
    elif args.dm_selection == "nondm":
        mask &= dmfrac <= args.dm_mask

    if args.dm_selection == "dm" and args.dm_normalisation == "dm-events":
        norm_events = len(np.unique(event[mask]))
    else:
        norm_events = n_events

    hist, xedges, yedges = np.histogram2d(
        x[mask],
        y[mask],
        bins=(np.linspace(1.0, 6.0, 80), np.linspace(-4.0, 6.0, 80)),
    )
    bin_area = np.diff(xedges)[:, None] * np.diff(yedges)[None, :]
    if norm_events > 0:
        density = hist / (norm_events * bin_area)
    else:
        density = hist

    return density, xedges, yedges


def parse_args():
    parser = argparse.ArgumentParser(
        description="Plot a scan of Lund planes for a two-variable scan."
    )
    parser.add_argument("sample_dir", type=Path)
    parser.add_argument(
        "--sample-prefix",
        default=None,
        help="Dataset prefix to select, for example hvZ2. Defaults to all .lunddat files.",
    )
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
        help="Jet pT threshold used when --jets pt-min. Defaults to each matching .eventdat header.",
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
    parser.add_argument(
        "--dm-normalisation",
        choices=("total-events", "dm-events"),
        default="total-events",
        help="For --dm-selection dm, normalise by all events or by events with selected DM entries.",
    )
    return parser.parse_args()


args = parse_args()

if sum(SCAN_VARIABLES) != 2:
    raise ValueError("SCAN_VARIABLES must contain exactly two entries set to 1.")

pattern = "*.lunddat"
if args.sample_prefix is not None:
    pattern = f"{args.sample_prefix}_*.lunddat"

datasets = []
scan_names = None
for input_file in sorted(args.sample_dir.glob(pattern)):
    metadata, event_jet_pt_min, n_events = read_event_metadata(input_file)
    names = list(metadata)
    selected_names = [name for name, scan in zip(names, SCAN_VARIABLES) if scan]
    if scan_names is None:
        scan_names = selected_names
    elif selected_names != scan_names:
        raise ValueError("All event metadata headers must use the same variable order.")

    jet_pt_min = args.jet_pt_min
    if jet_pt_min is None:
        jet_pt_min = event_jet_pt_min

    density, xedges, yedges = lund_density(input_file, jet_pt_min, n_events, args)
    datasets.append(
        {
            "file": input_file,
            "metadata": metadata,
            "density": density,
            "xedges": xedges,
            "yedges": yedges,
        }
    )

if not datasets:
    raise ValueError(f"No .lunddat files found in {args.sample_dir}.")

x_scan, y_scan = scan_names
x_values = sorted({dataset["metadata"][x_scan] for dataset in datasets})
y_values = sorted({dataset["metadata"][y_scan] for dataset in datasets}, reverse=True)
density_max = max(np.max(dataset["density"]) for dataset in datasets)

fig, axes = plt.subplots(
    len(y_values),
    len(x_values),
    figsize=(3.0 * len(x_values), 2.7 * len(y_values)),
    squeeze=False,
    sharex=True,
    sharey=True,
)

mesh = None
for row, y_value in enumerate(y_values):
    for col, x_value in enumerate(x_values):
        ax = axes[row, col]
        match = [
            dataset
            for dataset in datasets
            if dataset["metadata"][x_scan] == x_value
            and dataset["metadata"][y_scan] == y_value
        ]
        if not match:
            ax.axis("off")
            continue

        dataset = match[0]
        mesh = ax.pcolormesh(
            dataset["xedges"],
            dataset["yedges"],
            dataset["density"].T,
            cmap="magma",
            vmin=0,
            vmax=density_max,
        )
        ax.set_title(
            f"{x_scan}={x_value:g}, {y_scan}={y_value:g}",
            fontsize=9,
        )

        if row == len(y_values) - 1:
            ax.set_xlabel(r"$\ln(1/\Delta)$")
        if col == 0:
            ax.set_ylabel(r"$\ln(k_t/\mathrm{GeV})$")

if mesh is not None:
    colorbar = fig.colorbar(mesh, ax=axes, location="right", pad=0.02)
    colorbar.set_label(r"Lund density per event per Lund-plane area")
    colorbar.ax.yaxis.set_major_formatter(FormatStrFormatter("%.2g"))

x_label = AXIS_LABELS.get(x_scan, x_scan)
y_label = AXIS_LABELS.get(y_scan, y_scan)
fig.supxlabel(x_label)
fig.supylabel(y_label)

name_parts = [
    args.sample_dir.name,
    "lundscan",
    x_scan.replace("'", ""),
    "vs",
    y_scan.replace("'", ""),
    f"jets-{args.jets}",
]
if args.jets == "pt-min" and args.jet_pt_min is not None:
    name_parts.append(f"jetptmin-{format_value(args.jet_pt_min)}")
if args.dm_selection != "all":
    name_parts.append(f"{args.dm_selection}-mask-{format_value(args.dm_mask)}")
else:
    name_parts.append("dm-all")
if args.dm_selection == "dm":
    name_parts.append(f"norm-{args.dm_normalisation}")

output_file = args.sample_dir / ("_".join(name_parts) + ".pdf")
fig.suptitle(f"Lund-plane scan: {x_label} vs {y_label}", fontsize=12)
fig.tight_layout(rect=(0.04, 0.04, 0.92, 0.95))
fig.savefig(output_file)
