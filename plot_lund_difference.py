import sys
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

if len(sys.argv) != 3:
    print(f"Usage: python {sys.argv[0]} file1.dat file2.dat")
    sys.exit(1)

file1 = sys.argv[1]
file2 = sys.argv[2]

def count_events(fname):
    last_event = -1

    with open(fname) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue

            event = int(line.split()[0])
            last_event = event

    return last_event + 1

def check_event_numbers(file1, file2):
    n_events1 = count_events(file1)
    n_events2 = count_events(file2)

    if n_events1 != n_events2:
        print(f"Warning: {file1} has {n_events1} events, but {file2} has {n_events2} events.")

    return min(n_events1, n_events2)

n_events = check_event_numbers(file1, file2)

xbins = np.linspace(0.8, 8.0, 80)
ybins = np.linspace(-4.0,8.0, 80)

def load_hist(fname, max_events):
    Delta = []
    kt = []

    with open(fname) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue

            cols = line.split()
            event = int(cols[0])

            if event >= max_events:
                break

            Delta.append(float(cols[2]))
            kt.append(float(cols[3]))
    
    ln_oneoverDelta = np.log(1.0 / np.array(Delta))
    ln_kt = np.log(np.array(kt))

    h, _, _ = np.histogram2d(ln_oneoverDelta, ln_kt, bins=(xbins, ybins))
    return h

h1 = load_hist(file1, n_events)
h2 = load_hist(file2, n_events)

h1 /= n_events
h2 /= n_events

diff = h2 - h1
vmax = np.abs(diff).max()

ratio = diff/np.sqrt(h1 + h2 + 1.0e-12)  # Avoid division by zero
ratiomax = np.abs(ratio).max()

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

ax.set_title(f"({name2} - {name1}) / {n_events} events")
ax.set_xlabel(r"$\ln(1/\Delta)$")
ax.set_ylabel(r"$\ln(k_t/\mathrm{GeV})$")
fig.colorbar(im1, ax=ax, label="Normalised entry difference")

plt.tight_layout()

outfile = f"{name2}_minus_{name1}.pdf"
plt.savefig(outfile)
print(f"Saved {outfile}")