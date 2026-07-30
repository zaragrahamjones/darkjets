#!/usr/bin/env python3

import argparse
import sys
import textwrap
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

input_files = [Path(arg) for arg in sys.argv[1:]]
output_file = input_files[0].with_suffix("").with_name(
    "_".join(input_file.stem for input_file in input_files) + "_met.pdf")

mets = [np.loadtxt(input_file, comments="#", usecols=(3), unpack=True)
        for input_file in input_files]

bins = np.linspace(0, max(np.max(met) for met in mets), 51)

plt.figure(figsize=(7, 5))
for input_file, met in zip(input_files, mets):
    plt.hist(met, bins=bins, histtype="step", linewidth=1.8,
             label=input_file.stem)
plt.xlabel("MET [GeV]")
plt.ylabel("log( Number of events )")
plt.yscale("log")
plt.legend()
plt.tight_layout()
plt.title(textwrap.fill("Missing transverse energy at detector level", 60))
plt.savefig(output_file, dpi=200)
