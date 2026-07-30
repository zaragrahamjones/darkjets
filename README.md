# Dark Jets Hidden Valley Workflow

A small Pythia 8/FastJet analysis program for Hidden Valley dark-jet studies, plus scripts for scanning model parameters and plotting Lund-plane densities. This will be fed into a ML workflow for training a neural network to distinguish dark jets from QCD jets.

main executable:
```bash
./hiddenvalley
```

## Running

Basic usage:

```bash
./hiddenvalley -e 10000 -o hvZdefault hvZ.cmnd
```

With an override card:

```bash
./hiddenvalley -e 10k -o hvZ2_Lambda75_width5_mq2p5 hvZ.cmnd hvZ_var.cmnd
```

Command-line options:

| Option | Meaning |
|---|---|
| `-e` | Number of events. Suffixes `k`, `m`, and `g` are accepted, e.g. `10k`. |
| `-o` | Output stem. Files are written as `<stem>.eventdat`, `<stem>.jetdat`, and `<stem>.lunddat`. |
| `-s` | If omitted, no explicit seed is set. `-s i` uses seed `i` for integer `i>0`, else random. |

The executable is intentionally quiet on successful runs, in preparation for batch submission. Check log files for errors.

## Scan Script

`scan_hvZ.sh` scans over Z' width, Hidden Valley confinement scale, and dark quark mass. It writes a small override card, runs `hiddenvalley`, and names files using the scanned values.

## Output Files

Each run produces three data products, and one log file, with the same stem.

### `.eventdat`

Event-level missing-momentum information, and event metadata.

### `.jetdat`

Jet-level metadata. This includes jets with zero Lund declusterings, that do not get captures in the lunddat file.

### `.lunddat`

Cluster-level Lund-plane data. 

## Lund Density

The plotting scripts use averaged Lund density from eqn 2.3 of http://arxiv.org/abs/1807.04758 :
```text
rho(Delta, kt) = (1 / N_jet) d n_emission / (d ln kt d ln(1/Delta))
```

In finite bins this is implemented as:
```text
rho_bin = N_bin / (N_selected_jets * bin_area)
```

This is not a probability and does not have to lie between 0 and 1. It is the average number of declusterings per selected jet per unit Lund-plane area.

Jets with `nLund = 0` contribute zero emissions but are still counted in `N_selected_jets` when they pass the requested jet and DM masks, which is why `.jetdat` is needed.

## Plotting One Sample

Use:

```bash
python3 plot_lund_options.py sample.lunddat
```

Jet selection:
```bash
--jets all
--jets hardest
--jets pt-min --jet-pt-min 50
```

If `--jets pt-min` is used without `--jet-pt-min`, the threshold is read from the matching `.eventdat` header.

DM selection:
```bash
--dm-selection all
--dm-selection dm --dm-mask 0.75
--dm-selection nondm --dm-mask 0.75
```

The script reads:
```text
sample.lunddat
sample.jetdat
sample.eventdat
```

The plot title and output filename include non-default jet and DM choices.

## Plotting A Difference

Use:

```bash
python3 plot_lund_difference.py sample1.lunddat sample2.lunddat
```

The script independently builds a Lund density for each sample and plots:

```text
rho(sample2) - rho(sample1)
```

Independent jet masks:
```bash
--jets1 all|hardest|pt-min
--jet-pt-min1 50
--jets2 all|hardest|pt-min
--jet-pt-min2 50
```

Independent DM masks:
```bash
--mask1 none|dm|nondm
--dm-mask1 0.75
--mask2 none|dm|nondm
--dm-mask2 0.75
```

Filename and title mask labels are only added when a mask differs from the default.

## Notes

- Old `.lunddat` files made before the `.jetdat` split are not compatible with the current plotting scripts.
- Do not silently replace missing hardest jets with jet 1 in plotting. A hardest jet with `nLund = 0` is a real selected jet with no resolved declustering
  entries, and should stay in the jet-level denominator when selected.
- The `oldstuff` folder contains a few old scripts and notes that are not part of the current workflow, but may be useful for reference, especially given the streamlined `hiddenvalley.cc` interface.