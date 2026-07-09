===================================
INITIALIZATION:
===================================

darkjetplane is based on main507 in the PYTHIA examples, with the addition of a lund.dat output file for the dark jet analysis. The .cmnd file is used to set up the Pythia parameters for the simulation. The .cc file contains the main program that initializes Pythia, runs the event loop, and performs jet clustering using FastJet. The Makefiles are directly copied from the example folder too and Makefile.inc has been updated to include darkjetplane as a target for compilation.

plot_lund.py is a simple script to plot the lund.dat output file. It uses matplotlib to create a scatter plot of the log(1/Delta) vs log(kt) for the jets in the events.

===================================
COMPARISON WITH SM MODELS:
===================================

Within the .cmnd file, Zp mostly changes the rate and width, not a clean "amount of DM".

A clean SM comparison model would be something that gives the same _visible_ final state: a QCD jet recoiling against missing energy (neutrinos, anything else?). 
This could be achieved by:
- a Z' with the same mass and width, but decaying to SM particles instead of DM,
- pp-> Z -> QCD jets + invisible particles, and adjusting Z'=Z
