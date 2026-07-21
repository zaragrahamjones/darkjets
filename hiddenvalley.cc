// hiddenvalley.cc is a part of the PYTHIA event generator.
// Copyright (C) 2025 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Authors: Torbjorn Sjostrand <torbjorn.sjostrand@fysik.lu.se>

// Keywords: Hidden Valley;

// Test of Hidden Valley production in a few different channels.

#include "Pythia8/Pythia.h"
#include "Pythia8Plugins/FastJet3.h"
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"

#include <unistd.h>
#include <fstream>
#include <string>
#include <cctype>
#include <stdexcept>
#include <cmath>

int parseEvents(const std::string &s)
{
  char suffix = s.back();
  long multiplier = 1;
  std::string number = s;
  if (std::isalpha(static_cast<unsigned char>(suffix)))
  {
    number.pop_back();

    switch (std::tolower(static_cast<unsigned char>(suffix)))
    {
    case 'k':
      multiplier = 1000;
      break;
    case 'm':
      multiplier = 1000000;
      break;
    case 'g':
      multiplier = 1000000000;
      break;
    default:
      throw std::runtime_error("Unknown suffix: " + std::string(1, suffix));
    }
  }
  return std::stol(number) * multiplier;
}

//==========================================================================

using namespace Pythia8;

int main(int argc, char *argv[])
{
  //Generator.
  Pythia pythia;

  // Generator. Process selection from command-line arguments.
  int nEvents = 1000;
  int seed = -1;
  std::string outfile = "hv_lund.dat";
  std::string cmndfile = ""; // no default command file, must be specified

  int opt;
  while ((opt = getopt(argc, argv, "e:r:s:")) != -1)
  {
    switch (opt)
    {
    case 'e':
      nEvents = parseEvents(optarg);
      break;
    case 'r':
      outfile = optarg;
      break;
    case 's':
      seed = std::stoi(optarg);
      break;
    default:
      std::cerr << "Usage: " << argv[0]
                << " [-e events] [-r outfile] [-s seed] [cmndfile]\n";
      return 1;
    }
  }
  if (optind < argc)
    cmndfile = argv[optind];
  pythia.readFile(cmndfile);
  if (seed > 0) {
    pythia.readString("Random:setSeed = on");
    pythia.readString("Random:seed = " + std::to_string(seed));
  }
  int numberCount = std::max(1, nEvents / 10);
  pythia.readString("Next:numberCount = " + std::to_string(numberCount));

  // Production process via H, and decay to gv gv or gammav gammav. 
  // .cmnd cannot handle this, so we do it here.
  if (pythia.settings.flag("25:onMode") == 0)
  {
    if (pythia.settings.flag("Hiddenvalley:alphaOrder") == 1)
      pythia.readString("25:addChannel = 1 0.1 100 4900021 4900021");
    else
      pythia.readString("25:addChannel = 1 0.1 100 4900022 4900022");
  }

  // If Pythia fails to initialize, exit with error.
  if (!pythia.init()) return 1;

  // Event/process shorthand.
  Event &event = pythia.event;
  Event &process = pythia.process;

  // Fastjet analysis - select algorithm and parameters.
  double Rparam = 0.4;
  fastjet::Strategy strategy = fastjet::Best;
  fastjet::RecombinationScheme recombScheme = fastjet::E_scheme;
  fastjet::JetDefinition *jetDef = NULL;
  // jetDef = new fastjet::JetDefinition( fastjet::kt_algorithm, Rparam,
  //          recombScheme, strategy);
  jetDef = new fastjet::JetDefinition(fastjet::cambridge_algorithm, Rparam,
                                      recombScheme, strategy);

  // Fastjet input.
  std::vector<fastjet::PseudoJet> fjInputs;

  std::ofstream lundOut(outfile);
  lundOut << "# event jet Delta kt z psi m2 multiplicity \n";

  // Book histograms. also error counter.
  Hist nGluonv( "number of HV gluons",  100, -0.5, 99.5);
  Hist nGammav( "number of HV gammas",  100, -0.5, 99.5);
  Hist nHadronv("number of HV hadrons", 100, -0.5, 99.5);
  Hist pTj("dN/dpTj", 100, 0., 100.);
  Hist mRec("mRec", 100, 0., 1000.);
  int iErr = 0;

  // print all settings
  pythia.settings.listChanged();

  std::vector<int> hvIds = {
      4900001, 4900101, 4900102, 4900103,
      4900021, 4900022, 4900023,
      4900111, 4900113, 4900211, 4900213};
  pythia.particleData.list(hvIds);

  // Begin event loop. Generate event. Extra HV colour output.
  for (int iEvent = 0; iEvent < nEvents; ++iEvent) {
    if (!pythia.next()) {
      if (++iErr < 100) continue;
      else {
        cout << "Too many errors" << endl;
        break;
      }
    }
    if (iEvent == 0 && (process.hasHVcols() || event.hasHVcols())) {
      process.listHVcols();
      event.listHVcols();
    }

    // Number of "final" gv, gammav and hadronv in current event.
    int nGluonvNow  = 0;
    int nGammavNow  = 0;
    int nHadronvNow = 0;
    for (int i = 0; i < event.size(); ++i) {
      int idNow = event[i].idAbs();
      int idDau = event[ event[i].daughter1() ].idAbs();
      if      (idNow == 4900021 && idDau != 4900021) ++nGluonvNow;
      else if (idNow == 4900022 && idDau != 4900022) ++nGammavNow;
      else if (idNow >  4900110) ++nHadronvNow;
    }
    nGluonv.fill( nGluonvNow);
    nGammav.fill( nGammavNow);
    nHadronv.fill( nHadronvNow);

    // Invariant mass of DM system.
    Vec4 mRes = process[5].p() + process[6].p();
    mRec.fill(mRes.mCalc());

    // Keep track of missing ET.
    Vec4 missingETvec;
    fjInputs.clear();

    // Loop over event record to decide what to pass to FastJet.
    for (int i = 0; i < pythia.event.size(); ++i)
    {
      // Final state only.
      if (!pythia.event[i].isFinal())
        continue;

      // Detector-level jets should only use stable visible particles.
      if (!pythia.event[i].isVisible())
        continue;

      // Only |eta| < 3.6.
      if (std::abs(pythia.event[i].eta()) > 3.6)
        continue;

      // Visible momentum entering the jet definition.
      missingETvec += pythia.event[i].p();

      // Store as input to Fastjet.
      fjInputs.push_back(fastjet::PseudoJet(pythia.event[i].px(),
                                            pythia.event[i].py(), pythia.event[i].pz(), pythia.event[i].e()));
    }

    // Check that event contains analyzable particles.
    if (fjInputs.size() == 0)
    {
      cout << "Error: event with no final state particles" << endl;
      continue;
    }

    // Run Fastjet algorithm.
    vector<fastjet::PseudoJet> inclusiveJets, sortedJets;
    fastjet::ClusterSequence clustSeq(fjInputs, *jetDef);

    // Extract inclusive jets sorted by pT (note minimum pT of 20.0 GeV).
    inclusiveJets = clustSeq.inclusive_jets(20.0);
    sortedJets = sorted_by_pt(inclusiveJets);

    if (sortedJets.size() < 1)
    {
      cout << " No jets found in event " << iEvent << endl;
      continue;
    }
    // Extract Lund plane information for every selected jet. sortedJets is
    // ordered by pT, so iJet = 0 is the leading jet, iJet = 1 is subleading, etc.
    for (int iJet = 0; iJet < int(sortedJets.size()); ++iJet) {
      pTj.fill(sortedJets[iJet].pt());

      fastjet::PseudoJet j = sortedJets[iJet];
      fastjet::PseudoJet j1, j2;

      while (j.has_parents(j1, j2))
      {
        if (j1.pt() < j2.pt())
          std::swap(j1, j2);
        // In declustering language: j is the pre-branching object,
        // j1 is the harder post-branching branch, j2 is the softer emitted branch.

        double Delta = j1.delta_R(j2);
        if (Delta <= 0.)
          break;

        double z = j2.pt() / (j1.pt() + j2.pt());
        double kt = j2.pt() * Delta;
        double psi = std::atan2(j1.rap() - j2.rap(), j1.delta_phi_to(j2));
        double m2 = j.m2();
        double multiplicity = j1.constituents().size() + j2.constituents().size();

        lundOut << iEvent << " " << iJet << " "
                << Delta << " "
                << kt << " "
                << z << " "
                << psi << " "
                << m2 << " "
                << multiplicity << "\n";

        j = j1;
      }
    }

    // End of event loop. Print statistics and histograms.
    }
  pythia.stat();
  cout << nGluonv << nGammav << nHadronv << pTj;

  // Done.
  return 0;
}
