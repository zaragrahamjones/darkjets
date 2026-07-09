// darkjetplane.cc is a part of the PYTHIA event generator.
// Copyright (C) 2025 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Authors: Nishita Desai <nishita.desai@tifr.res.in>

// Keywords: jet finding; fastjet; BSM; dark matter

// This is a simple test program to study jets in Dark Matter production.

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

int parseEvents(const std::string& s) {
  char suffix = s.back();
  long multiplier = 1;
  std::string number = s;
  if (std::isalpha(static_cast<unsigned char>(suffix))) {
      number.pop_back();

      switch (std::tolower(static_cast<unsigned char>(suffix))) {
      case 'k': multiplier = 1000; break;
      case 'm': multiplier = 1000000; break;
      case 'g': multiplier = 1000000000; break;
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

  // Generator. Process selection.
  Pythia pythia;

  int nEvents = 1000;
  std::string outfile = "lund.dat";
  std::string cmndfile = "darkjetplane.cmnd";

  int opt;
  while ((opt = getopt(argc, argv, "e:r:")) != -1)
  {
    switch (opt)
    {
    case 'e':
      nEvents = parseEvents(optarg);
      break;
    case 'r':
      outfile = optarg;
      break;
    default:
      std::cerr << "Usage: " << argv[0]
                << " [-e events] [-r outfile] [cmndfile]\n";
      return 1;
    }
  }
  if (optind < argc)
    cmndfile = argv[optind];
  pythia.readFile(cmndfile);
  int numberCount = std::max(1, nEvents / 100);
  pythia.readString("Next:numberCount = " + std::to_string(numberCount));
  std::ofstream out(outfile);

  // If Pythia fails to initialize, exit with error.
  if (!pythia.init()) return 1;

  // Event (process) shorthand.
  Event& process = pythia.process;

  // Fastjet analysis - select algorithm and parameters.
  double Rparam = 0.4;
  fastjet::Strategy               strategy = fastjet::Best;
  fastjet::RecombinationScheme    recombScheme = fastjet::E_scheme;
  fastjet::JetDefinition         *jetDef = NULL;
  // jetDef = new fastjet::JetDefinition( fastjet::kt_algorithm, Rparam,
  //          recombScheme, strategy);
  jetDef = new fastjet::JetDefinition(fastjet::cambridge_algorithm, Rparam,
                                      recombScheme, strategy);

  // Fastjet input.
  std::vector <fastjet::PseudoJet> fjInputs;

  // Histograms. Error counter.
  Hist pTj("dN/dpTj", 100, 0., 100.);
  Hist mRec("mRec", 100, 0., 1000.);
  int iErr = 0;

  std::ofstream lundOut(outfile);
  lundOut << "# event jet log1overDelta logkt z Delta kt\n";

  // print all settings
  pythia.settings.listAll();
  // pythia.particleData.listAll();
  pythia.particleData.list(55);
  pythia.particleData.list(52);

  // Begin event loop. Generate event. Skip if error.
  for (int iEvent = 0; iEvent < nEvents; ++iEvent) {
    if (!pythia.next()) {
      if (++iErr < 100) continue;
      else {
        cout << "Too many errors" << endl;
        break;
      }
    }

    // Invariant mass of DM system.
    Vec4 mRes = process[5].p() + process[6].p();
    mRec.fill(mRes.mCalc());

    // Keep track of missing ET.
    Vec4 missingETvec;
    fjInputs.clear();

    // Loop over event record to decide what to pass to FastJet.
    for (int i = 0; i < pythia.event.size(); ++i) {
      // Final state only.
      if (!pythia.event[i].isFinal()) continue;

      // No neutrinos or DM.
      if ( pythia.event[i].idAbs() == 12 || pythia.event[i].idAbs() == 14
        || pythia.event[i].idAbs() == 16 || pythia.event[i].idAbs() == 52)
        continue;

      // Only |eta| < 3.6.
      if (abs(pythia.event[i].eta()) > 3.6) continue;

      // Missing ET.
      missingETvec += pythia.event[i].p();

      // Store as input to Fastjet.
      fjInputs.push_back( fastjet::PseudoJet( pythia.event[i].px(),
        pythia.event[i].py(), pythia.event[i].pz(), pythia.event[i].e() ) );
    }

    // Check that event contains analyzable particles.
    if (fjInputs.size() == 0) {
      cout << "Error: event with no final state particles" << endl;
      continue;
    }

    // Run Fastjet algorithm.
    vector <fastjet::PseudoJet> inclusiveJets, sortedJets;
    fastjet::ClusterSequence clustSeq(fjInputs, *jetDef);

    // Extract inclusive jets sorted by pT (note minimum pT of 20.0 GeV).
    inclusiveJets = clustSeq.inclusive_jets(20.0);
    sortedJets    = sorted_by_pt(inclusiveJets);

    if(sortedJets.size() < 1) {
      // cout << " No jets found in event " << iEvent << endl;
      continue;
    }
    pTj.fill( sortedJets[0].pt() );

    // Loop over the hardest jet and extract Lund plane information.
    fastjet::PseudoJet j = sortedJets[0];
    fastjet::PseudoJet j1, j2;

    while (j.has_parents(j1, j2))
    {
      if (j1.pt() < j2.pt())
        std::swap(j1, j2);

      double Delta = j1.delta_R(j2);
      if (Delta <= 0.)
        break;

      double z = j2.pt() / (j1.pt() + j2.pt());
      double kt = j2.pt() * Delta;

      lundOut << iEvent << " 0 "
              << log(1. / Delta) << " "
              << log(kt) << " "
              << z << " "
              << Delta << " "
              << kt << "\n";

      j = j1;
    }

  // End of event loop. Statistics. Histogram.
  }
  pythia.stat();
  cout << pTj;

  // Done.
  return 0;
}
