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
#include <sstream>
#include <string>
#include <cctype>
#include <stdexcept>
#include <cmath>
#include <vector>
#include <iomanip>
#include <algorithm>

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

struct RunOptions
{
  int nEvents = 10000;
  int seed = -1;
  std::string outfile = "hv";
  std::vector<std::string> cmndfiles;
};

RunOptions parseCommandLine(int argc, char *argv[])
{
  RunOptions options;

  int opt;
  while ((opt = getopt(argc, argv, "e:o:s:")) != -1)
  {
    switch (opt)
    {
    case 'e':
      options.nEvents = parseEvents(optarg);
      break;
    case 'o':
      options.outfile = optarg;
      break;
    case 's':
      options.seed = std::stoi(optarg);
      break;
    default:
      std::cerr << "Usage: " << argv[0]
                << " [-e events] [-o outfiles label] [-s seed]"
                << " [cmndfiles...]\n";
      throw std::runtime_error("Invalid command-line arguments");
    }
  }

  if (optind < argc)
    options.cmndfiles = std::vector<std::string>(argv + optind, argv + argc);

  return options;
}

//==========================================================================
using namespace Pythia8;

bool fromDarkHadron(int i, const Event& event)
{
  std::vector<int> stack(1, i);
  std::vector<bool> seen(event.size(), false);
  while (!stack.empty()) {
    int iNow = stack.back();
    stack.pop_back();
    if (iNow <= 0 || iNow >= event.size() || seen[iNow])
      continue;
    seen[iNow] = true;
    int id = event[iNow].idAbs();
    if (id >= 4900111 && id <= 4900213)
      return true;
    int mother1 = event[iNow].mother1();
    int mother2 = event[iNow].mother2();
    if (mother1 > 0)
      stack.push_back(mother1);
    if (mother2 > 0 && mother2 != mother1)
      stack.push_back(mother2);
  }
  return false;
}

void darkHadronFractions(const fastjet::PseudoJet& jet, const Event& event,
                         double& eFrac, double& ptFrac)
{
  std::vector<fastjet::PseudoJet> constituents = jet.constituents();
  if (constituents.size() == 0) {
    eFrac = 0.;
    ptFrac = 0.;
    return;
  }
  double eDark = 0.;
  double ptDark = 0.;
  double eConstituents = 0.;
  double ptConstituents = 0.;
  for (int iConst = 0; iConst < int(constituents.size()); ++iConst){
    eConstituents += constituents[iConst].E();
    ptConstituents += constituents[iConst].pt();
    if (fromDarkHadron(constituents[iConst].user_index(), event)) {
      eDark += constituents[iConst].E();
      ptDark += constituents[iConst].pt();
    }
  }
  eFrac = eDark / eConstituents;
  ptFrac = ptDark / ptConstituents;
}

//==========================================================================

int main(int argc, char *argv[])
{
  // Parse command line options.
  RunOptions options;
  try {
    options = parseCommandLine(argc, argv);
  } catch (const std::runtime_error&) {
    return 1;
  }
  int nEvents = options.nEvents;
  int seed = options.seed;
  std::string outfile = options.outfile;
  std::vector<std::string> cmndfiles = options.cmndfiles;

  // Suppress Pythia setup output.
  stringstream coutBuf;
  streambuf *oldCout = cout.rdbuf(coutBuf.rdbuf());

  //Generator.
  Pythia pythia;
  pythia.settings.addFlag("Main:writeLog", false);

  pythia.settings.addParm("JetAnalysis:etaMax", 3.6, false, false, 0., 10.); 
  pythia.settings.addParm("JetAnalysis:jetPtMin", 20.0, false, false, 0., 1000.);

  pythia.settings.addFlag("HiddenValley:setMesonMassesFromQv", false);
  pythia.settings.addParm("HiddenValley:ratioPiRho", 2.5, true, false, 0., 0.);

  pythia.settings.addFlag("HiddenValley:useCouplings", false);
  pythia.settings.addParm("HiddenValley:gSM", 1, false, false, 0., 0.);
  pythia.settings.addParm("HiddenValley:gHV", 2.4, false, false, 0., 0.);

  for (int iCmnd = 0; iCmnd < (int)cmndfiles.size(); ++iCmnd)
    if (!cmndfiles[iCmnd].empty())
      pythia.readFile(cmndfiles[iCmnd]);
  if (seed > 0) {
    pythia.readString("Random:setSeed = on");
    pythia.readString("Random:seed = " + std::to_string(seed));
  }

  // Detector-level parameters for jet analysis.
  double etaMax = pythia.settings.parm("JetAnalysis:etaMax");
  double jetPtMin = pythia.settings.parm("JetAnalysis:jetPtMin");

  // other variables for hidden valley showering
  if (pythia.settings.flag("HiddenValley:useCouplings"))
  {
    pythia.readString("HiddenValley:gSM = " + std::to_string(pythia.settings.parm("HiddenValley:gSM")));
    pythia.readString("HiddenValley:gHV = " + std::to_string(pythia.settings.parm("HiddenValley:gHV")));
    int nGauge = pythia.settings.mode("HiddenValley:Ngauge");
    int nFlav = pythia.settings.mode("HiddenValley:nFlav");
    double mZp = pythia.particleData.m0(4900023);
    double SMwidth = (3.0 * mZp * pythia.settings.parm("HiddenValley:gSM") * pythia.settings.parm("HiddenValley:gSM")) / (2.0 * 3.14159);
    double HVwidth = (nGauge * nFlav * mZp * pythia.settings.parm("HiddenValley:gHV") * pythia.settings.parm("HiddenValley:gHV")) / (12.0 * 3.14159);
    double Zwidth = SMwidth + HVwidth;
    pythia.readString("4900023:mWidth = " + std::to_string(Zwidth));
  }
  double Zwidth = pythia.particleData.mWidth(4900023);
  double Lambda = pythia.settings.parm("HiddenValley:Lambda");
  double mq = pythia.particleData.m0(4900101);

  if (pythia.settings.flag("HiddenValley:setMesonMassesFromQv") && mq > 0.) {
    double ratioPiRho = pythia.settings.parm("HiddenValley:ratioPiRho");
    double mPi = 8.0 * mq / (3.0 * ratioPiRho + 1.0);
    double mRho = ratioPiRho * mPi;
    pythia.readString("4900111:m0 = " + std::to_string(mPi));
    pythia.readString("4900113:m0 = " + std::to_string(mRho));
  }

  // Logfile initialization.
  ofstream logBuf;
  if (pythia.settings.flag("Main:writeLog"))
  {
    logBuf.open(outfile + ".log");
    cout.rdbuf(logBuf.rdbuf());
  }
  else
    cout.rdbuf(oldCout);
  cout << coutBuf.str();

  // If Pythia fails to initialize, exit with error.
  if (!pythia.init())
  {
    cout.rdbuf(oldCout);
    std::cerr << coutBuf.str();
    return 1;
  }

  // Event/process shorthand.
  Event &event = pythia.event;

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

  // lund plane output
  std::ofstream lundOut(outfile + ".lunddat");
  int precision = 6;
  int fracPrecision = 2;
  int lundWidth = std::max(precision + 8, 12);
  lundOut << "#evt"
          << std::setw(4) << "jet"
          << std::setw(lundWidth-4) << "clustE"
          << std::setw(lundWidth-3) << "clustPt"
          << std::setw(lundWidth) << "ln_1/Delta"
          << std::setw(lundWidth) << "ln_kt"
          << std::setw(lundWidth) << "z"
          << std::setw(lundWidth) << "psi"
          << std::setw(lundWidth) << "m2"
          << std::setw(5) << "mult" << "\n";
  lundOut << std::fixed << std::setprecision(precision);

  // jet level output
  std::ofstream jetOut(outfile + ".jetdat");
  jetOut << "#evt"
         << std::setw(4) << "jet"
         << std::setw(lundWidth-1) << "pt"
         << std::setw(lundWidth) << "eta"
         << std::setw(lundWidth) << "darkEfrac"
         << std::setw(lundWidth) << "darkPtfrac"
         << std::setw(7) << "nLund" << "\n";
  jetOut << std::fixed << std::setprecision(precision);

  // event level output file
  std::ofstream eventOut(outfile + ".eventdat");
  int evtprecision = 6;
  int evtWidth = std::max(evtprecision + 8, 12);
  eventOut << "#etacut jetPtMin Z'width Lambda m_qv \n";
  eventOut << "#  " << etaMax << "    " << jetPtMin << "        " << Zwidth << "     " << Lambda << "     " << mq << "\n";
  eventOut << "#evt"
          << std::setw(evtWidth-3) << "metPx"
          << std::setw(evtWidth) << "metPy"
          << std::setw(evtWidth) << "metDet"
          << std::setw(evtWidth) << "metTruth"
          << "\n";
  eventOut << std::fixed << std::setprecision(evtprecision);

  int iErr = 0;

  // Begin event loop. Generate event.
  for (int iEvent = 0; iEvent < nEvents; ++iEvent) {
    if (!pythia.next()) {
      if (++iErr < 100) continue;
      else break;
    }

    // Keep track of visible momentum.
    Vec4 visibleDetMomvec;
    Vec4 invisibleTruthMomvec;
    fjInputs.clear();

    // Loop over event record to decide what to pass to FastJet.
    for (int i = 0; i < pythia.event.size(); ++i)
    {
      // Final state only.
      if (!pythia.event[i].isFinal())
        continue;

      if (!pythia.event[i].isVisible())
        invisibleTruthMomvec += pythia.event[i].p();

      // Detector-level jets should only use stable visible particles.
      if (!pythia.event[i].isVisible())
        continue;

      // Only |eta| < 3.6.
      if (std::abs(pythia.event[i].eta()) > etaMax)
        continue;

      // Visible momentum entering the jet definition 
      visibleDetMomvec += pythia.event[i].p();

      // Store as input to Fastjet.
      fastjet::PseudoJet particle(pythia.event[i].px(), pythia.event[i].py(),
                                  pythia.event[i].pz(), pythia.event[i].e());
      particle.set_user_index(i);
      fjInputs.push_back(particle);
    }

    // Check that event contains analyzable particles.
    if (fjInputs.size() == 0)
      continue;

    // Run Fastjet algorithm.
    vector<fastjet::PseudoJet> inclusiveJets, sortedJets;
    fastjet::ClusterSequence clustSeq(fjInputs, *jetDef);

    // Extract inclusive jets sorted by pT (note minimum pT of 20.0 GeV).
    inclusiveJets = clustSeq.inclusive_jets(jetPtMin);
    sortedJets = sorted_by_pt(inclusiveJets);

    if (sortedJets.size() < 1)
      continue;

    double metPx = -visibleDetMomvec.px();
    double metPy = -visibleDetMomvec.py();
    double met = std::sqrt(metPx * metPx + metPy * metPy);
    eventOut << iEvent
             << std::setw(evtWidth) << metPx
             << std::setw(evtWidth) << metPy
             << std::setw(evtWidth) << met
             << std::setw(evtWidth) << invisibleTruthMomvec.pT()
             << "\n";

    for (int iJet = 0; iJet < int(sortedJets.size()); ++iJet) {
      double darkHadronEnergyFrac = 0.;
      double darkHadronPtFrac = 0.;
      darkHadronFractions(sortedJets[iJet], event, darkHadronEnergyFrac,
                          darkHadronPtFrac);

      fastjet::PseudoJet j = sortedJets[iJet];
      fastjet::PseudoJet j1, j2;
      int nLundEntries = 0;

      while (j.has_parents(j1, j2))
      {
        if (j1.pt() < j2.pt())
          std::swap(j1, j2);
        // In declustering language: j is the pre-branching object,
        // j1 is the harder post-branching branch, j2 is the softer emitted branch.
        double clusterE = j.e();
        double clusterPt = j.pt();
        double Delta = j1.delta_R(j2);
        double kt = j2.pt() * Delta;
        if (Delta <= 0. || kt <= 0. || !std::isfinite(Delta) || !std::isfinite(kt))
          break;
        double ln_oneoverDelta = std::log(1.0 / Delta);
        double ln_kt = std::log(kt);
        double z = j2.pt() / (j1.pt() + j2.pt());
        double psi = std::atan2(j1.rap() - j2.rap(), j1.delta_phi_to(j2));
        double m2 = j.m2();
        int multiplicity = j.constituents().size();

        lundOut << iEvent
                << std::setw(6) << iJet
                << std::setprecision(fracPrecision)
                << std::setw(lundWidth-3) << clusterE
                << std::setw(lundWidth-3) << clusterPt
                << std::setprecision(evtprecision)
                << std::setw(lundWidth) << ln_oneoverDelta
                << std::setw(lundWidth) << ln_kt
                << std::setw(lundWidth) << z
                << std::setw(lundWidth) << psi
                << std::setw(lundWidth) << m2
                << std::setw(5) << multiplicity 
                << "\n";

        ++nLundEntries;
        j = j1;
      }

      jetOut << iEvent
             << std::setw(6) << iJet
             << std::setw(lundWidth) << sortedJets[iJet].pt()
             << std::setw(lundWidth) << sortedJets[iJet].eta()
             << std::setprecision(fracPrecision)
             << std::setw(lundWidth) << darkHadronEnergyFrac
             << std::setw(lundWidth) << darkHadronPtFrac
             << std::setprecision(precision)
             << std::setw(7) << nLundEntries
             << "\n";
    }

    // End of event loop.
    }
  // Done.
  return 0;
}
