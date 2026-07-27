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
#include "Pythia8Plugins/InputParser.h"

#include <chrono>
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
#include <map>

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
  int nEvents = 1000;
  int seed = -1;
  int eventRecordEvent = -1;
  std::string outfile = "hv";
  std::vector<std::string> cmndfiles;
};

RunOptions parseCommandLine(int argc, char *argv[])
{
  RunOptions options;

  int opt;
  while ((opt = getopt(argc, argv, "e:o:s:r:")) != -1)
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
    case 'r':
      options.eventRecordEvent = std::stoi(optarg);
      break;
    default:
      std::cerr << "Usage: " << argv[0]
                << " [-e events] [-o outfiles label] [-s seed]"
                << " [-r event# record] [cmndfiles...]\n";
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

bool fromDarkParton(int i, const Event &event)
{
  std::vector<int> stack(1, i);
  std::vector<bool> seen(event.size(), false);
  while (!stack.empty())
  {
    int iNow = stack.back();
    stack.pop_back();
    if (iNow <= 0 || iNow >= event.size() || seen[iNow])
      continue;
    seen[iNow] = true;
    int id = event[iNow].idAbs();
    if (id >= 4900101 && id <= 4900108 || id == 4900021)
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

std::vector<double> darkHadronFractions(const fastjet::PseudoJet& jet,
                                  const Event& event)
{
  std::vector<fastjet::PseudoJet> constituents = jet.constituents();
  if (constituents.size() == 0)
    return {0., 0., 0.};
  int nDark = 0;
  double eDark = 0.;
  double ptDark = 0.;
  double eConstituents = 0.;
  double ptConstituents = 0.;
  for (int iConst = 0; iConst < int(constituents.size()); ++iConst){
    eConstituents += constituents[iConst].E();
    ptConstituents += constituents[iConst].pt();
    if (fromDarkHadron(constituents[iConst].user_index(), event)) {
      ++nDark;
      eDark += constituents[iConst].E();
      ptDark += constituents[iConst].pt();
    }
  }
  return {double(nDark) / double(constituents.size()), eDark / eConstituents, ptDark / ptConstituents};
}

void printEventBranch(int i, const Event& event, std::ostream& out,
                      const std::map<int, int>& jetIndexOfParticle, int depth,
                      const std::string& branch)
{
  for (int iDepth = 0; iDepth < depth; ++iDepth)
    out << "  ";

  if (i <= 0 || i >= event.size()) {
    out << branch << " unknown\n";
    return;
  }

  int daughter1 = event[i].daughter1();
  int daughter2 = event[i].daughter2();
  out << branch << " " << event[i].name();
  if (daughter1 <= 0) {
    std::map<int, int>::const_iterator it = jetIndexOfParticle.find(i);
    if (it != jetIndexOfParticle.end())
      out << " [jet " << it->second << "]";
    else
      out << " [no jet]";
    out << "\n";
    return;
  }
  out << "\n";

  std::vector<int> daughters;
  int lastDaughter = (daughter2 > daughter1) ? daughter2 : daughter1;
  for (int iDaughter = daughter1; iDaughter <= lastDaughter; ++iDaughter)
    if (iDaughter > 0 && iDaughter < event.size())
      daughters.push_back(iDaughter);

  std::sort(daughters.begin(), daughters.end(),
            [&event](int a, int b) { return event[a].pT() > event[b].pT(); });

  for (int iDaughter = 0; iDaughter < int(daughters.size()); ++iDaughter) {
    std::string daughterBranch = (iDaughter == 0) ? "hard" : "soft";
    printEventBranch(daughters[iDaughter], event, out, jetIndexOfParticle,
                     depth + 1, daughterBranch);
  }
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
  int eventRecordEvent = options.eventRecordEvent;
  std::string outfile = options.outfile;
  std::vector<std::string> cmndfiles = options.cmndfiles;

  // Catch all Pythia output until the log setting has been read.
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
  int numberCount = std::max(1, nEvents / 10);
  pythia.readString("Next:numberCount = " + std::to_string(numberCount));

  // Detector-level parameters for jet analysis.
  double etaMax = pythia.settings.parm("JetAnalysis:etaMax");
  double jetPtMin = pythia.settings.parm("JetAnalysis:jetPtMin");

  // Coupling constants input to width for hidden valley showering
  if (pythia.settings.flag("HiddenValley:useCouplings"))
  {
    pythia.readString("HiddenValley:gSM = " + std::to_string(pythia.settings.parm("HiddenValley:gSM")));
    pythia.readString("HiddenValley:gHV = " + std::to_string(pythia.settings.parm("HiddenValley:gHV")));
    int nGauge = pythia.settings.mode("HiddenValley:Ngauge");
    int nFlav = pythia.settings.mode("HiddenValley:nFlav");
    double mZp = pythia.particleData.m0(4900023);
    double SMwidth = (3.0 *mZp* pythia.settings.parm("HiddenValley:gSM") * pythia.settings.parm("HiddenValley:gSM")) / (2.0 * 3.14159);
    double HVwidth = (nGauge * nFlav *mZp* pythia.settings.parm("HiddenValley:gHV") * pythia.settings.parm("HiddenValley:gHV")) / (12.0 * 3.14159);
    double width = SMwidth + HVwidth;
    pythia.readString("4900023:mWidth = " + std::to_string(width));
  }

  // Production process via H, and decay to gv gv or gammav gammav. 
  // .cmnd cannot handle this, so we do it here.
  if (pythia.settings.flag("HiggsSM:all"))
  {
    if (pythia.settings.mode("Hiddenvalley:alphaOrder") == 1)
      pythia.readString("25:addChannel = 1 0.1 100 4900021 4900021");
    else
      pythia.readString("25:addChannel = 1 0.1 100 4900022 4900022");
  }

  double mQv = pythia.particleData.m0(4900101);
  if (pythia.settings.flag("HiddenValley:setMesonMassesFromQv") && mQv > 0.) {
    double ratioPiRho = pythia.settings.parm("HiddenValley:ratioPiRho");
    double mPi = 8.0 / (3.0 * ratioPiRho + 1.0) * mQv;
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

  std::ofstream lundOut(outfile + ".dat");
  int precision = 6;
  int fracPrecision = 2;
  int realWidth = std::max(precision + 8, 12);
  lundOut << "#evt"
          << std::setw(4) << "jet"
          << std::setw(realWidth-3) << "darkEfrac"
          << std::setw(realWidth-4) << "clustE"
          << std::setw(realWidth-3) << "darkPtfrac"
          << std::setw(realWidth-3) << "clustPt"
          << std::setw(realWidth) << "ln_1/Delta"
          << std::setw(realWidth) << "ln_kt"
          << std::setw(realWidth) << "z"
          << std::setw(realWidth) << "psi"
          << std::setw(realWidth) << "m2"
          << std::setw(5) << "mult" << "\n";
  lundOut << std::fixed << std::setprecision(precision);

  // Book histograms. also error counter.
  Hist nGluonv( "number of HV gluons",  100, -0.5, 99.5);
  Hist nGammav( "number of HV gammas",  100, -0.5, 99.5);
  Hist nHadronv("number of HV hadrons", 100, -0.5, 99.5);
  Hist pTj("dN/dpTj", 100, 0., 100.);
  Hist mRec("mRec", 100, 0., 1000.);
  int iErr = 0;

  // print all settings
  pythia.settings.listChanged();

  // for nflav <= 3
  std::vector<int> hvConnectors = {4900001, 4900002, 4900003, 4900023};
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
      if (std::abs(pythia.event[i].eta()) > etaMax)
        continue;

      // Visible momentum entering the jet definition.
      missingETvec += pythia.event[i].p();

      // Store as input to Fastjet.
      fastjet::PseudoJet particle(pythia.event[i].px(), pythia.event[i].py(),
                                  pythia.event[i].pz(), pythia.event[i].e());
      particle.set_user_index(i);
      fjInputs.push_back(particle);
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
    inclusiveJets = clustSeq.inclusive_jets(jetPtMin);
    sortedJets = sorted_by_pt(inclusiveJets);

    if (sortedJets.size() < 1)
    {
      cout << " No jets found in event " << iEvent << endl;
      continue;
    }
    // Extract Lund plane information for every selected jet. sortedJets is
    // ordered by pT, so iJet = 0 is the leading jet, iJet = 1 is subleading, etc.
    if (eventRecordEvent >= 0 && iEvent == eventRecordEvent) {
      std::ofstream diagramOut(outfile + "_event"
                               + std::to_string(eventRecordEvent)
                               + "_eventrecord.dat");
      std::map<int, int> jetIndexOfParticle;
      for (int iJet = 0; iJet < int(sortedJets.size()); ++iJet) {
        std::vector<fastjet::PseudoJet> constituents = sortedJets[iJet].constituents();
        for (int iConst = 0; iConst < int(constituents.size()); ++iConst)
          jetIndexOfParticle[constituents[iConst].user_index()] = iJet;
      }
      for (int i = 0; i < event.size(); ++i)
        if (std::find(hvConnectors.begin(), hvConnectors.end(), event[i].idAbs()) != hvConnectors.end()
            && event[i].status() == -22)
        {
          printEventBranch(i, event, diagramOut, jetIndexOfParticle, 0, "root");
        }
    }

    for (int iJet = 0; iJet < int(sortedJets.size()); ++iJet) {
      double fullJetPt = sortedJets[iJet].pt();
      pTj.fill(fullJetPt);

      double darkHadronEnergyFrac = darkHadronFractions(sortedJets[iJet], event)[1];
      double darkHadronPtFrac = darkHadronFractions(sortedJets[iJet], event)[2];

      fastjet::PseudoJet j = sortedJets[iJet];
      fastjet::PseudoJet j1, j2;

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
                << std::setw(realWidth-3) << darkHadronEnergyFrac
                << std::setw(realWidth-3) << clusterE
                << std::setw(realWidth-3) << darkHadronPtFrac
                << std::setw(realWidth-3) << clusterPt
                << std::setprecision(precision)
                << std::setw(realWidth) << ln_oneoverDelta
                << std::setw(realWidth) << ln_kt
                << std::setw(realWidth) << z
                << std::setw(realWidth) << psi
                << std::setw(realWidth) << m2
                << std::setw(5) << multiplicity << "\n";

        j = j1;
      }
    }

    // End of event loop. Print statistics and histograms.
    }
  pythia.stat();
  cout << nGluonv << nGammav << nHadronv << pTj;

  if (pythia.settings.flag("Main:writeLog"))
    cout.rdbuf(oldCout);

  // Done.
  return 0;
}
