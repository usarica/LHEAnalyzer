#ifndef OPTIONPARSER_H
#define OPTIONPARSER_H

#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>
#include <utility>
#include <cassert>
#include "TString.h"
#include "melaHelpers.h"
#include "ParticleComparators.h"
#include "TopComparators.h"
#include "HiggsComparators.h"


class OptionParser{
protected:
  void extractSkippedEvents(std::string const& rawoption);
  void extractGlobalRecordSet(std::string const& rawoption);

  void configureMela();
  void deconfigureMela();
  void extractMelaGenProdId(std::string const& rawoption);

  Bool_t checkListVariable(std::vector<std::string> const& list, std::string const& var)const;

  std::vector<std::string> rawOptions;

  Double_t mPOLE;
  Double_t wPOLE;
  Double_t wPOLEStandard;
  Int_t erg_tev;
  Int_t includeGenInfo;
  Int_t includeRecoInfo;
  Int_t removeDaughterMasses;
  Int_t computeDecayAngles;
  Int_t computeVBFAngles;
  Int_t computeVHAngles;
  Int_t computeTTHAngles;
  Int_t fileLevel;
  Int_t pythiaStep;
  Int_t isGenHZZ;
  Int_t isRecoHZZ;
  Int_t genDecayMode;
  Int_t recoDecayMode;
  Int_t recoSelBehaviour;
  Int_t recoSmearBehaviour;
  HiggsComparators::CandidateSelection genHiggsCandidateSelectionScheme;
  HiggsComparators::CandidateSelection recoHiggsCandidateSelectionScheme;

  Bool_t recastGenTopologyToLOQCDVH;
  Bool_t recastGenTopologyToLOQCDVBF;

  Double_t jetDeltaRIso;
  std::string jetAlgo;

  std::string indir;
  std::string outdir;
  std::string coutput;
  std::string tmpDir;
  std::vector<std::string> filename;
  std::vector<std::string> excludedBranch;
  Int_t maxEvents;

  // Mela probabilities to include, has to be in abbreviated form (eg. "All", "None", "p0plus", "g1", "g1_prime2" etc.)
  std::pair<TVar::Production, TVar::MatrixElement> sampleProductionId;
  std::vector<std::string> includeGenDecayProb;
  std::vector<std::string> includeRecoDecayProb;
  std::vector<std::string> includeGenProdProb;
  std::vector<std::string> includeRecoProdProb;

  std::vector<std::pair<Int_t, Int_t>> eventSkipRanges;
  std::vector<std::pair<std::string, std::string>> globalRecordSet;

public:
  OptionParser(int argc, char** argv);
  ~OptionParser(){ deconfigureMela(); };

  void analyze();
  void splitOption(std::string const& rawoption, std::string& wish, std::string& value, char delimiter='=');
  void splitOptionRecursive(std::string const& rawoption, std::vector<std::string>& splitoptions, char delimiter=',');
  void interpretOption(std::string const& wish, std::string const& value);
  void printOptionsHelp();
  void printOptionSummary();

  Bool_t processGenInfo(){ bool doProcess=true; if (includeGenInfo==0) doProcess=false; return doProcess; }
  Bool_t processRecoInfo(){ bool doProcess=true; if (includeRecoInfo==0) doProcess=false; return doProcess; }
  Bool_t isAnExcludedBranch(std::string const& branchname);
  Int_t analysisLevel(){ return fileLevel; }
  Int_t pythiaType(){ return pythiaStep; }
  std::string jetAlgorithm(){ return jetAlgo; }
  Int_t doGenHZZdecay(){ return isGenHZZ; }
  Int_t doRecoHZZdecay(){ return isRecoHZZ; }
  Int_t genDecayProducts(){ return genDecayMode; }
  Int_t recoDecayProducts(){ return recoDecayMode; }
  Int_t recoSelectionMode(){ return recoSelBehaviour; }
  Int_t recoSmearingMode(){ return recoSmearBehaviour; }
  HiggsComparators::CandidateSelection getHiggsCandidateSelectionScheme(bool isGen=false){ if (isGen) return genHiggsCandidateSelectionScheme; else return recoHiggsCandidateSelectionScheme; }
  std::string inputDir(){ return indir; }
  std::string outputDir(){ return outdir; }
  std::string getTempDir(){ return tmpDir; }
  std::string outputFilename(){ return coutput; }
  std::vector<std::string> inputfiles(){ return filename; }
  Int_t maxEventsToProcess(){ return maxEvents; };
  std::vector<std::pair<Int_t, Int_t>> getSkippedEvents(){ return eventSkipRanges; };

  Bool_t hasGlobalRecord(){ return !globalRecordSet.empty(); }
  std::vector<std::pair<std::string, std::string>> getGlobalRecordSet(){ return globalRecordSet; }

  // MELA-related options
  Double_t mH(){ return mPOLE; }
  Double_t GammaH(){ return wPOLE; }
  Int_t sqrts(){ return erg_tev; }
  Bool_t initializeMELABranches(){ return (!includeGenDecayProb.empty() && !includeRecoDecayProb.empty() && !includeGenProdProb.empty() && !includeRecoProdProb.empty()); }
  Bool_t doRemoveLepMasses(){ bool doProcess=true; if (removeDaughterMasses==0) doProcess=false; return doProcess; }
  Bool_t doComputeDecayAngles(){ bool doProcess=true; if (computeDecayAngles==0) doProcess=false; return doProcess; }
  Bool_t doComputeVBFAngles(){ bool doProcess=true; if (computeVBFAngles==0) doProcess=false; return doProcess; }
  Bool_t doComputeVHAngles(){ bool doProcess=true; if (computeVHAngles==0) doProcess=false; return doProcess; }
  Bool_t doComputeTTHAngles(){ bool doProcess=true; if (computeTTHAngles==0) doProcess=false; return doProcess; }
  Int_t computeVHAnglesOrder(){ return computeVHAngles; }
  Bool_t hasGenDecayME(std::string const& str);
  Bool_t hasRecoDecayME(std::string const& str);
  Bool_t hasGenProdME(std::string const& str);
  Bool_t hasRecoProdME(std::string const& str);

  const Bool_t& doRecastGenTopologyToLOQCDVH() const{ return recastGenTopologyToLOQCDVH; }
  const Bool_t& doRecastGenTopologyToLOQCDVBF() const{ return recastGenTopologyToLOQCDVBF; }

  std::pair<TVar::Production, TVar::MatrixElement> getSampleProductionId(){ return sampleProductionId; }

};

#endif
