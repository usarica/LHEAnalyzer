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
#include "HiggsComparators.h"

using namespace std;

class OptionParser{
protected:
  void extractSkippedEvents(const string& rawoption);
  void extractGlobalRecordSet(const string& rawoption);

  void configureMela();
  void deconfigureMela();
  void extractMelaGenProdId(string rawoption);

  Bool_t checkListVariable(const vector<string>& list, const string& var)const;

  vector<string> rawOptions;

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
  string jetAlgo;

  string indir;
  string outdir;
  string coutput;
  string tmpDir;
  vector<string> filename;
  vector<string> excludedBranch;
  Int_t maxEvents;

  // Mela probabilities to include, has to be in abbreviated form (eg. "All", "None", "p0plus", "g1", "g1_prime2" etc.)
  pair<TVar::Production, TVar::MatrixElement> sampleProductionId;
  vector<string> includeGenDecayProb;
  vector<string> includeRecoDecayProb;
  vector<string> includeGenProdProb;
  vector<string> includeRecoProdProb;

  vector < pair<Int_t, Int_t> > eventSkipRanges;
  vector < pair<string, string> > globalRecordSet;

public:
  OptionParser(int argc, char** argv);
  ~OptionParser(){ deconfigureMela(); };

  void analyze();
  void splitOption(const string& rawoption, string& wish, string& value, char delimiter='=');
  void splitOptionRecursive(const string& rawoption, vector<string>& splitoptions, char delimiter=',');
  void interpretOption(const string& wish, string value);
  void printOptionsHelp();
  void printOptionSummary();

  Bool_t processGenInfo(){ bool doProcess=true; if (includeGenInfo==0) doProcess=false; return doProcess; }
  Bool_t processRecoInfo(){ bool doProcess=true; if (includeRecoInfo==0) doProcess=false; return doProcess; }
  Bool_t isAnExcludedBranch(const string& branchname);
  Int_t analysisLevel(){ return fileLevel; }
  Int_t pythiaType(){ return pythiaStep; }
  string jetAlgorithm(){ return jetAlgo; }
  Int_t doGenHZZdecay(){ return isGenHZZ; }
  Int_t doRecoHZZdecay(){ return isRecoHZZ; }
  Int_t genDecayProducts(){ return genDecayMode; }
  Int_t recoDecayProducts(){ return recoDecayMode; }
  Int_t recoSelectionMode(){ return recoSelBehaviour; }
  Int_t recoSmearingMode(){ return recoSmearBehaviour; }
  HiggsComparators::CandidateSelection getHiggsCandidateSelectionScheme(bool isGen=false){ if (isGen) return genHiggsCandidateSelectionScheme; else return recoHiggsCandidateSelectionScheme; }
  string inputDir(){ return indir; }
  string outputDir(){ return outdir; }
  string getTempDir(){ return tmpDir; }
  string outputFilename(){ return coutput; }
  vector<string> inputfiles(){ return filename; }
  Int_t maxEventsToProcess(){ return maxEvents; };
  vector < pair<Int_t, Int_t> > getSkippedEvents(){ return eventSkipRanges; };

  Bool_t hasGlobalRecord(){ return !globalRecordSet.empty(); }
  vector < pair<string, string> > getGlobalRecordSet(){ return globalRecordSet; }

  // MELA-related options
  Double_t mH(){ return mPOLE; }
  Double_t GammaH(){ return wPOLE; }
  Int_t sqrts(){ return erg_tev; }
  Bool_t initializeMELABranches(){ return (!includeGenDecayProb.empty() && !includeRecoDecayProb.empty() && !includeGenProdProb.empty() && !includeRecoProdProb.empty()); }
  Bool_t doRemoveLepMasses(){ bool doProcess=true; if (removeDaughterMasses==0) doProcess=false; return doProcess; }
  Bool_t doComputeDecayAngles(){ bool doProcess=true; if (computeDecayAngles==0) doProcess=false; return doProcess; }
  Bool_t doComputeVBFAngles(){ bool doProcess=true; if (computeVBFAngles==0) doProcess=false; return doProcess; }
  Bool_t doComputeVHAngles(){ bool doProcess=true; if (computeVHAngles==0) doProcess=false; return doProcess; }
  Int_t computeVHAnglesOrder(){ return computeVHAngles; }
  Bool_t hasGenDecayME(const string& str);
  Bool_t hasRecoDecayME(const string& str);
  Bool_t hasGenProdME(const string& str);
  Bool_t hasRecoProdME(const string& str);

  const Bool_t& doRecastGenTopologyToLOQCDVH() const{ return recastGenTopologyToLOQCDVH; }
  const Bool_t& doRecastGenTopologyToLOQCDVBF() const{ return recastGenTopologyToLOQCDVBF; }

  pair<TVar::Production, TVar::MatrixElement> getSampleProductionId(){ return sampleProductionId; }

};

#endif
