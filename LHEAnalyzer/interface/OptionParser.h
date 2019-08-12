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
#include "MELAEvent.h"
#include "melaHelpers.h"
#include "ParticleComparators.h"
#include "TopComparators.h"
#include "HiggsComparators.h"


class OptionParser{
protected:
  void extractSkippedEvents(std::string const& rawoption);
  void extractGlobalRecordSet(std::string const& rawoption);

  void configureMela()const;
  void deconfigureMela()const;
  void extractMElines();
  void extractXsec();
  static void extractMElines(std::string const& sfile, std::vector<std::string>& llist);

  Bool_t checkListVariable(std::vector<std::string> const& list, std::string const& var)const;

  std::vector<std::string> rawOptions;

  Double_t mPOLE;
  Double_t wPOLE;
  Double_t wPOLEStandard;

  Float_t sampleXsec;
  Float_t sampleXsecErr;

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
  MELAEvent::CandidateVVMode isGenHZZ;
  MELAEvent::CandidateVVMode isRecoHZZ;
  Int_t genDecayMode;
  Int_t recoDecayMode;
  Int_t recoSelBehaviour;
  Int_t recoSmearBehaviour;

  HiggsComparators::CandidateSelection genHiggsCandidateSelectionScheme;
  HiggsComparators::CandidateSelection recoHiggsCandidateSelectionScheme;

  Double_t jetDeltaRIso;
  std::string jetAlgo;

  std::string indir;
  std::string outdir;
  std::string coutput;
  std::string tmpDir;
  std::vector<std::string> filename;
  std::vector<std::string> excludedBranch;
  Int_t maxEvents;

  // MELA probabilities to compute
  std::string lheMEfile;
  std::string recoMEfile;
  std::vector<std::string> lheMElist;
  std::vector<std::string> recoMElist;

  std::vector<std::pair<Int_t, Int_t>> eventSkipRanges;
  std::vector<std::pair<std::string, std::string>> globalRecordSet;

public:
  OptionParser(int argc, char** argv);
  ~OptionParser(){ deconfigureMela(); };

  void analyze();
  void interpretOption(std::string const& wish, std::string const& value);
  void printOptionsHelp(bool command_fail)const;

  Bool_t processGenInfo()const{ return (includeGenInfo!=0); }
  Bool_t processRecoInfo()const{ return (includeRecoInfo!=0); }
  Bool_t isAnExcludedBranch(std::string const& branchname) const;
  Int_t analysisLevel()const{ return fileLevel; }
  Int_t pythiaType()const{ return pythiaStep; }
  std::string jetAlgorithm()const{ return jetAlgo; }
  MELAEvent::CandidateVVMode doGenHZZdecay()const{ return isGenHZZ; }
  MELAEvent::CandidateVVMode doRecoHZZdecay()const{ return isRecoHZZ; }
  Int_t genDecayProducts()const{ return genDecayMode; }
  Int_t recoDecayProducts()const{ return recoDecayMode; }
  Int_t recoSelectionMode()const{ return recoSelBehaviour; }
  Int_t recoSmearingMode()const{ return recoSmearBehaviour; }
  HiggsComparators::CandidateSelection getHiggsCandidateSelectionScheme(bool isGen)const{ return (isGen ? genHiggsCandidateSelectionScheme : recoHiggsCandidateSelectionScheme); }
  std::string inputDir()const{ return indir; }
  std::string outputDir()const{ return outdir; }
  std::string getTempDir()const{ return tmpDir; }
  std::string outputFilename()const{ return coutput; }
  std::vector<std::string> const& inputfiles()const{ return filename; }
  std::vector<std::string> const& getLHEMEList()const{ return lheMElist; }
  std::vector<std::string> const& getRecoMEList()const{ return recoMElist; }

  Int_t maxEventsToProcess()const{ return maxEvents; };
  std::vector<std::pair<Int_t, Int_t>> getSkippedEvents()const{ return eventSkipRanges; };

  Bool_t hasGlobalRecord()const{ return !globalRecordSet.empty(); }
  std::vector<std::pair<std::string, std::string>> getGlobalRecordSet()const{ return globalRecordSet; }

  // MELA-related options
  Double_t mH()const{ return mPOLE; }
  Double_t GammaH()const{ return wPOLE; }
  Float_t get_xsec()const{ return sampleXsec; }
  Float_t get_xsecerr()const{ return sampleXsecErr; }

  Int_t sqrts()const{ return erg_tev; }
  Bool_t doMEComputations()const{ return !(lheMElist.empty() && recoMElist.empty()); }
  Bool_t doRemoveLepMasses()const{ return (removeDaughterMasses!=0); }
  Bool_t doComputeDecayAngles()const{ return (computeDecayAngles!=0); }
  Bool_t doComputeVBFAngles()const{ return (computeVBFAngles!=0); }
  Bool_t doComputeVHAngles()const{ return (computeVHAngles!=0); }
  Bool_t doComputeTTHAngles()const{ return (computeTTHAngles!=0); }
  Int_t computeVHAnglesOrder()const{ return computeVHAngles; }

};

#endif
