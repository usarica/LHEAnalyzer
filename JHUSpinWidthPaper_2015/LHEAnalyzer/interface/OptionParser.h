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
#include "PDGHelpers.h"
#include "ParticleComparators.h"
#include "HiggsComparators.h"

using namespace std;

class OptionParser{
public:
  OptionParser(int argc, char** argv);
  ~OptionParser(){};

  void analyze();
  void splitOption(string rawoption, string& wish, string& value);
  void interpretOption(string wish, string value);
  void printOptionsHelp();
  void printOptionSummary();

  Double_t mH(){ return mPOLE; }
  Double_t GammaH(){ return wPOLE; }
  Int_t sqrts(){ return erg_tev; }
  Bool_t processGenInfo(){ bool doProcess=true; if (includeGenInfo==0) doProcess=false; return doProcess; }
  Bool_t processRecoInfo(){ bool doProcess=true; if (includeRecoInfo==0) doProcess=false; return doProcess; }
  Bool_t doRemoveLepMasses(){ bool doProcess=true; if (removeDaughterMasses==0) doProcess=false; return doProcess; }
  Int_t analysisLevel(){ return fileLevel; }
  Bool_t doGenHZZdecay(){
    bool doHZZ=true;
    if (isGenHZZ==0) doHZZ=false;
    return doHZZ;
  }
  Bool_t doRecoHZZdecay(){
    bool doHZZ=true;
    if (isRecoHZZ==0) doHZZ=false;
    return doHZZ;
  }
  Int_t genDecayProducts(){ return genDecayMode; }
  Int_t recoDecayProducts(){ return recoDecayMode; }
  HiggsComparators::CandidateSelection getHiggsCandidateSelectionScheme(bool isGen=false){ if (isGen) return genHiggsCandidateSelectionScheme; else return recoHiggsCandidateSelectionScheme; }
  string inputDir(){ return indir; }
  string outputDir(){ return outdir; }
  string getTempDir(){ return tmpDir; }
  string outputFilename(){ return coutput; }
  vector<string> inputfiles(){ return filename; }

protected:
  vector<string> rawOptions;

  Double_t mPOLE;
  Double_t wPOLE;
  Int_t erg_tev;
  Int_t includeGenInfo;
  Int_t includeRecoInfo;
  Int_t removeDaughterMasses;
  Int_t fileLevel;
  Int_t isGenHZZ;
  Int_t isRecoHZZ;
  Int_t genDecayMode;
  Int_t recoDecayMode;
  HiggsComparators::CandidateSelection genHiggsCandidateSelectionScheme;
  HiggsComparators::CandidateSelection recoHiggsCandidateSelectionScheme;
  string indir;
  string outdir;
  string coutput;
  string tmpDir;

  vector<string> filename;
};

#endif
