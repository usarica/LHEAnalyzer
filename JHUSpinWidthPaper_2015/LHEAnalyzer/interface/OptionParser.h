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
  Bool_t doHZZdecay(){ bool doHZZ=true; if (isHZZ==0) doHZZ=false; return doHZZ; }
  Int_t decayProducts(){ return decayMode; }
  string inputDir(){ return indir; }
  string outputDir(){ return outdir; }
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
  Int_t isHZZ;
  Int_t decayMode;
  string indir;
  string outdir;
  string coutput;

  vector<string> filename;
};

#endif
