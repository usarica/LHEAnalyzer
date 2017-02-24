#include "TString.h"

TString user_dir = "/afs/cern.ch/work/u/usarica/scratch-ZZdevelopment/CJLST_4l/CMSSW_8_0_26_patch1/src/ZZAnalysis/AnalysisStep/test/prod/PROD/";

enum{
  nBkgSamples=0,
  //nBkgSamples=1,
  nSigSamples=4,
  //nSamples=5
  nSamples=4
};
pair<TString, TString> strSamples[nSamples]={
  //pair<TString, TString>("DYJetsToLL_M50_Chunk0/ZZ4lAnalysis.root", "DY"),
  pair<TString, TString>("ggH115_Chunk0/ZZ4lAnalysis.root", "ggH115"),
  pair<TString, TString>("ggH125_Chunk0/ZZ4lAnalysis.root", "ggH125"),
  pair<TString, TString>("ggH150_Chunk0/ZZ4lAnalysis.root", "ggH150"),
  pair<TString, TString>("ggH750_Chunk0/ZZ4lAnalysis.root", "ggH750")

};

enum{
  nChannels=3
};

TString channame[nChannels]={ "4mu", "4e", "2e2mu" };
TString chanlabel[nChannels]={ "4#mu", "4e", "2e2#mu" };

