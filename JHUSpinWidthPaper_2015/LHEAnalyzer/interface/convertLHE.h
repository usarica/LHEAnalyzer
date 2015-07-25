#ifndef CONVERT_LHE_H
#define CONVERT_LHE_H

#include <iostream>
#include <fstream>
#include <cstdio>
#include <string>
#include <iomanip>
#include <vector>
#include "LHEParticleSmear.h"
#include "HVVTree.h"
#include "HiggsComparators.h"


using namespace std;

class convertLHE{
public:
  convertLHE(OptionParser* options_);
  ~convertLHE(){};
  void run();

protected:
  void configure(OptionParser* options_); // Set output file, tree
  void finalizeRun();
  vector<Particle*> readEvent(ifstream& input_lhe, double& weight);

  OptionParser* options;
  vector<string> filename;
  TFile* foutput;
  HVVTree* tree;
};
#endif
