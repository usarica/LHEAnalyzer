#ifndef CONVERTER_H
#define CONVERTER_H

#include <iostream>
#include <fstream>
#include <cstdio>
#include <string>
#include <iomanip>
#include <vector>
#include "LHEParticleSmear.h"
#include "HVVTree.h"


class converter{
public:
  converter(OptionParser* options_);
  virtual ~converter(){ finalizeRun(); }

  virtual void run()=0;

protected:
  virtual void configure(); // Set output file, tree
  virtual void finalizeRun();

  OptionParser* options;
  std::vector<std::string> filename; // Input files
  TFile* foutput;
  HVVTree* tree;

};
#endif
