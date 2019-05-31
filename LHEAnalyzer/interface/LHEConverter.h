#ifndef LHECONVERTER_H
#define LHECONVERTER_H

#include "converter.h"

using namespace std;

class LHEConverter : public converter{
public:
  LHEConverter(OptionParser* options_);
  ~LHEConverter(){};
  void run();

protected:
  void configure(); // Set output file, tree
  void finalizeRun();
  vector<MELAParticle*> readEvent(ifstream& input_lhe, int& fline, double& weight);
};
#endif
