#ifndef CONVERT_LHE_H
#define CONVERT_LHE_H

#include "converter.h"

using namespace std;

class convertLHE : public converter{
public:
  convertLHE(OptionParser* options_);
  ~convertLHE(){};
  void run();

protected:
  void configure(); // Set output file, tree
  void finalizeRun();
  vector<Particle*> readEvent(ifstream& input_lhe, double& weight);
};
#endif
