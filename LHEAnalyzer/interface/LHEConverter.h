#ifndef LHECONVERTER_H
#define LHECONVERTER_H

#include "converter.h"


class LHEConverter : public converter{
public:
  LHEConverter(OptionParser* options_);
  ~LHEConverter(){};
  void run();

protected:
  void configure(); // Set output file, tree
  void finalizeRun();
  std::vector<MELAParticle*> readEvent(std::ifstream& input_lhe, int& fline, double& weight);

};
#endif
