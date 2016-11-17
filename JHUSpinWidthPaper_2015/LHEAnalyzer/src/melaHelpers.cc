#include "../interface/melaHelpers.h"


namespace melaHelpers{
  Mela* melaHandle=0;
  Float_t genPOLEWidth=4.07e-3;
  Float_t standardPOLEWidth=4.07e-3;
}

using namespace PDGHelpers;

void melaHelpers::setSamplePoleWidth(Float_t width_){ genPOLEWidth = width_; }
void melaHelpers::setStandardPoleWidth(Float_t width_){ standardPOLEWidth = width_; }


Float_t melaHelpers::melaBranchMEInterpreter(const ZZCandidate* cand, string& branchname){
  // Do absolutely nothing right now!
  return 0.;
}

