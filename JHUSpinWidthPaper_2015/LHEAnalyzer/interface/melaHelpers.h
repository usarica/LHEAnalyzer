#ifndef MELAHELPERS_H
#define MELAHELPERS_H

#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <utility>
#include <algorithm>
#include <ZZMatrixElement/MELA/interface/Mela.h>
#include <ZZMatrixElement/MELA/src/computeAngles.h>
#include "Event.h"

namespace melaHelpers{
  extern Mela* melaHandle;
  extern Float_t genPOLEWidth;
  extern Float_t standardPOLEWidth;

  void setSamplePoleWidth(Float_t width_);
  void setStandardPoleWidth(Float_t width_);

  Float_t melaBranchMEInterpreter(const ZZCandidate* cand, string& branchname);

  void computeVHAngles(
    TLorentzVector thep4Z1, TLorentzVector thep4Z2, TLorentzVector thep4H,
    double& costheta1, double& costheta2, double& Phi, double& costhetastar, double& Phi1
    );
}

#endif
