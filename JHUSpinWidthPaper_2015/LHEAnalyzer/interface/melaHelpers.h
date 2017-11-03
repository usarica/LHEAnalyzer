#ifndef MELAHELPERS_H
#define MELAHELPERS_H

#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <utility>
#include <algorithm>
#include <ZZMatrixElement/MELA/interface/Mela.h>
#include <MelaAnalytics/CandidateLOCaster/interface/MELACandidateRecaster.h>
#include "Event.h"
#include "TLorentzRotation.h"

namespace melaHelpers{
  extern Mela* melaHandle;
  extern Float_t genPOLEWidth;
  extern Float_t standardPOLEWidth;

  void setSamplePoleWidth(Float_t width_);
  void setStandardPoleWidth(Float_t width_);

  Float_t melaBranchMEInterpreter(const ZZCandidate* cand, string& branchname);
}

#endif
