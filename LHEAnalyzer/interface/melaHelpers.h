#ifndef MELAHELPERS_H
#define MELAHELPERS_H

#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <utility>
#include <algorithm>
#include "Mela.h"
#include "MELACandidateRecaster.h"
#include "MELAEvent.h"
#include "TLorentzRotation.h"

namespace melaHelpers{
  extern Mela* melaHandle;
  extern Float_t genPOLEWidth;
  extern Float_t standardPOLEWidth;

  void setSamplePoleWidth(Float_t width_);
  void setStandardPoleWidth(Float_t width_);

  TVar::Production getFirstAssociatedLeptonicVProduction(MELACandidate const* cand);
  TVar::Production getFirstAssociatedHadronicVProduction(MELACandidate const* cand);

  Float_t melaBranchMEInterpreter(std::string& branchname);

}

#endif
