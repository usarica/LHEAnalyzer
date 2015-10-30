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
#include "TLorentzRotation.h"
#include "Event.h"

namespace melaHelpers{
  extern Mela* melaHandle;
  extern Float_t genPOLEWidth;
  extern Float_t standardPOLEWidth;

  void setSamplePoleWidth(Float_t width_);
  void setStandardPoleWidth(Float_t width_);

  Float_t melaBranchMEInterpreter(const ZZCandidate* cand, string& branchname);


  void constrainedRemoveLeptonMass(TLorentzVector& p1, TLorentzVector& p2);

  void computeAngles(
    TLorentzVector p4M11, int Z1_lept1Id,
    TLorentzVector p4M12, int Z1_lept2Id,
    TLorentzVector p4M21, int Z2_lept1Id,
    TLorentzVector p4M22, int Z2_lept2Id,
    float& costhetastar,
    float& costheta1,
    float& costheta2,
    float& Phi,
    float& Phi1);

  void computeVBFangles(
    float& costhetastar,
    float& costheta1,
    float& costheta2,
    float& Phi,
    float& Phi1,
    float& Q2V1,
    float& Q2V2,
    TLorentzVector p4M11, int Z1_lept1Id,
    TLorentzVector p4M12, int Z1_lept2Id,
    TLorentzVector p4M21, int Z2_lept1Id,
    TLorentzVector p4M22, int Z2_lept2Id,
    TLorentzVector jet1, int jet1Id,
    TLorentzVector jet2, int jet2Id,
    TLorentzVector* injet1=0, int injet1Id=0,
    TLorentzVector* injet2=0, int injet2Id=0);
  void computeVHangles(
    float& costhetastar,
    float& costheta1,
    float& costheta2,
    float& Phi,
    float& Phi1,
    TLorentzVector p4M11, int Z1_lept1Id,
    TLorentzVector p4M12, int Z1_lept2Id,
    TLorentzVector p4M21, int Z2_lept1Id,
    TLorentzVector p4M22, int Z2_lept2Id,
    TLorentzVector jet1, int jet1Id,
    TLorentzVector jet2, int jet2Id,
    TLorentzVector* injet1=0, int injet1Id=0,
    TLorentzVector* injet2=0, int injet2Id=0);
}

#endif
