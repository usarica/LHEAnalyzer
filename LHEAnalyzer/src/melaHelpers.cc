#include "melaHelpers.h"


namespace melaHelpers{
  Mela* melaHandle=nullptr;
  Float_t genPOLEWidth=4.07e-3;
  Float_t standardPOLEWidth=4.07e-3;
}

using namespace PDGHelpers;

void melaHelpers::setSamplePoleWidth(Float_t width_){ genPOLEWidth = width_; }
void melaHelpers::setStandardPoleWidth(Float_t width_){ standardPOLEWidth = width_; }

TVar::Production melaHelpers::getFirstAssociatedLeptonicVProduction(MELACandidate const* cand){
  TVar::Production res=TVar::nProductions;
  if (!cand) return res;
  for (const MELAParticle* tmpV:cand->getAssociatedSortedVs()){
    if (tmpV->getDaughter(0) && PDGHelpers::isALepton(tmpV->getDaughter(0)->id)){
      if (PDGHelpers::isAWBoson(tmpV->id)) res=TVar::Lep_WH;
      else if (PDGHelpers::isAZBoson(tmpV->id)) res=TVar::Lep_ZH;
    }
    else if (PDGHelpers::isAPhoton(tmpV->id)) res=TVar::GammaH;
    if (res!=TVar::nProductions) break;
  }
  return res;
}
TVar::Production melaHelpers::getFirstAssociatedHadronicVProduction(MELACandidate const* cand){
  TVar::Production res=TVar::nProductions;
  if (!cand) return res;
  for (const MELAParticle* tmpV:cand->getAssociatedSortedVs()){
    if (tmpV->getDaughter(0) && (PDGHelpers::isAJet(tmpV->getDaughter(0)->id) && !PDGHelpers::isAGluon(tmpV->getDaughter(0)->id))){
      if (PDGHelpers::isAWBoson(tmpV->id)) res=TVar::Had_WH;
      else if (PDGHelpers::isAZBoson(tmpV->id)) res=TVar::Had_ZH;
    }
    if (res!=TVar::nProductions) break;
  }
  return res;
}

Float_t melaHelpers::melaBranchMEInterpreter(std::string& branchname){
  // Do absolutely nothing right now!
  //MELACandidate* cand=melaHandle->getCurrentCandidate(),
  return 0.;
}

