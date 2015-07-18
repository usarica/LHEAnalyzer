#ifndef HVVTREE_H
#define HVVTREE_H

#include "BaseTree.h"
#include "ZZCandidate.h"
#include <ZZMatrixElement/MELA/interface/Mela.h>
#include <ZZMatrixElement/MELA/src/computeAngles.h>

class HVVTree : public BaseTree{
public:
  HVVTree():BaseTree(){};
  HVVTree(string treename) : BaseTree(treename){};
  HVVTree(string treename, string treetitle) : BaseTree(treename, treetitle){};

  void bookAllBranches();

  void fillEventVariables(Float_t weight, Int_t passSelection);
  void fillMotherInfo(Particle* mother);

  void fillCandidate(ZZCandidate* pH, bool isGen=false);
  void fillCandidateDaughters(ZZCandidate* pH, bool isGen=false);
  void fillDaughterProducts(ZZCandidate* pH, bool isGen=false);
  void fillAssociatedInfo(ZZCandidate* pH, bool isGen=false);

//  void calculateDecayAngles(TLorentzVector thep4M11, TLorentzVector thep4M12, TLorentzVector thep4M21, TLorentzVector thep4M22, Float_t& costheta1, Float_t& costheta2, Float_t& phi, Float_t& costhetastar, Float_t& phistar1);
  void fillDecayAngles(ZZCandidate* pH, bool isGen=false);
//  void fillProductionAngles(ZZCandidate* pH, bool isGen=false);

};

#endif
