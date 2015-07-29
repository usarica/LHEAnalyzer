#ifndef HVVTREE_H
#define HVVTREE_H

#include "BaseTree.h"
#include "OptionParser.h"

class HVVTree : public BaseTree{
public:
  HVVTree():BaseTree(), options(0){}
  HVVTree(string treename) : BaseTree(treename), options(0){}
  HVVTree(string treename, string treetitle) : BaseTree(treename, treetitle), options(0){}
  HVVTree(string treename, TFile* fin) : BaseTree(treename, fin), options(0){}

  void setOptions(OptionParser* options_){ options=options_; }

  bool reserveBranch(string branchname, BaseTree::BranchTypes branchtype, bool doSetAddress);
  void bookAllBranches(bool doSetAddress);

  void fillEventVariables(Float_t weight, Int_t passSelection);
  void fillMotherInfo(Particle* mother);

  void fillCandidate(ZZCandidate* pH, bool isGen=false);
  void fillCandidateDaughters(ZZCandidate* pH, bool isGen=false);
  void fillDaughterProducts(ZZCandidate* pH, bool isGen=false);
  void fillAssociatedInfo(ZZCandidate* pH, bool isGen=false);

//  void calculateDecayAngles(TLorentzVector thep4M11, TLorentzVector thep4M12, TLorentzVector thep4M21, TLorentzVector thep4M22, Float_t& costheta1, Float_t& costheta2, Float_t& phi, Float_t& costhetastar, Float_t& phistar1);
  void fillDecayAngles(ZZCandidate* pH, bool isGen=false);
//  void fillProductionAngles(ZZCandidate* pH, bool isGen=false);

protected:
  OptionParser* options;
};

#endif
