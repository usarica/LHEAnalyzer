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
  vector<string> getMELABranchList()const;

  bool reserveBranch(string branchname, const BaseTree::BranchTypes& branchtype, const bool& doSetAddress);
  void bookAllBranches(const bool& doSetAddress);

  void fillEventVariables(const Float_t& weight, const Int_t& passSelection);
  void fillMotherInfo(const MELAParticle* mother);

  void fillCandidate(MELACandidate* pH, bool isGen=false);
  void fillCandidateDaughters(MELACandidate* pH, bool isGen=false);
  void fillDaughterProducts(MELACandidate* pH, bool isGen=false);
  void fillAssociatedInfo(MELACandidate* pH, bool isGen=false);

  void fillDecayAngles(bool isGen=false);
  void fillVBFProductionAngles(bool isGen=false);
  void fillVHProductionAngles(bool isGen=false);
  void fillMELAProbabilities(bool isGen);

protected:
  void bookPtEtaPhiMassIdBranches(const string& owner, const BaseTree::BranchTypes& btype, const bool& doSetAddress, const bool& addId, const bool& usePz, bool isGen);
  void getPtEtaPhiMIdBranches(vector<string>& blist, const string& owner, const bool& addId, const bool& usePz, bool isGen);

  void bookAngularBranches(const bool& doSetAddress);
  void getAngularBranches(vector<string>& blist, const Int_t& prodFlag /* 0: Decay, 1: VBF, 2: VH */, bool isGen);

  void bookMELABranches(const bool& doSetAddress);

  vector<string> constructMELABranchList(const bool& doSetAddress);
  void setupMELASignalMECases(vector<string>& accumulatedlist, TVar::Production prod, TVar::MatrixElement me, bool isGen, bool isProdME, bool doSetAddress);
  vector<string> getMELASignalMEBranches(TVar::Production prod, TVar::MatrixElement me, vector<string> gList, vector<int> gCountRe, vector<int> gCountIm, bool isGen, bool isProdME, bool doSetAddress);

  OptionParser* options;
  vector<string> melaProbBranches;
};

#endif
