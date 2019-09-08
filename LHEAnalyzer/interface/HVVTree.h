#ifndef HVVTREE_H
#define HVVTREE_H

#include "BaseTree.h"
#include "OptionParser.h"
#include "GMECHelperFunctions.h"


using namespace BranchHelpers;


class HVVTree : public BaseTree{
protected:
  OptionParser* options;

  std::vector<MELAOptionParser*> recome_originalopts; // Be careful: Only for reading
  std::vector<MELAOptionParser*> recome_copyopts;
  std::vector<MELAHypothesis*> recome_units;
  std::vector<MELAHypothesis*> recome_aliased_units;
  std::vector<MELAComputation*> recome_computers;
  std::vector<MELACluster*> recome_clusters;
  std::vector<MELABranch*> recome_branches;

  std::vector<MELAOptionParser*> lheme_originalopts; // Be careful: Only for reading
  std::vector<MELAOptionParser*> lheme_copyopts;
  std::vector<MELAHypothesis*> lheme_units;
  std::vector<MELAHypothesis*> lheme_aliased_units;
  std::vector<MELAComputation*> lheme_computers;
  std::vector<MELACluster*> lheme_clusters;
  std::vector<MELABranch*> lheme_branches;

public:
  HVVTree();
  HVVTree(std::string treename);
  HVVTree(std::string treename, std::string treetitle);
  HVVTree(std::string treename, TFile* fin);
  virtual ~HVVTree();

  void setOptions(OptionParser* options_){ options=options_; }

  std::vector<MELABranch*>* getRecoMELABranches(){ return &recome_branches; }
  std::vector<MELABranch*>* getLHEMELABranches(){ return &lheme_branches; }
  std::vector<MELABranch*> const* getRecoMELABranches() const{ return &recome_branches; }
  std::vector<MELABranch*> const* getLHEMELABranches() const{ return &lheme_branches; }

  bool reserveBranch(std::string branchname, const BaseTree::BranchTypes& branchtype, const bool& doSetAddress);
  void bookAllBranches(const bool& doSetAddress);

  void fillXsec(const Float_t& val, const Float_t& err);
  void fillEventVariables(const Float_t& weight, const Int_t& passSelection);
  void fillMotherInfo(const MELAParticle* mother);

  void fillCandidate(MELACandidate* pH, bool isGen=false);
  void fillCandidateDaughters(MELACandidate* pH, bool isGen=false);
  void fillDaughterProducts(MELACandidate* pH, bool isGen=false);
  void fillAssociatedInfo(MELACandidate* pH, bool isGen=false);

  void fillDecayAngles(bool isGen=false);
  void fillVBFProductionAngles(bool isGen=false);
  void fillVHProductionAngles(bool isGen=false);
  void fillTTHProductionAngles(bool isGen=false);
  void fillMELAProbabilities(bool isGen);

protected:
  void bookPtEtaPhiMassIdBranches(const std::string& owner, const BaseTree::BranchTypes& btype, const bool& doSetAddress, const bool& addId, const bool& usePz, bool isGen);
  void bookMotherParticleBranches(const BaseTree::BranchTypes& btype, const bool& doSetAddress);
  void getPtEtaPhiMIdBranches(std::vector<std::string>& blist, const std::string& owner, const bool& addId, const bool& usePz, bool isGen);

  void bookDisplacementBranches(const std::string& owner, const BaseTree::BranchTypes& btype, const bool& doSetAddress, bool isGen);
  void getDisplacementBranches(std::vector<std::string>& blist, const std::string& owner, bool isGen);

  void bookAngularBranches(const bool& doSetAddress);
  void getAngularBranches(std::vector<std::string>& blist, const Int_t& prodFlag /* 0: Decay, 1: VBF, 2: VH */, bool isGen);

  void buildMELABranches(bool doSetAddress);
  void bookMELABranches(MELAOptionParser* me_opt, MELAComputation* computer, bool doCopy);
  void clearMELABranches();

  void computeMELABranches(bool isGen);
  void pushMELABranches(bool isGen);

  void updateMELAClusters_Common(const std::string clustertype, bool isGen);
  void updateMELAClusters_J1JEC(const std::string clustertype, bool isGen);
  void updateMELAClusters_J2JEC(const std::string clustertype, bool isGen);
  void updateMELAClusters_LepWH(const std::string clustertype, bool isGen);
  void updateMELAClusters_LepZH(const std::string clustertype, bool isGen);
  void updateMELAClusters_NoInitialQ(const std::string clustertype, bool isGen);
  void updateMELAClusters_NoInitialG(const std::string clustertype, bool isGen);
  void updateMELAClusters_NoAssociatedG(const std::string clustertype, bool isGen);
  void updateMELAClusters_NoInitialGNoAssociatedG(const std::string clustertype, bool isGen);
  void updateMELAClusters_BestLOAssociatedZ(const std::string clustertype, bool isGen);
  void updateMELAClusters_BestLOAssociatedW(const std::string clustertype, bool isGen);
  void updateMELAClusters_BestLOAssociatedVBF(const std::string clustertype, bool isGen);
  void updateMELAClusters_BestNLOVHApproximation(const std::string clustertype, bool isGen);
  void updateMELAClusters_BestNLOVBFApproximation(const std::string clustertype, bool isGen);

};

#endif
