#include "HVVTree.h"


using namespace std;


HVVTree::HVVTree() :BaseTree(), options(0){}
HVVTree::HVVTree(std::string treename) : BaseTree(treename), options(0){}
HVVTree::HVVTree(std::string treename, std::string treetitle) : BaseTree(treename, treetitle), options(0){}
HVVTree::HVVTree(std::string treename, TFile* fin) : BaseTree(treename, fin), options(0){}
HVVTree::~HVVTree(){
  // Delete MELA branches
  clearMELABranches();
}

bool HVVTree::reserveBranch(string branchname, const BaseTree::BranchTypes& branchtype, const bool& doSetAddress){
  bool isAvailable = true;
  if (options->isAnExcludedBranch(branchname)){
    isAvailable=false;
    if (doSetAddress && hvvtree->GetBranchStatus(branchname.c_str())) hvvtree->SetBranchStatus(branchname.c_str(), 0);
  }
  else if (
    (doSetAddress && !hvvtree->GetBranchStatus(branchname.c_str()))
    ||
    (!doSetAddress && hvvtree->GetBranchStatus(branchname.c_str()))
    ) isAvailable=false; // If setAddress to a non-existing branch or branch to an existing address
  if (isAvailable) bookBranch(branchname, branchtype, doSetAddress);
  return isAvailable;
}
void HVVTree::bookAllBranches(const bool& doSetAddress){
  if (!options){
    cerr << "HVVTree::bookAllBranches -> No options are set for the HVVTree!" << endl;
    return;
  }

  reserveBranch("MC_weight", BaseTree::bFloat, doSetAddress);
  reserveBranch("xsec", BaseTree::bFloat, doSetAddress);
  reserveBranch("xsecerr", BaseTree::bFloat, doSetAddress);
  if (options->processGenInfo()){
    bookMotherParticleBranches(BaseTree::bVectorFloat, doSetAddress);
    bookPtEtaPhiMassIdBranches("H", BaseTree::bFloat, doSetAddress, false, true, true);
    bookPtEtaPhiMassIdBranches("Z1", BaseTree::bFloat, doSetAddress, false, false, true);
    bookPtEtaPhiMassIdBranches("Z2", BaseTree::bFloat, doSetAddress, false, false, true);
    bookPtEtaPhiMassIdBranches("Za", BaseTree::bFloat, doSetAddress, false, false, true);
    bookPtEtaPhiMassIdBranches("Zb", BaseTree::bFloat, doSetAddress, false, false, true);

    bookPtEtaPhiMassIdBranches("AssociatedParticle", BaseTree::bVectorFloat, doSetAddress, true, false, true);
    reserveBranch("GenDijetMass", BaseTree::bFloat, doSetAddress);
    reserveBranch("GenDileptonMass", BaseTree::bFloat, doSetAddress);
    reserveBranch("GenDijetVVMass", BaseTree::bFloat, doSetAddress);
    reserveBranch("GenDileptonVVMass", BaseTree::bFloat, doSetAddress);
    reserveBranch("GenDRjet", BaseTree::bFloat, doSetAddress);
    reserveBranch("GenDRlepton", BaseTree::bFloat, doSetAddress);

    reserveBranch("NGenAssociatedVs", BaseTree::bInt, doSetAddress);
    bookPtEtaPhiMassIdBranches("AssociatedV", BaseTree::bVectorFloat, doSetAddress, true, false, true);
    reserveBranch("GenAssociatedV_Particle1Index", BaseTree::bVectorInt, doSetAddress);
    reserveBranch("GenAssociatedV_Particle2Index", BaseTree::bVectorInt, doSetAddress);

    reserveBranch("NGenAssociatedTops", BaseTree::bInt, doSetAddress);
    bookPtEtaPhiMassIdBranches("AssociatedTop", BaseTree::bVectorFloat, doSetAddress, true, false, true);
    reserveBranch("GenAssociatedTop_PartnerParticleIndex", BaseTree::bVectorInt, doSetAddress);
    reserveBranch("GenAssociatedTop_WFermionIndex", BaseTree::bVectorInt, doSetAddress);
    reserveBranch("GenAssociatedTop_WAntifermionIndex", BaseTree::bVectorInt, doSetAddress);

    bookPtEtaPhiMassIdBranches("CandDau", BaseTree::bVectorFloat, doSetAddress, true, false, true);
  }
  if (options->processRecoInfo()){
    reserveBranch("isSelected", BaseTree::bInt, doSetAddress);

    bookPtEtaPhiMassIdBranches("ZZ", BaseTree::bFloat, doSetAddress, false, true, false);
    bookPtEtaPhiMassIdBranches("Z1", BaseTree::bFloat, doSetAddress, false, false, false);
    bookPtEtaPhiMassIdBranches("Z2", BaseTree::bFloat, doSetAddress, false, false, false);
    bookPtEtaPhiMassIdBranches("Za", BaseTree::bFloat, doSetAddress, false, false, false);
    bookPtEtaPhiMassIdBranches("Zb", BaseTree::bFloat, doSetAddress, false, false, false);

    bookPtEtaPhiMassIdBranches("AssociatedParticle", BaseTree::bVectorFloat, doSetAddress, true, false, false);
    reserveBranch("DijetMass", BaseTree::bFloat, doSetAddress);
    reserveBranch("DileptonMass", BaseTree::bFloat, doSetAddress);
    reserveBranch("DijetVVMass", BaseTree::bFloat, doSetAddress);
    reserveBranch("DileptonVVMass", BaseTree::bFloat, doSetAddress);
    reserveBranch("DRjet", BaseTree::bFloat, doSetAddress);
    reserveBranch("DRlepton", BaseTree::bFloat, doSetAddress);

    reserveBranch("NAssociatedVs", BaseTree::bInt, doSetAddress);
    bookPtEtaPhiMassIdBranches("AssociatedV", BaseTree::bVectorFloat, doSetAddress, true, false, false);
    reserveBranch("AssociatedV_Particle1Index", BaseTree::bVectorInt, doSetAddress);
    reserveBranch("AssociatedV_Particle2Index", BaseTree::bVectorInt, doSetAddress);

    reserveBranch("NAssociatedTops", BaseTree::bInt, doSetAddress);
    bookPtEtaPhiMassIdBranches("AssociatedTop", BaseTree::bVectorFloat, doSetAddress, true, false, false);
    reserveBranch("AssociatedTop_PartnerParticleIndex", BaseTree::bVectorInt, doSetAddress);
    reserveBranch("AssociatedTop_WFermionIndex", BaseTree::bVectorInt, doSetAddress);
    reserveBranch("AssociatedTop_WAntifermionIndex", BaseTree::bVectorInt, doSetAddress);

    bookPtEtaPhiMassIdBranches("CandDau", BaseTree::bVectorFloat, doSetAddress, true, false, false);
  }

  bookAngularBranches(doSetAddress);
  buildMELABranches(doSetAddress);
  actuateBranches(doSetAddress);
}


void HVVTree::bookPtEtaPhiMassIdBranches(const string& owner, const BaseTree::BranchTypes& btype, const bool& doSetAddress, const bool& addId, const bool& usePz, bool isGen){
  bool const btype_is_vector = (btype==BaseTree::bVectorFloat || btype==BaseTree::bVectorDouble);
  vector<string> tmpBranchList;
  getPtEtaPhiMIdBranches(tmpBranchList, owner, addId, usePz, isGen);
  for (unsigned int b=0; b<tmpBranchList.size(); b++){
    string branchname = tmpBranchList.at(b);
    if(!addId || branchname.find("Id")==string::npos) reserveBranch(tmpBranchList.at(b), btype, doSetAddress);
    else{
      BaseTree::BranchTypes bInttype = BaseTree::bInt;
      if (btype_is_vector) bInttype = BaseTree::bVectorInt;
      reserveBranch(tmpBranchList.at(b), bInttype, doSetAddress);
    }
  }
}
void HVVTree::bookMotherParticleBranches(const BaseTree::BranchTypes& btype, const bool& doSetAddress){
  const string strGen = "Gen";
  const string strOwner="Mother";
  bool const btype_is_vector = (btype==BaseTree::bVectorFloat || btype==BaseTree::bVectorDouble);
  reserveBranch((strGen+strOwner+"Pz"), btype, doSetAddress);
  reserveBranch((strGen+strOwner+"E"), btype, doSetAddress);
  reserveBranch((strGen+strOwner+"Id"), (btype_is_vector ? BaseTree::bVectorInt : BaseTree::bInt), doSetAddress);
}
void HVVTree::getPtEtaPhiMIdBranches(vector<string>& blist, const string& owner, const bool& addId, const bool& usePz, bool isGen){
  string strGen = "Gen";
  vector<string> strtmp;

  strtmp.push_back("Pt");
  if (usePz) strtmp.push_back("Pz");
  else strtmp.push_back("Eta");
  strtmp.push_back("Phi");
  strtmp.push_back("Mass");
  if (addId) strtmp.push_back("Id");

  for (string& varname:strtmp){
    varname.insert(0, owner);
    if (isGen) varname.insert(0, strGen);
    blist.push_back(varname);
  }
}
void HVVTree::bookAngularBranches(const bool& doSetAddress){
  vector<string> tmpBranchList;
  if (options->processGenInfo()){
    if (options->doComputeDecayAngles() || doSetAddress) getAngularBranches(tmpBranchList, 0, true);
    if (options->doComputeVBFAngles() || doSetAddress) getAngularBranches(tmpBranchList, 1, true);
    if (options->doComputeVHAngles() || doSetAddress){
      if (options->computeVHAnglesOrder()!=3) getAngularBranches(tmpBranchList, 2, true);
      else{
        getAngularBranches(tmpBranchList, 3, true);
        getAngularBranches(tmpBranchList, 4, true);
      }
    }
    if (options->doComputeTTHAngles() || doSetAddress) getAngularBranches(tmpBranchList, 5, true);
  }
  if (options->processRecoInfo()){
    if (options->doComputeDecayAngles() || doSetAddress) getAngularBranches(tmpBranchList, 0, false);
    if (options->doComputeVBFAngles() || doSetAddress) getAngularBranches(tmpBranchList, 1, false);
    if (options->doComputeVHAngles() || doSetAddress){
      if (options->computeVHAnglesOrder()!=3) getAngularBranches(tmpBranchList, 2, false);
      else{
        getAngularBranches(tmpBranchList, 3, false);
        getAngularBranches(tmpBranchList, 4, false);
      }
    }
    if (options->doComputeTTHAngles() || doSetAddress) getAngularBranches(tmpBranchList, 5, false);
  }
  for (unsigned int b=0; b<tmpBranchList.size(); b++){
    reserveBranch(tmpBranchList.at(b), BaseTree::bFloat, doSetAddress);
  }
}
void HVVTree::getAngularBranches(vector<string>& blist, const Int_t& prodFlag /* 0: Decay, 1: VBF, 2-4: VH, 5: TTH */, bool isGen){
  string strGen = "Gen";
  vector<string> strtmp;
  if (prodFlag!=5){
    strtmp.push_back("costhetastar");
    strtmp.push_back("costheta1");
    strtmp.push_back("costheta2");
    strtmp.push_back("Phi");
    strtmp.push_back("Phi1");
    if (prodFlag==1){
      strtmp.push_back("Q1");
      strtmp.push_back("Q2");
    }
  }
  else{
    // Masses
    strtmp.push_back("mTop1");
    strtmp.push_back("mW1");
    strtmp.push_back("mTop2");
    strtmp.push_back("mW2");

    // ttH system
    strtmp.push_back("costhetastar");
    strtmp.push_back("costheta1");
    strtmp.push_back("costheta2");
    strtmp.push_back("Phi");
    strtmp.push_back("Phi1");

    // tt system
    strtmp.push_back("costhetabb");
    strtmp.push_back("costhetaWW");
    strtmp.push_back("Phibb");
    strtmp.push_back("Phi1bb");

    // Wplus system
    strtmp.push_back("costhetaWplus");
    strtmp.push_back("PhiWplus");

    // Wminus system
    strtmp.push_back("costhetaWminus");
    strtmp.push_back("PhiWminus");
  }
  for (unsigned int b=0; b<strtmp.size(); b++){
    string varname = strtmp.at(b);
    if (isGen) varname.insert(0, strGen);
    if (prodFlag==1) varname.append("_VBF");
    else if (prodFlag==2) varname.append("_VH");
    else if (prodFlag==3) varname.append("_VHhadronic");
    else if (prodFlag==4) varname.append("_VHleptonic");
    else if (prodFlag==5) varname.append("_ttH");
    blist.push_back(varname);
  }
}

void HVVTree::fillMotherInfo(const MELAParticle* mother){
  if (options && options->processGenInfo() && mother){
    setVal("GenMotherPz", mother->z());
    setVal("GenMotherE", mother->t());
    setVal("GenMotherId", mother->id);
  }
}

void HVVTree::fillCandidate(MELACandidate* pH, bool isGen){
  if (!options) return;
  if ((!options->processGenInfo() && isGen) || (!options->processRecoInfo() && !isGen)) return;
  if (!pH) return;

  string varname;
  string strcore = "ZZ";
  if (isGen) strcore = "GenH";
  varname = strcore + "Mass"; setVal(varname, (pH ? pH->m() : 0.));
  varname = strcore + "Pt"; setVal(varname, (pH ? pH->pt() : 0));
  varname = strcore + "Pz"; setVal(varname, (pH ? pH->z() : 0));
  varname = strcore + "Phi"; setVal(varname, (pH ? pH->phi() : 0));

  fillCandidateDaughters(pH, isGen);
  fillDaughterProducts(pH, isGen);
  fillAssociatedInfo(pH, isGen);

  if (melaHelpers::melaHandle && pH){
    melaHelpers::melaHandle->setCurrentCandidate(pH);

    if (options->doComputeDecayAngles()) fillDecayAngles(isGen);
    if (options->doComputeVBFAngles()) fillVBFProductionAngles(isGen);
    if (options->doComputeVHAngles()) fillVHProductionAngles(isGen);
    if (options->doComputeTTHAngles()) fillTTHProductionAngles(isGen);
    fillMELAProbabilities(isGen); // Do it at the last step

    melaHelpers::melaHandle->resetInputEvent();
  }
}
void HVVTree::fillCandidateDaughters(MELACandidate* pH, bool isGen){
  string varname;
  string strcore;

  MELAParticle* pV1=(pH ? pH->getSortedV(0) : 0);
  MELAParticle* pV2=(pH ? pH->getSortedV(1) : 0);
  if (pH){
    if (pV1 && pV1->getMother(0)!=pH) pV1=0;
    if (pV2 && pV2->getMother(0)!=pH) pV2=0;
  }

  strcore = "Z1";
  if (isGen) strcore.insert(0, "Gen");
  varname = strcore + "Mass"; setVal(varname, (pV1 ? pV1->m() : 0));
  varname = strcore + "Pt"; setVal(varname, (pV1 ? pV1->pt() : 0));
  varname = strcore + "Eta"; setVal(varname, (pV1 ? (pV1->t()>0 ? pV1->eta() : 0) : 0));
  varname = strcore + "Phi"; setVal(varname, (pV1 ? (pV1->pt()>0 ? pV1->phi() : 0) : 0));
  strcore = "Z2";
  if (isGen) strcore.insert(0, "Gen");
  varname = strcore + "Mass"; setVal(varname, (pV2 ? pV2->m() : 0));
  varname = strcore + "Pt"; setVal(varname, (pV2 ? pV2->pt() : 0));
  varname = strcore + "Eta"; setVal(varname, (pV2 ? (pV2->t()>0 ? pV2->eta() : 0) : 0));
  varname = strcore + "Phi"; setVal(varname, (pV2 ? (pV2->pt()>0 ? pV2->phi() : 0) : 0));

  TLorentzVector nullVector(0, 0, 0, 0);
  TLorentzVector pZ1alt = (pH ? pH->getAlternativeVMomentum(0) : nullVector);
  TLorentzVector pZ2alt = (pH ? pH->getAlternativeVMomentum(1) : nullVector);

  strcore = "Za";
  if (isGen) strcore.insert(0, "Gen");
  varname = strcore + "Mass"; setVal(varname, (pH ? pZ1alt.M() : 0));
  varname = strcore + "Pt"; setVal(varname, (pH ? pZ1alt.Pt() : 0));
  varname = strcore + "Eta"; setVal(varname, (pH ? (pZ1alt.T()>0 ? pZ1alt.Eta() : 0) : 0));
  varname = strcore + "Phi"; setVal(varname, (pH ? (pZ1alt.Pt()>0 ? pZ1alt.Phi() : 0) : 0));
  strcore = "Zb";
  if (isGen) strcore.insert(0, "Gen");
  varname = strcore + "Mass"; setVal(varname, (pH ? pZ2alt.M() : 0));
  varname = strcore + "Pt"; setVal(varname, (pH ? pZ2alt.Pt() : 0));
  varname = strcore + "Eta"; setVal(varname, (pH ? (pZ2alt.T()>0 ? pZ2alt.Eta() : 0) : 0));
  varname = strcore + "Phi"; setVal(varname, (pH ? (pZ2alt.Pt()>0 ? pZ2alt.Phi() : 0) : 0));
}
void HVVTree::fillDaughterProducts(MELACandidate* pH, bool isGen){
  string varname;
  string strcore = "CandDau";
  if (isGen) strcore.insert(0, "Gen");

  int nDau = std::min((pH ? pH->getNSortedVs() : 0), 2);
  for (int v=0; v<nDau; v++){
    MELAParticle* intermediateV = pH->getSortedV(v);
    if (intermediateV && intermediateV->getMother(0)!=pH){ intermediateV = nullptr; nDau--; }
    if (v==nDau) break;
    int nVDau = (intermediateV ? intermediateV->getNDaughters() : 0);

    for (int d=0; d<2; d++){
      bool isNew = false;
      MELAParticle* lepton = (intermediateV ? intermediateV->getDaughter(d) : nullptr);

      if (!lepton){
        isNew=true;
        TLorentzVector pDaughter(0, 0, 0, 0);
        int idDaughter = -9000;
        if (intermediateV && nVDau==0 && d==0){
          pDaughter = intermediateV->p4;
          idDaughter = intermediateV->id;
        }
        else if (nDau==0 && v==0 && d==0){
          pDaughter = pH->p4;
          idDaughter = pH->id;
        }
        if (idDaughter!=-9000) lepton = new MELAParticle(idDaughter, pDaughter);
      }

      if (lepton){
        varname = strcore + "Mass"; setVal(varname, lepton->m());
        varname = strcore + "Pt"; setVal(varname, lepton->pt());
        varname = strcore + "Eta"; setVal(varname, (lepton->pt()>0. ? lepton->eta() : 0));
        varname = strcore + "Phi"; setVal(varname, (lepton->pt()>0. ? lepton->phi() : 0));
        varname = strcore + "Id"; setVal(varname, lepton->id);
      }
      if (isNew){ delete lepton; lepton=nullptr; }
    }
  }
}
void HVVTree::fillAssociatedInfo(MELACandidate* pH, bool isGen){
  vector<MELAParticle*> AssociatedParticle;
  vector<MELAParticle*> tmpAssociatedParticle;
  Float_t DijetMass=-1;
  Float_t DileptonMass=-1;
  Float_t DijetVVMass=-1;
  Float_t DileptonVVMass=-1;
  Float_t dRjet=0;
  Float_t dRlep=0;

  Int_t NAssociatedVs=0;
  vector<MELAParticle*> AssociatedV;
  vectorInt AssociatedV_Particle1Index;
  vectorInt AssociatedV_Particle2Index;

  Int_t NAssociatedTops=0;
  vector<MELATopCandidate_t*> AssociatedTop;
  vectorInt AssociatedTop_PartnerParticleIndex;
  vectorInt AssociatedTop_WFermionIndex;
  vectorInt AssociatedTop_WAntifermionIndex;

  if (pH){
    for (MELAParticle* apart:pH->getAssociatedJets()){ if (apart->passSelection) tmpAssociatedParticle.push_back(apart); }
    for (MELAParticle* apart:pH->getAssociatedLeptons()){ if (apart->passSelection) tmpAssociatedParticle.push_back(apart); }
    for (MELAParticle* apart:pH->getAssociatedNeutrinos()){ if (apart->passSelection) tmpAssociatedParticle.push_back(apart); }
    for (MELAParticle* apart:pH->getAssociatedPhotons()){ if (apart->passSelection) tmpAssociatedParticle.push_back(apart); }
  }

  while (!tmpAssociatedParticle.empty()){ // Re-sort all associated particles by leading pT (categories are individually sorted, but mixing categories loses this sorting)
    MELAParticle* tmpPart=nullptr;
    size_t pos=0;
    for (size_t el=0; el<tmpAssociatedParticle.size(); el++){
      if (!tmpPart){
        tmpPart = tmpAssociatedParticle.at(el); pos=el;
      }
      else if (tmpPart->pt()<tmpAssociatedParticle.at(el)->pt()){
        tmpPart = tmpAssociatedParticle.at(el); pos=el;
      }// Safer to do in two steps
    }
    AssociatedParticle.push_back(tmpPart);
    tmpAssociatedParticle.erase(tmpAssociatedParticle.begin()+pos);
  }

  pair<MELAParticle*, MELAParticle*> leadingPtJetPair(nullptr, nullptr);
  pair<MELAParticle*, MELAParticle*> leadingPtLeptonPair(nullptr, nullptr);
  for (MELAParticle* part:AssociatedParticle){
    if (PDGHelpers::isAJet(part->id)){
      if (!leadingPtJetPair.first) leadingPtJetPair.first = part;
      else if (!leadingPtJetPair.second) leadingPtJetPair.second = part;
    }
    if (PDGHelpers::isALepton(part->id)){
      if (!leadingPtLeptonPair.first) leadingPtLeptonPair.first = part;
      else if (!leadingPtLeptonPair.second) leadingPtLeptonPair.second = part;
    }
  }


  NAssociatedVs = (pH ? pH->getAssociatedSortedVs().size() : 0);
  if (pH){
    for (MELAParticle* pAV:pH->getAssociatedSortedVs()){
      bool doSkip = false;
      MELAParticle* avd1 = pAV->getDaughter(0);
      MELAParticle* avd2 = pAV->getDaughter(1);

      if (!pAV->passSelection || (avd1 && !avd1->passSelection) || (avd2 && !avd2->passSelection)) doSkip = true;
      if (doSkip) continue;

      AssociatedV.push_back(pAV);
      for (unsigned int aa=0; aa<AssociatedParticle.size(); aa++){
        if (avd1==AssociatedParticle.at(aa)) AssociatedV_Particle1Index.push_back(aa);
        else if (avd2==AssociatedParticle.at(aa)) AssociatedV_Particle2Index.push_back(aa);
      }
    }
  }

  NAssociatedTops = (pH ? pH->getNAssociatedTops() : 0);
  if (pH){
    for (MELATopCandidate_t* pTop:pH->getAssociatedTops()){
      bool doSkip = false;
      MELAParticle* tpp = pTop->getPartnerParticle();
      MELAParticle* tWf = pTop->getWFermion();
      MELAParticle* tWfb = pTop->getWAntifermion();

      if (!pTop->passSelection || (tpp && !tpp->passSelection) || (tWf && !tWf->passSelection) || (tWfb && !tWfb->passSelection)) doSkip = true;
      if (doSkip) continue;

      AssociatedTop.push_back(pTop);

      int i_tpp=-1;
      int i_tWf=-1;
      int i_tWfb=-1;
      for (unsigned int aa=0; aa<AssociatedParticle.size(); aa++){
        if (tpp && tpp==AssociatedParticle.at(aa)) i_tpp=aa;
        else if (tWf && tWf==AssociatedParticle.at(aa)) i_tWf=aa;
        else if (tWfb && tWfb==AssociatedParticle.at(aa)) i_tWfb=aa;
      }
      AssociatedTop_PartnerParticleIndex.push_back(i_tpp);
      AssociatedTop_WFermionIndex.push_back(i_tWf);
      AssociatedTop_WAntifermionIndex.push_back(i_tWfb);
    }
  }

  string varname;
  string strcore;

  if (leadingPtJetPair.first && leadingPtJetPair.second){
    DijetMass = (leadingPtJetPair.first->p4+leadingPtJetPair.second->p4).M();
    varname = "DijetMass";
    if (isGen) varname.insert(0, "Gen");
    setVal(varname, DijetMass);
    DijetVVMass = (pH->p4+leadingPtJetPair.first->p4+leadingPtJetPair.second->p4).M();
    varname = "DijetVVMass";
    if (isGen) varname.insert(0, "Gen");
    setVal(varname, DijetVVMass);
    dRjet = leadingPtJetPair.first->deltaR(leadingPtJetPair.second->p4);
    varname = "DRjet";
    if (isGen) varname.insert(0, "Gen");
    setVal(varname, dRjet);
  }
  if (leadingPtLeptonPair.first && leadingPtLeptonPair.second){
    DileptonMass = (leadingPtLeptonPair.first->p4+leadingPtLeptonPair.second->p4).M();
    varname = "DileptonMass";
    if (isGen) varname.insert(0, "Gen");
    setVal(varname, DileptonMass);
    DileptonVVMass = (pH->p4+leadingPtLeptonPair.first->p4+leadingPtLeptonPair.second->p4).M();
    varname = "DileptonVVMass";
    if (isGen) varname.insert(0, "Gen");
    setVal(varname, DileptonVVMass);
    dRlep = leadingPtLeptonPair.first->deltaR(leadingPtLeptonPair.second->p4);
    varname = "DRlepton";
    if (isGen) varname.insert(0, "Gen");
    setVal(varname, dRlep);
  }

  varname = "NAssociatedVs";
  if (isGen) varname.insert(1, "Gen");
  setVal(varname, NAssociatedVs);

  strcore = "AssociatedParticle";
  if (isGen) strcore.insert(0, "Gen");
  for (unsigned int aa=0; aa<AssociatedParticle.size(); aa++){
    varname = strcore + "Mass"; setVal(varname, AssociatedParticle.at(aa)->m());
    varname = strcore + "Pt"; setVal(varname, AssociatedParticle.at(aa)->pt());
    varname = strcore + "Eta"; setVal(varname, AssociatedParticle.at(aa)->eta());
    varname = strcore + "Phi"; setVal(varname, AssociatedParticle.at(aa)->phi());
    varname = strcore + "Id"; setVal(varname, AssociatedParticle.at(aa)->id);
  }

  strcore = "AssociatedV";
  if (isGen) strcore.insert(0, "Gen");
  for (int aa=0; aa<NAssociatedVs; aa++){
    varname = strcore + "Mass"; setVal(varname, AssociatedV.at(aa)->m());
    varname = strcore + "Pt"; setVal(varname, AssociatedV.at(aa)->pt());
    varname = strcore + "Eta"; setVal(varname, AssociatedV.at(aa)->eta());
    varname = strcore + "Phi"; setVal(varname, AssociatedV.at(aa)->phi());
    varname = strcore + "Id"; setVal(varname, AssociatedV.at(aa)->id);
    varname = strcore + "_Particle1Index"; setVal(varname, AssociatedV_Particle1Index.at(aa));
    varname = strcore + "_Particle2Index"; setVal(varname, AssociatedV_Particle2Index.at(aa));
  }

  strcore = "AssociatedTop";
  if (isGen) strcore.insert(0, "Gen");
  for (int aa=0; aa<NAssociatedTops; aa++){
    varname = strcore + "Mass"; setVal(varname, AssociatedTop.at(aa)->m());
    varname = strcore + "Pt"; setVal(varname, AssociatedTop.at(aa)->pt());
    varname = strcore + "Eta"; setVal(varname, AssociatedTop.at(aa)->eta());
    varname = strcore + "Phi"; setVal(varname, AssociatedTop.at(aa)->phi());
    varname = strcore + "Id"; setVal(varname, AssociatedTop.at(aa)->id);
    varname = strcore + "_PartnerParticleIndex"; setVal(varname, AssociatedTop_PartnerParticleIndex.at(aa));
    varname = strcore + "_WFermionIndex"; setVal(varname, AssociatedTop_WFermionIndex.at(aa));
    varname = strcore + "_WAntifermionIndex"; setVal(varname, AssociatedTop_WAntifermionIndex.at(aa));
  }
}

void HVVTree::fillDecayAngles(bool isGen){
  Float_t qH, m1, m2, costheta1=0, costheta2=0, Phi=0, costhetastar=0, Phi1=0;

  melaHelpers::melaHandle->computeDecayAngles(
    qH,
    m1,
    m2,
    costheta1,
    costheta2,
    Phi,
    costhetastar,
    Phi1
  );

  vector<string> varlist;
  getAngularBranches(varlist, 0, isGen);
  for (unsigned int b=0; b<varlist.size(); b++){
    string varname = varlist.at(b);
    if (varname.find("costhetastar")!=string::npos) setVal(varname, costhetastar);
    else if (varname.find("costheta1")!=string::npos) setVal(varname, costheta1);
    else if (varname.find("costheta2")!=string::npos) setVal(varname, costheta2);
    else if (varname.find("Phi1")!=string::npos) setVal(varname, Phi1);
    else if (varname.find("Phi")!=string::npos) setVal(varname, Phi);
    else cerr << "HVVTree::fillDecayAngles -> ERROR: Branch " << varname << " is invalid!" << endl;
  }
}
void HVVTree::fillVBFProductionAngles(bool isGen){
  Float_t costheta1=0, costheta2=0, Phi=0, costhetastar=0, Phi1=0;
  Float_t q1sq=-1, q2sq=-1;

  if (melaHelpers::melaHandle->getCurrentCandidate()){
    const std::vector<MELAParticle*>& ajets = melaHelpers::melaHandle->getCurrentCandidate()->getAssociatedJets();

    unsigned int nselajets=0;
    for (auto const* jet:ajets) if (jet && jet->passSelection) nselajets++;
    if (nselajets>=2){
      melaHelpers::melaHandle->computeVBFAngles(
        q1sq,
        q2sq,
        costheta1,
        costheta2,
        Phi,
        costhetastar,
        Phi1
      );
      q1sq = (q1sq>=0 ? sqrt(q1sq) : -sqrt(-q1sq));
      q2sq = (q2sq>=0 ? sqrt(q2sq) : -sqrt(-q2sq));
    }
  }
  vector<string> varlist;
  getAngularBranches(varlist, 1, isGen);
  for (unsigned int b=0; b<varlist.size(); b++){
    string varname = varlist.at(b);
    if (varname.find("costhetastar_VBF")!=string::npos) setVal(varname, costhetastar);
    else if (varname.find("costheta1_VBF")!=string::npos) setVal(varname, costheta1);
    else if (varname.find("costheta2_VBF")!=string::npos) setVal(varname, costheta2);
    else if (varname.find("Phi1_VBF")!=string::npos) setVal(varname, Phi1);
    else if (varname.find("Phi_VBF")!=string::npos) setVal(varname, Phi);
    else if (varname.find("Q1_VBF")!=string::npos) setVal(varname, q1sq);
    else if (varname.find("Q2_VBF")!=string::npos) setVal(varname, q2sq);
    else cerr << "HVVTree::fillVBFProductionAngles -> ERROR: Branch " << varname << " is invalid!" << endl;
  }
}
void HVVTree::fillVHProductionAngles(bool isGen){
  Float_t mVstar_hadronic=-1, mV_hadronic=-1, costheta1_hadronic=0, costheta2_hadronic=0, Phi_hadronic=0, costhetastar_hadronic=0, Phi1_hadronic=0;
  Float_t mVstar_leptonic=-1, mV_leptonic=-1, costheta1_leptonic=0, costheta2_leptonic=0, Phi_leptonic=0, costhetastar_leptonic=0, Phi1_leptonic=0;


  if (options->computeVHAnglesOrder()!=2){ // Vjj
    TVar::Production Vprod = melaHelpers::getFirstAssociatedHadronicVProduction(melaHelpers::melaHandle->getCurrentCandidate());
    if (Vprod!=TVar::nProductions){
      melaHelpers::melaHandle->setProcess(TVar::HSMHiggs, TVar::JHUGen, Vprod);
      melaHelpers::melaHandle->computeVHAngles(
        mVstar_hadronic,
        mV_hadronic,
        costheta1_hadronic,
        costheta2_hadronic,
        Phi_hadronic,
        costhetastar_hadronic,
        Phi1_hadronic
      );
    }
  }
  if (options->computeVHAnglesOrder()!=1){ // Vll
    TVar::Production Vprod = melaHelpers::getFirstAssociatedLeptonicVProduction(melaHelpers::melaHandle->getCurrentCandidate());
    if (Vprod!=TVar::nProductions){
      melaHelpers::melaHandle->setProcess(TVar::HSMHiggs, TVar::JHUGen, Vprod);
      melaHelpers::melaHandle->computeVHAngles(
        mVstar_leptonic,
        mV_leptonic,
        costheta1_leptonic,
        costheta2_leptonic,
        Phi_leptonic,
        costhetastar_leptonic,
        Phi1_leptonic
      );
    }
  }

  vector<string> varlist;
  if (options->computeVHAnglesOrder()!=3) getAngularBranches(varlist, 2, isGen);
  else{
    getAngularBranches(varlist, 3, isGen);
    getAngularBranches(varlist, 4, isGen);
  }
  for (unsigned int b=0; b<varlist.size(); b++){
    string varname = varlist.at(b);
    if (options->computeVHAnglesOrder()==1){
      if (varname.find("costhetastar_VH")!=string::npos) setVal(varname, costhetastar_hadronic);
      else if (varname.find("costheta1_VH")!=string::npos) setVal(varname, costheta1_hadronic);
      else if (varname.find("costheta2_VH")!=string::npos) setVal(varname, costheta2_hadronic);
      else if (varname.find("Phi1_VH")!=string::npos) setVal(varname, Phi1_hadronic);
      else if (varname.find("Phi_VH")!=string::npos) setVal(varname, Phi_hadronic);
      else if (varname.find("mVstar_VH")!=string::npos) setVal(varname, mVstar_hadronic);
      else if (varname.find("mV_VH")!=string::npos) setVal(varname, mV_hadronic);
      else cerr << "HVVTree::fillVHProductionAngles -> ERROR: Branch " << varname << " is invalid!" << endl;
    }
    if (options->computeVHAnglesOrder()==2){
      if (varname.find("costhetastar_VH")!=string::npos) setVal(varname, costhetastar_leptonic);
      else if (varname.find("costheta1_VH")!=string::npos) setVal(varname, costheta1_leptonic);
      else if (varname.find("costheta2_VH")!=string::npos) setVal(varname, costheta2_leptonic);
      else if (varname.find("Phi1_VH")!=string::npos) setVal(varname, Phi1_leptonic);
      else if (varname.find("Phi_VH")!=string::npos) setVal(varname, Phi_leptonic);
      else if (varname.find("mVstar_VH")!=string::npos) setVal(varname, mVstar_leptonic);
      else if (varname.find("mV_VH")!=string::npos) setVal(varname, mV_leptonic);
      else cerr << "HVVTree::fillVHProductionAngles -> ERROR: Branch " << varname << " is invalid!" << endl;
    }
    if (options->computeVHAnglesOrder()==3){
      if (varname.find("costhetastar_VHhadronic")!=string::npos) setVal(varname, costhetastar_hadronic);
      else if (varname.find("costheta1_VHhadronic")!=string::npos) setVal(varname, costheta1_hadronic);
      else if (varname.find("costheta2_VHhadronic")!=string::npos) setVal(varname, costheta2_hadronic);
      else if (varname.find("Phi1_VHhadronic")!=string::npos) setVal(varname, Phi1_hadronic);
      else if (varname.find("Phi_VHhadronic")!=string::npos) setVal(varname, Phi_hadronic);
      else if (varname.find("mVstar_VHhadronic")!=string::npos) setVal(varname, mVstar_hadronic);
      else if (varname.find("mV_VHhadronic")!=string::npos) setVal(varname, mV_hadronic);
      else if (varname.find("costhetastar_VHleptonic")!=string::npos) setVal(varname, costhetastar_leptonic);
      else if (varname.find("costheta1_VHleptonic")!=string::npos) setVal(varname, costheta1_leptonic);
      else if (varname.find("costheta2_VHleptonic")!=string::npos) setVal(varname, costheta2_leptonic);
      else if (varname.find("Phi1_VHleptonic")!=string::npos) setVal(varname, Phi1_leptonic);
      else if (varname.find("Phi_VHleptonic")!=string::npos) setVal(varname, Phi_leptonic);
      else if (varname.find("mVstar_VHleptonic")!=string::npos) setVal(varname, mVstar_leptonic);
      else if (varname.find("mV_VHleptonic")!=string::npos) setVal(varname, mV_leptonic);
      else cerr << "HVVTree::fillVHProductionAngles -> ERROR: Branch " << varname << " is invalid!" << endl;
    }
  }
}
void HVVTree::fillTTHProductionAngles(bool isGen){
  // Masses
  float mT1=0;
  float mW1=0;
  float mT2=0;
  float mW2=0;

  // ttH system
  float hs=0;
  float h1=0;
  float h2=0;
  float Phi=0;
  float Phi1=0;

  // tt system
  float hbb=0;
  float hWW=0;
  float Phibb=0;
  float Phi1bb=0;

  // Wplus system
  float hWplusf=0;
  float PhiWplusf=0;

  // Wminus system
  float hWminusf=0;
  float PhiWminusf=0;


  if (melaHelpers::melaHandle->getCurrentCandidate()){
    melaHelpers::melaHandle->computeTTHAngles(
      1,

      mT1, mW1,
      mT2, mW2,

      // ttH system
      h1, h2, Phi, hs, Phi1,

      // tt system
      hbb, hWW, Phibb, Phi1bb,

      // Wplus system
      hWplusf, PhiWplusf,

      // Wminus system
      hWminusf, PhiWminusf
    );
  }
  vector<string> varlist;
  getAngularBranches(varlist, 5, isGen);
  for (unsigned int b=0; b<varlist.size(); b++){
    string varname = varlist.at(b);
    if (varname.find("costhetastar_ttH")!=string::npos) setVal(varname, hs);
    else if (varname.find("costheta1_ttH")!=string::npos) setVal(varname, h1);
    else if (varname.find("costheta2_ttH")!=string::npos) setVal(varname, h2);
    else if (varname.find("Phi1_ttH")!=string::npos) setVal(varname, Phi1);
    else if (varname.find("Phi_ttH")!=string::npos) setVal(varname, Phi);

    else if (varname.find("costhetabb_ttH")!=string::npos) setVal(varname, hbb);
    else if (varname.find("costhetaWW_ttH")!=string::npos) setVal(varname, hWW);
    else if (varname.find("Phi1bb_ttH")!=string::npos) setVal(varname, Phi1bb);
    else if (varname.find("Phibb_ttH")!=string::npos) setVal(varname, Phibb);

    else if (varname.find("costhetaWplus_ttH")!=string::npos) setVal(varname, hWplusf);
    else if (varname.find("PhiWplus_ttH")!=string::npos) setVal(varname, PhiWplusf);

    else if (varname.find("costhetaWminus_ttH")!=string::npos) setVal(varname, hWminusf);
    else if (varname.find("PhiWminus_ttH")!=string::npos) setVal(varname, PhiWminusf);

    else if (varname.find("mTop1_ttH")!=string::npos) setVal(varname, mT1);
    else if (varname.find("mW1_ttH")!=string::npos) setVal(varname, mW1);
    else if (varname.find("mTop2_ttH")!=string::npos) setVal(varname, mT2);
    else if (varname.find("mW2_ttH")!=string::npos) setVal(varname, mW2);
    else cerr << "HVVTree::fillTTHProductionAngles -> ERROR: Branch " << varname << " is invalid!" << endl;
  }
}
void HVVTree::fillMELAProbabilities(bool isGen){
  computeMELABranches(isGen);
  pushMELABranches(isGen);
}


void HVVTree::fillXsec(const Float_t& val, const Float_t& err){
  setVal("xsec", val);
  setVal("xsecerr", err);
}
void HVVTree::fillEventVariables(const Float_t& weight, const Int_t& passSelection){
  setVal("MC_weight", weight);
  if (options && options->processRecoInfo()) setVal("isSelected", passSelection);
}


void HVVTree::buildMELABranches(bool doSetAddress){
  std::vector<std::string> const& lheMElist = options->getLHEMEList();
  std::vector<std::string> const& recoMElist = options->getRecoMEList();

  if (!melaHelpers::melaHandle){
    // If the MELA object exists, there is no need to prompt binding because the branches will be recomputed.
    /***************************/
    /***** LHE ME BRANCHES *****/
    /***************************/
    for (auto const& strLHEME:lheMElist){
      MELAOptionParser* me_opt = new MELAOptionParser(strLHEME);
      if (strLHEME.find("Copy")!=string::npos) lheme_copyopts.push_back(me_opt);
      else lheme_originalopts.push_back(me_opt);
    }
    // Resolve original options
    for (MELAOptionParser* me_opt:lheme_originalopts) bookMELABranches(me_opt, nullptr, false);
    // Resolve copy options
    for (MELAOptionParser* me_opt:lheme_copyopts){
      MELAOptionParser* original_opt=nullptr;
      // Find the original options
      for (MELAOptionParser* tmp_opt:lheme_originalopts){
        if (me_opt->testCopyAlias(tmp_opt->getAlias())){
          original_opt = tmp_opt;
          break;
        }
      }
      if (!original_opt) continue;
      else me_opt->pickOriginalOptions(original_opt);
      bookMELABranches(me_opt, nullptr, false);
    }

    /****************************/
    /***** RECO ME BRANCHES *****/
    /****************************/
    for (auto const& strRecoME:recoMElist){
      MELAOptionParser* me_opt = new MELAOptionParser(strRecoME);
      if (strRecoME.find("Copy")!=string::npos) recome_copyopts.push_back(me_opt);
      else recome_originalopts.push_back(me_opt);
    }
    // Resolve original options
    for (MELAOptionParser* me_opt:recome_originalopts) bookMELABranches(me_opt, nullptr, false);
    // Resolve copy options
    for (MELAOptionParser* me_opt:recome_copyopts){
      MELAOptionParser* original_opt=nullptr;
      // Find the original options
      for (MELAOptionParser* tmp_opt:recome_originalopts){
        if (me_opt->testCopyAlias(tmp_opt->getAlias())){
          original_opt = tmp_opt;
          break;
        }
      }
      if (!original_opt) continue;
      else me_opt->pickOriginalOptions(original_opt);
      bookMELABranches(me_opt, nullptr, false);
    }

    // Iterate over the MELABranches and do some address set/unsetting.
    for (MELABranch* br:recome_branches){
      TString const& brname = br->bname;
      TBranch* root_br = br->getBranch();
      root_br->ResetAddress();
      this->reserveBranch(brname.Data(), BaseTree::bFloat, doSetAddress);
    }
    for (MELABranch* br:lheme_branches){
      TString const& brname = br->bname;
      TBranch* root_br = br->getBranch();
      root_br->ResetAddress();
      this->reserveBranch(brname.Data(), BaseTree::bFloat, doSetAddress);
    }
    clearMELABranches();
  }
  else if (!doSetAddress){
    /***************************/
    /***** LHE ME BRANCHES *****/
    /***************************/
    for (auto const& strLHEME:lheMElist){
      MELAOptionParser* me_opt;
      // First find out if the option has a copy specification
      // These copy options will be evaulated in a separate loop
      if (strLHEME.find("Copy")!=string::npos){
        me_opt = new MELAOptionParser(strLHEME);
        lheme_copyopts.push_back(me_opt);
        continue;
      }

      // Create a hypothesis for each option
      MELAHypothesis* me_hypo = new MELAHypothesis(melaHelpers::melaHandle, strLHEME);
      lheme_units.push_back(me_hypo);

      me_opt = me_hypo->getOption();
      if (me_opt->isAliased()) lheme_aliased_units.push_back(me_hypo);

      // Create a computation for each hypothesis
      MELAComputation* me_computer = new MELAComputation(me_hypo);
      lheme_computers.push_back(me_computer);

      // Add the computation to a named cluster to keep track of JECUp/JECDn, or for best-pWH_SM Lep_WH computations
      GMECHelperFunctions::addToMELACluster(me_computer, lheme_clusters);

      this->bookMELABranches(me_opt, me_computer, false);
    }
    // Resolve copy options
    for (MELAOptionParser* me_opt:lheme_copyopts){
      MELAHypothesis* original_hypo=nullptr;
      MELAOptionParser* original_opt=nullptr;
      // Find the original options
      for (auto* me_aliased_unit:lheme_aliased_units){
        if (me_opt->testCopyAlias(me_aliased_unit->getOption()->getAlias())){
          original_hypo = me_aliased_unit;
          original_opt = original_hypo->getOption();
          break;
        }
      }
      if (!original_opt) continue;
      else me_opt->pickOriginalOptions(original_opt);
      // Create a new computation for the copy options
      MELAComputation* me_computer = new MELAComputation(original_hypo);
      me_computer->setOption(me_opt);
      lheme_computers.push_back(me_computer);

      // The rest is the same story...
      // Add the computation to a named cluster to keep track of JECUp/JECDn, or for best-pWH_SM Lep_WH computations
      GMECHelperFunctions::addToMELACluster(me_computer, lheme_clusters);

      // Create the necessary branches for each computation
      // Notice that no tree is passed, so no TBranches are created.
      this->bookMELABranches(me_opt, me_computer, true);
    }
    // Loop over the computations to add any contingencies to aliased hypotheses
    for (auto& me_computer:lheme_computers) me_computer->addContingencies(lheme_aliased_units);

    /****************************/
    /***** RECO ME BRANCHES *****/
    /****************************/
    for (auto const& strRecoME:recoMElist){
      MELAOptionParser* me_opt;
      // First find out if the option has a copy specification
      // These copy options will be evaulated in a separate loop
      if (strRecoME.find("Copy")!=string::npos){
        me_opt = new MELAOptionParser(strRecoME);
        recome_copyopts.push_back(me_opt);
        continue;
      }

      // Create a hypothesis for each option
      MELAHypothesis* me_hypo = new MELAHypothesis(melaHelpers::melaHandle, strRecoME);
      recome_units.push_back(me_hypo);

      me_opt = me_hypo->getOption();
      if (me_opt->isAliased()) recome_aliased_units.push_back(me_hypo);

      // Create a computation for each hypothesis
      MELAComputation* me_computer = new MELAComputation(me_hypo);
      recome_computers.push_back(me_computer);

      // Add the computation to a named cluster to keep track of JECUp/JECDn, or for best-pWH_SM Lep_WH computations
      GMECHelperFunctions::addToMELACluster(me_computer, recome_clusters);

      this->bookMELABranches(me_opt, me_computer, false);
    }
    // Resolve copy options
    for (MELAOptionParser* me_opt:recome_copyopts){
      MELAHypothesis* original_hypo=nullptr;
      MELAOptionParser* original_opt=nullptr;
      // Find the original options
      for (auto* me_aliased_unit:recome_aliased_units){
        if (me_opt->testCopyAlias(me_aliased_unit->getOption()->getAlias())){
          original_hypo = me_aliased_unit;
          original_opt = original_hypo->getOption();
          break;
        }
      }
      if (!original_opt) continue;
      else me_opt->pickOriginalOptions(original_opt);
      // Create a new computation for the copy options
      MELAComputation* me_computer = new MELAComputation(original_hypo);
      me_computer->setOption(me_opt);
      recome_computers.push_back(me_computer);

      // The rest is the same story...
      // Add the computation to a named cluster to keep track of JECUp/JECDn, or for best-pWH_SM Lep_WH computations
      GMECHelperFunctions::addToMELACluster(me_computer, recome_clusters);

      // Create the necessary branches for each computation
      // Notice that no tree is passed, so no TBranches are created.
      this->bookMELABranches(me_opt, me_computer, true);
    }
    // Loop over the computations to add any contingencies to aliased hypotheses
    for (auto& me_computer:recome_computers) me_computer->addContingencies(recome_aliased_units);

    if (!lheme_clusters.empty()){
      cout << "HVVTree::buildMELABranches: LHE ME clusters:" << endl;
      for (auto const* me_cluster:lheme_clusters){
        cout << "\t- Cluster " << me_cluster->getName() << " has " << me_cluster->getComputations()->size() << " computations registered." << endl;
      }
    }
    if (!recome_clusters.empty()){
      cout << "HVVTree::buildMELABranches: Reco ME clusters:" << endl;
      for (auto const* me_cluster:recome_clusters){
        cout << "\t- Cluster " << me_cluster->getName() << " has " << me_cluster->getComputations()->size() << " computations registered." << endl;
      }
    }
  }
}
void HVVTree::bookMELABranches(MELAOptionParser* me_opt, MELAComputation* computer, bool doCopy){
  if (!me_opt){
    cerr << "HVVTree::bookMELABranches: Did not receive a valid me_opt. Something went wrong." << endl;
    throw std::exception();
  }

  std::vector<MELABranch*>* me_branches = (me_opt->isGen() ? &(this->lheme_branches) : &(this->recome_branches));

  vector<TTree*> trees;
  if (!doCopy){
    trees.push_back(hvvtree);
    //if (me_opt->isGen()) trees.push_back(_failedTree);
  }
  else trees.push_back(nullptr); // Copy branches should not contain a tree reference

  if (me_opt->doBranch()){
    for (auto* tree:trees){
      string basename = me_opt->getName();
      if (me_opt->isGen()) basename = string("Gen_") + basename;
      MELABranch* tmpbranch;
      Float_t defVal=1.;
      if (me_opt->hasPAux()){
        tmpbranch = new MELABranch(
          tree, TString((string("pAux_") + basename).c_str()),
          defVal, computer
        );
        me_branches->push_back(tmpbranch);
      }
      if (me_opt->hasPConst()){
        tmpbranch = new MELABranch(
          tree, TString((string("pConst_") + basename).c_str()),
          defVal, computer
        );
        me_branches->push_back(tmpbranch);
      }
      defVal = me_opt->getDefaultME();
      tmpbranch = new MELABranch(
        tree, TString((string("p_") + basename).c_str()),
        defVal, computer
      );
      me_branches->push_back(tmpbranch);
      cout << "HVVTree::bookMELABranches: Constructed branch with base name " << basename << endl;
    }
  }
}
void HVVTree::clearMELABranches(){
#define CLEAR_MELA_BRANCHES_CMD(thelist) \
for (auto*& v:thelist) delete v; \
thelist.clear();

  CLEAR_MELA_BRANCHES_CMD(lheme_branches);
  CLEAR_MELA_BRANCHES_CMD(lheme_clusters);
  CLEAR_MELA_BRANCHES_CMD(lheme_computers);
  CLEAR_MELA_BRANCHES_CMD(lheme_copyopts);
  CLEAR_MELA_BRANCHES_CMD(lheme_originalopts);
  // Do not delete me_aliased_units. They are deleted together with me_units.
  CLEAR_MELA_BRANCHES_CMD(lheme_units);

  CLEAR_MELA_BRANCHES_CMD(recome_branches);
  CLEAR_MELA_BRANCHES_CMD(recome_clusters);
  CLEAR_MELA_BRANCHES_CMD(recome_computers);
  CLEAR_MELA_BRANCHES_CMD(recome_copyopts);
  CLEAR_MELA_BRANCHES_CMD(recome_originalopts);
  // Do not delete me_aliased_units. They are deleted together with me_units.
  CLEAR_MELA_BRANCHES_CMD(recome_units);

#undef CLEAR_MELA_BRANCHES_CMD
}


void HVVTree::computeMELABranches(bool isGen){
  if (!melaHelpers::melaHandle) return;
  // Sequantial computation
  updateMELAClusters_Common("Common", isGen);
  updateMELAClusters_J1JEC("J1JECNominal", isGen); updateMELAClusters_J1JEC("J1JECUp", isGen); updateMELAClusters_J1JEC("J1JECDn", isGen);
  updateMELAClusters_J2JEC("J2JECNominal", isGen); updateMELAClusters_J2JEC("J2JECUp", isGen); updateMELAClusters_J2JEC("J2JECDn", isGen);
  updateMELAClusters_LepWH("LepWH", isGen);
  updateMELAClusters_LepZH("LepZH", isGen);
  updateMELAClusters_NoInitialQ("NoInitialQ", isGen);
  updateMELAClusters_NoInitialG("NoInitialG", isGen);
  updateMELAClusters_BestLOAssociatedZ("BestLOAssociatedZ", isGen);
  updateMELAClusters_BestLOAssociatedW("BestLOAssociatedW", isGen);
  updateMELAClusters_BestLOAssociatedVBF("BestLOAssociatedVBF", isGen);
  updateMELAClusters_BestNLOVHApproximation("BestNLOZHApproximation", isGen);
  updateMELAClusters_BestNLOVHApproximation("BestNLOWHApproximation", isGen);
  updateMELAClusters_BestNLOVBFApproximation("BestNLOVBFApproximation", isGen);
  updateMELAClusters_NoAssociatedG("NoAssociatedG", isGen);
  updateMELAClusters_NoInitialGNoAssociatedG("NoInitialGNoAssociatedG", isGen);
  // Reverse sequence
  updateMELAClusters_NoInitialGNoAssociatedG("NoInitialGNoAssociatedGLast", isGen);
  updateMELAClusters_NoAssociatedG("NoAssociatedGLast", isGen);
  updateMELAClusters_BestNLOVBFApproximation("BestNLOVBFApproximationLast", isGen);
  updateMELAClusters_BestNLOVHApproximation("BestNLOWHApproximationLast", isGen);
  updateMELAClusters_BestNLOVHApproximation("BestNLOZHApproximationLast", isGen);
  updateMELAClusters_BestLOAssociatedVBF("BestLOAssociatedVBFLast", isGen);
  updateMELAClusters_BestLOAssociatedW("BestLOAssociatedWLast", isGen);
  updateMELAClusters_BestLOAssociatedZ("BestLOAssociatedZLast", isGen);
  updateMELAClusters_NoInitialG("NoInitialGLast", isGen);
  updateMELAClusters_NoInitialQ("NoInitialQLast", isGen);
  updateMELAClusters_LepZH("LepZHLast", isGen);
  updateMELAClusters_LepWH("LepWHLast", isGen);
  updateMELAClusters_J2JEC("J2JECNominalLast", isGen); updateMELAClusters_J2JEC("J2JECUpLast", isGen); updateMELAClusters_J2JEC("J2JECDnLast", isGen);
  updateMELAClusters_J1JEC("J1JECNominalLast", isGen); updateMELAClusters_J1JEC("J1JECUpLast", isGen); updateMELAClusters_J1JEC("J1JECDnLast", isGen);
  updateMELAClusters_Common("CommonLast", isGen);
}
void HVVTree::pushMELABranches(bool isGen){
  std::vector<MELABranch*>& me_branches = (isGen ? lheme_branches : recome_branches);
  std::vector<MELACluster*>& me_clusters = (isGen ? lheme_clusters : recome_clusters);
  // Pull + push...
  for (MELABranch* v:me_branches) v->setVal();
  // ...then reset
  for (MELACluster* v:me_clusters) v->reset();
}
// Common ME computations that do not manipulate the LHE candidate
void HVVTree::updateMELAClusters_Common(const string clustertype, bool isGen){
  MELACandidate* melaCand = melaHelpers::melaHandle->getCurrentCandidate();
  if (!melaCand) return;

  std::vector<MELACluster*>& me_clusters = (isGen ? lheme_clusters : recome_clusters);
  for (MELACluster* theCluster:me_clusters){
    if (theCluster->getName()==clustertype){
      // Re-compute all related hypotheses first...
      theCluster->computeAll();
      // ...then update the cluster
      theCluster->update();
    }
  }
}
// Common ME computations for JECNominal, Up and Down variations, case where ME requires 1 jet
void HVVTree::updateMELAClusters_J1JEC(const string clustertype, bool isGen){
  constexpr int njecnum = 3;
  const int jecnumsel = -1
    + int(clustertype=="J1JECNominal")*1
    + int(clustertype=="J1JECUp")*2
    + int(clustertype=="J1JECDn")*3;
  if (jecnumsel<0) return;

  //cout << "Begin HVVTree::updateMELAClusters_J1JEC(" << clustertype << "," << isGen << ")" << endl;

  // First determine if any of the candidates has only one jet
  bool doSkip=true;
  int nMelaStored = melaHelpers::melaHandle->getNCandidates();
  bool setJECcand = (nMelaStored == njecnum);
  if (!(setJECcand || nMelaStored==0)) return;

  for (int jecnum=0; jecnum<njecnum; jecnum++){
    if (setJECcand) melaHelpers::melaHandle->setCurrentCandidateFromIndex(jecnum);
    else if (jecnum>0) break;

    MELACandidate* melaCand = melaHelpers::melaHandle->getCurrentCandidate();
    if (!melaCand) continue;

    unsigned int nGoodJets=melaCand->getNAssociatedJets();
    doSkip = doSkip && (nGoodJets!=1);
  }
  if (doSkip) return; // If none of the candidates have exactly 1 jet, skip the computations

  std::vector<MELACluster*>& me_clusters = (isGen ? lheme_clusters : recome_clusters);

  for (int jecnum=0; jecnum<njecnum; jecnum++){
    if (jecnum!=jecnumsel) continue;

    if (setJECcand) melaHelpers::melaHandle->setCurrentCandidateFromIndex(jecnum);
    else if (jecnum>0) break;

    MELACandidate* melaCand = melaHelpers::melaHandle->getCurrentCandidate();
    if (!melaCand) continue;

    const unsigned int nGoodJets=std::min(1, melaCand->getNAssociatedJets());
    for (unsigned int firstjet = 0; firstjet < nGoodJets; firstjet++){ // Loop over first jet

      for (unsigned int disableJet=0; disableJet<nGoodJets; disableJet++) melaCand->getAssociatedJet(disableJet)->setSelected(disableJet==firstjet); // Disable the other jets

      for (MELACluster* theCluster:me_clusters){
        if (theCluster->getName()==clustertype){
          // Re-compute all related hypotheses first...
          theCluster->computeAll();
          // ...then force an update the cluster
          theCluster->forceUpdate(); // MELAComputation::testMaximizationCache still prevents update for a second time if there are no contingencies.
        }
      } // End loop over clusters

    } // End loop over first jet
    // Turn associated jets back on
    for (unsigned int disableJet=0; disableJet<nGoodJets; disableJet++) melaCand->getAssociatedJet(disableJet)->setSelected(true); // Turn all jets back on
  } // End jecnum loop

  //cout << "End HVVTree::updateMELAClusters_J1JEC(" << clustertype << "," << isGen << ")" << endl;
}
// Common ME computations for JECNominal, Up and Down variations, case where ME requires 2 jets
void HVVTree::updateMELAClusters_J2JEC(const string clustertype, bool isGen){
  constexpr int njecnum = 3;
  const int jecnumsel = -1
    + int(clustertype=="J2JECNominal")*1
    + int(clustertype=="J2JECUp")*2
    + int(clustertype=="J2JECDn")*3;
  if (jecnumsel<0) return;

  //cout << "Begin HVVTree::updateMELAClusters_J2JEC(" << clustertype << "," << isGen << ")" << endl;

  int nMelaStored = melaHelpers::melaHandle->getNCandidates();
  bool setJECcand = (nMelaStored == njecnum);
  if (!(setJECcand || nMelaStored==0)) return;

  std::vector<MELACluster*>& me_clusters = (isGen ? lheme_clusters : recome_clusters);

  for (int jecnum=0; jecnum<njecnum; jecnum++){
    if (jecnum!=jecnumsel) continue;

    if (setJECcand) melaHelpers::melaHandle->setCurrentCandidateFromIndex(jecnum);
    else if (jecnum>0) break;

    MELACandidate* melaCand = melaHelpers::melaHandle->getCurrentCandidate();
    if (!melaCand) continue;

    const unsigned int nGoodJets=melaCand->getNAssociatedJets();
    for (unsigned int firstjet = 0; firstjet < nGoodJets; firstjet++){ // Loop over first jet
      for (unsigned int secondjet = firstjet+1; secondjet < nGoodJets; secondjet++){ // Loop over second jet
        // Disable jets and tops
        for (unsigned int disableJet=0; disableJet<nGoodJets; disableJet++) melaCand->getAssociatedJet(disableJet)->setSelected((disableJet==firstjet || disableJet==secondjet)); // Disable the other jets
        unsigned int nDisabledStableTops=0;
        for (MELATopCandidate_t* einTop:melaCand->getAssociatedTops()){
          if (einTop->getNDaughters()==3) einTop->setSelected(false); // All unstable tops are disabled in the loop for jets (where "jet"=="stable top") since we are looping over jecnum
          else{
            einTop->setSelected((nDisabledStableTops==firstjet || nDisabledStableTops==secondjet)); // Disable the other stable tops
            nDisabledStableTops++;
          }
        }

        for (MELACluster* theCluster:me_clusters){
          if (theCluster->getName()==clustertype){
            // Re-compute all related hypotheses first...
            theCluster->computeAll();
            // ...then force an update the cluster
            theCluster->forceUpdate(); // MELAComputation::testMaximizationCache still prevents update for a second time if there are no contingencies.
          }
        } // End loop over clusters

      } // End loop over second jet
    } // End loop over first jet
    // Turn associated jets/tops back on
    for (unsigned int disableJet=0; disableJet<nGoodJets; disableJet++) melaCand->getAssociatedJet(disableJet)->setSelected(true); // Turn all jets back on
    for (MELATopCandidate_t* einTop:melaCand->getAssociatedTops()) einTop->setSelected(true); // Turn all tops back on
  } // End jecnum loop

  //cout << "End HVVTree::updateMELAClusters_J2JEC(" << clustertype << "," << isGen << ")" << endl;
}
// Common ME computations for leptonic WH: Loops over possible fake neutrinos
void HVVTree::updateMELAClusters_LepWH(const string clustertype, bool isGen){
  int nMelaStored = melaHelpers::melaHandle->getNCandidates();
  bool setLepHypoCand = (nMelaStored == 1);
  if (!(setLepHypoCand || nMelaStored==0)) return;

  if (setLepHypoCand) melaHelpers::melaHandle->setCurrentCandidateFromIndex(0);
  MELACandidate* melaCand = melaHelpers::melaHandle->getCurrentCandidate();
  if (!melaCand) return;

  std::vector<MELACluster*>& me_clusters = (isGen ? lheme_clusters : recome_clusters);

  int nNeutrinos = melaCand->getNAssociatedNeutrinos();
  for (int inu=0; inu<nNeutrinos; inu++){
    // Notice: Looping over Ws does not make much sense unless you have more than one lepton since the fake neutrino is already calculated from the available lepton with W mass constraint.
    // Such a loop over Ws only makes sense if there are more than one lepton in the event, but in that case, it still does not make sense to cross-match neutrinos and leptons.
    for (int disableNu=0; disableNu<nNeutrinos; disableNu++) melaCand->getAssociatedNeutrino(disableNu)->setSelected(disableNu==inu); // Disable all neutrinos other than index==inu
    for (MELACluster* theCluster:me_clusters){
      if (theCluster->getName()==clustertype){
        // Re-compute all related hypotheses first...
        theCluster->computeAll();
        // ...then force an update the cluster
        theCluster->forceUpdate(); // MELAComputation::testMaximizationCache still prevents update for a second time if there are no contingencies.
      }
    } // End loop over clusters
  } // End loop over possible neutrinos
  // Re-enable all neutrinos
  for (int disableNu=0; disableNu<nNeutrinos; disableNu++) melaCand->getAssociatedNeutrino(disableNu)->setSelected(true);
}
// Common ME computations for leptonic ZH: Picks best Z3
void HVVTree::updateMELAClusters_LepZH(const string clustertype, bool isGen){
  int nMelaStored = melaHelpers::melaHandle->getNCandidates();
  bool setLepHypoCand = (nMelaStored == 1);
  if (!(setLepHypoCand || nMelaStored==0)) return;

  if (setLepHypoCand) melaHelpers::melaHandle->setCurrentCandidateFromIndex(0);
  MELACandidate* melaCand = melaHelpers::melaHandle->getCurrentCandidate();
  if (!melaCand) return;

  std::vector<MELACluster*>& me_clusters = (isGen ? lheme_clusters : recome_clusters);

  std::vector<MELAParticle*> associatedVs = melaCand->getAssociatedSortedVs();
  const unsigned int nSortedVs = associatedVs.size();
  constexpr unsigned int iSortedVstart=0;
  double dZmass=-1;
  int chosenZ=-1;
  // Choose the Z by mass closest to mZ (~equivalent to ordering by best SM ME but would be equally valid for BSM MEs as well)
  for (unsigned int iV=iSortedVstart; iV<nSortedVs; iV++){
    MELAParticle const* associatedV = associatedVs.at(iV);
    if (!PDGHelpers::isAZBoson(associatedV->id)) continue;
    if (!PDGHelpers::isALepton(associatedV->getDaughter(0)->id)) continue;
    if (chosenZ<0 || fabs(associatedV->m()-PDGHelpers::Zmass)<dZmass){ dZmass=fabs(associatedV->m()-PDGHelpers::Zmass); chosenZ=(int) iV; }
  }
  if (chosenZ>=0){
    // Disable every associated Z boson (and not its daughters!) unless it is the chosen one
    for (unsigned int disableV=iSortedVstart; disableV<nSortedVs; disableV++){
      bool flag=(((int) disableV)==chosenZ);
      MELAParticle* einV = associatedVs.at(disableV);
      if (PDGHelpers::isAZBoson(einV->id)) einV->setSelected(flag);
    }

    for (MELACluster* theCluster:me_clusters){
      if (theCluster->getName()==clustertype){
        // Re-compute all related hypotheses first...
        theCluster->computeAll();
        // ...then force an update the cluster
        theCluster->forceUpdate(); // MELAComputation::testMaximizationCache still prevents update for a second time if there are no contingencies.
      }
    } // End loop over clusters

  } // End if chosenZ>=0
  // Re-enable every associated Z boson and its daughters unless it is the chosen one
  for (MELAParticle* einV:associatedVs){ if (PDGHelpers::isAZBoson(einV->id)) einV->setSelected(true); }
}
// ME computations that require no quark initial state
void HVVTree::updateMELAClusters_NoInitialQ(const string clustertype, bool isGen){
  MELACandidate* melaCand = melaHelpers::melaHandle->getCurrentCandidate();
  if (!melaCand) return;

  // Manipulate the candidate
  // Assign 0 to the id of quark mothers
  std::vector<int> motherIds;
  for (int imot=0; imot<melaCand->getNMothers(); imot++){
    motherIds.push_back(melaCand->getMother(imot)->id);
    if (PDGHelpers::isAQuark(melaCand->getMother(imot)->id)) melaCand->getMother(imot)->id = 0;
  }

  std::vector<MELACluster*>& me_clusters = (isGen ? lheme_clusters : recome_clusters);
  for (MELACluster* theCluster:me_clusters){
    if (theCluster->getName()==clustertype){
      // Re-compute all related hypotheses first...
      theCluster->computeAll();
      // ...then update the cluster
      theCluster->update();
    }
  }

  // Restore the candidate properties
  for (int imot=0; imot<melaCand->getNMothers(); imot++) melaCand->getMother(imot)->id = motherIds.at(imot); // Restore all mother ids
}
// ME computations that require no gluon initial state
void HVVTree::updateMELAClusters_NoInitialG(const string clustertype, bool isGen){
  MELACandidate* melaCand = melaHelpers::melaHandle->getCurrentCandidate();
  if (!melaCand) return;

  // Manipulate the candidate
  // Assign 0 to the id of gluon mothers
  std::vector<int> motherIds;
  for (int imot=0; imot<melaCand->getNMothers(); imot++){
    motherIds.push_back(melaCand->getMother(imot)->id);
    if (PDGHelpers::isAGluon(melaCand->getMother(imot)->id)) melaCand->getMother(imot)->id = 0;
  }

  std::vector<MELACluster*>& me_clusters = (isGen ? lheme_clusters : recome_clusters);
  for (MELACluster* theCluster:me_clusters){
    if (theCluster->getName()==clustertype){
      // Re-compute all related hypotheses first...
      theCluster->computeAll();
      // ...then update the cluster
      theCluster->update();
    }
  }

  // Restore the candidate properties
  for (int imot=0; imot<melaCand->getNMothers(); imot++) melaCand->getMother(imot)->id = motherIds.at(imot); // Restore all mother ids
}
// ME computations that require no gluons as associated particles
void HVVTree::updateMELAClusters_NoAssociatedG(const string clustertype, bool isGen){
  MELACandidate* melaCand = melaHelpers::melaHandle->getCurrentCandidate();
  if (!melaCand) return;

  // Manipulate the candidate
  // Assign 0 to the id of gluon mothers
  std::vector<int> ajetIds;
  for (int ijet=0; ijet<melaCand->getNAssociatedJets(); ijet++){
    ajetIds.push_back(melaCand->getAssociatedJet(ijet)->id);
    if (PDGHelpers::isAGluon(melaCand->getAssociatedJet(ijet)->id)) melaCand->getAssociatedJet(ijet)->id = 0;
  }

  std::vector<MELACluster*>& me_clusters = (isGen ? lheme_clusters : recome_clusters);
  for (MELACluster* theCluster:me_clusters){
    if (theCluster->getName()==clustertype){
      // Re-compute all related hypotheses first...
      theCluster->computeAll();
      // ...then update the cluster
      theCluster->update();
    }
  }

  // Restore the candidate properties
  for (int ijet=0; ijet<melaCand->getNAssociatedJets(); ijet++) melaCand->getAssociatedJet(ijet)->id = ajetIds.at(ijet); // Restore all jets
}
// ME computations that require no gluon initial state and no gluons as associated particles
void HVVTree::updateMELAClusters_NoInitialGNoAssociatedG(const string clustertype, bool isGen){
  MELACandidate* melaCand = melaHelpers::melaHandle->getCurrentCandidate();
  if (!melaCand) return;

  // Manipulate the candidate
  // Assign 0 to the id of gluon mothers
  std::vector<int> motherIds;
  std::vector<int> ajetIds;
  for (int imot=0; imot<melaCand->getNMothers(); imot++){
    motherIds.push_back(melaCand->getMother(imot)->id);
    if (PDGHelpers::isAGluon(melaCand->getMother(imot)->id)) melaCand->getMother(imot)->id = 0;
  }
  for (int ijet=0; ijet<melaCand->getNAssociatedJets(); ijet++){
    ajetIds.push_back(melaCand->getAssociatedJet(ijet)->id);
    if (PDGHelpers::isAGluon(melaCand->getAssociatedJet(ijet)->id)) melaCand->getAssociatedJet(ijet)->id = 0;
  }

  std::vector<MELACluster*>& me_clusters = (isGen ? lheme_clusters : recome_clusters);
  for (MELACluster* theCluster:me_clusters){
    if (theCluster->getName()==clustertype){
      // Re-compute all related hypotheses first...
      theCluster->computeAll();
      // ...then update the cluster
      theCluster->update();
    }
  }

  // Restore the candidate properties
  for (int imot=0; imot<melaCand->getNMothers(); imot++) melaCand->getMother(imot)->id = motherIds.at(imot); // Restore all mother ids
  for (int ijet=0; ijet<melaCand->getNAssociatedJets(); ijet++) melaCand->getAssociatedJet(ijet)->id = ajetIds.at(ijet); // Restore all jets
}
// ME computations that require best Z, W or VBF topology at LO (no gluons)
void HVVTree::updateMELAClusters_BestLOAssociatedZ(const string clustertype, bool isGen){
  MELACandidate* melaCand = melaHelpers::melaHandle->getCurrentCandidate();
  if (!melaCand) return;

  // Manipulate the candidate
  // Assign 0 to the id of gluon mothers
  std::vector<int> motherIds;
  std::vector<int> ajetIds;
  for (int imot=0; imot<melaCand->getNMothers(); imot++){
    motherIds.push_back(melaCand->getMother(imot)->id);
    if (PDGHelpers::isAGluon(melaCand->getMother(imot)->id)) melaCand->getMother(imot)->id = 0;
  }
  for (int ijet=0; ijet<melaCand->getNAssociatedJets(); ijet++){
    ajetIds.push_back(melaCand->getAssociatedJet(ijet)->id);
    if (PDGHelpers::isAGluon(melaCand->getAssociatedJet(ijet)->id)) melaCand->getAssociatedJet(ijet)->id = 0;
  }

  std::vector<MELAParticle*> associatedVs; // Vector of Zs to loop over
  for (MELAParticle* Vtmp:melaCand->getAssociatedSortedVs()){
    if (Vtmp && PDGHelpers::isAZBoson(Vtmp->id) && Vtmp->getNDaughters()>=1){
      bool passSelection=true;
      for (MELAParticle* dauVtmp:Vtmp->getDaughters()) passSelection &= dauVtmp->passSelection;
      if (!passSelection) continue;

      associatedVs.push_back(Vtmp);
    }
  }

  // Give precedence to leptonic V decays
  bool hasALepV=false;
  for (MELAParticle* Vtmp:associatedVs){
    const int& Vtmp_dauid = Vtmp->getDaughter(0)->id;
    if (
      PDGHelpers::isALepton(Vtmp_dauid)
      ||
      PDGHelpers::isANeutrino(Vtmp_dauid)
      ){
      hasALepV=true;
      break;
    }
  }

  std::vector<MELAParticle*> modifiedVs;
  MELAParticle* bestVbyMass=nullptr;
  float bestVMassDiff=-1;
  for (MELAParticle* Vtmp:associatedVs){
    if (
      hasALepV &&
      PDGHelpers::isAJet(Vtmp->getDaughter(0)->id)
      ){
      for (MELAParticle* dauVtmp:Vtmp->getDaughters()) dauVtmp->setSelected(false);
      modifiedVs.push_back(Vtmp);
    }
    else if (!bestVbyMass || fabs(Vtmp->m()-PDGHelpers::Zmass)<bestVMassDiff){
      bestVMassDiff=fabs(Vtmp->m()-PDGHelpers::Zmass);
      bestVbyMass = Vtmp;
    }
  }

  std::vector<MELACluster*>& me_clusters = (isGen ? lheme_clusters : recome_clusters);
  for (MELACluster* theCluster:me_clusters){
    if (theCluster->getName()==clustertype){
      // Re-compute all related hypotheses first...
      theCluster->computeAll();
      // ...then update the cluster
      theCluster->update();
    }
  }

  // Restore the candidate properties
  for (int imot=0; imot<melaCand->getNMothers(); imot++) melaCand->getMother(imot)->id = motherIds.at(imot); // Restore all mother ids
  for (int ijet=0; ijet<melaCand->getNAssociatedJets(); ijet++) melaCand->getAssociatedJet(ijet)->id = ajetIds.at(ijet); // Restore all jets
  for (MELAParticle* Vtmp:modifiedVs){ for (MELAParticle* dauVtmp:Vtmp->getDaughters()) dauVtmp->setSelected(true); }
}
void HVVTree::updateMELAClusters_BestLOAssociatedW(const string clustertype, bool isGen){
  MELACandidate* melaCand = melaHelpers::melaHandle->getCurrentCandidate();
  if (!melaCand) return;

  // Manipulate the candidate
  // Assign 0 to the id of gluon mothers
  std::vector<int> motherIds;
  std::vector<int> ajetIds;
  for (int imot=0; imot<melaCand->getNMothers(); imot++){
    motherIds.push_back(melaCand->getMother(imot)->id);
    if (PDGHelpers::isAGluon(melaCand->getMother(imot)->id)) melaCand->getMother(imot)->id = 0;
  }
  for (int ijet=0; ijet<melaCand->getNAssociatedJets(); ijet++){
    ajetIds.push_back(melaCand->getAssociatedJet(ijet)->id);
    if (PDGHelpers::isAGluon(melaCand->getAssociatedJet(ijet)->id)) melaCand->getAssociatedJet(ijet)->id = 0;
  }

  std::vector<MELAParticle*> associatedVs; // Vector of Ws to loop over
  for (MELAParticle* Vtmp:melaCand->getAssociatedSortedVs()){
    if (Vtmp && PDGHelpers::isAWBoson(Vtmp->id) && Vtmp->getNDaughters()>=1){
      bool passSelection=true;
      for (MELAParticle* dauVtmp:Vtmp->getDaughters()) passSelection &= dauVtmp->passSelection;
      if (!passSelection) continue;

      associatedVs.push_back(Vtmp);
    }
  }
  // Give precedence to leptonic V decays
  bool hasALepV=false;
  for (MELAParticle* Vtmp:associatedVs){
    const int& Vtmp_dauid = Vtmp->getDaughter(0)->id;
    if (
      PDGHelpers::isALepton(Vtmp_dauid)
      ||
      PDGHelpers::isANeutrino(Vtmp_dauid)
      ){
      hasALepV=true;
      break;
    }
  }

  std::vector<MELAParticle*> modifiedVs;
  MELAParticle* bestVbyMass=nullptr;
  float bestVMassDiff=-1;
  for (MELAParticle* Vtmp:associatedVs){
    if (
      hasALepV &&
      PDGHelpers::isAJet(Vtmp->getDaughter(0)->id)
      ){
      for (MELAParticle* dauVtmp:Vtmp->getDaughters()) dauVtmp->setSelected(false);
      modifiedVs.push_back(Vtmp);
    }
    else if (!bestVbyMass || fabs(Vtmp->m()-PDGHelpers::Wmass)<bestVMassDiff){
      bestVMassDiff = fabs(Vtmp->m()-PDGHelpers::Wmass);
      bestVbyMass = Vtmp;
    }
  }

  std::vector<MELACluster*>& me_clusters = (isGen ? lheme_clusters : recome_clusters);
  for (MELACluster* theCluster:me_clusters){
    if (theCluster->getName()==clustertype){
      // Re-compute all related hypotheses first...
      theCluster->computeAll();
      // ...then update the cluster
      theCluster->update();
    }
  }

  // Restore the candidate properties
  for (int imot=0; imot<melaCand->getNMothers(); imot++) melaCand->getMother(imot)->id = motherIds.at(imot); // Restore all mother ids
  for (int ijet=0; ijet<melaCand->getNAssociatedJets(); ijet++) melaCand->getAssociatedJet(ijet)->id = ajetIds.at(ijet); // Restore all jets
  for (MELAParticle* Vtmp:modifiedVs){ for (MELAParticle* dauVtmp:Vtmp->getDaughters()) dauVtmp->setSelected(true); }
}
void HVVTree::updateMELAClusters_BestLOAssociatedVBF(const string clustertype, bool isGen){
  // Same as updateMELAClusters_NoInitialGNoAssociatedG, but keep a separate function for future studies
  MELACandidate* melaCand = melaHelpers::melaHandle->getCurrentCandidate();
  if (!melaCand) return;

  // Manipulate the candidate
  // Assign 0 to the id of gluon mothers
  std::vector<int> motherIds;
  std::vector<int> ajetIds;
  for (int imot=0; imot<melaCand->getNMothers(); imot++){
    motherIds.push_back(melaCand->getMother(imot)->id);
    if (PDGHelpers::isAGluon(melaCand->getMother(imot)->id)) melaCand->getMother(imot)->id = 0;
  }
  for (int ijet=0; ijet<melaCand->getNAssociatedJets(); ijet++){
    ajetIds.push_back(melaCand->getAssociatedJet(ijet)->id);
    if (PDGHelpers::isAGluon(melaCand->getAssociatedJet(ijet)->id)) melaCand->getAssociatedJet(ijet)->id = 0;
  }

  std::vector<MELACluster*>& me_clusters = (isGen ? lheme_clusters : recome_clusters);
  for (MELACluster* theCluster:me_clusters){
    if (theCluster->getName()==clustertype){
      // Re-compute all related hypotheses first...
      theCluster->computeAll();
      // ...then update the cluster
      theCluster->update();
    }
  }

  // Restore the candidate properties
  for (int imot=0; imot<melaCand->getNMothers(); imot++) melaCand->getMother(imot)->id = motherIds.at(imot); // Restore all mother ids
  for (int ijet=0; ijet<melaCand->getNAssociatedJets(); ijet++) melaCand->getAssociatedJet(ijet)->id = ajetIds.at(ijet); // Restore all jets
}
// ME computations that can approximate the NLO QCD (-/+ MiNLO extra jet) phase space to LO QCD in signal VBF or VH
// Use these for POWHEG samples
// MELACandidateRecaster has very specific use cases, so do not use these functions for other cases.
void HVVTree::updateMELAClusters_BestNLOVHApproximation(const string clustertype, bool isGen){
  MELACandidate* melaCand = melaHelpers::melaHandle->getCurrentCandidate();
  if (!melaCand) return;

  // Check if any clusters request this computation
  bool clustersRequest=false;
  std::vector<MELACluster*>& me_clusters = (isGen ? lheme_clusters : recome_clusters);
  for (MELACluster* theCluster:me_clusters){
    if (theCluster->getName()==clustertype){
      clustersRequest=true;
      break;
    }
  }
  if (!clustersRequest) return;

  // Need one recaster for each of ZH and WH, so distinguish by the cluster name
  TVar::Production candScheme;
  if (clustertype.find("BestNLOZHApproximation")!=string::npos) candScheme = TVar::Had_ZH;
  else if (clustertype.find("BestNLOWHApproximation")!=string::npos) candScheme = TVar::Had_WH;
  else return;

  MELACandidateRecaster recaster(candScheme);
  MELACandidate* candModified=nullptr;
  MELAParticle* bestAV = MELACandidateRecaster::getBestAssociatedV(melaCand, candScheme);
  if (bestAV){
    recaster.copyCandidate(melaCand, candModified);
    recaster.deduceLOVHTopology(candModified);
    melaHelpers::melaHandle->setCurrentCandidate(candModified);
  }
  else return; // No associated Vs found. The algorithm won't work.

  for (MELACluster* theCluster:me_clusters){
    if (theCluster->getName()==clustertype){
      // Re-compute all related hypotheses first...
      theCluster->computeAll();
      // ...then update the cluster
      theCluster->update();
    }
  }

  delete candModified;
  melaHelpers::melaHandle->setCurrentCandidate(melaCand); // Go back to the original candidate
}
void HVVTree::updateMELAClusters_BestNLOVBFApproximation(const string clustertype, bool isGen){
  MELACandidate* melaCand = melaHelpers::melaHandle->getCurrentCandidate();
  if (!melaCand) return;

  // Check if any clusters request this computation
  bool clustersRequest=false;
  std::vector<MELACluster*>& me_clusters = (isGen ? lheme_clusters : recome_clusters);
  for (MELACluster* theCluster:me_clusters){
    if (theCluster->getName()==clustertype){
      clustersRequest=true;
      break;
    }
  }
  if (!clustersRequest) return;

  // Need one recaster for VBF
  TVar::Production candScheme;
  if (clustertype.find("BestNLOVBFApproximation")!=string::npos) candScheme = TVar::JJVBF;
  else return;

  MELACandidateRecaster recaster(candScheme);
  MELACandidate* candModified=nullptr;
  recaster.copyCandidate(melaCand, candModified);
  recaster.reduceJJtoQuarks(candModified);
  melaHelpers::melaHandle->setCurrentCandidate(candModified);

  for (MELACluster* theCluster:me_clusters){
    if (theCluster->getName()==clustertype){
      // Re-compute all related hypotheses first...
      theCluster->computeAll();
      // ...then update the cluster
      theCluster->update();
    }
  }

  delete candModified;
  melaHelpers::melaHandle->setCurrentCandidate(melaCand); // Go back to the original candidate
}
