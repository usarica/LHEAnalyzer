#include "HVVTree.h"


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
  if (options->initializeMELABranches() || doSetAddress) bookMELABranches(doSetAddress);
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


void HVVTree::bookMELABranches(const bool& doSetAddress){
  vector<string> tmpBranchList = constructMELABranchList(doSetAddress);
  for (unsigned int b=0; b<tmpBranchList.size(); b++){
    bool isReserved = reserveBranch(tmpBranchList.at(b), BaseTree::bFloat, doSetAddress);
    if (isReserved){
      melaProbBranches.push_back(tmpBranchList.at(b));
    }
  }
}

vector<string> HVVTree::constructMELABranchList(const bool& doSetAddress){
  vector<string> blist;
  return blist;
}
void HVVTree::setupMELASignalMECases(vector<string>& accumulatedlist, TVar::Production prod, TVar::MatrixElement me, bool isGen, bool isProdME, bool doSetAddress){
  return;
}
vector<string> HVVTree::getMELASignalMEBranches(TVar::Production prod, TVar::MatrixElement me, vector<string> gList, vector<int> gCountRe, vector<int> gCountIm, bool isGen, bool isProdME, bool doSetAddress){
  vector<string> blist;
  return blist;
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
  if (pH==0) return;

  MELACandidate* pHactive=pH;
  MELACandidateRecaster* recaster=0;
  if (isGen){
    if (options->doRecastGenTopologyToLOQCDVH()){
      recaster = new MELACandidateRecaster(options->getSampleProductionId().first);
      MELACandidate* candModified=nullptr;
      MELAParticle* bestAV = MELACandidateRecaster::getBestAssociatedV(pHactive, options->getSampleProductionId().first);
      if (bestAV){
        recaster->copyCandidate(pHactive, candModified);
        recaster->deduceLOVHTopology(candModified);
        pHactive = candModified;
      }
    }
    else if (options->doRecastGenTopologyToLOQCDVBF()){
      recaster = new MELACandidateRecaster(TVar::JJVBF);
      MELACandidate* candModified=nullptr;
      recaster->copyCandidate(pHactive, candModified);
      recaster->reduceJJtoQuarks(candModified);
      pHactive = candModified;
    }
  }

  string varname;
  string strcore = "ZZ";
  if (isGen) strcore = "GenH";
  varname = strcore + "Mass"; setVal(varname, (pHactive ? pHactive->m() : 0.));
  varname = strcore + "Pt"; setVal(varname, (pHactive ? pHactive->pt() : 0));
  varname = strcore + "Pz"; setVal(varname, (pHactive ? pHactive->z() : 0));
  varname = strcore + "Phi"; setVal(varname, (pHactive ? pHactive->phi() : 0));

  fillCandidateDaughters(pHactive, isGen);
  fillDaughterProducts(pHactive, isGen);
  fillAssociatedInfo(pHactive, isGen);

  if (melaHelpers::melaHandle && pHactive){
    melaHelpers::melaHandle->setCurrentCandidate(pHactive);

    if (options->doComputeDecayAngles()) fillDecayAngles(isGen);
    if (options->doComputeVBFAngles()) fillVBFProductionAngles(isGen);
    if (options->doComputeVHAngles()) fillVHProductionAngles(isGen);
    if (options->doComputeTTHAngles()) fillTTHProductionAngles(isGen);
    if (melaProbBranches.size()>0) fillMELAProbabilities(isGen); // Do it at the last step

    melaHelpers::melaHandle->resetInputEvent();
  }

  if (pHactive!=pH) delete pHactive;
  if (recaster) delete recaster;
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
      int iLep = 2*v+d+1;

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
    for (int aa=0; aa<pH->getNAssociatedJets(); aa++){
      MELAParticle* apart = pH->getAssociatedJet(aa);
      tmpAssociatedParticle.push_back(apart);
    }
    for (int aa=0; aa<pH->getNAssociatedLeptons(); aa++){
      MELAParticle* apart = pH->getAssociatedLepton(aa);
      if(!PDGHelpers::isANeutrino(apart->id)) tmpAssociatedParticle.push_back(apart);
    }
    for (int aa=0; aa<pH->getNAssociatedNeutrinos(); aa++){
      MELAParticle* apart = pH->getAssociatedNeutrino(aa);
      tmpAssociatedParticle.push_back(apart);
    }
    for (int aa=0; aa<pH->getNAssociatedPhotons(); aa++){
      MELAParticle* apart = pH->getAssociatedPhoton(aa);
      tmpAssociatedParticle.push_back(apart);
    }
  }

  while (tmpAssociatedParticle.size()>0){ // Re-sort all associated particles by leading pT (categories are individually sorted, but mixing categories loses this sorting)
    MELAParticle* tmpPart=0;
    int pos=0;
    for (unsigned int el=0; el<tmpAssociatedParticle.size(); el++){
      if (tmpPart==0){
        tmpPart = tmpAssociatedParticle.at(el); pos=el;
      }
      else if (tmpPart->pt()<tmpAssociatedParticle.at(el)->pt()){
        tmpPart = tmpAssociatedParticle.at(el); pos=el;
      }// Safer to do in two steps
    }
    AssociatedParticle.push_back(tmpPart);
    tmpAssociatedParticle.erase(tmpAssociatedParticle.begin()+pos);
  }

  NAssociatedVs = (pH ? pH->getNSortedVs()-2 : 0);
  for (int av=2; av<NAssociatedVs+2; av++){
    MELAParticle* pAV = pH->getSortedV(av);
    AssociatedV.push_back(pAV);
    MELAParticle* avd1 = pAV->getDaughter(0);
    MELAParticle* avd2 = pAV->getDaughter(1);

    for (unsigned int aa=0; aa<AssociatedParticle.size(); aa++){
      if (avd1==AssociatedParticle.at(aa)) AssociatedV_Particle1Index.push_back(aa);
      else if (avd2==AssociatedParticle.at(aa)) AssociatedV_Particle2Index.push_back(aa);
    }
  }

  NAssociatedTops = (pH ? pH->getNAssociatedTops() : 0);
  for (int it=0; it<NAssociatedTops; it++){
    MELATopCandidate_t* pTop = pH->getAssociatedTop(it);
    AssociatedTop.push_back(pTop);
    MELAParticle* tpp = pTop->getPartnerParticle();
    MELAParticle* tWf = pTop->getWFermion();
    MELAParticle* tWfb = pTop->getWAntifermion();

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

  string varname;
  string strcore;

  if (pH->getNAssociatedJets()>1){
    DijetMass = (pH->getAssociatedJet(0)->p4+pH->getAssociatedJet(1)->p4).M();
    varname = "DijetMass";
    if (isGen) varname.insert(0, "Gen");
    setVal(varname, DijetMass);
    DijetVVMass = (pH->p4+pH->getAssociatedJet(0)->p4+pH->getAssociatedJet(1)->p4).M();
    varname = "DijetVVMass";
    if (isGen) varname.insert(0, "Gen");
    setVal(varname, DijetVVMass);
    dRjet = pH->getAssociatedJet(0)->deltaR(pH->getAssociatedJet(1)->p4);
    varname = "DRjet";
    if (isGen) varname.insert(0, "Gen");
    setVal(varname, dRjet);
  }
  if (pH->getNAssociatedLeptons()>1){
    DileptonMass = (pH->getAssociatedLepton(0)->p4+pH->getAssociatedLepton(1)->p4).M();
    varname = "DileptonMass";
    if (isGen) varname.insert(0, "Gen");
    setVal(varname, DileptonMass);
    DileptonVVMass = (pH->p4+pH->getAssociatedLepton(0)->p4+pH->getAssociatedLepton(1)->p4).M();
    varname = "DileptonVVMass";
    if (isGen) varname.insert(0, "Gen");
    setVal(varname, DileptonVVMass);
    dRlep = pH->getAssociatedLepton(0)->deltaR(pH->getAssociatedLepton(1)->p4);
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
      hs, h1, h2, Phi, Phi1,

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
  for (unsigned int b=0; b<melaProbBranches.size(); b++){
    string branchname = melaProbBranches.at(b);
    if ((isGen && branchname.find("Gen")==string::npos) || (!isGen && branchname.find("Gen")!=string::npos)) continue;
    Float_t prob = melaHelpers::melaBranchMEInterpreter(branchname);
    setVal(branchname, prob);
  }
}


void HVVTree::fillEventVariables(const Float_t& weight, const Int_t& passSelection){
  setVal("MC_weight", weight);
  if (options && options->processRecoInfo()) setVal("isSelected", passSelection);
}

