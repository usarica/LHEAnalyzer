#include "../interface/HVVTree.h"


bool HVVTree::reserveBranch(string branchname, BaseTree::BranchTypes branchtype, bool doSetAddress){
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
void HVVTree::bookAllBranches(bool doSetAddress){
  if (!options){
    cerr << "HVVTree::bookAllBranches -> No options are set for the HVVTree!" << endl;
    return;
  }

  reserveBranch("MC_weight", BaseTree::bFloat, doSetAddress);
  if (options->processGenInfo()){
    reserveBranch("genFinalState", BaseTree::bInt, doSetAddress);

    bookPtEtaPhiMassIdBranches("Mother", BaseTree::bVectorDouble, doSetAddress, true, true, true);
    bookPtEtaPhiMassIdBranches("H", BaseTree::bFloat, doSetAddress, false, true, true);
    bookPtEtaPhiMassIdBranches("Z1", BaseTree::bFloat, doSetAddress, false, false, true);
    bookPtEtaPhiMassIdBranches("Z2", BaseTree::bFloat, doSetAddress, false, false, true);
    bookPtEtaPhiMassIdBranches("Za", BaseTree::bFloat, doSetAddress, false, false, true);
    bookPtEtaPhiMassIdBranches("Zb", BaseTree::bFloat, doSetAddress, false, false, true);

    bookPtEtaPhiMassIdBranches("AssociatedParticle", BaseTree::bVectorDouble, doSetAddress, true, false, true);
    reserveBranch("GenDijetMass", BaseTree::bFloat, doSetAddress);
    reserveBranch("GenDileptonMass", BaseTree::bFloat, doSetAddress);
    reserveBranch("GenDijetVVMass", BaseTree::bFloat, doSetAddress);
    reserveBranch("GenDileptonVVMass", BaseTree::bFloat, doSetAddress);
    reserveBranch("GenDRjet", BaseTree::bFloat, doSetAddress);
    reserveBranch("GenDRlepton", BaseTree::bFloat, doSetAddress);

    reserveBranch("NGenAssociatedVs", BaseTree::bInt, doSetAddress);
    bookPtEtaPhiMassIdBranches("AssociatedV", BaseTree::bVectorDouble, doSetAddress, true, false, true);
    reserveBranch("GenAssociatedV_Particle1Index", BaseTree::bVectorInt, doSetAddress);
    reserveBranch("GenAssociatedV_Particle2Index", BaseTree::bVectorInt, doSetAddress);

    bookPtEtaPhiMassIdBranches("Lep1", BaseTree::bFloat, doSetAddress, true, false, true);
    bookPtEtaPhiMassIdBranches("Lep2", BaseTree::bFloat, doSetAddress, true, false, true);
    bookPtEtaPhiMassIdBranches("Lep3", BaseTree::bFloat, doSetAddress, true, false, true);
    bookPtEtaPhiMassIdBranches("Lep4", BaseTree::bFloat, doSetAddress, true, false, true);
  }
  if (options->processRecoInfo()){
    reserveBranch("isSelected", BaseTree::bInt, doSetAddress);

    bookPtEtaPhiMassIdBranches("ZZ", BaseTree::bFloat, doSetAddress, false, true, false);
    bookPtEtaPhiMassIdBranches("Z1", BaseTree::bFloat, doSetAddress, false, false, false);
    bookPtEtaPhiMassIdBranches("Z2", BaseTree::bFloat, doSetAddress, false, false, false);
    bookPtEtaPhiMassIdBranches("Za", BaseTree::bFloat, doSetAddress, false, false, false);
    bookPtEtaPhiMassIdBranches("Zb", BaseTree::bFloat, doSetAddress, false, false, false);

    bookPtEtaPhiMassIdBranches("AssociatedParticle", BaseTree::bVectorDouble, doSetAddress, true, false, false);
    reserveBranch("DijetMass", BaseTree::bFloat, doSetAddress);
    reserveBranch("DileptonMass", BaseTree::bFloat, doSetAddress);
    reserveBranch("DijetVVMass", BaseTree::bFloat, doSetAddress);
    reserveBranch("DileptonVVMass", BaseTree::bFloat, doSetAddress);
    reserveBranch("DRjet", BaseTree::bFloat, doSetAddress);
    reserveBranch("DRlepton", BaseTree::bFloat, doSetAddress);

    reserveBranch("NAssociatedVs", BaseTree::bInt, doSetAddress);
    bookPtEtaPhiMassIdBranches("AssociatedV", BaseTree::bVectorDouble, doSetAddress, true, false, false);
    reserveBranch("AssociatedV_Particle1Index", BaseTree::bVectorInt, doSetAddress);
    reserveBranch("AssociatedV_Particle2Index", BaseTree::bVectorInt, doSetAddress);

    bookPtEtaPhiMassIdBranches("Lep1", BaseTree::bFloat, doSetAddress, true, false, false);
    bookPtEtaPhiMassIdBranches("Lep2", BaseTree::bFloat, doSetAddress, true, false, false);
    bookPtEtaPhiMassIdBranches("Lep3", BaseTree::bFloat, doSetAddress, true, false, false);
    bookPtEtaPhiMassIdBranches("Lep4", BaseTree::bFloat, doSetAddress, true, false, false);
  }
  bookAngularBranches(doSetAddress);
  if (options->initializeMELA() || doSetAddress) bookMELABranches(doSetAddress);
  actuateBranches(doSetAddress);
}


void HVVTree::bookPtEtaPhiMassIdBranches(string owner, BaseTree::BranchTypes btype, bool doSetAddress, bool addId, bool usePz, bool isGen){
  vector<string> tmpBranchList;
  getPtEtaPhiMIdBranches(tmpBranchList, owner, addId, usePz, isGen);
  for (unsigned int b=0; b<tmpBranchList.size(); b++){
    string branchname = tmpBranchList.at(b);
    if(!addId || branchname.find("Id")==string::npos) reserveBranch(tmpBranchList.at(b), btype, doSetAddress);
    else{
      BaseTree::BranchTypes bInttype = BaseTree::bInt;
      if (btype==BaseTree::bVectorDouble) bInttype = BaseTree::bVectorInt;
      reserveBranch(tmpBranchList.at(b), bInttype, doSetAddress);
    }
  }
}
void HVVTree::getPtEtaPhiMIdBranches(vector<string>& blist, string owner, bool addId, bool usePz, bool isGen){
  string strGen = "Gen";
  vector<string> strtmp;

  strtmp.push_back("Pt");
  if (usePz) strtmp.push_back("Pz");
  else strtmp.push_back("Eta");
  strtmp.push_back("Phi");
  strtmp.push_back("Mass");
  if(addId) strtmp.push_back("Id");

  for (unsigned int b=0; b<strtmp.size(); b++){
    string varname = strtmp.at(b);
    varname.insert(0, owner);
    if (isGen) varname.insert(0, strGen);
    blist.push_back(varname);
  }
}
void HVVTree::bookAngularBranches(bool doSetAddress){
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
  }
  for (unsigned int b=0; b<tmpBranchList.size(); b++){
    reserveBranch(tmpBranchList.at(b), BaseTree::bFloat, doSetAddress);
  }
}
void HVVTree::getAngularBranches(vector<string>& blist, Int_t prodFlag /* 0: Decay, 1: VBF, 2: VH */, bool isGen){
  string strGen = "Gen";
  vector<string> strtmp;
  if (prodFlag==0){
    strtmp.push_back("costhetastar");
    strtmp.push_back("helcosthetaZ1");
    strtmp.push_back("helcosthetaZ2");
    strtmp.push_back("helphi");
    strtmp.push_back("phistarZ1");
  }
  else{
    strtmp.push_back("costhetastar");
    strtmp.push_back("helcosthetaV1");
    strtmp.push_back("helcosthetaV2");
    strtmp.push_back("helphi");
    strtmp.push_back("phistarV1");
    if (prodFlag==1){
      strtmp.push_back("Q_V1");
      strtmp.push_back("Q_V2");
    }
  }
  for (unsigned int b=0; b<strtmp.size(); b++){
    string varname = strtmp.at(b);
    if (isGen) varname.insert(0, strGen);
    if (prodFlag==1 && varname.find("Q_V")==string::npos) varname.append("_VBF");
    else if (prodFlag==2) varname.append("_VH");
    else if (prodFlag==3) varname.append("_VHhadronic");
    else if (prodFlag==4) varname.append("_VHleptonic");
    blist.push_back(varname);
  }
}


void HVVTree::bookMELABranches(bool doSetAddress){
  vector<string> tmpBranchList = constructMELABranchList(doSetAddress);
  for (unsigned int b=0; b<tmpBranchList.size(); b++){
    bool isReserved = reserveBranch(tmpBranchList.at(b), BaseTree::bFloat, doSetAddress);
    if (isReserved){
      melaProbBranches.push_back(tmpBranchList.at(b));
    }
  }
}
vector<string> HVVTree::constructMELABranchList(bool doSetAddress){
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


void HVVTree::fillMotherInfo(Particle* mother){
  if (options!=0 && options->processGenInfo() && mother!=0){
    setVal("GenMotherMass", mother->m());
    setVal("GenMotherPt", mother->pt());
    setVal("GenMotherPz", mother->z());
    setVal("GenMotherPhi", mother->phi());
    setVal("GenMotherId", mother->id);
  }
}


void HVVTree::fillCandidate(ZZCandidate* pH, bool isGen){
  if (!options) return;
  if ((!options->processGenInfo() && isGen) || (!options->processRecoInfo() && !isGen)) return;
  if (pH==0) return;

  string varname;
  string strcore = "ZZ";
  if (isGen) strcore = "GenH";

  varname = strcore + "Mass"; setVal(varname, (pH!=0 ? pH->m() : 0.));
  varname = strcore + "Pt"; setVal(varname, (pH!=0 ? pH->pt() : 0));
  varname = strcore + "Pz"; setVal(varname, (pH!=0 ? pH->z() : 0));
  varname = strcore + "Phi"; setVal(varname, (pH!=0 ? pH->phi() : 0));

  fillCandidateDaughters(pH, isGen);
  fillDaughterProducts(pH, isGen);
  fillAssociatedInfo(pH, isGen);
  if (options->doComputeDecayAngles()) fillDecayAngles(pH, isGen);
  if (options->doComputeVBFAngles()) fillVBFProductionAngles(pH, isGen);
  if (options->doComputeVHAngles()) fillVHProductionAngles(pH, isGen);

  if (melaProbBranches.size()>0) fillMELAProbabilities(pH, isGen); // Do it at the last step
}
void HVVTree::fillCandidateDaughters(ZZCandidate* pH, bool isGen){
  string varname;
  string strcore;

  Particle* pV1=(pH!=0 ? pH->getSortedV(0) : 0);
  Particle* pV2=(pH!=0 ? pH->getSortedV(1) : 0);
  if (pH!=0){
    if (pV1!=0 && pV1->getMother(0)!=pH) pV1=0;
    if (pV2!=0 && pV2->getMother(0)!=pH) pV2=0;
  }

  strcore = "Z1";
  if (isGen) strcore.insert(0, "Gen");
  varname = strcore + "Mass"; setVal(varname, (pV1!=0 ? pV1->m() : 0));
  varname = strcore + "Pt"; setVal(varname, (pV1!=0 ? pV1->pt() : 0));
  varname = strcore + "Eta"; setVal(varname, (pV1!=0 ? (pV1->t()>0 ? pV1->eta() : 0) : 0));
  varname = strcore + "Phi"; setVal(varname, (pV1!=0 ? (pV1->pt()>0 ? pV1->phi() : 0) : 0));
  strcore = "Z2";
  if (isGen) strcore.insert(0, "Gen");
  varname = strcore + "Mass"; setVal(varname, (pV2!=0 ? pV2->m() : 0));
  varname = strcore + "Pt"; setVal(varname, (pV2!=0 ? pV2->pt() : 0));
  varname = strcore + "Eta"; setVal(varname, (pV2!=0 ? (pV2->t()>0 ? pV2->eta() : 0) : 0));
  varname = strcore + "Phi"; setVal(varname, (pV2!=0 ? (pV2->pt()>0 ? pV2->phi() : 0) : 0));

  TLorentzVector nullVector(0, 0, 0, 0);
  TLorentzVector pZ1alt = (pH!=0 ? pH->getAlternativeVMomentum(0) : nullVector);
  TLorentzVector pZ2alt = (pH!=0 ? pH->getAlternativeVMomentum(1) : nullVector);

  strcore = "Za";
  if (isGen) strcore.insert(0, "Gen");
  varname = strcore + "Mass"; setVal(varname, (pH!=0 ? pZ1alt.M() : 0));
  varname = strcore + "Pt"; setVal(varname, (pH!=0 ? pZ1alt.Pt() : 0));
  varname = strcore + "Eta"; setVal(varname, (pH!=0 ? (pZ1alt.T()>0 ? pZ1alt.Eta() : 0) : 0));
  varname = strcore + "Phi"; setVal(varname, (pH!=0 ? (pZ1alt.Pt()>0 ? pZ1alt.Phi() : 0) : 0));
  strcore = "Zb";
  if (isGen) strcore.insert(0, "Gen");
  varname = strcore + "Mass"; setVal(varname, (pH!=0 ? pZ2alt.M() : 0));
  varname = strcore + "Pt"; setVal(varname, (pH!=0 ? pZ2alt.Pt() : 0));
  varname = strcore + "Eta"; setVal(varname, (pH!=0 ? (pZ2alt.T()>0 ? pZ2alt.Eta() : 0) : 0));
  varname = strcore + "Phi"; setVal(varname, (pH!=0 ? (pZ2alt.Pt()>0 ? pZ2alt.Phi() : 0) : 0));
}
void HVVTree::fillDaughterProducts(ZZCandidate* pH, bool isGen){
  string varname;
  string strcore = "Lep";
  if (isGen) strcore.insert(0, "Gen");

  int nDau = std::min((pH!=0 ? pH->getNSortedVs() : 0), 2);
  for (int v=0; v<nDau; v++){
    Particle* intermediateV = pH->getSortedV(v);
    if (intermediateV!=0 && intermediateV->getMother(0)!=pH){ intermediateV = 0; nDau--; }
    if (v==nDau) break;
    int nVDau = (intermediateV!=0 ? intermediateV->getNDaughters() : 0);

    for (int d=0; d<2; d++){
      int iLep = 2*v+d+1;
      char cILep[2];
      sprintf(cILep, "%i", iLep);
      string strILep = string(cILep);

      bool isNew = false;
      Particle* lepton = (intermediateV!=0 ? intermediateV->getDaughter(d) : 0);

      if (lepton==0){
        isNew=true;
        TLorentzVector pDaughter(0, 0, 0, 0);
        int idDaughter = 0;
        if (intermediateV!=0 && nVDau==0 && d==0){
          pDaughter = intermediateV->p4;
          idDaughter = intermediateV->id;
        }
        else if (nDau==0 && v==0 && d==0){
          pDaughter = pH->p4;
          idDaughter = pH->id;
        }
        lepton = new Particle(idDaughter, pDaughter);
      }

      varname = strcore + strILep + "Mass"; setVal(varname, (lepton!=0 ? lepton->m() : 0));
      varname = strcore + strILep + "Pt"; setVal(varname, (lepton!=0 ? lepton->pt() : 0));
      varname = strcore + strILep + "Eta"; setVal(varname, (lepton!=0 ? (lepton->t()>0 ? lepton->eta() : 0) : 0));
      varname = strcore + strILep + "Phi"; setVal(varname, (lepton!=0 ? (lepton->pt()>0 ? lepton->phi() : 0) : 0));
      varname = strcore + strILep + "Id"; setVal(varname, (lepton!=0 ? lepton->id : 0));

      if (isNew){ delete lepton; lepton=0; }
    }
  }

  if (isGen){
    Int_t genFinalState=-1;
    if (nDau>=2 && pH->getSortedV(0)->getNDaughters()>=2 && pH->getSortedV(1)->getNDaughters()>=2){
      if (PDGHelpers::isAZBoson(pH->getSortedV(0)->id) && PDGHelpers::isAZBoson(pH->getSortedV(1)->id)){
        // 4l
        if (std::abs(pH->getSortedV(0)->getDaughter(0)->id)==13 && std::abs(pH->getSortedV(1)->getDaughter(0)->id)==13) genFinalState=0; // 4mu
        else if (std::abs(pH->getSortedV(0)->getDaughter(0)->id)==11 && std::abs(pH->getSortedV(1)->getDaughter(0)->id)==11) genFinalState=1; // 4e
        else if ((std::abs(pH->getSortedV(0)->getDaughter(0)->id)==11 && std::abs(pH->getSortedV(1)->getDaughter(0)->id)==13) || (std::abs(pH->getSortedV(0)->getDaughter(0)->id)==13 && std::abs(pH->getSortedV(1)->getDaughter(0)->id)==11)) genFinalState=2; // 2e2mu
        else if (std::abs(pH->getSortedV(0)->getDaughter(0)->id)==15 && std::abs(pH->getSortedV(1)->getDaughter(0)->id)==15) genFinalState=3; // 4tau
        else if ((std::abs(pH->getSortedV(0)->getDaughter(0)->id)==15 && std::abs(pH->getSortedV(1)->getDaughter(0)->id)==13) || (std::abs(pH->getSortedV(0)->getDaughter(0)->id)==13 && std::abs(pH->getSortedV(1)->getDaughter(0)->id)==15)) genFinalState=4; // 2mu2tau
        else if ((std::abs(pH->getSortedV(0)->getDaughter(0)->id)==11 && std::abs(pH->getSortedV(1)->getDaughter(0)->id)==15) || (std::abs(pH->getSortedV(0)->getDaughter(0)->id)==15 && std::abs(pH->getSortedV(1)->getDaughter(0)->id)==11)) genFinalState=5; // 2e2tau
        // 4nu, 4q
        else if (std::abs(pH->getSortedV(0)->getDaughter(0)->id)==std::abs(pH->getSortedV(1)->getDaughter(0)->id)) genFinalState=0;
        else genFinalState=2;
      }
      else genFinalState=2; // WW has no interference
    }
    setVal("genFinalState", genFinalState);
  }
}
void HVVTree::fillAssociatedInfo(ZZCandidate* pH, bool isGen){
  vector<Particle*> AssociatedParticle;
  vector<Particle*> tmpAssociatedParticle;
  Float_t DijetMass=-1;
  Float_t DileptonMass=-1;
  Float_t DijetVVMass=-1;
  Float_t DileptonVVMass=-1;
  Float_t dRjet=0;
  Float_t dRlep=0;

  Int_t NAssociatedVs=0;
  vector<Particle*> AssociatedV;
  vectorInt AssociatedV_Particle1Index;
  vectorInt AssociatedV_Particle2Index;

  if (pH!=0){
    for (int aa=0; aa<pH->getNAssociatedJets(); aa++){
      Particle* apart = pH->getAssociatedJet(aa);
      tmpAssociatedParticle.push_back(apart);
    }
    for (int aa=0; aa<pH->getNAssociatedLeptons(); aa++){
      Particle* apart = pH->getAssociatedLepton(aa);
      if(!PDGHelpers::isANeutrino(apart->id)) tmpAssociatedParticle.push_back(apart);
    }
    for (int aa=0; aa<pH->getNAssociatedNeutrinos(); aa++){
      Particle* apart = pH->getAssociatedNeutrino(aa);
      tmpAssociatedParticle.push_back(apart);
    }
    for (int aa=0; aa<pH->getNAssociatedPhotons(); aa++){
      Particle* apart = pH->getAssociatedPhoton(aa);
      tmpAssociatedParticle.push_back(apart);
    }
  }

  while (tmpAssociatedParticle.size()>0){ // Re-sort all associated particles by leading pT (categories are individually sorted, but mixing categories loses this sorting)
    Particle* tmpPart=0;
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

  NAssociatedVs = (pH!=0 ? pH->getNSortedVs()-2 : 0);
  for (int av=2; av<NAssociatedVs+2; av++){
    Particle* pAV = pH->getSortedV(av);
    AssociatedV.push_back(pAV);
    Particle* avd1 = pAV->getDaughter(0);
    Particle* avd2 = pAV->getDaughter(1);

    for (unsigned int aa=0; aa<AssociatedParticle.size(); aa++){
      if (avd1==AssociatedParticle.at(aa)) AssociatedV_Particle1Index.push_back(aa);
      else if (avd2==AssociatedParticle.at(aa)) AssociatedV_Particle2Index.push_back(aa);
    }
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
}

void HVVTree::fillDecayAngles(ZZCandidate* pH, bool isGen){
  Float_t helcosthetaZ1=0, helcosthetaZ2=0, helphi=0, costhetastar=0, phistarZ1=0;
  TLorentzVector nullVector(0, 0, 0, 0);

  if (pH!=0 && pH->getNSortedVs()>=2 && pH->getSortedV(0)->getNDaughters()>=1 && pH->getSortedV(1)->getNDaughters()>=1){
    Particle* dau[2][2]={ { 0 } };
    for (int vv=0; vv<2; vv++){
      Particle* Vi = pH->getSortedV(vv);
      for (int dd=0; dd<Vi->getNDaughters(); dd++){
        dau[vv][dd] = Vi->getDaughter(dd);
      }
    }
    TUtil::computeAngles(
      (dau[0][0]!=0 ? dau[0][0]->p4 : nullVector), (dau[0][0]!=0 ? dau[0][0]->id : 0),
      (dau[0][1]!=0 ? dau[0][1]->p4 : nullVector), (dau[0][1]!=0 ? dau[0][1]->id : 0),
      (dau[1][0]!=0 ? dau[1][0]->p4 : nullVector), (dau[1][0]!=0 ? dau[1][0]->id : 0),
      (dau[1][1]!=0 ? dau[1][1]->p4 : nullVector), (dau[1][1]!=0 ? dau[1][1]->id : 0),
      costhetastar, helcosthetaZ1, helcosthetaZ2, helphi, phistarZ1
      );
  }
  // Protect against NaN
  if (!(costhetastar==costhetastar)) costhetastar=0;
  if (!(helcosthetaZ1==helcosthetaZ1)) helcosthetaZ1=0;
  if (!(helcosthetaZ2==helcosthetaZ2)) helcosthetaZ2=0;
  if (!(helphi==helphi)) helphi=0;
  if (!(phistarZ1==phistarZ1)) phistarZ1=0;

  vector<string> varlist;
  getAngularBranches(varlist, 0, isGen);
  for (unsigned int b=0; b<varlist.size(); b++){
    string varname = varlist.at(b);
    if (varname.find("costhetastar")!=string::npos) setVal(varname, costhetastar);
    else if (varname.find("helcosthetaZ1")!=string::npos) setVal(varname, helcosthetaZ1);
    else if (varname.find("helcosthetaZ2")!=string::npos) setVal(varname, helcosthetaZ2);
    else if (varname.find("helphi")!=string::npos) setVal(varname, helphi);
    else if (varname.find("phistarZ1")!=string::npos) setVal(varname, phistarZ1);
    else cerr << "HVVTree::fillDecayAngles -> ERROR: Branch " << varname << " is invalid!" << endl;
  }
}
void HVVTree::fillVBFProductionAngles(ZZCandidate* pH, bool isGen){
  Float_t helcosthetaV1=0, helcosthetaV2=0, helphi=0, costhetastar=0, phistarV1=0;
  Float_t q1sq=0, q2sq=0;

  if (pH!=0){
    TLorentzVector nullFourVector(0, 0, 0, 0);

    Particle* intermediateV[2]={ 0 };
    Particle* daughter[2][2]={ { 0 } };
    for (int v=0; v<2; v++){
      intermediateV[v] = pH->getSortedV(v);
      if (intermediateV[v]!=0 && intermediateV[v]->getMother(0)!=pH) intermediateV[v] = 0;
      if (intermediateV[v]!=0){
        for (int d=0; d<2; d++) daughter[v][d] = intermediateV[v]->getDaughter(d);
      }
    }

    TLorentzVector pDaughter[2][2];
    int idDaughter[2][2];
    for (int v=0; v<2; v++){
      for (int d=0; d<2; d++){
        pDaughter[v][d] = nullFourVector;
        idDaughter[v][d] = 0;
      }
    }
    for (int v=0; v<2; v++){
      for (int d=0; d<2; d++){
        if (daughter[v][d]!=0){
          pDaughter[v][d] = daughter[v][d]->p4;
          idDaughter[v][d] = daughter[v][d]->id;
        }
        else if (intermediateV[v]!=0 && daughter[v][d]==0 && daughter[v][1-d]==0){
          if (d==0){
            pDaughter[v][d] = intermediateV[v]->p4;
            idDaughter[v][d] = intermediateV[v]->id;
          }
        }
        else if (intermediateV[v]==0 && intermediateV[1-v]==0){
          if (v==0){
            pDaughter[v][d] = pH->p4;
            idDaughter[v][d] = pH->id;
          }
        }
      }
    }

    if (pH->getNAssociatedJets()>=1){
      TLorentzVector jet1, jet2;
      int jet1Id, jet2Id;

      TLorentzVector mother1;
      TLorentzVector mother2;
      int mother1Id=0;
      int mother2Id=0;
      TLorentzVector* mother1ref=0;
      TLorentzVector* mother2ref=0;

      jet1 = pH->getAssociatedJet(0)->p4;
      jet1Id = pH->getAssociatedJet(0)->id;
      if (pH->getNAssociatedJets()>=2){
        jet2 = pH->getAssociatedJet(1)->p4;
        jet2Id = pH->getAssociatedJet(1)->id;
      }
      else{
        TLorentzVector others = pH->p4;
        TUtil::computeFakeJet(jet1, others, jet2);
        jet2Id=0;
      }
      if (isGen && pH->getNMothers()>=2){
        mother1 = pH->getMother(0)->p4;
        mother2 = pH->getMother(1)->p4;
        mother1Id = pH->getMother(0)->id;
        mother2Id = pH->getMother(1)->id;
        mother1ref = &mother1;
        mother2ref = &mother2;
      }
      TUtil::computeVBFangles(
        costhetastar, helcosthetaV1, helcosthetaV2, helphi, phistarV1, q1sq, q2sq,
        pDaughter[0][0], idDaughter[0][0],
        pDaughter[0][1], idDaughter[0][1],
        pDaughter[1][0], idDaughter[1][0],
        pDaughter[1][1], idDaughter[1][1],
        jet1, jet1Id,
        jet2, jet2Id,
        mother1ref, mother1Id,
        mother2ref, mother2Id
        );
    }
  }
  // Protect against NaN
  if (!(costhetastar==costhetastar)) costhetastar=0;
  if (!(helcosthetaV1==helcosthetaV1)) helcosthetaV1=0;
  if (!(helcosthetaV2==helcosthetaV2)) helcosthetaV2=0;
  if (!(helphi==helphi)) helphi=0;
  if (!(phistarV1==phistarV1)) phistarV1=0;
  if (!(q1sq==q1sq)) q1sq=0;
  if (!(q2sq==q2sq)) q2sq=0;
  q1sq = (q1sq>=0 ? sqrt(q1sq) : -sqrt(-q1sq));
  q2sq = (q2sq>=0 ? sqrt(q2sq) : -sqrt(-q2sq));

  vector<string> varlist;
  getAngularBranches(varlist, 1, isGen);
  for (unsigned int b=0; b<varlist.size(); b++){
    string varname = varlist.at(b);
    if (varname.find("costhetastar_VBF")!=string::npos) setVal(varname, costhetastar);
    else if (varname.find("helcosthetaV1_VBF")!=string::npos) setVal(varname, helcosthetaV1);
    else if (varname.find("helcosthetaV2_VBF")!=string::npos) setVal(varname, helcosthetaV2);
    else if (varname.find("helphi_VBF")!=string::npos) setVal(varname, helphi);
    else if (varname.find("phistarV1_VBF")!=string::npos) setVal(varname, phistarV1);
    else if (varname.find("Q_V1")!=string::npos) setVal(varname, q1sq);
    else if (varname.find("Q_V2")!=string::npos) setVal(varname, q2sq);
    else cerr << "HVVTree::fillVBFProductionAngles -> ERROR: Branch " << varname << " is invalid!" << endl;
  }
}
void HVVTree::fillVHProductionAngles(ZZCandidate* pH, bool isGen){
  Float_t helcosthetaV1_hadronic=0, helcosthetaV2_hadronic=0, helphi_hadronic=0, costhetastar_hadronic=0, phistarV1_hadronic=0;
  Float_t helcosthetaV1_leptonic=0, helcosthetaV2_leptonic=0, helphi_leptonic=0, costhetastar_leptonic=0, phistarV1_leptonic=0;

  if (pH!=0){
    TLorentzVector nullFourVector(0, 0, 0, 0);

    TLorentzVector mother1;
    TLorentzVector mother2;
    int mother1Id=0;
    int mother2Id=0;
    TLorentzVector* mother1ref=0;
    TLorentzVector* mother2ref=0;
    if (isGen && pH->getNMothers()>=2){
      mother1 = pH->getMother(0)->p4;
      mother2 = pH->getMother(1)->p4;
      mother1Id = pH->getMother(0)->id;
      mother2Id = pH->getMother(1)->id;
      mother1ref = &mother1;
      mother2ref = &mother2;
    }

    Particle* intermediateV[2]={ 0 };
    Particle* daughter[2][2]={ { 0 } };
    for (int v=0; v<2; v++){
      intermediateV[v] = pH->getSortedV(v);
      if (intermediateV[v]!=0 && intermediateV[v]->getMother(0)!=pH) intermediateV[v] = 0;
      if (intermediateV[v]!=0){
        for (int d=0; d<2; d++) daughter[v][d] = intermediateV[v]->getDaughter(d);
      }
    }

    TLorentzVector pDaughter[2][2];
    int idDaughter[2][2];
    for (int v=0; v<2; v++){
      for (int d=0; d<2; d++){
        pDaughter[v][d] = nullFourVector;
        idDaughter[v][d] = 0;
      }
    }
    for (int v=0; v<2; v++){
      for (int d=0; d<2; d++){
        if (daughter[v][d]!=0){
          pDaughter[v][d] = daughter[v][d]->p4;
          idDaughter[v][d] = daughter[v][d]->id;
        }
        else if (intermediateV[v]!=0 && daughter[v][d]==0 && daughter[v][1-d]==0){
          if (d==0){
            pDaughter[v][d] = intermediateV[v]->p4;
            idDaughter[v][d] = intermediateV[v]->id;
          }
        }
        else if (intermediateV[v]==0 && intermediateV[1-v]==0){
          if (v==0){
            pDaughter[v][d] = pH->p4;
            idDaughter[v][d] = pH->id;
          }
        }
      }
    }

    if (pH->getNAssociatedJets()>=1 && options->computeVHAnglesOrder()!=2){
      TLorentzVector jet1, jet2;
      int jet1Id, jet2Id;

      jet1 = pH->getAssociatedJet(0)->p4;
      jet1Id = pH->getAssociatedJet(0)->id;
      if (pH->getNAssociatedJets()>=2){
        jet2 = pH->getAssociatedJet(1)->p4;
        jet2Id = pH->getAssociatedJet(1)->id;
      }
      else{
        TLorentzVector others = pH->p4;
        TUtil::computeFakeJet(jet1, others, jet2);
        jet2Id=0;
      }

      TUtil::computeVHangles(
        costhetastar_hadronic, helcosthetaV1_hadronic, helcosthetaV2_hadronic, helphi_hadronic, phistarV1_hadronic,
        pDaughter[0][0], idDaughter[0][0],
        pDaughter[0][1], idDaughter[0][1],
        pDaughter[1][0], idDaughter[1][0],
        pDaughter[1][1], idDaughter[1][1],
        jet1, jet1Id,
        jet2, jet2Id,
        mother1ref, mother1Id,
        mother2ref, mother2Id
        );
    }
    if (pH->getNAssociatedLeptons()>=1 && options->computeVHAnglesOrder()!=1){
      TLorentzVector jet1, jet2;
      int jet1Id, jet2Id;

      jet1 = pH->getAssociatedLepton(0)->p4;
      jet1Id = pH->getAssociatedLepton(0)->id;
      if (pH->getNAssociatedLeptons()>=2){
        jet2 = pH->getAssociatedLepton(1)->p4;
        jet2Id = pH->getAssociatedLepton(1)->id;
      }
      else{
        TLorentzVector others = pH->p4;
        TUtil::computeFakeJet(jet1, others, jet2);
        jet2Id=0;
      }

      TUtil::computeVHangles(
        costhetastar_leptonic, helcosthetaV1_leptonic, helcosthetaV2_leptonic, helphi_leptonic, phistarV1_leptonic,
        pDaughter[0][0], idDaughter[0][0],
        pDaughter[0][1], idDaughter[0][1],
        pDaughter[1][0], idDaughter[1][0],
        pDaughter[1][1], idDaughter[1][1],
        jet1, jet1Id,
        jet2, jet2Id,
        mother1ref, mother1Id,
        mother2ref, mother2Id
        );
    }
  }
  // Protect against NaN
  if (!(costhetastar_hadronic==costhetastar_hadronic)) costhetastar_hadronic=0;
  if (!(helcosthetaV1_hadronic==helcosthetaV1_hadronic)) helcosthetaV1_hadronic=0;
  if (!(helcosthetaV2_hadronic==helcosthetaV2_hadronic)) helcosthetaV2_hadronic=0;
  if (!(helphi_hadronic==helphi_hadronic)) helphi_hadronic=0;
  if (!(phistarV1_hadronic==phistarV1_hadronic)) phistarV1_hadronic=0;
  if (!(costhetastar_leptonic==costhetastar_leptonic)) costhetastar_leptonic=0;
  if (!(helcosthetaV1_leptonic==helcosthetaV1_leptonic)) helcosthetaV1_leptonic=0;
  if (!(helcosthetaV2_leptonic==helcosthetaV2_leptonic)) helcosthetaV2_leptonic=0;
  if (!(helphi_leptonic==helphi_leptonic)) helphi_leptonic=0;
  if (!(phistarV1_leptonic==phistarV1_leptonic)) phistarV1_leptonic=0;

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
      else if (varname.find("helcosthetaV1_VH")!=string::npos) setVal(varname, helcosthetaV1_hadronic);
      else if (varname.find("helcosthetaV2_VH")!=string::npos) setVal(varname, helcosthetaV2_hadronic);
      else if (varname.find("helphi_VH")!=string::npos) setVal(varname, helphi_hadronic);
      else if (varname.find("phistarV1_VH")!=string::npos) setVal(varname, phistarV1_hadronic);
      else cerr << "HVVTree::fillVHProductionAngles -> ERROR: Branch " << varname << " is invalid!" << endl;
    }
    if (options->computeVHAnglesOrder()==2){
      if (varname.find("costhetastar_VH")!=string::npos) setVal(varname, costhetastar_leptonic);
      else if (varname.find("helcosthetaV1_VH")!=string::npos) setVal(varname, helcosthetaV1_leptonic);
      else if (varname.find("helcosthetaV2_VH")!=string::npos) setVal(varname, helcosthetaV2_leptonic);
      else if (varname.find("helphi_VH")!=string::npos) setVal(varname, helphi_leptonic);
      else if (varname.find("phistarV1_VH")!=string::npos) setVal(varname, phistarV1_leptonic);
      else cerr << "HVVTree::fillVHProductionAngles -> ERROR: Branch " << varname << " is invalid!" << endl;
    }
    if (options->computeVHAnglesOrder()==3){
      if (varname.find("costhetastar_VHhadronic")!=string::npos) setVal(varname, costhetastar_hadronic);
      else if (varname.find("helcosthetaV1_VHhadronic")!=string::npos) setVal(varname, helcosthetaV1_hadronic);
      else if (varname.find("helcosthetaV2_VHhadronic")!=string::npos) setVal(varname, helcosthetaV2_hadronic);
      else if (varname.find("helphi_VHhadronic")!=string::npos) setVal(varname, helphi_hadronic);
      else if (varname.find("phistarV1_VHhadronic")!=string::npos) setVal(varname, phistarV1_hadronic);
      else if (varname.find("costhetastar_VHleptonic")!=string::npos) setVal(varname, costhetastar_leptonic);
      else if (varname.find("helcosthetaV1_VHleptonic")!=string::npos) setVal(varname, helcosthetaV1_leptonic);
      else if (varname.find("helcosthetaV2_VHleptonic")!=string::npos) setVal(varname, helcosthetaV2_leptonic);
      else if (varname.find("helphi_VHleptonic")!=string::npos) setVal(varname, helphi_leptonic);
      else if (varname.find("phistarV1_VHleptonic")!=string::npos) setVal(varname, phistarV1_leptonic);
      else cerr << "HVVTree::fillVHProductionAngles -> ERROR: Branch " << varname << " is invalid!" << endl;
    }
  }
}


void HVVTree::fillMELAProbabilities(ZZCandidate* pH, bool isGen){
  if (pH==0) return;

  for (unsigned int b=0; b<melaProbBranches.size(); b++){
    string branchname = melaProbBranches.at(b);
    if ((isGen && branchname.find("Gen")==string::npos) || (!isGen && branchname.find("Gen")!=string::npos)) continue;
    Float_t prob = melaHelpers::melaBranchMEInterpreter(pH, branchname);
    if (prob!=0) setVal(branchname, prob);
  }
}


void HVVTree::fillEventVariables(Float_t weight, Int_t passSelection){
  setVal("MC_weight", weight);
  if (options!=0 && options->processRecoInfo()) setVal("isSelected", passSelection);
}

