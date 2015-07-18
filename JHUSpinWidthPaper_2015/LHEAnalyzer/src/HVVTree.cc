#include "../interface/HVVTree.h"


void HVVTree::bookAllBranches(){
  bookBranch("MC_weight", BranchTypes::bFloat);
  bookBranch("isSelected", BranchTypes::bInt);
  bookBranch("genFinalState", BranchTypes::bInt);

  bookBranch("GenHMass", BranchTypes::bFloat);
  bookBranch("GenHPt", BranchTypes::bFloat);
  bookBranch("GenHPz", BranchTypes::bFloat);
  bookBranch("GenHPhi", BranchTypes::bFloat);

  bookBranch("GenZ1Mass", BranchTypes::bFloat);
  bookBranch("GenZ1Pt", BranchTypes::bFloat);
  bookBranch("GenZ1Phi", BranchTypes::bFloat);
  bookBranch("GenZ1Eta", BranchTypes::bFloat);

  bookBranch("GenZ2Mass", BranchTypes::bFloat);
  bookBranch("GenZ2Pt", BranchTypes::bFloat);
  bookBranch("GenZ2Phi", BranchTypes::bFloat);
  bookBranch("GenZ2Eta", BranchTypes::bFloat);

  bookBranch("GenZaMass", BranchTypes::bFloat);
  bookBranch("GenZaPt", BranchTypes::bFloat);
  bookBranch("GenZaPhi", BranchTypes::bFloat);
  bookBranch("GenZaEta", BranchTypes::bFloat);

  bookBranch("GenZbMass", BranchTypes::bFloat);
  bookBranch("GenZbPt", BranchTypes::bFloat);
  bookBranch("GenZbPhi", BranchTypes::bFloat);
  bookBranch("GenZbEta", BranchTypes::bFloat);

  bookBranch("GenMotherMass", BranchTypes::bVectorDouble);
  bookBranch("GenMotherPt", BranchTypes::bVectorDouble);
  bookBranch("GenMotherPz", BranchTypes::bVectorDouble);
  bookBranch("GenMotherPhi", BranchTypes::bVectorDouble);
  bookBranch("GenMotherId", BranchTypes::bVectorInt);

  bookBranch("GenNAssociatedVs", BranchTypes::bInt);
  bookBranch("GenAssociatedParticleMass", BranchTypes::bVectorDouble);
  bookBranch("GenAssociatedParticlePt", BranchTypes::bVectorDouble);
  bookBranch("GenAssociatedParticleEta", BranchTypes::bVectorDouble);
  bookBranch("GenAssociatedParticlePhi", BranchTypes::bVectorDouble);
  bookBranch("GenAssociatedParticleId", BranchTypes::bVectorInt);
  bookBranch("GenAssociatedVMass", BranchTypes::bVectorDouble);
  bookBranch("GenAssociatedVPt", BranchTypes::bVectorDouble);
  bookBranch("GenAssociatedVEta", BranchTypes::bVectorDouble);
  bookBranch("GenAssociatedVPhi", BranchTypes::bVectorDouble);
  bookBranch("GenAssociatedVId", BranchTypes::bVectorInt);
  bookBranch("GenAssociatedV_Particle1Index", BranchTypes::bVectorInt);
  bookBranch("GenAssociatedV_Particle2Index", BranchTypes::bVectorInt);

  bookBranch("GenhelcosthetaZ1", BranchTypes::bFloat);
  bookBranch("GenhelcosthetaZ2", BranchTypes::bFloat);
  bookBranch("Genhelphi", BranchTypes::bFloat);
  bookBranch("Gencosthetastar", BranchTypes::bFloat);
  bookBranch("GenphistarZ1", BranchTypes::bFloat);

  bookBranch("GenLep1Mass", BranchTypes::bFloat);
  bookBranch("GenLep2Mass", BranchTypes::bFloat);
  bookBranch("GenLep3Mass", BranchTypes::bFloat);
  bookBranch("GenLep4Mass", BranchTypes::bFloat);
  bookBranch("GenLep1Pt", BranchTypes::bFloat);
  bookBranch("GenLep2Pt", BranchTypes::bFloat);
  bookBranch("GenLep3Pt", BranchTypes::bFloat);
  bookBranch("GenLep4Pt", BranchTypes::bFloat);
  bookBranch("GenLep1Eta", BranchTypes::bFloat);
  bookBranch("GenLep2Eta", BranchTypes::bFloat);
  bookBranch("GenLep3Eta", BranchTypes::bFloat);
  bookBranch("GenLep4Eta", BranchTypes::bFloat);
  bookBranch("GenLep1Phi", BranchTypes::bFloat);
  bookBranch("GenLep2Phi", BranchTypes::bFloat);
  bookBranch("GenLep3Phi", BranchTypes::bFloat);
  bookBranch("GenLep4Phi", BranchTypes::bFloat);
  bookBranch("GenLep1Id", BranchTypes::bInt);
  bookBranch("GenLep2Id", BranchTypes::bInt);
  bookBranch("GenLep3Id", BranchTypes::bInt);
  bookBranch("GenLep4Id", BranchTypes::bInt);


  bookBranch("ZZMass", BranchTypes::bFloat);
  bookBranch("ZZPt", BranchTypes::bFloat);
  bookBranch("ZZPz", BranchTypes::bFloat);
  bookBranch("ZZPhi", BranchTypes::bFloat);

  bookBranch("Z1Mass", BranchTypes::bFloat);
  bookBranch("Z1Pt", BranchTypes::bFloat);
  bookBranch("Z1Phi", BranchTypes::bFloat);
  bookBranch("Z1Eta", BranchTypes::bFloat);

  bookBranch("Z2Mass", BranchTypes::bFloat);
  bookBranch("Z2Pt", BranchTypes::bFloat);
  bookBranch("Z2Phi", BranchTypes::bFloat);
  bookBranch("Z2Eta", BranchTypes::bFloat);

  bookBranch("ZaMass", BranchTypes::bFloat);
  bookBranch("ZaPt", BranchTypes::bFloat);
  bookBranch("ZaPhi", BranchTypes::bFloat);
  bookBranch("ZaEta", BranchTypes::bFloat);

  bookBranch("ZbMass", BranchTypes::bFloat);
  bookBranch("ZbPt", BranchTypes::bFloat);
  bookBranch("ZbPhi", BranchTypes::bFloat);
  bookBranch("ZbEta", BranchTypes::bFloat);

  bookBranch("NAssociatedVs", BranchTypes::bInt);
  bookBranch("AssociatedParticleMass", BranchTypes::bVectorDouble);
  bookBranch("AssociatedParticlePt", BranchTypes::bVectorDouble);
  bookBranch("AssociatedParticleEta", BranchTypes::bVectorDouble);
  bookBranch("AssociatedParticlePhi", BranchTypes::bVectorDouble);
  bookBranch("AssociatedParticleId", BranchTypes::bVectorInt);
  bookBranch("AssociatedVMass", BranchTypes::bVectorDouble);
  bookBranch("AssociatedVPt", BranchTypes::bVectorDouble);
  bookBranch("AssociatedVEta", BranchTypes::bVectorDouble);
  bookBranch("AssociatedVPhi", BranchTypes::bVectorDouble);
  bookBranch("AssociatedVId", BranchTypes::bVectorInt);
  bookBranch("AssociatedV_Particle1Index", BranchTypes::bVectorInt);
  bookBranch("AssociatedV_Particle2Index", BranchTypes::bVectorInt);

  bookBranch("helcosthetaZ1", BranchTypes::bFloat);
  bookBranch("helcosthetaZ2", BranchTypes::bFloat);
  bookBranch("helphi", BranchTypes::bFloat);
  bookBranch("costhetastar", BranchTypes::bFloat);
  bookBranch("phistarZ1", BranchTypes::bFloat);

  bookBranch("Lep1Mass", BranchTypes::bFloat);
  bookBranch("Lep2Mass", BranchTypes::bFloat);
  bookBranch("Lep3Mass", BranchTypes::bFloat);
  bookBranch("Lep4Mass", BranchTypes::bFloat);
  bookBranch("Lep1Pt", BranchTypes::bFloat);
  bookBranch("Lep2Pt", BranchTypes::bFloat);
  bookBranch("Lep3Pt", BranchTypes::bFloat);
  bookBranch("Lep4Pt", BranchTypes::bFloat);
  bookBranch("Lep1Eta", BranchTypes::bFloat);
  bookBranch("Lep2Eta", BranchTypes::bFloat);
  bookBranch("Lep3Eta", BranchTypes::bFloat);
  bookBranch("Lep4Eta", BranchTypes::bFloat);
  bookBranch("Lep1Phi", BranchTypes::bFloat);
  bookBranch("Lep2Phi", BranchTypes::bFloat);
  bookBranch("Lep3Phi", BranchTypes::bFloat);
  bookBranch("Lep4Phi", BranchTypes::bFloat);
  bookBranch("Lep1Id", BranchTypes::bInt);
  bookBranch("Lep2Id", BranchTypes::bInt);
  bookBranch("Lep3Id", BranchTypes::bInt);
  bookBranch("Lep4Id", BranchTypes::bInt);
}

void HVVTree::fillMotherInfo(Particle* mother){
  setVal("GenMotherMass", mother->m());
  setVal("GenMotherPt", mother->pt());
  setVal("GenMotherPz", mother->z());
  setVal("GenMotherPhi", mother->phi());
  setVal("GenMotherId", mother->id);
}


void HVVTree::fillCandidate(ZZCandidate* pH, bool isGen){
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
  fillDecayAngles(pH, isGen);
}
void HVVTree::fillCandidateDaughters(ZZCandidate* pH, bool isGen){
  string varname;
  string strcore;

  Particle* pV1=(pH!=0 ? pH->getSortedV(0) : 0);
  Particle* pV2=(pH!=0 ? pH->getSortedV(1) : 0);

  strcore = "Z1";
  if (isGen) strcore.insert(0, "Gen");
  varname = strcore + "Mass"; setVal(varname, (pV1!=0 ? pV1->m() : 0));
  varname = strcore + "Pt"; setVal(varname, (pV1!=0 ? pV1->pt() : 0));
  varname = strcore + "Eta"; setVal(varname, (pV1!=0 ? pV1->eta() : 0));
  varname = strcore + "Phi"; setVal(varname, (pV1!=0 ? pV1->phi() : 0));
  strcore = "Z2";
  if (isGen) strcore.insert(0, "Gen");
  varname = strcore + "Mass"; setVal(varname, (pV2!=0 ? pV2->m() : 0));
  varname = strcore + "Pt"; setVal(varname, (pV2!=0 ? pV2->pt() : 0));
  varname = strcore + "Eta"; setVal(varname, (pV2!=0 ? pV2->eta() : 0));
  varname = strcore + "Phi"; setVal(varname, (pV2!=0 ? pV2->phi() : 0));

  TLorentzVector nullVector(0, 0, 0, 0);
  TLorentzVector pZ1alt = (pH!=0 ? pH->getAlternativeVMomentum(0) : nullVector);
  TLorentzVector pZ2alt = (pH!=0 ? pH->getAlternativeVMomentum(1) : nullVector);

  strcore = "Za";
  if (isGen) strcore.insert(0, "Gen");
  varname = strcore + "Mass"; setVal(varname, (pH!=0 ? pZ1alt.M() : 0));
  varname = strcore + "Pt"; setVal(varname, (pH!=0 ? pZ1alt.Pt() : 0));
  varname = strcore + "Eta"; setVal(varname, (pH!=0 ? pZ1alt.Eta() : 0));
  varname = strcore + "Phi"; setVal(varname, (pH!=0 ? pZ1alt.Phi() : 0));
  strcore = "Zb";
  if (isGen) strcore.insert(0, "Gen");
  varname = strcore + "Mass"; setVal(varname, (pH!=0 ? pZ2alt.M() : 0));
  varname = strcore + "Pt"; setVal(varname, (pH!=0 ? pZ2alt.Pt() : 0));
  varname = strcore + "Eta"; setVal(varname, (pH!=0 ? pZ2alt.Eta() : 0));
  varname = strcore + "Phi"; setVal(varname, (pH!=0 ? pZ2alt.Phi() : 0));
}
void HVVTree::fillDaughterProducts(ZZCandidate* pH, bool isGen){
  string varname;
  string strcore = "Lep";
  if (isGen) strcore.insert(0, "Gen");

  for (int v=0; v<2; v++){
    for (int d=0; d<2; d++){
      int iLep = 2*v+d+1;
      char cILep[2];
      sprintf(cILep, "%i", iLep);
      string strILep = string(cILep);

      Particle* lepton = (pH!=0 ? pH->getSortedV(v)->getDaughter(d) : 0);

      varname = strcore + strILep + "Mass"; setVal(varname, (lepton!=0 ? lepton->m() : 0));
      varname = strcore + strILep + "Pt"; setVal(varname, (lepton!=0 ? lepton->pt() : 0));
      varname = strcore + strILep + "Eta"; setVal(varname, (lepton!=0 ? lepton->eta() : 0));
      varname = strcore + strILep + "Phi"; setVal(varname, (lepton!=0 ? lepton->phi() : 0));
      varname = strcore + strILep + "Id"; setVal(varname, (lepton!=0 ? lepton->id : 0));
    }
  }

  if (isGen){
    Int_t genFinalState=-1;
    if (pH!=0){
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
  int NAssociatedVs=0;
  vector<Particle*> AssociatedParticle;
  vector<Particle*> tmpAssociatedParticle;

  vector<Particle*> AssociatedV;
  vector<int> AssociatedV_Particle1Index;
  vector<int> AssociatedV_Particle2Index;

  if (pH!=0){
    for (int aa=0; aa<pH->getNAssociatedJets(); aa++){
      Particle* apart = pH->getAssociatedJet(aa);
      tmpAssociatedParticle.push_back(apart);
    }
    for (int aa=0; aa<pH->getNAssociatedLeptons(); aa++){
      Particle* apart = pH->getAssociatedLepton(aa);
      tmpAssociatedParticle.push_back(apart);
    }
    for (int aa=0; aa<pH->getNAssociatedNeutrinos(); aa++){
      Particle* apart = pH->getAssociatedNeutrino(aa);
      tmpAssociatedParticle.push_back(apart);
    }
  }

  while (tmpAssociatedParticle.size()>0){ // Re-sort all associated particles by leading pT (categories are individually soreted, but mixing categories loses this sorting)
    Particle* tmpPart=0;
    int pos=0;
    for (int el=0; el<tmpAssociatedParticle.size(); el++){
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

    for (int aa=0; aa<AssociatedParticle.size(); aa++){
      if (avd1==AssociatedParticle.at(aa)) AssociatedV_Particle1Index.push_back(aa);
      else if (avd2==AssociatedParticle.at(aa)) AssociatedV_Particle2Index.push_back(aa);
    }
  }

  string varname;
  string strcore;

  varname = "NAssociatedVs";
  if (isGen) varname.insert(0, "Gen");
  setVal(varname, NAssociatedVs);

  strcore = "AssociatedParticle";
  if (isGen) strcore.insert(0, "Gen");
  for (int aa=0; aa<AssociatedParticle.size(); aa++){
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
  if (pH!=0) mela::computeAngles(
    pH->getSortedV(0)->getDaughter(0)->p4, pH->getSortedV(0)->getDaughter(0)->id,
    pH->getSortedV(0)->getDaughter(1)->p4, pH->getSortedV(0)->getDaughter(1)->id,
    pH->getSortedV(1)->getDaughter(0)->p4, pH->getSortedV(1)->getDaughter(0)->id,
    pH->getSortedV(1)->getDaughter(1)->p4, pH->getSortedV(1)->getDaughter(1)->id,
    costhetastar, helcosthetaZ1, helcosthetaZ2, helphi, phistarZ1
    );
  // Protect against NaN
  if (!(costhetastar==costhetastar)) costhetastar=0;
  if (!(helcosthetaZ1==helcosthetaZ1)) helcosthetaZ1=0;
  if (!(helcosthetaZ2==helcosthetaZ2)) helcosthetaZ2=0;
  if (!(helphi==helphi)) helphi=0;
  if (!(phistarZ1==phistarZ1)) phistarZ1=0;

  string varname;
  varname = "helcosthetaZ1"; if (isGen) varname.insert(0, "Gen"); setVal(varname, helcosthetaZ1);
  varname = "helcosthetaZ2"; if (isGen) varname.insert(0, "Gen"); setVal(varname, helcosthetaZ2);
  varname = "helphi"; if (isGen) varname.insert(0, "Gen"); setVal(varname, helphi);
  varname = "costhetastar"; if (isGen) varname.insert(0, "Gen"); setVal(varname, costhetastar);
  varname = "phistarZ1"; if (isGen) varname.insert(0, "Gen"); setVal(varname, phistarZ1);
}
//  void HVVTree::fillProductionAngles(Particle* pH, Particle* pV1, Particle* pV2, bool isGen=false);

void HVVTree::fillEventVariables(Float_t weight, Int_t passSelection){
  setVal("MC_weight", weight);
  setVal("isSelected", passSelection);
}

/*
void HVVTree::calculateDecayAngles(TLorentzVector thep4M11, TLorentzVector thep4M12, TLorentzVector thep4M21, TLorentzVector thep4M22, float& costheta1, float& costheta2, float& phi, float& costhetastar, float& phistar1){

  TLorentzVector thep4Z1 = thep4M11+thep4M12;
  TLorentzVector thep4Z2 = thep4M21+thep4M22;
  TLorentzVector thep4H = thep4Z1+thep4Z2;

  float norm;

  TVector3 boostX = -(thep4H.BoostVector());
  TLorentzVector thep4Z1inXFrame(thep4Z1);
  TLorentzVector thep4Z2inXFrame(thep4Z2);
  thep4Z1inXFrame.Boost(boostX);
  thep4Z2inXFrame.Boost(boostX);
  TVector3 theZ1X_p3 = TVector3(thep4Z1inXFrame.X(), thep4Z1inXFrame.Y(), thep4Z1inXFrame.Z());
  TVector3 theZ2X_p3 = TVector3(thep4Z2inXFrame.X(), thep4Z2inXFrame.Y(), thep4Z2inXFrame.Z());


  ///////////////////////////////////////////////
  // check for z1/z2 convention, redefine all 4 vectors with convention
  ///////////////////////////////////////////////	
  TLorentzVector p4H, p4Z1, p4M11, p4M12, p4Z2, p4M21, p4M22;
  p4H = thep4H;

  p4Z1 = thep4Z1; p4M11 = thep4M11; p4M12 = thep4M12;
  p4Z2 = thep4Z2; p4M21 = thep4M21; p4M22 = thep4M22;
  costhetastar = theZ1X_p3.CosTheta();

  // now helicity angles................................
  // ...................................................
  TVector3 boostZ1 = -(p4Z1.BoostVector());
  TLorentzVector p4Z2Z1(p4Z2);
  p4Z2Z1.Boost(boostZ1);
  // find the decay axis
  TVector3 unitx_1(-p4Z2Z1.X(), -p4Z2Z1.Y(), -p4Z2Z1.Z());
  norm = 1/(unitx_1.Mag());
  unitx_1*=norm;
  // boost daughters of z2
  TLorentzVector p4M21Z1(p4M21);
  TLorentzVector p4M22Z1(p4M22);
  p4M21Z1.Boost(boostZ1);
  p4M22Z1.Boost(boostZ1);
  // create z and y axes
  TVector3 p4M21Z1_p3(p4M21Z1.X(), p4M21Z1.Y(), p4M21Z1.Z());
  TVector3 p4M22Z1_p3(p4M22Z1.X(), p4M22Z1.Y(), p4M22Z1.Z());
  TVector3 unitz_1 = p4M21Z1_p3.Cross(p4M22Z1_p3);
  norm = 1/(unitz_1.Mag());
  unitz_1 *= norm;
  TVector3 unity_1 = unitz_1.Cross(unitx_1);

  // calculate theta1
  TLorentzVector p4M11Z1(p4M11);
  p4M11Z1.Boost(boostZ1);
  TVector3 p3M11(p4M11Z1.X(), p4M11Z1.Y(), p4M11Z1.Z());
  TVector3 unitM11 = p3M11.Unit();
  float x_m11 = unitM11.Dot(unitx_1); float y_m11 = unitM11.Dot(unity_1); float z_m11 = unitM11.Dot(unitz_1);
  TVector3 M11_Z1frame(y_m11, z_m11, x_m11);
  costheta1 = M11_Z1frame.CosTheta();

  //////-----------------------old way of calculating phi---------------/////////
  phi = M11_Z1frame.Phi();

  // set axes for other system
  TVector3 boostZ2 = -(p4Z2.BoostVector());
  TLorentzVector p4Z1Z2(p4Z1);
  p4Z1Z2.Boost(boostZ2);
  TVector3 unitx_2(-p4Z1Z2.X(), -p4Z1Z2.Y(), -p4Z1Z2.Z());
  norm = 1/(unitx_2.Mag());
  unitx_2*=norm;
  // boost daughters of z2
  TLorentzVector p4M11Z2(p4M11);
  TLorentzVector p4M12Z2(p4M12);
  p4M11Z2.Boost(boostZ2);
  p4M12Z2.Boost(boostZ2);
  TVector3 p4M11Z2_p3(p4M11Z2.X(), p4M11Z2.Y(), p4M11Z2.Z());
  TVector3 p4M12Z2_p3(p4M12Z2.X(), p4M12Z2.Y(), p4M12Z2.Z());
  TVector3 unitz_2 = p4M11Z2_p3.Cross(p4M12Z2_p3);
  norm = 1/(unitz_2.Mag());
  unitz_2*=norm;
  TVector3 unity_2 = unitz_2.Cross(unitx_2);
  // calcuate theta2
  TLorentzVector p4M21Z2(p4M21);
  p4M21Z2.Boost(boostZ2);
  TVector3 p3M21(p4M21Z2.X(), p4M21Z2.Y(), p4M21Z2.Z());
  TVector3 unitM21 = p3M21.Unit();
  float x_m21 = unitM21.Dot(unitx_2); float y_m21 = unitM21.Dot(unity_2); float z_m21 = unitM21.Dot(unitz_2);
  TVector3 M21_Z2frame(y_m21, z_m21, x_m21);
  costheta2 = M21_Z2frame.CosTheta();

  // calculate phi
  // calculating phi_n
  TLorentzVector n_p4Z1inXFrame(p4Z1);
  TLorentzVector n_p4M11inXFrame(p4M11);
  n_p4Z1inXFrame.Boost(boostX);
  n_p4M11inXFrame.Boost(boostX);
  TVector3 n_p4Z1inXFrame_unit = n_p4Z1inXFrame.Vect().Unit();
  TVector3 n_p4M11inXFrame_unit = n_p4M11inXFrame.Vect().Unit();
  TVector3 n_unitz_1(n_p4Z1inXFrame_unit);
  //// y-axis is defined by neg lepton cross z-axis
  //// the subtle part is here...
  TVector3 n_unity_1 = n_unitz_1.Cross(n_p4M11inXFrame_unit);
  TVector3 n_unitx_1 = n_unity_1.Cross(n_unitz_1);

  TLorentzVector n_p4M21inXFrame(p4M21);
  n_p4M21inXFrame.Boost(boostX);
  TVector3 n_p4M21inXFrame_unit = n_p4M21inXFrame.Vect().Unit();
  //rotate into other plane
  TVector3 n_p4M21inXFrame_unitprime(n_p4M21inXFrame_unit.Dot(n_unitx_1), n_p4M21inXFrame_unit.Dot(n_unity_1), n_p4M21inXFrame_unit.Dot(n_unitz_1));

  ///////-----------------new way of calculating phi-----------------///////
  // float phi_n =  n_p4M21inXFrame_unitprime.Phi();
  /// and then calculate phistar1
  TVector3 n_p4PartoninXFrame_unit(0.0, 0.0, 1.0);
  TVector3 n_p4PartoninXFrame_unitprime(n_p4PartoninXFrame_unit.Dot(n_unitx_1), n_p4PartoninXFrame_unit.Dot(n_unity_1), n_p4PartoninXFrame_unit.Dot(n_unitz_1));
  // negative sign is for arrow convention in paper
  phistar1 = (n_p4PartoninXFrame_unitprime.Phi());
}
*/
