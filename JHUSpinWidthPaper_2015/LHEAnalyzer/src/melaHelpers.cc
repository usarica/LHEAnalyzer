#include "../interface/melaHelpers.h"


namespace melaHelpers{
  Mela* melaHandle=0;
  Float_t genPOLEWidth=4.07e-3;
  Float_t standardPOLEWidth=4.07e-3;
}

using namespace PDGHelpers;

void melaHelpers::setSamplePoleWidth(Float_t width_){ genPOLEWidth = width_; }
void melaHelpers::setStandardPoleWidth(Float_t width_){ standardPOLEWidth = width_; }


Float_t melaHelpers::melaBranchMEInterpreter(const ZZCandidate* cand, string& branchname){
  Float_t result=0;

  // Initialize all variables that could be used
  TVar::MatrixElement myME = TVar::JHUGen;
  TVar::Production myProduction = TVar::ZZGG;
  TVar::Process myProcess = TVar::SelfDefine_spin0;
  TVar::SuperMelaSyst superMELASyst = TVar::SMSyst_None;
  TVar::EventScaleScheme rScaleScheme = TVar::DefaultScaleScheme;
  TVar::LeptonInterference leptonInterfScheme = TVar::DefaultLeptonInterf;
  bool hasSuperMELA=false;

  double selfDHggcoupl[SIZE_HGG][2] ={ { 0 } };
  double selfDHvvcoupl[SIZE_HVV_VBF][2] ={ { 0 } };
  double selfDHwwcoupl[SIZE_HWW_VBF][2] ={ { 0 } };
  double selfDTotalHvvcoupl[SIZE_HVV][2] ={ { 0 } };
  double selfDTotalHwwcoupl[SIZE_HVV][2] ={ { 0 } };
  double ggvvcoupl[2]={ 0, 0 };
  double mePoleScale = 1.;

  // Here we go...
  if (branchname.find("VAJHU")!=string::npos || branchname.find("m4l")!=string::npos) myME = TVar::JHUGen;
  else if (branchname.find("VAMCFM")!=string::npos) { myME = TVar::MCFM; myProcess = TVar::HSMHiggs; }

  if (branchname.find("vbf")!=string::npos) myProduction = TVar::JJVBF;
  else if (branchname.find("hjj")!=string::npos) myProduction = TVar::JJGG;
  else if (branchname.find("hj")!=string::npos){ myProduction = TVar::JH; myProcess = TVar::HSMHiggs; }
  else if (branchname.find("wh")!=string::npos) myProduction = TVar::WH;
  else if (branchname.find("zh")!=string::npos) myProduction = TVar::ZH;
  else if (branchname.find("tth")!=string::npos) myProduction = TVar::ttH;
  else if (branchname.find("bbh")!=string::npos) myProduction = TVar::bbH;
  else if (branchname.find("m4l")!=string::npos) { myProduction = TVar::ZZGG; hasSuperMELA=true; } // Important to quote it here before bkg_VAMCFM
  else if (branchname.find("bkg_s")!=string::npos) myProduction = TVar::ZZQQB_S;
  else if (branchname.find("bkg_tu")!=string::npos) myProduction = TVar::ZZQQB_TU;
  else if (branchname.find("bkg")!=string::npos) myProduction = TVar::ZZQQB;

  if (branchname.find("p0plus_m4l")!=string::npos) myProcess = TVar::HSMHiggs;
  else if (branchname.find("bkg")!=string::npos) myProcess = TVar::bkgZZ;
  else if (branchname.find("ggzz_VAMCFM")!=string::npos || branchname.find("VVzz_VAMCFM")!=string::npos) myProcess = TVar::bkgZZ;
  else if ((branchname.find("ggzz")!=string::npos || branchname.find("VVzz")!=string::npos) && branchname.find("VAMCFM")!=string::npos) myProcess = TVar::bkgZZ_SMHiggs;

  if (branchname.find("ScaleUp")!=string::npos) superMELASyst = TVar::SMSyst_ScaleUp;
  else if (branchname.find("ScaleDown")!=string::npos) superMELASyst = TVar::SMSyst_ScaleDown;
  else if (branchname.find("ResUp")!=string::npos) superMELASyst = TVar::SMSyst_ResUp;
  else if (branchname.find("ResDown")!=string::npos) superMELASyst = TVar::SMSyst_ResDown;
  melaHelpers::melaHandle->setProcess(myProcess, myME, myProduction);
  // ...phew!
  //cout << "melaHelpers::melaBranchMEInterpreter: " << branchname << " received. ";
  //cout << "Production: " << myProduction << ", " << "ME: " << myME << ", " << "Process: " << myProcess << ", " << "superMELASyst: " << superMELASyst << endl;

  // Set MELA flags
  if (branchname.find("Gen")!=string::npos){
    leptonInterfScheme = TVar::InterfOn;
    melaHelpers::melaHandle->setMelaHiggsWidth(melaHelpers::standardPOLEWidth);
    if (branchname.find("VAMCFM")!=string::npos && branchname.find("GHscaled")!=string::npos){
      melaHelpers::melaHandle->setMelaHiggsWidth(melaHelpers::genPOLEWidth);
      mePoleScale = melaHelpers::genPOLEWidth/melaHelpers::standardPOLEWidth; // pGen for MCFM needs to be scaled consistently.
    }
  }
  melaHelpers::melaHandle->setMelaLeptonInterference(leptonInterfScheme);

  TLorentzVector nullFourVector(0,0,0,0);
  TLorentzVector pHiggs = cand->p4;
  int Higgs_id = cand->id; // == 25

  vector<TLorentzVector> Higgs_daughters;
  vector<int> Higgs_daughter_ids;
  for (int v=0; v<2; v++){
    Particle* pV = cand->getSortedV(v);
    if (pV==0) continue;
    for (int d=0; d<2; d++){
      Particle* pVD = pV->getDaughter(d);
      if (pVD==0) continue;
      int iVD = 2*v+d;
      TLorentzVector mom(pVD->x(), pVD->y(), pVD->z(), pVD->t());
      Higgs_daughters.push_back(mom);
      Higgs_daughter_ids.push_back(pVD->id);
    }
  }

  if (Higgs_daughter_ids.size()<4 && (myME==TVar::MCFM || myProduction==TVar::ZZGG)) return result;

  if (hasSuperMELA){ // Exit route for SuperMELA
    TVar::LeptonFlavor superFlavor=TVar::Flavor_Dummy;
    if (abs(Higgs_daughter_ids.at(0))==abs(Higgs_daughter_ids.at(1)) &&
      abs(Higgs_daughter_ids.at(0))==abs(Higgs_daughter_ids.at(2)) &&
      abs(Higgs_daughter_ids.at(0))==abs(Higgs_daughter_ids.at(3))){
      if (abs(Higgs_daughter_ids.at(0))==11) superFlavor=TVar::Flavor_4e;
      else superFlavor=TVar::Flavor_4mu;
    }
    else superFlavor=TVar::Flavor_2e2mu;
    if (superFlavor!=TVar::Flavor_Dummy)
      melaHelpers::melaHandle->computePM4l(
      cand->m(),
      superFlavor,
      superMELASyst,
      result
      );
    return result;
  }

  // Find channel flavor
  int flavor=-1;
  if (abs(Higgs_daughter_ids.at(1))==abs(Higgs_daughter_ids.at(3)) &&
    abs(Higgs_daughter_ids.at(0))==abs(Higgs_daughter_ids.at(2))) flavor = 1;
  else flavor=3;

  // Special treatment for VH involves calling the routine with a single Higgs
  if ((myProduction==TVar::WH || myProduction==TVar::ZH) && myME==TVar::JHUGen){ // Avoid factor of 4
    Higgs_daughters.clear();
    Higgs_daughter_ids.clear();
    for (int d=0; d<3; d++) { Higgs_daughters.push_back(nullFourVector); Higgs_daughter_ids.push_back(0); }
    Higgs_daughters.push_back(pHiggs); Higgs_daughter_ids.push_back(Higgs_id);
  }

  vector<TLorentzVector> V_daughters;
  vector<int> V_daughter_ids;
  // Find the associated production leptons/jets
  if ((myProduction==TVar::WH || myProduction==TVar::ZH)){ // VH
    if (cand->getNSortedVs()<3) return result;
    else{
      for (int av=2; av<cand->getNSortedVs(); av++){ // Vs are already constructed after the jets, leptons etc. are ordered. Looking at their daughters therefore also orders the jets, leptons etc. by default.
        Particle* pV = cand->getSortedV(av);
        if (pV==0) continue;
        Particle* pVD1 = pV->getDaughter(0);
        Particle* pVD2 = pV->getDaughter(1);
        if (
          (
          (
          ((myProduction==TVar::WH && isAWBoson(pV->id)) || (myProduction==TVar::ZH && isAZBoson(pV->id))) && (isALepton(pVD1->id) || isANeutrino(pVD1->id))
          ) && branchname.find("leptonic")!=string::npos // Leptonic VH
          )
          ||
          (
          (
          (myProduction==TVar::WH && ((isAWBoson(pV->id) && isAQuark(pVD1->id)) || (pV->id==0 && pVD1->id==0)))
          ||
          (myProduction==TVar::ZH && ((isAZBoson(pV->id) && isAQuark(pVD1->id)) || (pV->id==0 && pVD1->id==0)))
          ) && branchname.find("hadronic")!=string::npos // Hadronic VH
          )
          ){
          TLorentzVector mom1(pVD1->x(), pVD1->y(), pVD1->z(), pVD1->t());
          V_daughters.push_back(mom1);
          V_daughter_ids.push_back(pVD1->id);
          TLorentzVector mom2(pVD2->x(), pVD2->y(), pVD2->z(), pVD2->t());
          V_daughters.push_back(mom2);
          V_daughter_ids.push_back(pVD2->id);
          break;
        }
      }
      if (V_daughter_ids.size()!=2) return result;
    }
  }
  else{ // Hadronic 2-jet MEs
    if (cand->getNAssociatedJets()>0){ // Uniform treatment for 1-jet cases
      if (cand->getNAssociatedJets()<2 && myProduction!=TVar::JJVBF && myProduction!=TVar::JH) return result; // Escape for non-implemented MEs
      int nrealjets = min(2, cand->getNAssociatedJets());
      for (int av=0; av<nrealjets; av++){ // Jets are already ordered by pT.
        Particle* pJ = cand->getAssociatedJet(av);
        TLorentzVector mom(pJ->x(), pJ->y(), pJ->z(), pJ->t());
        V_daughters.push_back(mom);
        V_daughter_ids.push_back(pJ->id);
      }
      if (V_daughter_ids.size()==1){
        TLorentzVector mom(0, 0, 0, 0);
        V_daughters.push_back(mom);
        V_daughter_ids.push_back(0);
      }
    }
  }

  vector<string> gList[2];
  vector<pair<int, double>> gCoef;
  if (myProcess!=TVar::bkgZZ && !hasSuperMELA){
    if (myProduction==TVar::ttH || myProduction==TVar::bbH || myProduction==TVar::JH || myProduction==TVar::JJGG || (myProduction==TVar::ZZGG && branchname.find("prod")!=string::npos)){
      gList[0].push_back("g2"); gList[1].push_back("0plus");
      gList[0].push_back("g4"); gList[1].push_back("0minus");
      if (myProduction==TVar::ttH){ // ttH
        gCoef.push_back(pair<int, double>(0, 1));
        gCoef.push_back(pair<int, double>(2, 1.593));
      }
      else if (myProduction==TVar::bbH){ // bbH, the same as ttH for now
        gCoef.push_back(pair<int, double>(0, 1));
        gCoef.push_back(pair<int, double>(2, 1.593));
      }
      else if (myProduction==TVar::JH){ // 1-jet ggH, process == HSMHiggs, so selfDefine is not enabled.
        gCoef.push_back(pair<int, double>(0, 0));
        gCoef.push_back(pair<int, double>(2, 0));
      }
      else if (myProduction==TVar::JJGG){ // 2-jets ggH
        double coeffScale = 1.;
        if (myME==TVar::JHUGen) coeffScale = sqrt(1.8e-5);
        gCoef.push_back(pair<int, double>(0, 1.*coeffScale));
        gCoef.push_back(pair<int, double>(2, sqrt(1.0017)*coeffScale));
      }
      else if (myProduction==TVar::ZZGG && branchname.find("prod")!=string::npos){ // 0-jet ggH prod. ME, the same as 2-jets ggH for now in absolute scale
        gCoef.push_back(pair<int, double>(0, 1));
        gCoef.push_back(pair<int, double>(2, sqrt(1.0017)));
      }
    }
    else{
      gList[0].push_back("g1"); gList[1].push_back("0plus"); gCoef.push_back(pair<int, double>(0, 1));
      gList[0].push_back("g2"); gList[1].push_back("0hplus");
      gList[0].push_back("g4"); gList[1].push_back("0minus");
      gList[0].push_back("g1_prime2"); gList[1].push_back("0_g1prime2");
      if (myProduction==TVar::WH){ // WH, to be revised
        double coeffScale = 1.;
        if (myME==TVar::MCFM && branchname.find("Gen")!=string::npos && branchname.find("GHscaled")!=string::npos) coeffScale = pow(mePoleScale, 0.25);
        gCoef.push_back(pair<int, double>(0, 1.*coeffScale));
        gCoef.push_back(pair<int, double>(1, 1.*coeffScale));
        gCoef.push_back(pair<int, double>(3, 0.123659*coeffScale));
        if (myME==TVar::MCFM) gCoef.push_back(pair<int, double>(11, 1.*coeffScale));
        else gCoef.push_back(pair<int, double>(5, 1.*coeffScale));
      }
      else if (myProduction==TVar::ZH){ // ZH, to be revised
        double coeffScale = 1.;
        if (myME==TVar::MCFM && branchname.find("Gen")!=string::npos && branchname.find("GHscaled")!=string::npos) coeffScale = pow(mePoleScale, 0.25);
        gCoef.push_back(pair<int, double>(0, 1.*coeffScale));
        gCoef.push_back(pair<int, double>(1, 1.*coeffScale));
        gCoef.push_back(pair<int, double>(3, 0.143276*coeffScale));
        if (myME==TVar::MCFM) gCoef.push_back(pair<int, double>(11, 1.*coeffScale));
        else gCoef.push_back(pair<int, double>(5, 1.*coeffScale));
      }
      else if (myProduction==TVar::JJVBF){ // VBF
        double coeffScale = 1.;
        if (myME==TVar::MCFM && branchname.find("Gen")!=string::npos && branchname.find("GHscaled")!=string::npos) coeffScale = pow(mePoleScale, 0.25);
        gCoef.push_back(pair<int, double>(0, 1.*coeffScale));
        gCoef.push_back(pair<int, double>(1, 0.270955*coeffScale));
        gCoef.push_back(pair<int, double>(3, 0.29724129*coeffScale));
        if (myME==TVar::MCFM) gCoef.push_back(pair<int, double>(11, 2132.143*coeffScale));
        else gCoef.push_back(pair<int, double>(5, 2132.143*coeffScale));
      }
      else if (myProduction==TVar::ZZGG){ // 0-jet ggH dec. ME
        double coeffScale = 1.;
        if (myME==TVar::MCFM && branchname.find("Gen")!=string::npos && branchname.find("GHscaled")!=string::npos) coeffScale = sqrt(mePoleScale);
        gCoef.push_back(pair<int, double>(0, 1.*coeffScale));
        gCoef.push_back(pair<int, double>(1, 1.65684*coeffScale));
        gCoef.push_back(pair<int, double>(3, 2.55052*coeffScale));
        if (branchname.find("Gen")==string::npos) gCoef.push_back(pair<int, double>(11, -12100.42*coeffScale));
        else gCoef.push_back(pair<int, double>(11, 12100.42*coeffScale)); // For a mostly positive discriminant
      }
    }
  
    int sgList = gList[0].size();
    bool** gFind;
    if (sgList>0){
      gFind = new bool*[sgList];
      for (int gg=0; gg<sgList; gg++){
        gFind[gg] = new bool[2];
        for (int im=0; im<2; im++) gFind[gg][im] = false;
      }
      for (int gg=0; gg<sgList; gg++){
        for (int im=1; im>=0; im--){
          for (int tt=0; tt<2; tt++){
            string chvar = gList[tt].at(gg);
            if (im==1) chvar.append("_pi2");
            if (branchname.find(chvar)!=string::npos) gFind[gg][im]=true; // Does not care if it also finds g1g1_pi2, for example. It is the responsibility of other functions to pass the correct g's.
          }
        }
      }
    }

    if (myProduction==TVar::ttH || myProduction==TVar::bbH || myProduction==TVar::JJGG || (myProduction==TVar::ZZGG && branchname.find("prod")!=string::npos)){
      for (int gg=0; gg<sgList; gg++){
        for (int im=0; im<2; im++){
          if (gFind[gg][im]) selfDHggcoupl[gCoef.at(gg).first][im] = gCoef.at(gg).second;
        }
      }
    }
    else if (myProduction==TVar::JH) {} // Do nothing
    else if (myProduction==TVar::JJVBF || myProduction==TVar::WH || myProduction==TVar::ZH || myProduction==TVar::ZZGG){
      for (int gg=0; gg<sgList; gg++){
        for (int im=0; im<2; im++){
          if (gFind[gg][im]){
            // Account for array size differences between MCFM and JHUGen associated production
            if (myME!=TVar::MCFM && myProduction!=TVar::ZZGG){
              selfDHvvcoupl[gCoef.at(gg).first][im] = gCoef.at(gg).second;
              selfDHwwcoupl[gCoef.at(gg).first][im] = gCoef.at(gg).second;
            }
            else{
              selfDTotalHvvcoupl[gCoef.at(gg).first][im] = gCoef.at(gg).second;
              selfDTotalHwwcoupl[gCoef.at(gg).first][im] = gCoef.at(gg).second;
            }
          }
        }
      }
    }
    // The arrays's task is complete.
    if (sgList>0){
      for (int gg=0; gg<sgList; gg++) delete[] gFind[gg];
      delete[] gFind;
    }
  }

  // Note: No implementation of tt/bbH yet!
  if ((myProduction==TVar::JJVBF || myProduction==TVar::JJGG || myProduction==TVar::JH) && myME==TVar::JHUGen){
    Float_t tmpME = 0, auxME = 1;
    melaHelpers::melaHandle->computeProdP(
      V_daughters.at(0), V_daughter_ids.at(0),
      V_daughters.at(1), V_daughter_ids.at(1),
      pHiggs, Higgs_id,
      nullFourVector, 0,
      selfDHggcoupl,
      selfDHvvcoupl,
      selfDHwwcoupl,
      tmpME
      );
    melaHelpers::melaHandle->get_PAux(auxME);
    result = tmpME*auxME;
  }
  else if ((myProduction==TVar::WH || myProduction==TVar::ZH) && myME==TVar::JHUGen){
    Float_t tmpME = 0;
    
    // Unfortunately, cannot use vector::data
    TLorentzVector myjets[2] ={ V_daughters.at(0), V_daughters.at(1) };
    TLorentzVector daughters[4] ={ Higgs_daughters.at(0), Higgs_daughters.at(1), Higgs_daughters.at(2), Higgs_daughters.at(3) };
    int daughterids[4] ={ Higgs_daughter_ids.at(0), Higgs_daughter_ids.at(1), Higgs_daughter_ids.at(2), Higgs_daughter_ids.at(3) };
    int jetids[2] ={ V_daughter_ids.at(0), V_daughter_ids.at(1) };

    melaHelpers::melaHandle->computeProdP(
      myjets,
      daughters,
      jetids,
      daughterids,
      false,
      selfDHvvcoupl,
      tmpME);
    result = tmpME;
  }
  else if ((myProduction==TVar::ZZGG || myProduction==TVar::ZZQQB) && !hasSuperMELA){
    Float_t tmpME = 0;

    Float_t helcosthetaZ1=0, helcosthetaZ2=0, helphi=0, costhetastar=0, phistarZ1=0;
    if (cand!=0) mela::computeAngles(
      cand->getSortedV(0)->getDaughter(0)->p4, cand->getSortedV(0)->getDaughter(0)->id,
      cand->getSortedV(0)->getDaughter(1)->p4, cand->getSortedV(0)->getDaughter(1)->id,
      cand->getSortedV(1)->getDaughter(0)->p4, cand->getSortedV(1)->getDaughter(0)->id,
      cand->getSortedV(1)->getDaughter(1)->p4, cand->getSortedV(1)->getDaughter(1)->id,
      costhetastar, helcosthetaZ1, helcosthetaZ2, helphi, phistarZ1
      );
    // Protect against NaN
    if (!(costhetastar==costhetastar)) costhetastar=0;
    if (!(helcosthetaZ1==helcosthetaZ1)) helcosthetaZ1=0;
    if (!(helcosthetaZ2==helcosthetaZ2)) helcosthetaZ2=0;
    if (!(helphi==helphi)) helphi=0;
    if (!(phistarZ1==phistarZ1)) phistarZ1=0;

    if (myProduction==TVar::ZZGG && myProcess!=TVar::bkgZZ){
      melaHelpers::melaHandle->computeP(
        (float)cand->m(), (float)cand->getSortedV(0)->m(), (float)cand->getSortedV(1)->m(),
        costhetastar, helcosthetaZ1, helcosthetaZ2, helphi, phistarZ1,
        flavor,
        selfDTotalHvvcoupl,
        tmpME
        );
    }
    else{ // qqZZ or ggZZ bkg.
      melaHelpers::melaHandle->computeP(
        (float)cand->m(), (float)cand->getSortedV(0)->m(), (float)cand->getSortedV(1)->m(),
        costhetastar, helcosthetaZ1, helcosthetaZ2, helphi, phistarZ1,
        flavor,
        tmpME,
        branchname.find("wconst")!=string::npos
        );
    }
    result = tmpME;
  }
  return result;
}


void melaHelpers::computeVHAngles(TLorentzVector thep4Z1, TLorentzVector thep4Z2, TLorentzVector thep4H, double& costheta1, double& costheta2, double& Phi, double& costhetastar, double& Phi1){

  TLorentzVector p4H, p4Z1, p4M11, p4M12, p4Z2, p4M21, p4M22;

  p4M11 = thep4Z1;
  p4M12 = thep4Z2;
  p4Z2 = thep4H;
  p4M21 = p4Z2;
  p4M22 = p4Z2;
  p4Z1 = p4M11+p4M12;
  p4H = p4Z1+p4Z2;

  //// costhetastar
  TVector3 boostX = -(thep4H.BoostVector());
  TLorentzVector thep4Z1inXFrame(p4Z1);
  TLorentzVector thep4Z2inXFrame(p4Z2);
  thep4Z1inXFrame.Boost(boostX);
  thep4Z2inXFrame.Boost(boostX);
  TVector3 theZ1X_p3 = TVector3(thep4Z1inXFrame.X(), thep4Z1inXFrame.Y(), thep4Z1inXFrame.Z());
  TVector3 theZ2X_p3 = TVector3(thep4Z2inXFrame.X(), thep4Z2inXFrame.Y(), thep4Z2inXFrame.Z());
  costhetastar = theZ1X_p3.CosTheta();

  //// --------------------------- costheta1
  TVector3 boostV1 = -(thep4Z1.BoostVector());
  TLorentzVector p4M11_BV1(p4M11);
  TLorentzVector p4M12_BV1(p4M12);
  TLorentzVector p4M21_BV1(p4M21);
  TLorentzVector p4M22_BV1(p4M22);
  p4M11_BV1.Boost(boostV1);
  p4M12_BV1.Boost(boostV1);
  p4M21_BV1.Boost(boostV1);
  p4M22_BV1.Boost(boostV1);

  TLorentzVector p4V2_BV1 = p4M21_BV1 + p4M22_BV1;
  //// costheta1
  costheta1 = -p4V2_BV1.Vect().Dot(p4M11_BV1.Vect())/p4V2_BV1.Vect().Mag()/p4M11_BV1.Vect().Mag();

  //// --------------------------- costheta2
  TVector3 boostV2 = -(thep4Z2.BoostVector());
  TLorentzVector p4M11_BV2(p4M11);
  TLorentzVector p4M12_BV2(p4M12);
  TLorentzVector p4M21_BV2(p4M21);
  TLorentzVector p4M22_BV2(p4M22);
  p4M11_BV2.Boost(boostV2);
  p4M12_BV2.Boost(boostV2);
  p4M21_BV2.Boost(boostV2);
  p4M22_BV2.Boost(boostV2);

  TLorentzVector p4V1_BV2 = p4M11_BV2 + p4M12_BV2;
  //// costheta2
  costheta2 = -p4V1_BV2.Vect().Dot(p4M21_BV2.Vect())/p4V1_BV2.Vect().Mag()/p4M21_BV2.Vect().Mag();

  //// --------------------------- Phi and Phi1
  //    TVector3 boostX = -(thep4H.BoostVector());
  TLorentzVector p4M11_BX(p4M11);
  TLorentzVector p4M12_BX(p4M12);
  TLorentzVector p4M21_BX(p4M21);
  TLorentzVector p4M22_BX(p4M22);

  p4M11_BX.Boost(boostX);
  p4M12_BX.Boost(boostX);
  p4M21_BX.Boost(boostX);
  p4M22_BX.Boost(boostX);

  TVector3 tmp1 = p4M11_BX.Vect().Cross(p4M12_BX.Vect());
  TVector3 tmp2 = p4M21_BX.Vect().Cross(p4M22_BX.Vect());

  TVector3 normal1_BX(tmp1.X()/tmp1.Mag(), tmp1.Y()/tmp1.Mag(), tmp1.Z()/tmp1.Mag());
  TVector3 normal2_BX(tmp2.X()/tmp2.Mag(), tmp2.Y()/tmp2.Mag(), tmp2.Z()/tmp2.Mag());

  //// Phi
  TLorentzVector p4Z1_BX = p4M11_BX + p4M12_BX;
  double tmpSgnPhi = p4Z1_BX.Vect().Dot(normal1_BX.Cross(normal2_BX));
  double sgnPhi = tmpSgnPhi/fabs(tmpSgnPhi);
  Phi = sgnPhi * acos(-1.*normal1_BX.Dot(normal2_BX));


  //////////////

  TVector3 beamAxis(0, 0, 1);
  TVector3 tmp3 = (p4M11_BX + p4M12_BX).Vect();

  TVector3 p3V1_BX(tmp3.X()/tmp3.Mag(), tmp3.Y()/tmp3.Mag(), tmp3.Z()/tmp3.Mag());
  TVector3 tmp4 = beamAxis.Cross(p3V1_BX);
  TVector3 normalSC_BX(tmp4.X()/tmp4.Mag(), tmp4.Y()/tmp4.Mag(), tmp4.Z()/tmp4.Mag());

  //// Phi1
  double tmpSgnPhi1 = p4Z1_BX.Vect().Dot(normal1_BX.Cross(normalSC_BX));
  double sgnPhi1 = tmpSgnPhi1/fabs(tmpSgnPhi1);
  Phi1 = sgnPhi1 * acos(normal1_BX.Dot(normalSC_BX));
}

