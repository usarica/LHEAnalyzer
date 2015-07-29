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
  else if (branchname.find("m4l")!=string::npos) myProduction = TVar::ZZGG; // Important to quote it here before bkg_VAMCFM
  else if (branchname.find("bkg_s")!=string::npos) myProduction = TVar::ZZQQB_S;
  else if (branchname.find("bkg_tu")!=string::npos) myProduction = TVar::ZZQQB_TU;
  else if (branchname.find("bkg")!=string::npos) myProduction = TVar::ZZQQB;

  if (branchname.find("p0plus_m4l")!=string::npos) myProcess = TVar::HSMHiggs;
  else if (branchname.find("bkg")!=string::npos) myProcess = TVar::bkgZZ;
  else if (branchname.find("ggzz_VAMCFM")!=string::npos) myProcess = TVar::bkgZZ;
  else if (branchname.find("ggzz")!=string::npos && branchname.find("VAMCFM")!=string::npos) myProcess = TVar::bkgZZ_SMHiggs;

  if (branchname.find("ScaleUp")!=string::npos) superMELASyst = TVar::SMSyst_ScaleUp;
  else if (branchname.find("ScaleDown")!=string::npos) superMELASyst = TVar::SMSyst_ScaleDown;
  else if (branchname.find("ResUp")!=string::npos) superMELASyst = TVar::SMSyst_ResUp;
  else if (branchname.find("ResDown")!=string::npos) superMELASyst = TVar::SMSyst_ResDown;
  melaHandle->setProcess(myProcess, myME, myProduction);
  // ...phew!

  // Set MELA flags
  if (branchname.find("Gen")!=string::npos){
    leptonInterfScheme = TVar::InterfOn;
    melaHandle->setMelaHiggsWidth(genPOLEWidth);
//    if (branchname.find("VAMCFM")!=string::npos) mePoleScale = wPOLE/wPOLEStandard; // pGen for MCFM needs to be scaled consistently.
  }
  melaHandle->setMelaLeptonInterference(leptonInterfScheme);

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
  if (myProduction==TVar::WH || myProduction==TVar::ZH){ // Avoid factor of 4
    Higgs_daughters.clear();
    for (int d=0; d<3; d++) { Higgs_daughters.push_back(nullFourVector); Higgs_daughter_ids.push_back(0); }
    Higgs_daughters.push_back(pHiggs); Higgs_daughter_ids.push_back(Higgs_id);
  }

  vector<TLorentzVector> V_daughters;
  vector<int> V_daughter_ids;
  if ((myProduction==TVar::WH || myProduction==TVar::ZH) && branchname.find("leptonic")!=string::npos){ // Leptonic VH
    if(cand->getNAssociatedLeptons()<2) return result;
    else{
      for (int av=2; av<cand->getNSortedVs(); av++){ // Vs are already constructed after the jets, leptons etc. are ordered. Looking at their daughters therefore also orders the jets, leptons etc. by default.
        Particle* pV = cand->getSortedV(av);
        if (pV==0) continue;
        Particle* pVD1 = pV->getDaughter(0);
        Particle* pVD2 = pV->getDaughter(1);
        if (
          ((myProduction==TVar::WH && isAWBoson(pV->id)) && (isALepton(pVD1->id) || isANeutrino(pVD1->id)))
          ||
          ((myProduction==TVar::ZH && isAZBoson(pV->id)) && (isALepton(pVD1->id) || isANeutrino(pVD1->id)))
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



  // Note: No implementation of tt/bbH yet!
  if (myProduction==TVar::JJVBF || myProduction==TVar::JJGG || myProduction==TVar::JH){
    Float_t tmpME = 0, auxME = 1;
    melaHandle->computeProdP(
      V_daughters.at(0), V_daughter_ids.at(0),
      V_daughters.at(1), V_daughter_ids.at(1),
      pHiggs, Higgs_id,
      nullFourVector, 0,
      selfDHggcoupl,
      selfDHvvcoupl,
      selfDHwwcoupl,
      tmpME
      );
    melaHandle->get_PAux(auxME);
    result = tmpME*auxME;
  }
  else if ((myProduction==TVar::WH || myProduction==TVar::ZH) && myME==TVar::JHUGen){
    Float_t tmpME = 0;
    
    // Unfortunately, cannot use vector::data
    TLorentzVector myjets[2] ={ V_daughters.at(0), V_daughters.at(1) };
    TLorentzVector daughters[4] ={ Higgs_daughters.at(0), Higgs_daughters.at(1), Higgs_daughters.at(2), Higgs_daughters.at(3) };
    int daughterids[4] ={ Higgs_daughter_ids.at(0), Higgs_daughter_ids.at(1), Higgs_daughter_ids.at(2), Higgs_daughter_ids.at(3) };
    int jetids[2] ={ V_daughter_ids.at(0), V_daughter_ids.at(1) };

    melaHandle->computeProdP(
      myjets,
      daughters,
      jetids,
      daughterids,
      false,
      selfDHvvcoupl,
      tmpME);
    result = tmpME;
  }
  else if (myProduction==TVar::ZZGG){
    Float_t tmpME = 0;

    // Will need to get angles and then compute...

    result = tmpME;
  }
  return result;
}

