#include <iostream>
#include <iomanip>
#include <cmath>
#include <string>
#include <vector>
#include <fstream>
#include <cstdlib>
#include <utility>
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TString.h"
#include "TMath.h"
#include "TF1.h"
#include "TMatrixD.h"
#include "TVectorD.h"
#include "TSpline.h"
#include "TSystem.h"
#include "TGraph.h"
#include "TGraphAsymmErrors.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "./data/QuantFuncMathCore.h"
#include "TLorentzVector.h"
#include "TPaveText.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TColor.h"
#include "TDirectory.h"
#include "./data/ZZ4l_125_6_Samples.h"

using namespace std;
using namespace ROOT::Math;


const int nSIPCuts = 6; // <4 old SIP vs new SIP schemes
TString cutNames[nSIPCuts] ={
  "PVSIPcut",
  "NoSIPcut",
  "Z1SIPcut",
  "4lchi2cut",
  "Z1SIP_4lchi2cut",
  "CombinedCut"
};
const int nCR=2;
double Z1masswidth[2]={ 10, 5 };
double Z1pole = 91.1876;

TString cutLabel[nSIPCuts] ={
  "Old cut",
  "No cut",
  "Z1-SIP",
  "chi**2",
  "New cut",
  "New + old"
};
TString strZ1Category_label[2]={
  "Z1->mumu",
  "Z1->ee"
};
TString strZ2Category_label[6]={
  " + mumu (OS)",
  " + emu (OS)",
  " + ee (OS)",
  " + ee (SS)",
  " + emu (SS)",
  " + mumu (SS)"
};
TString strLep3Category_label[2]={
  " + mu",
  " + e"
};
TString strLepIsoCut_label[3] ={
  "Tight",
  "2p2f",
  "3p1f"
};
TString Z1masswidth_label[3]={
  "40<m1<120",
  "|m1-mZ|<10",
  "|m(3l)-mZ|<5"
};


void restructureSIPYields_CR(int erg_tev, int iCR=0, double m4l_low = 105.6, double m4l_high = 140.6);

TString produceCRname(int iCR, int iwgt, int isocut, int iM1cut, int catZ1, int catZ2, int icut, char* cappend);

void plot1DNewOld(int iCR, int iwgt, int isocut, int iM1cut, int catZ1, int catZ2, char* cappend, int erg_tev, TFile* frecord, double m4l_low = 105.6, double m4l_high = 140.6);
void plot2DNewOld(int iCR, int iwgt, int isocut, int iM1cut, int catZ1, int catZ2, char* cappend, int erg_tev, TFile* frecord, double m4l_low = 105.6, double m4l_high = 140.6);
void plot2DNewOldDoubleRatio(int iCR, int iwgt, int isocut, int iM1cut, int catZ1, int catZ2, char* cappend, int erg_tev, TFile* frecord, double m4l_low = 105.6, double m4l_high = 140.6);
double applyCRselection(
  int option[7],

  float m4l_low, float m4l_high,

  int CRflag,
  float Lep1combRelIsoPF, float Lep2combRelIsoPF, float Lep3combRelIsoPF, float Lep4combRelIsoPF,
  bool Lep1isID, bool Lep2isID, bool Lep3isID, bool Lep4isID,
  short Z1ids,
  short Z2ids,
  float PFMET,

  float Lep_Z1SIP[4],
  float KalmanCandVtx_chi2,
  float LepSIP[4],
  int LepID[4],

  float Z1Mass, float ZZMass,

  float ZXfake_weight,
  float ZXfake_weight_SS,
  float ZXfake_weight_OS[6]
  ){

  int iCR = option[0];
  int icut = option[1];
  int catZ1 = option[2];
  int catZ2 = option[3];
  int iwgt = option[4];
  int isocut = option[5];
  int iM1cut = option[6];

  int Lep1ID = LepID[0];
  int Lep2ID = LepID[1];
  int Lep3ID = LepID[2];
  int Lep4ID = LepID[3];

  float Lep1_Z1SIP = Lep_Z1SIP[0];
  float Lep2_Z1SIP = Lep_Z1SIP[1];
  float Lep3_Z1SIP = Lep_Z1SIP[2];
  float Lep4_Z1SIP = Lep_Z1SIP[3];

  float Lep1SIP = LepSIP[0];
  float Lep2SIP = LepSIP[1];
  float Lep3SIP = LepSIP[2];
  float Lep4SIP = LepSIP[3];

  float strWgt =1;
  bool strnewSIP[nSIPCuts] ={
    TMath::Max(Lep1SIP, TMath::Max(Lep2SIP, TMath::Max(Lep3SIP, Lep4SIP)))<4,
    true,
    fabs(Lep1_Z1SIP)<4 && fabs(Lep2_Z1SIP)<4 && fabs(Lep3_Z1SIP)<5 && fabs(Lep4_Z1SIP)<5,
    KalmanCandVtx_chi2<30,
    fabs(Lep1_Z1SIP)<4 && fabs(Lep2_Z1SIP)<4 && fabs(Lep3_Z1SIP)<5 && fabs(Lep4_Z1SIP)<5 && KalmanCandVtx_chi2<30,
    fabs(Lep1_Z1SIP)<4 && fabs(Lep2_Z1SIP)<4 && fabs(Lep3_Z1SIP)<5 && fabs(Lep4_Z1SIP)<5 && KalmanCandVtx_chi2<30 && TMath::Max(Lep1SIP, TMath::Max(Lep2SIP, TMath::Max(Lep3SIP, Lep4SIP)))<4
  };
  bool strnewSIP_CRZL[nSIPCuts] ={
    TMath::Max(Lep1SIP, TMath::Max(Lep2SIP, Lep3SIP))<4,
    true,
    fabs(Lep1_Z1SIP)<4 && fabs(Lep2_Z1SIP)<4 && fabs(Lep3_Z1SIP)<5,
    KalmanCandVtx_chi2<30,
    fabs(Lep1_Z1SIP)<4 && fabs(Lep2_Z1SIP)<4 && fabs(Lep3_Z1SIP)<5 && KalmanCandVtx_chi2<30,
    fabs(Lep1_Z1SIP)<4 && fabs(Lep2_Z1SIP)<4 && fabs(Lep3_Z1SIP)<5 && KalmanCandVtx_chi2<30 && TMath::Max(Lep1SIP, TMath::Max(Lep2SIP, Lep3SIP))<4
  };
  bool strZ1Category[2]={
    Z1ids == -169,
    Z1ids == -121
  };
  bool strZ2Category[6]={
    (CRflag==6 || CRflag==10), // DO require these, otherwise you will need mass cut, iso cut, ID as well, and has to be the best Z
    Z2ids == -143,
    (CRflag==8 || CRflag==12),
    (CRflag==7 || CRflag==11),
    Z2ids == 143,
    (CRflag==5 || CRflag==9)
  };
  bool strLep3Category[2]={
    fabs(Lep3ID) == 13,
    fabs(Lep3ID) == 11
  };
  bool strLepIsoCut[nCR]={
    TMath::Max(Lep1combRelIsoPF, TMath::Max(Lep2combRelIsoPF, TMath::Max(Lep3combRelIsoPF, Lep4combRelIsoPF)))<0.4 && Lep3isID && Lep4isID,
    Lep3combRelIsoPF<0.4 && Lep3isID  // Lep1,2 already have the cut; what is seen on ntuples is the relIso before FSR
  };
  bool strLepIsoCut_pf[2]={
    (Lep3combRelIsoPF>=0.4 || !Lep3isID)  && (Lep4combRelIsoPF>=0.4 || !Lep4isID), // 2p2f
    (!(Lep3combRelIsoPF<0.4 && Lep4combRelIsoPF<0.4 && Lep3isID && Lep4isID) && !((Lep3combRelIsoPF>=0.4 || !Lep3isID)  && (Lep4combRelIsoPF>=0.4 || !Lep4isID)))  // 3p1f
  };

  bool strcut = strnewSIP[icut];
  if (iCR==1) strcut = strnewSIP_CRZL[icut];

  if (iCR!=1) strcut = strcut && ZZMass>=m4l_low && ZZMass<m4l_high;
  strcut = strcut && strZ1Category[catZ1];

  if (iCR!=1) strcut = strcut && strZ2Category[catZ2];
  else strcut = strcut && strLep3Category[catZ2];

  if (isocut==0 && iCR==1){
    strcut = strcut;
  }
  else if (isocut==1){
    strcut = strcut && strLepIsoCut[iCR];
  }
  else if (isocut>1 && iCR!=1){
    strcut = strcut && strLepIsoCut_pf[isocut-2];
  }
  if (iM1cut==1){
    strcut = strcut && fabs(Z1Mass-Z1pole)<Z1masswidth[iM1cut-1];
  }
  else if (iM1cut==2){
    strcut = strcut && fabs(ZZMass-Z1pole)<Z1masswidth[iM1cut-1];
  }
  if (iCR==1){
    strcut = strcut && PFMET<25;
  }

  if (iwgt==0) strWgt=1;
  else{
    if (catZ2>=3) strWgt = ZXfake_weight_SS;
    else strWgt = ZXfake_weight_OS[icut];
  }

  double result = (strcut ? strWgt : 0);
  return result;
}


TString getstr_hOS(int iZ, int iL, int iEB, int iSel){
  TString hname = "ZL_Z";
  if (iZ==0) hname.Append("mumu_");
  else hname.Append("ee_");
  if (iL==0) hname.Append("mu_");
  else hname.Append("e_");
  hname.Append("Unweighted_");
  if (iSel==0) hname.Append(cutNames[iSel]);
  else hname.Append(cutNames[4]);
  hname.Append("_Tight_ZPeak_pTeta_leading_");
  if (iEB==0) hname.Append("Barrel_ratio");
  else hname.Append("Endcap_ratio");
  return hname;
}


double applyOSSFweight(
  TFile* fOS,
  int CRflag,
  int icut,
  short Z1ids, short Z2ids,
  float LepcombRelIsoPF,
  bool LepisID,

  int LepID,
  float LepPt,
  float LepEta,

  float& OSweighterror
  ){
  float result = 0;
  float OSSFweight=0;

  bool strZ1Category[2]={
    Z1ids == -169, // Zmumu
    Z1ids == -121 // Zee
  };
  bool strZ2Category[6]={
    (CRflag==6 || CRflag==10), // mumu (OS)
    Z2ids == -143, // emu (OS)
    (CRflag==8 || CRflag==12), // ee (OS)
    (CRflag==7 || CRflag==11), // ee (SS)
    Z2ids == 143, // emu (SS)
    (CRflag==5 || CRflag==9) // mumu (SS)
  };

  int pass = 1;
  int histoindex_1=-1;
  int histoindex_2=-1;
  int histoindex_3=-1;
  if (!(strZ1Category[0] || strZ1Category[1])) pass *= 0;
  if (!(strZ2Category[0] || strZ2Category[1] || strZ2Category[2])) pass *= 0;
  if (pass==0) return result;
  if (LepcombRelIsoPF>=0.4 || !LepisID){
    if (strZ1Category[0]) histoindex_1=0;
    if (strZ1Category[1]) histoindex_1=1;
    if (abs(LepID)==13 && LepEta<1.2){
      histoindex_2=0;
      histoindex_3=0;
    }
    if (abs(LepID)==13 && LepEta>=1.2){
      histoindex_2=0;
      histoindex_3=1;
    }
    if (abs(LepID)==11 && LepEta<1.45){
      histoindex_2=1;
      histoindex_3=0;
    }
    if (abs(LepID)==11 && LepEta>=1.45){
      histoindex_2=1;
      histoindex_3=1;
    }
  }
  if (histoindex_1>=0 && histoindex_2>=0 && histoindex_3>=0){
    TH1F* hOS = (TH1F*)fOS->Get(getstr_hOS(histoindex_1, histoindex_2, histoindex_3, icut));
    OSSFweight = hOS->GetBinContent(hOS->GetXaxis()->FindBin(LepPt));
    OSweighterror = hOS->GetBinError(hOS->GetXaxis()->FindBin(LepPt));
    delete hOS;
  }
  result = OSSFweight / (1.-OSSFweight);

  OSweighterror /= pow(1.-OSSFweight, 2);
  return result;
}


float getZXfake_weight_OS(
  int iCR,
  int CRflag, int ZXcut, TFile* ZXWeightTables_OS, TFile* fOS,

  short Z1id,
  int LepId[2],
  float LepPt[2],
  float LepEta[2],
  float LepcombRelIsoPF[2],
  bool LepisID[2],

  float ZZMass,
  float Z1Mass,
  float Z2Mass

  )
{
  float OSSFweight[2] ={ 0 };
  float OSSFweighterror[2] ={ 0 };
  float totalOS = 1;
  float totalOSerror=0;

  int nZ2Leps = 2;
  if (iCR==1) nZ2Leps--;


  short Z1ID = Z1id;
  short Z2ID = LepId[0];
  if (nZ2Leps==2) Z2ID *= LepId[1];

  for (int ilep=0; ilep<nZ2Leps; ilep++) {
    OSSFweight[ilep] = applyOSSFweight(
      fOS,
      CRflag,
      ZXcut,
      Z1ID, Z2ID,
      LepcombRelIsoPF[ilep],
      LepisID[ilep],

      LepId[ilep],
      LepPt[ilep],
      LepEta[ilep],

      OSSFweighterror[ilep]
      );
  } // for loop on Z2 legs

  if (OSSFweight[0]>0){
    totalOS *= OSSFweight[0];
    totalOSerror += pow(OSSFweighterror[0]/OSSFweight[0], 2);
  }
  if (OSSFweight[1]>0){
    totalOS *= OSSFweight[1];
    totalOSerror += pow(OSSFweighterror[1]/OSSFweight[1], 2);
  }

  if (OSSFweight[0]==0 && OSSFweight[1]==0) totalOS=0;
  else if ((OSSFweight[0]==0 || OSSFweight[1]==0) && nZ2Leps==2){ // 3p1f OS: SPECIAL CASE FOR SHAPE REWEIGHTING!!!
    float masses[3] ={
      ZZMass,
      Z1Mass,
      Z2Mass
    };

    bool pass=false;

    TString hname = "OSfakeweight__Z";
    if (Z1ID==-169){
      hname.Append("mumu_"); pass=true;
    }
    else if (Z1ID==-121){
      hname.Append("ee_"); pass=true;
    }

    if (pass){
      if (CRflag==6 || CRflag==10){
        hname.Append("mumuOS_");
      }
      else if (Z2ID == -143){
        hname.Append("emuOS_");
      }
      else if (CRflag==8 || CRflag==12){
        hname.Append("eeOS_");
      }
      else pass = false;
    }
    if (ZXcut<6) hname += cutNames[ZXcut];
    else pass=false;
    hname.Append("_3p1f");
    if (pass){
      TH3F* hCorr=0;
      if (ZXWeightTables_OS!=0 && !ZXWeightTables_OS->IsZombie()){
        hCorr = (TH3F*)ZXWeightTables_OS->Get(hname);
        if (hCorr==0) cout << "No such histogram as " << hname << endl;
        int binx = hCorr->GetXaxis()->FindBin(masses[0]);
        int biny = hCorr->GetYaxis()->FindBin(masses[1]);
        int binz = hCorr->GetZaxis()->FindBin(masses[2]);
        Float_t corrVal = hCorr->GetBinContent(binx, biny, binz);
        Float_t corrErr = hCorr->GetBinError(binx, biny, binz);
        totalOS*=corrVal;
        totalOSerror += pow(corrErr/corrVal, 2);
        delete hCorr;
      }
    }
  }
  totalOSerror = sqrt(totalOSerror) * totalOS;
  if (fabs(totalOS)>0.025){
    totalOS = pow(0.025, 2) / totalOS;
    totalOSerror = pow(0.025 / totalOS, 2)*totalOSerror;
  }

  return totalOS;
}


void restructureSIPYields_OwnFakeRates(int erg_tev=8, double m4l_low = 105.6, double m4l_high = 140.6){
  for (int iCR=0; iCR<1; iCR++) restructureSIPYields_CR(erg_tev, iCR, m4l_low, m4l_high);
}

void restructureSIPYields_CR(int erg_tev, int iCR, double m4l_low, double m4l_high){
  if (iCR>2) return;
  char TREE_NAME[]="SelectedTree";

  TString erg_dir;
  erg_dir.Form("LHC_%iTeV", erg_tev);
  TString comstring;
  comstring.Form("%iTeV", erg_tev);
  //	int EnergyIndex = 0; if (erg_tev == 8) EnergyIndex = 1;
  int sampleindex = kAllSamples-1 + iCR;

  TChain* tc = new TChain(TREE_NAME);

  TString cinput_common_noSIP = user_dir_hep + "No_SIP/" + erg_dir + "/CR/";
  TString cinput = cinput_common_noSIP;
  cinput = cinput + sample_FullSim[sampleindex] + ".root";
  tc->Add(cinput);
  if (tc->GetEntries() == 0){
    cout << "Could not find any files, aborting..." << endl; return;
  }

  string coutput_common = user_dir_hep + "Analysis/Auxiliary/";

  TString fOSfakeinput = "CRdistributions_";
  fOSfakeinput.Append("ZL_");
  fOSfakeinput += comstring;
  fOSfakeinput.Append("_m4l_105.6_140.6.root");
  fOSfakeinput.Prepend(coutput_common);
  TFile* fweight = new TFile(fOSfakeinput, "read");

  TString cfakeaftercorr = "OSOwnfakeweights_AfterCorrection_";
  cfakeaftercorr += comstring;
  cfakeaftercorr.Prepend(coutput_common);
  cfakeaftercorr.Append(".root");
  TFile* ZXWeightTables_OS = new TFile(cfakeaftercorr, "read");

  TString SSinputname = "OSoverSS_";
  SSinputname.Append("MCCR_");
  SSinputname += comstring;
  SSinputname.Append(".root");
  SSinputname.Prepend(coutput_common.c_str());
  TFile* fSSinput = new TFile(SSinputname, "read");
  TH2F* hRatio[nSIPCuts];
  for (int sipcut=0; sipcut<nSIPCuts; sipcut++){
    TString chratio = "hCR_MC_OSSSRatios_";
    chratio.Append(cutLabel[sipcut]);
    hRatio[sipcut] = (TH2F*) fSSinput->Get(chratio);
  }

  TFile* frecord;
  TString fname = "CRdistributions_OwnFakeRates_";
  if (iCR==0) fname.Append("ZLL_");
  else if (iCR==1) fname.Append("ZL_");
  else if (iCR==2) fname.Append("QQBZZ_");
  fname += comstring;
  fname.Append(Form("_m4l_%.1f_%.1f%s", m4l_low, m4l_high, ".root"));
  fname.Prepend(coutput_common.c_str());
  frecord = new TFile(fname, "recreate");

  TFile* foutput[3];
  TH1F* hYield[3][2]={ { 0 } };
  TH1F* hCount[3][2]={ { 0 } };
  string strtout = coutput_common;
  strtout = strtout + "logCR_OwnFakeRates_";
  char ctout[1000];
  sprintf(ctout, "%s%.1f_%.1f_%iTeV%s", strtout.c_str(), m4l_low, m4l_high, erg_tev, ".log");
  ofstream tout(ctout, ios::out);

  if (iCR==0){
    for (int ff=0; ff<3; ff++){
      TString OUTPUT_NAME = "LifetimeKD_RelativeSIPYields_OwnFakeRates_";
      OUTPUT_NAME.Append(Form("%s", "CR"));
      OUTPUT_NAME.Append(Form("_%s_", user_folder[ff]));
      OUTPUT_NAME.Append(comstring);
      OUTPUT_NAME.Append(Form("_m4l%.1f_%.1f", m4l_low, m4l_high));
      OUTPUT_NAME.Append(".root");
      TString coutput = coutput_common + OUTPUT_NAME;
      foutput[ff] = new TFile(coutput, "recreate");
      for (int sos=0; sos<2; sos++){
        TString chyield = "hCR";
        if (sos==0) chyield += "_OS";
        else chyield += "_SS";
        hYield[ff][sos] = new TH1F(chyield, "Yield of any SIP", nSIPCuts - 1, 0, nSIPCuts - 1);
        TString chcount = "hCR_Unweighted";
        if (sos==0) chcount += "_OS";
        else chcount += "_SS";
        hCount[ff][sos] = new TH1F(chcount, "Count of any SIP", nSIPCuts - 1, 0, nSIPCuts - 1);
      }
    }
  }

  frecord->cd();

  TString hname = "htemp";

  float ZXfake_weight;
  float ZXfake_weight_SS;
  float ZXfake_weight_OS[6];
  int CRflag;
  float Lep1combRelIsoPF, Lep2combRelIsoPF, Lep3combRelIsoPF, Lep4combRelIsoPF;
  bool Lep1isID, Lep2isID, Lep3isID, Lep4isID;
  short Z1ids;
  short Z2ids;
  float PFMET;

  float Lep1_Z1SIP, Lep2_Z1SIP, Lep3_Z1SIP, Lep4_Z1SIP;
  float Lep1SIP, Lep2SIP, Lep3SIP, Lep4SIP;
  int Lep1ID, Lep2ID, Lep3ID, Lep4ID;

  float Z1Mass, Z2Mass, ZZMass, ZZPt, ZZEta, ZZPhi;
  float KalmanCandVtx_x, KalmanCandVtx_y, KalmanCandVtx_z, KalmanCandVtx_chi2;
  float OfflinePrimaryVtx_x, OfflinePrimaryVtx_y, OfflinePrimaryVtx_z;

  float Lep1Pt, Lep2Pt, Lep3Pt, Lep4Pt;
  float Lep1Eta, Lep2Eta, Lep3Eta, Lep4Eta;

  tc->SetBranchAddress("ZXfake_weight", &ZXfake_weight);
  tc->SetBranchAddress("ZXfake_weight_SS", &ZXfake_weight_SS);
  tc->SetBranchAddress("CRflag", &CRflag);

  tc->SetBranchAddress("Z1ids", &Z1ids);
  tc->SetBranchAddress("Z2ids", &Z2ids);
  tc->SetBranchAddress("Lep1ID", &Lep1ID);
  tc->SetBranchAddress("Lep1isID", &Lep1isID);
  tc->SetBranchAddress("Lep1combRelIsoPF", &Lep1combRelIsoPF);
  tc->SetBranchAddress("Lep2ID", &Lep2ID);
  tc->SetBranchAddress("Lep2isID", &Lep2isID);
  tc->SetBranchAddress("Lep2combRelIsoPF", &Lep2combRelIsoPF);
  tc->SetBranchAddress("Lep3ID", &Lep3ID);
  tc->SetBranchAddress("Lep3isID", &Lep3isID);
  tc->SetBranchAddress("Lep3combRelIsoPF", &Lep3combRelIsoPF);

  tc->SetBranchAddress("PFMET", &PFMET);
  tc->SetBranchAddress("Z1Mass", &Z1Mass);
  tc->SetBranchAddress("Z2Mass", &Z2Mass);
  tc->SetBranchAddress("ZZMass", &ZZMass);
  tc->SetBranchAddress("ZZPt", &ZZPt);
  tc->SetBranchAddress("ZZEta", &ZZEta);
  tc->SetBranchAddress("ZZPhi", &ZZPhi);
  tc->SetBranchAddress("OfflinePrimaryVtx_x", &OfflinePrimaryVtx_x);
  tc->SetBranchAddress("OfflinePrimaryVtx_y", &OfflinePrimaryVtx_y);
  tc->SetBranchAddress("OfflinePrimaryVtx_z", &OfflinePrimaryVtx_z);
  tc->SetBranchAddress("KalmanCandVtx_x", &KalmanCandVtx_x);
  tc->SetBranchAddress("KalmanCandVtx_y", &KalmanCandVtx_y);
  tc->SetBranchAddress("KalmanCandVtx_z", &KalmanCandVtx_z);

  tc->SetBranchAddress("Lep1Pt", &Lep1Pt);
  tc->SetBranchAddress("Lep2Pt", &Lep2Pt);
  tc->SetBranchAddress("Lep3Pt", &Lep3Pt);
  tc->SetBranchAddress("Lep1Eta", &Lep1Eta);
  tc->SetBranchAddress("Lep2Eta", &Lep2Eta);
  tc->SetBranchAddress("Lep3Eta", &Lep3Eta);


  tc->SetBranchAddress("KalmanCandVtx_chi2", &KalmanCandVtx_chi2);
  tc->SetBranchAddress("Lep1_Z1SIP", &Lep1_Z1SIP);
  tc->SetBranchAddress("Lep2_Z1SIP", &Lep2_Z1SIP);
  tc->SetBranchAddress("Lep3_Z1SIP", &Lep3_Z1SIP);
  tc->SetBranchAddress("Lep1SIP", &Lep1SIP);
  tc->SetBranchAddress("Lep2SIP", &Lep2SIP);
  tc->SetBranchAddress("Lep3SIP", &Lep3SIP);

  if (iCR!=1){
    tc->SetBranchAddress("Lep4ID", &Lep4ID);
    tc->SetBranchAddress("Lep4isID", &Lep4isID);
    tc->SetBranchAddress("Lep4combRelIsoPF", &Lep4combRelIsoPF);
    tc->SetBranchAddress("Lep4_Z1SIP", &Lep4_Z1SIP);
    tc->SetBranchAddress("Lep4SIP", &Lep4SIP);
    tc->SetBranchAddress("Lep4Pt", &Lep4Pt);
    tc->SetBranchAddress("Lep4Eta", &Lep4Eta);
  }


  for (int iwgt=0; iwgt<2; iwgt++){
    if (iCR==1 && iwgt==1) continue;
    for (int isocut=0; isocut<4; isocut++){
      if (iCR==1 && isocut>=2) continue;
      for (int iM1cut=0; iM1cut<3; iM1cut++){
        if (iM1cut==2 && iCR!=1) continue;
        if (iwgt==0) tout << "Un-weighted";
        else tout << "Weighted";
        if (isocut>=1) tout << ", " << strLepIsoCut_label[isocut-1];
        tout << ", " << Z1masswidth_label[iM1cut] << endl;
        tout << "\t"; for (int icut=0; icut<nSIPCuts; icut++) tout << cutLabel[icut] << '\t' << cutLabel[icut] << " Error\t";
        tout << endl;
        for (int catZ1=0; catZ1<2; catZ1++){
          for (int catZ2=0; catZ2<(iCR==0 ? 6 : 2); catZ2++){
            tout << strZ1Category_label[catZ1] << " ";
            if (iCR!=1) tout << strZ2Category_label[catZ2];
            else tout << strLep3Category_label[catZ2];
            for (int icut=0; icut<nSIPCuts; icut++){
              gStyle->SetTitleFont(62, "t");
              gROOT->SetStyle(gStyle->GetName());
              gROOT->ForceStyle();

              TString chsip = hname;
              TH1F* hsip = new TH1F(chsip, "", 1, 0, 1000);
              hsip->Sumw2();

              double ptbins[9]={ 0, 5, 7, 10, 20, 30, 40, 50, 80 };
              double etabins[3]={ 0, 1.2, 2.5 };
              if ((iCR==1 && catZ2==1) || (iCR!=1 && (catZ2==2 || catZ2==3))) etabins[1]=1.45;

              TH1F* hTxy = new TH1F(produceCRname(iCR, iwgt, isocut, iM1cut, catZ1, catZ2, icut, "Txy"), "", 60, -1400, 1400);
              hTxy->Sumw2();
              hTxy->SetXTitle("T_{xy} (#mum)");
              hTxy->SetYTitle("Events / 46.7 #mum");
              TH1F* hmZZ = new TH1F(produceCRname(iCR, iwgt, isocut, iM1cut, catZ1, catZ2, icut, "mZZ"), "", (iCR!=1 ? 5 : 10), m4l_low, m4l_high);
              hmZZ->Sumw2();
              if (iCR!=1){
                hmZZ->SetXTitle("m_{4l} (GeV)");
                hmZZ->SetYTitle("Events / 7 GeV");
              }
              else{
                hmZZ->SetXTitle("m_{3l} (GeV)");
                hmZZ->SetYTitle("Events / 3.5 GeV");
              }
              TH2F* hpteta = new TH2F(produceCRname(iCR, iwgt, isocut, iM1cut, catZ1, catZ2, icut, "pTeta_leading"), "", 8, ptbins, 2, etabins);
              hpteta->Sumw2();
              if (iCR!=1) hpteta->SetXTitle("Leading p^{3,4}_{T} (GeV)");
              else hpteta->SetXTitle("p^{3}_{T} (GeV)");
              hpteta->SetYTitle("|#eta|");
              hpteta->SetZTitle("Events / bin");

              TH1F* hmZ1;
              TH1F* hmZ2;
              TH2F* hptetasub;
              if (iM1cut==0){
                hmZ1 = new TH1F(produceCRname(iCR, iwgt, isocut, iM1cut, catZ1, catZ2, icut, "mZ1"), "", (iCR!=1 ? 10 : 20), 40, 120);
                hmZ1->Sumw2();
                hmZ1->SetXTitle("m_{1} (GeV)");
                if (iCR==0) hmZ1->SetYTitle("Events / 8 GeV");
                else if (iCR==1) hmZ1->SetYTitle("Events / 4 GeV");
              }
              if (iCR!=1){
                hmZ2 = new TH1F(produceCRname(iCR, iwgt, isocut, iM1cut, catZ1, catZ2, icut, "mZ2"), "", 17, 0, 136);
                hmZ2->Sumw2();
                hmZ2->SetXTitle("m_{2} (GeV)");
                hmZ2->SetYTitle("Events / 8 GeV");

                hptetasub = new TH2F(produceCRname(iCR, iwgt, isocut, iM1cut, catZ1, catZ2, icut, "pTeta_subleading"), "", 8, ptbins, 2, etabins);
                hptetasub->Sumw2();
                hptetasub->SetXTitle("Sub-leading p^{3,4}_{T} (GeV)");
                hptetasub->SetYTitle("|#eta|");
                hptetasub->SetZTitle("Events / bin");
              }

              for (int ev=0; ev<tc->GetEntries(); ev++){
                tc->GetEntry(ev);

                int option[7]={
                  iCR, icut, catZ1, catZ2, iwgt, isocut, iM1cut
                };
                int LepID[4]={
                  Lep1ID, Lep2ID, Lep3ID, Lep4ID
                };
                float Lep_Z1SIP[4]={
                  Lep1_Z1SIP, Lep2_Z1SIP, Lep3_Z1SIP, Lep4_Z1SIP
                };
                float LepSIP[4]={
                  Lep1SIP, Lep2SIP, Lep3SIP, Lep4SIP
                };

                int LepID_Z2[2]={
                  Lep3ID, Lep4ID
                };
                float LepPt_Z2[2]={
                  Lep3Pt, Lep4Pt
                };
                float LepEta_Z2[2]={
                  Lep3Eta, Lep4Eta
                };
                float LepcombRelIsoPF_Z2[2]={
                  Lep3combRelIsoPF, Lep4combRelIsoPF
                };
                bool LepisID_Z2[2]={
                  Lep3isID, Lep4isID
                };

                for (int zzx=0; zzx<6; zzx++){
                  ZXfake_weight_OS[zzx] = getZXfake_weight_OS(
                    iCR,
                    CRflag, zzx, ZXWeightTables_OS, fweight,

                    Z1ids,
                    LepID_Z2,
                    LepPt_Z2,
                    LepEta_Z2,
                    LepcombRelIsoPF_Z2,
                    LepisID_Z2,

                    ZZMass,
                    Z1Mass,
                    Z2Mass

                    );
                }
                float wgt = applyCRselection(
                  option,

                  m4l_low, m4l_high,

                  CRflag,
                  Lep1combRelIsoPF, Lep2combRelIsoPF, Lep3combRelIsoPF, Lep4combRelIsoPF,
                  Lep1isID, Lep2isID, Lep3isID, Lep4isID,
                  Z1ids,
                  Z2ids,
                  PFMET,

                  Lep_Z1SIP,
                  KalmanCandVtx_chi2,
                  LepSIP,
                  LepID,

                  Z1Mass, ZZMass,

                  ZXfake_weight,
                  ZXfake_weight_SS,
                  ZXfake_weight_OS

                  );

                if (wgt==0) continue;

                int OSoverSS_biny = -1;
                if (catZ1==0 && catZ2>=3) OSoverSS_biny = 6-catZ2;
                else if (catZ1==1 && catZ2>=3) OSoverSS_biny = 6-catZ2+3;
                float SSwgt_OSoverSS = 1;
                if (OSoverSS_biny>0) SSwgt_OSoverSS = hRatio[icut]->GetBinContent(hRatio[icut]->GetXaxis()->FindBin(ZZMass), OSoverSS_biny);
                if(iwgt>0) wgt *= SSwgt_OSoverSS;
/*
                cout << "wgt: " << wgt << '\t'
                  << "icut: " << option[1] << '\t'
                  << "CRflag: " << CRflag << '\t'
                  << "SIP: " << max(Lep1SIP, max(Lep2SIP, max(Lep3SIP, Lep4SIP))) << endl;
*/
                float strdraw = Z1Mass;

                float strdraw_Txy = ((KalmanCandVtx_x - OfflinePrimaryVtx_x)*cos(ZZPhi) + (KalmanCandVtx_x - OfflinePrimaryVtx_x)*sin(ZZPhi))*10000.0*ZZMass / ZZPt;
                float strdraw_mZZ = ZZMass;
                float strdraw_mZ1 = Z1Mass;
                float strdraw_mZ2 = Z2Mass;
                float strdraw_pteta_ld[2] ={
                  fabs(Lep3Eta), Lep3Pt
                };
                float strdraw_pteta_subld[2] ={
                  fabs(Lep4Eta), Lep4Pt
                };
                if (iCR!=1 && Lep3Pt<=Lep4Pt){
                  strdraw_pteta_ld[0] = fabs(Lep4Eta);
                  strdraw_pteta_ld[1] = fabs(Lep4Pt);
                  strdraw_pteta_subld[0] = fabs(Lep3Eta);
                  strdraw_pteta_subld[1] = fabs(Lep3Pt);
                }
                hsip->Fill(strdraw, wgt);
                hTxy->Fill(strdraw_Txy, wgt);
                hmZZ->Fill(strdraw_mZZ, wgt);
                if (iM1cut==0) hmZ1->Fill(strdraw_mZ1, wgt);
                hpteta->Fill(strdraw_pteta_ld[1], strdraw_pteta_ld[0], wgt);
                if (iCR!=1){
                  hmZ2->Fill(strdraw_mZ2, wgt);
                  hptetasub->Fill(strdraw_pteta_subld[1], strdraw_pteta_subld[0], wgt);
                }
              }

              frecord->WriteTObject(hTxy);
              delete hTxy;
              frecord->WriteTObject(hmZZ);
              delete hmZZ;
              if (iM1cut==0){
                frecord->WriteTObject(hmZ1);
                delete hmZ1;
              }
              if (iCR!=1){
                frecord->WriteTObject(hmZ2);
                delete hmZ2;

                frecord->WriteTObject(hptetasub);
                delete hptetasub;
              }
              frecord->WriteTObject(hpteta);
              delete hpteta;

              tout << '\t' << hsip->Integral();

              if (iCR==0){
                if (icut<nSIPCuts-1 && !(catZ2==1 || catZ2==4) && isocut!=0 && iM1cut==0){
                  int channeltype=2; // 2e2mu + 2mu2e
                  if (catZ1==0 && (catZ2==0 || catZ2==5)) channeltype=0; // 4mu
                  else if (catZ1==1 && (catZ2==2 || catZ2==3)) channeltype=1; // 4e
                  else if ((catZ1==1 && (catZ2==0 || catZ2==5))
                    || (catZ1==0 && (catZ2==2 || catZ2==3))) channeltype=2; // 2e2mu
                  else{
                    cout << "Invalid cut" << endl; channeltype=-1;
                  }
                  int osss = (catZ2<3?0:1);

                  if (channeltype>=0){
                    double integral, integralerror;
                    integral = hsip->IntegralAndError(1, hsip->GetNbinsX(), integralerror);
                    if (isocut==1 && catZ2<3){
                      integral=0; integralerror=0;
                    }
                    if (iwgt>0){
                      hYield[channeltype][osss]->AddBinContent(icut+1, integral);
                      hYield[channeltype][osss]->SetBinError(icut+1, sqrt(pow(hYield[channeltype][osss]->GetBinError(icut+1), 2)  + pow(integralerror, 2)));
                    }
                    else{
                      hCount[channeltype][osss]->AddBinContent(icut+1, integral);
                      hCount[channeltype][osss]->SetBinError(icut+1, sqrt(pow(hCount[channeltype][osss]->GetBinError(icut+1), 2)  + pow(integralerror, 2)));
                    }
                  }
                }
              }

              delete hsip;
            }
            tout << endl;
            plot1DNewOld(iCR, iwgt, isocut, iM1cut, catZ1, catZ2, "Txy", erg_tev, frecord, m4l_low, m4l_high);
            plot1DNewOld(iCR, iwgt, isocut, iM1cut, catZ1, catZ2, "mZZ", erg_tev, frecord, m4l_low, m4l_high);
            if (iM1cut==0) plot1DNewOld(iCR, iwgt, isocut, iM1cut, catZ1, catZ2, "mZ1", erg_tev, frecord, m4l_low, m4l_high);
            if (iCR!=1) plot1DNewOld(iCR, iwgt, isocut, iM1cut, catZ1, catZ2, "mZ2", erg_tev, frecord, m4l_low, m4l_high);
            plot2DNewOld(iCR, iwgt, isocut, iM1cut, catZ1, catZ2, "pTeta_leading", erg_tev, frecord, m4l_low, m4l_high);
            if (iCR!=1) plot2DNewOld(iCR, iwgt, isocut, iM1cut, catZ1, catZ2, "pTeta_subleading", erg_tev, frecord, m4l_low, m4l_high);
          }
        }
        tout << endl;
      }
      tout << endl;
    }
    tout << endl;

    for (int isocut=1; isocut<4; isocut++){
      if (iCR==1 && isocut>=2) continue;
      for (int catZ1=0; catZ1<2; catZ1++){
        for (int catZ2=0; catZ2<(iCR==0 ? 6 : 2); catZ2++){
          for (int iM1cut=0; iM1cut<3; iM1cut++){
            if (iM1cut==2 && iCR!=1) continue;
            plot2DNewOldDoubleRatio(iCR, iwgt, isocut, iM1cut, catZ1, catZ2, "pTeta_leading", erg_tev, frecord, m4l_low, m4l_high);
            if (iCR!=1) plot2DNewOldDoubleRatio(iCR, iwgt, isocut, iM1cut, catZ1, catZ2, "pTeta_subleading", erg_tev, frecord, m4l_low, m4l_high);
          }
        }
      }
    }

  }
  delete tc;

  if (iCR==0){
    for (int ff=0; ff<3; ff++){
      for (int bin=1; bin<=nSIPCuts-1; bin++){
        double count_OS = hCount[ff][0]->GetBinContent(bin);
        double yield_OS = hYield[ff][0]->GetBinContent(bin);
        double yielderror_OS = hYield[ff][0]->GetBinError(bin);

        double count_SS = hCount[ff][1]->GetBinContent(bin);
        double yield_SS = hYield[ff][1]->GetBinContent(bin);
        double yielderror_SS = hYield[ff][1]->GetBinError(bin);

        yield_OS *= count_OS / (count_OS+count_SS);
        yield_SS *= count_SS / (count_OS+count_SS);
        yielderror_OS *= count_OS / (count_OS+count_SS);
        yielderror_SS *= count_SS / (count_OS+count_SS);
        hYield[ff][0]->SetBinContent(bin, yield_OS);
        hYield[ff][0]->SetBinError(bin, yielderror_OS);
        hYield[ff][1]->SetBinContent(bin, yield_SS);
        hYield[ff][1]->SetBinError(bin, yielderror_SS);
      }
      for (int ssos=0; ssos<2; ssos++){
        //      double nocutabs = hYield[ff][ssos]->GetBinContent(2);
        //      for (int bin=1; bin<=nSIPCuts-1; bin++) hYield[ff][ssos]->SetBinContent(bin, hYield[ff][ssos]->GetBinContent(bin) / nocutabs);
        tout << "Yield scaling ratio for channel " << ff << ": " << hYield[ff][ssos]->GetBinContent(nSIPCuts-1)/hYield[ff][ssos]->GetBinContent(1) << endl;
        foutput[ff]->cd();
        foutput[ff]->WriteTObject(hYield[ff][ssos]);
        foutput[ff]->WriteTObject(hCount[ff][ssos]);
        delete hYield[ff][ssos];
        delete hCount[ff][ssos];
      }
      foutput[ff]->Close();
    }
  }
  tout.close();
  frecord->Close();

  for (int sipcut=0; sipcut<nSIPCuts; sipcut++) delete hRatio[sipcut];
  fSSinput->Close();
  ZXWeightTables_OS->Close();
  fweight->Close();
}


void count_OSSS_CR_MC(int erg_tev){
  char TREE_NAME[]="SelectedTree";

  TString erg_dir;
  erg_dir.Form("LHC_%iTeV", erg_tev);
  TString comstring;
  comstring.Form("%iTeV", erg_tev);
  char cerg[10];
  sprintf(cerg,"%iTeV", erg_tev);

  TChain* tc = new TChain(TREE_NAME);

  TString cinput_common_noSIP = user_dir_hep + "No_SIP/" + erg_dir + "/CR/";
  for (int iCR=0; iCR<nCRZLLMC; iCR++){
    TString cinput = cinput_common_noSIP;
    cinput = cinput + sample_CR_MC[iCR] + "_CRZLLTree.root";
    tc->Add(cinput);
  }
  if (tc->GetEntries() == 0){
    cout << "Could not find any files, aborting..." << endl; return;
  }

  string coutput_common = user_dir_hep + "Analysis/Auxiliary/";

  TFile* frecord;
  TString fname = "OSoverSS_";
  fname.Append("MCCR_");
  fname += comstring;
  fname.Append(".root");
  fname.Prepend(coutput_common.c_str());
  frecord = new TFile(fname, "recreate");

  TH2F* hYield[nSIPCuts];
  TH2F* hRatio[nSIPCuts];
  string strtout = coutput_common;
  strtout = strtout + "logMCCR_OSoverSS_";
  strtout = strtout + cerg;
  strtout = strtout + ".log";
  ofstream tout(strtout.c_str(), ios::out);

  const int nmZZbins=1;
  double mZZbins[nmZZbins+1]={ 100, 3000 };
//  const int nmZZbins=5;
//  double mZZbins[nmZZbins+1]={ 100, 105.6, 140.6, 170, 800, 3000 };
  for (int sipcut=0; sipcut<nSIPCuts; sipcut++){
    TString chyield = "hCR_MC_OSSSYields_";
    chyield.Append(cutLabel[sipcut]);
    TString chratio = "hCR_MC_OSSSRatios_";
    chratio.Append(cutLabel[sipcut]);
    hYield[sipcut] = new TH2F(chyield, "Yield of any SIP", nmZZbins, mZZbins, 12, 0, 12);
    hRatio[sipcut] = new TH2F(chratio, "Ratio of any SIP", nmZZbins, mZZbins, 6, 0, 6);
  }

  frecord->cd();

  TString hname = "htemp";

  float MC_weight;
  float ZXfake_weight;
  float ZXfake_weight_SS;
  float ZXfake_weight_OS[6];
  int CRflag;
  float Lep1combRelIsoPF, Lep2combRelIsoPF, Lep3combRelIsoPF, Lep4combRelIsoPF;
  bool Lep1isID, Lep2isID, Lep3isID, Lep4isID;
  short Z1ids;
  short Z2ids;
  float PFMET;

  float Lep1_Z1SIP, Lep2_Z1SIP, Lep3_Z1SIP, Lep4_Z1SIP;
  float Lep1SIP, Lep2SIP, Lep3SIP, Lep4SIP;
  int Lep1ID, Lep2ID, Lep3ID, Lep4ID;

  float Z1Mass, Z2Mass, ZZMass, ZZPt, ZZEta, ZZPhi;
  float KalmanCandVtx_x, KalmanCandVtx_y, KalmanCandVtx_z, KalmanCandVtx_chi2;
  float OfflinePrimaryVtx_x, OfflinePrimaryVtx_y, OfflinePrimaryVtx_z;

  float Lep1Pt, Lep2Pt, Lep3Pt, Lep4Pt;
  float Lep1Eta, Lep2Eta, Lep3Eta, Lep4Eta;

  tc->SetBranchAddress("MC_weight", &MC_weight);
  tc->SetBranchAddress("ZXfake_weight", &ZXfake_weight);
  tc->SetBranchAddress("ZXfake_weight_OS", ZXfake_weight_OS);
  tc->SetBranchAddress("ZXfake_weight_SS", &ZXfake_weight_SS);
  tc->SetBranchAddress("CRflag", &CRflag);

  tc->SetBranchAddress("Z1ids", &Z1ids);
  tc->SetBranchAddress("Z2ids", &Z2ids);
  tc->SetBranchAddress("Lep1ID", &Lep1ID);
  tc->SetBranchAddress("Lep1isID", &Lep1isID);
  tc->SetBranchAddress("Lep1combRelIsoPF", &Lep1combRelIsoPF);
  tc->SetBranchAddress("Lep2ID", &Lep2ID);
  tc->SetBranchAddress("Lep2isID", &Lep2isID);
  tc->SetBranchAddress("Lep2combRelIsoPF", &Lep2combRelIsoPF);
  tc->SetBranchAddress("Lep3ID", &Lep3ID);
  tc->SetBranchAddress("Lep3isID", &Lep3isID);
  tc->SetBranchAddress("Lep3combRelIsoPF", &Lep3combRelIsoPF);

  tc->SetBranchAddress("PFMET", &PFMET);
  tc->SetBranchAddress("Z1Mass", &Z1Mass);
  tc->SetBranchAddress("Z2Mass", &Z2Mass);
  tc->SetBranchAddress("ZZMass", &ZZMass);
  tc->SetBranchAddress("ZZPt", &ZZPt);
  tc->SetBranchAddress("ZZEta", &ZZEta);
  tc->SetBranchAddress("ZZPhi", &ZZPhi);
  tc->SetBranchAddress("OfflinePrimaryVtx_x", &OfflinePrimaryVtx_x);
  tc->SetBranchAddress("OfflinePrimaryVtx_y", &OfflinePrimaryVtx_y);
  tc->SetBranchAddress("OfflinePrimaryVtx_z", &OfflinePrimaryVtx_z);
  tc->SetBranchAddress("KalmanCandVtx_x", &KalmanCandVtx_x);
  tc->SetBranchAddress("KalmanCandVtx_y", &KalmanCandVtx_y);
  tc->SetBranchAddress("KalmanCandVtx_z", &KalmanCandVtx_z);

  tc->SetBranchAddress("Lep1Pt", &Lep1Pt);
  tc->SetBranchAddress("Lep2Pt", &Lep2Pt);
  tc->SetBranchAddress("Lep3Pt", &Lep3Pt);
  tc->SetBranchAddress("Lep1Eta", &Lep1Eta);
  tc->SetBranchAddress("Lep2Eta", &Lep2Eta);
  tc->SetBranchAddress("Lep3Eta", &Lep3Eta);


  tc->SetBranchAddress("KalmanCandVtx_chi2", &KalmanCandVtx_chi2);
  tc->SetBranchAddress("Lep1_Z1SIP", &Lep1_Z1SIP);
  tc->SetBranchAddress("Lep2_Z1SIP", &Lep2_Z1SIP);
  tc->SetBranchAddress("Lep3_Z1SIP", &Lep3_Z1SIP);
  tc->SetBranchAddress("Lep1SIP", &Lep1SIP);
  tc->SetBranchAddress("Lep2SIP", &Lep2SIP);
  tc->SetBranchAddress("Lep3SIP", &Lep3SIP);

  tc->SetBranchAddress("Lep4ID", &Lep4ID);
  tc->SetBranchAddress("Lep4isID", &Lep4isID);
  tc->SetBranchAddress("Lep4combRelIsoPF", &Lep4combRelIsoPF);
  tc->SetBranchAddress("Lep4_Z1SIP", &Lep4_Z1SIP);
  tc->SetBranchAddress("Lep4SIP", &Lep4SIP);
  tc->SetBranchAddress("Lep4Pt", &Lep4Pt);
  tc->SetBranchAddress("Lep4Eta", &Lep4Eta);


  for (int iwgt=0; iwgt<2; iwgt++){
    for (int isocut=0; isocut<4; isocut++){
      for (int iM1cut=0; iM1cut<1; iM1cut++){
        if (iwgt==0) tout << "Un-weighted";
        else tout << "Weighted";
        if (isocut>=1) tout << ", " << strLepIsoCut_label[isocut-1];
        tout << ", " << Z1masswidth_label[iM1cut] << endl;
        tout << "\t"; for (int icut=0; icut<nSIPCuts; icut++) tout << cutLabel[icut] << '\t';
        tout << endl;
        for (int catZ1=0; catZ1<2; catZ1++){
          for (int catZ2=0; catZ2<6; catZ2++){
            tout << strZ1Category_label[catZ1] << " ";
            tout << strZ2Category_label[catZ2];
            for (int icut=0; icut<nSIPCuts; icut++){
              TString chsip = hname;
              TH1F* hsip = new TH1F(chsip, "", nmZZbins, mZZbins);
              hsip->Sumw2();

              for (int ev=0; ev<tc->GetEntries(); ev++){
                tc->GetEntry(ev);

                int option[7]={
                  0, icut, catZ1, catZ2, iwgt, isocut, iM1cut
                };
                int LepID[4]={
                  Lep1ID, Lep2ID, Lep3ID, Lep4ID
                };
                float Lep_Z1SIP[4]={
                  Lep1_Z1SIP, Lep2_Z1SIP, Lep3_Z1SIP, Lep4_Z1SIP
                };
                float LepSIP[4]={
                  Lep1SIP, Lep2SIP, Lep3SIP, Lep4SIP
                };

                float wgt = applyCRselection(
                  option,

                  0, 10000,

                  CRflag,
                  Lep1combRelIsoPF, Lep2combRelIsoPF, Lep3combRelIsoPF, Lep4combRelIsoPF,
                  Lep1isID, Lep2isID, Lep3isID, Lep4isID,
                  Z1ids,
                  Z2ids,
                  PFMET,

                  Lep_Z1SIP,
                  KalmanCandVtx_chi2,
                  LepSIP,
                  LepID,

                  Z1Mass, ZZMass,

                  ZXfake_weight,
                  ZXfake_weight_SS,
                  ZXfake_weight_OS

                  );

                if (wgt==0) continue;
                wgt *= MC_weight;
                if (erg_tev==7) wgt*= 5.051;
                else if (erg_tev==8) wgt*= 19.712;

                float strdraw = ZZMass;

                hsip->Fill(strdraw, wgt);
              }

              tout << '\t' << hsip->Integral();

              if (isocut!=0 && iM1cut==0 && iwgt>0){
                int channeltype=catZ2+6*catZ1;

                if (!(isocut==1 && catZ2<3)){
                  for (int binx=1; binx<=hsip->GetNbinsX()+1; binx++){
                    double bincontent = hsip->GetBinContent(binx);
                    double binerror = hsip->GetBinError(binx);
                    hYield[icut]->SetBinContent(binx, channeltype+1, bincontent + hYield[icut]->GetBinContent(binx, channeltype+1));
                    hYield[icut]->SetBinError(binx, channeltype+1, sqrt(pow(hYield[icut]->GetBinError(binx, channeltype+1), 2)  + pow(binerror, 2)));
                  }
                }
              }

              delete hsip;
            }
            tout << endl;
          }
        }
        tout << endl;
      }
      tout << endl;
    }
    tout << endl;
  }
  delete tc;

  for (int icut=0; icut<nSIPCuts; icut++){
    for (int binx=1; binx<=hYield[icut]->GetNbinsX()+1; binx++){
      cout << "Binx: " << binx << endl;
      for (int biny=1; biny<=hYield[icut]->GetNbinsY(); biny++){
        if (biny>=4 && biny<=6) continue;
        if (biny>=10 && biny<=12) continue;
        int oppositebiny = 6-biny+1;
        if (biny>6) oppositebiny = 12-biny+7;
        int ratiobiny = biny;
        if (biny>6) ratiobiny = biny-3;


        double yield_OS = hYield[icut]->GetBinContent(binx, biny);
        double yielderror_OS = hYield[icut]->GetBinError(binx, biny);
        double yield_SS = hYield[icut]->GetBinContent(binx, oppositebiny);
        double yielderror_SS = hYield[icut]->GetBinError(binx, oppositebiny);
        double ratio=1;
        if (yield_SS!=0) ratio = yield_OS / yield_SS;
        double error = 0;
        if (yield_OS!=0) error += pow(yielderror_OS/yield_OS, 2);
        if (yield_SS!=0) error += pow(yielderror_SS/yield_SS, 2);
        error = ratio * sqrt(error);
        cout << "Matched " << biny << " to " << oppositebiny << ". Recording bin " << ratiobiny << endl;
        cout << yield_OS << '/' << yield_SS << " = " << ratio << " +- " << error << endl;

        hRatio[icut]->SetBinContent(binx, ratiobiny, ratio);
        hRatio[icut]->SetBinError(binx, ratiobiny, error);
      }
    }
    frecord->cd();
    frecord->WriteTObject(hYield[icut]);
    frecord->WriteTObject(hRatio[icut]);
    delete hYield[icut];
    delete hRatio[icut];
  }
  tout.close();
  frecord->Close();
}


TString produceCRname(int iCR, int iwgt, int isocut, int iM1cut, int catZ1, int catZ2, int icut, char* cappend){
  TString name;
  if (iCR==0) name = "ZLL";
  else name = "ZL";

  if (catZ1==0) name.Append("_Zmumu");
  else name.Append("_Zee");
  if (iCR==0){
    if (catZ2==0) name.Append("_mumuOS");
    if (catZ2==1) name.Append("_emuOS");
    if (catZ2==2) name.Append("_eeOS");
    if (catZ2==3) name.Append("_eeSS");
    if (catZ2==4) name.Append("_emuSS");
    if (catZ2==5) name.Append("_mumuSS");
  }
  else if (iCR==1){
    if (catZ2==0) name.Append("_mu");
    if (catZ2==1) name.Append("_e");
  }

  if (iwgt==0) name.Append("_Unweighted");
  else if (iwgt==1) name.Append("_Weighted");

  name.Append("_");
  name += cutNames[icut];

  if (isocut==0) name.Append("_AllLoose");
  else{ name.Append("_"); name += strLepIsoCut_label[isocut-1]; }

  if (iM1cut==1) name.Append("_ZPeak");
  else if (iM1cut==2) name.Append("_ZGPeak");

  name.Append("_");
  name.Append(cappend);

  return name;
}

void plot1DNewOld(int iCR, int iwgt, int isocut, int iM1cut, int catZ1, int catZ2, char* cappend, int erg_tev, TFile* frecord, double m4l_low, double m4l_high){
  TString comstring;
  comstring.Form("%iTeV", erg_tev);
  TString coutput_common = user_dir_hep + "Analysis/Plots/CR_FakeRates_own/Distributions/";
  coutput_common.Append(comstring);
  coutput_common.Append("/");
  if (iCR==0){
    coutput_common = coutput_common + "ZLL/";
    coutput_common.Append(Form("m4l_%.1f_%.1f/", m4l_low, m4l_high));
  }
  else coutput_common = coutput_common +  "ZL/";
  if (isocut==0) coutput_common.Append("AllLoose");
  else coutput_common.Append(strLepIsoCut_label[isocut-1]);
  coutput_common.Append("/");
  if (iM1cut==0) coutput_common.Append("m1_All/");
  else if (iM1cut==1) coutput_common.Append("m1_ZPeak/");
  else if (iM1cut==2) coutput_common.Append("m1_ZGPeak/");

  TString canvasdir = coutput_common;
  canvasdir.Append(Form("%s/", cappend));

  TString mkdirCommand = "mkdir -p ";
  mkdirCommand.Append(canvasdir);
  gSystem->Exec(mkdirCommand);

  TString hname[2] ={
    produceCRname(iCR, iwgt, isocut, iM1cut, catZ1, catZ2, 4, cappend),
    produceCRname(iCR, iwgt, isocut, iM1cut, catZ1, catZ2, 0, cappend)
  };
  TString hlabel[2] ={
    "New cut", "SIP_{3D}<4"
  };
  TH1F* htemp[2] ={
    (TH1F*)frecord->Get(hname[0]),
    (TH1F*)frecord->Get(hname[1])
  };

  if (iCR==1){
    htemp[0]->Scale(0.001);
    htemp[1]->Scale(0.001);
    htemp[0]->SetYTitle(Form("%s / 1000", htemp[0]->GetYaxis()->GetTitle()));
    htemp[1]->SetYTitle(Form("%s / 1000", htemp[1]->GetYaxis()->GetTitle()));
  }

  TString identifier_label = "#frac{";
  identifier_label = identifier_label + strZ1Category_label[catZ1];
  if (iCR!=1) identifier_label = identifier_label + strZ2Category_label[catZ2];
  else identifier_label = identifier_label + strLep3Category_label[catZ2];
  identifier_label = identifier_label + "}{";
  if (isocut>1) identifier_label = identifier_label + strLepIsoCut_label[isocut-1];
  else if (isocut==1 && iCR==0) identifier_label = identifier_label + "4p";
  else if (isocut==1 && iCR==1) identifier_label = identifier_label + "3p";
  if (isocut!=0) identifier_label = identifier_label + ", ";
  identifier_label = identifier_label + Z1masswidth_label[iM1cut] + "}";

  TString canvasname = produceCRname(iCR, iwgt, isocut, iM1cut, catZ1, catZ2, 5, cappend);
  canvasname = canvasname + "_" + comstring;

  frecord->cd();

  TCanvas* c2D = new TCanvas(canvasname, "", 8, 30, 800, 800);
  c2D->cd();
  gStyle->SetOptStat(0);
  c2D->SetFillColor(0);
  c2D->SetBorderMode(0);
  c2D->SetBorderSize(2);
  c2D->SetTickx(1);
  c2D->SetTicky(1);
  c2D->SetLeftMargin(0.17);
  c2D->SetRightMargin(0.05);
  c2D->SetTopMargin(0.07);
  c2D->SetBottomMargin(0.13);
  c2D->SetFrameFillStyle(0);
  c2D->SetFrameBorderMode(0);
  c2D->SetFrameFillStyle(0);
  c2D->SetFrameBorderMode(0);

  TLegend *l2D = new TLegend(0.76, 0.74, 0.98, 0.90);
  l2D->SetBorderSize(0);
  l2D->SetTextFont(42);
  l2D->SetTextSize(0.03);
  l2D->SetLineColor(1);
  l2D->SetLineStyle(1);
  l2D->SetLineWidth(1);
  l2D->SetFillColor(0);
  l2D->SetFillStyle(0);

  TPaveText* pt = new TPaveText(0.15, 0.93, 0.85, 1, "brNDC");
  pt->SetBorderSize(0);
  pt->SetFillStyle(0);
  pt->SetTextAlign(12);
  pt->SetTextFont(42);
  pt->SetTextSize(0.045);
  TText* text = pt->AddText(0.025, 0.45, "#font[61]{CMS}");
  text->SetTextSize(0.044);
  text = pt->AddText(0.165, 0.42, "#font[52]{Preliminary}");
  text->SetTextSize(0.0315);
  if (erg_tev==8) text = pt->AddText(0.837, 0.45, "#font[42]{19.7 fb^{-1} (8 TeV)}");
  else if (erg_tev==7) text = pt->AddText(0.837, 0.45, "#font[42]{ 5.1 fb^{-1} (7 TeV)}");
  text->SetTextSize(0.0315);

  htemp[0]->GetYaxis()->SetRangeUser(0, (TMath::Max(htemp[0]->GetMaximum() + htemp[0]->GetBinError(htemp[0]->GetMaximumBin()), htemp[1]->GetMaximum() + htemp[1]->GetBinError(htemp[1]->GetMaximumBin()))) * 1.5);
  htemp[0]->GetXaxis()->SetNdivisions(505);
  htemp[0]->GetXaxis()->SetLabelFont(42);
  htemp[0]->GetXaxis()->SetLabelOffset(0.007);
  htemp[0]->GetXaxis()->SetLabelSize(0.04);
  htemp[0]->GetXaxis()->SetTitleSize(0.06);
  htemp[0]->GetXaxis()->SetTitleOffset(0.9);
  htemp[0]->GetXaxis()->SetTitleFont(42);
  htemp[0]->GetYaxis()->SetNdivisions(505);
  htemp[0]->GetYaxis()->SetLabelFont(42);
  htemp[0]->GetYaxis()->SetLabelOffset(0.007);
  htemp[0]->GetYaxis()->SetLabelSize(0.04);
  htemp[0]->GetYaxis()->SetTitleSize(0.06);
  htemp[0]->GetYaxis()->SetTitleOffset(1.1);
  htemp[0]->GetYaxis()->SetTitleFont(42);

  htemp[0]->SetLineWidth(2);
  htemp[1]->SetLineWidth(2);
  htemp[0]->SetLineColor(kBlue);
  htemp[1]->SetLineColor(kRed);

  l2D->AddEntry(htemp[0], hlabel[0], "lf");
  l2D->AddEntry(htemp[1], hlabel[1], "lf");

  htemp[0]->Draw("hist");
  htemp[1]->Draw("histsame");

  TH1F* hclone[2]={
    (TH1F*)htemp[0]->Clone("h1"),
    (TH1F*)htemp[1]->Clone("h2")
  };

  hclone[0]->SetFillColor(kBlue);
  hclone[1]->SetFillColor(kRed);
  hclone[0]->SetFillStyle(3004);
  hclone[1]->SetFillStyle(3005);
  hclone[0]->Draw("e2same");
  hclone[1]->Draw("e2same");
  l2D->Draw("same");
  pt->Draw();

  TPaveText *pt10 = new TPaveText(0.21, 0.84, 0.40, 0.92, "brNDC");
  pt10->SetBorderSize(0);
  pt10->SetTextAlign(12);
  pt10->SetTextSize(0.03);
  pt10->SetFillStyle(0);
  pt10->SetTextFont(42);
  TText* text10 = pt10->AddText(0.01, 0.01, identifier_label);
  pt10->Draw();

  frecord->WriteTObject(c2D);

  canvasname.Prepend(canvasdir);
  TString canvasname_pdf = canvasname;
  TString canvasname_eps = canvasname;
  TString canvasname_png = canvasname;
  TString canvasname_root = canvasname;
  TString canvasname_c = canvasname;
  canvasname_pdf.Append(".pdf");
  canvasname_eps.Append(".eps");
  canvasname_png.Append(".png");
  canvasname_root.Append(".root");
  canvasname_c.Append(".C");
  c2D->SaveAs(canvasname_pdf);
  c2D->SaveAs(canvasname_eps);
  c2D->SaveAs(canvasname_png);
  c2D->SaveAs(canvasname_root);
  c2D->SaveAs(canvasname_c);

  delete hclone[0];
  delete hclone[1];
  c2D->Close();
}

void plot2DNewOld(int iCR, int iwgt, int isocut, int iM1cut, int catZ1, int catZ2, char* cappend, int erg_tev, TFile* frecord, double m4l_low, double m4l_high){
  TString comstring;
  comstring.Form("%iTeV", erg_tev);
  TString coutput_common = user_dir_hep + "Analysis/Plots/CR_FakeRates_own/Distributions/";
  coutput_common.Append(comstring);
  coutput_common.Append("/");
  if (iCR==0){
    coutput_common = coutput_common + "ZLL/";
    coutput_common.Append(Form("m4l_%.1f_%.1f/", m4l_low, m4l_high));
  }
  else coutput_common = coutput_common +  "ZL/";
  if (isocut==0) coutput_common.Append("AllLoose");
  else coutput_common.Append(strLepIsoCut_label[isocut-1]);
  coutput_common.Append("/");
  if (iM1cut==0) coutput_common.Append("m1_All/");
  else if (iM1cut==1) coutput_common.Append("m1_ZPeak/");
  else if (iM1cut==2) coutput_common.Append("m1_ZGPeak/");

  TString canvasdir = coutput_common;
  canvasdir.Append(Form("%s/", cappend));

  TString mkdirCommand = "mkdir -p ";
  mkdirCommand.Append(canvasdir);
  gSystem->Exec(mkdirCommand);

  TString hname[2] ={
    produceCRname(iCR, iwgt, isocut, iM1cut, catZ1, catZ2, 4, cappend),
    produceCRname(iCR, iwgt, isocut, iM1cut, catZ1, catZ2, 0, cappend)
  };
  TString hlabel[4] ={
    "New cut (All)", "SIP_{3D}<4 (All)",
    "New cut (Endcap)", "SIP_{3D}<4 (Endcap)"
  };
  TH2F* htemp[2] ={
    (TH2F*)frecord->Get(hname[0]),
    (TH2F*)frecord->Get(hname[1])
  };

  if (iCR==1){
    htemp[0]->Scale(0.001);
    htemp[1]->Scale(0.001);
    htemp[0]->GetZaxis()->SetTitle(Form("%s / 1000", htemp[0]->GetZaxis()->GetTitle()));
    htemp[1]->GetZaxis()->SetTitle(Form("%s / 1000", htemp[1]->GetZaxis()->GetTitle()));
  }

  TH1F* htemp_1D[4] ={
    (TH1F*)htemp[0]->ProjectionX("hnew0_12", 0, -1, "e"),
    (TH1F*)htemp[1]->ProjectionX("hnew1_12", 0, -1, "e"),
    (TH1F*)htemp[0]->ProjectionX("hnew0_2", 0, -1, "e"),
    (TH1F*)htemp[1]->ProjectionX("hnew1_2", 0, -1, "e")
  };
  TH1F* hclone[4];
  for (int hh=0; hh<2; hh++){
    for (int binx=0; binx<=htemp_1D[hh]->GetNbinsX()+1; binx++){
      double inclusive = htemp_1D[hh]->GetBinContent(binx);
      double inclusive_err = htemp_1D[hh]->GetBinError(binx);
      double exclusive = htemp[hh]->GetBinContent(binx, 1);
      double exclusive_err = htemp[hh]->GetBinError(binx, 1);
      double bincontent = inclusive - exclusive;
      double binerror = sqrt(pow(inclusive_err, 2) - pow(exclusive_err, 2));
      htemp_1D[2+hh]->SetBinContent(binx, bincontent);
      htemp_1D[2+hh]->SetBinError(binx, binerror);
    }
  }

  TString identifier_label = "#frac{";
  identifier_label = identifier_label + strZ1Category_label[catZ1];
  if (iCR!=1) identifier_label = identifier_label + strZ2Category_label[catZ2];
  else identifier_label = identifier_label + strLep3Category_label[catZ2];
  identifier_label = identifier_label + "}{";
  if (isocut>1) identifier_label = identifier_label + strLepIsoCut_label[isocut-1];
  else if (isocut==1 && iCR==0) identifier_label = identifier_label + "4p";
  else if (isocut==1 && iCR==1) identifier_label = identifier_label + "3p";
  if (isocut!=0) identifier_label = identifier_label + ", ";
  identifier_label = identifier_label + Z1masswidth_label[iM1cut] + "}";

  TString canvasname = produceCRname(iCR, iwgt, isocut, iM1cut, catZ1, catZ2, 5, cappend);
  canvasname = canvasname + "_" + comstring;

  frecord->cd();

  TCanvas* c2D = new TCanvas(canvasname, "", 8, 30, 800, 800);
  c2D->cd();
  gStyle->SetOptStat(0);
  c2D->SetFillColor(0);
  c2D->SetBorderMode(0);
  c2D->SetBorderSize(2);
  c2D->SetTickx(1);
  c2D->SetTicky(1);
  c2D->SetLeftMargin(0.17);
  c2D->SetRightMargin(0.05);
  c2D->SetTopMargin(0.07);
  c2D->SetBottomMargin(0.13);
  c2D->SetFrameFillStyle(0);
  c2D->SetFrameBorderMode(0);
  c2D->SetFrameFillStyle(0);
  c2D->SetFrameBorderMode(0);

  TLegend *l2D = new TLegend(0.61, 0.74, 0.98, 0.90);
  l2D->SetBorderSize(0);
  l2D->SetTextFont(42);
  l2D->SetTextSize(0.03);
  l2D->SetLineColor(1);
  l2D->SetLineStyle(1);
  l2D->SetLineWidth(1);
  l2D->SetFillColor(0);
  l2D->SetFillStyle(0);

  TPaveText* pt = new TPaveText(0.15, 0.93, 0.85, 1, "brNDC");
  pt->SetBorderSize(0);
  pt->SetFillStyle(0);
  pt->SetTextAlign(12);
  pt->SetTextFont(42);
  pt->SetTextSize(0.045);
  TText* text = pt->AddText(0.025, 0.45, "#font[61]{CMS}");
  text->SetTextSize(0.044);
  text = pt->AddText(0.165, 0.42, "#font[52]{Preliminary}");
  text->SetTextSize(0.0315);
  if (erg_tev==8) text = pt->AddText(0.837, 0.45, "#font[42]{19.7 fb^{-1} (8 TeV)}");
  else if (erg_tev==7) text = pt->AddText(0.837, 0.45, "#font[42]{ 5.1 fb^{-1} (7 TeV)}");
  text->SetTextSize(0.0315);

  for (int hh=0; hh<2; hh++){
    htemp_1D[2*hh]->GetYaxis()->SetRangeUser(0, (TMath::Max(htemp_1D[2*hh]->GetMaximum() + htemp_1D[2*hh]->GetBinError(htemp_1D[2*hh]->GetMaximumBin()), htemp_1D[2*hh+1]->GetMaximum() + htemp_1D[2*hh+1]->GetBinError(htemp_1D[2*hh+1]->GetMaximumBin()))) * 1.5);
    htemp_1D[2*hh]->GetXaxis()->SetNdivisions(505);
    htemp_1D[2*hh]->GetXaxis()->SetLabelFont(42);
    htemp_1D[2*hh]->GetXaxis()->SetLabelOffset(0.007);
    htemp_1D[2*hh]->GetXaxis()->SetLabelSize(0.04);
    htemp_1D[2*hh]->GetXaxis()->SetTitleSize(0.06);
    htemp_1D[2*hh]->GetXaxis()->SetTitleOffset(0.9);
    htemp_1D[2*hh]->GetXaxis()->SetTitleFont(42);
    htemp_1D[2*hh]->GetYaxis()->SetNdivisions(505);
    htemp_1D[2*hh]->GetYaxis()->SetLabelFont(42);
    htemp_1D[2*hh]->GetYaxis()->SetLabelOffset(0.007);
    htemp_1D[2*hh]->GetYaxis()->SetLabelSize(0.04);
    htemp_1D[2*hh]->GetYaxis()->SetTitleSize(0.06);
    htemp_1D[2*hh]->GetYaxis()->SetTitleOffset(1.1);
    htemp_1D[2*hh]->GetYaxis()->SetTitleFont(42);

    htemp_1D[2*hh]->SetYTitle(htemp[0]->GetZaxis()->GetTitle());
    htemp_1D[2*hh+1]->SetYTitle(htemp[0]->GetZaxis()->GetTitle());

    htemp_1D[2*hh]->SetLineWidth(2);
    htemp_1D[2*hh+1]->SetLineWidth(2);
    htemp_1D[2*hh]->SetLineColor(kBlue);
    htemp_1D[2*hh+1]->SetLineColor(kRed);

    if (hh==1){
      htemp_1D[2*hh]->SetLineStyle(7);
      htemp_1D[2*hh+1]->SetLineStyle(7);
    }

    l2D->AddEntry(htemp_1D[2*hh], hlabel[2*hh], "lf");
    l2D->AddEntry(htemp_1D[2*hh+1], hlabel[2*hh+1], "lf");

    if (hh==0) htemp_1D[2*hh]->Draw("hist");
    else htemp_1D[2*hh]->Draw("histsame");
    htemp_1D[2*hh+1]->Draw("histsame");

    hclone[2*hh]=(TH1F*)htemp_1D[2*hh]->Clone(Form("h%i", 2*hh));
    hclone[2*hh+1]=(TH1F*)htemp_1D[2*hh+1]->Clone(Form("h%i", 2*hh+1));

    hclone[2*hh]->SetFillColor(kBlue);
    hclone[2*hh+1]->SetFillColor(kRed);
    hclone[2*hh]->SetFillStyle(3004);
    hclone[2*hh+1]->SetFillStyle(3005);
    hclone[2*hh]->Draw("e2same");
    hclone[2*hh+1]->Draw("e2same");
  }
  l2D->Draw("same");
  pt->Draw();

  TPaveText *pt10 = new TPaveText(0.21, 0.84, 0.40, 0.92, "brNDC");
  pt10->SetBorderSize(0);
  pt10->SetTextAlign(12);
  pt10->SetTextSize(0.03);
  pt10->SetFillStyle(0);
  pt10->SetTextFont(42);
  TText* text10 = pt10->AddText(0.01, 0.01, identifier_label);
  pt10->Draw();

  frecord->WriteTObject(c2D);

  canvasname.Prepend(canvasdir);
  TString canvasname_pdf = canvasname;
  TString canvasname_eps = canvasname;
  TString canvasname_png = canvasname;
  TString canvasname_root = canvasname;
  TString canvasname_c = canvasname;
  canvasname_pdf.Append(".pdf");
  canvasname_eps.Append(".eps");
  canvasname_png.Append(".png");
  canvasname_root.Append(".root");
  canvasname_c.Append(".C");
  c2D->SaveAs(canvasname_pdf);
  c2D->SaveAs(canvasname_eps);
  c2D->SaveAs(canvasname_png);
  c2D->SaveAs(canvasname_root);
  c2D->SaveAs(canvasname_c);

  for (int hh=0; hh<4; hh++){
    delete hclone[hh];
    delete htemp_1D[hh];
  }
  c2D->Close();
}

void plot2DNewOldDoubleRatio(int iCR, int iwgt, int isocut, int iM1cut, int catZ1, int catZ2, char* cappend, int erg_tev, TFile* frecord, double m4l_low, double m4l_high){
  TString comstring;
  comstring.Form("%iTeV", erg_tev);
  TString coutput_common = user_dir_hep + "Analysis/Plots/CR_FakeRates_own/DoubleRatios/";
  coutput_common.Append(comstring);
  coutput_common.Append("/");
  if (iCR==0){
    coutput_common = coutput_common + "ZLL/";
    coutput_common.Append(Form("m4l_%.1f_%.1f/", m4l_low, m4l_high));
  }
  else coutput_common = coutput_common +  "ZL/";
  if (isocut==0) return;
  else coutput_common.Append(strLepIsoCut_label[isocut-1]);
  coutput_common.Append("/");
  if (iM1cut==0) coutput_common.Append("m1_All/");
  else if (iM1cut==1) coutput_common.Append("m1_ZPeak/");
  else if (iM1cut==2) coutput_common.Append("m1_ZGPeak/");

  TString canvasdir = coutput_common;
  canvasdir.Append(Form("%s/", cappend));

  TString mkdirCommand = "mkdir -p ";
  mkdirCommand.Append(canvasdir);
  gSystem->Exec(mkdirCommand);

  TString hname[6] ={
    produceCRname(iCR, iwgt, isocut, iM1cut, catZ1, catZ2, 4, cappend),
    produceCRname(iCR, iwgt, isocut, iM1cut, catZ1, catZ2, 0, cappend),
    produceCRname(iCR, iwgt, isocut, iM1cut, catZ1, catZ2, 5, cappend),
    produceCRname(iCR, iwgt, 0, iM1cut, catZ1, catZ2, 4, cappend),
    produceCRname(iCR, iwgt, 0, iM1cut, catZ1, catZ2, 0, cappend),
    produceCRname(iCR, iwgt, 0, iM1cut, catZ1, catZ2, 5, cappend)
  };
  TString hlabel[2][2] ={
    { "New cut (Barrel)", "New cut (Endcap)" },
    { "SIP_{3D}<4 (Barrel)", "SIP_{3D}<4 (Endcap)" }
  };
  TH2F* htemp[6];
  for (int hh=0; hh<6; hh++){
    htemp[hh] = (TH2F*)frecord->Get(hname[hh]);
    if (iCR==1){
      htemp[hh]->Scale(0.001);
      htemp[hh]->GetZaxis()->SetTitle(Form("%s / 1000", htemp[hh]->GetZaxis()->GetTitle()));
    }
  }

  double maxplot=0;
  TH1F* htemp_1D[2][2];
  for (int hh=0; hh<2; hh++){
    for (int it=0; it<2; it++){
      TString hrationame = Form("%s_%s_ratio", htemp[hh]->GetName(), "Barrel");
      if (it==1) hrationame = Form("%s_%s_ratio", htemp[hh]->GetName(), "Endcap");
      htemp_1D[hh][it] = (TH1F*)htemp[hh]->ProjectionX(hrationame);
      htemp_1D[hh][it]->SetXTitle(htemp[hh]->GetXaxis()->GetTitle());
      htemp_1D[hh][it]->SetYTitle("Fake rate");
      if (iCR==0) htemp_1D[hh][it]->SetYTitle("Ratio");
      for (int bin=0; bin<=htemp[hh]->GetNbinsX()+1; bin++){
        double tightval[2] ={ htemp[hh]->GetBinContent(bin, it+1), htemp[hh]->GetBinError(bin, it+1) };
        double looseval[2] ={ htemp[hh+3]->GetBinContent(bin, it+1), htemp[hh+3]->GetBinError(bin, it+1) };

        double commontightval[2] ={ htemp[2]->GetBinContent(bin, it+1), htemp[2]->GetBinError(bin, it+1) };
        double commonlooseval[2] ={ htemp[5]->GetBinContent(bin, it+1), htemp[5]->GetBinError(bin, it+1) };

        double bincontent = 0;
        double binerror = 0;
        if (looseval[0]>0){
          bincontent = tightval[0] / looseval[0];
          double binerror_num = sqrt(pow(tightval[1]*(looseval[0] - tightval[0]), 2) + pow(tightval[0]*sqrt(pow(looseval[1], 2) - pow(tightval[1], 2)), 2));
          binerror = binerror_num / pow(looseval[0], 2);
        }
        htemp_1D[hh][it]->SetBinContent(bin, bincontent);
        htemp_1D[hh][it]->SetBinError(bin, binerror);

        if (maxplot<(bincontent + binerror))maxplot=(bincontent + binerror);
      }
    }
  }
  //  cout << maxplot << endl;
  TH1F* hclone[2][2];

  TString identifier_label = "#frac{";
  identifier_label = identifier_label + strZ1Category_label[catZ1];
  if (iCR!=1) identifier_label = identifier_label + strZ2Category_label[catZ2];
  else identifier_label = identifier_label + strLep3Category_label[catZ2];
  identifier_label = identifier_label + "}{";
  if (isocut>1) identifier_label = identifier_label + strLepIsoCut_label[isocut-1];
  else if (isocut==1 && iCR==0) identifier_label = identifier_label + "4p";
  else if (isocut==1 && iCR==1) identifier_label = identifier_label + "3p";
  if (isocut!=0) identifier_label = identifier_label + ", ";
  identifier_label = identifier_label + Z1masswidth_label[iM1cut] + "}";

  TString canvasname = produceCRname(iCR, iwgt, isocut, iM1cut, catZ1, catZ2, 5, cappend);
  canvasname = canvasname + "_" + comstring + "_RatioComparison";

  frecord->cd();

  TCanvas* c2D = new TCanvas(canvasname, "", 8, 30, 800, 800);
  c2D->cd();
  gStyle->SetOptStat(0);
  c2D->SetFillColor(0);
  c2D->SetBorderMode(0);
  c2D->SetBorderSize(2);
  c2D->SetTickx(1);
  c2D->SetTicky(1);
  c2D->SetLeftMargin(0.17);
  c2D->SetRightMargin(0.05);
  c2D->SetTopMargin(0.07);
  c2D->SetBottomMargin(0.13);
  c2D->SetFrameFillStyle(0);
  c2D->SetFrameBorderMode(0);
  c2D->SetFrameFillStyle(0);
  c2D->SetFrameBorderMode(0);

  TLegend *l2D = new TLegend(0.61, 0.74, 0.98, 0.90);
  l2D->SetBorderSize(0);
  l2D->SetTextFont(42);
  l2D->SetTextSize(0.03);
  l2D->SetLineColor(1);
  l2D->SetLineStyle(1);
  l2D->SetLineWidth(1);
  l2D->SetFillColor(0);
  l2D->SetFillStyle(0);

  TPaveText* pt = new TPaveText(0.15, 0.93, 0.85, 1, "brNDC");
  pt->SetBorderSize(0);
  pt->SetFillStyle(0);
  pt->SetTextAlign(12);
  pt->SetTextFont(42);
  pt->SetTextSize(0.045);
  TText* text = pt->AddText(0.025, 0.45, "#font[61]{CMS}");
  text->SetTextSize(0.044);
  text = pt->AddText(0.165, 0.42, "#font[52]{Preliminary}");
  text->SetTextSize(0.0315);
  if (erg_tev==8) text = pt->AddText(0.837, 0.45, "#font[42]{19.7 fb^{-1} (8 TeV)}");
  else if (erg_tev==7) text = pt->AddText(0.837, 0.45, "#font[42]{ 5.1 fb^{-1} (7 TeV)}");
  text->SetTextSize(0.0315);

  for (int hh=0; hh<2; hh++){
    htemp_1D[hh][0]->GetYaxis()->SetRangeUser(0, maxplot * 1.5);
    htemp_1D[hh][0]->GetXaxis()->SetNdivisions(505);
    htemp_1D[hh][0]->GetXaxis()->SetLabelFont(42);
    htemp_1D[hh][0]->GetXaxis()->SetLabelOffset(0.007);
    htemp_1D[hh][0]->GetXaxis()->SetLabelSize(0.04);
    htemp_1D[hh][0]->GetXaxis()->SetTitleSize(0.06);
    htemp_1D[hh][0]->GetXaxis()->SetTitleOffset(0.9);
    htemp_1D[hh][0]->GetXaxis()->SetTitleFont(42);
    htemp_1D[hh][0]->GetYaxis()->SetNdivisions(505);
    htemp_1D[hh][0]->GetYaxis()->SetLabelFont(42);
    htemp_1D[hh][0]->GetYaxis()->SetLabelOffset(0.007);
    htemp_1D[hh][0]->GetYaxis()->SetLabelSize(0.04);
    htemp_1D[hh][0]->GetYaxis()->SetTitleSize(0.06);
    htemp_1D[hh][0]->GetYaxis()->SetTitleOffset(1.1);
    htemp_1D[hh][0]->GetYaxis()->SetTitleFont(42);


    htemp_1D[hh][0]->SetLineWidth(2);
    htemp_1D[hh][1]->SetLineWidth(2);

    if (hh==0){
      htemp_1D[hh][0]->SetLineColor(kBlue);
      htemp_1D[hh][1]->SetLineColor(kBlue);
    }
    else if (hh==1){
      htemp_1D[hh][0]->SetLineColor(kRed);
      htemp_1D[hh][1]->SetLineColor(kRed);
    }
    htemp_1D[hh][1]->SetLineStyle(7);

    frecord->WriteTObject(htemp_1D[hh][0]);
    frecord->WriteTObject(htemp_1D[hh][1]);

    l2D->AddEntry(htemp_1D[hh][0], hlabel[hh][0], "lf");
    l2D->AddEntry(htemp_1D[hh][1], hlabel[hh][1], "lf");

    if (hh==0) htemp_1D[hh][0]->Draw("hist");
    else htemp_1D[hh][0]->Draw("histsame");
    htemp_1D[hh][1]->Draw("histsame");

    hclone[hh][0]=(TH1F*)htemp_1D[hh][0]->Clone(Form("h%i_0", hh));
    hclone[hh][1]=(TH1F*)htemp_1D[hh][1]->Clone(Form("h%i_1", hh+1));

    if (hh==0){
      hclone[hh][0]->SetFillColor(kBlue);
      hclone[hh][1]->SetFillColor(kBlue);
      hclone[hh][0]->SetFillStyle(3004);
      hclone[hh][1]->SetFillStyle(3006);
    }
    else if (hh==1){
      hclone[hh][0]->SetFillColor(kRed);
      hclone[hh][1]->SetFillColor(kRed);
      hclone[hh][0]->SetFillStyle(3005);
      hclone[hh][1]->SetFillStyle(3007);
    }
    hclone[hh][0]->Draw("e2same");
    hclone[hh][1]->Draw("e2same");
  }
  l2D->Draw("same");
  pt->Draw();

  TPaveText *pt10 = new TPaveText(0.21, 0.84, 0.40, 0.92, "brNDC");
  pt10->SetBorderSize(0);
  pt10->SetTextAlign(12);
  pt10->SetTextSize(0.03);
  pt10->SetFillStyle(0);
  pt10->SetTextFont(42);
  TText* text10 = pt10->AddText(0.01, 0.01, identifier_label);
  pt10->Draw();

  frecord->WriteTObject(c2D);

  canvasname.Prepend(canvasdir);
  TString canvasname_pdf = canvasname;
  TString canvasname_eps = canvasname;
  TString canvasname_png = canvasname;
  TString canvasname_root = canvasname;
  TString canvasname_c = canvasname;
  canvasname_pdf.Append(".pdf");
  canvasname_eps.Append(".eps");
  canvasname_png.Append(".png");
  canvasname_root.Append(".root");
  canvasname_c.Append(".C");
  c2D->SaveAs(canvasname_pdf);
  c2D->SaveAs(canvasname_eps);
  c2D->SaveAs(canvasname_png);
  c2D->SaveAs(canvasname_root);
  c2D->SaveAs(canvasname_c);

  for (int hh=0; hh<2; hh++){
    for (int it=0; it<2; it++){
      delete hclone[hh][it];
      delete htemp_1D[hh][it];
    }
  }
  c2D->Close();
}
