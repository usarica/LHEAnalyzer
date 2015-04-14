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


void computeOSOwnWeights_CR(int erg_tev=8, int iCR=0, double m4l_low = 100, double m4l_high = 3000);

void combineOSOwnWeights(int erg_tev=8, double m4l_low=100, double m4l_high=3000);

TString produceCRname(int iCR, int iwgt, int isocut, int iM1cut, int catZ1, int catZ2, int icut, char* cappend);

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

  float ZXfake_weight
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

  float strWgt[2] ={ 1, ZXfake_weight };
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

  double result = (strcut ? strWgt[iwgt] : 0);
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
  if (icut>=nSIPCuts) return 0;
  if (icut<0) return 0;

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
  if (histoindex_1>=0 && histoindex_2>=0 && histoindex_3>=0 ){
    TH1F* hOS = (TH1F*)fOS->Get(getstr_hOS(histoindex_1, histoindex_2, histoindex_3, icut));
    OSSFweight = hOS->GetBinContent(hOS->GetXaxis()->FindBin(LepPt));
    OSweighterror = hOS->GetBinError(hOS->GetXaxis()->FindBin(LepPt));
    delete hOS;
  }
  result = OSSFweight / (1.-OSSFweight);

  OSweighterror /= pow(1.-OSSFweight, 2);
  return result;
}

double getQQZZEWKCorrection(float GenHMass, int systematics = 0){
  double A = -0.0782143;
  double MA = 125.369;
  double EWKcorr = A*log(GenHMass / MA);
  double m4lThreshold = 2.0*(91.1876 - 2.4952*2.0);

  double B = 1.89547;
  double MB = 68.1253;
  double OB = 1.53646;
  double QCDcorr = OB - B*exp(-GenHMass / MB);

  if (EWKcorr < -1) EWKcorr = -1;
  if (QCDcorr < 0) QCDcorr = 0;
  if (EWKcorr>1.0){
    EWKcorr = 1.0;
    QCDcorr = 1.0;
  }

  double result=1;
  if (systematics == 0) result = (1.0 + EWKcorr);
  else if (systematics == 1) result = (1.0 + EWKcorr*(2.0 - QCDcorr));
  else if (systematics == -1) result = (1.0 + EWKcorr*QCDcorr);
  else result = 1;
  if (GenHMass<m4lThreshold) result=1;
  return result;
}


void computeOSOwnWeights(int erg_tev=8, double m4l_low = 100, double m4l_high = 3000){
  for (int iCR=0; iCR<=2; iCR++) computeOSOwnWeights_CR(erg_tev, iCR, m4l_low, m4l_high);
  combineOSOwnWeights(erg_tev, m4l_low, m4l_high);
}

void computeOSOwnWeights_CR(int erg_tev, int iCR, double m4l_low, double m4l_high){
//  if (iCR>2 || iCR==1) return;
  char TREE_NAME[]="SelectedTree";

  int EnergyIndex=0; if (erg_tev==8)EnergyIndex=1;
  double my_luminosity = luminosity[EnergyIndex];

  TString erg_dir;
  erg_dir.Form("LHC_%iTeV", erg_tev);
  TString comstring;
  comstring.Form("%iTeV", erg_tev);

  TChain* tc = new TChain(TREE_NAME);

  TString cinput_common_noSIP = user_dir_hep + "No_SIP/" + erg_dir + "/CR/";
  if (iCR<2){
    int sampleindex = kAllSamples-1 + iCR;
    TString cinput = cinput_common_noSIP;
    cinput = cinput + sample_FullSim[sampleindex] + ".root";
    tc->Add(cinput);
  }
  else{
    for (int smp=kGGSamples; smp<kQQBZZSamples; smp++){
      TString cinput = cinput_common_noSIP;
      cinput = cinput + sample_FullSim[smp] + "_CRZLLTree.root";
      tc->Add(cinput);
    }
  }

  if (tc->GetEntries() == 0){
    cout << "Could not find any files, aborting..." << endl; return;
  }
  else if (iCR==2) cout << "qqZZ CR: " << tc->GetEntries() << endl;

  TString coutput_common = user_dir_hep + "Analysis/Auxiliary/";

  TFile* frecord;
  TString fname = "OSOwnfakeweights_";
  if (iCR==0) fname.Append("ZLL_");
  else if (iCR==1) fname.Append("ZL_");
  else if (iCR==2) fname.Append("QQBZZ_");
  fname += comstring;
  fname.Append(Form("_m4l_%.1f_%.1f%s", m4l_low, m4l_high, ".root"));
  fname.Prepend(coutput_common);
  frecord = new TFile(fname, "recreate");

  TFile* foutput[3];

  TString fOSfakeinput = "CRdistributions_";
  fOSfakeinput.Append("ZL_");
  fOSfakeinput += comstring;
  fOSfakeinput.Append("_m4l_105.6_140.6.root");
  fOSfakeinput.Prepend(coutput_common);
  TFile* fweight = new TFile(fOSfakeinput,"read");

  frecord->cd();

  TString hname = "htemp";

  float GenHMass;
  float ZXfake_weight,MC_weight,MC_weight_QQBZZEWK;
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
  if (iCR==2){
    tc->SetBranchAddress("MC_weight", &MC_weight);
    tc->SetBranchAddress("GenHMass", &GenHMass);
  }

  for (int iwgt=0; iwgt<2; iwgt++){
    if (iCR==1 && iwgt==1) continue;
    for (int isocut=0; isocut<4; isocut++){
      if (iCR==1 && isocut>=2) continue;
      if (isocut<2) continue;

      for (int iM1cut=0; iM1cut<1; iM1cut++){
        if (iM1cut==2 && iCR!=1) continue;
        if (iwgt==0) cout << "Un-weighted";
        else cout << "Weighted";
        if (isocut>=1) cout << ", " << strLepIsoCut_label[isocut-1];
        cout << ", " << Z1masswidth_label[iM1cut] << endl;
        cout << "\t"; for (int icut=0; icut<nSIPCuts; icut++) cout << cutLabel[icut] << '\t' << cutLabel[icut] << " Error\t";
        cout << endl;
        for (int catZ1=0; catZ1<2; catZ1++){
          for (int catZ2=0; catZ2<(iCR!=1 ? 3 : 2); catZ2++){
//          for (int catZ2=0; catZ2<(iCR!=1 ? 6 : 2); catZ2++){
            cout << strZ1Category_label[catZ1] << " ";
            if (iCR!=1) cout << strZ2Category_label[catZ2];
            else cout << strLep3Category_label[catZ2];
            for (int icut=0; icut<nSIPCuts; icut++){
              gStyle->SetTitleFont(62, "t");
              gROOT->SetStyle(gStyle->GetName());
              gROOT->ForceStyle();

              TH1F* hsip = new TH1F(hname, "", 1, 0, 1000);
              hsip->Sumw2();

              double ptbins[9]={ 0, 5, 7, 10, 20, 30, 40, 50, 80 };
              double etabins[3]={ 0, 1.2, 2.5 };
              if ((iCR==1 && catZ2==1) || (iCR!=1 && (catZ2==2 || catZ2==3))) etabins[1]=1.45;

//              const int nbinsTxy = 60;
//              double Txybins[nbinsTxy+1];
//              for (int bin=0; bin<=nbinsTxy; bin++)Txybins[bin] = -1400.+2800.*bin/nbinsTxy;
              const int nbinsTxy = 7;
              double Txybins[nbinsTxy+1]={-2000,-500,-200,-50,50,200,500,2000};
/*              const int nbinsmZZ = 4;
              double mZZbins[nbinsmZZ+1]={ 100, 140.6, 170, 250, m4l_high };
              const int nbinsmZ1 = 10;
              double mZ1bins[nbinsmZ1+1];
              for (int bin=0; bin<=nbinsmZ1; bin++)mZ1bins[bin] = 0.+120.*bin/nbinsmZ1;
              const int nbinsmZ2 = 6;
              double mZ2bins[nbinsmZ2+1];
              for (int bin=0; bin<=nbinsmZ2; bin++)mZ2bins[bin] = 12.+132.*bin/nbinsmZ2;
*/
              const int nbinsmZZ = 4;
              double mZZbins[nbinsmZZ+1]={ 0, 100, 140, 170, 3000 };
              const int nbinsmZ1 = 5;
              double mZ1bins[nbinsmZ1+1]={ 0, 40, 62.5, 85, 97.5, 120 };
              const int nbinsmZ2 = 3;
              double mZ2bins[nbinsmZ2+1]={ 0, 35, 65, 150 };

              TH2F* hmZZTxyOS;
              TH2F* hmZZTxyOS_extra;
              TH3F* hmZZmZ1mZ2OS;
              TH3F* hmZZmZ1mZ2OS_extra;
              if (iCR!=1 && isocut>1 && catZ2<3){
                hmZZTxyOS = new TH2F(produceCRname(iCR, iwgt, isocut, iM1cut, catZ1, catZ2, icut, "mZZ_Txy_OS"), "", nbinsmZZ, mZZbins, nbinsTxy, Txybins);
                hmZZTxyOS->Sumw2();
                hmZZTxyOS->SetXTitle("m_{4l} (GeV)");
                hmZZTxyOS->SetYTitle("T_{xy} (#mum)");
                hmZZTxyOS->SetZTitle("Events / bin");
                hmZZTxyOS->SetOption("colz");

                hmZZmZ1mZ2OS = new TH3F(produceCRname(iCR, iwgt, isocut, iM1cut, catZ1, catZ2, icut, "mZZ_mZ1mZ2_OS"), "", nbinsmZZ, mZZbins, nbinsmZ1, mZ1bins, nbinsmZ2, mZ2bins);
                hmZZmZ1mZ2OS->Sumw2();
                hmZZmZ1mZ2OS->SetXTitle("m_{4l} (GeV)");
                hmZZmZ1mZ2OS->SetYTitle("m_{1} (GeV)");
                hmZZmZ1mZ2OS->SetZTitle("m_{2} (GeV)");
                hmZZmZ1mZ2OS->SetOption("box");
                if (isocut==2 && iwgt>0){
                  hmZZTxyOS_extra = new TH2F(produceCRname(iCR, iwgt, isocut, iM1cut, catZ1, catZ2, icut, "mZZ_Txy_OS_2p2fin3p1f"), "", nbinsmZZ, mZZbins, nbinsTxy, Txybins);
                  hmZZTxyOS_extra->Sumw2();
                  hmZZTxyOS_extra->SetXTitle("m_{4l} (GeV)");
                  hmZZTxyOS_extra->SetYTitle("T_{xy} (#mum)");
                  hmZZTxyOS_extra->SetZTitle("Events / bin");
                  hmZZTxyOS_extra->SetOption("colz");

                  hmZZmZ1mZ2OS_extra = new TH3F(produceCRname(iCR, iwgt, isocut, iM1cut, catZ1, catZ2, icut, "mZZ_mZ1mZ2_OS_2p2fin3p1f"), "", nbinsmZZ, mZZbins, nbinsmZ1, mZ1bins, nbinsmZ2, mZ2bins);
                  hmZZmZ1mZ2OS_extra->Sumw2();
                  hmZZmZ1mZ2OS_extra->SetXTitle("m_{4l} (GeV)");
                  hmZZmZ1mZ2OS_extra->SetYTitle("m_{1} (GeV)");
                  hmZZmZ1mZ2OS_extra->SetZTitle("m_{2} (GeV)");
                  hmZZmZ1mZ2OS_extra->SetOption("box");
                }
              }


              for (int ev=0; ev<tc->GetEntries(); ev++){
                MC_weight=1;
                MC_weight_QQBZZEWK=1;
                
                tc->GetEntry(ev);

                if (iCR==2){
                  MC_weight*=my_luminosity;
                  MC_weight_QQBZZEWK = getQQZZEWKCorrection(GenHMass);
                  MC_weight*=MC_weight_QQBZZEWK;
                  ZXfake_weight = ZXfake_weight*MC_weight;
                }

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

                  1
//                  ZXfake_weight
                  );
                if (wgt==0) continue;
                float wgterror=0;


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

                float OSSFweight[2] ={ 0 };
                float OSSFweighterror[2] ={ 0 };
                if (iCR!=1 && isocut>1 && catZ2<3 && wgt>0){

                  OSSFweight[0] = applyOSSFweight(
                    fweight,

                    CRflag,
                    icut,
                    Z1ids, Z2ids,
                    Lep3combRelIsoPF,
                    Lep3isID,

                    Lep3ID,
                    Lep3Pt,
                    Lep3Eta,

                    OSSFweighterror[0]
                    );

                  OSSFweight[1] = applyOSSFweight(
                    fweight,

                    CRflag,
                    icut,
                    Z1ids, Z2ids,
                    Lep4combRelIsoPF,
                    Lep4isID,

                    Lep4ID,
                    Lep4Pt,
                    Lep4Eta,

                    OSSFweighterror[1]
                    );

                  if (isocut==2 && (OSSFweight[0]==0 || OSSFweight[1]==0)) cout << "WARNING: 2p2f OSSF weights are not all positive!" << endl;
                  if (isocut==3 && (OSSFweight[0]==0 && OSSFweight[1]==0)) cout << "WARNING: 3p1f OSSF weights are all zero!" << endl;
                  if (isocut==3 && (OSSFweight[0]!=0 && OSSFweight[1]!=0)) cout << "WARNING: 3p1f OSSF weights are all non-zero!" << endl;

                  float totalOS = 1;
                  float sumOS = 0;

                  float totalOSerror=0;
                  float sumOSerror=0;

                  if (iwgt>0){
                    if (OSSFweight[0]>0){
                      totalOS *= OSSFweight[0];
                      totalOSerror += pow(OSSFweighterror[0]/OSSFweight[0], 2);
                    }
                    if (OSSFweight[1]>0){
                      totalOS *= OSSFweight[1];
                      totalOSerror += pow(OSSFweighterror[1]/OSSFweight[1], 2);
                    }
                    totalOSerror = sqrt(totalOSerror) * totalOS;
                    if (isocut==2){
                      sumOS = OSSFweight[0] + OSSFweight[1];
                      sumOSerror = sqrt(pow(OSSFweighterror[0], 2) + pow(OSSFweighterror[1], 2));
                    }
                    if (iCR==2){
                      totalOS*=MC_weight;
                      sumOS*=MC_weight;
                      totalOSerror*=MC_weight;
                      sumOSerror*=MC_weight;
                    }
                    wgt = totalOS;
                    wgterror = totalOSerror;
                  }
                  else{
                    if (iCR==2){
                      totalOS*=MC_weight;
                      sumOS*=MC_weight;
                    }
                  }

                  hmZZTxyOS->Fill(strdraw_mZZ, strdraw_Txy, totalOS);
                  hmZZmZ1mZ2OS->Fill(strdraw_mZZ, strdraw_mZ1, strdraw_mZ2, totalOS);

                  if (iwgt>0){
/*
                    int binx = hmZZTxyOS->GetXaxis()->FindBin(strdraw_mZZ);
                    int biny = hmZZTxyOS->GetYaxis()->FindBin(strdraw_Txy);
                    hmZZTxyOS->SetBinError(binx, biny, sqrt(pow(hmZZTxyOS->GetBinError(binx, biny), 2)+pow(totalOSerror, 2)));

                    binx = hmZZmZ1mZ2OS->GetXaxis()->FindBin(strdraw_mZZ);
                    biny = hmZZmZ1mZ2OS->GetYaxis()->FindBin(strdraw_mZ1);
                    int binz = hmZZmZ1mZ2OS->GetZaxis()->FindBin(strdraw_mZ2);
                    hmZZmZ1mZ2OS->SetBinError(binx, biny, binz, sqrt(pow(hmZZmZ1mZ2OS->GetBinError(binx, biny, binz), 2)+pow(totalOSerror, 2)));
*/
                    if (isocut==2){
                      hmZZTxyOS_extra->Fill(strdraw_mZZ, strdraw_Txy, sumOS);
                      hmZZmZ1mZ2OS_extra->Fill(strdraw_mZZ, strdraw_mZ1, strdraw_mZ2, sumOS);
/*
                      int binx = hmZZTxyOS_extra->GetXaxis()->FindBin(strdraw_mZZ);
                      int biny = hmZZTxyOS_extra->GetYaxis()->FindBin(strdraw_Txy);
                      hmZZTxyOS_extra->SetBinError(binx, biny, sqrt(pow(hmZZTxyOS_extra->GetBinError(binx, biny), 2)+pow(sumOSerror, 2)));

                      binx = hmZZmZ1mZ2OS_extra->GetXaxis()->FindBin(strdraw_mZZ);
                      biny = hmZZmZ1mZ2OS_extra->GetYaxis()->FindBin(strdraw_mZ1);
                      int binz = hmZZmZ1mZ2OS_extra->GetZaxis()->FindBin(strdraw_mZ2);
                      hmZZmZ1mZ2OS_extra->SetBinError(binx, biny, binz, sqrt(pow(hmZZmZ1mZ2OS_extra->GetBinError(binx, biny, binz), 2)+pow(sumOSerror, 2)));
*/
                    }
                  }

                }

                hsip->Fill(strdraw, wgt);
                if (iwgt>0){
                  int binx = hsip->GetXaxis()->FindBin(strdraw);
                  hsip->SetBinError(binx, sqrt(pow(hsip->GetBinError(binx), 2)+pow(wgterror, 2)));
                }

              }

              if (iCR!=1 && isocut>1 && catZ2<3){
                frecord->WriteTObject(hmZZTxyOS);
                delete hmZZTxyOS;
                frecord->WriteTObject(hmZZmZ1mZ2OS);
                delete hmZZmZ1mZ2OS;
                if (isocut==2 && iwgt>0){
                  frecord->WriteTObject(hmZZTxyOS_extra);
                  delete hmZZTxyOS_extra;
                  frecord->WriteTObject(hmZZmZ1mZ2OS_extra);
                  delete hmZZmZ1mZ2OS_extra;
                }
              }

              double totalerror = 0;
              double totalintegral = hsip->IntegralAndError(1,hsip->GetNbinsX(),totalerror);

              cout << '\t' << totalintegral << '\t' << totalerror;

              delete hsip;
            }
            cout << endl;
          }
        }
        cout << endl;
      }
      cout << endl;
    }
    cout << endl;
  }
  delete tc;

  fweight->Close();
  frecord->Close();
}

void combineOSOwnWeights(int erg_tev, double m4l_low, double m4l_high){
  TString erg_dir;
  erg_dir.Form("LHC_%iTeV", erg_tev);
  TString comstring;
  comstring.Form("%iTeV", erg_tev);

  TString cinput_common = user_dir_hep + "Analysis/Auxiliary/";
  TString coutput_common = user_dir_hep + "Analysis/Auxiliary/";

  TFile* finput[2];
  for (int iCR=0; iCR<2; iCR++){
    TString fname = "OSOwnfakeweights_";
    if (iCR==0) fname.Append("ZLL_");
    else fname.Append("QQBZZ_");
    fname += comstring;
    fname.Append(Form("_m4l_%.1f_%.1f%s", m4l_low, m4l_high, ".root"));
    fname.Prepend(cinput_common);
    finput[iCR] = new TFile(fname, "read");
  }

  TString coutput = "OSOwnfakeweights_AfterCorrection_";
  coutput += comstring;
  coutput.Prepend(coutput_common);
  coutput.Append(".root");
  TFile* foutput = new TFile(coutput, "recreate");

  TH3F* hOSweight_ZLL_unweighted[nSIPCuts][3*2];
  TH3F* hOSweight_ZLL_2p2fin3p1f[nSIPCuts][3*2];
  TH3F* hOSweight_ZLL_weighted[nSIPCuts][3*2];
  TH3F* hOSweight_qqZZ[nSIPCuts][3*2];
  TH3F* hOSweight_3p1fFinal[nSIPCuts][3*2];

  int iM1cut=0;
  int isocut=3;
  for (int icut=0; icut<nSIPCuts; icut++){
    for (int catZ1=0; catZ1<2; catZ1++){
      for (int catZ2=0; catZ2<3; catZ2++){
        int index_catZ = 3*catZ1+catZ2;
        TString strOSweight_ZLL_unweighted = produceCRname(0, 0, isocut, iM1cut, catZ1, catZ2, icut, "mZZ_mZ1mZ2_OS");
        TString strOSweight_ZLL_2p2fin3p1f = produceCRname(0, 1, 2, iM1cut, catZ1, catZ2, icut, "mZZ_mZ1mZ2_OS_2p2fin3p1f");
        TString strOSweight_ZLL_weighted = produceCRname(0, 1, isocut, iM1cut, catZ1, catZ2, icut, "mZZ_mZ1mZ2_OS");
        TString strOSweight_qqZZ = produceCRname(2, 0, isocut, iM1cut, catZ1, catZ2, icut, "mZZ_mZ1mZ2_OS");

        cout << strOSweight_ZLL_unweighted << endl;
        hOSweight_ZLL_unweighted[icut][index_catZ] = (TH3F*)finput[0]->Get(strOSweight_ZLL_unweighted);
        cout << strOSweight_ZLL_2p2fin3p1f << endl;
        hOSweight_ZLL_2p2fin3p1f[icut][index_catZ] = (TH3F*)finput[0]->Get(strOSweight_ZLL_2p2fin3p1f);
        cout << strOSweight_ZLL_weighted << endl;
        hOSweight_ZLL_weighted[icut][index_catZ] = (TH3F*)finput[0]->Get(strOSweight_ZLL_weighted);
        cout << strOSweight_qqZZ << endl;
        hOSweight_qqZZ[icut][index_catZ] = (TH3F*)finput[1]->Get(strOSweight_qqZZ);
/*
        cout << "Nbins hOSweight_ZLL_unweighted: " << hOSweight_ZLL_unweighted[icut][index_catZ]->GetNbinsX()*hOSweight_ZLL_unweighted[icut][index_catZ]->GetNbinsY()*hOSweight_ZLL_unweighted[icut][index_catZ]->GetNbinsZ() << endl;
        cout << "Nbins hOSweight_ZLL_2p2fin3p1f: " << hOSweight_ZLL_2p2fin3p1f[icut][index_catZ]->GetNbinsX()*hOSweight_ZLL_2p2fin3p1f[icut][index_catZ]->GetNbinsY()*hOSweight_ZLL_2p2fin3p1f[icut][index_catZ]->GetNbinsZ() << endl;
        cout << "Nbins hOSweight_ZLL_weighted: " << hOSweight_ZLL_weighted[icut][index_catZ]->GetNbinsX()*hOSweight_ZLL_weighted[icut][index_catZ]->GetNbinsY()*hOSweight_ZLL_weighted[icut][index_catZ]->GetNbinsZ() << endl;
        cout << "Nbins hOSweight_qqZZ: " << hOSweight_qqZZ[icut][index_catZ]->GetNbinsX()*hOSweight_qqZZ[icut][index_catZ]->GetNbinsY()*hOSweight_qqZZ[icut][index_catZ]->GetNbinsZ() << endl;
*/
        TString strOSweight_3p1fFinal = "OSfakeweight_";
        if (catZ1==0) strOSweight_3p1fFinal.Append("_Zmumu");
        else strOSweight_3p1fFinal.Append("_Zee");
        if (catZ2==0) strOSweight_3p1fFinal.Append("_mumuOS");
        else if (catZ2==1) strOSweight_3p1fFinal.Append("_emuOS");
        else if (catZ2==2) strOSweight_3p1fFinal.Append("_eeOS");
        strOSweight_3p1fFinal.Append("_");
        strOSweight_3p1fFinal += cutNames[icut];
        strOSweight_3p1fFinal.Append("_");
        strOSweight_3p1fFinal += strLepIsoCut_label[isocut-1];

        hOSweight_3p1fFinal[icut][index_catZ] = (TH3F*)hOSweight_ZLL_weighted[icut][index_catZ]->Clone(strOSweight_3p1fFinal);
        hOSweight_3p1fFinal[icut][index_catZ]->Sumw2();
        for (int binz=0; binz<=hOSweight_3p1fFinal[icut][index_catZ]->GetNbinsZ()+1; binz++){
          for (int biny=0; biny<=hOSweight_3p1fFinal[icut][index_catZ]->GetNbinsY()+1; biny++){
            for (int binx=0; binx<=hOSweight_3p1fFinal[icut][index_catZ]->GetNbinsX()+1; binx++){
              hOSweight_3p1fFinal[icut][index_catZ]->SetBinContent(binx, biny, binz, 1);
              hOSweight_3p1fFinal[icut][index_catZ]->SetBinError(binx, biny, binz, 0);
            }
          }
        }

        for (int binz=0; binz<=hOSweight_3p1fFinal[icut][index_catZ]->GetNbinsZ()+1; binz++){
          for (int biny=0; biny<=hOSweight_3p1fFinal[icut][index_catZ]->GetNbinsY()+1; biny++){
            for (int binx=0; binx<=hOSweight_3p1fFinal[icut][index_catZ]->GetNbinsX()+1; binx++){
              hOSweight_ZLL_unweighted[icut][index_catZ]->SetBinError(binx, biny, binz, 0);
              hOSweight_ZLL_2p2fin3p1f[icut][index_catZ]->SetBinError(binx, biny, binz, 0);
            }
          }
        }

        hOSweight_qqZZ[icut][index_catZ]->Divide(hOSweight_ZLL_unweighted[icut][index_catZ]);
        hOSweight_ZLL_2p2fin3p1f[icut][index_catZ]->Divide(hOSweight_ZLL_unweighted[icut][index_catZ]);
        hOSweight_3p1fFinal[icut][index_catZ]->Add(hOSweight_qqZZ[icut][index_catZ], -1);
        hOSweight_3p1fFinal[icut][index_catZ]->Add(hOSweight_ZLL_2p2fin3p1f[icut][index_catZ], -1);
        int nzeros=0;
        for (int binz=0; binz<=hOSweight_3p1fFinal[icut][index_catZ]->GetNbinsZ()+1; binz++){
          for (int biny=0; biny<=hOSweight_3p1fFinal[icut][index_catZ]->GetNbinsY()+1; biny++){
            for (int binx=0; binx<=hOSweight_3p1fFinal[icut][index_catZ]->GetNbinsX()+1; binx++){
              double bincontent = hOSweight_3p1fFinal[icut][index_catZ]->GetBinContent(binx, biny, binz);
              double binerror = hOSweight_3p1fFinal[icut][index_catZ]->GetBinError(binx, biny, binz);

              if (binerror>=fabs(bincontent)){
                cout << "Bin error larger than content: " << strOSweight_3p1fFinal << '\t' << bincontent << '\t' << binerror
                  << '\t' << hOSweight_3p1fFinal[icut][index_catZ]->GetXaxis()->GetBinCenter(binx)
                  << '\t' << hOSweight_3p1fFinal[icut][index_catZ]->GetYaxis()->GetBinCenter(biny)
                  << '\t' << hOSweight_3p1fFinal[icut][index_catZ]->GetZaxis()->GetBinCenter(binz)
                  << endl;
              }

              if (bincontent<0){
                cout << "Bin content <0: " << strOSweight_3p1fFinal << '\t' << bincontent
                  << '\t' << hOSweight_3p1fFinal[icut][index_catZ]->GetXaxis()->GetBinCenter(binx)
                  << '\t' << hOSweight_3p1fFinal[icut][index_catZ]->GetYaxis()->GetBinCenter(biny)
                  << '\t' << hOSweight_3p1fFinal[icut][index_catZ]->GetZaxis()->GetBinCenter(binz)
                  << endl;
                hOSweight_3p1fFinal[icut][index_catZ]->SetBinContent(binx, biny, binz,0);
                hOSweight_3p1fFinal[icut][index_catZ]->SetBinError(binx, biny, binz, (fabs(bincontent)>binerror?bincontent:binerror));
              }
              if (bincontent<1 && bincontent>0)nzeros++;
            }
          }
        }

        cout << "Nvalid: " << nzeros << '/' << hOSweight_3p1fFinal[icut][index_catZ]->GetNbinsX()*hOSweight_3p1fFinal[icut][index_catZ]->GetNbinsY()*hOSweight_3p1fFinal[icut][index_catZ]->GetNbinsZ() << endl;

        hOSweight_ZLL_weighted[icut][index_catZ]->Multiply(hOSweight_3p1fFinal[icut][index_catZ]);
        hOSweight_ZLL_weighted[icut][index_catZ]->SetName(Form("%s_corrected", hOSweight_ZLL_weighted[icut][index_catZ]->GetName()));
        foutput->WriteTObject(hOSweight_3p1fFinal[icut][index_catZ]);
        foutput->WriteTObject(hOSweight_ZLL_weighted[icut][index_catZ]);
      }
    }
  }

  TString fOSfakeinput = "CRdistributions_";
  fOSfakeinput.Append("ZL_");
  fOSfakeinput += comstring;
  fOSfakeinput.Append("_m4l_105.6_140.6.root");
  fOSfakeinput.Prepend(coutput_common);
  TFile* fweight = new TFile(fOSfakeinput, "read");
  for (int iZ=0; iZ<2; iZ++){
    for (int iL=0; iL<2; iL++){
      for (int iEB=0; iEB<2; iEB++){
        for (int icut=0; icut<nSIPCuts; icut++){
          TH1F* hOS = (TH1F*)fweight->Get(getstr_hOS(iZ, iL, iEB, icut));
          foutput->WriteTObject(hOS);
          delete hOS;
        }
      }
    }
  }

  fweight->Close();
  foutput->Close();
  for (int iCR=0; iCR<2; iCR++) finput[iCR]->Close();
}

TString produceCRname(int iCR, int iwgt, int isocut, int iM1cut, int catZ1, int catZ2, int icut, char* cappend){
  TString name;
  if (iCR==0) name = "ZLL";
  else if (iCR==1) name = "ZL";
  else name = "QQBZZ";

  if (catZ1==0) name.Append("_Zmumu");
  else name.Append("_Zee");
  if (iCR!=1){
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
