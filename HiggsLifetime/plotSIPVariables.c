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

string processName[8] ={
  "CTau0", "CTau100", "CTau500", "CTau1000",
  "qqZZ", "ggZZ", "CR", "data"
};
TString processTitle[8] ={
  "Higgs c#tau_{H}=0 #mum", "Higgs c#tau_{H}=100 #mum", "Higgs c#tau_{H}=500 #mum", "Higgs c#tau_{H}=1000 #mum",
  "q#bar{q}#rightarrow4l bkg.", "gg#rightarrow4l bkg.", "Z+X bkg.", "Observed"
};
enum{
  vnOPVTracks,
  vsigmaDirRatio,
  vdiffDotTransverse_DirSubPV,
  vPVSIP,
  vLep1_Z1SIP,
  vLep2_Z1SIP,
  vLep3_Z1SIP,
  vLep4_Z1SIP,
  v4l_chi2,
  vTxy,
  vTxy_BS,
  vDbkg,
  nVariables
};
TString varNames[nVariables] ={
  "nOPVTracks",
  "sigmaDirRatio",
  "diffDotTransverse_DirSubPV",

  "PVSIP",
  "Lep1_Z1SIP",
  "Lep2_Z1SIP",
  "Lep3_Z1SIP",
  "Lep4_Z1SIP",

  "4l_chi2",
  "Txy",
  "Txy_BS",
  "Dbkg"
};
TString varTitles[nVariables] ={
  "N^{PV}_{tracks}",
  "#bar{#sigma}^{2}_{PV} / #bar{#sigma}^{2}_{4l}",
  "(r^{reco}_{#lower[-0.3]{PV}} - r^{#lower[-0.15]{true}}_{PV}) . #hat{p}^{#lower[0.2]{reco}}_{#lower[-0.2]{T}} (#mum)",

  "#it{SIP}_{3D}",
  "#it{SIP}^{Z1}_{1} (Signed)",
  "#it{SIP}^{Z1}_{2} (Signed)",
  "#it{SIP}^{Z1}_{3} (Signed)",
  "#it{SIP}^{Z1}_{4} (Signed)",

  "#chi^{2}_{4l}",
  "T_{xy} (#mum)",
  "c#Deltat (#mum)",
  "D_{bkg}"
};
float KD[nVariables] ={ 0 };
float KD_ranges[nVariables][2] ={
  { 0, 100 },
  { 0, 6 },
  { -1000, 1000 },

  { 0, 30 },
  { -4, 4 },
  { -4, 4 },
  { -5, 5 },
  { -5, 5 },

  { 0, 40 },

  { -1000, 1000 },
  { -1000, 1000 },
  { 0, 1 }

};
int KD_nbins[nVariables] ={
  25,
  60,
  80,

  30,
  30,
  30,
  40,
  40,

  40,

  50,
  50,
  20
};
const int nVariables_2D=13;
int KD_pairing[nVariables_2D][2] ={
  { vsigmaDirRatio, vnOPVTracks },
  { vdiffDotTransverse_DirSubPV, vnOPVTracks },
  { vLep1_Z1SIP, vPVSIP },
  { vLep2_Z1SIP, vPVSIP },
  { vLep3_Z1SIP, vPVSIP },
  { vLep4_Z1SIP, vPVSIP },
  { vLep1_Z1SIP, v4l_chi2 },
  { vLep2_Z1SIP, v4l_chi2 },
  { vLep3_Z1SIP, v4l_chi2 },
  { vLep4_Z1SIP, v4l_chi2 },
  { v4l_chi2, vPVSIP },

  { vLep2_Z1SIP, vLep1_Z1SIP },
  { vLep4_Z1SIP, vLep3_Z1SIP }
};

const int nMZZ_ranges=4;
float systZZMass_range[nMZZ_ranges][2] ={
  { 105.6, 140.6 }, { 70, 105.6 }, { 140.6, 170 }, { 170, 800 }
};

const int nSIPCuts = 2; // <4 old SIP vs new SIP schemes
const int minSIPCuts = 1; // 0: <4 old SIP, 1: new SIP schemes
const int newSIPCuts = 5; // 0: No SIP; 1: Z1SIP, 2: KalmanChi2, 3: Z1SIP + KalmanChi2, 4: Combined new+old
TString cutNames[newSIPCuts + 1] ={
  "PVSIPcut",
  "NoSIPcut",
  "Z1SIPcut",
  "4lchi2cut",
  "Z1SIP_4lchi2cut",
  "Combinedcut"
};
TString cutLabel[newSIPCuts + 1] ={
  "Old cut",
  "No cut",
  "Z1-SIP",
  "chi**2",
  "New cut",
  "New + old"
};
TString ctitle_SIPCut[newSIPCuts + 1] ={
  "#it{SIP}_{3D}<4",
  "No vertex cut",
  "#it{SIP}^{Z1}_{1,2}<4, #it{SIP}^{Z1}_{3,4}<5",
  "#chi^{2}_{4l}<30",
  "#it{SIP}^{Z1}_{1,2}<4, #it{SIP}^{Z1}_{3,4}<5, #chi^{2}_{4l}<30",
  "Combined cut"
};


void generic_Histo2DPlotter(TFile* foutput, TString canvasDir, TH2* hSlice, int erg_tev, float lumi, int channel, bool isCRonly=false);
void symmetrize_PromptTemplates(TH1F* hrepair);


void create_SIPCut_YieldRatios_single(int iProcess, float m4l_low=105.6, float m4l_high=140.6){
  if (iProcess==6){ cout << "ZX estimation is separated." << endl; return; }

  char TREE_NAME[]="SelectedTree";
  int processMap[7] ={
    0, 1, 2, 3,
    kGGSamples, // qqZZ full range samples; 6 of them
    kGGOLDSamples, // MCFM gg, 3 of them
    (kAllSamples-1)
  };
  TString strnewSIP[newSIPCuts] ={
    Form("MC_weight*( ZZMass>=%.1f && ZZMass<%.1f )", m4l_low, m4l_high),
    Form("MC_weight*( ZZMass>=%.1f && ZZMass<%.1f && abs(Lep1_Z1SIP)<4 && abs(Lep2_Z1SIP)<4 && abs(Lep3_Z1SIP)<5 && abs(Lep4_Z1SIP)<5 )", m4l_low, m4l_high),
    Form("MC_weight*( ZZMass>=%.1f && ZZMass<%.1f && KalmanCandVtx_chi2<30 )", m4l_low, m4l_high),
    Form("MC_weight*( ZZMass>=%.1f && ZZMass<%.1f && abs(Lep1_Z1SIP)<4 && abs(Lep2_Z1SIP)<4 && abs(Lep3_Z1SIP)<5 && abs(Lep4_Z1SIP)<5 && KalmanCandVtx_chi2<30 )", m4l_low, m4l_high),
    Form("MC_weight*( ZZMass>=%.1f && ZZMass<%.1f && abs(Lep1_Z1SIP)<4 && abs(Lep2_Z1SIP)<4 && abs(Lep3_Z1SIP)<5 && abs(Lep4_Z1SIP)<5 && KalmanCandVtx_chi2<30 && max(max(max(Lep3SIP,Lep4SIP),Lep2SIP),Lep1SIP)<4 )", m4l_low, m4l_high)
  };
  TString strnewSIP_CR[newSIPCuts] ={
    Form("ZXfake_weight*( ZZMass>=%.1f && ZZMass<%.1f )", m4l_low, m4l_high),
    Form("ZXfake_weight*( ZZMass>=%.1f && ZZMass<%.1f && abs(Lep1_Z1SIP)<4 && abs(Lep2_Z1SIP)<4 && abs(Lep3_Z1SIP)<5 && abs(Lep4_Z1SIP)<5 )", m4l_low, m4l_high),
    Form("ZXfake_weight*( ZZMass>=%.1f && ZZMass<%.1f && KalmanCandVtx_chi2<30 )", m4l_low, m4l_high),
    Form("ZXfake_weight*( ZZMass>=%.1f && ZZMass<%.1f && abs(Lep1_Z1SIP)<4 && abs(Lep2_Z1SIP)<4 && abs(Lep3_Z1SIP)<5 && abs(Lep4_Z1SIP)<5 && KalmanCandVtx_chi2<30 )", m4l_low, m4l_high),
    Form("ZXfake_weight*( ZZMass>=%.1f && ZZMass<%.1f && abs(Lep1_Z1SIP)<4 && abs(Lep2_Z1SIP)<4 && abs(Lep3_Z1SIP)<5 && abs(Lep4_Z1SIP)<5 && KalmanCandVtx_chi2<30 && max(max(max(Lep3SIP,Lep4SIP),Lep2SIP),Lep1SIP)<4 )", m4l_low, m4l_high)
  };

  for (int erg_tev = 7; erg_tev < 9; erg_tev++){
    TString erg_dir;
    erg_dir.Form("LHC_%iTeV", erg_tev);
    TString comstring;
    comstring.Form("%iTeV", erg_tev);
    int EnergyIndex = 0; if (erg_tev == 8) EnergyIndex = 1;

    for (int folder = 0; folder < 3; folder++){
      TString OUTPUT_NAME = "LifetimeKD_RelativeSIPYields_";
      OUTPUT_NAME.Append(Form("%s", processName[iProcess].c_str()));
      if (iProcess != 6){
        OUTPUT_NAME.Append(Form("_%s_", user_folder[folder]));
        OUTPUT_NAME.Append(comstring);
      }
      OUTPUT_NAME.Append(Form("_m4l%.1f_%.1f%s", m4l_low, m4l_high, ".root"));
      if (iProcess == 6 && (folder>0||erg_tev==7)) continue;

      TString cinput_common_withSIP = user_dir_hep + erg_dir + "/";
      TString cinput_common_noSIP = user_dir_hep + "No_SIP/" + erg_dir + "/";
      cinput_common_withSIP.Append(Form("%s/", user_folder[folder]));
      cinput_common_noSIP.Append(Form("%s/", user_folder[folder]));
      TString cinput_common = cinput_common_noSIP;
      TString coutput_common = user_dir_hep + "Analysis/Auxiliary/";
      TString coutput = coutput_common + OUTPUT_NAME;
      TFile* foutput = new TFile(coutput, "recreate");

      TChain* tc = new TChain(TREE_NAME);
      TString cinput;
      int smp = processMap[iProcess];
      bool treesExist = true;
      if (iProcess == 0){
        cinput = cinput_common + sample_FullSim[smp] + ".root";
        tc->Add(cinput);
      }
      else if (iProcess < 4 && iProcess>1){
        cinput = cinput_common + sample_FullSim[smp] + ".root";
        tc->Add(cinput);
      }
      else if (iProcess==1){
        cinput = cinput_common + sample_FullSim[smp] + ".root";
        tc->Add(cinput);
      }
      else if (iProcess == 4){
        smp = processMap[iProcess];
        for (int it = 0; it < 6; it++){ // qqZZ has 6 samples, gg has 2 samples
          cout << smp << '\t' << sample_FullSim[smp] << endl;
          cinput = cinput_common + sample_FullSim[smp] + ".root";
          tc->Add(cinput);
          smp++;
        }
      }
      else if (iProcess == 5){
        smp = processMap[iProcess];
        for (int it = 0; it < 3; it++){ // qqZZ has 6 samples, gg has 2 samples
          cout << smp << '\t' << sample_FullSim[smp] << endl;
          cinput = cinput_common + sample_FullSim[smp] + ".root";
          tc->Add(cinput);
          smp++;
        }
      }
      else if (iProcess == 6){
        for (int ee = 7; ee < 9; ee++){
          cinput = user_dir_hep;
          cinput.Append("No_SIP/");
          cinput.Append(Form("LHC_%iTeV/", ee));
          cinput = cinput + "CR/" + sample_FullSim[smp] + ".root";
          tc->Add(cinput);
        }
      }
      if (tc->GetEntries() == 0){ cout << "Could not find the file, aborting..." << endl; treesExist = false; }
      if (!treesExist){
        foutput->Close();
        TString rmCommand = "rm ";
        rmCommand.Append(coutput);
        gSystem->Exec(rmCommand);
        continue;
      }
      TString hname = Form("h%s", processName[iProcess].c_str());
      TH1F* hsip = new TH1F(hname, "Yield ratio to SIP (Old)<4", newSIPCuts + 1, 0, newSIPCuts + 1);
      hsip->Sumw2();
      TString strdraw = "0>>"; strdraw.Append(hname);
      TString strcut = strnewSIP[0];
      strcut = strcut + "*( max(max(max(Lep3SIP,Lep4SIP),Lep2SIP),Lep1SIP)<4 )";
//    if (iProcess>=4) strcut = "MC_weight*(ZZMass>=105.6 && ZZMass<140.6 && max(max(max(Lep3SIP,Lep4SIP),Lep2SIP),Lep1SIP)<4 )";
      if (iProcess==6) strcut = "ZXfake_weight*(ZZMass>=105.6 && ZZMass<140.6 && max(max(max(Lep3SIP,Lep4SIP),Lep2SIP),Lep1SIP)<4 )";
//    cout << "Drawing: " << strdraw << " with cut " << strcut << endl;
      tc->Draw(strdraw, strcut);
      for (int binx = 2; binx <= newSIPCuts + 1; binx++){
        strdraw = Form("%i>>+", binx - 1); strdraw.Append(hname);
        strcut = strnewSIP[binx - 2];
        if (iProcess == 6) strcut = strnewSIP_CR[binx - 2];
        //        cout << "Drawing: " << strdraw << " with cut " << strcut << endl;
        tc->Draw(strdraw, strcut);
      }
      double nosipval = hsip->GetBinContent(2);
      double nosiperr = hsip->GetBinError(2);

      double commonsipval = hsip->GetBinContent(6);
      double commonsiperr = hsip->GetBinError(6);
      double newsipval = hsip->GetBinContent(5);
      double newsiperr = hsip->GetBinError(5);
      double oldsipval = hsip->GetBinContent(1);
      double oldsiperr = hsip->GetBinError(1);
      double oldnotnewval = oldsipval-commonsipval;
      double oldnotnewerr = pow(oldsiperr, 2)-pow(commonsiperr, 2);
      double newnotoldval = newsipval-commonsipval;
      double newnotolderr = pow(newsiperr, 2)-pow(commonsiperr, 2);
      double term1 = newnotoldval/oldsipval;
      double err_term1 = 0;
      if (newnotoldval!=0) err_term1 += pow(newnotolderr/newnotoldval, 2);
      if (oldsipval!=0) err_term1 += pow(oldsiperr/oldsipval, 2);
      err_term1 = sqrt(err_term1)*term1;
      double term2 = oldnotnewval/commonsipval;
      double err_term2 = 0;
      if (oldnotnewval!=0) err_term2 += pow(oldnotnewerr/oldnotnewval, 2);
      if (commonsipval!=0) err_term2 += pow(commonsiperr/commonsipval, 2);
      err_term2 = sqrt(err_term2)*term2;
      term2 = 1./(1.+term2);
      err_term2 = err_term2*pow(term2, 2);
      double ratio = term1+term2;
      double ratio_err = sqrt(pow(err_term2, 2)+pow(err_term1, 2));
      cout << processName[iProcess] << '\t' << erg_tev << " TeV\t" << user_folder[folder] << ": " << ratio << " +- " << ratio_err << endl;

      for (int binx=1; binx<=hsip->GetNbinsX(); binx++){
        double val = hsip->GetBinContent(binx);
        double err = hsip->GetBinError(binx);

        double new_err = sqrt(pow(err*(nosipval-val), 2) + (pow(nosiperr, 2)-pow(err, 2))*pow(val, 2))/pow(nosipval, 2);
        double new_val = val/nosipval;
        hsip->SetBinContent(binx, new_val);
        hsip->SetBinError(binx, new_err);
      }
      foutput->WriteTObject(hsip);


      delete tc;
      delete hsip;
      foutput->Close();
    }
  }
}

void create_SIPCut_YieldRatios(){
  float yield_range[nMZZ_ranges][2] ={
    { 105.6, 140.6 }, { 70., 105.6 }, { 140.6, 170. }, { 170., 3000. }
  };
  for (int p=0; p<6; p++){
    for (int rr=0; rr<nMZZ_ranges; rr++){
      if (rr==1 && (p==6 || p==5)) continue;
      if (rr!=0 && p<=3) continue;
      create_SIPCut_YieldRatios_single(p, yield_range[rr][0], yield_range[rr][1]);
    }
  }
}

float compute_subtractedPV(
  float CandVtx_x,
  float CandVtx_y,
  float CandVtx_z,
  float CandVtx_cov_xx,
  float CandVtx_cov_xy,
  float CandVtx_cov_xz,
  float CandVtx_cov_yy,
  float CandVtx_cov_yz,
  float CandVtx_cov_zz,
  float PrimaryVtx_x,
  float PrimaryVtx_y,
  float PrimaryVtx_z,
  float PrimaryVtx_cov_xx,
  float PrimaryVtx_cov_xy,
  float PrimaryVtx_cov_xz,
  float PrimaryVtx_cov_yy,
  float PrimaryVtx_cov_yz,
  float PrimaryVtx_cov_zz,
  float& SubtractedVtx_x,
  float& SubtractedVtx_y,
  float& SubtractedVtx_z
  );


bool applyZXloose(
  int CRflag,
  float Lep1combRelIsoPF,
  float Lep2combRelIsoPF,
  float Lep3combRelIsoPF,
  float Lep4combRelIsoPF,
  bool Lep3isID,
  bool Lep4isID
  ){
  bool admit=false;
  if (
    (CRflag==6 || CRflag==10) ||
    (CRflag==8 || CRflag==12) ||
    (CRflag==7 || CRflag==11) ||
    (CRflag==5 || CRflag==9)) admit = true;
  else return admit;

  if (TMath::Max(Lep1combRelIsoPF, TMath::Max(Lep2combRelIsoPF, TMath::Max(Lep3combRelIsoPF, Lep4combRelIsoPF)))<0.4 && Lep3isID && Lep4isID &&
    (
    (CRflag==6 || CRflag==10) || (CRflag==8 || CRflag==12)
    )
    ) return false;
  else return true;
}
bool ZXchannelselection(
  int folder,
  short Z1ids,
  int CRflag){
  bool admit=false;
  if (folder==0 && (Z1ids==-169 && (CRflag==6 || CRflag==10 || CRflag==5 || CRflag==9))) admit=true;
  if (folder==1 && (Z1ids==-121 && (CRflag==8 || CRflag==12 || CRflag==7 || CRflag==11))) admit=true;
  if (folder==2 &&
    (
    (Z1ids==-169 && (CRflag==8 || CRflag==12 || CRflag==7 || CRflag==11)) ||
    (Z1ids==-121 && (CRflag==6 || CRflag==10 || CRflag==5 || CRflag==9))
    )
    ) admit=true;
  return admit;
}



void produce_SIPVariables(int iProcess, int iRegion = 0){ // iRegion: 0 for SR, 1-3 for SB 1-3
  gROOT->ProcessLine(".x tdrstyle.cc");
  char TREE_NAME[] = "SelectedTree";
  int processMap[7] ={
    0, 1, 2, 3,
    kGGSamples, // qqZZ full range samples; 6 of them
    kGGOLDSamples, // MCFM gg, 3 of them
    (kAllSamples - 1)
  };
  float rangeZZMass[2] ={ systZZMass_range[iRegion][0], systZZMass_range[iRegion][1] };

  float ZXfake_weight_SS;
  float ZXfake_weight_OS[6];
  float MC_weight;
  float MC_weight_noxsec;
  float MC_weight_QQZZEWK = 1;
  float MC_weight_Kfactor = 1;
  float MC_weight_QQBGGProper[4];
  float ZZMass, ZZPt, ZZEta, ZZPhi;
  float Lep1_Z1SIP, Lep2_Z1SIP, Lep3_Z1SIP, Lep4_Z1SIP;
  float Lep1SIP, Lep2SIP, Lep3SIP, Lep4SIP;
  int Lep1ID, Lep2ID, Lep3ID, Lep4ID;
  float KalmanCandVtx_x, KalmanCandVtx_y, KalmanCandVtx_z, KalmanCandVtx_chi2;
  float KalmanCandVtx_cov_xx;
  float KalmanCandVtx_cov_xy;
  float KalmanCandVtx_cov_xz;
  float KalmanCandVtx_cov_yy;
  float KalmanCandVtx_cov_yz;
  float KalmanCandVtx_cov_zz;
  float Z1CandVtx_x, Z1CandVtx_y, Z1CandVtx_z;
  float Z1CandVtx_cov_xx;
  float Z1CandVtx_cov_xy;
  float Z1CandVtx_cov_xz;
  float Z1CandVtx_cov_yy;
  float Z1CandVtx_cov_yz;
  float Z1CandVtx_cov_zz;
  float OfflinePrimaryVtx_x, OfflinePrimaryVtx_y, OfflinePrimaryVtx_z, OfflinePrimaryVtx_ndof;
  float OfflinePrimaryVtx_cov_xx;
  float OfflinePrimaryVtx_cov_xy;
  float OfflinePrimaryVtx_cov_xz;
  float OfflinePrimaryVtx_cov_yy;
  float OfflinePrimaryVtx_cov_yz;
  float OfflinePrimaryVtx_cov_zz;

  int RunNumber, LumiNumber, RunNumber_Ref, LumiNumber_Ref;
  float BeamPosX, BeamPosY, BeamPosZ, BeamPosXErr, BeamPosYErr, BeamPosZErr;

  float D_bkg, p0plus_VAJHU, bkg_VAMCFM, p0plus_m4l, bkg_m4l;

  float GenPrimaryVtx_x, GenPrimaryVtx_y, GenPrimaryVtx_z;
  float GenIntVtx_x, GenIntVtx_y, GenIntVtx_z;
  float GenHMass, GenHPt, GenHPhi;

  float SubtractedPrimaryVtx_x, SubtractedPrimaryVtx_y, SubtractedPrimaryVtx_z;
  float Txy, Dxy, Txy_true, Dxy_true, delTxy, delDxy;
  float Txy_BS, Dxy_BS;
  float sigmaPV_xy, sigmaInt_xy;
  float CandVtx_x, CandVtx_y, CandVtx_z;
  float CandBSVtx_x, CandBSVtx_y, CandBSVtx_z;
  float TrueCandVtx_x, TrueCandVtx_y, TrueCandVtx_z;
  float sigmaDirectional_ratio;

  int CRflag=-1;
  float Lep1combRelIsoPF=0, Lep2combRelIsoPF=0, Lep3combRelIsoPF=0, Lep4combRelIsoPF=0;
  bool Lep3isID=true, Lep4isID=true;
  short Z1ids;

  TH1F* hvars[nVariables][newSIPCuts] ={ { 0 } };
  TH2F* hvars_2D[nVariables_2D][newSIPCuts] ={ { 0 } };
  TString hname_core = "h_";
  for (int v = 0; v < nVariables; v++){
    float binwidth = (KD_ranges[v][1] - KD_ranges[v][0]) / KD_nbins[v];
    int nbins = KD_nbins[v]+1;
    float xmin = KD_ranges[v][0];
    float xmax = KD_ranges[v][1] + binwidth;
    if (xmin < 0){
      xmin -= binwidth; nbins++;
    }

    for (int sc = 0; sc < newSIPCuts; sc++){
      TString hname = hname_core;
      hname.Append(varNames[v]);
      hname.Append("_");
      hname.Append(cutNames[sc + 1]);
      hvars[v][sc] = new TH1F(hname, "", nbins, xmin, xmax);
      hvars[v][sc]->Sumw2();
      hvars[v][sc]->SetXTitle(varTitles[v]);
      
      TString ytitle = "Events / ";
      int binwidth_int = (int)binwidth;
      float binwidth_intfloat = (float)binwidth_int;
      float binwidth10 = 10.*binwidth;
      int binwidth10_int = (int)binwidth10;
      float binwidth10_intfloat = (float)binwidth10_int;
      if (binwidth == binwidth_intfloat) ytitle.Append(Form("%.0f", binwidth));
      else if (binwidth<1 && binwidth10==binwidth10_intfloat) ytitle.Append(Form("%.1f", binwidth));
      else ytitle.Append(Form("%.2f", binwidth));
      if (v==vdiffDotTransverse_DirSubPV || v==vTxy || v==vTxy_BS) ytitle.Append(" #mum");
      if (v==vnOPVTracks) ytitle.Append(" tracks");
      hvars[v][sc]->SetYTitle(ytitle);

      hvars[v][sc]->GetXaxis()->SetNdivisions(505);
      hvars[v][sc]->GetXaxis()->SetLabelFont(42);
      hvars[v][sc]->GetXaxis()->SetLabelOffset(0.007);
      hvars[v][sc]->GetXaxis()->SetLabelSize(0.04);
      hvars[v][sc]->GetXaxis()->SetTitleSize(0.06);
      hvars[v][sc]->GetXaxis()->SetTitleOffset(0.9);
      hvars[v][sc]->GetXaxis()->SetTitleFont(42);
      hvars[v][sc]->GetYaxis()->SetNdivisions(505);
      hvars[v][sc]->GetYaxis()->SetLabelFont(42);
      hvars[v][sc]->GetYaxis()->SetLabelOffset(0.007);
      hvars[v][sc]->GetYaxis()->SetLabelSize(0.04);
      hvars[v][sc]->GetYaxis()->SetTitleSize(0.06);
      hvars[v][sc]->GetYaxis()->SetTitleOffset(1.1);
      hvars[v][sc]->GetYaxis()->SetTitleFont(42);
    }
  }
  for (int v = 0; v < nVariables_2D; v++){
    int vx=KD_pairing[v][0];
    float binwidth_x = (KD_ranges[vx][1] - KD_ranges[vx][0]) / KD_nbins[vx];
    int nbins_x = KD_nbins[vx]+1;
    float xmin = KD_ranges[vx][0];
    float xmax = KD_ranges[vx][1] + binwidth_x;
    if (xmin < 0){
      xmin -= binwidth_x; nbins_x++;
    }

    int vy=KD_pairing[v][1];
    float binwidth_y = (KD_ranges[vy][1] - KD_ranges[vy][0]) / KD_nbins[vy];
    int nbins_y = KD_nbins[vy]+1;
    float ymin = KD_ranges[vy][0];
    float ymax = KD_ranges[vy][1] + binwidth_y;
    if (ymin < 0){
      ymin -= binwidth_y; nbins_y++;
    }

    for (int sc = 0; sc < newSIPCuts; sc++){
      TString hname = hname_core;
      hname.Append(varNames[vy]);
      hname.Append("_vs_");
      hname.Append(varNames[vx]);
      hname.Append("_");
      hname.Append(cutNames[sc + 1]);
      hvars_2D[v][sc] = new TH2F(hname, "", nbins_x, xmin, xmax, nbins_y, ymin, ymax);
      hvars_2D[v][sc]->Sumw2();
      hvars_2D[v][sc]->SetXTitle(varTitles[vx]);
      hvars_2D[v][sc]->SetYTitle(varTitles[vy]);
      hvars_2D[v][sc]->GetXaxis()->SetNdivisions(505);
      hvars_2D[v][sc]->GetXaxis()->SetLabelFont(42);
      hvars_2D[v][sc]->GetXaxis()->SetLabelOffset(0.007);
      hvars_2D[v][sc]->GetXaxis()->SetLabelSize(0.04);
      hvars_2D[v][sc]->GetXaxis()->SetTitleSize(0.06);
      hvars_2D[v][sc]->GetXaxis()->SetTitleOffset(1.0);
      hvars_2D[v][sc]->GetXaxis()->SetTitleFont(42);
      hvars_2D[v][sc]->GetYaxis()->SetNdivisions(505);
      hvars_2D[v][sc]->GetYaxis()->SetLabelFont(42);
      hvars_2D[v][sc]->GetYaxis()->SetLabelOffset(0.007);
      hvars_2D[v][sc]->GetYaxis()->SetLabelSize(0.04);
      hvars_2D[v][sc]->GetYaxis()->SetTitleSize(0.06);
      hvars_2D[v][sc]->GetYaxis()->SetTitleOffset(1.1);
      hvars_2D[v][sc]->GetYaxis()->SetTitleFont(42);
    }
  }


  for (int erg_tev = 7; erg_tev < 9; erg_tev++){
    TString erg_dir;
    erg_dir.Form("LHC_%iTeV", erg_tev);
    TString comstring;
    comstring.Form("%iTeV", erg_tev);
    int EnergyIndex = 0; if (erg_tev == 8) EnergyIndex = 1;

    // BeamSpot info
    TChain* tBeam[2];
    TString strBeamSpot[2][2]={
      { "data_GR_R_44_V15C", "MC_START44_V13" },
      { "data_FT_53_V21_AN4", "MC_START53_V23" },
    };
    for (int bs=0; bs<2; bs++){
      tBeam[bs] = new TChain("BeamSpotRecord");
      TString cinput_BS = "./data/BeamSpotRecord_" + strBeamSpot[EnergyIndex][bs] + "_" + comstring + ".root";
      tBeam[bs]->Add(cinput_BS);
      tBeam[bs]->SetBranchAddress("RunNumber", &RunNumber_Ref);
      tBeam[bs]->SetBranchAddress("LumiNumber", &LumiNumber_Ref);
      tBeam[bs]->SetBranchAddress("BeamPosX", &BeamPosX);
      tBeam[bs]->SetBranchAddress("BeamPosY", &BeamPosY);
      tBeam[bs]->SetBranchAddress("BeamPosZ", &BeamPosZ);
      tBeam[bs]->SetBranchAddress("BeamPosXErr", &BeamPosXErr);
      tBeam[bs]->SetBranchAddress("BeamPosYErr", &BeamPosYErr);
      tBeam[bs]->SetBranchAddress("BeamPosZErr", &BeamPosZErr);
    }

    for (int folder = 0; folder < 3; folder++){
      TString OUTPUT_NAME = "LifetimeKD_SIPVariables_";
      OUTPUT_NAME.Append(Form("%s", processName[iProcess].c_str()));
      OUTPUT_NAME.Append(Form("_%s_", user_folder[folder]));
      OUTPUT_NAME.Append(comstring);
      if (iRegion == 0) OUTPUT_NAME.Append("_SR");
      else if (iRegion <= 3) OUTPUT_NAME.Append(Form("_SB%i", iRegion));
      else assert(0);
      OUTPUT_NAME.Append(".root");

      TString cinput_common_withSIP = user_dir_hep + erg_dir + "/";
      TString cinput_common_noSIP = user_dir_hep + "No_SIP/" + erg_dir + "/";
      if (iProcess < 7){
        cinput_common_withSIP.Append(Form("%s/", user_folder[folder]));
        cinput_common_noSIP.Append(Form("%s/", user_folder[folder]));
      }
      else{
        cinput_common_withSIP.Append(Form("%s/", user_folder[3]));
        cinput_common_noSIP.Append(Form("%s/", user_folder[3]));
      }
      TString cinput_common = cinput_common_noSIP;
      TString coutput_common = user_dir_hep + "Analysis/Auxiliary/";

      TString plotDir = coutput_common;
      plotDir.Append("../Plots/SIP/");
      if (iRegion==0) plotDir.Append("SignalRegion/");
      else plotDir.Append(Form("SidebandRegion%i/", iRegion));
      plotDir.Append(Form("%s/", processName[iProcess].c_str()));
      TString plotDir_1D=plotDir;
      TString plotDir_2D=plotDir;
      plotDir_1D.Append("1D/");
      plotDir_2D.Append("2D/");
      TString mkdirCommand = "mkdir -p ";
      TString mkdirCommand_1D = mkdirCommand;
      TString mkdirCommand_2D = mkdirCommand;
      mkdirCommand_1D.Append(plotDir_1D);
      gSystem->Exec(mkdirCommand_1D);
      mkdirCommand_2D.Append(plotDir_2D);
      gSystem->Exec(mkdirCommand_2D);

      string cinput_aux = user_dir_hep + "Analysis/Auxiliary/";
      TString SSinputname = "OSoverSS_";
      SSinputname.Append("MCCR_");
      SSinputname += comstring;
      SSinputname.Append(".root");
      SSinputname.Prepend(cinput_aux.c_str());
      TFile* fSSinput = new TFile(SSinputname, "read");
      TH2F* hRatio[newSIPCuts + 1];
      for (int hr=0; hr<newSIPCuts + 1; hr++){
        TString chratio = "hCR_MC_OSSSRatios_";
        chratio.Append(cutLabel[hr]);
        hRatio[hr] = (TH2F*)fSSinput->Get(chratio);
        cout << "Obtained ratio histogram " << hRatio[hr]->GetName() << endl;
      }
      double ratio_targetsignalyield[8][newSIPCuts + 1][4] ={ { { 0 } } };
      double targetCRcount[2][newSIPCuts + 1][4] ={ { { 0 } } };
      float m4l_lowhigh[4][2]={
        { 105.6, 140.6 },
        { 70., 105.6 },
        { 140.6, 170. },
        { 170., 3000. }
      };
      for (int p = 0; p <= 6; p++){
        for (int ir=0; ir<4; ir++){
          TString INPUTYIELD_NAME = "LifetimeKD_RelativeSIPYields_";
          INPUTYIELD_NAME.Append(Form("%s", processName[p].c_str()));
          INPUTYIELD_NAME.Append(Form("_%s_", user_folder[folder]));
          INPUTYIELD_NAME.Append(comstring);
          INPUTYIELD_NAME.Append(Form("_m4l%.1f_%.1f", m4l_lowhigh[ir][0], m4l_lowhigh[ir][1]));
          INPUTYIELD_NAME.Append(".root");
          TString cinput_yield_common = user_dir_hep + "Analysis/Auxiliary/";
          TString cinput_yield = cinput_yield_common + INPUTYIELD_NAME;
          TFile* finput = new TFile(cinput_yield, "read");
          if (finput==0 || finput->IsZombie()){
            cout << cinput_yield << " not found." << endl; continue;
          }
          TH1F* htemp;
          if (p!=6){
            htemp = (TH1F*)finput->Get(Form("h%s", processName[p].c_str()));
            for (int b = 1; b <= newSIPCuts + 1; b++) ratio_targetsignalyield[p][b-1][ir] = htemp->GetBinContent(b);
            delete htemp;
          }
          else{
            htemp = (TH1F*)finput->Get(Form("h%s_OS", processName[p].c_str()));
            for (int b = 1; b <= newSIPCuts + 1; b++) ratio_targetsignalyield[p][b-1][ir] = htemp->GetBinContent(b);
            delete htemp;
            htemp = (TH1F*)finput->Get(Form("h%s_SS", processName[p].c_str()));
            for (int b = 1; b <= newSIPCuts + 1; b++) ratio_targetsignalyield[p][b-1][ir] += htemp->GetBinContent(b);
            delete htemp;

            htemp = (TH1F*)finput->Get(Form("h%s_Unweighted_OS", processName[p].c_str()));
            for (int b = 1; b <= newSIPCuts + 1; b++) targetCRcount[0][b-1][ir] = htemp->GetBinContent(b);
            delete htemp;
            htemp = (TH1F*)finput->Get(Form("h%s_Unweighted_SS", processName[p].c_str()));
            for (int b = 1; b <= newSIPCuts + 1; b++) targetCRcount[1][b-1][ir] = htemp->GetBinContent(b);
            delete htemp;
          }
          finput->Close();
        }
      }
      for (int b = 1; b <= newSIPCuts + 1; b++){
        ratio_targetsignalyield[5][b-1][1] = ratio_targetsignalyield[5][b-1][0];
        ratio_targetsignalyield[6][b-1][1] = ratio_targetsignalyield[6][b-1][0];
        targetCRcount[0][b-1][1] = targetCRcount[0][b-1][0];
        targetCRcount[1][b-1][1] = targetCRcount[1][b-1][0];
      }

      TString coutput = coutput_common + OUTPUT_NAME;
      TFile* foutput = new TFile(coutput, "recreate");

      TChain* tc;
      TChain* tc_extra=0;
      TString cinput;
      int smp = processMap[iProcess];
      bool treesExist = true;
      tc = new TChain(TREE_NAME);

      if (iProcess < 4){
        cinput = cinput_common + sample_FullSim[smp] + ".root";
        tc->Add(cinput);
      }
      else if (iProcess == 4 || iProcess == 5){
        int iProcess_new = iProcess;
        if (iProcess==5 && iRegion==1) iProcess_new--;
        smp = processMap[iProcess_new];
        for (int it = 0; it < (iProcess_new == 4 ? 6 : 3); it++){ // qqZZ has 6 samples, gg has 3 samples
          cinput = cinput_common + sample_FullSim[smp] + "_Reprocessed.root";
          if (iProcess==5 && iRegion==1) cout << "Region 1 gg: " << cinput << endl;;
          tc->Add(cinput);
          smp++;
        }
      }
      else if (iProcess == 6){
        cinput = user_dir_hep;
        cinput.Append("No_SIP/");
        cinput.Append(erg_dir);
        cinput = cinput + "/CR/" + sample_FullSim[smp] + ".root";
        tc->Add(cinput);
      }
      else if (iProcess == 7){
        cinput = cinput_common + data_files[folder] + ".root";
        tc->Add(cinput);
      }
      if (iRegion==0 && (iProcess==4 || iProcess==5)){
        tc_extra = new TChain(TREE_NAME);
        if (iProcess==4 || iProcess==5){
          for (int it = kQQBZZSamples; it < kQQBZZSamples_Dedicated; it++){ // qq/ggZZ extra for signal region
            cinput = cinput_common + sample_FullSim[it] + "_Reprocessed.root";
            tc_extra->Add(cinput);
            smp++;
          }
        }
        if (iProcess==5){
          for (int it = kGGHSamples; it < kQQBZZSamples; it++){ // ggZZ extra for signal region
            if (it >= kGGOLDSamples && it < kGGMCFMSamples) continue;
            cinput = cinput_common + sample_FullSim[it] + "_Reprocessed.root";
            tc_extra->Add(cinput);
            smp++;
          }
        }
      }

      if (tc->GetEntries() == 0){ cout << "Could not find the file, aborting..." << endl; treesExist = false; }
      else{
        if (tc->GetBranchStatus("MC_weight_QQBGGProper")) tc->SetBranchAddress("MC_weight_QQBGGProper", MC_weight_QQBGGProper);
        if (tc->GetBranchStatus("MC_weight_QQZZEWK")) tc->SetBranchAddress("MC_weight_QQZZEWK", &MC_weight_QQZZEWK);
        if (tc->GetBranchStatus("MC_weight_Kfactor")) tc->SetBranchAddress("MC_weight_Kfactor", &MC_weight_Kfactor);

        if (iProcess < 6){
          tc->SetBranchAddress("MC_weight", &MC_weight);
          tc->SetBranchAddress("MC_weight_noxsec", &MC_weight_noxsec);
        }
        else if (iProcess == 6){
//          if (tc->GetBranchStatus("ZXfake_weightProper")) tc->SetBranchAddress("ZXfake_weightProper", &MC_weight);
//          else{
//            if (tc->GetBranchStatus("ZXfake_weight_SS")){
              tc->SetBranchAddress("ZXfake_weight_OS", ZXfake_weight_OS);
              if (tc->GetBranchStatus("ZXfake_weight_SS")) tc->SetBranchAddress("ZXfake_weight_SS", &ZXfake_weight_SS);
              else cout << "NOT SET" << endl;
//            }
//            else tc->SetBranchAddress("ZXfake_weight", &MC_weight);
            tc->SetBranchAddress("CRflag", &CRflag);

            tc->SetBranchAddress("Lep1combRelIsoPF", &Lep1combRelIsoPF);
            tc->SetBranchAddress("Lep2combRelIsoPF", &Lep2combRelIsoPF);
            tc->SetBranchAddress("Lep3combRelIsoPF", &Lep3combRelIsoPF);
            tc->SetBranchAddress("Lep4combRelIsoPF", &Lep4combRelIsoPF);

            tc->SetBranchAddress("Lep3isID", &Lep3isID);
            tc->SetBranchAddress("Lep4isID", &Lep4isID);

            tc->SetBranchAddress("Z1ids", &Z1ids);
            tc->SetBranchAddress("Lep3ID", &Lep3ID);
            tc->SetBranchAddress("Lep4ID", &Lep4ID);
//          }
        }
        if (iProcess>=6){
          tc->SetBranchAddress("RunNumber", &RunNumber);
          tc->SetBranchAddress("LumiNumber", &LumiNumber);
        }

        tc->SetBranchAddress("ZZMass", &ZZMass);
        tc->SetBranchAddress("ZZPt", &ZZPt);
        tc->SetBranchAddress("ZZEta", &ZZEta);
        tc->SetBranchAddress("ZZPhi", &ZZPhi);
        tc->SetBranchAddress("Lep1SIP", &Lep1SIP);
        tc->SetBranchAddress("Lep2SIP", &Lep2SIP);
        tc->SetBranchAddress("Lep3SIP", &Lep3SIP);
        tc->SetBranchAddress("Lep4SIP", &Lep4SIP);
        tc->SetBranchAddress("Lep1_Z1SIP", &Lep1_Z1SIP);
        tc->SetBranchAddress("Lep2_Z1SIP", &Lep2_Z1SIP);
        tc->SetBranchAddress("Lep3_Z1SIP", &Lep3_Z1SIP);
        tc->SetBranchAddress("Lep4_Z1SIP", &Lep4_Z1SIP);
        tc->SetBranchAddress("OfflinePrimaryVtx_x", &OfflinePrimaryVtx_x);
        tc->SetBranchAddress("OfflinePrimaryVtx_y", &OfflinePrimaryVtx_y);
        tc->SetBranchAddress("OfflinePrimaryVtx_z", &OfflinePrimaryVtx_z);
        tc->SetBranchAddress("OfflinePrimaryVtx_cov_xx", &OfflinePrimaryVtx_cov_xx);
        tc->SetBranchAddress("OfflinePrimaryVtx_cov_xy", &OfflinePrimaryVtx_cov_xy);
        tc->SetBranchAddress("OfflinePrimaryVtx_cov_xz", &OfflinePrimaryVtx_cov_xz);
        tc->SetBranchAddress("OfflinePrimaryVtx_cov_yy", &OfflinePrimaryVtx_cov_yy);
        tc->SetBranchAddress("OfflinePrimaryVtx_cov_yz", &OfflinePrimaryVtx_cov_yz);
        tc->SetBranchAddress("OfflinePrimaryVtx_cov_zz", &OfflinePrimaryVtx_cov_zz);
        tc->SetBranchAddress("OfflinePrimaryVtx_ndof", &OfflinePrimaryVtx_ndof);
        tc->SetBranchAddress("KalmanCandVtx_x", &KalmanCandVtx_x);
        tc->SetBranchAddress("KalmanCandVtx_y", &KalmanCandVtx_y);
        tc->SetBranchAddress("KalmanCandVtx_z", &KalmanCandVtx_z);
        tc->SetBranchAddress("KalmanCandVtx_chi2", &KalmanCandVtx_chi2);
        tc->SetBranchAddress("KalmanCandVtx_cov_xx", &KalmanCandVtx_cov_xx);
        tc->SetBranchAddress("KalmanCandVtx_cov_xy", &KalmanCandVtx_cov_xy);
        tc->SetBranchAddress("KalmanCandVtx_cov_xz", &KalmanCandVtx_cov_xz);
        tc->SetBranchAddress("KalmanCandVtx_cov_yy", &KalmanCandVtx_cov_yy);
        tc->SetBranchAddress("KalmanCandVtx_cov_yz", &KalmanCandVtx_cov_yz);
        tc->SetBranchAddress("KalmanCandVtx_cov_zz", &KalmanCandVtx_cov_zz);
        tc->SetBranchAddress("Z1CandVtx_x", &Z1CandVtx_x);
        tc->SetBranchAddress("Z1CandVtx_y", &Z1CandVtx_y);
        tc->SetBranchAddress("Z1CandVtx_z", &Z1CandVtx_z);
        tc->SetBranchAddress("Z1CandVtx_cov_xx", &Z1CandVtx_cov_xx);
        tc->SetBranchAddress("Z1CandVtx_cov_xy", &Z1CandVtx_cov_xy);
        tc->SetBranchAddress("Z1CandVtx_cov_xz", &Z1CandVtx_cov_xz);
        tc->SetBranchAddress("Z1CandVtx_cov_yy", &Z1CandVtx_cov_yy);
        tc->SetBranchAddress("Z1CandVtx_cov_yz", &Z1CandVtx_cov_yz);
        tc->SetBranchAddress("Z1CandVtx_cov_zz", &Z1CandVtx_cov_zz);

        if (tc->GetBranchStatus("GenPrimaryVtx_x")){
          tc->SetBranchAddress("GenPrimaryVtx_x", &GenPrimaryVtx_x);
          tc->SetBranchAddress("GenPrimaryVtx_y", &GenPrimaryVtx_y);
          tc->SetBranchAddress("GenPrimaryVtx_z", &GenPrimaryVtx_z);
          tc->SetBranchAddress("GenIntVtx_x", &GenIntVtx_x);
          tc->SetBranchAddress("GenIntVtx_y", &GenIntVtx_y);
          tc->SetBranchAddress("GenIntVtx_z", &GenIntVtx_z);
          tc->SetBranchAddress("GenHPt", &GenHPt);
          tc->SetBranchAddress("GenHPhi", &GenHPhi);
          tc->SetBranchAddress("GenHMass", &GenHMass);
        }

        if (tc->GetBranchStatus("D_bkg")){
          tc->SetBranchAddress("D_bkg", &D_bkg);
        }
        else if (tc->GetBranchStatus("p0plus_VAJHU")){
          tc->SetBranchAddress("p0plus_VAJHU", &p0plus_VAJHU);
          tc->SetBranchAddress("bkg_VAMCFM", &bkg_VAMCFM);
          tc->SetBranchAddress("p0plus_m4l", &p0plus_m4l);
          tc->SetBranchAddress("bkg_m4l", &bkg_m4l);
        }

      }
      if (tc_extra!=0){
        if (tc_extra->GetEntries() == 0){ cout << "Could not find the file, aborting..." << endl; treesExist = false; }
        else{
          tc_extra->SetBranchAddress("MC_weight", &MC_weight);
          tc_extra->SetBranchAddress("MC_weight_noxsec", &MC_weight_noxsec);
          if (tc_extra->GetBranchStatus("MC_weight_QQBGGProper")) tc_extra->SetBranchAddress("MC_weight_QQBGGProper", MC_weight_QQBGGProper);
          if (tc_extra->GetBranchStatus("MC_weight_QQZZEWK")) tc_extra->SetBranchAddress("MC_weight_QQZZEWK", &MC_weight_QQBGGProper);
          if (tc_extra->GetBranchStatus("MC_weight_Kfactor")) tc_extra->SetBranchAddress("MC_weight_Kfactor", &MC_weight_Kfactor);

          tc_extra->SetBranchAddress("ZZMass", &ZZMass);
          tc_extra->SetBranchAddress("ZZPt", &ZZPt);
          tc_extra->SetBranchAddress("ZZEta", &ZZEta);
          tc_extra->SetBranchAddress("ZZPhi", &ZZPhi);
          tc_extra->SetBranchAddress("Lep1SIP", &Lep1SIP);
          tc_extra->SetBranchAddress("Lep2SIP", &Lep2SIP);
          tc_extra->SetBranchAddress("Lep3SIP", &Lep3SIP);
          tc_extra->SetBranchAddress("Lep4SIP", &Lep4SIP);
          tc_extra->SetBranchAddress("Lep1_Z1SIP", &Lep1_Z1SIP);
          tc_extra->SetBranchAddress("Lep2_Z1SIP", &Lep2_Z1SIP);
          tc_extra->SetBranchAddress("Lep3_Z1SIP", &Lep3_Z1SIP);
          tc_extra->SetBranchAddress("Lep4_Z1SIP", &Lep4_Z1SIP);
          tc_extra->SetBranchAddress("OfflinePrimaryVtx_x", &OfflinePrimaryVtx_x);
          tc_extra->SetBranchAddress("OfflinePrimaryVtx_y", &OfflinePrimaryVtx_y);
          tc_extra->SetBranchAddress("OfflinePrimaryVtx_z", &OfflinePrimaryVtx_z);
          tc_extra->SetBranchAddress("OfflinePrimaryVtx_cov_xx", &OfflinePrimaryVtx_cov_xx);
          tc_extra->SetBranchAddress("OfflinePrimaryVtx_cov_xy", &OfflinePrimaryVtx_cov_xy);
          tc_extra->SetBranchAddress("OfflinePrimaryVtx_cov_xz", &OfflinePrimaryVtx_cov_xz);
          tc_extra->SetBranchAddress("OfflinePrimaryVtx_cov_yy", &OfflinePrimaryVtx_cov_yy);
          tc_extra->SetBranchAddress("OfflinePrimaryVtx_cov_yz", &OfflinePrimaryVtx_cov_yz);
          tc_extra->SetBranchAddress("OfflinePrimaryVtx_cov_zz", &OfflinePrimaryVtx_cov_zz);
          tc_extra->SetBranchAddress("OfflinePrimaryVtx_ndof", &OfflinePrimaryVtx_ndof);
          tc_extra->SetBranchAddress("KalmanCandVtx_x", &KalmanCandVtx_x);
          tc_extra->SetBranchAddress("KalmanCandVtx_y", &KalmanCandVtx_y);
          tc_extra->SetBranchAddress("KalmanCandVtx_z", &KalmanCandVtx_z);
          tc_extra->SetBranchAddress("KalmanCandVtx_chi2", &KalmanCandVtx_chi2);
          tc_extra->SetBranchAddress("KalmanCandVtx_cov_xx", &KalmanCandVtx_cov_xx);
          tc_extra->SetBranchAddress("KalmanCandVtx_cov_xy", &KalmanCandVtx_cov_xy);
          tc_extra->SetBranchAddress("KalmanCandVtx_cov_xz", &KalmanCandVtx_cov_xz);
          tc_extra->SetBranchAddress("KalmanCandVtx_cov_yy", &KalmanCandVtx_cov_yy);
          tc_extra->SetBranchAddress("KalmanCandVtx_cov_yz", &KalmanCandVtx_cov_yz);
          tc_extra->SetBranchAddress("KalmanCandVtx_cov_zz", &KalmanCandVtx_cov_zz);
          tc_extra->SetBranchAddress("Z1CandVtx_x", &Z1CandVtx_x);
          tc_extra->SetBranchAddress("Z1CandVtx_y", &Z1CandVtx_y);
          tc_extra->SetBranchAddress("Z1CandVtx_z", &Z1CandVtx_z);
          tc_extra->SetBranchAddress("Z1CandVtx_cov_xx", &Z1CandVtx_cov_xx);
          tc_extra->SetBranchAddress("Z1CandVtx_cov_xy", &Z1CandVtx_cov_xy);
          tc_extra->SetBranchAddress("Z1CandVtx_cov_xz", &Z1CandVtx_cov_xz);
          tc_extra->SetBranchAddress("Z1CandVtx_cov_yy", &Z1CandVtx_cov_yy);
          tc_extra->SetBranchAddress("Z1CandVtx_cov_yz", &Z1CandVtx_cov_yz);
          tc_extra->SetBranchAddress("Z1CandVtx_cov_zz", &Z1CandVtx_cov_zz);

          if (tc_extra->GetBranchStatus("GenPrimaryVtx_x")){
            tc_extra->SetBranchAddress("GenPrimaryVtx_x", &GenPrimaryVtx_x);
            tc_extra->SetBranchAddress("GenPrimaryVtx_y", &GenPrimaryVtx_y);
            tc_extra->SetBranchAddress("GenPrimaryVtx_z", &GenPrimaryVtx_z);
            tc_extra->SetBranchAddress("GenIntVtx_x", &GenIntVtx_x);
            tc_extra->SetBranchAddress("GenIntVtx_y", &GenIntVtx_y);
            tc_extra->SetBranchAddress("GenIntVtx_z", &GenIntVtx_z);
            tc_extra->SetBranchAddress("GenHPt", &GenHPt);
            tc_extra->SetBranchAddress("GenHPhi", &GenHPhi);
            tc_extra->SetBranchAddress("GenHMass", &GenHMass);
          }

          if (tc_extra->GetBranchStatus("D_bkg")){
            tc_extra->SetBranchAddress("D_bkg", &D_bkg);
          }
          else if (tc_extra->GetBranchStatus("p0plus_VAJHU")){
            tc_extra->SetBranchAddress("p0plus_VAJHU", &p0plus_VAJHU);
            tc_extra->SetBranchAddress("bkg_VAMCFM", &bkg_VAMCFM);
            tc_extra->SetBranchAddress("p0plus_m4l", &p0plus_m4l);
            tc_extra->SetBranchAddress("bkg_m4l", &bkg_m4l);
          }

        }
      }

      if (!treesExist){
        foutput->Close();
        TString rmCommand = "rm ";
        rmCommand.Append(coutput);
        gSystem->Exec(rmCommand);
        continue;
      }

      for (int sc = 0; sc < newSIPCuts; sc++){
        for (int v = 0; v < nVariables; v++){
          hvars[v][sc]->Reset("ICESM");
          hvars[v][sc]->Sumw2();
        }
        for (int v = 0; v < nVariables_2D; v++){
          hvars_2D[v][sc]->Reset("ICESM");
          hvars_2D[v][sc]->Sumw2();
        }
      }

      double sum_offshell[newSIPCuts] ={ 0 };
      double sum_80_100[newSIPCuts] ={ 0 };
      double sum_signal[newSIPCuts] ={ 0 };

      for (int ev = 0; ev < tc->GetEntries(); ev++){
        RunNumber = 1;
        LumiNumber = 1;

        MC_weight = 1;
        MC_weight_noxsec = 1;
        MC_weight_QQZZEWK = 1;
        MC_weight_Kfactor = 1;

        GenPrimaryVtx_x=0;
        GenPrimaryVtx_y=0;
        GenPrimaryVtx_z=0;
        GenIntVtx_x=0;
        GenIntVtx_y=0;
        GenIntVtx_z=0;

        GenHMass=125.6;
        GenHPt=20;
        GenHPhi=0;
        Lep1_Z1SIP=0; Lep2_Z1SIP=0; Lep3_Z1SIP=0; Lep4_Z1SIP=0; KalmanCandVtx_chi2=0;
        Lep1SIP=0; Lep2SIP=0; Lep3SIP=0; Lep4SIP=0;
        CRflag=-1;

        ZXfake_weight_SS=1;

        tc->GetEntry(ev);

        float oldSIPdef = max(max(max(Lep1SIP, Lep2SIP), Lep3SIP), Lep4SIP);

        if (iProcess==4 && iRegion==0) MC_weight = MC_weight_noxsec*MC_weight_QQBGGProper[1];
        if (iProcess==4) MC_weight *= MC_weight_QQZZEWK;

        if (iProcess==5 && iRegion==0) MC_weight = MC_weight_noxsec;
        if (iProcess==5 && iRegion<=1) MC_weight *= MC_weight_QQBGGProper[0];
        if (iProcess==5) MC_weight *= MC_weight_Kfactor;

        if (iProcess==6){
          if (
            !applyZXloose(CRflag,
            Lep1combRelIsoPF,
            Lep2combRelIsoPF,
            Lep3combRelIsoPF,
            Lep4combRelIsoPF,
            Lep3isID,
            Lep4isID)
            ) continue;
          if (!ZXchannelselection(folder, Z1ids, CRflag)) continue;
        }

        int index_BS = 1;
        bool matchFound = false;
        if (iProcess>=6) index_BS = 0;
        for (int evBS=0; evBS < tBeam[index_BS]->GetEntries(); evBS++){
          tBeam[index_BS]->GetEntry(evBS);
          if (RunNumber == RunNumber_Ref && LumiNumber == LumiNumber_Ref){
            matchFound = true;
            break;
          }
          else continue;
        }
        if (!matchFound) cerr << "No beamspot matching possible!" << endl;

        CandBSVtx_x = KalmanCandVtx_x - BeamPosX;
        CandBSVtx_y = KalmanCandVtx_y - BeamPosY;
        CandBSVtx_z = KalmanCandVtx_z - BeamPosZ;
        Dxy_BS = (CandBSVtx_x*cos(ZZPhi) + CandBSVtx_y*sin(ZZPhi));
        Txy_BS = Dxy_BS*ZZMass / ZZPt;

        TrueCandVtx_x = GenIntVtx_x - GenPrimaryVtx_x;
        TrueCandVtx_y = GenIntVtx_y - GenPrimaryVtx_y;
        TrueCandVtx_z = GenIntVtx_z - GenPrimaryVtx_z;
        Dxy_true = (TrueCandVtx_x*cos(GenHPhi) + TrueCandVtx_y*sin(GenHPhi));
        Txy_true = Dxy_true*GenHMass / GenHPt;

        CandVtx_x = KalmanCandVtx_x - OfflinePrimaryVtx_x;
        CandVtx_y = KalmanCandVtx_y - OfflinePrimaryVtx_y;
        CandVtx_z = KalmanCandVtx_z - OfflinePrimaryVtx_z;
        Dxy = (CandVtx_x*cos(ZZPhi) + CandVtx_y*sin(ZZPhi));
        Txy = Dxy*ZZMass / ZZPt;

        sigmaDirectional_ratio = compute_subtractedPV(
          KalmanCandVtx_x,
          KalmanCandVtx_y,
          KalmanCandVtx_z,
          KalmanCandVtx_cov_xx,
          KalmanCandVtx_cov_xy,
          KalmanCandVtx_cov_xz,
          KalmanCandVtx_cov_yy,
          KalmanCandVtx_cov_yz,
          KalmanCandVtx_cov_zz,
          OfflinePrimaryVtx_x,
          OfflinePrimaryVtx_y,
          OfflinePrimaryVtx_z,
          OfflinePrimaryVtx_cov_xx,
          OfflinePrimaryVtx_cov_xy,
          OfflinePrimaryVtx_cov_xz,
          OfflinePrimaryVtx_cov_yy,
          OfflinePrimaryVtx_cov_yz,
          OfflinePrimaryVtx_cov_zz,
          SubtractedPrimaryVtx_x,
          SubtractedPrimaryVtx_y,
          SubtractedPrimaryVtx_z
          );

        if (!tc->GetBranchStatus("D_bkg")) D_bkg = p0plus_VAJHU*p0plus_m4l / (p0plus_VAJHU*p0plus_m4l + bkg_VAMCFM*bkg_m4l);

        KD[vnOPVTracks] = (OfflinePrimaryVtx_ndof+3.)/2.;
        KD[vsigmaDirRatio] = sigmaDirectional_ratio;
        KD[vdiffDotTransverse_DirSubPV] = ((SubtractedPrimaryVtx_x - GenPrimaryVtx_x)*cos(ZZPhi) + (SubtractedPrimaryVtx_y - GenPrimaryVtx_y)*sin(ZZPhi))*10000.;
        KD[vPVSIP] = oldSIPdef;
        KD[vLep1_Z1SIP] = Lep1_Z1SIP;
        KD[vLep2_Z1SIP] = Lep2_Z1SIP;
        KD[vLep3_Z1SIP] = Lep3_Z1SIP;
        KD[vLep4_Z1SIP] = Lep4_Z1SIP;
        KD[v4l_chi2] = KalmanCandVtx_chi2;
        KD[vTxy] = Txy*10000.;
        KD[vTxy_BS] = Txy_BS*10000.;
        KD[vDbkg] = D_bkg;


        for (int sc = 0; sc < newSIPCuts-1; sc++){
          if ((sc == 1 || sc >= 3) && (fabs(Lep1_Z1SIP) >= 4 || fabs(Lep2_Z1SIP) >= 4 || fabs(Lep3_Z1SIP) >= 5 || fabs(Lep4_Z1SIP) >= 5)) continue;
          if ((sc == 2 || sc >= 3) && KalmanCandVtx_chi2 >= 30) continue;
          if (sc == 4  && oldSIPdef>=4) continue;

          if (iProcess==6){
            if (
              (CRflag==6 || CRflag==10) ||
              (CRflag==8 || CRflag==12)
              ){
              MC_weight = ZXfake_weight_OS[sc+1] * (targetCRcount[0][sc+1][iRegion] / (targetCRcount[0][sc+1][iRegion] + targetCRcount[1][sc+1][iRegion]));
            }
            else if (
              (CRflag==7 || CRflag==11) ||
              (CRflag==5 || CRflag==9)
              ){
              if (tc->GetBranchStatus("ZXfake_weight_SS")){
                int OSoverSS_biny = -1;
                if ((CRflag==5 || CRflag==9) && Z1ids==-169) OSoverSS_biny = 1; // 4mu
                if ((CRflag==7 || CRflag==11) && Z1ids==-169) OSoverSS_biny = 3; // 2mu2e
                if ((CRflag==5 || CRflag==9) && Z1ids==-121) OSoverSS_biny = 4; // 2e2mu
                if ((CRflag==7 || CRflag==11) && Z1ids==-121) OSoverSS_biny = 6; // 4e
                float SSwgt_OSoverSS = 1;
                if (OSoverSS_biny>0) SSwgt_OSoverSS = hRatio[sc+1]->GetBinContent(hRatio[sc+1]->GetXaxis()->FindBin(ZZMass), OSoverSS_biny);
//                if (sc==3 && iProcess==6 && iRegion==1 && folder==0 && ZZMass>=110 && ZZMass<130) cout << ZZMass << '\t' << SSwgt_OSoverSS<< '\t' << ZXfake_weight_SS << '\t' << (targetCRcount[1][sc+1][iRegion] / (targetCRcount[0][sc+1][iRegion] + targetCRcount[1][sc+1][iRegion])) << endl;

                MC_weight = ZXfake_weight_SS*SSwgt_OSoverSS;
                MC_weight *= (targetCRcount[1][sc+1][iRegion] / (targetCRcount[0][sc+1][iRegion] + targetCRcount[1][sc+1][iRegion]));
              }
            }
          }

          if ((ZZMass >= 220 && ZZMass < 1600)) sum_offshell[sc] += MC_weight;
          if ((ZZMass >= 80 + (iProcess==6 ? 30 : 0) && ZZMass < 100 + (iProcess==6 ? 30 : 0))){
            sum_80_100[sc] += MC_weight;
//            if (sc==3 && iProcess==6 && iRegion==1) cout << ZZMass << '\t' << KalmanCandVtx_chi2<< '\t' << CRflag << '\t' << MC_weight << endl;
          }
          if ((ZZMass >= 105.6 && ZZMass < 140.6)) sum_signal[sc] += MC_weight;

          if ((ZZMass >= rangeZZMass[1] || ZZMass < rangeZZMass[0]) && !(iRegion==1 && iProcess==6)) continue;
          if ((ZZMass >= rangeZZMass[1]+30 || ZZMass < rangeZZMass[0]+30) && (iRegion==1 && iProcess==6)) continue;
          float fillvar_2D[nVariables_2D][2] ={ { 0 } };
          for (int v = 0; v < nVariables; v++){
            float fillvar = KD[v];
            if (fillvar<hvars[v][sc]->GetXaxis()->GetXmin()) fillvar = hvars[v][sc]->GetXaxis()->GetXmin();
            if (fillvar>=hvars[v][sc]->GetXaxis()->GetXmax()) fillvar = hvars[v][sc]->GetXaxis()->GetXmax() - (hvars[v][sc]->GetXaxis()->GetBinWidth(1))*0.5;

//            if (v==0 && sc==3 && iProcess==6 && iRegion==1) cout << ZZMass << '\t' << KalmanCandVtx_chi2<< '\t' << CRflag << '\t' << MC_weight << endl;
            hvars[v][sc]->Fill(fillvar, MC_weight);

            for (int k = 0; k < nVariables_2D; k++){
              if (v == KD_pairing[k][0]) fillvar_2D[k][0] = fillvar;
              if (v == KD_pairing[k][1]) fillvar_2D[k][1] = fillvar;
            }
          }
          for (int v = 0; v < nVariables_2D; v++) hvars_2D[v][sc]->Fill(fillvar_2D[v][0], fillvar_2D[v][1], MC_weight);
        }
      }

      if (iRegion==0 && tc_extra!=0){
        cout << "Running over extra trees..." << endl;
        for (int ev = 0; ev < tc_extra->GetEntries(); ev++){
          RunNumber = 1;
          LumiNumber = 1;

          MC_weight = 1;
          MC_weight_noxsec = 1;
          MC_weight_QQZZEWK = 1;
          MC_weight_Kfactor = 1;

          GenPrimaryVtx_x=0;
          GenPrimaryVtx_y=0;
          GenPrimaryVtx_z=0;
          GenIntVtx_x=0;
          GenIntVtx_y=0;
          GenIntVtx_z=0;

          GenHMass=125.6;
          GenHPt=20;
          GenHPhi=0;
          Lep1_Z1SIP=0; Lep2_Z1SIP=0; Lep3_Z1SIP=0; Lep4_Z1SIP=0; KalmanCandVtx_chi2=0;
          Lep1SIP=0; Lep2SIP=0; Lep3SIP=0; Lep4SIP=0;

          tc_extra->GetEntry(ev);

          float oldSIPdef = max(max(max(Lep1SIP, Lep2SIP), Lep3SIP), Lep4SIP);

          MC_weight = MC_weight_noxsec;
          if (iProcess==4) MC_weight *= MC_weight_QQZZEWK*MC_weight_QQBGGProper[1];
          if (iProcess==5) MC_weight *= MC_weight_Kfactor*MC_weight_QQBGGProper[0];

          int index_BS = 1;
          bool matchFound = false;
          for (int evBS=0; evBS < tBeam[index_BS]->GetEntries(); evBS++){
            tBeam[index_BS]->GetEntry(evBS);
            if (RunNumber == RunNumber_Ref && LumiNumber == LumiNumber_Ref){
              matchFound = true;
              break;
            }
            else continue;
          }
          if (!matchFound) cerr << "No beamspot matching possible!" << endl;

          CandBSVtx_x = KalmanCandVtx_x - BeamPosX;
          CandBSVtx_y = KalmanCandVtx_y - BeamPosY;
          CandBSVtx_z = KalmanCandVtx_z - BeamPosZ;
          Dxy_BS = (CandBSVtx_x*cos(ZZPhi) + CandBSVtx_y*sin(ZZPhi));
          Txy_BS = Dxy_BS*ZZMass / ZZPt;

          TrueCandVtx_x = GenIntVtx_x - GenPrimaryVtx_x;
          TrueCandVtx_y = GenIntVtx_y - GenPrimaryVtx_y;
          TrueCandVtx_z = GenIntVtx_z - GenPrimaryVtx_z;
          Dxy_true = (TrueCandVtx_x*cos(GenHPhi) + TrueCandVtx_y*sin(GenHPhi));
          Txy_true = Dxy_true*GenHMass / GenHPt;

          CandVtx_x = KalmanCandVtx_x - OfflinePrimaryVtx_x;
          CandVtx_y = KalmanCandVtx_y - OfflinePrimaryVtx_y;
          CandVtx_z = KalmanCandVtx_z - OfflinePrimaryVtx_z;
          Dxy = (CandVtx_x*cos(ZZPhi) + CandVtx_y*sin(ZZPhi));
          Txy = Dxy*ZZMass / ZZPt;

          sigmaDirectional_ratio = compute_subtractedPV(
            KalmanCandVtx_x,
            KalmanCandVtx_y,
            KalmanCandVtx_z,
            KalmanCandVtx_cov_xx,
            KalmanCandVtx_cov_xy,
            KalmanCandVtx_cov_xz,
            KalmanCandVtx_cov_yy,
            KalmanCandVtx_cov_yz,
            KalmanCandVtx_cov_zz,
            OfflinePrimaryVtx_x,
            OfflinePrimaryVtx_y,
            OfflinePrimaryVtx_z,
            OfflinePrimaryVtx_cov_xx,
            OfflinePrimaryVtx_cov_xy,
            OfflinePrimaryVtx_cov_xz,
            OfflinePrimaryVtx_cov_yy,
            OfflinePrimaryVtx_cov_yz,
            OfflinePrimaryVtx_cov_zz,
            SubtractedPrimaryVtx_x,
            SubtractedPrimaryVtx_y,
            SubtractedPrimaryVtx_z
            );

          if (!tc->GetBranchStatus("D_bkg")) D_bkg = p0plus_VAJHU*p0plus_m4l / (p0plus_VAJHU*p0plus_m4l + bkg_VAMCFM*bkg_m4l);

          KD[vnOPVTracks] = (OfflinePrimaryVtx_ndof+3.)/2.;
          KD[vsigmaDirRatio] = sigmaDirectional_ratio;
          KD[vdiffDotTransverse_DirSubPV] = ((SubtractedPrimaryVtx_x - GenPrimaryVtx_x)*cos(ZZPhi) + (SubtractedPrimaryVtx_y - GenPrimaryVtx_y)*sin(ZZPhi))*10000.;
          KD[vPVSIP] = oldSIPdef;
          KD[vLep1_Z1SIP] = Lep1_Z1SIP;
          KD[vLep2_Z1SIP] = Lep2_Z1SIP;
          KD[vLep3_Z1SIP] = Lep3_Z1SIP;
          KD[vLep4_Z1SIP] = Lep4_Z1SIP;
          KD[v4l_chi2] = KalmanCandVtx_chi2;
          KD[vTxy] = Txy*10000.;
          KD[vTxy_BS] = Txy_BS*10000.;
          KD[vDbkg] = D_bkg;


          for (int sc = 0; sc < newSIPCuts-1; sc++){
            if ((sc == 1 || sc >= 3) && (fabs(Lep1_Z1SIP) >= 4 || fabs(Lep2_Z1SIP) >= 4 || fabs(Lep3_Z1SIP) >= 5 || fabs(Lep4_Z1SIP) >= 5)) continue;
            if ((sc == 2 || sc >= 3) && KalmanCandVtx_chi2 >= 30) continue;
            if (sc == 4  && oldSIPdef>=4) continue;

//            if ((ZZMass >= 220 && ZZMass < 1600)) sum_offshell[sc] += MC_weight;
//            if ((ZZMass >= 80 + (iProcess==6 ? 20 : 0) && ZZMass < 100 + (iProcess==6 ? 20 : 0))) sum_80_100[sc] += MC_weight;
            if ((ZZMass >= 105.6 && ZZMass < 140.6)) sum_signal[sc] += MC_weight;

            if (ZZMass >= rangeZZMass[1] || ZZMass < rangeZZMass[0]) continue;
            float fillvar_2D[nVariables_2D][2] ={ { 0 } };
            for (int v = 0; v < nVariables; v++){
              float fillvar = KD[v];
              if (fillvar<hvars[v][sc]->GetXaxis()->GetXmin()) fillvar = hvars[v][sc]->GetXaxis()->GetXmin();
              if (fillvar>=hvars[v][sc]->GetXaxis()->GetXmax()) fillvar = hvars[v][sc]->GetXaxis()->GetXmax() - (hvars[v][sc]->GetXaxis()->GetBinWidth(1))*0.5;

              hvars[v][sc]->Fill(fillvar, MC_weight);

              for (int k = 0; k < nVariables_2D; k++){
                if (v == KD_pairing[k][0]) fillvar_2D[k][0] = fillvar;
                if (v == KD_pairing[k][1]) fillvar_2D[k][1] = fillvar;
              }
            }
            for (int v = 0; v < nVariables_2D; v++) hvars_2D[v][sc]->Fill(fillvar_2D[v][0], fillvar_2D[v][1], MC_weight);
          }
        }
      }
      cout << "Beginning scaling " << endl;
      for (int sc = 0; sc < newSIPCuts-1; sc++){
        double inputyield = 0, targetyield = 0, yieldScale = 1;
        if (iRegion == 0){ // SR
          inputyield = sum_signal[sc];
          cout << ratio_targetsignalyield[iProcess][sc+1][iRegion]/ratio_targetsignalyield[iProcess][0][iRegion] << endl;
          if (iProcess<=3) targetyield = yield_signal_higgs[EnergyIndex][folder]*ratio_targetsignalyield[iProcess][sc+1][iRegion]/ratio_targetsignalyield[0][0][iRegion];
          else if (iProcess==4) targetyield = yield_signal_qqzz[EnergyIndex][folder]*ratio_targetsignalyield[iProcess][sc+1][iRegion]/ratio_targetsignalyield[iProcess][0][iRegion];
          else if (iProcess==5) targetyield = yield_signal_ggzz[EnergyIndex][folder]*ratio_targetsignalyield[iProcess][sc+1][iRegion]/ratio_targetsignalyield[iProcess][0][iRegion];
          else if (iProcess==6) targetyield = yield_signal_zx[EnergyIndex][folder]*ratio_targetsignalyield[iProcess][sc+1][iRegion]/ratio_targetsignalyield[iProcess][0][iRegion];
        }
        else if (iRegion == 1){ // SB1
          inputyield = sum_80_100[sc];
          cout << ratio_targetsignalyield[iProcess][sc+1][0]/ratio_targetsignalyield[iProcess][0][0] << endl;
          if (iProcess<=3) targetyield = 0;
          else if (iProcess==4) targetyield = inputyield*luminosity[EnergyIndex];
          else if (iProcess==5) targetyield = yield_80_100_ggzz[EnergyIndex][folder]*ratio_targetsignalyield[iProcess][sc+1][0]/ratio_targetsignalyield[iProcess][0][0];
          else if (iProcess==6) targetyield = yield_80_100_zx[EnergyIndex][folder]*ratio_targetsignalyield[iProcess][sc+1][0]/ratio_targetsignalyield[iProcess][0][0];
          if (sc==3 && iProcess==6) cout << sum_80_100[sc] << '\t' << hvars[0][sc]->Integral() << '\t' << targetyield << '\t' << ratio_targetsignalyield[iProcess][sc+1][0]/ratio_targetsignalyield[iProcess][0][0] << endl;
        }
        else if (iRegion >= 2){ // SB2,3
          inputyield = sum_offshell[sc];
          cout << ratio_targetsignalyield[iProcess][sc+1][iRegion]/ratio_targetsignalyield[iProcess][0][iRegion] << endl;
          if (iProcess<=3) targetyield = 0;
          else if (iProcess==4) targetyield = inputyield*luminosity[EnergyIndex];
          else if (iProcess==5) targetyield = yield_offshell_ggzz[EnergyIndex][folder]*ratio_targetsignalyield[iProcess][sc+1][iRegion]/ratio_targetsignalyield[iProcess][0][iRegion];
          else if (iProcess==6) targetyield = yield_offshell_zx[EnergyIndex][folder]*ratio_targetsignalyield[iProcess][sc+1][iRegion]/ratio_targetsignalyield[iProcess][0][iRegion];
        }
        if (iProcess<7){
          yieldScale = (inputyield!=0 ? (targetyield/inputyield) : 0);
          for (int v = 0; v < nVariables; v++){
            hvars[v][sc]->Scale(yieldScale);
          }
          for (int v = 0; v < nVariables_2D; v++){
            hvars_2D[v][sc]->Scale(yieldScale);
          }
        }
      }
      // 1D recording
      for (int sc = 0; sc < newSIPCuts-1; sc++){
        for (int v = 0; v < nVariables; v++){
          if (v >= vLep1_Z1SIP && v <= vLep4_Z1SIP && (sc == 1 || sc == 3)) continue;
          if (v == v4l_chi2 && sc >= 2) continue;
          if ((v == vLep1_Z1SIP || v == vLep3_Z1SIP) && (sc != 1 && sc != 3)){
            hvars[v][sc]->Add(hvars[v + 1][sc], 1);
            hvars[v][sc]->Scale(0.5);
            if (v == vLep1_Z1SIP) hvars[v][sc]->SetXTitle("#it{SIP}^{Z1}_{1,2} (Signed)");
            else if (v == vLep3_Z1SIP) hvars[v][sc]->SetXTitle("#it{SIP}^{Z1}_{3,4} (Signed)");
          }
          else if ((v == vLep2_Z1SIP || v == vLep4_Z1SIP) && (sc != 1 && sc != 3)) continue;
          if (v == vdiffDotTransverse_DirSubPV && iProcess >= 6) continue;
          foutput->WriteTObject(hvars[v][sc]);
        }
      }
      // 1D plotting
      for (int v = 0; v < nVariables; v++){
        bool canSkip[newSIPCuts] ={ 0 };
        int nDrawn=0;

        for (int sc = 0; sc < newSIPCuts-1; sc++){
          if (v >= vLep1_Z1SIP && v <= vLep4_Z1SIP && (sc == 1 || sc == 3)) canSkip[sc]=true;
          if (v == v4l_chi2 && sc >= 2) canSkip[sc]=true;
          if ((v == vLep2_Z1SIP || v == vLep4_Z1SIP) && (sc != 1 && sc != 3)) canSkip[sc]=true;
          if (v == vdiffDotTransverse_DirSubPV && iProcess >= 6) canSkip[sc]=true;
          if (!canSkip[sc]) nDrawn++;
        }
        if (canSkip[0]) continue;

        double maxplot = hvars[v][0]->GetMaximum();
        hvars[v][0]->GetYaxis()->SetRangeUser(0, maxplot*1.65);

        TPaveText* pt = new TPaveText(0.15, 0.93, 0.85, 1, "brNDC");
        pt->SetBorderSize(0);
        pt->SetFillStyle(0);
        pt->SetTextAlign(12);
        pt->SetTextFont(42);
        pt->SetTextSize(0.045);
        TText* text = pt->AddText(0.025, 0.45, "#font[61]{CMS}");
        text->SetTextSize(0.044);
        text = pt->AddText(0.165, 0.42, "#font[52]{Unpublished}");
        text->SetTextSize(0.0315);
        TString cErgTev = Form("#font[42]{%.1f fb^{-1} (%i TeV)}", luminosity[EnergyIndex], erg_tev);
        if (folder == 0) cErgTev.Prepend("4#mu ");
        if (folder == 1) cErgTev.Prepend("4e ");
        if (folder == 2) cErgTev.Prepend("2e+2#mu ");
        text = pt->AddText((folder<2 ? 0.790 : 0.725), 0.45, cErgTev);
        text->SetTextSize(0.0315);

        hvars[v][0]->SetLineColor(kBlack);
        hvars[v][0]->SetLineWidth(2);
        hvars[v][1]->SetLineColor(kBlue);
        hvars[v][1]->SetLineWidth(2);
        hvars[v][2]->SetLineColor(kGreen+2);
        hvars[v][2]->SetLineWidth(2);
        hvars[v][3]->SetLineColor(kRed);
        hvars[v][3]->SetLineWidth(2);

        foutput->cd();
        TString appendName;

        if (iRegion == 0) appendName = Form("_%s_SignalRegion_%s_%iTeV", processName[iProcess].c_str(), user_folder[folder], erg_tev);
        else appendName = Form("_%s_SidebandRegion%i_%s_%iTeV", processName[iProcess].c_str(), iRegion, user_folder[folder], erg_tev);
        TString canvasname = "cCompareSIPCuts_";
        canvasname.Append(varNames[v]);
        canvasname.Append(appendName);
        TCanvas* cc = new TCanvas(canvasname, "", 8, 30, 800, 800);
        cc->cd();
        gStyle->SetOptStat(0);
        cc->SetFillColor(0);
        cc->SetBorderMode(0);
        cc->SetBorderSize(2);
        cc->SetTickx(1);
        cc->SetTicky(1);
        cc->SetLeftMargin(0.17);
        cc->SetRightMargin(0.05);
        cc->SetTopMargin(0.07);
        cc->SetBottomMargin(0.13);
        cc->SetFrameFillStyle(0);
        cc->SetFrameBorderMode(0);
        cc->SetFrameFillStyle(0);
        cc->SetFrameBorderMode(0);

        TLegend *ll;
        if (nDrawn>2) ll = new TLegend(0.20, 0.60, 0.45, 0.90);
        else if (canSkip[1]) ll = new TLegend(0.20, 0.735, 0.45, 0.90); // if 4l_chi cut is plotted
        else if (canSkip[2]) ll = new TLegend(0.20, 0.735, 0.45, 0.90); // if Z!SIP cut is plotted
        ll->SetBorderSize(0);
        ll->SetTextFont(42);
        ll->SetTextSize(0.03);
        ll->SetLineColor(1);
        ll->SetLineStyle(1);
        ll->SetLineWidth(1);
        ll->SetFillColor(0);
        ll->SetFillStyle(0);

        for (int sc = 0; sc < newSIPCuts-1; sc++){
          if (canSkip[sc]) continue;
          ll->AddEntry(hvars[v][sc], ctitle_SIPCut[sc+1], "l");
          if (sc==0) hvars[v][sc]->Draw("hist");
          else hvars[v][sc]->Draw("histsame");
        }
        ll->Draw("same");
        pt->Draw();

        TPaveText *pt10 = new TPaveText(0.67, 0.86, 0.90, 0.90, "brNDC");
        pt10->SetBorderSize(0);
        pt10->SetTextAlign(12);
        pt10->SetTextSize(0.03);
        pt10->SetFillStyle(0);
        pt10->SetTextFont(42);
        TText* text10;
        if (iRegion==0) text10 = pt10->AddText(0.01, 0.01, "Signal Region");
        else text10 = pt10->AddText(0.01, 0.01, Form("Sideband %i", iRegion));
        pt10->Draw();

        TPaveText *pt20 = new TPaveText(0.67, 0.82, 0.90, 0.86, "brNDC");
        pt20->SetBorderSize(0);
        pt20->SetTextAlign(12);
        pt20->SetTextSize(0.03);
        pt20->SetFillStyle(0);
        pt20->SetTextFont(42);
        TText* text20 = pt20->AddText(0.01, 0.01, processTitle[iProcess]);
        pt20->Draw();

        cc->RedrawAxis();
        cc->Modified();
        cc->Update();
        foutput->WriteTObject(cc);
        canvasname.Prepend(plotDir_1D);
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
        cc->SaveAs(canvasname_pdf);
        cc->SaveAs(canvasname_eps);
        cc->SaveAs(canvasname_png);
        cc->SaveAs(canvasname_root);
        cc->SaveAs(canvasname_c);
        delete ll;
        cc->Close();
        delete pt;
      }
      // 2D recording + plotting
      for (int sc = 0; sc < newSIPCuts-1; sc++){
        for (int v = 0; v < nVariables_2D; v++){
          bool notRecordable = false;
          for (int m = 0; m < 2; m++){
            if (KD_pairing[v][m] >= vLep1_Z1SIP && KD_pairing[v][m] <= vLep4_Z1SIP && (sc == 1 || sc == 3)) notRecordable = true;
            if (KD_pairing[v][m] == v4l_chi2 && sc >= 2) notRecordable = true;
            if (KD_pairing[v][m] == vdiffDotTransverse_DirSubPV && iProcess >= 6) notRecordable = true;
          }
          if (notRecordable) continue;
          if ((KD_pairing[v][0] == vLep1_Z1SIP || KD_pairing[v][0] == vLep3_Z1SIP) && (KD_pairing[v][1] == vPVSIP || KD_pairing[v][1] == v4l_chi2) && (sc != 1 && sc != 3)){
            hvars_2D[v][sc]->Add(hvars_2D[v + 1][sc], 1);
            hvars_2D[v][sc]->Scale(0.5);
            if (KD_pairing[v][0] == vLep1_Z1SIP) hvars_2D[v][sc]->SetXTitle("#it{SIP}^{Z1}_{1,2} (Signed)");
            else if (KD_pairing[v][0] == vLep3_Z1SIP) hvars_2D[v][sc]->SetXTitle("#it{SIP}^{Z1}_{3,4} (Signed)");
          }
          else if ((KD_pairing[v][0] == vLep2_Z1SIP || KD_pairing[v][0] == vLep4_Z1SIP) && (KD_pairing[v][1] == vPVSIP || KD_pairing[v][1] == v4l_chi2) && (sc != 1 && sc != 3)) continue;
          foutput->WriteTObject(hvars_2D[v][sc]);
          string defaultName = hvars_2D[v][sc]->GetName();
          if (iRegion == 0) hvars_2D[v][sc]->SetName(Form("%s_%s_SignalRegion", processName[iProcess].c_str(), defaultName.c_str()));
          else hvars_2D[v][sc]->SetName(Form("%s_%s_SidebandRegion%i", processName[iProcess].c_str(), defaultName.c_str(), iRegion));
          generic_Histo2DPlotter(foutput, plotDir_2D, hvars_2D[v][sc], erg_tev, luminosity[EnergyIndex], folder, (iProcess==6));
          hvars_2D[v][sc]->SetName(defaultName.c_str());
        }
      }
      if (tc_extra!=0) delete tc_extra;
      delete tc;
      foutput->Close();
      for (int sc = 0; sc < newSIPCuts+1; sc++) delete hRatio[sc];
      fSSinput->Close();
    }

    for (int bs=0; bs<2; bs++) delete tBeam[bs];
  }
}

void extract_SIPSyst(){
  const int nProcesses = 4;
  string processName[nProcesses] ={
    "qqZZ", "ggZZ", "CR", "data"
  };
  TH1F* hvars[nProcesses][3*2][newSIPCuts] ={ { { 0 } } };

  TString cinput_common = user_dir_hep + "Analysis/Auxiliary/";

  for (int iProcess=0; iProcess<nProcesses; iProcess++){
    for (int erg_tev = 7; erg_tev < 9; erg_tev++){
      TString erg_dir;
      erg_dir.Form("LHC_%iTeV", erg_tev);
      TString comstring;
      comstring.Form("%iTeV", erg_tev);
      int EnergyIndex = 0; if (erg_tev == 8) EnergyIndex = 1;

      for (int folder = 0; folder < 3; folder++){
        int v = 3*EnergyIndex + folder;

        TString INPUT_NAME = "LifetimeKD_SIPVariables_";
        INPUT_NAME.Append(Form("%s", processName[iProcess].c_str()));
        INPUT_NAME.Append(Form("_%s_", user_folder[folder]));
        INPUT_NAME.Append(comstring);
        bool combineSidebands=true;
        INPUT_NAME.Append(Form("_SB%i", 1));
        INPUT_NAME.Append(".root");

        TString cinput = cinput_common + INPUT_NAME;
        TFile* finput = new TFile(cinput, "read");

        cout << "Read " << cinput << endl;
        if (finput!=0){
          if (!finput->IsZombie()){
            TString hname_core = "h_";
            for (int sc = 0; sc < newSIPCuts; sc++){
              TString hname = hname_core;
              hname.Append(varNames[vPVSIP]);
              hname.Append("_");
              hname.Append(cutNames[sc + 1]);
              TH1F* htemp = (TH1F*)finput->Get(hname);
              if (htemp!=0){
                htemp->Sumw2();
                //                  cout << htemp->GetName() << endl;
                int h_index = iProcess;
                TString tempname = htemp->GetName();
                tempname.Append(processName[iProcess]);
                htemp->SetName(tempname);

                if (hvars[h_index][v][sc]!=0) hvars[h_index][v][sc]->Add(htemp);
                else{
                  gROOT->cd();
                  hvars[h_index][v][sc] = (TH1F*)htemp->Clone(Form("%s_clone", htemp->GetName()));
                }
                delete htemp;
              }
            }
            finput->Close();
          }
        }
        if (combineSidebands){
          INPUT_NAME = "LifetimeKD_SIPVariables_";
          INPUT_NAME.Append(Form("%s", processName[iProcess].c_str()));
          INPUT_NAME.Append(Form("_%s_", user_folder[folder]));
          INPUT_NAME.Append(comstring);
          INPUT_NAME.Append(Form("_SB%i", 3));
          INPUT_NAME.Append(".root");

          cinput = cinput_common + INPUT_NAME;
          finput = new TFile(cinput, "read");

          cout << "Read " << cinput << endl;
          if (finput!=0){
            if (!finput->IsZombie()){
              TString hname_core = "h_";
              for (int sc = 0; sc < newSIPCuts; sc++){
                TString hname = hname_core;
                hname.Append(varNames[vPVSIP]);
                hname.Append("_");
                hname.Append(cutNames[sc + 1]);
                TH1F* htemp = (TH1F*)finput->Get(hname);
                if (htemp!=0){
                  //                  cout << htemp->GetName() << endl;
                  int h_index = iProcess;
                  TString tempname = htemp->GetName();
                  tempname.Append(processName[iProcess]);
                  htemp->SetName(tempname);
                  if (hvars[h_index][v][sc]!=0) hvars[h_index][v][sc]->Add(htemp);
                  else{
                    gROOT->cd();
                    hvars[h_index][v][sc] = (TH1F*)htemp->Clone(Form("%s_clone", htemp->GetName()));
                  }
                  delete htemp;
                }
              }
              finput->Close();
            }
          }
        }
      }
    }
  }

  gROOT->cd();

  cout << "sqrt(s)\tChannel\tProcess\tNo cut\tSIP_3D cut\tNew cut\tNo cut error\tSIP_3D cut error\tNew cut error" << endl;
  for (int erg_tev = 7; erg_tev < 9; erg_tev++){
    int EnergyIndex = 0; if (erg_tev == 8) EnergyIndex = 1;
    for (int folder = 0; folder < 3; folder++){
      int v = 3*EnergyIndex + folder;
      for (int hh=0; hh<nProcesses; hh++){
        int bin_sip3d_4 = hvars[hh][v][0]->GetXaxis()->FindBin(4.001);
        if (bin_sip3d_4!=5) cout << "bin_sip3d_4 = " << bin_sip3d_4 << endl;
        double integral_nocut=0, error_nocut=0;
        double integral_newcut=0, error_newcut=0;
        double integral_sip3dcut=0, error_sip3dcut=0;
        integral_nocut = hvars[hh][v][0]->IntegralAndError(1, hvars[hh][v][0]->GetNbinsX(), error_nocut);
        integral_newcut = hvars[hh][v][3]->IntegralAndError(1, hvars[hh][v][3]->GetNbinsX(), error_newcut);
        integral_sip3dcut = hvars[hh][v][0]->IntegralAndError(1, bin_sip3d_4-1, error_sip3dcut);

        cout << erg_tev << '\t' << user_folder[folder] << '\t' << processName[hh] << '\t'
          << integral_nocut << '\t'
          << integral_sip3dcut << '\t'
          << integral_newcut << '\t'
          << error_nocut << '\t'
          << error_sip3dcut << '\t'
          << error_newcut
          << endl;
      }
    }
  }

  for (int hh=0; hh<nProcesses; hh++){
    for (int erg_tev = 7; erg_tev < 9; erg_tev++){
      int EnergyIndex = 0; if (erg_tev == 8) EnergyIndex = 1;
      for (int folder = 0; folder < 3; folder++){
        int v = 3*EnergyIndex + folder;
        for (int sc = 0; sc < newSIPCuts; sc++){
          if (hvars[hh][v][sc]!=0) delete hvars[hh][v][sc];
        }
      }
    }
  }
}

void combineAllChannels_perProcess_ratio(TH1F* hvars[newSIPCuts-1], bool canSkip[newSIPCuts-1], int nDrawn, int iProcess, int v, int iRegion, TString plotDir_1D);

void combineAllChannels_perProcess(TH1F* hvars[newSIPCuts-1], int iProcess, int v, int iRegion, TString plotDir){
  bool canSkip[newSIPCuts-1] ={ 0 };
  int nDrawn=0;
  for (int sc=0; sc<newSIPCuts-1; sc++){
    if (v >= vLep1_Z1SIP && v <= vLep4_Z1SIP && (sc == 1 || sc == 3)) canSkip[sc]=true;
    if (v == v4l_chi2 && sc >= 2) canSkip[sc]=true;
    if ((v == vLep2_Z1SIP || v == vLep4_Z1SIP) && (sc != 1 && sc != 3)) canSkip[sc]=true;
    if (v == vdiffDotTransverse_DirSubPV && iProcess >= 6) canSkip[sc]=true;
    if (hvars[sc]==0) canSkip[sc]=true;
    if (!canSkip[sc]) nDrawn++;
    else continue;
    cout << "Obtained " << hvars[sc]->GetName() << endl;
  }
  if (canSkip[0] || nDrawn<=1) return;
  else cout << "Proceeding... ";

  TString plotDir_1D=plotDir;
  plotDir_1D.Append(Form("%s/", processName[iProcess].c_str()));
  plotDir_1D.Append("1D/");
  TString mkdirCommand = "mkdir -p ";
  TString mkdirCommand_1D = mkdirCommand;
  mkdirCommand_1D.Append(plotDir_1D);
  gSystem->Exec(mkdirCommand_1D);

  double maxplot = hvars[0]->GetMaximum();
  hvars[0]->GetYaxis()->SetRangeUser(0, maxplot*1.65);

  TPaveText* pt = new TPaveText(0.15, 0.93, 0.85, 1, "brNDC");
  pt->SetBorderSize(0);
  pt->SetFillStyle(0);
  pt->SetTextAlign(12);
  pt->SetTextFont(42);
  pt->SetTextSize(0.045);
  TText* text = pt->AddText(0.025, 0.45, "#font[61]{CMS}");
  text->SetTextSize(0.044);
  text = pt->AddText(0.165, 0.42, "#font[52]{Unpublished}");
  text->SetTextSize(0.0315);
  TString cErgTev = "#font[42]{19.7 fb^{-1} (8 TeV) + 5.1 fb^{-1} (7 TeV)}";
  text = pt->AddText(0.537, 0.45, cErgTev);
  text->SetTextSize(0.0315);

  hvars[0]->SetLineColor(kBlack);
  hvars[0]->SetLineWidth(2);
  if (!canSkip[1]){
    hvars[1]->SetLineColor(kBlue);
    hvars[1]->SetLineWidth(2);
  }
  if (!canSkip[2]){
    hvars[2]->SetLineColor(kGreen+2);
    hvars[2]->SetLineWidth(2);
  }
  if (!canSkip[3]){
    hvars[3]->SetLineColor(kRed);
    hvars[3]->SetLineWidth(2);
  }
  TString appendName;
  if (iRegion == 0) appendName = Form("_%s_SignalRegion_AllTeV", processName[iProcess].c_str());
  else if (iRegion<4) appendName = Form("_%s_SidebandRegion%i_AllTeV", processName[iProcess].c_str(), iRegion);
  else if (iRegion==4) appendName = Form("_%s_SidebandRegion13_AllTeV", processName[iProcess].c_str());
  TString canvasname = "cCompareSIPCuts_";
  canvasname.Append(varNames[v]);
  canvasname.Append(appendName);
  TCanvas* cc = new TCanvas(canvasname, "", 8, 30, 800, 800);
  cc->cd();
  gStyle->SetOptStat(0);
  cc->SetFillColor(0);
  cc->SetBorderMode(0);
  cc->SetBorderSize(2);
  cc->SetTickx(1);
  cc->SetTicky(1);
  cc->SetLeftMargin(0.17);
  cc->SetRightMargin(0.05);
  cc->SetTopMargin(0.07);
  cc->SetBottomMargin(0.13);
  cc->SetFrameFillStyle(0);
  cc->SetFrameBorderMode(0);
  cc->SetFrameFillStyle(0);
  cc->SetFrameBorderMode(0);

  TLegend *ll;
  if (nDrawn>2) ll = new TLegend(0.20, 0.60, 0.45, 0.90);
  else if (canSkip[1]) ll = new TLegend(0.20, 0.735, 0.45, 0.90); // if 4l_chi cut is plotted
  else if (canSkip[2]) ll = new TLegend(0.20, 0.735, 0.45, 0.90); // if Z1-SIP cut is plotted
  ll->SetBorderSize(0);
  ll->SetTextFont(42);
  ll->SetTextSize(0.03);
  ll->SetLineColor(1);
  ll->SetLineStyle(1);
  ll->SetLineWidth(1);
  ll->SetFillColor(0);
  ll->SetFillStyle(0);

  for (int sc = 0; sc < newSIPCuts-1; sc++){
    if (canSkip[sc]) continue;
    else cout << "Plotting " << sc;
    ll->AddEntry(hvars[sc], ctitle_SIPCut[sc+1], "l");
    if (sc==0) hvars[sc]->Draw("hist");
    else hvars[sc]->Draw("histsame");
  }
  ll->Draw("same");
  pt->Draw();

  TPaveText *pt10 = new TPaveText(0.67, 0.86, 0.90, 0.90, "brNDC");
  pt10->SetBorderSize(0);
  pt10->SetTextAlign(12);
  pt10->SetTextSize(0.03);
  pt10->SetFillStyle(0);
  pt10->SetTextFont(42);
  TText* text10;
  if (iRegion==0) text10 = pt10->AddText(0.01, 0.01, "Signal Region");
  else if (iRegion<4) text10 = pt10->AddText(0.01, 0.01, Form("Sideband %i", iRegion));
  else if (iRegion==4) text10 = pt10->AddText(0.01, 0.01, "Sideband 1&3");
  pt10->Draw();

  TPaveText *pt20 = new TPaveText(0.67, 0.82, 0.90, 0.86, "brNDC");
  pt20->SetBorderSize(0);
  pt20->SetTextAlign(12);
  pt20->SetTextSize(0.03);
  pt20->SetFillStyle(0);
  pt20->SetTextFont(42);
  TText* text20 = pt20->AddText(0.01, 0.01, processTitle[iProcess]);
  pt20->Draw();

  cc->RedrawAxis();
  cc->Modified();
  cc->Update();
  canvasname.Prepend(plotDir_1D);
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
  cc->SaveAs(canvasname_pdf);
  cc->SaveAs(canvasname_eps);
  cc->SaveAs(canvasname_png);
  cc->SaveAs(canvasname_root);
  cc->SaveAs(canvasname_c);
  delete pt20;
  delete pt10;
  delete ll;
  cc->Close();
  delete pt;

  combineAllChannels_perProcess_ratio(hvars, canSkip, nDrawn, iProcess, v, iRegion, plotDir_1D);
}

void combineAllChannels_perProcess_ratio(TH1F* hvars[newSIPCuts-1], bool canSkip[newSIPCuts-1], int nDrawn, int iProcess, int v, int iRegion, TString plotDir_1D){
  double maxplot = 0, minplot=99;
  for (int sc=1; sc<newSIPCuts-1; sc++){
    if (canSkip[sc]) continue;
    for (int bin=1; bin<=hvars[0]->GetNbinsX(); bin++){
      double all_val = hvars[0]->GetBinContent(bin);
      double all_err = pow(hvars[0]->GetBinError(bin), 2);
      double val = hvars[sc]->GetBinContent(bin);
      double err = pow(hvars[sc]->GetBinError(bin), 2);

      double ratio = 0;
      if(all_val>0) ratio = val/all_val;
      double ratio_error = 0;
      if(all_val>0) ratio_error = sqrt(max(pow(val, 2)*(all_err-err) + pow(all_val-val, 2)*err,0.))/pow(all_val, 2);
      hvars[sc]->SetBinContent(bin, ratio);
      hvars[sc]->SetBinError(bin, ratio_error);

      if (ratio+ratio_error>maxplot) maxplot = ratio+ratio_error;
      if (ratio-ratio_error<minplot) minplot = ratio-ratio_error;
    }
  }
  for (int sc=1; sc<newSIPCuts-1; sc++){
    if (canSkip[sc]) continue;
    hvars[sc]->GetYaxis()->SetRangeUser(minplot*0.70, maxplot*1.65);
  }

  TPaveText* pt = new TPaveText(0.15, 0.93, 0.85, 1, "brNDC");
  pt->SetBorderSize(0);
  pt->SetFillStyle(0);
  pt->SetTextAlign(12);
  pt->SetTextFont(42);
  pt->SetTextSize(0.045);
  TText* text = pt->AddText(0.025, 0.45, "#font[61]{CMS}");
  text->SetTextSize(0.044);
  text = pt->AddText(0.165, 0.42, "#font[52]{Unpublished}");
  text->SetTextSize(0.0315);
  TString cErgTev = "#font[42]{19.7 fb^{-1} (8 TeV) + 5.1 fb^{-1} (7 TeV)}";
  text = pt->AddText(0.537, 0.45, cErgTev);
  text->SetTextSize(0.0315);

  TString appendName;
  if (iRegion == 0) appendName = Form("_%s_SignalRegion_AllTeV", processName[iProcess].c_str());
  else if (iRegion<4) appendName = Form("_%s_SidebandRegion%i_AllTeV", processName[iProcess].c_str(), iRegion);
  else if (iRegion==4) appendName = Form("_%s_SidebandRegion13_AllTeV", processName[iProcess].c_str());
  TString canvasname = "cCompareSIPCuts_Ratio_";
  canvasname.Append(varNames[v]);
  canvasname.Append(appendName);
  TCanvas* cc = new TCanvas(canvasname, "", 8, 30, 800, 800);
  cc->cd();
  gStyle->SetOptStat(0);
  cc->SetFillColor(0);
  cc->SetBorderMode(0);
  cc->SetBorderSize(2);
  cc->SetTickx(1);
  cc->SetTicky(1);
  cc->SetLeftMargin(0.17);
  cc->SetRightMargin(0.05);
  cc->SetTopMargin(0.07);
  cc->SetBottomMargin(0.13);
  cc->SetFrameFillStyle(0);
  cc->SetFrameBorderMode(0);
  cc->SetFrameFillStyle(0);
  cc->SetFrameBorderMode(0);

  TLegend *ll;
  if (nDrawn>2) ll = new TLegend(0.20, 0.60, 0.45, 0.90);
  else if (canSkip[1]) ll = new TLegend(0.20, 0.735, 0.45, 0.90); // if 4l_chi cut is plotted
  else if (canSkip[2]) ll = new TLegend(0.20, 0.735, 0.45, 0.90); // if Z1-SIP cut is plotted
  ll->SetBorderSize(0);
  ll->SetTextFont(42);
  ll->SetTextSize(0.03);
  ll->SetLineColor(1);
  ll->SetLineStyle(1);
  ll->SetLineWidth(1);
  ll->SetFillColor(0);
  ll->SetFillStyle(0);

  bool initialized=false;
  for (int sc = 1; sc < newSIPCuts-1; sc++){
    if (canSkip[sc]) continue;

    hvars[sc]->GetYaxis()->SetTitle("Ratio");
    hvars[sc]->SetMarkerStyle(1);
    hvars[sc]->SetMarkerSize(1);
    hvars[sc]->SetMarkerColor(hvars[sc]->GetLineColor());
//    hvars[sc]->SetFillStyle(3005);
//    hvars[sc]->SetFillColor(hvars[sc]->GetLineColor());
    if (!initialized){
      hvars[sc]->Draw("histe1");
      initialized=true;
    }
    else hvars[sc]->Draw("histe1same");
    ll->AddEntry(hvars[sc], ctitle_SIPCut[sc+1], "l");
  }
  ll->Draw("same");
  pt->Draw();

  TPaveText *pt10 = new TPaveText(0.67, 0.86, 0.90, 0.90, "brNDC");
  pt10->SetBorderSize(0);
  pt10->SetTextAlign(12);
  pt10->SetTextSize(0.03);
  pt10->SetFillStyle(0);
  pt10->SetTextFont(42);
  TText* text10;
  if (iRegion==0) text10 = pt10->AddText(0.01, 0.01, "Signal Region");
  else if (iRegion<4) text10 = pt10->AddText(0.01, 0.01, Form("Sideband %i", iRegion));
  else if (iRegion==4) text10 = pt10->AddText(0.01, 0.01, "Sideband 1&3");
  pt10->Draw();

  TPaveText *pt20 = new TPaveText(0.67, 0.82, 0.90, 0.86, "brNDC");
  pt20->SetBorderSize(0);
  pt20->SetTextAlign(12);
  pt20->SetTextSize(0.03);
  pt20->SetFillStyle(0);
  pt20->SetTextFont(42);
  TText* text20 = pt20->AddText(0.01, 0.01, processTitle[iProcess]);
  pt20->Draw();

  cc->RedrawAxis();
  cc->Modified();
  cc->Update();
  canvasname.Prepend(plotDir_1D);
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
  cc->SaveAs(canvasname_pdf);
  cc->SaveAs(canvasname_eps);
  cc->SaveAs(canvasname_png);
  cc->SaveAs(canvasname_root);
  cc->SaveAs(canvasname_c);
  delete pt20;
  delete pt10;
  delete ll;
  cc->Close();
  delete pt;
}

void plotCumulative_SIPVariables(int iRegion=0){
  gROOT->ProcessLine(".x tdrstyle.cc");
  float rangeZZMass[2] ={ systZZMass_range[iRegion][0], systZZMass_range[iRegion][1] };
  TH1F* hvars[8][nVariables][newSIPCuts] ={ { { 0 } } };


  TString cinput_common = user_dir_hep + "Analysis/Auxiliary/";

  for (int iProcess=0; iProcess<8; iProcess++){
    if (iRegion!=0 && iProcess<=3)continue;

    for (int erg_tev = 7; erg_tev < 9; erg_tev++){
      TString erg_dir;
      erg_dir.Form("LHC_%iTeV", erg_tev);
      TString comstring;
      comstring.Form("%iTeV", erg_tev);
      int EnergyIndex = 0; if (erg_tev == 8) EnergyIndex = 1;

      for (int folder = 0; folder < 3; folder++){
        TString INPUT_NAME = "LifetimeKD_SIPVariables_";
        INPUT_NAME.Append(Form("%s", processName[iProcess].c_str()));
        INPUT_NAME.Append(Form("_%s_", user_folder[folder]));
        INPUT_NAME.Append(comstring);
        bool combineSidebands=false;
        if (iRegion == 0) INPUT_NAME.Append("_SR");
        else if (iRegion <= 3) INPUT_NAME.Append(Form("_SB%i", iRegion));
        else if (iRegion==4){ INPUT_NAME.Append(Form("_SB%i", 1)); combineSidebands=true; }
        else assert(0);
        INPUT_NAME.Append(".root");

        TString cinput = cinput_common + INPUT_NAME;
        TFile* finput = new TFile(cinput, "read");

        cout << "Read " << cinput << endl;
        if (finput!=0){
          if (!finput->IsZombie()){
            TString hname_core = "h_";
            for (int v = 0; v < nVariables; v++){
              for (int sc = 0; sc < newSIPCuts; sc++){
                TString hname = hname_core;
                hname.Append(varNames[v]);
                hname.Append("_");
                hname.Append(cutNames[sc + 1]);
                TH1F* htemp = (TH1F*)finput->Get(hname);
                if (htemp!=0){
                  //                  cout << htemp->GetName() << endl;
                  int h_index = iProcess;
                  TString tempname = htemp->GetName();
                  tempname.Append(processName[iProcess]);
                  htemp->SetName(tempname);
                  if (h_index==6 && v==0 && sc==3) cout << htemp->Integral() << endl;

                  if (hvars[h_index][v][sc]!=0) hvars[h_index][v][sc]->Add(htemp);
                  else{
                    gROOT->cd();
                    hvars[h_index][v][sc] = (TH1F*)htemp->Clone(Form("%s_clone", htemp->GetName()));
                  }
                  delete htemp;
                }
              }
            }
            finput->Close();
          }
        }
        if (combineSidebands){
          INPUT_NAME = "LifetimeKD_SIPVariables_";
          INPUT_NAME.Append(Form("%s", processName[iProcess].c_str()));
          INPUT_NAME.Append(Form("_%s_", user_folder[folder]));
          INPUT_NAME.Append(comstring);
          INPUT_NAME.Append(Form("_SB%i", 3));
          INPUT_NAME.Append(".root");

          cinput = cinput_common + INPUT_NAME;
          finput = new TFile(cinput, "read");

          cout << "Read " << cinput << endl;
          if (finput!=0){
            if (!finput->IsZombie()){
              TString hname_core = "h_";
              for (int v = 0; v < nVariables; v++){
                for (int sc = 0; sc < newSIPCuts; sc++){
                  TString hname = hname_core;
                  hname.Append(varNames[v]);
                  hname.Append("_");
                  hname.Append(cutNames[sc + 1]);
                  TH1F* htemp = (TH1F*)finput->Get(hname);
                  if (htemp!=0){
                    //                  cout << htemp->GetName() << endl;
                    int h_index = iProcess;
                    TString tempname = htemp->GetName();
                    tempname.Append(processName[iProcess]);
                    htemp->SetName(tempname);
                    if (hvars[h_index][v][sc]!=0) hvars[h_index][v][sc]->Add(htemp);
                    else{
                      gROOT->cd();
                      hvars[h_index][v][sc] = (TH1F*)htemp->Clone(Form("%s_clone", htemp->GetName()));
                    }
                    delete htemp;
                  }
                }
              }
              finput->Close();
            }
          }
        }

      }
    }
  }
  TString plotDir = cinput_common;
  plotDir.Append("../Plots/SIP/");
  if (iRegion==0) plotDir.Append("SignalRegion/");
  else if (iRegion<4) plotDir.Append(Form("SidebandRegion%i/", iRegion));
  else if (iRegion==4) plotDir.Append("SidebandRegion13/");
  TString plotDir_1D=plotDir;
  plotDir_1D.Append("Cumulative/");
  plotDir_1D.Append("1D/");
  TString mkdirCommand = "mkdir -p ";
  TString mkdirCommand_1D = mkdirCommand;
  mkdirCommand_1D.Append(plotDir_1D);
  gSystem->Exec(mkdirCommand_1D);

  for (int v = 0; v < nVariables; v++){
    for (int pp=0; pp<8; pp++){
      TH1F* hcvars[newSIPCuts-1]={ 0 };
      bool processDNE=true;
      int nPlotted=0;
      for (int sc=0; sc<newSIPCuts-1; sc++){
        bool canSkip = false;
        if (v >= vLep1_Z1SIP && v <= vLep4_Z1SIP && (sc == 1 || sc == 3)) canSkip=true;
        if (v == v4l_chi2 && sc >= 2) canSkip=true;
        if ((v == vLep2_Z1SIP || v == vLep4_Z1SIP) && (sc != 1 && sc != 3)) canSkip=true;
        if (v == vdiffDotTransverse_DirSubPV) canSkip=true;
        if (v == vdiffDotTransverse_DirSubPV && pp >= 6) canSkip=true;
        if (canSkip) continue;
        if (hvars[pp][v][sc]!=0){
          hcvars[sc]=(TH1F*)hvars[pp][v][sc]->Clone(Form("%s_copy", hvars[pp][v][sc]->GetName()));
          cout << "Admitting " << hcvars[sc]->GetName() << endl;
          processDNE=false;
          nPlotted++;
        }
      }
      if (!processDNE && nPlotted>1) combineAllChannels_perProcess(hcvars, pp, v, iRegion, plotDir);
      for (int sc=0; sc<newSIPCuts-1; sc++){
        bool canSkip = false;
        if (v >= vLep1_Z1SIP && v <= vLep4_Z1SIP && (sc == 1 || sc == 3)) canSkip=true;
        if (v == v4l_chi2 && sc >= 2) canSkip=true;
        if ((v == vLep2_Z1SIP || v == vLep4_Z1SIP) && (sc != 1 && sc != 3)) canSkip=true;
        if (v == vdiffDotTransverse_DirSubPV) canSkip=true;
        if (v == vdiffDotTransverse_DirSubPV && pp >= 6) canSkip=true;
        if (canSkip) continue;
        if (hcvars[sc]!=0){
          cout << "Deleting " << hcvars[sc]->GetName() << endl;
          delete hcvars[sc];
        }
      }
    }

    for (int sc = 0; sc < newSIPCuts-1; sc++){
      bool canSkip = false;
      if (v >= vLep1_Z1SIP && v <= vLep4_Z1SIP && (sc == 1 || sc == 3)) canSkip=true;
      if (v == v4l_chi2 && sc >= 2) canSkip=true;
      if ((v == vLep2_Z1SIP || v == vLep4_Z1SIP) && (sc != 1 && sc != 3)) canSkip=true;
      if (v == vdiffDotTransverse_DirSubPV) canSkip=true;
      if (canSkip) continue;

      double max_plot = 0;

      TPaveText* pt = new TPaveText(0.15, 0.93, 0.85, 1, "brNDC");
      pt->SetBorderSize(0);
      pt->SetFillStyle(0);
      pt->SetTextAlign(12);
      pt->SetTextFont(42);
      pt->SetTextSize(0.045);
      TText* text = pt->AddText(0.025, 0.45, "#font[61]{CMS}");
      text->SetTextSize(0.044);
      text = pt->AddText(0.165, 0.42, "#font[52]{Unpublished}");
      text->SetTextSize(0.0315);
      TString cErgTev = "#font[42]{19.7 fb^{-1} (8 TeV) + 5.1 fb^{-1} (7 TeV)}";
      text = pt->AddText(0.537, 0.45, cErgTev);
      text->SetTextSize(0.0315);

      hvars[7][v][sc]->SetLineColor(kBlack);
      hvars[7][v][sc]->SetLineWidth(2);

      for (int iProcess=0; iProcess<7; iProcess++){
        if (hvars[iProcess][v][sc]==0) continue;
        hvars[iProcess][v][sc]->GetXaxis()->SetLabelFont(42);
        hvars[iProcess][v][sc]->GetXaxis()->SetLabelOffset(0.007);
        hvars[iProcess][v][sc]->GetXaxis()->SetLabelSize(0.04);
        hvars[iProcess][v][sc]->GetXaxis()->SetTitleSize(0.06);
        hvars[iProcess][v][sc]->GetXaxis()->SetTitleOffset(0.9);
        hvars[iProcess][v][sc]->GetXaxis()->SetTitleFont(42);
        hvars[iProcess][v][sc]->GetYaxis()->SetNdivisions(505);
        hvars[iProcess][v][sc]->GetYaxis()->SetLabelFont(42);
        hvars[iProcess][v][sc]->GetYaxis()->SetLabelOffset(0.007);
        hvars[iProcess][v][sc]->GetYaxis()->SetLabelSize(0.04);
        hvars[iProcess][v][sc]->GetYaxis()->SetTitleSize(0.06);
        hvars[iProcess][v][sc]->GetYaxis()->SetTitleOffset(1.1);
        hvars[iProcess][v][sc]->GetYaxis()->SetTitleFont(42);

        if ((v==vTxy_BS || v==vTxy) && iProcess<6) symmetrize_PromptTemplates(hvars[iProcess][v][sc]);
      }

      int ndata = 0;
      double xx_data[100];
      double xu_data[100];
      double xd_data[100];
      double yy_data[100];
      double yu_data[100];
      double yd_data[100];
      double rr_data[100];
      double ru_data[100];
      double rd_data[100];
      const double quant = (1.0 - 0.6827) / 2.0;
      double integral_data = 1;
      //      if (normOption==1) integral_data = hvars[7][v][sc]->Integral();

      double maxplot=0;
      for (int bin = 1; bin <= hvars[7][v][sc]->GetNbinsX(); bin++){
        double bincenter = hvars[7][v][sc]->GetBinCenter(bin);
        double bincontent = hvars[7][v][sc]->GetBinContent(bin);

        if (bincontent > 0){
          xx_data[ndata] = bincenter;
          yy_data[ndata] = bincontent / integral_data;
          xu_data[ndata] = 0;
          xd_data[ndata] = 0;
          yu_data[ndata] = (ROOT::Math::chisquared_quantile_c(quant, 2 * (bincontent + 1)) / 2. - bincontent) / integral_data;
          yd_data[ndata] = ((bincontent == 0) ? 0 : (bincontent - ROOT::Math::chisquared_quantile_c(1 - quant, 2 * bincontent) / 2.)) / integral_data;

          double high_data = yy_data[ndata] + yu_data[ndata];
          if (high_data > maxplot) maxplot = high_data;
          ndata++;
        }
      }
      gROOT->cd();
      TGraphAsymmErrors* tgdata = new TGraphAsymmErrors(ndata, xx_data, yy_data, xd_data, xu_data, yd_data, yu_data);
      tgdata->SetName("tgdata");
      tgdata->SetMarkerSize(1.2);
      tgdata->SetMarkerStyle(20);
      tgdata->SetMarkerColor(kBlack);
      tgdata->SetLineColor(kBlack);
      tgdata->SetLineWidth(1);

      gROOT->cd();
      TString appendName;
      appendName = "_";
      appendName += cutNames[sc+1];
      if (iRegion == 0) appendName += "_SignalRegion_AllTeV";
      else if (iRegion<4) appendName += Form("_SidebandRegion%i_AllTeV", iRegion);
      else if (iRegion==4) appendName += "_SidebandRegions13_AllTeV";
      TString canvasname = "cCompareSIPCuts_";
      canvasname.Append(varNames[v]);
      canvasname.Append(appendName);
      TCanvas* cc = new TCanvas(canvasname, "", 8, 30, 800, 800);
      cc->cd();
      gStyle->SetOptStat(0);
      cc->SetFillColor(0);
      cc->SetBorderMode(0);
      cc->SetBorderSize(2);
      cc->SetTickx(1);
      cc->SetTicky(1);
      cc->SetLeftMargin(0.17);
      cc->SetRightMargin(0.05);
      cc->SetTopMargin(0.07);
      cc->SetBottomMargin(0.13);
      cc->SetFrameFillStyle(0);
      cc->SetFrameBorderMode(0);
      cc->SetFrameFillStyle(0);
      cc->SetFrameBorderMode(0);

      TLegend *ll;
      ll = new TLegend(0.22, (iRegion!=0 ? 0.66 : 0.54), 0.60, 0.90);
      ll->SetBorderSize(0);
      ll->SetTextFont(42);
      ll->SetTextSize(0.04);
      ll->SetLineColor(1);
      ll->SetLineStyle(1);
      ll->SetLineWidth(1);
      ll->SetFillColor(0);
      ll->SetFillStyle(0);

      ll->AddEntry(tgdata, "Observed", "ep");
      if (sc==3) cout << "New cut \t | ZX norm: " << hvars[6][v][sc]->Integral() << "\tggZZ norm: " << hvars[5][v][sc]->Integral() << "\tqqZZ norm: " << hvars[4][v][sc]->Integral() << endl;
      if (sc==0) cout << "No cut \t | ZX norm: " << hvars[6][v][sc]->Integral() << "\tggZZ norm: " << hvars[5][v][sc]->Integral() << "\tqqZZ norm: " << hvars[4][v][sc]->Integral() << endl;

      hvars[5][v][sc]->Add(hvars[6][v][sc]);
      hvars[4][v][sc]->Add(hvars[5][v][sc]);
      if (hvars[3][v][sc]!=0 && hvars[0][v][sc]!=0) hvars[3][v][sc]->Scale(hvars[0][v][sc]->Integral()/hvars[3][v][sc]->Integral());
      if (hvars[0][v][sc]!=0){
        hvars[0][v][sc]->SetLineWidth(2);
        hvars[0][v][sc]->SetLineColor(kRed);
        cout << "Higgs norm: " << hvars[0][v][sc]->Integral() << endl;
        hvars[0][v][sc]->Add(hvars[4][v][sc]);
//        hvars[0][v][sc]->Draw("histsame");
        ll->AddEntry(hvars[0][v][sc], processTitle[0], "l");
      }
      if (hvars[2][v][sc]!=0){
        hvars[2][v][sc]->SetLineWidth(2);
        hvars[2][v][sc]->SetLineStyle(7);
        hvars[2][v][sc]->SetLineColor(kViolet);
        cout << "Higgs (BSM) norm: " << hvars[2][v][sc]->Integral() << endl;
        hvars[2][v][sc]->Add(hvars[4][v][sc]);
//        hvars[2][v][sc]->Draw("histsame");
        ll->AddEntry(hvars[2][v][sc], processTitle[2], "l");
      }
      hvars[4][v][sc]->SetLineWidth(2);
      hvars[5][v][sc]->SetLineWidth(2);
      hvars[6][v][sc]->SetLineWidth(2);
      hvars[4][v][sc]->SetLineColor(kBlack);
      hvars[5][v][sc]->SetLineColor(kBlack);
      hvars[6][v][sc]->SetLineColor(kBlack);
      hvars[4][v][sc]->SetFillColor(kAzure-9);
      hvars[5][v][sc]->SetFillColor(kAzure-2);
      hvars[6][v][sc]->SetFillColor(TColor::GetColor("#669966"));
      hvars[4][v][sc]->SetFillStyle(1001);
      hvars[5][v][sc]->SetFillStyle(1001);
      hvars[6][v][sc]->SetFillStyle(1001);
      for (int pp=4; pp<=6; pp++) ll->AddEntry(hvars[pp][v][sc], processTitle[pp], "f");

      double yaxmax = maxplot;
      if (iRegion==0 && v==v4l_chi2) yaxmax *= 2.;
      else if (iRegion==4 && v==v4l_chi2) yaxmax *= 1.7;
      else if (iRegion==4) yaxmax *= 1.3;
      else if (iRegion==0) yaxmax *= 1.5;
      else yaxmax *= 1.7;
      hvars[4][v][sc]->GetYaxis()->SetRangeUser(0, yaxmax);
      if (v==vDbkg) hvars[4][v][sc]->GetXaxis()->SetRangeUser(KD_ranges[v][0], KD_ranges[v][1]-0.001);
      hvars[4][v][sc]->Draw("hist");
      hvars[5][v][sc]->Draw("histsame");
      hvars[6][v][sc]->Draw("histsame");
      if (hvars[0][v][sc]!=0) hvars[0][v][sc]->Draw("histsame");
      if (hvars[2][v][sc]!=0) hvars[2][v][sc]->Draw("histsame");
      tgdata->Draw("e1psame");
      ll->Draw("same");
      pt->Draw();

      TPaveText *pt10 = new TPaveText((sc==newSIPCuts-2 ? 0.67 : 0.63), (sc==newSIPCuts-2 ? 0.85 : 0.80), 0.90, (sc==newSIPCuts-2 ? 0.90 : 0.85), "brNDC");
      pt10->SetBorderSize(0);
      pt10->SetTextAlign(12);
      pt10->SetTextSize(0.04);
      pt10->SetFillStyle(0);
      pt10->SetTextFont(42);
      TText* text10;
      if (iRegion==0) text10 = pt10->AddText(0.01, 0.01, "Signal region");
      else if (iRegion<4) text10 = pt10->AddText(0.01, 0.01, "Control region sideband");
      else if (iRegion==4) text10 = pt10->AddText(0.01, 0.01, "Sideband m_{4l}");
      if(iRegion!=0) pt10->Draw();

      TPaveText *pt20 = new TPaveText((sc==newSIPCuts-2 ? 0.67 : 0.63), 0.85, 0.90, 0.90, "brNDC");
      pt20->SetBorderSize(0);
      pt20->SetTextAlign(12);
      pt20->SetTextSize(0.04);
      pt20->SetFillStyle(0);
      pt20->SetTextFont(42);
      TText* text20;
      text20 = pt20->AddText(0.01, 0.01, ctitle_SIPCut[sc+1]);
      if (!(iRegion==0 && (v==vDbkg || v==vTxy || v==vTxy_BS)) && sc!=newSIPCuts-2) pt20->Draw();

      cc->RedrawAxis();
      cc->Modified();
      cc->Update();
      canvasname.Prepend(plotDir_1D);
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
      cc->SaveAs(canvasname_pdf);
      cc->SaveAs(canvasname_eps);
      cc->SaveAs(canvasname_png);
      cc->SaveAs(canvasname_root);
      cc->SaveAs(canvasname_c);
      delete pt20;
      delete pt10;
      delete ll;
      cc->Close();
      delete tgdata;
      delete pt;
    }
  }

  for (int hh=0; hh<8; hh++){
    for (int v = 0; v < nVariables; v++){
      for (int sc = 0; sc < newSIPCuts; sc++){
        if (hvars[hh][v][sc]!=0) delete hvars[hh][v][sc];
      }
    }
  }
}

void generic_Histo2DPlotter(TFile* foutput, TString canvasDir, TH2* hSlice, int erg_tev, float lumi, int channel, bool isCRonly){
  gStyle->SetOptStat(0);
  gStyle->SetPalette(1);
  gStyle->SetPadRightMargin(0.20);
  gROOT->ForceStyle();
  foutput->cd();
  gStyle->SetOptStat(0);
  TString canvasname_2D;
  if (!isCRonly) canvasname_2D = Form("c%s_%s_%iTeV", hSlice->GetName(), user_folder[channel], erg_tev);
  else canvasname_2D = Form("c%s_%iTeV", hSlice->GetName(), erg_tev);
  TCanvas* c = new TCanvas(canvasname_2D, "", 1000, 800);
  gStyle->SetOptStat(0);
  c->SetFillColor(0);
  c->SetBorderMode(0);
  c->SetBorderSize(2);
  c->SetTickx(1);
  c->SetTicky(1);
  c->cd();
  c->SetFrameFillStyle(0);
  c->SetFrameBorderMode(0);
  c->SetFrameFillStyle(0);
  c->SetFrameBorderMode(0);
  c->cd();
  hSlice->SetTitle("");
  hSlice->GetXaxis()->SetLabelSize(0.04);
  hSlice->GetYaxis()->SetLabelSize(0.04);
  hSlice->GetXaxis()->SetTitleOffset(0.95);
  hSlice->GetYaxis()->SetTitleOffset(1.1);
  hSlice->Draw("colz");
  TPaveText* pt = new TPaveText(0.15, 0.955, 0.92, 1, "brNDC");
  pt->SetBorderSize(0);
  pt->SetFillStyle(0);
  pt->SetTextAlign(12);
  pt->SetTextFont(42);
  pt->SetTextSize(0.04);
  TText* text = pt->AddText(0.02, 0.45, "#font[61]{CMS}");
  text->SetTextSize(0.044);
  text = pt->AddText(0.12, 0.42, "#font[52]{Unpublished}");
  text->SetTextSize(0.0315);
  TString cErgTev = Form("#font[42]{%.1f fb^{-1} (%i TeV)}", lumi, erg_tev);
  if (!isCRonly){
    if (channel == 0) cErgTev.Prepend("4#mu ");
    if (channel == 1) cErgTev.Prepend("4e ");
    if (channel == 2) cErgTev.Prepend("2e+2#mu ");
  }
  if (!isCRonly) text = pt->AddText((channel<2 ? 0.575 : 0.545), 0.45, cErgTev);
  else text = pt->AddText(0.625, 0.45, cErgTev);
  //	text = pt->AddText(0.41,0.45,"#font[42]{19.7 fb^{-1} (8 TeV) + 5.1 fb^{-1} (7 TeV)}");
  text->SetTextSize(0.0315);
  pt->Draw();
  TPaveText* pt2 = new TPaveText(0.985, 0.72, 0.99, 0.99, "brNDC");
  pt2->SetBorderSize(0);
  pt2->SetFillStyle(0);
  pt2->SetTextAlign(21);
  pt2->SetTextFont(42);
  pt2->SetTextSize(0.04);
  TText* text2 = pt2->AddText(0.01, 0.3, "#font[42]{Events / bin}");
  text2->SetTextSize(0.056);
  text2->SetTextAngle(90);
  pt2->Draw();
  foutput->cd();
  foutput->WriteTObject(c);

  canvasname_2D.Prepend(canvasDir);
  TString canvasname_2D_pdf = canvasname_2D;
  TString canvasname_2D_eps = canvasname_2D;
  TString canvasname_2D_png = canvasname_2D;
  TString canvasname_2D_root = canvasname_2D;
  TString canvasname_2D_c = canvasname_2D;
  canvasname_2D_pdf.Append(".pdf");
  canvasname_2D_eps.Append(".eps");
  canvasname_2D_png.Append(".png");
  canvasname_2D_root.Append(".root");
  canvasname_2D_c.Append(".C");
  c->SaveAs(canvasname_2D_pdf);
  c->SaveAs(canvasname_2D_eps);
  c->SaveAs(canvasname_2D_png);
  c->SaveAs(canvasname_2D_root);
  c->SaveAs(canvasname_2D_c);

  delete pt; delete pt2;
  c->Close();
}

float compute_subtractedPV(
  float CandVtx_x,
  float CandVtx_y,
  float CandVtx_z,
  float CandVtx_cov_xx,
  float CandVtx_cov_xy,
  float CandVtx_cov_xz,
  float CandVtx_cov_yy,
  float CandVtx_cov_yz,
  float CandVtx_cov_zz,
  float PrimaryVtx_x,
  float PrimaryVtx_y,
  float PrimaryVtx_z,
  float PrimaryVtx_cov_xx,
  float PrimaryVtx_cov_xy,
  float PrimaryVtx_cov_xz,
  float PrimaryVtx_cov_yy,
  float PrimaryVtx_cov_yz,
  float PrimaryVtx_cov_zz,
  float& SubtractedVtx_x,
  float& SubtractedVtx_y,
  float& SubtractedVtx_z
  ){
  double cand_vtx[3] ={
    CandVtx_x,
    CandVtx_y,
    CandVtx_z
  };
  double opv_vtx[3] ={
    PrimaryVtx_x,
    PrimaryVtx_y,
    PrimaryVtx_z
  };
  double cand_vtx_T[3] ={
    CandVtx_x,
    CandVtx_y,
    0
  };
  double opv_vtx_T[3] ={
    PrimaryVtx_x,
    PrimaryVtx_y,
    0
  };
  double cand_cov[9] ={
    CandVtx_cov_xx,
    CandVtx_cov_xy,
    CandVtx_cov_xz,
    CandVtx_cov_xy,
    CandVtx_cov_yy,
    CandVtx_cov_yz,
    CandVtx_cov_xz,
    CandVtx_cov_yz,
    CandVtx_cov_zz
  };
  double opv_cov[9] ={
    PrimaryVtx_cov_xx,
    PrimaryVtx_cov_xy,
    PrimaryVtx_cov_xz,
    PrimaryVtx_cov_xy,
    PrimaryVtx_cov_yy,
    PrimaryVtx_cov_yz,
    PrimaryVtx_cov_xz,
    PrimaryVtx_cov_yz,
    PrimaryVtx_cov_zz
  };

  TMatrixD candCov(3, 3);
  TMatrixD opvCov(3, 3);
  TVectorD candVtx(3, cand_vtx);
  TVectorD opvVtx(3, opv_vtx);
  candCov.SetMatrixArray(cand_cov);
  opvCov.SetMatrixArray(opv_cov);

  TVectorD candVtxTransverse(3, cand_vtx_T);
  TVectorD opvVtxTransverse(3, opv_vtx_T);

  TVectorD orient = candVtxTransverse - opvVtxTransverse;
  double norm_orient = orient.Norm2Sqr();
  if (norm_orient != 0) orient *= 1. / sqrt(norm_orient);
  double orienterr_opv = opvCov.Similarity(orient);
  double orienterr_cand = candCov.Similarity(orient);
  double errratio = orienterr_opv / orienterr_cand;
  double errratio_return = errratio;
  double errscale = 1. - errratio;
  bool zeroScale = false;
  bool negativeScale = false;
  if (errscale == 0) zeroScale = true;
  else if (errscale < 0) negativeScale = true;

  candVtx *= -errratio;
  if (zeroScale) errratio = 0.9;
  TVectorD result = (opvVtx + candVtx);
  result *= (1. / errscale);

  SubtractedVtx_x = result[0];
  SubtractedVtx_y = result[1];
  SubtractedVtx_z = result[2];

  return ((float)errratio_return);
}


void symmetrize_PromptTemplates(TH1F* hrepair){
  double oldInt = hrepair->Integral();
  int nbins = hrepair->GetNbinsX();
  int left_minbin = -1, left_maxbin=-1, right_minbin=-1, right_maxbin=-1;
  if (nbins % 2 == 0){
    left_minbin = 1; left_maxbin=nbins/2; right_minbin=left_maxbin+1; right_maxbin=nbins;
  }
  else{
    left_minbin = 1; left_maxbin=(nbins-1)/2; right_minbin=left_maxbin+2; right_maxbin=nbins;
  }
  for (int bin=right_minbin; bin<=right_maxbin; bin++){
    int bin2 = left_maxbin - (bin-right_minbin);
    double bincontent1 = hrepair->GetBinContent(bin);
    double bincontent2 = hrepair->GetBinContent(bin2);
    double bincontent = (bincontent1 + bincontent2)/2;

    double binerror1 = hrepair->GetBinError(bin);
    double binerror2 = hrepair->GetBinError(bin2);
    double binerror = (binerror1 + binerror2)/2;

    hrepair->SetBinContent(bin, bincontent);
    hrepair->SetBinContent(bin2, bincontent);
    hrepair->SetBinError(bin, binerror);
    hrepair->SetBinError(bin2, binerror);
  }
  double newInt = hrepair->Integral();
  newInt /= oldInt;
  hrepair->Scale(newInt);
}


void plotSIPVariables(int recreateYields=0){
  if (recreateYields==1) create_SIPCut_YieldRatios();
  for (int p=0; p<=7; p++){
    for (int r=0; r<4; r++){
      produce_SIPVariables(p, r);
    }
  }
  for (int r=0; r<4; r++){
    plotCumulative_SIPVariables(r);
  }
}

