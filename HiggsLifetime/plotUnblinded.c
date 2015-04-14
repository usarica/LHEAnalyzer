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

const int nProcesses = 24;
string processName[nProcesses] ={
  "CTau0", "CTau100", "CTau500", "CTau1000",
  "CTau0_VBF", "CTau100_VBF", "CTau500_VBF", "CTau1000_VBF",
  "CTau0_WH", "CTau100_WH", "CTau500_WH", "CTau1000_WH",
  "CTau0_ZH", "CTau100_ZH", "CTau500_ZH", "CTau1000_ZH",
  "CTau0_ttH", "CTau100_ttH", "CTau500_ttH", "CTau1000_ttH",
  "qqZZ", "ggZZ", "CR", "data"
};
TString processTitle[8] ={
  "SM signal", "c#tau_{H}=100 #mum", "c#tau_{H}=500 #mum", "c#tau_{H}=1000 #mum",
  "q#bar{q}#rightarrow4l bkg.", "gg#rightarrow4l bkg.", "Z+X", "Observed"
};
TString productionName[5]={ "ggH", "VBFH", "WH", "ZH", "ttH" };
TString strHiggsProductionFit[4] ={ "VBFH125", "WH125", "ZH125", "ttH125" };
TString templatesMainDir = "/Analysis/CMSSW_6_1_1/src/HiggsAnalysis/HZZ4l_Combination/CreateDatacards/templates2D/HCTau_NewSIP_BS/";

enum{
  vTxy_BS,
  vDbkg,
  vpT4l,
  vm4l,
  nVariables
};
TString varNames[nVariables] ={
  "Txy_BS",
  "Dbkg",
  "pT",
  "mZZ"
};
TString varTitles[nVariables] ={
  "c#Deltat (#mum)",
  "D_{bkg}",
  "p_{T} (GeV)",
  "m_{4l} (GeV)"
};
float KD[nVariables] ={ 0 };
float KD_ranges[nVariables][2] ={
  { -1000, 1000 },
  { 0, 1 },
  { 0, 200 },
  { 105.6, 140.6 }
};
int KD_nbins[nVariables] ={
  50,
  20,
  20,
  14
};
const int nVariables_2D=2;
int KD_pairing[nVariables_2D][2] ={
  { vTxy_BS, vDbkg },
  { vm4l, vpT4l }
};

const int nMZZ_ranges=1;
float systZZMass_range[nMZZ_ranges][2] ={
  { 105.6, 140.6 }
};

TString cutName = "Z1SIP_4lchi2cut";
TString cutLabel = "New cut";
TString ctitle_SIPCut = "#it{SIP}^{Z1}_{1,2}<4, #it{SIP}^{Z1}_{3,4}<5, #chi^{2}_{4l}<30";

TString templateMainDir = "/Analysis/CMSSW_6_1_1/src/HiggsAnalysis/HZZ4l_Combination/CreateDatacards/templates2D/HCTau_NewSIP_BS/";

void generic_Histo2DPlotter(TFile* foutput, TString canvasDir, TH2* hSlice, int erg_tev, float lumi, int channel, bool isCRonly=false);


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
void symmetrize_PromptTemplates(TH1F* hrepair);



void plotUnblinded_single(int iProcess, int isEnriched=0){
  int iRegion = 0;
  gROOT->ProcessLine(".x tdrstyle.cc");
  char TREE_NAME[] = "SelectedTree";
  int processMap[nProcesses-1] ={
    0, 1, 2, 3,
    0, 1, 2, 3,
    0, 1, 2, 3,
    0, 1, 2, 3,
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

  float D_bkg, D_bkg_kin, p0plus_VAJHU, bkg_VAMCFM, p0plus_m4l, bkg_m4l;

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

  TH1F* hvars[nVariables] ={ 0 };
  TH2F* hvars_2D[nVariables_2D] ={ 0 };
  TString hname_core = "h_";
  for (int v = 0; v < nVariables; v++){
    float binwidth = (KD_ranges[v][1] - KD_ranges[v][0]) / KD_nbins[v];
    int nbins = KD_nbins[v];
    float xmin = KD_ranges[v][0];
    float xmax = KD_ranges[v][1];
    if (isEnriched && v==vTxy_BS){
      xmin *= 0.5;
      xmax *= 0.5;
      nbins /= 2;
    }
    if (!(v==vDbkg || v==vm4l)){
      xmax += binwidth;
      nbins++;
    }
    if (xmin < 0 && v!=vm4l){
      xmin -= binwidth; nbins++;
    }

    TString hname = hname_core;
    hname.Append(varNames[v]);
    hname.Append("_");
    hname.Append(cutName);
    hvars[v] = new TH1F(hname, "", nbins, xmin, xmax);
    hvars[v]->Sumw2();
    hvars[v]->SetXTitle(varTitles[v]);
    hvars[v]->SetYTitle("Events / bin");
    hvars[v]->GetXaxis()->SetNdivisions(505);
    hvars[v]->GetXaxis()->SetLabelFont(42);
    hvars[v]->GetXaxis()->SetLabelOffset(0.007);
    hvars[v]->GetXaxis()->SetLabelSize(0.04);
    hvars[v]->GetXaxis()->SetTitleSize(0.06);
    hvars[v]->GetXaxis()->SetTitleOffset(0.9);
    hvars[v]->GetXaxis()->SetTitleFont(42);
    hvars[v]->GetYaxis()->SetNdivisions(505);
    hvars[v]->GetYaxis()->SetLabelFont(42);
    hvars[v]->GetYaxis()->SetLabelOffset(0.007);
    hvars[v]->GetYaxis()->SetLabelSize(0.04);
    hvars[v]->GetYaxis()->SetTitleSize(0.06);
    hvars[v]->GetYaxis()->SetTitleOffset(1.1);
    hvars[v]->GetYaxis()->SetTitleFont(42);
  }
  for (int v = 0; v < nVariables_2D; v++){
    int vx=KD_pairing[v][0];
    float binwidth_x = (KD_ranges[vx][1] - KD_ranges[vx][0]) / KD_nbins[vx];
    int nbins_x = KD_nbins[vx];
    float xmin = KD_ranges[vx][0];
    float xmax = KD_ranges[vx][1];
    if (!(vx==vDbkg || vx==vm4l)){
      xmax += binwidth_x;
      nbins_x++;
    }
    if (xmin < 0 && vx!=vm4l){
      xmin -= binwidth_x; nbins_x++;
    }

    int vy=KD_pairing[v][1];
    float binwidth_y = (KD_ranges[vy][1] - KD_ranges[vy][0]) / KD_nbins[vy];
    int nbins_y = KD_nbins[vy];
    float ymin = KD_ranges[vy][0];
    float ymax = KD_ranges[vy][1];
    if (vy!=vDbkg){
      ymax += binwidth_y;
      nbins_y++;
    }
    if (ymin < 0){
      ymin -= binwidth_y; nbins_y++;
    }

    TString hname = hname_core;
    hname.Append(varNames[vy]);
    hname.Append("_vs_");
    hname.Append(varNames[vx]);
    hname.Append("_");
    hname.Append(cutName);
    hvars_2D[v] = new TH2F(hname, "", nbins_x, xmin, xmax, nbins_y, ymin, ymax);
    hvars_2D[v]->Sumw2();
    hvars_2D[v]->SetXTitle(varTitles[vx]);
    hvars_2D[v]->SetYTitle(varTitles[vy]);
    hvars_2D[v]->GetXaxis()->SetNdivisions(505);
    hvars_2D[v]->GetXaxis()->SetLabelFont(42);
    hvars_2D[v]->GetXaxis()->SetLabelOffset(0.007);
    hvars_2D[v]->GetXaxis()->SetLabelSize(0.04);
    hvars_2D[v]->GetXaxis()->SetTitleSize(0.06);
    hvars_2D[v]->GetXaxis()->SetTitleOffset(1.0);
    hvars_2D[v]->GetXaxis()->SetTitleFont(42);
    hvars_2D[v]->GetYaxis()->SetNdivisions(505);
    hvars_2D[v]->GetYaxis()->SetLabelFont(42);
    hvars_2D[v]->GetYaxis()->SetLabelOffset(0.007);
    hvars_2D[v]->GetYaxis()->SetLabelSize(0.04);
    hvars_2D[v]->GetYaxis()->SetTitleSize(0.06);
    hvars_2D[v]->GetYaxis()->SetTitleOffset(1.1);
    hvars_2D[v]->GetYaxis()->SetTitleFont(42);
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
      TString OUTPUT_NAME = "LifetimeKD_Unblinded_";
      OUTPUT_NAME.Append(Form("%s", processName[iProcess].c_str()));
      OUTPUT_NAME.Append(Form("_%s_", user_folder[folder]));
      OUTPUT_NAME.Append(comstring);
      if (isEnriched==1) OUTPUT_NAME.Append("_SignalEnriched");
      OUTPUT_NAME.Append(".root");

      TString cinput_common_noSIP = user_dir_hep + "No_SIP/" + erg_dir + "/";
      if (iProcess < nProcesses-1){
        cinput_common_noSIP.Append(Form("%s/", user_folder[folder]));
      }
      else{
        cinput_common_noSIP.Append(Form("%s/", user_folder[3]));
      }
      TString cinput_common = cinput_common_noSIP;
      TString coutput_common = user_dir_hep + "Analysis/Auxiliary/";

      TString plotDir = coutput_common;
      plotDir.Append("../Plots/Unblinded/");
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
      TH2F* hRatio;

      TString chratio = "hCR_MC_OSSSRatios_";
      chratio.Append(cutLabel);
      hRatio = (TH2F*)fSSinput->Get(chratio);
      cout << "Obtained ratio histogram " << hRatio->GetName() << endl;

      double ratio_targetsignalyield[2] = { 0 }; // New vs Old
      double targetCRcount[2] = { 0 }; // OS vs SS for CR
      float m4l_lowhigh[2] = { 105.6, 140.6 };
      int pCode = iProcess;
      if (iProcess<20) pCode = 0;
      TString INPUTYIELD_NAME = "LifetimeKD_RelativeSIPYields_";
      INPUTYIELD_NAME.Append(Form("%s", processName[pCode].c_str()));
      INPUTYIELD_NAME.Append(Form("_%s_", user_folder[folder]));
      INPUTYIELD_NAME.Append(comstring);
      INPUTYIELD_NAME.Append(Form("_m4l%.1f_%.1f", m4l_lowhigh[0], m4l_lowhigh[1]));
      INPUTYIELD_NAME.Append(".root");
      TString cinput_yield_common = user_dir_hep + "Analysis/Auxiliary/";
      TString cinput_yield = cinput_yield_common + INPUTYIELD_NAME;
      TFile* finput = new TFile(cinput_yield, "read");
      TH1F* htemp;
      if (finput==0 || finput->IsZombie()){
        cout << cinput_yield << " not found." << endl;
      }
      else{
        if (pCode!=22){
          htemp = (TH1F*)finput->Get(Form("h%s", processName[pCode].c_str()));
          ratio_targetsignalyield[1] = htemp->GetBinContent(5);
          ratio_targetsignalyield[0] = htemp->GetBinContent(1);
          delete htemp;
        }
        else{
          htemp = (TH1F*)finput->Get(Form("h%s_OS", processName[pCode].c_str()));
          ratio_targetsignalyield[1] = htemp->GetBinContent(5);
          ratio_targetsignalyield[0] = htemp->GetBinContent(1);
          delete htemp;
          htemp = (TH1F*)finput->Get(Form("h%s_SS", processName[pCode].c_str()));
          ratio_targetsignalyield[1] += htemp->GetBinContent(5);
          ratio_targetsignalyield[0] += htemp->GetBinContent(1);
          delete htemp;

          htemp = (TH1F*)finput->Get(Form("h%s_Unweighted_OS", processName[pCode].c_str()));
          targetCRcount[0] = htemp->GetBinContent(5);
          delete htemp;
          htemp = (TH1F*)finput->Get(Form("h%s_Unweighted_SS", processName[pCode].c_str()));
          targetCRcount[1] = htemp->GetBinContent(5);
          delete htemp;
        }
        finput->Close();
      }
      double cutscale = 1;
      if(iProcess<23) cutscale = ratio_targetsignalyield[1]/ratio_targetsignalyield[0];
      cout << "Cut scale: " << cutscale << endl;

      TString coutput = coutput_common + OUTPUT_NAME;
      TFile* foutput = new TFile(coutput, "recreate");

      TString templatesDir = templatesMainDir;
      templatesDir.Prepend(user_dir_hep);
      templatesDir.Append(comstring);
      templatesDir.Append("/");
      if (folder==2) templatesDir.Append("2e2mu");
      else templatesDir.Append(user_folder[folder]);
      templatesDir.Append("_templates_");
      TH1F* htemplate_Dbkg=0;
      TString templateFile_Dbkg = templatesDir;
      if (iProcess<20){
        templateFile_Dbkg.Append("SignalScaleResSyst.root");
        finput = new TFile(templateFile_Dbkg, "read");
        TString hname = "T_2D_ScaleResNominal";
        htemp = (TH1F*)finput->Get(hname);
        hname.Append("_template");
        foutput->cd();
        htemplate_Dbkg = (TH1F*)htemp->Clone(hname);
        finput->cd();
        delete htemp;
        finput->Close();
        foutput->cd();
        htemplate_Dbkg->Scale(1./htemplate_Dbkg->Integral());
      }
      else if(iProcess<23){
        templateFile_Dbkg.Append("Merged_bkg.root");
        finput = new TFile(templateFile_Dbkg, "read");
        TString hname;
        if (iProcess<22) hname = Form("template_%s_Dbkg", processName[iProcess].c_str());
        else hname = "template_ZX_Dbkg_Nominal";
        htemp = (TH1F*)finput->Get(hname);
        hname.Append("_template");
        foutput->cd();
        htemplate_Dbkg = (TH1F*)htemp->Clone(hname);
        finput->cd();
        delete htemp;
        finput->Close();
        foutput->cd();
        htemplate_Dbkg->Scale(1./htemplate_Dbkg->Integral());
      }
      double signalEnrichedScale = 1;
      if (isEnriched==1 && htemplate_Dbkg!=0) signalEnrichedScale = htemplate_Dbkg->Integral(26, 50);
      cout << "Enriching scale: " << signalEnrichedScale << endl;
      cutscale *= signalEnrichedScale;

      TFile* finput_pTrewgt = new TFile(Form("./data/compareProductionModes_pT_%s%s", comstring.Data(), ".root"), "read");
      TF1* fit_pT = 0;
      if (iProcess>=4 && iProcess<20){
        int iProd = (iProcess-pCode)/4;
        if (iProd>0){
          TString fitname = strHiggsProductionFit[iProd-1];
          fitname.Append("_pToverMZZ_ratio_fit");
          cout << "Retrieving " << fitname;
          fit_pT = (TF1*)finput_pTrewgt->Get(fitname);
          cout << ", address: " << fit_pT << endl;
        }
      }
      foutput->cd();

      TChain* tc;
      TChain* tc_extra=0;
      TString cinput;
      int smp = processMap[iProcess];
      bool treesExist = true;
      tc = new TChain(TREE_NAME);

      if (iProcess < 20){
        cinput = cinput_common + sample_FullSim[smp] + ".root";
        tc->Add(cinput);
      }
      else if (iProcess == 20 || iProcess == 21){
        int iProcess_new = iProcess;
        smp = processMap[iProcess_new];
        for (int it = 0; it < (iProcess_new == 20 ? 6 : 3); it++){ // qqZZ has 6 samples, gg has 3 samples
          cinput = cinput_common + sample_FullSim[smp] + "_Reprocessed.root";
          tc->Add(cinput);
          smp++;
        }
      }
      else if (iProcess == 22){
        cinput = user_dir_hep;
        cinput.Append("No_SIP/");
        cinput.Append(erg_dir);
        cinput = cinput + "/CR/" + sample_FullSim[smp] + ".root";
        tc->Add(cinput);
      }
      else if (iProcess == 23){
        cinput = cinput_common + data_files[folder] + ".root";
        tc->Add(cinput);
      }
      if (iRegion==0 && (iProcess==20 || iProcess==21)){
        tc_extra = new TChain(TREE_NAME);
        if (iProcess==20 || iProcess==21){
          for (int it = kQQBZZSamples; it < kQQBZZSamples_Dedicated; it++){ // qq/ggZZ extra for signal region
            cinput = cinput_common + sample_FullSim[it] + "_Reprocessed.root";
            tc_extra->Add(cinput);
            smp++;
          }
        }
        if (iProcess==21){
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

        if (iProcess < 22){
          tc->SetBranchAddress("MC_weight", &MC_weight);
          tc->SetBranchAddress("MC_weight_noxsec", &MC_weight_noxsec);
        }
        else if (iProcess == 22){
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
        if (iProcess>=22){
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
          tc->SetBranchAddress("D_bkg_kin", &D_bkg_kin);
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
            tc_extra->SetBranchAddress("D_bkg_kin", &D_bkg_kin);
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

      for (int v = 0; v < nVariables; v++){
        hvars[v]->Reset("ICESM");
        hvars[v]->SetStats(0);
        hvars[v]->Sumw2();
      }
      for (int v = 0; v < nVariables_2D; v++){
        hvars_2D[v]->Reset("ICESM");
        hvars_2D[v]->SetStats(0);
        hvars_2D[v]->Sumw2();
      }

      double sum_signal = 0;

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

        if (iProcess==20) MC_weight = MC_weight_noxsec*MC_weight_QQBGGProper[1];
        if (iProcess==20) MC_weight *= MC_weight_QQZZEWK;

        if (iProcess==21) MC_weight = MC_weight_noxsec;
        if (iProcess==21) MC_weight *= MC_weight_QQBGGProper[0];
        if (iProcess==21) MC_weight *= MC_weight_Kfactor;

        // LEFT HERE
        if (iProcess==22){
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
        if (fabs(Lep1_Z1SIP) >= 4 || fabs(Lep2_Z1SIP) >= 4 || fabs(Lep3_Z1SIP) >= 5 || fabs(Lep4_Z1SIP) >= 5 || KalmanCandVtx_chi2 >= 30) continue;
        if (ZZMass >= rangeZZMass[1] || ZZMass < rangeZZMass[0]) continue;
        if (!tc->GetBranchStatus("D_bkg")){
          D_bkg = p0plus_VAJHU*p0plus_m4l / (p0plus_VAJHU*p0plus_m4l + bkg_VAMCFM*bkg_m4l);
          D_bkg_kin = p0plus_VAJHU / (p0plus_VAJHU + bkg_VAMCFM);
        }
        if (isEnriched==1 && D_bkg<=0.5) continue;

        int index_BS = 1;
        bool matchFound = false;
        if (iProcess>=22) index_BS = 0;
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

        KD[vTxy_BS] = Txy_BS*10000.;
        KD[vDbkg] = D_bkg;
        KD[vm4l] = ZZMass;
        KD[vpT4l] = ZZPt;

        if (iProcess==22){
          if (
            (CRflag==6 || CRflag==10) ||
            (CRflag==8 || CRflag==12)
            ){
            MC_weight = ZXfake_weight_OS[4] * (targetCRcount[0] / (targetCRcount[0] + targetCRcount[1]));
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
              if (OSoverSS_biny>0) SSwgt_OSoverSS = hRatio->GetBinContent(hRatio->GetXaxis()->FindBin(ZZMass), OSoverSS_biny);
              //                if (sc==3 && iProcess==6 && iRegion==1 && folder==0 && ZZMass>=110 && ZZMass<130) cout << ZZMass << '\t' << SSwgt_OSoverSS<< '\t' << ZXfake_weight_SS << '\t' << (targetCRcount[1][sc+1][iRegion] / (targetCRcount[0][sc+1][iRegion] + targetCRcount[1][sc+1][iRegion])) << endl;

              MC_weight = ZXfake_weight_SS*SSwgt_OSoverSS;
              MC_weight *= (targetCRcount[1] / (targetCRcount[0] + targetCRcount[1]));
            }
          }
        }

        if (iProcess>=4 && iProcess<20 && fit_pT!=0){
          double pTRatio = GenHPt/GenHMass;
          if (ev==0) cout << pTRatio << '\t';
          double pTscale = fit_pT->Eval(pTRatio);
          if (ev==0) cout << pTscale << endl;
          MC_weight *= pTscale;
        }

        sum_signal += MC_weight;

        float fillvar_2D[nVariables_2D][2] ={ { 0 } };
        for (int v = 0; v < nVariables; v++){
          float fillvar = KD[v];
          if (fillvar<hvars[v]->GetXaxis()->GetXmin()) fillvar = hvars[v]->GetXaxis()->GetXmin();
          if (fillvar>=hvars[v]->GetXaxis()->GetXmax()) fillvar = hvars[v]->GetXaxis()->GetXmax() - (hvars[v]->GetXaxis()->GetBinWidth(1))*0.5;

          hvars[v]->Fill(fillvar, MC_weight);

          for (int k = 0; k < nVariables_2D; k++){
            if (v == KD_pairing[k][0]) fillvar_2D[k][0] = fillvar;
            if (v == KD_pairing[k][1]) fillvar_2D[k][1] = fillvar;
          }
        }
        for (int v = 0; v < nVariables_2D; v++) hvars_2D[v]->Fill(fillvar_2D[v][0], fillvar_2D[v][1], MC_weight);
      }

      if (tc_extra!=0){
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

          if (fabs(Lep1_Z1SIP) >= 4 || fabs(Lep2_Z1SIP) >= 4 || fabs(Lep3_Z1SIP) >= 5 || fabs(Lep4_Z1SIP) >= 5 || KalmanCandVtx_chi2 >= 30) continue;
          if (ZZMass >= rangeZZMass[1] || ZZMass < rangeZZMass[0]) continue;
          if (!tc_extra->GetBranchStatus("D_bkg")){
            D_bkg = p0plus_VAJHU*p0plus_m4l / (p0plus_VAJHU*p0plus_m4l + bkg_VAMCFM*bkg_m4l);
            D_bkg_kin = p0plus_VAJHU / (p0plus_VAJHU + bkg_VAMCFM);
          }
          if (isEnriched==1 && D_bkg<=0.5) continue;

          MC_weight = MC_weight_noxsec;
          if (iProcess==20) MC_weight *= MC_weight_QQZZEWK*MC_weight_QQBGGProper[1];
          if (iProcess==21) MC_weight *= MC_weight_Kfactor*MC_weight_QQBGGProper[0];

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

          KD[vTxy_BS] = Txy_BS*10000.;
          KD[vDbkg] = D_bkg;
          KD[vm4l] = ZZMass;
          KD[vpT4l] = ZZPt;

          sum_signal += MC_weight;

          float fillvar_2D[nVariables_2D][2] ={ { 0 } };
          for (int v = 0; v < nVariables; v++){
            float fillvar = KD[v];
            if (fillvar<hvars[v]->GetXaxis()->GetXmin()) fillvar = hvars[v]->GetXaxis()->GetXmin();
            if (fillvar>=hvars[v]->GetXaxis()->GetXmax()) fillvar = hvars[v]->GetXaxis()->GetXmax() - (hvars[v]->GetXaxis()->GetBinWidth(1))*0.5;

            hvars[v]->Fill(fillvar, MC_weight);

            for (int k = 0; k < nVariables_2D; k++){
              if (v == KD_pairing[k][0]) fillvar_2D[k][0] = fillvar;
              if (v == KD_pairing[k][1]) fillvar_2D[k][1] = fillvar;
            }
          }
          for (int v = 0; v < nVariables_2D; v++) hvars_2D[v]->Fill(fillvar_2D[v][0], fillvar_2D[v][1], MC_weight);
        }
      }
      for (int v = 0; v < nVariables_2D; v++){
        if (iProcess==23) continue;
        TH1F* hproj;
        if (vTxy_BS == KD_pairing[v][0]){
          hproj = (TH1F*) hvars_2D[v]->ProjectionX(Form("%s%S", hvars_2D[v]->GetName(), "tempproj"));
          double initialIntegral = hvars_2D[v]->Integral();
          for (int biny = 1; biny<=hvars_2D[v]->GetNbinsY(); biny++){
            double sliceIntegral = hvars_2D[v]->Integral(1, hvars_2D[v]->Integral(1, hvars_2D[v]->GetNbinsX(), biny, biny));
            double sliceScale = sliceIntegral / initialIntegral;
            for (int binx=1; binx<=hvars_2D[v]->GetNbinsX(); binx++){
              hvars_2D[v]->SetBinContent(binx, biny, sliceScale*hproj->GetBinContent(binx));
              hvars_2D[v]->SetBinError(binx, biny, sliceScale*hproj->GetBinError(binx));
            }
          }
          double finalIntegral = hvars_2D[v]->Integral();
          double integralScale = finalIntegral / initialIntegral;
          hvars_2D[v]->Scale(integralScale);
        }
      }

      cout << "Beginning scaling " << endl;

      double inputyield = 0, targetyield = 0, yieldScale = 1;
      inputyield = sum_signal;
      if (iProcess<4) targetyield = yield_signal_ggh[EnergyIndex][folder]*cutscale;
      else if (iProcess<8) targetyield = yield_signal_vbfh[EnergyIndex][folder]*cutscale;
      else if (iProcess<12) targetyield = yield_signal_wh[EnergyIndex][folder]*cutscale;
      else if (iProcess<16) targetyield = yield_signal_zh[EnergyIndex][folder]*cutscale;
      else if (iProcess<20) targetyield = yield_signal_tth[EnergyIndex][folder]*cutscale;
      else if (iProcess==20) targetyield = yield_signal_qqzz[EnergyIndex][folder]*cutscale;
      else if (iProcess==21) targetyield = yield_signal_ggzz[EnergyIndex][folder]*cutscale;
      else if (iProcess==22) targetyield = yield_signal_zx[EnergyIndex][folder]*cutscale;
      else if (iProcess==23) targetyield = inputyield;

      if (iProcess<23){
        yieldScale = (inputyield!=0 ? (targetyield/inputyield) : 0);
        for (int v = 0; v < nVariables; v++){
          hvars[v]->Scale(yieldScale);
        }
        for (int v = 0; v < nVariables_2D; v++){
          hvars_2D[v]->Scale(yieldScale);
        }
      }
      cout << "Overall yield scale: " << yieldScale << " == " << targetyield << " / " << inputyield << endl;

      // 1D recording
      for (int v = 0; v < nVariables; v++){
        foutput->WriteTObject(hvars[v]);
      }
      for (int v = 0; v < nVariables_2D; v++){
        bool notRecordable = false;
        if (notRecordable) continue;
        foutput->WriteTObject(hvars_2D[v]);
        string defaultName = hvars_2D[v]->GetName();
        if (isEnriched==0) hvars_2D[v]->SetName(Form("%s_%s_SR", processName[iProcess].c_str(), defaultName.c_str()));
        else hvars_2D[v]->SetName(Form("%s_%s_SR_SignalEnriched", processName[iProcess].c_str(), defaultName.c_str()));
        generic_Histo2DPlotter(foutput, plotDir_2D, hvars_2D[v], erg_tev, luminosity[EnergyIndex], folder, (iProcess==22));
        hvars_2D[v]->SetName(defaultName.c_str());
      }
      if (tc_extra!=0) delete tc_extra;
      delete tc;
      if (fit_pT!=0) delete fit_pT;
      finput_pTrewgt->Close();
      delete htemplate_Dbkg;
      foutput->Close();
      delete hRatio;
      fSSinput->Close();
    }

    for (int bs=0; bs<2; bs++) delete tBeam[bs];
  }
  for (int v = 0; v < nVariables_2D; v++){
    delete hvars_2D[v];
  }
  for (int v = 0; v < nVariables; v++){
    delete hvars[v];
  }
}

void plotCumulative_Unblinded(int isEnriched=0){
  int iRegion=0;

  gROOT->ProcessLine(".x tdrstyle.cc");
  float rangeZZMass[2] ={ systZZMass_range[iRegion][0], systZZMass_range[iRegion][1] };
  TH1F* hvars[nProcesses][nVariables] ={ { 0 } };

  TString cinput_common = user_dir_hep + "Analysis/Auxiliary/";

  for (int iProcess=0; iProcess<nProcesses; iProcess++){
    for (int erg_tev = 7; erg_tev < 9; erg_tev++){
      TString erg_dir;
      erg_dir.Form("LHC_%iTeV", erg_tev);
      TString comstring;
      comstring.Form("%iTeV", erg_tev);
      int EnergyIndex = 0; if (erg_tev == 8) EnergyIndex = 1;

      for (int folder = 0; folder < 3; folder++){
        TString INPUT_NAME = "LifetimeKD_Unblinded_";
        INPUT_NAME.Append(Form("%s", processName[iProcess].c_str()));
        INPUT_NAME.Append(Form("_%s_", user_folder[folder]));
        INPUT_NAME.Append(comstring);
        if (isEnriched==1) INPUT_NAME.Append("_SignalEnriched");
        INPUT_NAME.Append(".root");

        TString cinput = cinput_common + INPUT_NAME;
        TFile* finput = new TFile(cinput, "read");

        cout << "Read " << cinput << endl;
        if (finput!=0){
          if (!finput->IsZombie()){
            TString hname_core = "h_";
            for (int v = 0; v < nVariables; v++){
              TString hname = hname_core;
              hname.Append(varNames[v]);
              hname.Append("_");
              hname.Append(cutName);
              TH1F* htemp = (TH1F*)finput->Get(hname);
              if (htemp!=0){
                int h_index = iProcess;
                TString tempname = htemp->GetName();
                tempname.Append(processName[iProcess]);
                htemp->SetName(tempname);

                if (hvars[h_index][v]!=0) hvars[h_index][v]->Add(htemp);
                else{
                  gROOT->cd();
                  hvars[h_index][v] = (TH1F*)htemp->Clone(Form("%s_clone", htemp->GetName()));
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
  TString plotDir = cinput_common;
  plotDir.Append("../Plots/Unblinded/");
  TString plotDir_1D=plotDir;
  plotDir_1D.Append("Cumulative/");
  plotDir_1D.Append("1D/");
  TString mkdirCommand = "mkdir -p ";
  TString mkdirCommand_1D = mkdirCommand;
  mkdirCommand_1D.Append(plotDir_1D);
  gSystem->Exec(mkdirCommand_1D);

  for (int v = 0; v < nVariables; v++){
    bool canSkip = false;

    double max_plot = 0;

    for (int iProcess=0; iProcess<nProcesses; iProcess++){
      if (hvars[iProcess][v]==0) continue;
      hvars[iProcess][v]->GetXaxis()->SetLabelFont(42);
      hvars[iProcess][v]->GetXaxis()->SetLabelOffset(0.007);
      hvars[iProcess][v]->GetXaxis()->SetLabelSize(0.04);
      hvars[iProcess][v]->GetXaxis()->SetTitleSize(0.06);
      hvars[iProcess][v]->GetXaxis()->SetTitleOffset(0.9);
      hvars[iProcess][v]->GetXaxis()->SetTitleFont(42);
      hvars[iProcess][v]->GetYaxis()->SetNdivisions(505);
      hvars[iProcess][v]->GetYaxis()->SetLabelFont(42);
      hvars[iProcess][v]->GetYaxis()->SetLabelOffset(0.007);
      hvars[iProcess][v]->GetYaxis()->SetLabelSize(0.04);
      hvars[iProcess][v]->GetYaxis()->SetTitleSize(0.06);
      hvars[iProcess][v]->GetYaxis()->SetTitleOffset(1.1);
      hvars[iProcess][v]->GetYaxis()->SetTitleFont(42);

      if (v==vTxy_BS && !(iProcess>=22 || (iProcess<20 && iProcess % 4 !=0))) symmetrize_PromptTemplates(hvars[iProcess][v]);
    }

    TPaveText* pt = new TPaveText(0.15, 0.93, 0.85, 1, "brNDC");
    pt->SetBorderSize(0);
    pt->SetFillStyle(0);
    pt->SetTextAlign(12);
    pt->SetTextFont(42);
    pt->SetTextSize(0.045);
    TText* text = pt->AddText(0.025, 0.45, "#font[61]{CMS}");
    text->SetTextSize(0.044);
    if (isEnriched==0 && !(v==vDbkg) ){
      text = pt->AddText(0.165, 0.42, "#font[52]{Unpublished}");
      text->SetTextSize(0.0315);
    }
    TString cErgTev = "#font[42]{19.7 fb^{-1} (8 TeV) + 5.1 fb^{-1} (7 TeV)}";
    text = pt->AddText(0.537, 0.45, cErgTev);
    text->SetTextSize(0.0315);

    hvars[23][v]->SetLineColor(kBlack);
    hvars[23][v]->SetLineWidth(2);

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

    double maxplot=0;
    for (int bin = 1; bin <= hvars[23][v]->GetNbinsX(); bin++){
      double bincenter = hvars[23][v]->GetBinCenter(bin);
      double bincontent = hvars[23][v]->GetBinContent(bin);

      if (bincontent >= 0){
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
    appendName += cutName;
    appendName += "_SR_AllTeV";
    if (isEnriched==1) appendName.Append("_SignalEnriched");
    TString canvasname = "cCanvas_Unblinded_";
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

    int plotBSM = 0;
    if (!(v==vDbkg||v==vm4l||v==vpT4l)) plotBSM = 1;
    if (v==vpT4l) plotBSM = 3;
//    if (v==vTxy_BS) plotBSM = 3;
    double nPlotted_base = 5;
    double nPlotted = nPlotted_base;
    if (plotBSM>0 && plotBSM<=2) nPlotted += 1.;
    if (plotBSM>2) nPlotted += 2.;

    TLegend *ll;
    float lxmin = 0.22, lxwidth = 0.38;
    float lymax = 0.9, lywidth = 0.3;
    if (plotBSM>0) lywidth *= (nPlotted/nPlotted_base);
    if (v==vpT4l || v==vDbkg) lxmin += 0.39;
    float lxmax = lxmin + lxwidth;
    float lymin = lymax - lywidth;
    ll = new TLegend(lxmin, lymin, lxmax, lymax);
    ll->SetBorderSize(0);
    ll->SetTextFont(42);
    ll->SetTextSize(0.04);
    ll->SetLineColor(1);
    ll->SetLineStyle(1);
    ll->SetLineWidth(1);
    ll->SetFillColor(0);
    ll->SetFillStyle(0);

    ll->AddEntry(tgdata, "Observed", "ep");
    cout << "New cut \t | ZX norm: " << hvars[22][v]->Integral() << "\tggZZ norm: " << hvars[21][v]->Integral() << "\tqqZZ norm: " << hvars[20][v]->Integral() << endl;

    hvars[21][v]->Add(hvars[22][v]);
    hvars[20][v]->Add(hvars[21][v]);

    TH1F* hBSM = 0;
    TH1F* hVBFH = 0;
    TH1F* hTTH = 0;
    TH1F* hSM = 0;
    int useBSM = 1;

    if (hvars[0][v]!=0){
      hSM = (TH1F*)hvars[0][v]->Clone("hSM");
      for (int iProd=1; iProd<=4; iProd++) hSM->Add(hvars[iProd*4+0][v]);
    }
    if (hvars[useBSM][v]!=0){
      hBSM = (TH1F*)hvars[useBSM][v]->Clone("hBSM");
      for (int iProd=1; iProd<=4; iProd++) hBSM->Add(hvars[iProd*4+useBSM][v]);
    }
    if (hvars[4][v]!=0){
      hVBFH = (TH1F*)hvars[4][v]->Clone("hVBFH");
      hVBFH->Scale(hSM->Integral() / hVBFH->Integral());
    }
    if (hvars[16][v]!=0){
      hTTH = (TH1F*)hvars[16][v]->Clone("hTTH");
      hTTH->Scale(hSM->Integral() / hTTH->Integral());
    }

    if (hSM!=0){
      hSM->SetLineWidth(2);
      hSM->SetLineColor(kRed);
      hSM->SetFillColor(0);
      hSM->SetFillStyle(1001);
      cout << "Higgs norm: " << hSM->Integral() << endl;
      hSM->Add(hvars[20][v]);
      ll->AddEntry(hSM, processTitle[0], "f");
    }
    if (hBSM!=0 && plotBSM==1){
      hBSM->SetLineWidth(2);
      hBSM->SetLineStyle(7);
      hBSM->SetLineColor(kViolet);
      hBSM->SetFillColor(0);
      hBSM->SetFillStyle(1001);
      cout << "Higgs (BSM) norm: " << hBSM->Integral() << endl;
      hBSM->Add(hvars[20][v]);
      ll->AddEntry(hBSM, processTitle[useBSM], "f");
    }
    if (hVBFH!=0 && plotBSM==3){
      hVBFH->SetLineWidth(2);
      hVBFH->SetLineStyle(7);
      hVBFH->SetLineColor(kViolet);
      hVBFH->SetFillColor(0);
      hVBFH->SetFillStyle(1001);
      hVBFH->Add(hvars[20][v]);
      ll->AddEntry(hVBFH, "VBF signal", "f");
    }
    if (hTTH!=0 && (plotBSM==2 || plotBSM==3)){
      hTTH->SetLineWidth(2);
      hTTH->SetLineStyle(5);
      hTTH->SetLineColor(kGreen+3);
      hTTH->SetFillColor(0);
      hTTH->SetFillStyle(1001);
      hTTH->Add(hvars[20][v]);
      ll->AddEntry(hTTH, "t#bar{t}H signal", "f");
    }
    hvars[20][v]->SetLineWidth(2);
    hvars[21][v]->SetLineWidth(2);
    hvars[22][v]->SetLineWidth(2);
    hvars[20][v]->SetLineColor(kBlack);
    hvars[21][v]->SetLineColor(kBlack);
    hvars[22][v]->SetLineColor(kBlack);
    hvars[20][v]->SetFillColor(kAzure-9);
    hvars[21][v]->SetFillColor(kAzure-2);
    hvars[22][v]->SetFillColor(TColor::GetColor("#669966"));
    hvars[20][v]->SetFillStyle(1001);
    hvars[21][v]->SetFillStyle(1001);
    hvars[22][v]->SetFillStyle(1001);
    for (int pp=20; pp<=22; pp++) ll->AddEntry(hvars[pp][v], processTitle[pp-16], "f");

    hSM->GetXaxis()->SetNdivisions(505);
    hSM->GetYaxis()->SetNdivisions(505);
    TString xtitle = varTitles[v];
    TString ytitle = "Events / ";
    double binwidth = hSM->GetBinWidth(1);
    if (binwidth>=1){
      int intWidth = (int)binwidth;
      if (binwidth == (double)intWidth) ytitle.Append(Form("%.0f", binwidth));
      else ytitle.Append(Form("%.1f", binwidth));
    }
    else if (binwidth>=0.1) ytitle.Append(Form("%.1f", binwidth));
    else if (binwidth>=0.01) ytitle.Append(Form("%.2f", binwidth));
    else ytitle.Append(Form("%.3f", binwidth));
    if (v==vm4l || v==vpT4l) ytitle.Append(" GeV");
    if (v==vTxy_BS) ytitle.Append(" #mum");
    hSM->GetXaxis()->SetTitle(xtitle);
    hSM->GetYaxis()->SetTitle(ytitle);
    if (isEnriched==1 && v==vTxy_BS) hSM->GetYaxis()->SetRangeUser(0, maxplot*1.15);
    else hSM->GetYaxis()->SetRangeUser(0, maxplot*1.1);
//    hSM->GetYaxis()->SetRangeUser(0, maxplot*1.1);

    hSM->Draw("hist");
    if (hBSM!=0 && plotBSM==1) hBSM->Draw("histsame");
    if (hVBFH!=0 && plotBSM==3) hVBFH->Draw("histsame");
    if (hTTH!=0 && (plotBSM==2||plotBSM==3)) hTTH->Draw("histsame");
    hvars[20][v]->Draw("histsame");
    hvars[21][v]->Draw("histsame");
    hvars[22][v]->Draw("histsame");
    tgdata->Draw("e1psame");
    ll->Draw("same");
    pt->Draw();

    float pt_xmin = 0.71, pt_xwidth = 0.21;
    float pt_ymax = 0.92, pt_ywidth = 0.08;
    if (v==vpT4l){
      pt_xmin = 0.25;
    }
    float pt_xmax = pt_xmin + pt_xwidth;
    float pt_ymin = pt_ymax - pt_ywidth;
    TPaveText *pt10 = new TPaveText(pt_xmin, pt_ymin, pt_xmax, pt_ymax, "brNDC");
    pt10->SetBorderSize(0);
    pt10->SetTextAlign(12);
    pt10->SetTextSize(0.04);
    pt10->SetFillStyle(0);
    pt10->SetTextFont(42);
    TText* text10;
    if (isEnriched==1){
      text10 = pt10->AddText(0.01, 0.01, "D_{bkg} > 0.5");
//      pt10->Draw();
    }

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

    delete pt10;

    delete hBSM;
    delete hTTH;
    delete hVBFH;
    delete hSM;

    delete ll;
    cc->Close();
    delete tgdata;
    delete pt;
  }

  for (int hh=0; hh<nProcesses; hh++){
    for (int v = 0; v < nVariables; v++){
      if (hvars[hh][v]!=0) delete hvars[hh][v];
    }
  }
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

/*
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
*/
