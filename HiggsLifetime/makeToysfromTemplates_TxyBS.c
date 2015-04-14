#include <iostream>
#include <iomanip>
#include <cmath>
#include <string>
#include <vector>
#include <fstream>
#include <cstdlib>
#include <algorithm>
#include <cstdlib>
#include "TFile.h"
#include "TList.h"
#include "TRandom3.h"
#include "TTree.h"
#include "TChain.h"
#include "TString.h"
#include "TMath.h"
#include "TF1.h"
#include "TSpline.h"
#include "TSystem.h"
#include "TGraph.h"
#include "TGraphAsymmErrors.h"
#include "TMatrixD.h"
#include "TVectorD.h"
#include "TH1D.h"
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
#include "RooRealVar.h"
#include "RooArgSet.h"
#include "RooDataSet.h"
#include "RooCategory.h"
#include "RooPlot.h"
#include "RooMsgService.h"
#include "RooFitResult.h"
#include "./Pdfs/RooRealFlooredSumPdf.h"
#include "./data/ZZ4l_125_6_Samples.h"

using namespace std;
using namespace ROOT::Math;
using namespace RooFit;


const double ctau_limits[2] ={ 0, 1000 };
const int nctaus = 100;
const double cm_to_microns = 10000.;
const double KD_scale = 1250.;
const int N_Kinematics=1;
const int N_KDs=2;
double Kinematics_ranges[N_Kinematics][2] ={
  { 105.6, 140.6 }
};
double KD_ranges[N_KDs][2] ={
  { 0, 1 },
  { -1, 1 }
};
double Kinematics_bins[N_Kinematics] ={
  1000
};
double KD_bins[N_KDs] ={
  50,
  50
};


string processName[8] ={
  "CTau0", "CTau100", "CTau500", "CTau1000",
  "qqZZ", "ggZZ", "CR", "data"
};

const int nHiggsProd = 5;
TString strHiggsProductionFit[nHiggsProd-1] ={ "VBFH125", "WH125", "ZH125", "ttH125" };
TString strHiggsProduction[nHiggsProd-1] ={ "VBFH", "WH", "ZH", "ttH" };
TString strHiggsProduction_label[nHiggsProd] ={ "ggH", "VBF", "WH", "ZH", "t#bar{t}H" };

const int newSIPCuts = 4; // 0: No SIP; 1: Z1SIP, 2: KalmanChi2, 3: Z1SIP + KalmanChi2
TString cutNames[newSIPCuts + 1] ={
  "PVSIPcut",
  "NoSIPcut",
  "Z1SIPcut",
  "4lchi2cut",
  "Z1SIP_4lchi2cut"
};

TString templateDir = "/scratch0/hep/usarical/HiggsLifetime/Analysis/CMSSW_6_1_1/src/HiggsAnalysis/HZZ4l_Combination/CreateDatacards/templates2D/HCTau_NewSIP_BS/";

bool testSelection(
  float Lep1SIP, float Lep2SIP, float Lep3SIP, float Lep4SIP,
  float Lep1_Z1SIP, float Lep2_Z1SIP, float Lep3_Z1SIP, float Lep4_Z1SIP,
  float KalmanCandVtx_chi2,
  int scheme=4
  ){
  bool pass=true;
  switch (scheme){
  case 0:
    if (Lep1SIP>=4 || Lep2SIP>=4 || Lep3SIP>=4 || Lep4SIP>=4) pass = false;
    break;
  case 1:
    break;
  case 2:
    if (fabs(Lep1_Z1SIP) >= 4 || fabs(Lep2_Z1SIP) >= 4 || fabs(Lep3_Z1SIP) >= 5 || fabs(Lep4_Z1SIP) >= 5) pass = false;
    break;
  case 3:
    if (KalmanCandVtx_chi2>=30) pass = false;
    break;
  case 4:
    if (fabs(Lep1_Z1SIP) >= 4 || fabs(Lep2_Z1SIP) >= 4 || fabs(Lep3_Z1SIP) >= 5 || fabs(Lep4_Z1SIP) >= 5 || KalmanCandVtx_chi2>=30) pass = false;
    break;
  default:
    break;
  }
  return pass;
}
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
  int CRflag
  ){
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


float compute_delDxy(
  float ZZPhi,
  float CandVtx_cov_xx,
  float CandVtx_cov_xy,
  float CandVtx_cov_xz,
  float CandVtx_cov_yy,
  float CandVtx_cov_yz,
  float CandVtx_cov_zz,
  float PrimaryVtx_cov_xx,
  float PrimaryVtx_cov_xy,
  float PrimaryVtx_cov_xz,
  float PrimaryVtx_cov_yy,
  float PrimaryVtx_cov_yz,
  float PrimaryVtx_cov_zz
  );


void makeCombineSignalAsimovTrees_Signal_single(int folder, int erg_tev, int iProd, int BSMSample=0, int BSMNormScheme=0, int scheme=4);
void makeCombineSignalAsimovTrees_Bkg_single(int folder, int erg_tev, int scheme=4);
void assembleCombineEmbeddedToys(double target_ctau = 0, int BSMSample=0, int iProd = 0, double BSMNormScheme=0, int scheme=4);
void compareEmbeddedAsimov_SignalTemplates(double target_ctau, int iProd, int scheme=4);
int randomModulo(int i) { return rand()%i; }
TH1F* transform_tanh_varBin(TH1F* hinput);
void symmetrize_PromptTemplates(TH1F* hrepair);

void makeCombineAsimovTrees_NoSIP(int isBkg = 0, int selection=4){
  for (int CoM=7; CoM<9; ++CoM){
    for (int channel=0; channel<3; channel++){
      if (isBkg==0 || isBkg>1){
        for (int iProd=0; iProd<nHiggsProd; iProd++){
          for (int BSMSample=0; BSMSample<=kfLambda1_05; BSMSample++){
            for (int BSMNormScheme=0; BSMNormScheme<=1; BSMNormScheme++){
              makeCombineSignalAsimovTrees_Signal_single(channel, CoM, iProd, BSMSample, BSMNormScheme, selection);
            }
          }
        }
      }
      if(isBkg>=1) makeCombineSignalAsimovTrees_Bkg_single(channel, CoM, selection);
    }
  }
  if (abs(isBkg)>1){
    for (int iProd=0; iProd<nHiggsProd; iProd++){
      for (int BSMSample=0; BSMSample<=kfLambda1_05; BSMSample++){
        for (int BSMNormScheme=0; BSMNormScheme<=1; BSMNormScheme++){
          assembleCombineEmbeddedToys(0, BSMSample, iProd, BSMNormScheme, selection);
          if (iProd==0 && BSMSample==0 && BSMNormScheme==0){
            assembleCombineEmbeddedToys(100, BSMSample, iProd, BSMNormScheme, selection);
            assembleCombineEmbeddedToys(500, BSMSample, iProd, BSMNormScheme, selection);
            assembleCombineEmbeddedToys(1000, BSMSample, iProd, BSMNormScheme, selection);
          }
        }
      }
    }
  }
}

void makeCombineSignalAsimovTrees_Signal_single(int folder, int erg_tev, int iProd, int BSMSample, int BSMNormScheme, int scheme){
  if (iProd==0 && !(BSMSample<=kfLambda1_1 || BSMSample==kfLambda1_05)) return;
  else if (iProd==1 && BSMSample>=kNumFiles_VBFHVVBSM) return;
  else if (iProd==2 && BSMSample>=kNumFiles_WHVVBSM) return;
  else if (iProd==3 && (BSMSample>=kNumFiles_ZHVVBSM || BSMSample==kfg2_1_fg4_0)) return;
  else if (iProd==4 && (BSMSample>0)) return;
  if (BSMNormScheme!=0 && BSMSample==0) return;

  if (scheme==0) return;
  unsigned randomSeed = 3*erg_tev + folder;
  srand(randomSeed);

  TString comstring;
  comstring.Form("%iTeV", erg_tev);
  TString erg_dir;
  erg_dir.Form("LHC_%iTeV", erg_tev);
  int EnergyIndex=1;
  if (erg_tev==7) EnergyIndex=0;
  float ZZMass_PeakCut[2]={ 105.6, 140.6 };

  char TREE_NAME[]="SelectedTree";
  TString INPUT_NAME = "HZZ4lTree_powheg15jhuGenV3-CTau";
  TString OUTPUT_NAME = "AsimovTree_Sig_TxyBS_";
  if (iProd>0) OUTPUT_NAME.Append(Form("%s_", strHiggsProduction[iProd-1].Data()));
  if (BSMSample>0) OUTPUT_NAME.Append(Form("Hypo%i_", BSMSample));
  OUTPUT_NAME.Append(user_folder[folder]);
  OUTPUT_NAME.Append("_");
  OUTPUT_NAME.Append(comstring);
  OUTPUT_NAME.Append("_");
  OUTPUT_NAME = OUTPUT_NAME + cutNames[scheme];
  if (BSMNormScheme!=0) OUTPUT_NAME.Append("_SMLikeNorm");
  OUTPUT_NAME = OUTPUT_NAME + ".root";

  TString cinput_common = user_dir_hep + "No_SIP/" + erg_dir + "/" + user_folder[folder] + "/";
  TString coutput_common = user_dir_hep + "Analysis/AsimovToys/";

  // BSM HVV vertex block
  const int kNumFiles=kNumFiles_GGHVVBSM;
  const int kNumHypo=kNumGGHVVBSM;
  const int gapZZHypo=kfLambda1_03;
  const int gapZZHypo_2=kfLambda1_05_pL190;
  const int firstNonZZHypo=kfZG_1_fGG_0;
  TString shuffledTreeName = "HZZ4lTree_H125p6_ShuffledSignalBkg.root";
  int BSMFile=BSMSample;
  if (BSMSample<0 || BSMSample>=kNumHypo){
    cerr << "No hypothesis" << endl;
    assert(0);
  }
  else if (BSMSample>gapZZHypo && BSMSample<gapZZHypo_2){
    cout << "Passed gap hypothesis no. 1, " << gapZZHypo << ", while file range is still usable." << endl;
    BSMFile-=1;
  }
  else if (BSMSample>gapZZHypo_2 && BSMSample < (kNumFiles + 1)){
    cout << "Passed gap hypothesis no. 2, " << gapZZHypo << ", while file range is still usable." << endl;
    BSMFile-=2;
  }
  else if (BSMSample == gapZZHypo || BSMSample == gapZZHypo_2 || BSMSample >= (kNumFiles + 1)){
    BSMFile = 0;
    cout << "No dedicated sample is available." << endl;
    assert(0);
  }
  cout << "BSM Hypothesis is " << BSMSample << ", and  dedicated file is " << BSMFile << endl;

  double signalNormScale=1;
  if (BSMSample>0 && BSMNormScheme==0){
    char* cMCwgt[2] ={
      "MC_CV_weight",
      "MC_CV_4GeVcut_weight"
    };
    int chooseMCwgt=0;
    if (BSMSample>=firstNonZZHypo) chooseMCwgt=1;
    char cMCwgtchosen[30];
    sprintf(cMCwgtchosen, "%s", cMCwgt[chooseMCwgt]);
    char cMCwgtchosen_indexed[30];
    sprintf(cMCwgtchosen_indexed, "%s[%i]", cMCwgtchosen, BSMSample);
    cout << "Signal MC weight string: " << cMCwgtchosen_indexed << endl;
    char cBSMSampleAll_cut[100];
    sprintf(cBSMSampleAll_cut, "(ZZMass<%.1f && ZZMass>=%.1f)*%s", ZZMass_PeakCut[1], ZZMass_PeakCut[0], cMCwgtchosen_indexed);
    cout << "BSM cut for all events: " << cBSMSampleAll_cut << endl;
    //    char cBSMSampleDedicated_cut[100];
    //    sprintf(cBSMSampleDedicated_cut, "(ZZMass<%.1f && ZZMass>=%.1f && EventSample==%i)*%s", ZZMass_PeakCut[1], ZZMass_PeakCut[0], BSMFile, cMCwgtchosen_indexed);
    //    cout << "BSM cut for dedicated events: " << cBSMSampleDedicated_cut << endl;

    TChain *sig = new TChain("SelectedTree");
    TString cinput_ggBSM = cinput_common;
    cinput_ggBSM.Append(shuffledTreeName);
    sig->Add(cinput_ggBSM);

    TH1D *htempSM=new TH1D("htSM", "", 1, 0.0, 1000.);
    TH1D *htempSMCut=new TH1D("htSMCut", "", 1, 0.0, 1000.);
    TH1D *htempBSM = new TH1D("htBSM", "", 1, 0.0, 1000.);
    //    TH1D *htempBSMDedicated = new TH1D("htBSMDedicated", "", 1, 0.0, 1000.);

    sig->Draw("ZZMass>>htSM", Form("(ZZMass<%.1f && ZZMass>=%.1f)*%s[0]", ZZMass_PeakCut[1], ZZMass_PeakCut[0], cMCwgt[0]));
    sig->Draw("ZZMass>>htSMCut", Form("(ZZMass<%.1f && ZZMass>=%.1f)*%s[0]", ZZMass_PeakCut[1], ZZMass_PeakCut[0], cMCwgt[1]));
    sig->Draw("ZZMass>>htBSM", cBSMSampleAll_cut);
    //    sig->Draw("ZZMass>>htBSMDedicated", cBSMSampleDedicated_cut);

    double sumSM = htempSM->GetSumOfWeights();
    double sumSMCut = htempSMCut->GetSumOfWeights();
    double sumBSM = htempBSM->GetSumOfWeights();
    //    double sumBSMDedicated = htempBSMDedicated->GetSumOfWeights();//(EventSample==BSMFile)
    cout << "Rsignal for SM hypothesis all files (without 4 GeV cut): " << sumSM << endl;
    cout << "Rsignal for SM hypothesis all files (with 4 GeV cut): " << sumSMCut << endl;
    cout << "Rsignal for BSM hypothesis all files: " << sumBSM << endl;
    //    cout << "Rsignal for BSM hypothesis dedicated file: " << sumBSMDedicated << endl;
    if (chooseMCwgt==1) sumBSM *= (sumSM/sumSMCut);

    signalNormScale = (sumBSM/sumSM);
    cout << "Signal scale for BSM hypothesis: " << signalNormScale << endl;

    //    delete htempBSMDedicated;
    delete htempBSM;
    delete htempSMCut;
    delete htempSM;
    delete sig;
  }

  // Determine target yield
  TString INPUTYIELD_NAME = "LifetimeKD_RelativeSIPYields_";
  INPUTYIELD_NAME.Append(Form("%s", processName[0].c_str()));
  INPUTYIELD_NAME.Append(Form("_%s_", user_folder[folder]));
  INPUTYIELD_NAME.Append(comstring);
  INPUTYIELD_NAME.Append(Form("_m4l%.1f_%.1f", 105.6, 140.6));
  INPUTYIELD_NAME.Append(".root");
  TString cinput_yield_common = user_dir_hep + "Analysis/Auxiliary/";
  TString cinput_yield = cinput_yield_common + INPUTYIELD_NAME;
  TFile* finput = new TFile(cinput_yield, "read");
  TH1F* htemp = (TH1F*)finput->Get(Form("h%s", processName[0].c_str()));
  double ratio_targetsignalyield = htemp->GetBinContent(scheme+1)/htemp->GetBinContent(1);
  delete htemp;
  finput->Close();

  double sumXSEC = XSEC_ggH[EnergyIndex] + XSEC_VBFH[EnergyIndex] + XSEC_WH[EnergyIndex] + XSEC_ZH[EnergyIndex] + XSEC_ttH[EnergyIndex];
  //use ggH yields directly from SIP<4 case
  double targetYield = yield_signal_ggh[EnergyIndex][folder]*ratio_targetsignalyield * (sumXSEC/XSEC_ggH[EnergyIndex]);

  TFile* finput_pTrewgt = new TFile(Form("./data/compareProductionModes_pT_%iTeV%s", erg_tev, ".root"), "read");
  TF1* fit_pT = 0;
  if (iProd>0){
    TString fitname = strHiggsProductionFit[iProd-1];
    fitname.Append("_pToverMZZ_ratio_fit");
    fit_pT = (TF1*)finput_pTrewgt->Get(fitname);
    if (iProd==1) targetYield  = yield_signal_vbfh[EnergyIndex][folder]*ratio_targetsignalyield * (sumXSEC/XSEC_VBFH[EnergyIndex]);
    else if (iProd==2) targetYield  = yield_signal_wh[EnergyIndex][folder]*ratio_targetsignalyield * (sumXSEC/XSEC_WH[EnergyIndex]);
    else if (iProd==3) targetYield  = yield_signal_zh[EnergyIndex][folder]*ratio_targetsignalyield * (sumXSEC/XSEC_ZH[EnergyIndex]);
    else if (iProd==4) targetYield  = yield_signal_tth[EnergyIndex][folder]*ratio_targetsignalyield * (sumXSEC/XSEC_ttH[EnergyIndex]);
  }
  targetYield *= signalNormScale;
  cout << comstring << ' ' << user_folder[folder] << " yield: " << targetYield << endl;


  double templateWeight = 1.;
  int fileID;
  double ctauWeight=1;
  float GenHMass, GenHPt, GenHPhi;
  float ZZMass, ZZPt, ZZEta, ZZPhi;
  float MC_weight;
  float Lep1SIP, Lep2SIP, Lep3SIP, Lep4SIP;
  float Lep1_Z1SIP, Lep2_Z1SIP, Lep3_Z1SIP, Lep4_Z1SIP;
  float Z1CandVtx_x, Z1CandVtx_y, Z1CandVtx_z;

  float KalmanCandVtx_x, KalmanCandVtx_y, KalmanCandVtx_z, KalmanCandVtx_chi2;
  float OfflinePrimaryVtx_x, OfflinePrimaryVtx_y, OfflinePrimaryVtx_z;
  float GenPrimaryVtx_x, GenPrimaryVtx_y, GenPrimaryVtx_z;
  float GenIntVtx_x, GenIntVtx_y, GenIntVtx_z;
  float CandVtx_x, CandVtx_y, CandVtx_z;
  float TrueCandVtx_x, TrueCandVtx_y, TrueCandVtx_z;
  float CandBSVtx_x, CandBSVtx_y, CandBSVtx_z;
  float TrueCandBSVtx_x, TrueCandBSVtx_y, TrueCandBSVtx_z;

  int RunNumber, LumiNumber, RunNumber_Ref, LumiNumber_Ref;
  float BeamPosX, BeamPosY, BeamPosZ, BeamPosXErr, BeamPosYErr, BeamPosZErr;

  float KalmanCandVtx_cov_xx;
  float KalmanCandVtx_cov_xy;
  float KalmanCandVtx_cov_xz;
  float KalmanCandVtx_cov_yy;
  float KalmanCandVtx_cov_yz;
  float KalmanCandVtx_cov_zz;
  float OfflinePrimaryVtx_ndof;
  float OfflinePrimaryVtx_cov_xx;
  float OfflinePrimaryVtx_cov_xy;
  float OfflinePrimaryVtx_cov_xz;
  float OfflinePrimaryVtx_cov_yy;
  float OfflinePrimaryVtx_cov_yz;
  float OfflinePrimaryVtx_cov_zz;

  float p0plus_VAJHU;
  float p0plus_VAMCFM;
  float bkg_VAMCFM;
  float p0plus_m4l, p0plus_m4l_ScaleUp, p0plus_m4l_ScaleDown, p0plus_m4l_ResUp, p0plus_m4l_ResDown;
  float bkg_m4l, bkg_m4l_ScaleUp, bkg_m4l_ScaleDown, bkg_m4l_ResUp, bkg_m4l_ResDown;

  double KD=-99, KD_up=-99, KD_dn=-99, smd=-99, mzz;
  float Txy, Dxy, Txy_true, Dxy_true, delTxy, delDxy;
  float Txy_BS, Dxy_BS, Txy_BS_true, Dxy_BS_true;
  float delTxy_BS=-99;
  float sigmaPV_xy, sigmaInt_xy;

  float D_bkg=-99;
  float D_bkg_ScaleUp=-99;
  float D_bkg_ScaleDown=-99;
  float D_bkg_ResUp=-99;
  float D_bkg_ResDown=-99;

  double ctau_gridsize = (ctau_limits[1]-ctau_limits[0])/nctaus;
  int nbinsx = KD_bins[0];
  double xlow = -1, xhigh = 1;
  int nbinsy = KD_bins[1];
  double ylow=0, yhigh=1;

  TString coutput = coutput_common + OUTPUT_NAME;
  TFile* foutput = new TFile(coutput, "recreate");
  cout << "Starting to create " << coutput << endl;

  // BeamSpot info
  TString strBeamSpot[2]={ "MC_START44_V13", "MC_START53_V23" };
  TChain* tBeam = new TChain("BeamSpotRecord");
  TString cinput_BS = "./data/BeamSpotRecord_" + strBeamSpot[EnergyIndex] + "_" + comstring + ".root";
  tBeam->Add(cinput_BS);
  tBeam->SetBranchAddress("RunNumber", &RunNumber_Ref);
  tBeam->SetBranchAddress("LumiNumber", &LumiNumber_Ref);
  tBeam->SetBranchAddress("BeamPosX", &BeamPosX);
  tBeam->SetBranchAddress("BeamPosY", &BeamPosY);
  tBeam->SetBranchAddress("BeamPosZ", &BeamPosZ);
  tBeam->SetBranchAddress("BeamPosXErr", &BeamPosXErr);
  tBeam->SetBranchAddress("BeamPosYErr", &BeamPosYErr);
  tBeam->SetBranchAddress("BeamPosZErr", &BeamPosZErr);

  int ntrees = kNumSamples;
  if (BSMSample>0) ntrees=1;
  TChain* tree[kNumSamples] ={ 0 };
  for (int f=0; f<ntrees; f++){
    TString cinput = cinput_common;
    if (BSMSample==0 && (iProd==0 || iProd==4 || f>0)) cinput = cinput + sample_FullSim[f] + ".root";
    else if (iProd==0){
      cinput = cinput + sample_FullSim_GGHVVBSM[BSMFile] + ".root";
    }
    else if (iProd==1){
      cinput = cinput + sample_FullSim_VBFHVVBSM[BSMFile] + ".root";
    }
    else if (iProd==2){
      cinput = cinput + sample_FullSim_WHVVBSM[BSMFile] + ".root";
    }
    else if (iProd==3){
      cinput = cinput + sample_FullSim_ZHVVBSM[BSMFile] + ".root";
    }
    cout << "Attaching " << cinput << endl;
    tree[f] = new TChain(TREE_NAME);
    tree[f]->Add(cinput);

    tree[f]->SetBranchAddress("MC_weight", &MC_weight);
    tree[f]->SetBranchAddress("ZZMass", &ZZMass);
    tree[f]->SetBranchAddress("ZZPt", &ZZPt);
    tree[f]->SetBranchAddress("ZZEta", &ZZEta);
    tree[f]->SetBranchAddress("ZZPhi", &ZZPhi);
    tree[f]->SetBranchAddress("p0plus_VAJHU", &p0plus_VAJHU);
    tree[f]->SetBranchAddress("p0plus_VAMCFM", &p0plus_VAMCFM);
    tree[f]->SetBranchAddress("bkg_VAMCFM", &bkg_VAMCFM);
    tree[f]->SetBranchAddress("p0plus_m4l", &p0plus_m4l);
    tree[f]->SetBranchAddress("bkg_m4l", &bkg_m4l);
    tree[f]->SetBranchAddress("p0plus_m4l_ScaleUp", &p0plus_m4l_ScaleUp);
    tree[f]->SetBranchAddress("bkg_m4l_ScaleUp", &bkg_m4l_ScaleUp);
    tree[f]->SetBranchAddress("p0plus_m4l_ScaleDown", &p0plus_m4l_ScaleDown);
    tree[f]->SetBranchAddress("bkg_m4l_ScaleDown", &bkg_m4l_ScaleDown);
    tree[f]->SetBranchAddress("p0plus_m4l_ResUp", &p0plus_m4l_ResUp);
    tree[f]->SetBranchAddress("bkg_m4l_ResUp", &bkg_m4l_ResUp);
    tree[f]->SetBranchAddress("p0plus_m4l_ResDown", &p0plus_m4l_ResDown);
    tree[f]->SetBranchAddress("bkg_m4l_ResDown", &bkg_m4l_ResDown);

    tree[f]->SetBranchAddress("Lep1SIP", &Lep1SIP);
    tree[f]->SetBranchAddress("Lep2SIP", &Lep2SIP);
    tree[f]->SetBranchAddress("Lep3SIP", &Lep3SIP);
    tree[f]->SetBranchAddress("Lep4SIP", &Lep4SIP);
    tree[f]->SetBranchAddress("Lep1_Z1SIP", &Lep1_Z1SIP);
    tree[f]->SetBranchAddress("Lep2_Z1SIP", &Lep2_Z1SIP);
    tree[f]->SetBranchAddress("Lep3_Z1SIP", &Lep3_Z1SIP);
    tree[f]->SetBranchAddress("Lep4_Z1SIP", &Lep4_Z1SIP);
    tree[f]->SetBranchAddress("OfflinePrimaryVtx_x", &OfflinePrimaryVtx_x);
    tree[f]->SetBranchAddress("OfflinePrimaryVtx_y", &OfflinePrimaryVtx_y);
    tree[f]->SetBranchAddress("OfflinePrimaryVtx_z", &OfflinePrimaryVtx_z);
    tree[f]->SetBranchAddress("OfflinePrimaryVtx_cov_xx", &OfflinePrimaryVtx_cov_xx);
    tree[f]->SetBranchAddress("OfflinePrimaryVtx_cov_xy", &OfflinePrimaryVtx_cov_xy);
    tree[f]->SetBranchAddress("OfflinePrimaryVtx_cov_xz", &OfflinePrimaryVtx_cov_xz);
    tree[f]->SetBranchAddress("OfflinePrimaryVtx_cov_yy", &OfflinePrimaryVtx_cov_yy);
    tree[f]->SetBranchAddress("OfflinePrimaryVtx_cov_yz", &OfflinePrimaryVtx_cov_yz);
    tree[f]->SetBranchAddress("OfflinePrimaryVtx_cov_zz", &OfflinePrimaryVtx_cov_zz);
    tree[f]->SetBranchAddress("OfflinePrimaryVtx_ndof", &OfflinePrimaryVtx_ndof);
    tree[f]->SetBranchAddress("KalmanCandVtx_x", &KalmanCandVtx_x);
    tree[f]->SetBranchAddress("KalmanCandVtx_y", &KalmanCandVtx_y);
    tree[f]->SetBranchAddress("KalmanCandVtx_z", &KalmanCandVtx_z);
    tree[f]->SetBranchAddress("KalmanCandVtx_chi2", &KalmanCandVtx_chi2);
    tree[f]->SetBranchAddress("KalmanCandVtx_cov_xx", &KalmanCandVtx_cov_xx);
    tree[f]->SetBranchAddress("KalmanCandVtx_cov_xy", &KalmanCandVtx_cov_xy);
    tree[f]->SetBranchAddress("KalmanCandVtx_cov_xz", &KalmanCandVtx_cov_xz);
    tree[f]->SetBranchAddress("KalmanCandVtx_cov_yy", &KalmanCandVtx_cov_yy);
    tree[f]->SetBranchAddress("KalmanCandVtx_cov_yz", &KalmanCandVtx_cov_yz);
    tree[f]->SetBranchAddress("KalmanCandVtx_cov_zz", &KalmanCandVtx_cov_zz);
    tree[f]->SetBranchAddress("Z1CandVtx_x", &Z1CandVtx_x);
    tree[f]->SetBranchAddress("Z1CandVtx_y", &Z1CandVtx_y);
    tree[f]->SetBranchAddress("Z1CandVtx_z", &Z1CandVtx_z);

    if (tree[f]->GetBranchStatus("GenPrimaryVtx_x")){
      tree[f]->SetBranchAddress("GenPrimaryVtx_x", &GenPrimaryVtx_x);
      tree[f]->SetBranchAddress("GenPrimaryVtx_y", &GenPrimaryVtx_y);
      tree[f]->SetBranchAddress("GenPrimaryVtx_z", &GenPrimaryVtx_z);
      tree[f]->SetBranchAddress("GenIntVtx_x", &GenIntVtx_x);
      tree[f]->SetBranchAddress("GenIntVtx_y", &GenIntVtx_y);
      tree[f]->SetBranchAddress("GenIntVtx_z", &GenIntVtx_z);
      tree[f]->SetBranchAddress("GenHPt", &GenHPt);
      tree[f]->SetBranchAddress("GenHPhi", &GenHPhi);
      tree[f]->SetBranchAddress("GenHMass", &GenHMass);
    }
  }

  for (int genctau = 0; genctau <= nctaus; genctau++){
    double target_ctau = ctau_limits[0] + ctau_gridsize*genctau;
    bool skipTargetCTau = true;
    for (int f = 0; f < ntrees; f++){
      double initial_ctau = sample_CTau[f];
      if (initial_ctau == target_ctau){ skipTargetCTau = false; break; }
    }
    if (skipTargetCTau) continue;

    TH2F* hfill;
    hfill = new TH2F(Form("H_2D_CTau%.0f", target_ctau), "", nbinsx, xlow, xhigh, nbinsy, ylow, yhigh);
    hfill->SetTitle(Form("%.0f #mum %s at %i TeV", target_ctau, user_folder[folder], erg_tev));
    hfill->GetXaxis()->SetTitle(Form("tanh(T_{xy} / %.0f mm)", KD_scale));
    hfill->GetYaxis()->SetTitle("#it{D}_{bkg}");

    TTree* mytree = new TTree(Form("T_2D_CTau%.0f", target_ctau), "");
    mytree->Branch("_weight_", &templateWeight);
    mytree->Branch("CMS_zz4l_smd", &smd);
    mytree->Branch("CMS_zz4l_KD", &KD);
    mytree->Branch("CMS_zz4l_mass", &mzz);
    mytree->Branch("GenCTau", &Txy_true);
    mytree->Branch("fileID", &fileID);

    for (int f = 0; f < ntrees; f++){
      if (f == 0 && genctau>0) continue;
      if (f > 0 && genctau==0) continue;
      double initial_ctau = sample_CTau[f];
      if (initial_ctau!=target_ctau) continue;
      double effectiveInvCTau = (f>0 ? 1. / target_ctau - 1. / initial_ctau : 0);
      //			double xsec_ctau_wgt = hCounters->GetBinContent(f+1, 1) / hCounters->GetBinContent(f+1, 2);
      double xsec_ctau_wgt = (initial_ctau>0 ? initial_ctau / target_ctau : 1);
      cout << "CTau xsec reweighting from " << initial_ctau << " to " << target_ctau << " = " << xsec_ctau_wgt << endl;

      for (int ev = 0; ev < tree[f]->GetEntries(); ev++){
        RunNumber = 1;
        LumiNumber = 1;

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

        tree[f]->GetEntry(ev);
        if (!(ZZMass >= ZZMass_PeakCut[0] && ZZMass < ZZMass_PeakCut[1])) continue;
        bool passSelection = testSelection(
          Lep1SIP, Lep2SIP, Lep3SIP, Lep4SIP,
          Lep1_Z1SIP, Lep2_Z1SIP, Lep3_Z1SIP, Lep4_Z1SIP,
          KalmanCandVtx_chi2,
          scheme
          );
        if (!passSelection) continue;

        bool matchFound = false;
        for (int evBS=0; evBS < tBeam->GetEntries(); evBS++){
          tBeam->GetEntry(evBS);
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

        TrueCandBSVtx_x = GenIntVtx_x - BeamPosX;
        TrueCandBSVtx_y = GenIntVtx_y - BeamPosY;
        TrueCandBSVtx_z = GenIntVtx_z - BeamPosZ;

        Dxy_BS_true = (TrueCandBSVtx_x*cos(GenHPhi) + TrueCandBSVtx_y*sin(GenHPhi)) * cm_to_microns;
        Txy_BS_true = Dxy_BS_true*GenHMass / GenHPt;

        Dxy_BS = (CandBSVtx_x*cos(ZZPhi) + CandBSVtx_y*sin(ZZPhi)) * cm_to_microns;
        Txy_BS = Dxy_BS*ZZMass / ZZPt;

        CandVtx_x = KalmanCandVtx_x - OfflinePrimaryVtx_x;
        CandVtx_y = KalmanCandVtx_y - OfflinePrimaryVtx_y;
        CandVtx_z = KalmanCandVtx_z - OfflinePrimaryVtx_z;

        TrueCandVtx_x = GenIntVtx_x - GenPrimaryVtx_x;
        TrueCandVtx_y = GenIntVtx_y - GenPrimaryVtx_y;
        TrueCandVtx_z = GenIntVtx_z - GenPrimaryVtx_z;

        Dxy_true = (TrueCandVtx_x*cos(GenHPhi) + TrueCandVtx_y*sin(GenHPhi)) * cm_to_microns;
        Txy_true = Dxy_true*GenHMass / GenHPt;
        if (Txy_true<0) continue;

        D_bkg = p0plus_VAJHU*p0plus_m4l / (p0plus_VAJHU*p0plus_m4l + bkg_VAMCFM*bkg_m4l);
        D_bkg_ScaleUp = (p0plus_VAJHU*p0plus_m4l_ScaleUp) / (p0plus_VAJHU*p0plus_m4l_ScaleUp + bkg_VAMCFM*bkg_m4l_ScaleUp);
        D_bkg_ScaleDown = (p0plus_VAJHU*p0plus_m4l_ScaleDown) / (p0plus_VAJHU*p0plus_m4l_ScaleDown + bkg_VAMCFM*bkg_m4l_ScaleDown);
        D_bkg_ResUp = (p0plus_VAJHU*p0plus_m4l_ResUp) / (p0plus_VAJHU*p0plus_m4l_ResUp + bkg_VAMCFM*bkg_m4l_ResUp);
        D_bkg_ResDown = (p0plus_VAJHU*p0plus_m4l_ResDown) / (p0plus_VAJHU*p0plus_m4l_ResDown + bkg_VAMCFM*bkg_m4l_ResDown);

        mzz = ZZMass;
        smd = D_bkg;
        KD = Txy_BS;
        KD = tanh(KD / KD_scale); // Scaled to due to tanh function
        if (KD >= 1.) KD = 1. - 1.0e-4;
        if (KD <= -1.) KD = -1. + 1.0e-4;
        templateWeight = MC_weight;
        ctauWeight = exp(-Txy_true*effectiveInvCTau);
        templateWeight *= ctauWeight*xsec_ctau_wgt;
        if (fit_pT!=0 && BSMSample==0 && ((iProd>0 && iProd<4 && f>0) || iProd==4)) templateWeight *= fit_pT->Eval(GenHPt/GenHMass);

        if (smd<0) templateWeight = 0;
        if (templateWeight<=0) continue;
        hfill->Fill(KD, smd, templateWeight);

      }
    }
    double integral_rescale = 1;
    hfill->SetOption("colz");
    double integral_rewgt = hfill->Integral(0, nbinsx+1, 0, nbinsy+1);
    double integral_target = targetYield;
    integral_rescale = integral_target / integral_rewgt;
    hfill->Scale(integral_rescale);

    double sum_filled=0;
    for (int f = 0; f < ntrees; f++){
      if (f == 0 && genctau>0) continue;
      if (f > 0 && genctau==0) continue;
      double initial_ctau = sample_CTau[f];
      if (initial_ctau!=target_ctau) continue;
      double effectiveInvCTau = (f>0 ? 1. / target_ctau - 1. / initial_ctau : 0);
      //			double xsec_ctau_wgt = hCounters->GetBinContent(f+1, 1) / hCounters->GetBinContent(f+1, 2);
      double xsec_ctau_wgt = (initial_ctau>0 ? initial_ctau / target_ctau : 1);
      fileID = f;

      cout << "Starting to shuffle the entries... ";
      vector<int> treeEntries;
      for (int ev = 0; ev < tree[f]->GetEntries(); ev++) treeEntries.push_back(ev);
      random_shuffle(treeEntries.begin(), treeEntries.end(), randomModulo);
      cout << "done!" << endl;
      cout << treeEntries.size() << " / " << tree[f]->GetEntries() << endl;

      for (int ev = 0; ev < tree[f]->GetEntries(); ev++){
        RunNumber = 1;
        LumiNumber = 1;

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

        int entryNumber = ev;
        //        int entryNumber = (int)treeEntries[ev];
        //        cout << "Order " << ev << " -> " << entryNumber << endl;
        tree[f]->GetEntry(entryNumber);
        if (!(ZZMass >= ZZMass_PeakCut[0] && ZZMass < ZZMass_PeakCut[1])) continue;
        bool passSelection = testSelection(
          Lep1SIP, Lep2SIP, Lep3SIP, Lep4SIP,
          Lep1_Z1SIP, Lep2_Z1SIP, Lep3_Z1SIP, Lep4_Z1SIP,
          KalmanCandVtx_chi2,
          scheme
          );
        if (!passSelection) continue;

        bool matchFound = false;
        for (int evBS=0; evBS < tBeam->GetEntries(); evBS++){
          tBeam->GetEntry(evBS);
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

        TrueCandBSVtx_x = GenIntVtx_x - BeamPosX;
        TrueCandBSVtx_y = GenIntVtx_y - BeamPosY;
        TrueCandBSVtx_z = GenIntVtx_z - BeamPosZ;

        Dxy_BS_true = (TrueCandBSVtx_x*cos(GenHPhi) + TrueCandBSVtx_y*sin(GenHPhi)) * cm_to_microns;
        Txy_BS_true = Dxy_BS_true*GenHMass / GenHPt;

        Dxy_BS = (CandBSVtx_x*cos(ZZPhi) + CandBSVtx_y*sin(ZZPhi)) * cm_to_microns;
        Txy_BS = Dxy_BS*ZZMass / ZZPt;

        CandVtx_x = KalmanCandVtx_x - OfflinePrimaryVtx_x;
        CandVtx_y = KalmanCandVtx_y - OfflinePrimaryVtx_y;
        CandVtx_z = KalmanCandVtx_z - OfflinePrimaryVtx_z;

        TrueCandVtx_x = GenIntVtx_x - GenPrimaryVtx_x;
        TrueCandVtx_y = GenIntVtx_y - GenPrimaryVtx_y;
        TrueCandVtx_z = GenIntVtx_z - GenPrimaryVtx_z;

        Dxy_true = (TrueCandVtx_x*cos(GenHPhi) + TrueCandVtx_y*sin(GenHPhi)) * cm_to_microns;
        Txy_true = Dxy_true*GenHMass / GenHPt;
        if (Txy_true<0) continue;

        D_bkg = p0plus_VAJHU*p0plus_m4l / (p0plus_VAJHU*p0plus_m4l + bkg_VAMCFM*bkg_m4l);
        D_bkg_ScaleUp = (p0plus_VAJHU*p0plus_m4l_ScaleUp) / (p0plus_VAJHU*p0plus_m4l_ScaleUp + bkg_VAMCFM*bkg_m4l_ScaleUp);
        D_bkg_ScaleDown = (p0plus_VAJHU*p0plus_m4l_ScaleDown) / (p0plus_VAJHU*p0plus_m4l_ScaleDown + bkg_VAMCFM*bkg_m4l_ScaleDown);
        D_bkg_ResUp = (p0plus_VAJHU*p0plus_m4l_ResUp) / (p0plus_VAJHU*p0plus_m4l_ResUp + bkg_VAMCFM*bkg_m4l_ResUp);
        D_bkg_ResDown = (p0plus_VAJHU*p0plus_m4l_ResDown) / (p0plus_VAJHU*p0plus_m4l_ResDown + bkg_VAMCFM*bkg_m4l_ResDown);

        mzz = ZZMass;
        smd = D_bkg;
        KD = Txy_BS;
        KD = tanh(KD / KD_scale); // Scaled to due to tanh function
        if (KD >= 1.) KD = 1. - 1.0e-4;
        if (KD <= -1.) KD = -1. + 1.0e-4;

        templateWeight = MC_weight;
        ctauWeight = exp(-Txy_true*effectiveInvCTau);
        templateWeight *= ctauWeight*xsec_ctau_wgt;
        templateWeight *= integral_rescale;
        if (fit_pT!=0 && BSMSample==0 && ((iProd>0 && iProd<4 && f>0) || iProd==4)) templateWeight *= fit_pT->Eval(GenHPt/GenHMass);
        if (smd<0) templateWeight = 0;
        if (templateWeight<=0) continue;
        sum_filled += templateWeight;

        mytree->Fill();
      }
      cout << "Filled " << mytree->GetName() << ": " << mytree->GetEntries() << endl;
    }
    cout << "Final recorded yield: " << sum_filled << " | ";
    cout << endl;
    cout << "Filled the tree" << endl;
    foutput->WriteTObject(mytree);

    delete mytree;
    delete hfill;
  }
  for (int f=0; f<kNumSamples; f++){ if (tree[f]!=0) delete tree[f]; }

  delete tBeam;
  foutput->Close();
  cout << "Closed " << coutput << endl;
  if (fit_pT!=0) delete fit_pT;
  finput_pTrewgt->Close();
}

void makeCombineSignalAsimovTrees_Bkg_single(int folder, int erg_tev, int scheme){
  unsigned randomSeed = 3*erg_tev + folder;
  srand(randomSeed);

  TString comstring;
  comstring.Form("%iTeV", erg_tev);
  TString erg_dir;
  erg_dir.Form("LHC_%iTeV", erg_tev);
  int EnergyIndex=0; if (erg_tev==8) EnergyIndex=1;

  char TREE_NAME[]="SelectedTree";
  TString OUTPUT_NAME = "AsimovTree_background_TxyBS_";
  OUTPUT_NAME.Append(user_folder[folder]);
  OUTPUT_NAME.Append("_");
  OUTPUT_NAME.Append(comstring);
  OUTPUT_NAME =  OUTPUT_NAME + "_" + cutNames[scheme] + ".root";

  float ZZMass_PeakCut[2]={ 105.6, 140.6 };
  const int ntrees=3;
  double initialYield[ntrees]={
    yield_signal_qqzz[EnergyIndex][folder],
    yield_signal_ggzz[EnergyIndex][folder],
    yield_signal_zx[EnergyIndex][folder]
  };
  double targetYield[ntrees]={ 0 };
  double ratio_targetsignalyield[ntrees][5] ={ { 0 } };
  double targetCRcount[2][5] ={ { 0 } };

  // Determine target yield
  for (int iProcess = 4; iProcess<4+ntrees; iProcess++){
    TString INPUTYIELD_NAME = "LifetimeKD_RelativeSIPYields_";
    INPUTYIELD_NAME.Append(Form("%s", processName[iProcess].c_str()));
    INPUTYIELD_NAME.Append(Form("_%s_", user_folder[folder]));
    INPUTYIELD_NAME.Append(comstring);
    INPUTYIELD_NAME.Append(Form("_m4l%.1f_%.1f", 105.6, 140.6));
    INPUTYIELD_NAME.Append(".root");
    TString cinput_yield_common = user_dir_hep + "Analysis/Auxiliary/";
    TString cinput_yield = cinput_yield_common + INPUTYIELD_NAME;
    TFile* finput = new TFile(cinput_yield, "read");
    TH1F* htemp;
    if (iProcess!=6){
      htemp = (TH1F*)finput->Get(Form("h%s", processName[iProcess].c_str()));
      for (int b = 1; b <= 5; b++) ratio_targetsignalyield[iProcess-4][b-1] = htemp->GetBinContent(b);
      delete htemp;
    }
    else{
      htemp = (TH1F*)finput->Get(Form("h%s_OS", processName[iProcess].c_str()));
      for (int b = 1; b <= 5; b++) ratio_targetsignalyield[iProcess-4][b-1] = htemp->GetBinContent(b);
      delete htemp;
      htemp = (TH1F*)finput->Get(Form("h%s_SS", processName[iProcess].c_str()));
      for (int b = 1; b <= 5; b++) ratio_targetsignalyield[iProcess-4][b-1] += htemp->GetBinContent(b);
      delete htemp;

      htemp = (TH1F*)finput->Get(Form("h%s_Unweighted_OS", processName[iProcess].c_str()));
      for (int b = 1; b <= 5; b++) targetCRcount[0][b-1] = htemp->GetBinContent(b);
      delete htemp;
      htemp = (TH1F*)finput->Get(Form("h%s_Unweighted_SS", processName[iProcess].c_str()));
      for (int b = 1; b <= 5; b++) targetCRcount[1][b-1] = htemp->GetBinContent(b);
      delete htemp;
    }
    finput->Close();
    targetYield[iProcess-4] = initialYield[iProcess-4]*ratio_targetsignalyield[iProcess-4][scheme]/ratio_targetsignalyield[iProcess-4][0];
    cout << comstring << ' ' << user_folder[folder] << ' ' << processName[iProcess] <<  " yield: " << targetYield[iProcess-4] << endl;
  }
  double scale_yield[ntrees] ={ 1, 1, 1 };

  string cinput_aux = user_dir_hep + "Analysis/Auxiliary/";
  TString SSinputname = "OSoverSS_";
  SSinputname.Append("MCCR_");
  SSinputname += comstring;
  SSinputname.Append(".root");
  SSinputname.Prepend(cinput_aux.c_str());
  TFile* fSSinput = new TFile(SSinputname, "read");
  TH2F* hRatio;
  TString chratio = "hCR_MC_OSSSRatios_";
  if (scheme==0) chratio.Append("Old cut");
  else if (scheme==1) chratio.Append("No cut");
  else if (scheme==2) chratio.Append("Z1-SIP");
  else if (scheme==3) chratio.Append("chi**2");
  else if (scheme==4) chratio.Append("New cut");
  hRatio = (TH2F*)fSSinput->Get(chratio);
  cout << "Obtained ratio histogram " << hRatio->GetName() << endl;

  double templateWeight=1;
  float ZZMass, ZZPt, ZZEta, ZZPhi;
  float MC_weight;
  float MC_weight_noxsec;
  float ZXfake_weight_SS;
  float ZXfake_weight_OS[6];
  float MC_weight_QQZZEWK = 1;
  float MC_weight_Kfactor = 1;
  float MC_weight_QQBGGProper[4];
  float Lep1_Z1SIP, Lep2_Z1SIP, Lep3_Z1SIP, Lep4_Z1SIP;
  float Lep1SIP, Lep2SIP, Lep3SIP, Lep4SIP;
  int Lep1ID, Lep2ID, Lep3ID, Lep4ID;


  float Z1CandVtx_x, Z1CandVtx_y, Z1CandVtx_z;
  float KalmanCandVtx_x, KalmanCandVtx_y, KalmanCandVtx_z, KalmanCandVtx_chi2;
  float OfflinePrimaryVtx_x, OfflinePrimaryVtx_y, OfflinePrimaryVtx_z;
  float GenHMass, GenHPt, GenHPhi;
  float GenPrimaryVtx_x, GenPrimaryVtx_y, GenPrimaryVtx_z;
  float GenIntVtx_x, GenIntVtx_y, GenIntVtx_z;
  float CandVtx_x, CandVtx_y, CandVtx_z;
  float TrueCandVtx_x, TrueCandVtx_y, TrueCandVtx_z;
  float CandBSVtx_x, CandBSVtx_y, CandBSVtx_z;
  float TrueCandBSVtx_x, TrueCandBSVtx_y, TrueCandBSVtx_z;

  int RunNumber, LumiNumber, RunNumber_Ref, LumiNumber_Ref;
  float BeamPosX, BeamPosY, BeamPosZ, BeamPosXErr, BeamPosYErr, BeamPosZErr;

  int CRflag=-1;
  float Lep1combRelIsoPF=0, Lep2combRelIsoPF=0, Lep3combRelIsoPF=0, Lep4combRelIsoPF=0;
  bool Lep3isID=true, Lep4isID=true;
  short Z1ids;

  float KalmanCandVtx_cov_xx;
  float KalmanCandVtx_cov_xy;
  float KalmanCandVtx_cov_xz;
  float KalmanCandVtx_cov_yy;
  float KalmanCandVtx_cov_yz;
  float KalmanCandVtx_cov_zz;
  float OfflinePrimaryVtx_ndof;
  float OfflinePrimaryVtx_cov_xx;
  float OfflinePrimaryVtx_cov_xy;
  float OfflinePrimaryVtx_cov_xz;
  float OfflinePrimaryVtx_cov_yy;
  float OfflinePrimaryVtx_cov_yz;
  float OfflinePrimaryVtx_cov_zz;

  float p0plus_VAJHU;
  float p0plus_VAMCFM;
  float bkg_VAMCFM;
  float p0plus_m4l, p0plus_m4l_ScaleUp, p0plus_m4l_ScaleDown, p0plus_m4l_ResUp, p0plus_m4l_ResDown;
  float bkg_m4l, bkg_m4l_ScaleUp, bkg_m4l_ScaleDown, bkg_m4l_ResUp, bkg_m4l_ResDown;

  double KD=-99, smd=-99, mzz=-99;
  float Txy, Dxy, Txy_true, Dxy_true, delTxy, delDxy;
  float Txy_BS, Dxy_BS, Txy_BS_true, Dxy_BS_true;
  float delTxy_BS=-99;
  float sigmaPV_xy, sigmaInt_xy;
  float D_bkg=-99;

  TString cinput_common = user_dir_hep + "No_SIP/" + erg_dir + "/";
  TString coutput_common = user_dir_hep + "Analysis/AsimovToys/";
  TString coutput = coutput_common + OUTPUT_NAME;
  TFile* foutput = new TFile(coutput, "recreate");

  TChain* tBeam[2];
  TString strBeamSpot[2][2]={
    { "data_GR_R_44_V15C", "MC_START44_V13" },
    { "data_FT_53_V21_AN4", "MC_START53_V23" }
  };
  for (int bs=0; bs<2; bs++){
    tBeam[bs] = new TChain("BeamSpotRecord");
    if (bs==1){
      TString cinput_BS = "./data/BeamSpotRecord_" + strBeamSpot[EnergyIndex][bs] + "_" + comstring + ".root";
      tBeam[bs]->Add(cinput_BS);
    }
    else{
      for (int ee=0; ee<2; ee++){
        TString comstring_temp = "TeV";
        if (ee==0) comstring_temp.Prepend(Form("%i", 7));
        else if (ee==1) comstring_temp.Prepend(Form("%i", 8));
        TString cinput_BS = "./data/BeamSpotRecord_" + strBeamSpot[ee][bs] + "_" + comstring_temp + ".root";
        cout << "Beamspot record for data: " << cinput_BS << endl;
        tBeam[bs]->Add(cinput_BS);
      }
    }
    tBeam[bs]->SetBranchAddress("RunNumber", &RunNumber_Ref);
    tBeam[bs]->SetBranchAddress("LumiNumber", &LumiNumber_Ref);
    tBeam[bs]->SetBranchAddress("BeamPosX", &BeamPosX);
    tBeam[bs]->SetBranchAddress("BeamPosY", &BeamPosY);
    tBeam[bs]->SetBranchAddress("BeamPosZ", &BeamPosZ);
    tBeam[bs]->SetBranchAddress("BeamPosXErr", &BeamPosXErr);
    tBeam[bs]->SetBranchAddress("BeamPosYErr", &BeamPosYErr);
    tBeam[bs]->SetBranchAddress("BeamPosZErr", &BeamPosZErr);
  }

  TString cinput_qqzz_common = cinput_common + user_folder[folder] + "/";
  TString cinput_ggzz_common = cinput_qqzz_common;
  TString cinput_zx_common = cinput_common;

  const int nQQZZBkg = kQQBZZSamples-kGGSamples;
  const int nGGZZBkg = kGGMCFMSamples-kGGOLDSamples;
  const int nZXBkg = 2;
  const int bkgSize[ntrees]={ nQQZZBkg, nGGZZBkg, nZXBkg };
  TChain* tqqzz[nQQZZBkg];
  for (int smp = kGGSamples; smp < kQQBZZSamples; smp++){
    tqqzz[smp-kGGSamples] = new TChain(TREE_NAME);
    TString cinput_qqzz = cinput_qqzz_common + sample_FullSim[smp] + "_Reprocessed.root";
    tqqzz[smp-kGGSamples]->Add(cinput_qqzz);
  }
  TChain* tggzz[nGGZZBkg];
  for (int smp = kGGOLDSamples; smp < kGGMCFMSamples; smp++){
    tggzz[smp-kGGOLDSamples] = new TChain(TREE_NAME);
    TString cinput_ggzz = cinput_qqzz_common + sample_FullSim[smp] + "_Reprocessed.root";
    tggzz[smp-kGGOLDSamples]->Add(cinput_ggzz);
  }
  TChain* tzx[nZXBkg];
  for (int ee=0; ee<2; ee++){
    tzx[ee] = new TChain(TREE_NAME);
    int ee_tev = 7;
    if (ee==1) ee_tev=8;
    TString cinput_zx = user_dir_hep;
    cinput_zx.Append("No_SIP/");
    cinput_zx.Append(Form("LHC_%iTeV", ee_tev));
    cinput_zx = cinput_zx + "/CR/" + sample_FullSim[kAllSamples - 1] + ".root";
    tzx[ee]->Add(cinput_zx);
  }

  int nbinsx = KD_bins[0];
  double xlow = -1, xhigh = 1;
  double binwidthx = (xhigh-xlow)/nbinsx;
  int nbinsy = KD_bins[1];
  double ylow=0, yhigh=1;
  TH2F* hqqzz_full = new TH2F("h_template_qqZZ", Form("%i TeV %s q#bar{q}#rightarrow4l Background", erg_tev, user_folder[folder]), nbinsx, xlow, xhigh, nbinsy, ylow, yhigh);
  hqqzz_full->GetXaxis()->SetTitle(Form("tanh(T_{xy} / %.0f mm)", KD_scale));
  hqqzz_full->GetYaxis()->SetTitle("#it{D}_{bkg}");
  TH2F* hggzz_full = (TH2F*)hqqzz_full->Clone();
  TH2F* hzx_full = (TH2F*)hqqzz_full->Clone();
  TH2F* hmc_full = (TH2F*)hqqzz_full->Clone();
  hggzz_full->SetNameTitle("h_template_ggZZ", Form("%i TeV %s gg Background", erg_tev, user_folder[folder]));
  hzx_full->SetNameTitle("h_template_ZX", Form("%i TeV %s Z+X Background", erg_tev, user_folder[folder]));
  hmc_full->SetNameTitle("h_template_totalbkg", Form("%i TeV %s Total Background", erg_tev, user_folder[folder]));
  TH2F* hfull[ntrees+1] ={ hqqzz_full, hggzz_full, hzx_full, hmc_full };
  char* ctreefull[ntrees] ={ "qqZZ", "ggZZ", "ZX" };
  for (int tt = 0; tt < ntrees; tt++) hfull[tt]->Sumw2();
  TH1F* hproj[ntrees+1];

  TChain** tc[ntrees] ={ tqqzz, tggzz, tzx };
  TTree* mytree[ntrees];
  for (int tt = 0; tt < ntrees; tt++){
    const int sizeTree = bkgSize[tt];
    cout << "Tree " << tt << " size: " << sizeTree << endl;
    for (int ss=0; ss<sizeTree; ss++){
      if (tt != 2){
        //        tc[tt][ss]->SetBranchAddress("MC_weight_noxsec", &MC_weight_noxsec);
        tc[tt][ss]->SetBranchAddress("MC_weight", &MC_weight_noxsec);
        tc[tt][ss]->SetBranchAddress("MC_weight_Kfactor", &MC_weight_Kfactor);
        tc[tt][ss]->SetBranchAddress("MC_weight_QQBGGProper", MC_weight_QQBGGProper);
        if (tt == 0){
          tc[tt][ss]->SetBranchAddress("MC_weight_QQZZEWK", &MC_weight_QQZZEWK);
        }

        tc[tt][ss]->SetBranchAddress("D_bkg", &D_bkg);
      }
      else{
        tc[tt][ss]->SetBranchAddress("RunNumber", &RunNumber);
        tc[tt][ss]->SetBranchAddress("LumiNumber", &LumiNumber);

        tc[tt][ss]->SetBranchAddress("ZXfake_weight_OS", ZXfake_weight_OS);
        tc[tt][ss]->SetBranchAddress("ZXfake_weight_SS", &ZXfake_weight_SS);

        tc[tt][ss]->SetBranchAddress("CRflag", &CRflag);

        tc[tt][ss]->SetBranchAddress("Lep1combRelIsoPF", &Lep1combRelIsoPF);
        tc[tt][ss]->SetBranchAddress("Lep2combRelIsoPF", &Lep2combRelIsoPF);
        tc[tt][ss]->SetBranchAddress("Lep3combRelIsoPF", &Lep3combRelIsoPF);
        tc[tt][ss]->SetBranchAddress("Lep4combRelIsoPF", &Lep4combRelIsoPF);

        tc[tt][ss]->SetBranchAddress("Lep3isID", &Lep3isID);
        tc[tt][ss]->SetBranchAddress("Lep4isID", &Lep4isID);

        tc[tt][ss]->SetBranchAddress("Z1ids", &Z1ids);
        tc[tt][ss]->SetBranchAddress("Lep3ID", &Lep3ID);
        tc[tt][ss]->SetBranchAddress("Lep4ID", &Lep4ID);

        tc[tt][ss]->SetBranchAddress("p0plus_VAJHU", &p0plus_VAJHU);
        tc[tt][ss]->SetBranchAddress("p0plus_VAMCFM", &p0plus_VAMCFM);
        tc[tt][ss]->SetBranchAddress("bkg_VAMCFM", &bkg_VAMCFM);
        tc[tt][ss]->SetBranchAddress("p0plus_m4l", &p0plus_m4l);
        tc[tt][ss]->SetBranchAddress("bkg_m4l", &bkg_m4l);
        tc[tt][ss]->SetBranchAddress("p0plus_m4l_ScaleUp", &p0plus_m4l_ScaleUp);
        tc[tt][ss]->SetBranchAddress("bkg_m4l_ScaleUp", &bkg_m4l_ScaleUp);
        tc[tt][ss]->SetBranchAddress("p0plus_m4l_ScaleDown", &p0plus_m4l_ScaleDown);
        tc[tt][ss]->SetBranchAddress("bkg_m4l_ScaleDown", &bkg_m4l_ScaleDown);
        tc[tt][ss]->SetBranchAddress("p0plus_m4l_ResUp", &p0plus_m4l_ResUp);
        tc[tt][ss]->SetBranchAddress("bkg_m4l_ResUp", &bkg_m4l_ResUp);
        tc[tt][ss]->SetBranchAddress("p0plus_m4l_ResDown", &p0plus_m4l_ResDown);
        tc[tt][ss]->SetBranchAddress("bkg_m4l_ResDown", &bkg_m4l_ResDown);
      }

      tc[tt][ss]->SetBranchAddress("ZZMass", &ZZMass);
      tc[tt][ss]->SetBranchAddress("ZZPt", &ZZPt);
      tc[tt][ss]->SetBranchAddress("ZZEta", &ZZEta);
      tc[tt][ss]->SetBranchAddress("ZZPhi", &ZZPhi);

      tc[tt][ss]->SetBranchAddress("Lep1SIP", &Lep1SIP);
      tc[tt][ss]->SetBranchAddress("Lep2SIP", &Lep2SIP);
      tc[tt][ss]->SetBranchAddress("Lep3SIP", &Lep3SIP);
      tc[tt][ss]->SetBranchAddress("Lep4SIP", &Lep4SIP);
      tc[tt][ss]->SetBranchAddress("Lep1_Z1SIP", &Lep1_Z1SIP);
      tc[tt][ss]->SetBranchAddress("Lep2_Z1SIP", &Lep2_Z1SIP);
      tc[tt][ss]->SetBranchAddress("Lep3_Z1SIP", &Lep3_Z1SIP);
      tc[tt][ss]->SetBranchAddress("Lep4_Z1SIP", &Lep4_Z1SIP);
      tc[tt][ss]->SetBranchAddress("OfflinePrimaryVtx_x", &OfflinePrimaryVtx_x);
      tc[tt][ss]->SetBranchAddress("OfflinePrimaryVtx_y", &OfflinePrimaryVtx_y);
      tc[tt][ss]->SetBranchAddress("OfflinePrimaryVtx_z", &OfflinePrimaryVtx_z);
      tc[tt][ss]->SetBranchAddress("OfflinePrimaryVtx_cov_xx", &OfflinePrimaryVtx_cov_xx);
      tc[tt][ss]->SetBranchAddress("OfflinePrimaryVtx_cov_xy", &OfflinePrimaryVtx_cov_xy);
      tc[tt][ss]->SetBranchAddress("OfflinePrimaryVtx_cov_xz", &OfflinePrimaryVtx_cov_xz);
      tc[tt][ss]->SetBranchAddress("OfflinePrimaryVtx_cov_yy", &OfflinePrimaryVtx_cov_yy);
      tc[tt][ss]->SetBranchAddress("OfflinePrimaryVtx_cov_yz", &OfflinePrimaryVtx_cov_yz);
      tc[tt][ss]->SetBranchAddress("OfflinePrimaryVtx_cov_zz", &OfflinePrimaryVtx_cov_zz);
      tc[tt][ss]->SetBranchAddress("OfflinePrimaryVtx_ndof", &OfflinePrimaryVtx_ndof);
      tc[tt][ss]->SetBranchAddress("KalmanCandVtx_x", &KalmanCandVtx_x);
      tc[tt][ss]->SetBranchAddress("KalmanCandVtx_y", &KalmanCandVtx_y);
      tc[tt][ss]->SetBranchAddress("KalmanCandVtx_z", &KalmanCandVtx_z);
      tc[tt][ss]->SetBranchAddress("KalmanCandVtx_chi2", &KalmanCandVtx_chi2);
      tc[tt][ss]->SetBranchAddress("KalmanCandVtx_cov_xx", &KalmanCandVtx_cov_xx);
      tc[tt][ss]->SetBranchAddress("KalmanCandVtx_cov_xy", &KalmanCandVtx_cov_xy);
      tc[tt][ss]->SetBranchAddress("KalmanCandVtx_cov_xz", &KalmanCandVtx_cov_xz);
      tc[tt][ss]->SetBranchAddress("KalmanCandVtx_cov_yy", &KalmanCandVtx_cov_yy);
      tc[tt][ss]->SetBranchAddress("KalmanCandVtx_cov_yz", &KalmanCandVtx_cov_yz);
      tc[tt][ss]->SetBranchAddress("KalmanCandVtx_cov_zz", &KalmanCandVtx_cov_zz);
      tc[tt][ss]->SetBranchAddress("Z1CandVtx_x", &Z1CandVtx_x);
      tc[tt][ss]->SetBranchAddress("Z1CandVtx_y", &Z1CandVtx_y);
      tc[tt][ss]->SetBranchAddress("Z1CandVtx_z", &Z1CandVtx_z);
    }
    mytree[tt] = new TTree(Form("T_2D_%s", ctreefull[tt]), "");
    mytree[tt]->Branch("_weight_", &templateWeight);
    mytree[tt]->Branch("CMS_zz4l_smd", &smd);
    mytree[tt]->Branch("CMS_zz4l_KD", &KD);
    mytree[tt]->Branch("CMS_zz4l_mass", &mzz);
  }
  for (int tt = 0; tt < ntrees; tt++){
    const int sizeTree = bkgSize[tt];
    int totalTEntries[sizeTree];
    int nTraversedEntries = 0;
    for (int ss=0; ss<sizeTree; ss++){
      totalTEntries[ss] = tc[tt][ss]->GetEntries();
      cout << "NEntries for histogram: " << totalTEntries[ss] << endl;
      for (int ev = 0; ev < totalTEntries[ss]; ev++){
        RunNumber = 1;
        LumiNumber = 1;

        MC_weight_noxsec = 1;
        MC_weight_QQZZEWK = 1;
        MC_weight_Kfactor = 1;
        MC_weight_QQBGGProper[0]=1; MC_weight_QQBGGProper[1]=1;

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

        tc[tt][ss]->GetEntry(ev);
        nTraversedEntries++;

        bool passSelection = testSelection(
          Lep1SIP, Lep2SIP, Lep3SIP, Lep4SIP,
          Lep1_Z1SIP, Lep2_Z1SIP, Lep3_Z1SIP, Lep4_Z1SIP,
          KalmanCandVtx_chi2,
          scheme
          );
        if (!passSelection) continue;
        if (tt==2){
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
          if (
            (CRflag==6 || CRflag==10) ||
            (CRflag==8 || CRflag==12)
            ){
            MC_weight_noxsec = ZXfake_weight_OS[scheme] * (targetCRcount[0][scheme] / (targetCRcount[0][scheme] + targetCRcount[1][scheme]));
          }
          else if (
            (CRflag==7 || CRflag==11) ||
            (CRflag==5 || CRflag==9)
            ){
            int OSoverSS_biny = -1;
            if ((CRflag==5 || CRflag==9) && Z1ids==-169) OSoverSS_biny = 1; // 4mu
            if ((CRflag==7 || CRflag==11) && Z1ids==-169) OSoverSS_biny = 3; // 2mu2e
            if ((CRflag==5 || CRflag==9) && Z1ids==-121) OSoverSS_biny = 4; // 2e2mu
            if ((CRflag==7 || CRflag==11) && Z1ids==-121) OSoverSS_biny = 6; // 4e
            float SSwgt_OSoverSS = 1;
            if (OSoverSS_biny>0) SSwgt_OSoverSS = hRatio->GetBinContent(hRatio->GetXaxis()->FindBin(ZZMass), OSoverSS_biny);

            MC_weight_noxsec = ZXfake_weight_SS*SSwgt_OSoverSS;
            MC_weight_noxsec *= (targetCRcount[1][scheme] / (targetCRcount[0][scheme] + targetCRcount[1][scheme]));
          }
        }

        double wgt = MC_weight_noxsec;
        if (tt==0) wgt *= MC_weight_QQZZEWK;
        else if (tt==1) wgt *= MC_weight_Kfactor;

        if (!(ZZMass >= ZZMass_PeakCut[0] && ZZMass < ZZMass_PeakCut[1])) continue;

        int index_BS = 1;
        bool matchFound = false;
        if (tt==2) index_BS = 0;
        for (int evBS=0; evBS < tBeam[index_BS]->GetEntries(); evBS++){
          tBeam[index_BS]->GetEntry(evBS);
          if (RunNumber == RunNumber_Ref && LumiNumber == LumiNumber_Ref){
            matchFound = true;
            break;
          }
          else continue;
        }
        if (!matchFound) cerr << "No beamspot matching possible for tree " << tt << " / " << ss << "!" << endl;

        CandBSVtx_x = KalmanCandVtx_x - BeamPosX;
        CandBSVtx_y = KalmanCandVtx_y - BeamPosY;
        CandBSVtx_z = KalmanCandVtx_z - BeamPosZ;

        TrueCandBSVtx_x = GenIntVtx_x - BeamPosX;
        TrueCandBSVtx_y = GenIntVtx_y - BeamPosY;
        TrueCandBSVtx_z = GenIntVtx_z - BeamPosZ;

        Dxy_BS_true = (TrueCandBSVtx_x*cos(GenHPhi) + TrueCandBSVtx_y*sin(GenHPhi)) * cm_to_microns;
        Txy_BS_true = Dxy_BS_true*GenHMass / GenHPt;

        Dxy_BS = (CandBSVtx_x*cos(ZZPhi) + CandBSVtx_y*sin(ZZPhi)) * cm_to_microns;
        Txy_BS = Dxy_BS*ZZMass / ZZPt;

        CandVtx_x = KalmanCandVtx_x - OfflinePrimaryVtx_x;
        CandVtx_y = KalmanCandVtx_y - OfflinePrimaryVtx_y;
        CandVtx_z = KalmanCandVtx_z - OfflinePrimaryVtx_z;

        TrueCandVtx_x = GenIntVtx_x - GenPrimaryVtx_x;
        TrueCandVtx_y = GenIntVtx_y - GenPrimaryVtx_y;
        TrueCandVtx_z = GenIntVtx_z - GenPrimaryVtx_z;

        Dxy_true = (TrueCandVtx_x*cos(GenHPhi) + TrueCandVtx_y*sin(GenHPhi)) * cm_to_microns;
        Txy_true = Dxy_true*GenHMass / GenHPt;

        KD = Txy_BS;
        KD = tanh(KD / KD_scale); // Scaled to due to tanh function
        if (KD >= 1.) {
          KD = 1. - 1e-4;
        }
        if (KD <= -1.) {
          KD = -1. + 1e-4;
        }
        if (tt==2) D_bkg = p0plus_VAJHU*p0plus_m4l / (p0plus_VAJHU*p0plus_m4l + bkg_VAMCFM*bkg_m4l);
        if (D_bkg!=D_bkg) continue;
        smd = D_bkg;
        mzz = ZZMass;
        templateWeight = wgt;
        if (templateWeight<0) continue;
        hfull[tt]->Fill(KD, smd, templateWeight);
      }
    }
    cout << "nTraversed: " << nTraversedEntries << endl;
    nTraversedEntries=0;

    scale_yield[tt] = targetYield[tt]/hfull[tt]->Integral(0, hfull[tt]->GetNbinsX() + 1, 0, hfull[tt]->GetNbinsY() + 1);
    cout << "Scaling " << hfull[tt]->GetName() << " by " << scale_yield[tt] << " to " << targetYield[tt] << endl;
    hfull[tt]->Scale(scale_yield[tt]);

    double sum_filled=0;

    //    cout << "Starting to shuffle the entries... ";
    int treeEntries = 0;
    int* totalAccEv = new int[sizeTree];
    int* totalAccTEntries = new int[sizeTree];
    cout << "Tree set " << tt << " nEvents | nAccumulated: ";
    for (int ss=0; ss<sizeTree; ss++){
      cout << totalTEntries[ss] << " | ";
      treeEntries += totalTEntries[ss];
      totalAccEv[ss] = 0;
      totalAccTEntries[ss] = totalTEntries[ss];
      if (ss>0) totalAccTEntries[ss] += totalAccTEntries[ss-1];
      cout << totalAccTEntries[ss] << "\t";
    }
    cout << "\t" << treeEntries << endl;
    TRandom3 rand(ntrees*(3*EnergyIndex+folder)+tt);

    while (treeEntries>0){
      RunNumber = 1;
      LumiNumber = 1;

      MC_weight_noxsec = 1;
      MC_weight_QQZZEWK = 1;
      MC_weight_Kfactor = 1;
      MC_weight_QQBGGProper[0]=1; MC_weight_QQBGGProper[1]=1;

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

      int ev = rand.Integer(treeEntries);
      int treeCode = -1;
      for (int ss=0; ss<sizeTree; ss++){
        if (ev<totalAccTEntries[ss] && (totalTEntries[ss]-totalAccEv[ss])!=0){ treeCode=ss; break; }
      }
      if (treeCode<0) continue;
      int entryNumber = totalAccEv[treeCode];
      //      cout << "ev: " << ev << " treeCode: " << treeCode << " entryNumber: " << entryNumber << endl;
      if (tc[tt][treeCode]->GetEntries()<=entryNumber) cout << "ERROR!" << endl;
      tc[tt][treeCode]->GetEntry(entryNumber);

      totalAccEv[treeCode]++;
      for (int ss=treeCode; ss<sizeTree; ss++) totalAccTEntries[ss]--;
      treeEntries--;
      nTraversedEntries++;

      bool passSelection = testSelection(
        Lep1SIP, Lep2SIP, Lep3SIP, Lep4SIP,
        Lep1_Z1SIP, Lep2_Z1SIP, Lep3_Z1SIP, Lep4_Z1SIP,
        KalmanCandVtx_chi2,
        scheme
        );
      if (!passSelection) continue;
      if (tt==2){
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
        if (
          (CRflag==6 || CRflag==10) ||
          (CRflag==8 || CRflag==12)
          ){
          MC_weight_noxsec = ZXfake_weight_OS[scheme] * (targetCRcount[0][scheme] / (targetCRcount[0][scheme] + targetCRcount[1][scheme]));
        }
        else if (
          (CRflag==7 || CRflag==11) ||
          (CRflag==5 || CRflag==9)
          ){
          int OSoverSS_biny = -1;
          if ((CRflag==5 || CRflag==9) && Z1ids==-169) OSoverSS_biny = 1; // 4mu
          if ((CRflag==7 || CRflag==11) && Z1ids==-169) OSoverSS_biny = 3; // 2mu2e
          if ((CRflag==5 || CRflag==9) && Z1ids==-121) OSoverSS_biny = 4; // 2e2mu
          if ((CRflag==7 || CRflag==11) && Z1ids==-121) OSoverSS_biny = 6; // 4e
          float SSwgt_OSoverSS = 1;
          if (OSoverSS_biny>0) SSwgt_OSoverSS = hRatio->GetBinContent(hRatio->GetXaxis()->FindBin(ZZMass), OSoverSS_biny);

          MC_weight_noxsec = ZXfake_weight_SS*SSwgt_OSoverSS;
          MC_weight_noxsec *= (targetCRcount[1][scheme] / (targetCRcount[0][scheme] + targetCRcount[1][scheme]));
        }
      }

      if (!(ZZMass >= ZZMass_PeakCut[0] && ZZMass < ZZMass_PeakCut[1])) continue;

      int index_BS = 1;
      bool matchFound = false;
      if (tt==2) index_BS = 0;
      for (int evBS=0; evBS < tBeam[index_BS]->GetEntries(); evBS++){
        tBeam[index_BS]->GetEntry(evBS);
        if (RunNumber == RunNumber_Ref && LumiNumber == LumiNumber_Ref){
          matchFound = true;
          break;
        }
        else continue;
      }
      if (!matchFound) cerr << "No beamspot matching possible for tree " << tt << " / " << treeCode << "!" << endl;

      CandBSVtx_x = KalmanCandVtx_x - BeamPosX;
      CandBSVtx_y = KalmanCandVtx_y - BeamPosY;
      CandBSVtx_z = KalmanCandVtx_z - BeamPosZ;

      TrueCandBSVtx_x = GenIntVtx_x - BeamPosX;
      TrueCandBSVtx_y = GenIntVtx_y - BeamPosY;
      TrueCandBSVtx_z = GenIntVtx_z - BeamPosZ;

      Dxy_BS_true = (TrueCandBSVtx_x*cos(GenHPhi) + TrueCandBSVtx_y*sin(GenHPhi)) * cm_to_microns;
      Txy_BS_true = Dxy_BS_true*GenHMass / GenHPt;

      Dxy_BS = (CandBSVtx_x*cos(ZZPhi) + CandBSVtx_y*sin(ZZPhi)) * cm_to_microns;
      Txy_BS = Dxy_BS*ZZMass / ZZPt;

      CandVtx_x = KalmanCandVtx_x - OfflinePrimaryVtx_x;
      CandVtx_y = KalmanCandVtx_y - OfflinePrimaryVtx_y;
      CandVtx_z = KalmanCandVtx_z - OfflinePrimaryVtx_z;

      TrueCandVtx_x = GenIntVtx_x - GenPrimaryVtx_x;
      TrueCandVtx_y = GenIntVtx_y - GenPrimaryVtx_y;
      TrueCandVtx_z = GenIntVtx_z - GenPrimaryVtx_z;

      Dxy_true = (TrueCandVtx_x*cos(GenHPhi) + TrueCandVtx_y*sin(GenHPhi)) * cm_to_microns;
      Txy_true = Dxy_true*GenHMass / GenHPt;

      KD = Txy_BS;
      KD = tanh(KD / KD_scale); // Scaled to due to tanh function
      if (KD >= 1.) {
        KD = 1. - 1e-4;
      }
      if (KD <= -1.) {
        KD = -1. + 1e-4;
      }

      if (tt==2) D_bkg = p0plus_VAJHU*p0plus_m4l / (p0plus_VAJHU*p0plus_m4l + bkg_VAMCFM*bkg_m4l);
      if (D_bkg!=D_bkg) continue;

      double wgt = MC_weight_noxsec;
      if (tt==0) wgt *= MC_weight_QQZZEWK;
      else if (tt==1) wgt *= MC_weight_Kfactor;

      templateWeight = wgt*scale_yield[tt];

      smd = D_bkg;
      mzz = ZZMass;
      if (templateWeight<0) continue;

      sum_filled += templateWeight;
      mytree[tt]->Fill();
    }

    delete[] totalAccEv;
    delete[] totalAccTEntries;

    cout << mytree[tt]->GetName() << " events: " << sum_filled << endl;
    cout << "nTraversed final: " << nTraversedEntries << endl;
  }

  for (int tt = 0; tt < ntrees; tt++) hfull[ntrees]->Add(hfull[tt]);
  for (int tt = 0; tt < ntrees+1; tt++){
    TH1F* hProjKD = (TH1F*)hfull[tt]->ProjectionX();
    TH2F* hfull_conditional = (TH2F*)hfull[tt]->Clone(Form("%s_conditional", hfull[tt]->GetName()));
    hProjKD->SetTitle(Form("%s Projection", hfull[tt]->GetTitle()));
    hfull_conditional->SetTitle(Form("%s Conditional", hfull[tt]->GetTitle()));
    hproj[tt] = hProjKD;

    for (int biny = 0; biny <= hfull_conditional->GetNbinsY()+1; biny++){
      double integralX = hfull_conditional->Integral(0, hfull_conditional->GetNbinsX() + 1, biny, biny);
      if (integralX != 0){
        for (int binx = 0; binx <= hfull_conditional->GetNbinsX()+1; binx++) hfull_conditional->SetBinContent(binx, biny, (hfull_conditional->GetBinContent(binx, biny) / integralX));
      }
    }

    foutput->WriteTObject(hfull[tt]);
    foutput->WriteTObject(hfull_conditional);
    foutput->WriteTObject(hproj[tt]);
    if (tt<ntrees) foutput->WriteTObject(mytree[tt]);
    delete hfull_conditional;
    delete hfull[tt];
    if (tt < ntrees){
      delete mytree[tt];
      for (int ss=0; ss<bkgSize[tt]; ss++) delete tc[tt][ss];
    }
  }
  for (int bs=0; bs<2; bs++) delete tBeam[bs];

  foutput->Close();
  delete hRatio;
  fSSinput->Close();
}


void assembleCombineEmbeddedToys(double target_ctau, int BSMSample, int iProd, double BSMNormScheme, int scheme){
  if (iProd==0 && !(BSMSample<=kfLambda1_1 || BSMSample==kfLambda1_05)) return;
  else if (iProd==1 && BSMSample>=kNumFiles_VBFHVVBSM) return;
  else if (iProd==2 && BSMSample>=kNumFiles_WHVVBSM) return;
  else if (iProd==3 && (BSMSample>=kNumFiles_ZHVVBSM || BSMSample==kfg2_1_fg4_0)) return;
  else if (iProd==4 && (BSMSample>0)) return;
  if (BSMNormScheme!=0 && BSMSample==0) return;

  RooFit::Verbose(false);
  RooFit::PrintLevel(-1000);
  RooMsgService::instance().setStreamStatus(1, false);
  RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);

  const int kNumFiles=kNumFiles_GGHVVBSM;
  const int kNumHypo=kNumGGHVVBSM;
  const int gapZZHypo=kfLambda1_03;
  const int gapZZHypo_2=kfLambda1_05_pL190;
  const int firstNonZZHypo=kfZG_1_fGG_0;

  /***** RUN CONDITIONS *****/
  if (BSMNormScheme<0) assert(0);
  if (BSMSample>=kNumHypo) assert(0);
  /*** END RUN CONDITIONS ***/
  const int ntrees = 4;
  TString strProcess[ntrees] ={ "CTau", "qqZZ", "ggZZ", "ZX" };

  TString KDlist_Trees[N_KDs] ={
    "CMS_zz4l_smd",
    "CMS_zz4l_KD"
  };
  TString Kinematicslist_Trees[N_Kinematics] ={
    "CMS_zz4l_mass"
  };

  TString KDlist_Write[N_KDs] ={ // CHANGE THESE TO MATCH DATACARDS
    "CMS_zz4l_smd",
    "CMS_zz4l_KD"
  };
  TString Kinematicslist_Write[N_Kinematics] ={ // CHANGE THESE TO MATCH DATACARDS
    "CMS_zz4l_mass"
  };
  TString KDlist_WriteD[N_KDs] ={ // CHANGE THESE TO MATCH DATACARDS
    "CMS_zz4l_smd/D",
    "CMS_zz4l_KD/D"
  };
  TString Kinematicslist_WriteD[N_Kinematics] ={ // CHANGE THESE TO MATCH DATACARDS
    "CMS_zz4l_mass/D"
  };

  const int nCat_Tree=6; // Input categories: Should match input file
  const int nCat_Write=nCat_Tree; // Written categories: Whatever names and count you wish
  TString chan[nCat_Tree]={ "ch1", "ch2", "ch3", "ch4", "ch5", "ch6" }; // Input categories: Should match input file
  TString chan_Write[nCat_Write]={ "ch1", "ch2", "ch3", "ch4", "ch5", "ch6" }; // Written categories: Whatever names and count you wish


  double weight=0, weightToy=0;
  double KD_pass[N_KDs];
  double Kinematics_pass[N_Kinematics];

  RooDataSet *data[nCat_Write];
  RooDataSet *toydataset[nCat_Write];
  RooDataSet *data_asimov[nCat_Write];
  RooDataSet *toydataset_asimov[nCat_Write];

  TTree* tin[ntrees][nCat_Tree];
  TFile* finput[ntrees][nCat_Tree];
  double nEvents[ntrees][nCat_Tree]={ { 0 } };
  int nEntries[ntrees][nCat_Tree]={ { 0 } };
  vector<int> nEvents_Toy[ntrees][nCat_Tree];
  TTree* outTree[nCat_Write];

  TString OUTPUT_NAME = "AsimovToys_TxyBS_Hypo";
  OUTPUT_NAME.Append(Form("%i", BSMSample));
  OUTPUT_NAME.Append("_");
  OUTPUT_NAME.Append(Form("CTau%.0f", target_ctau));
  if (iProd>0) OUTPUT_NAME.Append(Form("_%s", strHiggsProduction[iProd-1].Data()));
  OUTPUT_NAME =  OUTPUT_NAME + "_" + cutNames[scheme];
  if (BSMNormScheme!=0) OUTPUT_NAME.Append("_SMLikeNorm");
  OUTPUT_NAME = OUTPUT_NAME + ".root";

  TString cinput_common = user_dir_hep + "Analysis/AsimovToys/";
  TString coutput_common = user_dir_hep + "Analysis/AsimovToys/";
  TString coutput = coutput_common + OUTPUT_NAME;
  TFile* foutput = new TFile(coutput, "recreate");

  /*** THE BOSS ***/
  RooCategory cat("CMS_channel", "CMS_channel");
  for (int idx=0; idx<nCat_Write; idx++){
    cat.defineType(chan_Write[idx], idx);
    cat.setLabel(chan_Write[idx]);
  };


  for (int idx=0; idx<nCat_Write; idx++){
    outTree[idx] = new TTree(chan_Write[idx], chan_Write[idx]);
    outTree[idx]->Branch("_weight_", &weight, "_weight_/D");
    for (int k=0; k<N_KDs; k++) outTree[idx]->Branch(KDlist_Write[k], (KD_pass+k), KDlist_WriteD[k]);
    for (int k=0; k<N_Kinematics; k++) outTree[idx]->Branch(Kinematicslist_Write[k], (Kinematics_pass+k), Kinematicslist_WriteD[k]);
  }

  for (int erg_tev=7; erg_tev<=8; erg_tev++){
    TString comstring;
    comstring.Form("%iTeV", erg_tev);
    TString erg_dir;
    erg_dir.Form("LHC_%iTeV", erg_tev);
    int EnergyIndex=0; if (erg_tev==8) EnergyIndex=1;

    for (int folder=0; folder<3; folder++){
      int idx = 3*EnergyIndex + folder;

      cout << "Channel " << chan_Write[idx] << " corresponds to " << comstring << " " << user_folder[folder] << endl;

      TString INPUT_NAME = "AsimovTree_background_TxyBS_";
      INPUT_NAME.Append(user_folder[folder]);
      INPUT_NAME.Append("_");
      INPUT_NAME.Append(comstring);
      INPUT_NAME =  INPUT_NAME + "_" + cutNames[scheme] + ".root";
      TString cinput = INPUT_NAME;
      cinput.Prepend(cinput_common);

      for (int tt=1; tt<ntrees; tt++){
        TString treename = strProcess[tt];
        treename.Prepend("T_2D_");
        finput[tt][idx] = new TFile(cinput, "read");
        tin[tt][idx] = (TTree*)finput[tt][idx]->Get(treename);
      }

      INPUT_NAME = "AsimovTree_Sig_TxyBS_";
      if (iProd>0) INPUT_NAME.Append(Form("%s_", strHiggsProduction[iProd-1].Data()));
      if (BSMSample>0) INPUT_NAME.Append(Form("Hypo%i_", BSMSample));
      INPUT_NAME.Append(user_folder[folder]);
      INPUT_NAME.Append("_");
      INPUT_NAME.Append(comstring);
      INPUT_NAME.Append("_");
      INPUT_NAME = INPUT_NAME + cutNames[scheme];
      if (BSMNormScheme!=0) INPUT_NAME.Append("_SMLikeNorm");
      INPUT_NAME = INPUT_NAME + ".root";


      cinput = INPUT_NAME;
      cinput.Prepend(cinput_common);
      for (int tt=0; tt<1; tt++){
        TString treename = strProcess[tt];
        treename.Append(Form("%.0f", target_ctau));
        treename.Prepend("T_2D_");
        finput[tt][idx] = new TFile(cinput, "read");
        tin[tt][idx] = (TTree*)finput[tt][idx]->Get(treename);
      }
      foutput->cd();

      cout << "Number of signal events received: " << tin[0][idx]->GetEntries() << endl;

      for (int tt=0; tt<ntrees; tt++){
        tin[tt][idx]->SetBranchAddress("_weight_", &weight);
        for (int k=0; k<N_KDs; k++) tin[tt][idx]->SetBranchAddress(KDlist_Write[k], (KD_pass+k));
        for (int k=0; k<N_Kinematics; k++) tin[tt][idx]->SetBranchAddress(Kinematicslist_Write[k], (Kinematics_pass+k));
        nEntries[tt][idx] = tin[tt][idx]->GetEntries();
        for (int ev=0; ev<nEntries[tt][idx]; ev++){
          tin[tt][idx]->GetEntry(ev); nEvents[tt][idx] += weight; outTree[idx]->Fill();
        }

        cout << "Number of events from process " << tt << " channel " << idx << ": " << nEvents[tt][idx] << endl;

        TRandom3 rand(nCat_Tree*tt+idx);
        int accEntries = 0;
        while (accEntries<nEntries[tt][idx]){
          int nToyInst = (int)rand.PoissonD(nEvents[tt][idx])+0.5;
          accEntries += nToyInst;
          if (accEntries<nEntries[tt][idx]){
            nEvents_Toy[tt][idx].push_back(nToyInst);
          }
        }
      }

    }
  }

  int nAllowedToys = -1;
  for (int idx=0; idx<nCat_Write; idx++){
    for (int tt=0; tt<ntrees; tt++){
      int nToysInst = nEvents_Toy[tt][idx].size();
      if (nAllowedToys<0 || nToysInst<nAllowedToys) nAllowedToys = nToysInst;
    }
  }
  cout << "Number of allowed toys is " << nAllowedToys << endl;
  TTree** tToys[nCat_Write];
  for (int idx=0; idx<nCat_Write; idx++){
    tToys[idx] = new TTree*[nAllowedToys];
    for (int iToy=0; iToy<nAllowedToys; iToy++){
      TString toyTreeName = Form("%s_Toy%i", chan_Write[idx].Data(), iToy);
      tToys[idx][iToy] = new TTree(toyTreeName, toyTreeName);
      tToys[idx][iToy]->Branch("_weight_", &weightToy, "_weight_/D");
      for (int k=0; k<N_KDs; k++) tToys[idx][iToy]->Branch(KDlist_Write[k], (KD_pass+k), KDlist_WriteD[k]);
      for (int k=0; k<N_Kinematics; k++) tToys[idx][iToy]->Branch(Kinematicslist_Write[k], (Kinematics_pass+k), Kinematicslist_WriteD[k]);
    }

    for (int tt=0; tt<ntrees; tt++){
      int accEvents=0;
      int iToy=0;
      double* sumWeights = new double[nAllowedToys];
      for (int it=0; it<nAllowedToys; it++) sumWeights[it] = 0;
      for (int ev=0; ev<nEntries[tt][idx]; ev++){
        tin[tt][idx]->GetEntry(ev);
        if (accEvents>=nEvents_Toy[tt][idx].at(iToy)){
          iToy++;
          accEvents=0;
        }
        if (iToy>=nAllowedToys) break;

        accEvents++;
        sumWeights[iToy] += weight;
      }

      iToy=0;
      accEvents=0;
      for (int ev=0; ev<nEntries[tt][idx]; ev++){
        tin[tt][idx]->GetEntry(ev);
        if (accEvents>=nEvents_Toy[tt][idx].at(iToy)){
          iToy++;
          accEvents=0;
        }
        if (iToy>=nAllowedToys) break;

        accEvents++;
        weightToy = weight;
        double event_scale = ((double)nEvents_Toy[tt][idx].at(iToy)) / sumWeights[iToy];
        weightToy *= event_scale;
        tToys[idx][iToy]->Fill();
      }

      delete[] sumWeights;
      delete tin[tt][idx];
      finput[tt][idx]->Close();
      foutput->cd();
    }
  }

  RooRealVar* rweightFit = new RooRealVar("_weight_", "_weight_", 0., 1.e20);
  RooRealVar* rCMS_zz4l_VarContainer[N_KDs+N_Kinematics] ={ 0 };
  RooArgSet rContainer(*rweightFit);
  for (int k=0; k<N_KDs; k++){
    rCMS_zz4l_VarContainer[k] = new RooRealVar(KDlist_Write[k], KDlist_Write[k], KD_ranges[k][0], KD_ranges[k][1]);
    rCMS_zz4l_VarContainer[k]->setBins(KD_bins[k]);
    rContainer.add(*rCMS_zz4l_VarContainer[k]);
  }
  for (int k=0; k<N_Kinematics; k++){
    rCMS_zz4l_VarContainer[k+N_KDs] = new RooRealVar(Kinematicslist_Write[k], Kinematicslist_Write[k], Kinematics_ranges[k][0], Kinematics_ranges[k][1]);
    rCMS_zz4l_VarContainer[k+N_KDs]->setBins(Kinematics_bins[k]);
    rContainer.add(*rCMS_zz4l_VarContainer[k+N_KDs]);
  }

  RooDataSet* total_asimov;
  for (int idx=0; idx<nCat_Write; idx++){
    cout << outTree[idx]->GetEntries() << "...";
    data_asimov[idx] = new RooDataSet(Form("data%d", idx), Form("data%d", idx), outTree[idx], rContainer, "", "_weight_");
    data_asimov[idx]->weightError(RooAbsData::None);

    RooArgSet toyArgs(rContainer);
    toyArgs.add(cat);
    toydataset_asimov[idx] = new RooDataSet(Form("toy_asimov%i", idx), Form("toy_asimov%i", idx), toyArgs, Index(cat), WeightVar("_weight_"), Import(chan_Write[idx], *data_asimov[idx]));

    //    toydataset_asimov[idx]->Print("v");
    if (idx==0) total_asimov = toydataset_asimov[idx];
    else total_asimov->append(*(toydataset_asimov[idx]));
  }
  total_asimov->SetName("toy_asimov");
  RooCategory* cata = dynamic_cast<RooCategory *>(total_asimov->get()->find("CMS_channel"));
  int na = cata->numBins((const char *)0);
  foutput->cd();
  total_asimov->Write("toy_asimov");
  total_asimov->Print("v");
  if (total_asimov->isWeighted()) cout << "Asimov dataset is weighted." << endl;
  //  double correlationVal = total_asimov->correlation(*rCMS_zz4l_VarContainer[2], *rCMS_zz4l_VarContainer[1]);
  //  cout << "Correlation from the asimov: " << correlationVal << endl;

  for (int kk=0; kk<N_KDs+N_Kinematics; kk++){
    TString cValPlot = "ValidateAsimov_";
    cValPlot.Append(rCMS_zz4l_VarContainer[kk]->GetName());
    TString cCanValPlot = "canvasValidateAsimov_";
    cCanValPlot.Append(rCMS_zz4l_VarContainer[kk]->GetName());
    TCanvas* can = new TCanvas(cCanValPlot, "", 800, 800);
    RooPlot* valPlot = new RooPlot(cValPlot, "", *rCMS_zz4l_VarContainer[kk], rCMS_zz4l_VarContainer[kk]->getMin(), rCMS_zz4l_VarContainer[kk]->getMax(), (kk<N_KDs ? KD_bins[kk] : Kinematics_bins[kk-N_KDs]));
    total_asimov->plotOn(valPlot);
    can->cd();
    valPlot->Draw();
    foutput->WriteTObject(can);
    delete valPlot;
    can->Close();
  }

  for (int idx=0; idx<nCat_Write; idx++){
    if (data_asimov[idx]!=0){
      delete data_asimov[idx]; data_asimov[idx]=0;
    }
    if (toydataset_asimov[idx]!=0){
      delete toydataset_asimov[idx]; toydataset_asimov[idx]=0;
    }
  }
  //  if (total_asimov!=0) delete total_asimov;
  cout << "\nRemoved asimov dataset" << endl;

  TDirectory *toysdir = foutput->mkdir("toys"); // This is where the toys will reside
  toysdir->cd();
  for (int iToy=0; iToy<nAllowedToys; iToy++){
    RooDataSet* total_toy;
    for (int idx=0; idx<nCat_Write; idx++){
      cout << tToys[idx][iToy]->GetEntries() << "...";
      data[idx] = new RooDataSet(Form("data%d", idx), Form("data%d", idx), tToys[idx][iToy], rContainer, "", "_weight_");
      data[idx]->weightError(RooAbsData::None);

      RooArgSet toyArgs(rContainer);
      toyArgs.add(cat);
      toydataset[idx] = new RooDataSet(Form("toy%i", idx), Form("toy%i", idx), toyArgs, Index(cat), WeightVar("_weight_"), Import(chan_Write[idx], *data[idx]));

      //      toydataset[idx]->Print("v");
      if (idx==0) total_toy = toydataset[idx];
      else total_toy->append(*(toydataset[idx]));
    }
    total_toy->SetName(Form("toy_%d", iToy));
    foutput->cd();
    toysdir->cd();
    total_toy->Write(Form("toy_%d", iToy));
    total_toy->Print("v");
    if (total_toy->isWeighted()) cout << "Toy dataset is weighted." << endl;

    for (int idx=0; idx<nCat_Write; idx++){
      if (data[idx]!=0){
        delete data[idx]; data[idx]=0;
      }
      if (toydataset[idx]!=0){
        delete toydataset[idx]; toydataset[idx]=0;
      }
    }
    //    if (total_toy!=0) delete total_toy;
    cout << "\nRemoved toy dataset " << iToy << endl;
  }


  for (int idx=0; idx<nCat_Write; idx++){
    for (int iToy=0; iToy<nAllowedToys; iToy++) delete tToys[idx][iToy];
    delete[] tToys[idx];
    delete outTree[idx];
  }
  foutput->Close();
}

void compareEmbeddedAsimov_SignalTemplates(double target_ctau, int iProd, int scheme){
  if (!(target_ctau==0 || target_ctau==100 || target_ctau==500 || target_ctau==1000)) return;
  if (target_ctau!=0 && iProd>0) return;
  bool isBkg = false;
  if ((iProd<0 && target_ctau>0) || (iProd>0 && target_ctau<0) || iProd<-3) return;
  else if (iProd<0 || target_ctau<0) isBkg = true;

  TString channame[3] ={ "4mu", "4e", "2e2mu" };
  TString yTitle = "Events / bin";
  TString xTitle = "c#Deltat (#mum)";

  TFile* finput = 0;

  TString OUTPUT_NAME;
  if (!isBkg){
    OUTPUT_NAME = "AsimovToys_SignalComparison_TxyBS_";
    OUTPUT_NAME.Append("_");
    OUTPUT_NAME.Append(Form("CTau%.0f", target_ctau));
    if (iProd>0) OUTPUT_NAME.Append(Form("_%s", strHiggsProduction[iProd-1].Data()));
    else OUTPUT_NAME.Append("_ggH");
    OUTPUT_NAME =  OUTPUT_NAME + "_" + cutNames[scheme];
    OUTPUT_NAME = OUTPUT_NAME + ".root";
  }
  else{
    OUTPUT_NAME = "AsimovToys_BackgroundComparison_TxyBS";
    if (iProd==-1) OUTPUT_NAME.Append("_ggZZ");
    else if (iProd==-2) OUTPUT_NAME.Append("_qqZZ");
    else if (iProd==-3) OUTPUT_NAME.Append("_ZX");
    OUTPUT_NAME =  OUTPUT_NAME + "_" + cutNames[scheme];
    OUTPUT_NAME = OUTPUT_NAME + ".root";

    cout << "Prod is a bkg. Creating " << OUTPUT_NAME << endl;
  }

  TString cinput_common = user_dir_hep + "Analysis/AsimovToys/";
  TString coutput_common = user_dir_hep + "Analysis/AsimovToys/";
  TString coutput = coutput_common + OUTPUT_NAME;
  TFile* foutput = new TFile(coutput, "recreate");

  TH1D* htemplate_Total = 0;
  TH1D* htemplate_tth_Total = 0;
  TH1F* hAsimov_Total[kNumGGHVVBSM] ={ 0 };

  for (int erg_tev=7; erg_tev<=8; erg_tev++){
    TString comstring;
    comstring.Form("%iTeV", erg_tev);
    TString erg_dir;
    erg_dir.Form("LHC_%iTeV", erg_tev);
    int EnergyIndex=1;
    if (erg_tev==7) EnergyIndex=0;

    for (int folder=0; folder<3; folder++){
      TString cinput_Txymain_templates = templateDir + comstring + "/" + channame[folder] + "_templates_";
      TString cinput = cinput_Txymain_templates;
      if (!isBkg) cinput.Append(Form("%.0f%s", target_ctau, "_Modified.root"));
      else cinput.Append("Merged_bkg.root");
      finput = new TFile(cinput, "read");

      if (isBkg && (finput==0 || finput->IsZombie())){
        cout << "Failed to open input bkg file " << cinput << endl;
      }
      if (isBkg) cout << "Attempted template file " << cinput << endl;
      TString hname;
      if (!isBkg){
        hname = "T_2D_TxyNominal_";
        if (iProd==0) hname.Append("ggH");
        else hname.Append(strHiggsProduction[iProd-1]);
      }
      else{
        if (iProd==-1) hname = "ggZZ";
        else if (iProd==-2) hname = "qqZZ";
        else if (iProd==-3) hname = "ZX";
        hname.Prepend("template_");
        hname.Append("_TxyNominal");
      }

      TH1F* hAsimov[kNumGGHVVBSM]={ 0 };
      TH1F* hAsimov_fill[kNumGGHVVBSM]={ 0 };

      TH1D* htemp = (TH1D*)finput->Get(hname);
      if (!htemp) cout << "htemp " << hname << " not found!" << endl;
      foutput->cd();
      TH1D* htemplate = (TH1D*)htemp->Clone(Form("%s_clone", htemp->GetName()));
      delete htemp;
      htemp = 0;
      TH1D* htemplate_tth = 0;
      if (!isBkg){
        htemp = (TH1D*)finput->Get("T_2D_TxyNominal_ttH");
        if (!htemp) cout << "htemp " << hname << " not found!" << endl;
        else{
          htemplate_tth = (TH1D*)htemp->Clone(Form("%s_clone", htemp->GetName()));
          delete htemp;
        }
      }
      finput->Close();

      TH1F* hstyle = (TH1F*)htemplate->Clone("hstyle");

      htemplate = (TH1D*)transform_tanh_varBin((TH1F*)htemplate);
      htemplate->SetName(Form("%s_%s_%s", htemplate->GetName(), comstring.Data(), channame[folder].Data()));
      foutput->WriteTObject(htemplate);
      cout << "Template integral: " << htemplate->Integral() << endl;

      htemplate->SetLineWidth(2);
      htemplate->SetLineColor(kBlack);
      htemplate->SetTitle("");
      if (htemplate_tth!=0){
        htemplate_tth = (TH1D*)transform_tanh_varBin((TH1F*) htemplate_tth);
        htemplate_tth->SetLineWidth(2);
        htemplate_tth->SetLineStyle(3);
        htemplate_tth->SetLineColor(kBlack);
        htemplate_tth->SetTitle("");
      }

      htemplate->GetXaxis()->SetTitle(xTitle);
      htemplate->GetYaxis()->SetTitle(yTitle);
      htemplate->GetXaxis()->SetNdivisions(505);
      htemplate->GetXaxis()->SetLabelFont(42);
      htemplate->GetXaxis()->SetLabelOffset(0.007);
      htemplate->GetXaxis()->SetLabelSize(0.04);
      htemplate->GetXaxis()->SetTitleSize(0.06);
      htemplate->GetXaxis()->SetTitleOffset(0.9);
      htemplate->GetXaxis()->SetTitleFont(42);
      htemplate->GetYaxis()->SetNdivisions(505);
      htemplate->GetYaxis()->SetLabelFont(42);
      htemplate->GetYaxis()->SetLabelOffset(0.007);
      htemplate->GetYaxis()->SetLabelSize(0.04);
      htemplate->GetYaxis()->SetTitleSize(0.06);
      htemplate->GetYaxis()->SetTitleOffset(1.1);
      htemplate->GetYaxis()->SetTitleFont(42);

      if (htemplate_Total!=0) htemplate_Total->Add(htemplate);
      else htemplate_Total = (TH1D*)htemplate->Clone("htemplate_total");
      htemplate->Scale(1./htemplate->Integral());
      if (htemplate_tth!=0){
        if (htemplate_tth_Total!=0) htemplate_tth_Total->Add(htemplate_tth);
        else htemplate_tth_Total = (TH1D*)htemplate_tth->Clone("htemplate_tth_total");
        htemplate_tth->Scale(1./htemplate_tth->Integral());
      }

      hstyle->Scale(1./hstyle->Integral());
      hstyle->SetLineWidth(2);
      hstyle->SetLineColor(kBlack);
      hstyle->SetTitle("");
      hstyle->GetXaxis()->SetTitle(xTitle);
      hstyle->GetYaxis()->SetTitle(yTitle);
      hstyle->GetXaxis()->SetNdivisions(505);
      hstyle->GetXaxis()->SetLabelFont(42);
      hstyle->GetXaxis()->SetLabelOffset(0.007);
      hstyle->GetXaxis()->SetLabelSize(0.04);
      hstyle->GetXaxis()->SetTitleSize(0.06);
      hstyle->GetXaxis()->SetTitleOffset(0.9);
      hstyle->GetXaxis()->SetTitleFont(42);
      hstyle->GetYaxis()->SetNdivisions(505);
      hstyle->GetYaxis()->SetLabelFont(42);
      hstyle->GetYaxis()->SetLabelOffset(0.007);
      hstyle->GetYaxis()->SetLabelSize(0.04);
      hstyle->GetYaxis()->SetTitleSize(0.06);
      hstyle->GetYaxis()->SetTitleOffset(1.1);
      hstyle->GetYaxis()->SetTitleFont(42);

      double max_plot = htemplate->GetMaximum();

      for (int BSMSample=0; BSMSample<kNumGGHVVBSM; BSMSample++){
        if (iProd==0 && !(BSMSample<=kfLambda1_1 || BSMSample==kfLambda1_05)) continue;
        else if (iProd==1 && BSMSample>=kNumFiles_VBFHVVBSM) continue;
        else if (iProd==2 && BSMSample>=kNumFiles_WHVVBSM) continue;
        else if (iProd==3 && (BSMSample>=kNumFiles_ZHVVBSM || BSMSample==kfg2_1_fg4_0)) continue;
        else if (iProd==4 && BSMSample>0) continue;
        if (isBkg && BSMSample>0) continue;

        TString INPUT_NAME;
        if (!isBkg){
          INPUT_NAME = "AsimovTree_Sig_TxyBS_";
          if (iProd>0) INPUT_NAME.Append(Form("%s_", strHiggsProduction[iProd-1].Data()));
          if (BSMSample>0) INPUT_NAME.Append(Form("Hypo%i_", BSMSample));
          INPUT_NAME.Append(user_folder[folder]);
          INPUT_NAME.Append("_");
          INPUT_NAME.Append(comstring);
          INPUT_NAME.Append("_");
          INPUT_NAME = INPUT_NAME + cutNames[scheme];
          INPUT_NAME = INPUT_NAME + ".root";
        }
        else{
          INPUT_NAME = "AsimovTree_background_TxyBS_";
          INPUT_NAME.Append(user_folder[folder]);
          INPUT_NAME.Append("_");
          INPUT_NAME.Append(comstring);
          INPUT_NAME =  INPUT_NAME + "_" + cutNames[scheme] + ".root";
        }
        cinput = INPUT_NAME;
        cinput.Prepend(cinput_common);

        finput = new TFile(cinput, "read");
        if (finput->IsZombie() || finput==0){
          cout << "Could not find " << cinput << endl;
          if (finput!=0) delete finput;
          continue;
        }
        TString treename;
        if (!isBkg) treename = Form("CTau%.0f", target_ctau);
        else if (iProd==-1) treename = "ggZZ";
        else if (iProd==-2) treename = "qqZZ";
        else if (iProd==-3) treename = "ZX";
        treename.Prepend("T_2D_");
        foutput->cd();
        TTree* tin = 0;
        tin = (TTree*)finput->Get(treename);
        if (tin==0){
          finput->Close();
          continue;
        }
        hAsimov[BSMSample] = (TH1F*)hstyle->Clone(Form("%s_Asimov_BSM%i", htemplate->GetName(), BSMSample));
        hAsimov[BSMSample]->Reset("ICESM");
        TString strdraw = hAsimov[BSMSample]->GetName();
        strdraw.Prepend("CMS_zz4l_KD>>");
        TString strweight = "_weight_";
        tin->Draw(strdraw, strweight);
        delete tin;
        finput->Close();
        hAsimov[BSMSample] = (TH1F*)transform_tanh_varBin(hAsimov[BSMSample]);

        hAsimov[BSMSample]->SetName(Form("%s_%s_%s", hAsimov[BSMSample]->GetName(), comstring.Data(), channame[folder].Data()));
        foutput->WriteTObject(hAsimov[BSMSample]);
        cout << "BSM " << BSMSample << " integral: " << hAsimov[BSMSample]->Integral() << endl;
        if (BSMSample==0 || BSMSample>kfLambda1_1) hAsimov[BSMSample]->SetLineStyle(7);
        if (BSMSample==kfg2_1_fg4_0 || BSMSample==kfg2_05_fg4_0) hAsimov[BSMSample]->SetLineColor(kBlue);
        else if (BSMSample==kfg2_0_fg4_1 || BSMSample==kfg2_0_fg4_05) hAsimov[BSMSample]->SetLineColor(kRed);
        else if (BSMSample==kfLambda1_1 || BSMSample==kfLambda1_05) hAsimov[BSMSample]->SetLineColor(kViolet);
        else if (BSMSample!=0) cout << "\nUnknown color scheme for BSM code " << BSMSample << endl;

        if (hAsimov_Total[BSMSample]!=0) hAsimov_Total[BSMSample]->Add(hAsimov[BSMSample]);
        else hAsimov_Total[BSMSample] = (TH1F*)hAsimov[BSMSample]->Clone(Form("%s%s", hAsimov[BSMSample]->GetName(), "_total"));
        hAsimov[BSMSample]->Scale(1./hAsimov[BSMSample]->Integral());

        hAsimov_fill[BSMSample] = (TH1F*)hAsimov[BSMSample]->Clone(Form("%s_Fill", hAsimov[BSMSample]->GetName()));
        hAsimov_fill[BSMSample]->SetFillColor(hAsimov[BSMSample]->GetLineColor());
        hAsimov_fill[BSMSample]->SetFillStyle(3002);

        max_plot = max(max_plot, hAsimov[BSMSample]->GetMaximum());
      }
      if (iProd>0 && !isBkg && target_ctau==0 && htemplate_tth!=0) max_plot = max(max_plot, htemplate_tth->GetMaximum());
      htemplate->GetYaxis()->SetRangeUser(0, max_plot*1.3);
      if (target_ctau==0 && !isBkg && iProd>0) htemplate->GetYaxis()->SetRangeUser(1e-5, 10.);
      else if (isBkg || (target_ctau==0 && !isBkg && iProd==0 && folder!=1)) htemplate->GetYaxis()->SetRangeUser(1e-5, 10.);
      else if (target_ctau==0 && !isBkg && iProd==0 && folder==1) htemplate->GetYaxis()->SetRangeUser(1e-5, 10.);

      TPaveText* pt = new TPaveText(0.15, 0.93, 0.85, 1, "brNDC");
      pt->SetBorderSize(0);
      pt->SetFillStyle(0);
      pt->SetTextAlign(12);
      pt->SetTextFont(42);
      pt->SetTextSize(0.045);
      TText* text = pt->AddText(0.025, 0.45, "#font[61]{CMS}");
      text->SetTextSize(0.044);
      text = pt->AddText(0.165, 0.42, "#font[52]{Simulation}");
      text->SetTextSize(0.0315);
      TString cErgTev = Form("#font[42]{%.1f fb^{-1} (%i TeV)}", luminosity[EnergyIndex], erg_tev);
      if (EnergyIndex==0) cErgTev.Append(" ");
      text = pt->AddText(0.837, 0.45, cErgTev);
      text->SetTextSize(0.0315);

      foutput->cd();
      TString canvasname_2D = "cCompare_EmbeddedAsimov_Templates_TxyBS_";
      if (iProd>0) canvasname_2D.Append(Form("%s_", strHiggsProduction[iProd-1].Data()));
      else if (iProd==0) canvasname_2D.Append("ggH_");
      else if (iProd==-1 && isBkg) canvasname_2D.Append("ggZZ_");
      else if (iProd==-2 && isBkg) canvasname_2D.Append("qqZZ_");
      else if (iProd==-3 && isBkg) canvasname_2D.Append("ZX_");
      if (!isBkg) canvasname_2D.Append(Form("CTau%.0f_", target_ctau));
      canvasname_2D.Append(comstring);
      canvasname_2D.Append("_");
      canvasname_2D.Append(user_folder[folder]);
      TCanvas* c2D = new TCanvas(canvasname_2D, "", 8, 30, 800, 800);
      c2D->cd();
      gStyle->SetOptStat(0);
      c2D->SetFillColor(0);
      c2D->SetBorderMode(0);
      c2D->SetBorderSize(2);
      c2D->SetTickx(1);
      c2D->SetTicky(1);
      if (target_ctau==0 || isBkg) c2D->SetLogy();
      c2D->SetLeftMargin(0.17);
      c2D->SetRightMargin(0.05);
      c2D->SetTopMargin(0.07);
      c2D->SetBottomMargin(0.13);
      c2D->SetFrameFillStyle(0);
      c2D->SetFrameBorderMode(0);
      c2D->SetFrameFillStyle(0);
      c2D->SetFrameBorderMode(0);

      TLegend *l2D;
      if (!isBkg) l2D = new TLegend(0.20, 0.57, 0.58, 0.90);
      else l2D = new TLegend(0.20, 0.70, 0.58, 0.90);
      l2D->SetBorderSize(0);
      l2D->SetTextFont(42);
      l2D->SetTextSize(0.03);
      l2D->SetLineColor(1);
      l2D->SetLineStyle(1);
      l2D->SetLineWidth(1);
      l2D->SetFillColor(0);
      l2D->SetFillStyle(0);

      l2D->AddEntry(htemplate, "Template", "l");
      htemplate->Draw("hist");
      if (iProd>0 && iProd<4 && !isBkg && target_ctau==0 && htemplate_tth!=0){
        l2D->AddEntry(htemplate_tth, "Template (t#bar{t}H)", "l");
        htemplate_tth->Draw("histsame");
      }
      for (int BSMSample=0; BSMSample<kNumGGHVVBSM; BSMSample++){
        if (hAsimov[BSMSample]!=0){
          TString embeddedLabel;
          if (!isBkg) embeddedLabel = sample_GGHVVBSM_label[BSMSample];
          else embeddedLabel = "Bkg.";
          if (!(isBkg&&iProd==-3)) embeddedLabel.Append(" simulation");
          else embeddedLabel.Append(" sample");
          l2D->AddEntry(hAsimov[BSMSample], embeddedLabel, "l");
          hAsimov[BSMSample]->Draw("histsame");
          hAsimov_fill[BSMSample]->Draw("e2same");
        }
      }

      l2D->Draw();
      pt->Draw();

      TString folder_id;
      if (folder==0) folder_id = "4#mu";
      else if (folder==1) folder_id = "4e";
      else if (folder==2) folder_id = "2e2#mu";

      TPaveText* pt10;
      if (!isBkg) pt10 = new TPaveText(0.65, 0.84, 0.90, 0.92, "brNDC");
      else pt10 = new TPaveText(0.76, 0.84, 0.90, 0.92, "brNDC");
      pt10->SetBorderSize(0);
      pt10->SetTextAlign(12);
      pt10->SetTextSize(0.03);
      pt10->SetFillStyle(0);
      pt10->SetTextFont(42);
      TText* text10 = pt10->AddText(0.01, 0.01, folder_id);
      pt10->Draw();

      TPaveText* pt20;
      if (!isBkg) pt20 = new TPaveText(0.65, 0.80, 0.90, 0.84, "brNDC");
      else pt20 = new TPaveText(0.76, 0.80, 0.90, 0.84, "brNDC");
      pt20->SetBorderSize(0);
      pt20->SetTextAlign(12);
      pt20->SetTextSize(0.03);
      pt20->SetFillStyle(0);
      pt20->SetTextFont(42);
      TString productionLabel;
      if (!isBkg){
        productionLabel = strHiggsProduction_label[iProd];
        productionLabel.Append(Form(", c#tau_{H} = %.0f #mum", target_ctau));
      }
      else{
        if (iProd==-1) productionLabel = "gg bkg.";
        else if (iProd==-2) productionLabel = "q#bar{q} bkg.";
        else if (iProd==-3) productionLabel = "Z+X bkg.";
      }
      TText* text20 = pt20->AddText(0.01, 0.01, productionLabel);
      pt20->Draw();

      c2D->RedrawAxis();
      c2D->Modified();
      c2D->Update();
      foutput->WriteTObject(c2D);

      TString canvasDir = coutput_common;
      canvasDir.Append("ValidationPlots/");
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
      c2D->SaveAs(canvasname_2D_pdf);
      c2D->SaveAs(canvasname_2D_eps);
      c2D->SaveAs(canvasname_2D_png);
      c2D->SaveAs(canvasname_2D_root);
      c2D->SaveAs(canvasname_2D_c);

      delete pt20;
      delete pt10;
      delete l2D;
      c2D->Close();
      delete pt;

      for (int BSMSample=0; BSMSample<kNumGGHVVBSM; BSMSample++){
        if (hAsimov[BSMSample]!=0){
          delete hAsimov[BSMSample];
          delete hAsimov_fill[BSMSample];
        }
      }
      delete htemplate_tth;
      delete hstyle;
      delete htemplate;
      foutput->cd();
    }
  }

  if (htemplate_Total==0) cout << "WARNING! template total empty" << endl;
  htemplate_Total->Scale(1./htemplate_Total->Integral());
  if (htemplate_tth_Total==0) cout << "WARNING! template_tth total empty" << endl;
  if (htemplate_tth_Total!=0) htemplate_tth_Total->Scale(1./htemplate_tth_Total->Integral());
  for (int BSMSample=0; BSMSample<kNumGGHVVBSM; BSMSample++){
    if (hAsimov_Total[BSMSample]!=0) hAsimov_Total[BSMSample]->Scale(1./hAsimov_Total[BSMSample]->Integral());
  }

  TPaveText* pt = new TPaveText(0.15, 0.93, 0.85, 1, "brNDC");
  pt->SetBorderSize(0);
  pt->SetFillStyle(0);
  pt->SetTextAlign(12);
  pt->SetTextFont(42);
  pt->SetTextSize(0.045);
  TText* text = pt->AddText(0.025, 0.45, "#font[61]{CMS}");
  text->SetTextSize(0.044);
  text = pt->AddText(0.165, 0.42, "#font[52]{Simulation}");
  text->SetTextSize(0.0315);
  TString cErgTev = "#font[42]{19.7 fb^{-1} (8 TeV) + 5.1 fb^{-1} (7 TeV)}";
  text = pt->AddText(0.537, 0.45, cErgTev);
  text->SetTextSize(0.0315);

  foutput->cd();
  TString canvasname_2D = "cCompare_EmbeddedAsimov_Templates_TxyBS_";
  if (iProd>0) canvasname_2D.Append(Form("%s_", strHiggsProduction[iProd-1].Data()));
  else if (iProd==0) canvasname_2D.Append("ggH_");
  else if (iProd==-1 && isBkg) canvasname_2D.Append("ggZZ_");
  else if (iProd==-2 && isBkg) canvasname_2D.Append("qqZZ_");
  else if (iProd==-3 && isBkg) canvasname_2D.Append("ZX_");
  if (!isBkg) canvasname_2D.Append(Form("CTau%.0f_", target_ctau));
  canvasname_2D.Append("AllTeV");
  TCanvas* c2D = new TCanvas(canvasname_2D, "", 8, 30, 800, 800);
  c2D->cd();
  gStyle->SetOptStat(0);
  c2D->SetFillColor(0);
  c2D->SetBorderMode(0);
  c2D->SetBorderSize(2);
  c2D->SetTickx(1);
  c2D->SetTicky(1);
  if (target_ctau==0 || isBkg) c2D->SetLogy();
  c2D->SetLeftMargin(0.17);
  c2D->SetRightMargin(0.05);
  c2D->SetTopMargin(0.07);
  c2D->SetBottomMargin(0.13);
  c2D->SetFrameFillStyle(0);
  c2D->SetFrameBorderMode(0);
  c2D->SetFrameFillStyle(0);
  c2D->SetFrameBorderMode(0);

  double xl1 = 0.20;
  double xl2 = 0.58;
  double yl2 = 0.90;
  double ylwidth = (yl2-0.57) / 5.;
  TLegend* l2D = new TLegend(xl1, yl2-ylwidth*5., xl2, yl2);
  l2D->SetBorderSize(0);
  l2D->SetTextFont(42);
  l2D->SetTextSize(0.03);
  l2D->SetLineColor(1);
  l2D->SetLineStyle(1);
  l2D->SetLineWidth(1);
  l2D->SetFillColor(0);
  l2D->SetFillStyle(0);

  double max_global=0;
  double min_global=99999;
  for (int BSMSample=0; BSMSample<kNumGGHVVBSM; BSMSample++){
    if (hAsimov_Total[BSMSample]!=0){
      hAsimov_Total[BSMSample]->SetLineWidth(2);
      hAsimov_Total[BSMSample]->SetTitle("");
      hAsimov_Total[BSMSample]->GetXaxis()->SetTitle(xTitle);
      hAsimov_Total[BSMSample]->GetYaxis()->SetTitle(yTitle);
      hAsimov_Total[BSMSample]->GetXaxis()->SetNdivisions(505);
      hAsimov_Total[BSMSample]->GetXaxis()->SetLabelFont(42);
      hAsimov_Total[BSMSample]->GetXaxis()->SetLabelOffset(0.007);
      hAsimov_Total[BSMSample]->GetXaxis()->SetLabelSize(0.04);
      hAsimov_Total[BSMSample]->GetXaxis()->SetTitleSize(0.06);
      hAsimov_Total[BSMSample]->GetXaxis()->SetTitleOffset(0.9);
      hAsimov_Total[BSMSample]->GetXaxis()->SetTitleFont(42);
      hAsimov_Total[BSMSample]->GetYaxis()->SetNdivisions(505);
      hAsimov_Total[BSMSample]->GetYaxis()->SetLabelFont(42);
      hAsimov_Total[BSMSample]->GetYaxis()->SetLabelOffset(0.007);
      hAsimov_Total[BSMSample]->GetYaxis()->SetLabelSize(0.04);
      hAsimov_Total[BSMSample]->GetYaxis()->SetTitleSize(0.06);
      hAsimov_Total[BSMSample]->GetYaxis()->SetTitleOffset(1.1);
      hAsimov_Total[BSMSample]->GetYaxis()->SetTitleFont(42);

      if (BSMSample==0) hAsimov_Total[BSMSample]->SetLineStyle(1);
      max_global = max(max_global, hAsimov_Total[BSMSample]->GetMaximum());
      min_global = min(min_global, hAsimov_Total[BSMSample]->GetMinimum());
    }
  }

  TString productionLabel;
  if (!isBkg) productionLabel = strHiggsProduction_label[iProd];
  else{
    if (iProd==-1) productionLabel = "gg bkg.";
    else if (iProd==-2) productionLabel = "q#bar{q} bkg.";
    else if (iProd==-3) productionLabel = "Z+X bkg.";
  }

  if (iProd>0 && !isBkg && target_ctau==0 && htemplate_tth_Total!=0){
    max_global = max(max_global, htemplate_tth_Total->GetMaximum());
    min_global = min(min_global, htemplate_tth_Total->GetMinimum());
  }
  if (min_global<=0) min_global = 7e-6;

  int ctrDrawn=0;

  //l2D->AddEntry(htemplate_Total, "Template", "l");
  //htemplate_Total->Draw("hist");
  for (int BSMSample=0; BSMSample<kNumGGHVVBSM; BSMSample++){
    if (hAsimov_Total[BSMSample]!=0){
      TString embeddedLabel;
      if (!isBkg) embeddedLabel = sample_GGHVVBSM_label[BSMSample];
      else embeddedLabel = productionLabel;
      if (!(isBkg&&iProd==-3)) embeddedLabel.Append(Form(" %s", productionLabel.Data()));
      l2D->AddEntry(hAsimov_Total[BSMSample], embeddedLabel, "l");
//      hAsimov_Total[BSMSample]->Draw("histsame");
      if (target_ctau==0) symmetrize_PromptTemplates(hAsimov_Total[BSMSample]);
      if (BSMSample==0){
        if (target_ctau==0 || isBkg) hAsimov_Total[BSMSample]->GetYaxis()->SetRangeUser(min_global*0.7, max_global*5.);
        hAsimov_Total[BSMSample]->Draw("hist");
      }
      else hAsimov_Total[BSMSample]->Draw("histsame");
      ctrDrawn++;
    }
  }
  if (hAsimov_Total[0]!=0) hAsimov_Total[0]->Draw("histsame");
  if (iProd>0 && iProd<4 && !isBkg && target_ctau==0 && htemplate_tth_Total!=0){
    symmetrize_PromptTemplates((TH1F*) htemplate_tth_Total);
    l2D->AddEntry(htemplate_tth_Total, "SM t#bar{t}H", "l");
    htemplate_tth_Total->Draw("histsame");
    ctrDrawn++;
  }

  double yl1 = yl2-ylwidth*ctrDrawn;
  l2D->SetY1(yl1);

  l2D->Draw();
  pt->Draw();

  TPaveText* pt20=0;
  if (!isBkg){
    pt20 = new TPaveText(0.75, 0.84, 0.90, 0.92, "brNDC");
    pt20->SetBorderSize(0);
    pt20->SetTextAlign(12);
    pt20->SetTextSize(0.03);
    pt20->SetFillStyle(0);
    pt20->SetTextFont(42);

    TString label = Form("c#tau_{H} = %.0f #mum", target_ctau);
    TText* text20 = pt20->AddText(0.01, 0.01, label);
    pt20->Draw();
  }

  c2D->RedrawAxis();
  c2D->Modified();
  c2D->Update();
  foutput->WriteTObject(c2D);

  TString canvasDir = coutput_common;
  canvasDir.Append("ValidationPlots/");
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
  c2D->SaveAs(canvasname_2D_pdf);
  c2D->SaveAs(canvasname_2D_eps);
  c2D->SaveAs(canvasname_2D_png);
  c2D->SaveAs(canvasname_2D_root);
  c2D->SaveAs(canvasname_2D_c);

  if (pt20!=0) delete pt20;
  delete l2D;
  c2D->Close();
  delete pt;





  delete htemplate_Total;
  delete htemplate_tth_Total;
  for (int BSMSample=0; BSMSample<kNumGGHVVBSM; BSMSample++){
    if (hAsimov_Total[BSMSample]!=0){
      delete hAsimov_Total[BSMSample];
    }
  }

  foutput->Close();
}

TH1F* transform_tanh_varBin(TH1F* hinput){
  const int nbinsx = hinput->GetNbinsX()+1;
  double xdef[nbinsx];
  double xnew[nbinsx];

  for (int bin=1; bin<=nbinsx; bin++){
    xdef[bin-1] = hinput->GetXaxis()->GetBinLowEdge(bin);
  }
  for (int bin=0; bin<nbinsx; bin++){
    if (bin!=0 && bin!=nbinsx-1) xnew[bin] = TMath::ATanH(xdef[bin])*KD_scale;
    else if (bin==0) xnew[bin] = TMath::ATanH(xdef[bin+1])*KD_scale-200.;
    else xnew[bin] = TMath::ATanH(xdef[bin-1])*KD_scale+200.;
  }

  TString hname_core = hinput->GetName();
  TString hname = hname_core;
  hname.Append("_transformed");
  TString htitle = hinput->GetTitle();
  TString xtitle = hinput->GetXaxis()->GetTitle();
  TString ytitle = hinput->GetYaxis()->GetTitle();

  TH1F* houtput = new TH1F(hname,htitle,nbinsx-1,xnew);
  for (int bin=1; bin<=nbinsx-1; bin++){
    houtput->SetBinContent(bin, hinput->GetBinContent(bin));
    houtput->SetBinError(bin, hinput->GetBinError(bin));
  }
  houtput->GetXaxis()->SetTitle(xtitle);
  houtput->GetYaxis()->SetTitle(ytitle);
  houtput->SetLineStyle(hinput->GetLineStyle());
  houtput->SetLineColor(hinput->GetLineColor());
  houtput->SetLineWidth(hinput->GetLineWidth());
  houtput->SetMarkerStyle(hinput->GetMarkerStyle());
  houtput->SetMarkerColor(hinput->GetMarkerColor());
  houtput->SetFillStyle(hinput->GetFillStyle());
  houtput->SetFillColor(hinput->GetFillColor());
  delete hinput;
  houtput->SetName(hname_core);
  return houtput;
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



float compute_delDxy(
  float ZZPhi,
  float CandVtx_cov_xx,
  float CandVtx_cov_xy,
  float CandVtx_cov_xz,
  float CandVtx_cov_yy,
  float CandVtx_cov_yz,
  float CandVtx_cov_zz,
  float PrimaryVtx_cov_xx,
  float PrimaryVtx_cov_xy,
  float PrimaryVtx_cov_xz,
  float PrimaryVtx_cov_yy,
  float PrimaryVtx_cov_yz,
  float PrimaryVtx_cov_zz
  ){

  double unit_pt[3] ={
    cos(ZZPhi), sin(ZZPhi), 0
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
  bool nonInvertible = false;
  if (cand_cov[0] * cand_cov[4] * cand_cov[8] == 0){ // Determinant of triangular matrices = product of diagonal entries
    nonInvertible = true;
    cout << "Cand vtx determinant is 0!" << endl;
  }
  if (opv_cov[0] * opv_cov[4] * opv_cov[8] == 0){ // Determinant of triangular matrices = product of diagonal entries
    nonInvertible = true;
    cout << "OPV determinant is 0!" << endl;
  }
  if (!nonInvertible){
    TMatrixD candCov(3, 3);
    TMatrixD opvCov(3, 3);
    TVectorD unitVector(3, unit_pt);
    candCov.SetMatrixArray(cand_cov);
    opvCov.SetMatrixArray(opv_cov);

    TMatrixD candInvCov = candCov;
    TMatrixD opvInvCov = opvCov;

    double testdet = 0;
    candInvCov.Invert(&testdet);
    if (testdet == 0){ cout << "CandVtx inversion unsuccessful" << endl; nonInvertible = true; }
    testdet = 0;
    opvInvCov.Invert(&testdet);
    if (testdet == 0){ cout << "OPV inversion unsuccessful" << endl; nonInvertible = true; }
    TMatrixD intInvCov = opvInvCov;
    intInvCov += candInvCov;
    TMatrixD intCov = intInvCov;
    testdet = 0;
    intCov.Invert(&testdet);
    if (testdet == 0){ cout << "Int. vtx inv. cov. inversion unsuccessful" << endl; nonInvertible = true; }

    if (!nonInvertible){
      TVectorD intVector = unitVector;
      intVector *= intCov;
      double dot_product = intVector*unitVector;
      return sqrt(dot_product);
    }
    else{
      return 0;
    }
  }
  else{
    return 0;
  }
}

