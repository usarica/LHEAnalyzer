#include <iostream>
#include <iomanip>
#include <cmath>
#include <string>
#include <vector>
#include <fstream>
#include <cstdlib>
#include "TFile.h"
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
#include "RooAbsPdf.h"
#include "RooRealSumPdf.h"
#include "RooRealVar.h"
#include "RooFormulaVar.h"
#include "RooDataHist.h"
#include "RooHistFunc.h"
#include "RooMsgService.h"
#include "RooFitResult.h"
#include "./Pdfs/RooRealFlooredSumPdf.h"
#include "./data/ZZ4l_125_6_Samples.h"

using namespace std;
using namespace ROOT::Math;
using namespace RooFit;


const double ctau_limits[2] = { 0, 1000 };
const int nctaus = 100;
const double cm_to_microns = 10000.;
//const double KD_scale = 800.;
const double KD_scale = 1250.;

string processName[8] ={
	"CTau0", "CTau100", "CTau500", "CTau1000",
	"qqZZ", "ggZZ", "CR", "data"
};

const int nHiggsProd = 5;
TString strHiggsProductionFit[nHiggsProd-1] ={ "VBFH125", "WH125", "ZH125", "ttH125" };
TString strHiggsProduction[nHiggsProd-1] ={ "VBFH", "WH", "ZH", "ttH" };

const int newSIPCuts = 4; // 0: No SIP; 1: Z1SIP, 2: KalmanChi2, 3: Z1SIP + KalmanChi2
TString cutNames[newSIPCuts + 1] ={
	"PVSIPcut",
	"NoSIPcut",
	"Z1SIPcut",
	"4lchi2cut",
	"Z1SIP_4lchi2cut"
};

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


void makeCombineSignalTemplateswithTrees_Signal_single(int folder, int erg_tev, int iProd, int scheme=4);
void makeCombineSignalTemplateswithTrees_Bkg_single(int folder, int erg_tev, int Systematics, int scheme=4);
void makeCombineSignalTemplateSyst_Signal_single(int folder, int erg_tev, int iProd, int scheme);
void regularizeKDSyst(TH1F* txy_syst[3]);


void makeCombineTemplateswithTrees_NoSIP(int recreateSyst=0, int selection=4){
//	const int kNumSyst=7;
//	int systematics[kNumSyst]={0,1,-1,2,-2,3,-3};
  if (recreateSyst==1){
    for (int CoM=7; CoM<9; ++CoM){
      for (int channel=0; channel<3; channel++){
        for (int iProd=0; iProd<nHiggsProd; iProd++) makeCombineSignalTemplateSyst_Signal_single(channel, CoM, iProd, selection);
      }
    }
  }
	for(int CoM=7;CoM<9;++CoM){
		for(int channel=0;channel<3;channel++){
      for (int iProd=0; iProd<nHiggsProd; iProd++) makeCombineSignalTemplateswithTrees_Signal_single(channel, CoM, iProd, selection);
			for(int syst=-1; syst<=1; syst++) makeCombineSignalTemplateswithTrees_Bkg_single(channel, CoM, syst, selection);
		}
	}
}

float get_KD_Syst(float val, float err, int kSystVar=0){
  if (kSystVar==0) return 1;
  else{
    if (val>=1){
      if (kSystVar>0) return val+err;
      else if (kSystVar<0) return 1./(val+err);
    }
    else{
      if (kSystVar<0) return val-err;
      else if (kSystVar>0) return 1./(val-err);
    }
  }

}

float compute_KDAlternatives(float KD, float KD_true, float systVar){
  return (KD-KD_true)*systVar+KD_true;
}

void makeCombineSignalTemplateswithTrees_Signal_single(int folder, int erg_tev, int iProd, int scheme){
  if (scheme==0) return;

  bool applySpecialFill = (iProd==4 && erg_tev==7 && folder==1);

  char TREE_NAME[]="SelectedTree";
	TString INPUT_NAME = "HZZ4lTree_powheg15jhuGenV3-CTau";
  TString OUTPUT_NAME = "_templates_Sig_TxyBS_";
  if (iProd>0) OUTPUT_NAME.Append(Form("%s_", strHiggsProduction[iProd-1].Data()));
  OUTPUT_NAME = OUTPUT_NAME + cutNames[scheme] + ".root";
  OUTPUT_NAME.Prepend(user_folder[folder]);
	TString comstring;
	comstring.Form("%iTeV", erg_tev);
	TString erg_dir;
	erg_dir.Form("LHC_%iTeV", erg_tev);

	const int nSyst=7;
	int kSyst=nSyst;
	const int kSystNominal=3;
	const int iSystRes=1;
	const int iSystScale=2;
	const int iSystKD=3;
	char* cSyst[nSyst] ={
		"_TxyDown", "_ScaleDown", "_ResDown",
		"",
		"_ResUp", "_ScaleUp", "_TxyUp"
	};

	int EnergyIndex=1;
	if (erg_tev==7) EnergyIndex=0;
	float ZZMass_PeakCut[2]={ 105.6, 140.6 };

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
  //use ggH yields directly from SIP<4 case
  double targetYield = yield_signal_ggh[EnergyIndex][folder]*ratio_targetsignalyield/luminosity[EnergyIndex];

  TFile* finput_pTrewgt = new TFile(Form("./data/compareProductionModes_pT_%iTeV%s", erg_tev, ".root"), "read");
  TF1* fit_pT = 0;
  if (iProd>0){
    TString fitname = strHiggsProductionFit[iProd-1];
    fitname.Append("_pToverMZZ_ratio_fit");
    fit_pT = (TF1*)finput_pTrewgt->Get(fitname);
    if (iProd==1) targetYield  = yield_signal_vbfh[EnergyIndex][folder]*ratio_targetsignalyield/luminosity[EnergyIndex];
    else if (iProd==2) targetYield  = yield_signal_wh[EnergyIndex][folder]*ratio_targetsignalyield/luminosity[EnergyIndex];
    else if (iProd==3) targetYield  = yield_signal_zh[EnergyIndex][folder]*ratio_targetsignalyield/luminosity[EnergyIndex];
    else if (iProd==4) targetYield  = yield_signal_tth[EnergyIndex][folder]*ratio_targetsignalyield/luminosity[EnergyIndex];
  }
  cout << comstring << ' ' << user_folder[folder] << " yield: " << targetYield*luminosity[EnergyIndex] << endl;


	float templateWeight[nSyst] ={ 1 };
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
  float KD=-99, KD_up=-99, KD_dn=-99;
  float Txy, Dxy, Txy_true, Dxy_true, delTxy, delDxy;
  float Txy_BS, Dxy_BS, Txy_BS_true, Dxy_BS_true;
  float delTxy_BS=-99;
  float sigmaPV_xy, sigmaInt_xy;

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

	float D_bkg=-99;
	float D_bkg_ScaleUp=-99;
	float D_bkg_ScaleDown=-99;
	float D_bkg_ResUp=-99;
	float D_bkg_ResDown=-99;

	double ctau_gridsize = (ctau_limits[1]-ctau_limits[0])/nctaus;
	//	int nbinsx = 400;
	//	double xlow = -10000, xhigh = 10000;
	int nbinsx = 50;
	double xlow = -1, xhigh = 1;
	int nbinsy = 50;
	double ylow=0, yhigh=1;

  TString cinput_common = user_dir_hep + "No_SIP/" + erg_dir + "/" + user_folder[folder] + "/";
  TString coutput_common = user_dir_hep + "Analysis/Templates/" + comstring + "/";

  TString INPUT_SYST_NAME = "_templateSystematics_Sig_TxyBS_";
  if (iProd>0) INPUT_SYST_NAME.Append(Form("%s_", strHiggsProduction[iProd-1].Data()));
  INPUT_SYST_NAME = INPUT_SYST_NAME + cutNames[scheme] + ".root";
  INPUT_SYST_NAME.Prepend(user_folder[folder]);
  INPUT_SYST_NAME.Prepend(coutput_common);
  TFile* finput_syst = new TFile(INPUT_SYST_NAME, "read");

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

	TChain* tree[kNumSamples] ={ 0 };
	for (int f=0; f<kNumSamples; f++){
		TString cinput = cinput_common + sample_FullSim[f] + ".root";;
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
		double Txy_BS_systWeight = 1.;
    bool hasSpecialFill = applySpecialFill;
    bool hasSpecialFill_v2 = (target_ctau==0);

    TH1F* txy_syst[3]={ 0 };
    if(genctau>0){
      txy_syst[0]=(TH1F*)finput_syst->Get(Form("fit_CTau%.0f__KD", target_ctau));
      txy_syst[1]=(TH1F*)finput_syst->Get(Form("fit_CTau%.0f__KD_up", target_ctau));
      txy_syst[2]=(TH1F*)finput_syst->Get(Form("fit_CTau%.0f__KD_dn", target_ctau));

//      cout << "Extracted Txy_BS syst histogram " << Form("fit_CTau%.0f__KD", target_ctau) << endl;
      for (int st=0; st<3; st++){
        if (txy_syst[st]==0) cout << "txy_syst[" << st << "] empty" << endl;
      }
    }
    else{
      txy_syst[0]=(TH1F*)finput_syst->Get(Form("H_2D_CTau%.0f_px", target_ctau));
      txy_syst[1]=(TH1F*)finput_syst->Get(Form("H_2D_CTau%.0f_TxyUp_px", target_ctau));
      txy_syst[2]=(TH1F*)finput_syst->Get(Form("H_2D_CTau%.0f_TxyDown_px", target_ctau));

//      cout << "Extracted Txy_BS syst histogram " << Form("H_2D_CTau%.0f_px", target_ctau) << endl;
      for (int st=0; st<3; st++){
        if (txy_syst[st]==0) cout << "txy_syst[" << st << "] empty" << endl;
        else cout << "Txy BS syst " << st << " norm: " << txy_syst[st]->Integral() << endl;
      }
    }
    regularizeKDSyst(txy_syst);

		TH2F* hfill[nSyst];
		for (int ss = 0; ss < nSyst; ss++){
			hfill[ss] = new TH2F(Form("H_2D_CTau%.0f%s", target_ctau, cSyst[ss]), "", nbinsx, xlow, xhigh, nbinsy, ylow, yhigh);
			hfill[ss]->SetTitle(Form("%.0f #mum %s at %i TeV", target_ctau, user_folder[folder], erg_tev));
      hfill[ss]->GetXaxis()->SetTitle(Form("tanh(T_{xy} / %.0f mm)", KD_scale));
			hfill[ss]->GetYaxis()->SetTitle("#it{D}_{bkg}");
		}

    TTree* mytree = new TTree(Form("T_2D_CTau%.0f", target_ctau), "");
		mytree->Branch("kSystematics", &kSyst);
		mytree->Branch("MC_weight", templateWeight, "MC_weight[kSystematics]/F");
		mytree->Branch("D_bkg", &D_bkg);
		mytree->Branch("D_bkg_ScaleUp", &D_bkg_ScaleUp);
		mytree->Branch("D_bkg_ScaleDown", &D_bkg_ScaleDown);
		mytree->Branch("D_bkg_ResUp", &D_bkg_ResUp);
		mytree->Branch("D_bkg_ResDown", &D_bkg_ResDown);
    mytree->Branch("KD", &KD);
//    mytree->Branch("KD_up", &KD_up);
//    mytree->Branch("KD_dn", &KD_dn);
    mytree->Branch("GenCTau", &Txy_true);
		mytree->Branch("fileID", &fileID);

		TH2F* hCounters = new TH2F(Form("hCounters_2D_CTau%.0f", target_ctau), "", kNumSamples, 0, kNumSamples, 2, 0, 2);
		hCounters->SetOption("colz");
		for (int f = 0; f < kNumSamples; f++){ // DO NOT APPLY CUTS ON PURPOSE
			if (f == 0 && genctau>0) continue;
			if (f > 0 && genctau==0) break;
			double initial_ctau = sample_CTau[f];
			double effectiveInvCTau = (f>0 ? 1. / target_ctau - 1. / initial_ctau : 0);

			for (int ev = 0; ev < tree[f]->GetEntries(); ev++){
				GenPrimaryVtx_x=0;
				GenPrimaryVtx_y=0;
				GenPrimaryVtx_z=0;
				GenIntVtx_x=0;
				GenIntVtx_y=0;
				GenIntVtx_z=0;

				GenHMass=125.6;
				GenHPt=20;
				GenHPhi=0;

				tree[f]->GetEntry(ev);

				TrueCandVtx_x = GenIntVtx_x - GenPrimaryVtx_x;
				TrueCandVtx_y = GenIntVtx_y - GenPrimaryVtx_y;
				TrueCandVtx_z = GenIntVtx_z - GenPrimaryVtx_z;
				Dxy_true = (TrueCandVtx_x*cos(GenHPhi) + TrueCandVtx_y*sin(GenHPhi)) * cm_to_microns;
				Txy_true = Dxy_true*GenHMass / GenHPt;
        if (Txy_true<0) continue;

				templateWeight[kSystNominal] = 1;
				ctauWeight = exp(-Txy_true*effectiveInvCTau);
				templateWeight[kSystNominal] *= ctauWeight;
        if (fit_pT!=0 && iProd>0) templateWeight[kSystNominal] *= fit_pT->Eval(GenHPt/GenHMass);

				hCounters->AddBinContent(hCounters->GetBin(f+1, 1), 1.);
				hCounters->AddBinContent(hCounters->GetBin(f+1, 2), templateWeight[kSystNominal]);
				templateWeight[kSystNominal]=1;
			}
		}
		foutput->WriteTObject(hCounters);

		for (int f = 0; f < kNumSamples; f++){
			if (f == 0 && genctau>0) continue;
			if (f > 0 && genctau==0) continue;
			double initial_ctau = sample_CTau[f];
			if (initial_ctau<target_ctau) continue;
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

				float smd[nSyst] ={
					D_bkg, D_bkg_ScaleDown, D_bkg_ResDown,
					D_bkg,
					D_bkg_ResUp, D_bkg_ScaleUp, D_bkg
				};
//        float KDSystwgt_up = get_KD_Syst(smearingVal, smearingErr, 1);
//        float KDSystwgt_dn = get_KD_Syst(smearingVal, smearingErr, -1);
//        float KDSystwgt_nominal = get_KD_Syst(smearingVal, smearingErr, 0);

//        KD = compute_KDAlternatives(Txy_BS, Txy_BS_true, KDSystwgt_nominal);
        KD = Txy_BS;
        KD = tanh(KD / KD_scale); // Scaled to due to tanh function
        if (KD >= 1.) KD = 1. - 1.0e-4;
        if (KD <= -1.) KD = -1. + 1.0e-4;
/*
        KD_up = compute_KDAlternatives(Txy_BS, Txy_BS_true, KDSystwgt_up);
        KD_up = tanh(KD_up / KD_scale); // Scaled to due to tanh function
        if (KD_up >= 1.) KD_up = 1. - 1.0e-10;
        if (KD_up <= -1.) KD_up = -1. + 1.0e-10;

        KD_dn = compute_KDAlternatives(Txy_BS, Txy_BS_true, KDSystwgt_dn);
        KD_dn = tanh(KD_dn / KD_scale); // Scaled to due to tanh function
        if (KD_dn >= 1.) KD_dn = 1. - 1.0e-10;
        if (KD_dn <= -1.) KD_dn = -1. + 1.0e-10;
*/
        for (int ss = 0; ss < nSyst; ss++){
					templateWeight[ss] = MC_weight;
					ctauWeight = exp(-Txy_true*effectiveInvCTau);
					templateWeight[ss] *= ctauWeight*xsec_ctau_wgt;
          if (fit_pT!=0 && iProd>0) templateWeight[ss] *= fit_pT->Eval(GenHPt/GenHMass);

					if (smd[ss]<0) templateWeight[ss] = 0;
          double KDSystWgt=1;
          int txy_syst_bin = txy_syst[0]->GetXaxis()->FindBin(Txy_BS);
          if (txy_syst_bin > txy_syst[0]->GetNbinsX()) txy_syst_bin = txy_syst[0]->GetNbinsX();
          if (txy_syst_bin < 1) txy_syst_bin = 1;
          if (ss == (kSystNominal + iSystKD)) KDSystWgt = txy_syst[1]->GetBinContent(txy_syst_bin);
          else if (ss == (kSystNominal - iSystKD)) KDSystWgt = txy_syst[2]->GetBinContent(txy_syst_bin);
          if (KDSystWgt>1.6) KDSystWgt=1.6;
          if (KDSystWgt<0.6) KDSystWgt=0.6;
          if (KDSystWgt!=KDSystWgt) KDSystWgt = 1;
          templateWeight[ss] *= KDSystWgt;
//          if ((ss == (kSystNominal + iSystKD) || ss == (kSystNominal - iSystKD)) && Txy_BS>1000 && target_ctau==600.) cout << Txy_BS << '\t' << KDSystWgt << endl;

//          if (ss == (kSystNominal + iSystKD)) hfill[ss]->Fill(KD_up, smd[ss], templateWeight[ss]);
//          else if (ss == (kSystNominal - iSystKD)) hfill[ss]->Fill(KD_dn, smd[ss], templateWeight[ss]);
//          else hfill[ss]->Fill(KD, smd[ss], templateWeight[ss]);
          hfill[ss]->Fill(KD, smd[ss], templateWeight[ss]);
        }
			}
		}
		double integral_rescale[nSyst] ={ 1 };
		for (int ss = 0; ss < nSyst; ss++){
			hfill[ss]->SetOption("colz");
			double integral_rewgt = hfill[ss]->Integral(0,nbinsx+1,0,nbinsy+1);
			double integral_target = targetYield;
			integral_rescale[ss] = integral_target / integral_rewgt;
			cout << "Template for cTau = " << target_ctau << " of " << user_folder[folder] << " @ " << comstring << " is rescaled with " << integral_rescale[ss] << " for systematic " << ss << "..." << endl;
			hfill[ss]->Scale(integral_rescale[ss]);
		}

    double sum_filled[nSyst]={ 0 };
    for (int f = 0; f < kNumSamples; f++){
			if (f == 0 && genctau>0) continue;
			if (f > 0 && genctau==0) continue;
			double initial_ctau = sample_CTau[f];
			if (initial_ctau<target_ctau) continue;
			double effectiveInvCTau = (f>0 ? 1. / target_ctau - 1. / initial_ctau : 0);
//			double xsec_ctau_wgt = hCounters->GetBinContent(f+1, 1) / hCounters->GetBinContent(f+1, 2);
			double xsec_ctau_wgt = (initial_ctau>0 ? initial_ctau / target_ctau : 1);
			fileID = f;

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

        float smd[nSyst] ={
          D_bkg, D_bkg_ScaleDown, D_bkg_ResDown,
          D_bkg,
          D_bkg_ResUp, D_bkg_ScaleUp, D_bkg
        };
        //        float KDSystwgt_up = get_KD_Syst(smearingVal, smearingErr, 1);
        //        float KDSystwgt_dn = get_KD_Syst(smearingVal, smearingErr, -1);
        //        float KDSystwgt_nominal = get_KD_Syst(smearingVal, smearingErr, 0);

//        KD = compute_KDAlternatives(Txy_BS, Txy_BS_true, KDSystwgt_nominal);
        KD = Txy_BS;
        KD = tanh(KD / KD_scale); // Scaled to due to tanh function
        if (KD >= 1.) KD = 1. - 1.0e-4;
        if (KD <= -1.) KD = -1. + 1.0e-4;
        /*
        KD_up = compute_KDAlternatives(Txy_BS, Txy_BS_true, KDSystwgt_up);
        KD_up = tanh(KD_up / KD_scale); // Scaled to due to tanh function
        if (KD_up >= 1.) KD_up = 1. - 1.0e-10;
        if (KD_up <= -1.) KD_up = -1. + 1.0e-10;

        KD_dn = compute_KDAlternatives(Txy_BS, Txy_BS_true, KDSystwgt_dn);
        KD_dn = tanh(KD_dn / KD_scale); // Scaled to due to tanh function
        if (KD_dn >= 1.) KD_dn = 1. - 1.0e-10;
        if (KD_dn <= -1.) KD_dn = -1. + 1.0e-10;
        */
        for (int ss = 0; ss < nSyst; ss++){
          templateWeight[ss] = MC_weight;
					ctauWeight = exp(-Txy_true*effectiveInvCTau);
					templateWeight[ss] *= ctauWeight*xsec_ctau_wgt;
					templateWeight[ss] *= integral_rescale[ss];
					if (smd[ss]<0) templateWeight[ss] = 0;

          double KDSystWgt=1;
          int txy_syst_bin = txy_syst[0]->GetXaxis()->FindBin(Txy_BS);
          if (txy_syst_bin > txy_syst[0]->GetNbinsX()) txy_syst_bin = txy_syst[0]->GetNbinsX();
          if (txy_syst_bin < 1) txy_syst_bin = 1;
          if (ss == (kSystNominal + iSystKD)) KDSystWgt = txy_syst[1]->GetBinContent(txy_syst_bin);
          else if (ss == (kSystNominal - iSystKD)) KDSystWgt = txy_syst[2]->GetBinContent(txy_syst_bin);
          if (KDSystWgt>1.6) KDSystWgt=1.6;
          if (KDSystWgt<0.6) KDSystWgt=0.6;
          if (KDSystWgt!=KDSystWgt) KDSystWgt = 1;

          templateWeight[ss] *= KDSystWgt;
          if (fit_pT!=0 && iProd>0) templateWeight[ss] *= fit_pT->Eval(GenHPt/GenHMass);
          sum_filled[ss] += templateWeight[ss];
          if (hasSpecialFill) templateWeight[ss]*=0.5;
          if (hasSpecialFill_v2) templateWeight[ss]*=0.5;
        }
        int fill_Iterator = 1;
        if (hasSpecialFill) fill_Iterator = 2;
        for (int itFill=0; itFill<fill_Iterator; itFill++){
          mytree->Fill();
          if (hasSpecialFill_v2){
            KD = -KD;
            mytree->Fill();
          }
        }
      }
		}
    cout << "Final recorded yield for ";
    for (int ss = 0; ss < nSyst; ss++){
      cout << "systematic " << ss << ": " << sum_filled[ss] << " | ";
    }
    cout << endl;
    cout << "Filled the tree" << endl;
    TH1F* hprojX_array[nSyst];
    TH2F* hconditional_array[nSyst];
    for (int ss = 0; ss < nSyst; ss++){
      TH2F* hfill_conditional;
      TH1F* hfill_projX;

      hfill_conditional = (TH2F*)hfill[ss]->Clone(Form("%s_Conditional", hfill[ss]->GetName()));
      hfill_conditional->SetTitle(Form("%s Conditional over tanh(T_{xy} / %.0f #mum)", hfill[ss]->GetTitle(), KD_scale));
      /*			for (int biny = 1; biny <= hfill_conditional->GetNbinsY(); biny++){
              double integralX = hfill_conditional->Integral(0, hfill_conditional->GetNbinsX() + 1, biny, biny);
              if (integralX != 0){
              for (int binx = 1; binx <= hfill_conditional->GetNbinsX(); binx++) hfill_conditional->SetBinContent(binx, biny, (hfill_conditional->GetBinContent(binx, biny) / integralX));
              }
              }
              */
      for (int binx = 1; binx <= hfill_conditional->GetNbinsX(); binx++){
        double integralX = hfill_conditional->Integral(binx, binx, 0, hfill_conditional->GetNbinsY() + 1);
        if (integralX != 0){
          for (int biny = 1; biny <= hfill_conditional->GetNbinsY(); biny++) hfill_conditional->SetBinContent(binx, biny, (hfill_conditional->GetBinContent(binx, biny) / integralX));
        }
      }

      hfill_projX = (TH1F*)hfill[ss]->ProjectionX();
      hfill_projX->SetNameTitle(Form("H_1DTxy_BS_CTau%.0f%s", target_ctau, cSyst[ss]),Form("%s Projection on tanh(T_{xy} / %.0f #mum)", hfill[ss]->GetTitle(), KD_scale));

      foutput->WriteTObject(hfill[ss]);
      foutput->WriteTObject(hfill_conditional);
      foutput->WriteTObject(hfill_projX);
      hprojX_array[ss] = hfill_projX;
      hconditional_array[ss] = hfill_conditional;
    }
    TH1F* hKDSystRatio[2];
    hKDSystRatio[0] = (TH1F*)hprojX_array[(kSystNominal + iSystKD)]->Clone(Form("%s_Ratio", hprojX_array[(kSystNominal + iSystKD)]->GetName()));
    hKDSystRatio[1] = (TH1F*)hprojX_array[(kSystNominal - iSystKD)]->Clone(Form("%s_Ratio", hprojX_array[(kSystNominal - iSystKD)]->GetName()));
    for (int rr=0; rr<2; rr++){
      hKDSystRatio[rr]->Divide(hprojX_array[kSystNominal]);
      foutput->WriteTObject(hKDSystRatio[rr]);
      delete hKDSystRatio[rr];
    }
    for (int ss = 0; ss < nSyst; ss++){
      delete hprojX_array[ss];
      delete hconditional_array[ss];
			delete hfill[ss];
		}
		foutput->WriteTObject(mytree);
		delete mytree;
    for (int st=0; st<3; st++) delete txy_syst[st];

	}
	for (int f=0; f<kNumSamples; f++){ if (tree[f]!=0) delete tree[f]; }

  delete tBeam;
	foutput->Close();
	cout << "Closed " << coutput << endl;
	finput_syst->Close();
  if (fit_pT!=0) delete fit_pT;
  finput_pTrewgt->Close();
}

TH2F** makeCombineSignalTemplateMCSyst_Signal_single(int folder, int erg_tev, int iProd, int scheme){
  if (scheme==0) return 0;
  gStyle->SetOptStat(0);
  TH2F** hArray = new TH2F*[nctaus+1];

  char TREE_NAME[]="SelectedTree";
  TString INPUT_NAME = "HZZ4lTree_powheg15jhuGenV3-CTau";
  TString comstring;
  comstring.Form("%iTeV", erg_tev);
  TString erg_dir;
  erg_dir.Form("LHC_%iTeV", erg_tev);

  int EnergyIndex=1;
  if (erg_tev==7) EnergyIndex=0;
  float ZZMass_PeakCut[2]={ 105.6, 140.6 };

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
  //use ggH yields directly from SIP<4 case
  double targetYield = yield_signal_ggh[EnergyIndex][folder]*ratio_targetsignalyield/luminosity[EnergyIndex];

  TFile* finput_pTrewgt = new TFile(Form("./data/compareProductionModes_pT_%iTeV%s", erg_tev, ".root"), "read");
  TF1* fit_pT = 0;
  if (iProd>0){
    TString fitname = strHiggsProductionFit[iProd-1];
    fitname.Append("_pToverMZZ_ratio_fit");
    fit_pT = (TF1*)finput_pTrewgt->Get(fitname);
    if (iProd==1) targetYield  = yield_signal_vbfh[EnergyIndex][folder]*ratio_targetsignalyield/luminosity[EnergyIndex];
    else if (iProd==2) targetYield  = yield_signal_wh[EnergyIndex][folder]*ratio_targetsignalyield/luminosity[EnergyIndex];
    else if (iProd==3) targetYield  = yield_signal_zh[EnergyIndex][folder]*ratio_targetsignalyield/luminosity[EnergyIndex];
    else if (iProd==4) targetYield  = yield_signal_tth[EnergyIndex][folder]*ratio_targetsignalyield/luminosity[EnergyIndex];
  }
  cout << comstring << ' ' << user_folder[folder] << " yield: " << targetYield*luminosity[EnergyIndex] << endl;
  
  float templateWeight = 1;
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
  float KD=-99, KD_up=-99, KD_dn=-99;
  float Txy, Dxy, Txy_true, Dxy_true, delTxy, delDxy;
  float Txy_BS, Dxy_BS, Txy_BS_true, Dxy_BS_true;
  float delTxy_BS=-99;
  float sigmaPV_xy, sigmaInt_xy;

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

  float D_bkg=-99;
  float D_bkg_ScaleUp=-99;
  float D_bkg_ScaleDown=-99;
  float D_bkg_ResUp=-99;
  float D_bkg_ResDown=-99;

  double ctau_gridsize = (ctau_limits[1]-ctau_limits[0])/nctaus;
  int nbinsx = 30;
  double xlow = -2340, xhigh = 2340;
  double binwidthx = (xhigh-xlow)/nbinsx;
  xlow-=binwidthx*1.5; xhigh+=binwidthx*1.5; nbinsx+=3;
  int nbinsy = 50;
  double ylow=0, yhigh=1;

  TString INPUT_SYST_NAME = "LifetimeKD_SystFits_DataMC_";
  INPUT_SYST_NAME.Append("newSIP_");
  INPUT_SYST_NAME.Append(user_folder[folder]);
  INPUT_SYST_NAME.Append(Form("_%s_", "Txy_BeamSpot"));
  INPUT_SYST_NAME.Append("AllTeV");
  INPUT_SYST_NAME.Append(".root");
  TString cinput_syst = user_dir_hep + "Analysis/ShapeSystematics/";
  cinput_syst.Append("Txy_BeamSpot/");
  cinput_syst.Append(INPUT_SYST_NAME);
  TFile* finput_syst = new TFile(cinput_syst, "read");
  TH1F* txy_syst = (TH1F*)finput_syst->Get("hCommonSmearing");
  float smearingVal = txy_syst->GetBinContent(1);
  float smearingErr = txy_syst->GetBinError(1);

  // BeamSpot info
  TString strBeamSpot[2]={ "MC_START44_V13" , "MC_START53_V23" };
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

  TString cinput_common = user_dir_hep + "No_SIP/" + erg_dir + "/" + user_folder[folder] + "/";
  TChain* tree[kNumSamples] ={ 0 };
  for (int f=0; f<kNumSamples; f++){
    TString cinput = cinput_common + sample_FullSim[f] + ".root";;
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
    double Txy_BS_systWeight = 1.;

    gROOT->cd();
    gStyle->SetOptStat(0);
    TH2F* hfill;
    hfill = new TH2F(Form("H_2D_CTau%.0f_MC", target_ctau), "", nbinsx, xlow, xhigh, nbinsy, ylow, yhigh);
    hfill->SetTitle("");
    hfill->GetXaxis()->SetTitle("T^{Beam}_{xy} (#mum)");
    hfill->GetYaxis()->SetTitle("#it{D}_{bkg}");
    hfill->GetZaxis()->SetTitle("Events / bin");

    for (int f = 0; f < kNumSamples; f++){
      if (f == 0 && genctau>0) continue;
      if (f > 0 && genctau==0) continue;
      double initial_ctau = sample_CTau[f];
      if (initial_ctau<target_ctau) continue;
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

        float smd = D_bkg;
        float KDSystwgt_nominal = get_KD_Syst(smearingVal, smearingErr, 0);
        KD = compute_KDAlternatives(Txy_BS, Txy_BS_true, KDSystwgt_nominal);

        templateWeight = MC_weight;
        ctauWeight = exp(-Txy_true*effectiveInvCTau);
        templateWeight *= ctauWeight*xsec_ctau_wgt;
        if (fit_pT!=0 && iProd>0) templateWeight *= fit_pT->Eval(GenHPt/GenHMass);
        hfill->Fill(KD, smd, templateWeight);
      }
    }
    double integral_rescale ={ 1 };
    hfill->SetOption("colz");
    double integral_rewgt = hfill->Integral(0, nbinsx+1, 0, nbinsy+1);
    double integral_target = targetYield;
    integral_rescale = integral_target / integral_rewgt;
    hfill->Scale(integral_rescale);
    for (int biny=0; biny<=hfill->GetNbinsY()+1; biny++){
      double underflow = hfill->GetBinContent(0, biny);
      double overflow = hfill->GetBinContent(hfill->GetNbinsX()+1, biny);
      hfill->AddBinContent(hfill->GetBin(1, biny), underflow);
      hfill->AddBinContent(hfill->GetBin(hfill->GetNbinsX(), biny), overflow);
      hfill->SetBinContent(0, biny, 0);
      hfill->SetBinContent(hfill->GetNbinsX()+1, biny, 0);
    }

    hfill->SetStats(0);
    hArray[genctau] = hfill;
  }

  for (int f=0; f<kNumSamples; f++){ if (tree[f]!=0) delete tree[f]; }
  delete tBeam;
  delete txy_syst;
  finput_syst->Close();
  if (fit_pT!=0) delete fit_pT;
  finput_pTrewgt->Close();
  return hArray;
}


void makeCombineSignalTemplateSyst_Signal_single(int folder, int erg_tev, int iProd, int scheme){
  if (scheme==0) return;
  gROOT->ProcessLine(".x tdrstyle.cc");
  gStyle->SetOptStat(0);
  char TREE_NAME[]="SelectedTree";
  TString INPUT_NAME = "HZZ4lTree_powheg15jhuGenV3-CTau";
  TString OUTPUT_NAME = "_templateSystematics_Sig_TxyBS_";
  if (iProd>0) OUTPUT_NAME.Append(Form("%s_", strHiggsProduction[iProd-1].Data()));
  OUTPUT_NAME = OUTPUT_NAME + cutNames[scheme] + ".root";
  OUTPUT_NAME.Prepend(user_folder[folder]);
  TString comstring;
  comstring.Form("%iTeV", erg_tev);
  TString erg_dir;
  erg_dir.Form("LHC_%iTeV", erg_tev);

  const int nSyst=3;
  int kSyst=nSyst;
  const int kSystNominal=1;
  const int iSystKD=1;
  char* cSyst[nSyst] ={
    "_TxyDown",
    "",
    "_TxyUp"
  };

  int EnergyIndex=1;
  if (erg_tev==7) EnergyIndex=0;
  float ZZMass_PeakCut[2]={ 105.6, 140.6 };

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
  //use ggH yields directly from SIP<4 case
  double targetYield = yield_signal_ggh[EnergyIndex][folder]*ratio_targetsignalyield/luminosity[EnergyIndex];

  TFile* finput_pTrewgt = new TFile(Form("./data/compareProductionModes_pT_%iTeV%s", erg_tev, ".root"), "read");
  TF1* fit_pT = 0;
  if (iProd>0){
    TString fitname = strHiggsProductionFit[iProd-1];
    fitname.Append("_pToverMZZ_ratio_fit");
    fit_pT = (TF1*)finput_pTrewgt->Get(fitname);
    if (iProd==1) targetYield  = yield_signal_vbfh[EnergyIndex][folder]*ratio_targetsignalyield/luminosity[EnergyIndex];
    else if (iProd==2) targetYield  = yield_signal_wh[EnergyIndex][folder]*ratio_targetsignalyield/luminosity[EnergyIndex];
    else if (iProd==3) targetYield  = yield_signal_zh[EnergyIndex][folder]*ratio_targetsignalyield/luminosity[EnergyIndex];
    else if (iProd==4) targetYield  = yield_signal_tth[EnergyIndex][folder]*ratio_targetsignalyield/luminosity[EnergyIndex];
  }
  cout << comstring << ' ' << user_folder[folder] << " yield: " << targetYield*luminosity[EnergyIndex] << endl;

  float templateWeight[nSyst] ={ 1 };
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
  float KD=-99, KD_up=-99, KD_dn=-99;
  float Txy, Dxy, Txy_true, Dxy_true, delTxy, delDxy;
  float Txy_BS, Dxy_BS, Txy_BS_true, Dxy_BS_true;
  float delTxy_BS=-99;
  float sigmaPV_xy, sigmaInt_xy;

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

  float D_bkg=-99;
  float D_bkg_ScaleUp=-99;
  float D_bkg_ScaleDown=-99;
  float D_bkg_ResUp=-99;
  float D_bkg_ResDown=-99;

  double ctau_gridsize = (ctau_limits[1]-ctau_limits[0])/(2*nctaus);
  int nbinsx = 30;
  double xlow = -2340, xhigh = 2340;
  double binwidthx = (xhigh-xlow)/nbinsx;
  xlow-=binwidthx*1.5; xhigh+=binwidthx*1.5; nbinsx+=3;
  int nbinsy = 50;
  double ylow=0, yhigh=1;

  TString INPUT_SYST_NAME = "LifetimeKD_SystFits_DataMC_";
  INPUT_SYST_NAME.Append("newSIP_");
  INPUT_SYST_NAME.Append(user_folder[folder]);
  INPUT_SYST_NAME.Append(Form("_%s_", "Txy_BeamSpot"));
  INPUT_SYST_NAME.Append("AllTeV");
  INPUT_SYST_NAME.Append(".root");
  TString cinput_syst = user_dir_hep + "Analysis/ShapeSystematics/";
  cinput_syst.Append("Txy_BeamSpot/");
  cinput_syst.Append(INPUT_SYST_NAME);
  TFile* finput_syst = new TFile(cinput_syst, "read");
  TH1F* txy_syst = (TH1F*)finput_syst->Get("hCommonSmearing");
  float smearingVal = txy_syst->GetBinContent(1);
  float smearingErr = txy_syst->GetBinError(1);

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

  TString cinput_common = user_dir_hep + "No_SIP/" + erg_dir + "/" + user_folder[folder] + "/";
  TString coutput_common = user_dir_hep + "Analysis/Templates/" + comstring + "/";
  TString coutput = coutput_common + OUTPUT_NAME;
  TFile* foutput = new TFile(coutput, "recreate");
  cout << "Starting to create " << coutput << endl;
  TH2F** hMC = makeCombineSignalTemplateMCSyst_Signal_single(folder, erg_tev, iProd, scheme);
  
  TChain* tree = 0;
  TString cinput = cinput_common + sample_FullSim[0] + ".root";;
  tree = new TChain(TREE_NAME);
  tree->Add(cinput);

  tree->SetBranchAddress("MC_weight", &MC_weight);
  tree->SetBranchAddress("ZZMass", &ZZMass);
  tree->SetBranchAddress("ZZPt", &ZZPt);
  tree->SetBranchAddress("ZZEta", &ZZEta);
  tree->SetBranchAddress("ZZPhi", &ZZPhi);
  tree->SetBranchAddress("p0plus_VAJHU", &p0plus_VAJHU);
  tree->SetBranchAddress("p0plus_VAMCFM", &p0plus_VAMCFM);
  tree->SetBranchAddress("bkg_VAMCFM", &bkg_VAMCFM);
  tree->SetBranchAddress("p0plus_m4l", &p0plus_m4l);
  tree->SetBranchAddress("bkg_m4l", &bkg_m4l);
  tree->SetBranchAddress("p0plus_m4l_ScaleUp", &p0plus_m4l_ScaleUp);
  tree->SetBranchAddress("bkg_m4l_ScaleUp", &bkg_m4l_ScaleUp);
  tree->SetBranchAddress("p0plus_m4l_ScaleDown", &p0plus_m4l_ScaleDown);
  tree->SetBranchAddress("bkg_m4l_ScaleDown", &bkg_m4l_ScaleDown);
  tree->SetBranchAddress("p0plus_m4l_ResUp", &p0plus_m4l_ResUp);
  tree->SetBranchAddress("bkg_m4l_ResUp", &bkg_m4l_ResUp);
  tree->SetBranchAddress("p0plus_m4l_ResDown", &p0plus_m4l_ResDown);
  tree->SetBranchAddress("bkg_m4l_ResDown", &bkg_m4l_ResDown);

  tree->SetBranchAddress("Lep1SIP", &Lep1SIP);
  tree->SetBranchAddress("Lep2SIP", &Lep2SIP);
  tree->SetBranchAddress("Lep3SIP", &Lep3SIP);
  tree->SetBranchAddress("Lep4SIP", &Lep4SIP);
  tree->SetBranchAddress("Lep1_Z1SIP", &Lep1_Z1SIP);
  tree->SetBranchAddress("Lep2_Z1SIP", &Lep2_Z1SIP);
  tree->SetBranchAddress("Lep3_Z1SIP", &Lep3_Z1SIP);
  tree->SetBranchAddress("Lep4_Z1SIP", &Lep4_Z1SIP);
  tree->SetBranchAddress("OfflinePrimaryVtx_x", &OfflinePrimaryVtx_x);
  tree->SetBranchAddress("OfflinePrimaryVtx_y", &OfflinePrimaryVtx_y);
  tree->SetBranchAddress("OfflinePrimaryVtx_z", &OfflinePrimaryVtx_z);
  tree->SetBranchAddress("OfflinePrimaryVtx_cov_xx", &OfflinePrimaryVtx_cov_xx);
  tree->SetBranchAddress("OfflinePrimaryVtx_cov_xy", &OfflinePrimaryVtx_cov_xy);
  tree->SetBranchAddress("OfflinePrimaryVtx_cov_xz", &OfflinePrimaryVtx_cov_xz);
  tree->SetBranchAddress("OfflinePrimaryVtx_cov_yy", &OfflinePrimaryVtx_cov_yy);
  tree->SetBranchAddress("OfflinePrimaryVtx_cov_yz", &OfflinePrimaryVtx_cov_yz);
  tree->SetBranchAddress("OfflinePrimaryVtx_cov_zz", &OfflinePrimaryVtx_cov_zz);
  tree->SetBranchAddress("OfflinePrimaryVtx_ndof", &OfflinePrimaryVtx_ndof);
  tree->SetBranchAddress("KalmanCandVtx_x", &KalmanCandVtx_x);
  tree->SetBranchAddress("KalmanCandVtx_y", &KalmanCandVtx_y);
  tree->SetBranchAddress("KalmanCandVtx_z", &KalmanCandVtx_z);
  tree->SetBranchAddress("KalmanCandVtx_chi2", &KalmanCandVtx_chi2);
  tree->SetBranchAddress("KalmanCandVtx_cov_xx", &KalmanCandVtx_cov_xx);
  tree->SetBranchAddress("KalmanCandVtx_cov_xy", &KalmanCandVtx_cov_xy);
  tree->SetBranchAddress("KalmanCandVtx_cov_xz", &KalmanCandVtx_cov_xz);
  tree->SetBranchAddress("KalmanCandVtx_cov_yy", &KalmanCandVtx_cov_yy);
  tree->SetBranchAddress("KalmanCandVtx_cov_yz", &KalmanCandVtx_cov_yz);
  tree->SetBranchAddress("KalmanCandVtx_cov_zz", &KalmanCandVtx_cov_zz);
  tree->SetBranchAddress("Z1CandVtx_x", &Z1CandVtx_x);
  tree->SetBranchAddress("Z1CandVtx_y", &Z1CandVtx_y);
  tree->SetBranchAddress("Z1CandVtx_z", &Z1CandVtx_z);

  if (tree->GetBranchStatus("GenPrimaryVtx_x")){
    tree->SetBranchAddress("GenPrimaryVtx_x", &GenPrimaryVtx_x);
    tree->SetBranchAddress("GenPrimaryVtx_y", &GenPrimaryVtx_y);
    tree->SetBranchAddress("GenPrimaryVtx_z", &GenPrimaryVtx_z);
    tree->SetBranchAddress("GenIntVtx_x", &GenIntVtx_x);
    tree->SetBranchAddress("GenIntVtx_y", &GenIntVtx_y);
    tree->SetBranchAddress("GenIntVtx_z", &GenIntVtx_z);
    tree->SetBranchAddress("GenHPt", &GenHPt);
    tree->SetBranchAddress("GenHPhi", &GenHPhi);
    tree->SetBranchAddress("GenHMass", &GenHMass);
  }

  TH1F* hProj_grandcontainer[2*nctaus+1][nSyst];
  for (int genctau = 0; genctau <= 2*nctaus; genctau++){
    double target_ctau = ctau_limits[0] + ctau_gridsize*genctau;
    double Txy_BS_systWeight = 1.;

    TH1F* hProj_container[nSyst];
    int genctau_2 = genctau/2;

    TH2F* hfill[nSyst];
    for (int ss = 0; ss < nSyst; ss++){
      hfill[ss] = new TH2F(Form("H_2D_CTau%.0f%s", target_ctau, cSyst[ss]), "", nbinsx, xlow, xhigh, nbinsy, ylow, yhigh);
      hfill[ss]->SetTitle(Form("%.0f #mum %s at %i TeV", target_ctau, user_folder[folder], erg_tev));
      hfill[ss]->GetXaxis()->SetTitle("T^{Beam}_{xy} (#mum)");
      hfill[ss]->GetYaxis()->SetTitle("#it{D}_{bkg}");
      hfill[ss]->GetZaxis()->SetTitle("Events / bin");
    }
    TH2F* hfill_conditional;
    TH1F* hfill_projX;

    for (int ev = 0; ev < tree->GetEntries(); ev++){
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

      tree->GetEntry(ev);
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

      float smd = D_bkg;

      float KDSystwgt_up = get_KD_Syst(smearingVal, smearingErr, 1);
      float KDSystwgt_dn = get_KD_Syst(smearingVal, smearingErr, -1);
      float KDSystwgt_nominal = get_KD_Syst(smearingVal, smearingErr, 0);

      KD = compute_KDAlternatives(Txy_BS, Txy_BS_true, KDSystwgt_nominal);
      //      KD = tanh(KD / KD_scale); // Scaled to due to tanh function
      //      if (KD >= 1.) KD = 1. - 1.0e-4;
      //      if (KD <= -1.) KD = -1. + 1.0e-4;

      KD_up = compute_KDAlternatives(Txy_BS, Txy_BS_true, KDSystwgt_up);
      //      KD_up = tanh(KD_up / KD_scale); // Scaled to due to tanh function
      //      if (KD_up >= 1.) KD_up = 1. - 1.0e-10;
      //      if (KD_up <= -1.) KD_up = -1. + 1.0e-10;

      KD_dn = compute_KDAlternatives(Txy_BS, Txy_BS_true, KDSystwgt_dn);
      //      KD_dn = tanh(KD_dn / KD_scale); // Scaled to due to tanh function
      //      if (KD_dn >= 1.) KD_dn = 1. - 1.0e-10;
      //      if (KD_dn <= -1.) KD_dn = -1. + 1.0e-10;

      for (int ss = 0; ss < nSyst; ss++){
        float KD_fill;
        if (ss == (kSystNominal + iSystKD)) KD_fill=KD_up;
        else if (ss == (kSystNominal - iSystKD)) KD_fill=KD_dn;
        else KD_fill=KD;

        int first_binX = hfill[ss]->GetXaxis()->FindBin(KD_fill);
        int biny = hfill[ss]->GetYaxis()->FindBin(smd);
        if (target_ctau>0){
          for (int binx=first_binX; binx<=hfill[ss]->GetNbinsX()+1; binx++){
            float xmin, xmax;
            templateWeight[ss] = MC_weight;
            if (fit_pT!=0 && iProd>0) templateWeight[ss] *= fit_pT->Eval(GenHPt/GenHMass);
            if (binx>0) xmin = max((double)KD_fill, hfill[ss]->GetXaxis()->GetBinLowEdge(binx));
            else xmin = KD_fill;
            xmin -= KD_fill;

            if (binx<hfill[ss]->GetNbinsX()+1){
              xmax = hfill[ss]->GetXaxis()->GetBinUpEdge(binx);
              xmax -= KD_fill;
              templateWeight[ss] *= (exp(-xmin/target_ctau) - exp(-xmax/target_ctau));
            }
            else{
              templateWeight[ss] *= exp(-xmin/target_ctau);
            }
            if (smd<0) templateWeight[ss] = 0;

            hfill[ss]->AddBinContent(hfill[ss]->GetBin(binx, biny), templateWeight[ss]);
          }
        }
        else{
          templateWeight[ss] = MC_weight;
          if (fit_pT!=0 && iProd>0) templateWeight[ss] *= fit_pT->Eval(GenHPt/GenHMass);
          hfill[ss]->Fill(KD_fill, smd, templateWeight[ss]);
        }
      }
    }
    double integral_rescale[nSyst] ={ 1 };
    for (int ss = 0; ss < nSyst; ss++){
      hfill[ss]->SetOption("colz");
      double integral_rewgt = hfill[ss]->Integral(0, nbinsx+1, 0, nbinsy+1);
      double integral_target = targetYield;
      integral_rescale[ss] = integral_target / integral_rewgt;
      cout << "Template for cTau = " << target_ctau << " of " << user_folder[folder] << " @ " << comstring << " is rescaled with " << integral_rescale[ss] << " for systematic " << ss << "..." << endl;
      hfill[ss]->Scale(integral_rescale[ss]);
      for (int biny=0; biny<=hfill[ss]->GetNbinsY()+1; biny++){
        double underflow = hfill[ss]->GetBinContent(0, biny);
        double overflow = hfill[ss]->GetBinContent(hfill[ss]->GetNbinsX()+1, biny);
        hfill[ss]->AddBinContent(hfill[ss]->GetBin(1, biny), underflow);
        hfill[ss]->AddBinContent(hfill[ss]->GetBin(hfill[ss]->GetNbinsX(), biny), overflow);
        hfill[ss]->SetBinContent(0, biny, 0);
        hfill[ss]->SetBinContent(hfill[ss]->GetNbinsX()+1, biny, 0);
      }
    }
    for (int ss = 0; ss < nSyst; ss++){
      hfill_conditional = (TH2F*)hfill[ss]->Clone(Form("%s_Conditional", hfill[ss]->GetName()));
      hfill_conditional->SetTitle("");
      for (int binx = 1; binx <= hfill_conditional->GetNbinsX(); binx++){
        double integralX = hfill_conditional->Integral(binx, binx, 0, hfill_conditional->GetNbinsY() + 1);
        if (integralX != 0){
          for (int biny = 1; biny <= hfill_conditional->GetNbinsY(); biny++) hfill_conditional->SetBinContent(binx, biny, (hfill_conditional->GetBinContent(binx, biny) / integralX));
        }
      }

      hfill_projX = (TH1F*)hfill[ss]->ProjectionX();
      hfill_projX->SetTitle("");
      hfill_projX->GetYaxis()->SetTitle("Events / bin");

      foutput->WriteTObject(hfill[ss]);
      foutput->WriteTObject(hfill_conditional);
      foutput->WriteTObject(hfill_projX);
      hProj_container[ss] = hfill_projX;
      hProj_grandcontainer[genctau][ss] = hfill_projX;
      //      delete hfill_projX;
      delete hfill_conditional;
      delete hfill[ss];
    }

    if (hMC==0) cout << "No hMC!" << endl;
    if ((genctau % 2)==0){
      TH1F* hMC_proj = (TH1F*)hMC[genctau_2]->ProjectionX();
      foutput->WriteTObject(hMC[genctau_2]);
      foutput->WriteTObject(hMC_proj);
      hMC_proj->SetStats(0);
      hMC[genctau_2]->SetStats(0);

      hMC_proj->Scale(luminosity[EnergyIndex]);
      hProj_container[kSystNominal + iSystKD]->Scale(luminosity[EnergyIndex]);
      hProj_container[kSystNominal - iSystKD]->Scale(luminosity[EnergyIndex]);
      hProj_container[kSystNominal]->Scale(luminosity[EnergyIndex]);

      TH1F* hFit=0;
      TH1F* hFit_up=0;
      TH1F* hFit_dn=0;
      double n1val=0, n2val=0;
      if (genctau>0){
        RooRealVar* frac1 = new RooRealVar(Form("frac1_CTau%.0f", target_ctau), "", 1, 0, 1);
        RooRealVar* frac2 = new RooRealVar(Form("frac2_CTau%.0f", target_ctau), "", 0, 0, 1);
        RooRealVar* xvar = new RooRealVar("KD", "", 0., xlow, xhigh);
        xvar->setBins(nbinsx);
        RooDataHist* data = new RooDataHist(Form("data_CTau%.0f", target_ctau), "", *xvar, hMC_proj);
        RooArgList fitArgs(*xvar);
        RooDataHist* h0 = new RooDataHist("h0", "", fitArgs, hProj_grandcontainer[0][kSystNominal]);
        RooDataHist* h1 = new RooDataHist("h1", "", fitArgs, hProj_container[kSystNominal]);
        RooDataHist* h2 = new RooDataHist("h2", "", fitArgs, hProj_grandcontainer[genctau_2][kSystNominal]);
        RooHistFunc* f0 = new RooHistFunc("f0", "", RooArgSet(*xvar), *h0);
        RooHistFunc* f1 = new RooHistFunc("f1", "", RooArgSet(*xvar), *h1);
        RooHistFunc* f2 = new RooHistFunc("f2", "", RooArgSet(*xvar), *h2);
        RooFormulaVar* norm0 = new RooFormulaVar("norm0", "(@0+@1)>1? 0. : (1.-@0-@1)", RooArgList(*frac1, *frac2));
        RooFormulaVar* norm1 = new RooFormulaVar("norm1", "(@0+@1)>1? 0. : @0", RooArgList(*frac1, *frac2));
        RooFormulaVar* norm2 = new RooFormulaVar("norm2", "(@0+@1)>1? 0. : @1", RooArgList(*frac1, *frac2));
        RooRealFlooredSumPdf* pdf = new RooRealFlooredSumPdf("pdf", "", RooArgList(*f0, *f1, *f2), RooArgList(*norm0, *norm1, *norm2));
//        pdf->fitTo(*data, SumW2Error(kFALSE));

        frac1->setVal(1);
        frac2->setVal(0);

        foutput->cd();
        hFit = (TH1F*)pdf->createHistogram(Form("fit_CTau%.0f", target_ctau), *xvar);
        hFit->SetTitle("");
        //      hFit->SetLineStyle(7);
        hFit->SetLineColor(kViolet);
        hFit->SetLineWidth(2);
        n1val = norm1->getVal();
        n2val = norm2->getVal();
        cout << "Normalization 1: " << n1val << " and Normalization 2: " << n2val << endl;

        delete pdf;
        delete norm2;
        delete norm1;
        delete norm0;
        delete f2;
        delete f1;
        delete f0;
        delete h2;
        delete h1;
        delete h0;
        delete data;
        delete xvar;
        delete frac2;
        delete frac1;
      }
      if (hFit!=0){
        hFit->Scale(hMC_proj->Integral()/hFit->Integral());

        hFit_up = (TH1F*)hProj_container[kSystNominal + iSystKD]->Clone(Form("%s_up", hFit->GetName()));
        hFit_up->Scale(n1val);
        hFit_up->Add(hProj_grandcontainer[genctau_2][kSystNominal + iSystKD], n2val);
        hFit_up->Add(hProj_grandcontainer[0][kSystNominal + iSystKD], 1.-n2val-n1val);

        hFit_dn = (TH1F*)hProj_container[kSystNominal - iSystKD]->Clone(Form("%s_dn", hFit->GetName()));
        hFit_dn->Scale(n1val);
        hFit_dn->Add(hProj_grandcontainer[genctau_2][kSystNominal - iSystKD], n2val);
        hFit_dn->Add(hProj_grandcontainer[0][kSystNominal - iSystKD], 1.-n2val-n1val);

        cout << "Scaling up by " << hFit->Integral()/hFit_up->Integral() << " and down by " << hFit->Integral()/hFit_dn->Integral() << endl;
        hFit_up->Scale(hFit->Integral()/hFit_up->Integral());
        hFit_dn->Scale(hFit->Integral()/hFit_dn->Integral());
      }
      else cout << "hfit na" << endl;

      double maxplot = max(hMC_proj->GetMaximum(), hProj_container[kSystNominal - iSystKD]->GetMaximum());
      if (hFit_dn!=0) maxplot = max(maxplot, hFit_dn->GetMaximum());
      hMC_proj->GetYaxis()->SetTitle("Events / bin");
      hMC_proj->GetYaxis()->SetRangeUser(0, maxplot*1.7);

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

      hMC_proj->SetLineColor(kBlack);
      hMC_proj->SetLineWidth(2);
      hProj_container[kSystNominal - iSystKD]->SetLineColor(kRed);
      hProj_container[kSystNominal - iSystKD]->SetLineWidth(2);
      hProj_container[kSystNominal - iSystKD]->SetLineStyle(7);
      hProj_container[kSystNominal + iSystKD]->SetLineColor(kGreen+2);
      hProj_container[kSystNominal + iSystKD]->SetLineWidth(2);
      hProj_container[kSystNominal + iSystKD]->SetLineStyle(7);
      hProj_container[kSystNominal]->SetLineColor(kBlue);
      hProj_container[kSystNominal]->SetLineWidth(2);

      if (hFit!=0){
        hFit_up->SetLineColor(hProj_container[kSystNominal + iSystKD]->GetLineColor());
        hFit_up->SetLineWidth(2);
        hFit_up->SetLineStyle(7);

        hFit_dn->SetLineColor(hProj_container[kSystNominal - iSystKD]->GetLineColor());
        hFit_dn->SetLineWidth(2);
        hFit_dn->SetLineStyle(7);
      }

      TString plotDir = coutput_common;
      plotDir.Append("ValidationPlots/");
      if (iProd>0) plotDir.Append(Form("%s/", strHiggsProduction[iProd-1].Data()));
      TString mkdirCommand = "mkdir -p ";
      mkdirCommand.Append(plotDir);
      gSystem->Exec(mkdirCommand);

      foutput->cd();
      TString appendName;
      appendName = Form("CTau%.0f_%s_%iTeV", target_ctau, user_folder[folder], erg_tev);
      TString canvasname = "cCompareLifetime_TxyBS_Syst_";
      if (iProd>0) canvasname.Append(Form("%s_", strHiggsProduction[iProd-1].Data()));
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
      ll = new TLegend(0.20, 0.60, 0.45, 0.90);
      ll->SetBorderSize(0);
      ll->SetTextFont(42);
      ll->SetTextSize(0.03);
      ll->SetLineColor(1);
      ll->SetLineStyle(1);
      ll->SetLineWidth(1);
      ll->SetFillColor(0);
      ll->SetFillStyle(0);

      ll->AddEntry(hMC_proj, "Higgs MC prediction", "l");
      if (genctau!=0) hMC_proj->GetXaxis()->SetRangeUser(-800, 2400);
      else hMC_proj->GetXaxis()->SetRangeUser(-1000, 1000);
      hMC_proj->Draw("hist");
      if (genctau==0){
        ll->AddEntry(hProj_container[kSystNominal + iSystKD], "Resolution up", "l");
        ll->AddEntry(hProj_container[kSystNominal - iSystKD], "Resolution down", "l");
        hProj_container[kSystNominal + iSystKD]->Draw("histsame");
        hProj_container[kSystNominal - iSystKD]->Draw("histsame");
      }
      else if (hFit!=0){
//        ll->AddEntry(hProj_container[kSystNominal], "Convoluted c#tau=0 #mum (nominal)", "l");
        ll->AddEntry(hFit, "Convolution (nominal)", "l");
        ll->AddEntry(hFit_up, "Convolution (res. up)", "l");
        ll->AddEntry(hFit_dn, "Convolution (res. down)", "l");

        hProj_container[kSystNominal]->Draw("histsame");
        hFit->Draw("histsame");
        hFit_up->Draw("histsame");
        hFit_dn->Draw("histsame");
      }
      ll->Draw("same");
      pt->Draw();


      TPaveText *pt20 = new TPaveText(0.75, 0.86, 0.90, 0.90, "brNDC");
      pt20->SetBorderSize(0);
      pt20->SetTextAlign(12);
      pt20->SetTextSize(0.03);
      pt20->SetFillStyle(0);
      pt20->SetTextFont(42);
      TText* text20 = pt20->AddText(0.01, 0.01, Form("c#tau_{H}=%.0f #mum", target_ctau));
      pt20->Draw();

      cc->RedrawAxis();
      cc->Modified();
      cc->Update();
      foutput->WriteTObject(cc);
      canvasname.Prepend(plotDir);
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
      delete ll;
      delete pt;
      cc->Close();

      delete hMC_proj;
      if (hFit!=0){
        foutput->WriteTObject(hFit);
        foutput->WriteTObject(hFit_up);
        foutput->WriteTObject(hFit_dn);
        delete hFit;
        delete hFit_up;
        delete hFit_dn;
      }
    }
  }
  for (int genctau = 0; genctau <= 2*nctaus; genctau++){
    for (int ss = 0; ss < nSyst; ss++) delete hProj_grandcontainer[genctau][ss];
  }
  if (tree!=0) delete tree;

  cout << "Deleting external MC" << endl;
  for (int tt=0; tt<=nctaus; tt++) delete hMC[tt];
  delete[] hMC;
  foutput->Close();
  cout << "Closed " << coutput << endl;
  delete tBeam;
  delete txy_syst;
  finput_syst->Close();
  if (fit_pT!=0) delete fit_pT;
  finput_pTrewgt->Close();
}


void makeCombineSignalTemplateswithTrees_Bkg_single(int folder, int erg_tev, int Systematics, int scheme){
/*	if (Systematics !=0){
		cout << "Txy_BS up and down are not yet supported" << endl;
		return;
	}
*/
	char TREE_NAME[]="SelectedTree";
	TString OUTPUT_NAME = "_templates_background_TxyBS_" + cutNames[scheme] + "_";
	OUTPUT_NAME.Prepend(user_folder[folder]);
	TString comstring;
	comstring.Form("%iTeV",erg_tev);
	TString erg_dir;
	erg_dir.Form("LHC_%iTeV",erg_tev);
	int EnergyIndex=0;if(erg_tev==8) EnergyIndex=1;

	if (Systematics == 0) OUTPUT_NAME += "Nominal.root";
	if (Systematics == 1) OUTPUT_NAME += "TxyUp.root";
	if (Systematics == -1) OUTPUT_NAME += "TxyDown.root";

  TString INPUT_SYST_NAME = "LifetimeKD_SystFits_DataMC_";
  INPUT_SYST_NAME.Append("newSIP_");
  INPUT_SYST_NAME.Append(user_folder[folder]);
  INPUT_SYST_NAME.Append(Form("_%s_", "Txy_BeamSpot"));
  INPUT_SYST_NAME.Append("AllTeV");
  INPUT_SYST_NAME.Append(".root");
  TString cinput_syst = user_dir_hep + "Analysis/ShapeSystematics/";
  cinput_syst.Append("Txy_BeamSpot/");
  cinput_syst.Append(INPUT_SYST_NAME);
  TFile* finput_syst = new TFile(cinput_syst, "read");
  TH1F* txy_syst = (TH1F*)finput_syst->Get("hCommonSmearing");
  float smearingVal = txy_syst->GetBinContent(1);
  float smearingErr = txy_syst->GetBinError(1);
  
  float ZZMass_PeakCut[2]={105.6,140.6};
	const int ntrees=3;
	double initialYield[ntrees]={
		yield_signal_qqzz[EnergyIndex][folder] / luminosity[EnergyIndex],
    yield_signal_ggzz[EnergyIndex][folder] / luminosity[EnergyIndex],
    yield_signal_zx[EnergyIndex][folder]
	};
  if (Systematics!=0) initialYield[2] = yield_offshell_zx[EnergyIndex][folder];
  double zx_offshell_sum = 0;
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
    if (iProcess!=6) cout << comstring << ' ' << user_folder[folder] << ' ' << processName[iProcess] <<  " yield: " << targetYield[iProcess-4]*luminosity[EnergyIndex] << endl;
    else cout << comstring << ' ' << user_folder[folder] << ' ' << processName[iProcess] <<  " yield: " << targetYield[iProcess-4] << endl;
  }
	double scale_yield[ntrees] = { 1, 1, 1 };

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

	float templateWeight=1;
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
  float KD=-99;
  float Txy, Dxy, Txy_true, Dxy_true, delTxy, delDxy;
  float Txy_BS, Dxy_BS, Txy_BS_true, Dxy_BS_true;
  float delTxy_BS=-99;
  float sigmaPV_xy, sigmaInt_xy;

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

	float D_bkg=-99;

	TString cinput_common = user_dir_hep + "No_SIP/" + erg_dir + "/";
	TString coutput_common = user_dir_hep + "Analysis/Templates/" + comstring + "/";
	TString coutput = coutput_common + OUTPUT_NAME;
	TFile* foutput = new TFile(coutput,"recreate");

  // BeamSpot info
  TChain* tBeam[2];
  TString strBeamSpot[2][2]={
    { "data_GR_R_44_V15C", "MC_START44_V13" },
    { "data_FT_53_V21_AN4", "MC_START53_V23" }
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

	TString cinput_qqzz_common = cinput_common + user_folder[folder] + "/";
	TString cinput_ggzz_common = cinput_qqzz_common;
	TString cinput_zx_common = cinput_common;

	TChain* tqqzz = new TChain(TREE_NAME);
	for (int smp = kGGSamples; smp < kQQBZZSamples_Dedicated; smp++){
		TString cinput_qqzz = cinput_qqzz_common + sample_FullSim[smp] + "_Reprocessed.root";
		tqqzz->Add(cinput_qqzz);
	}
	TChain* tggzz = new TChain(TREE_NAME);
	for (int smp = kGGHSamples; smp < kQQBZZSamples_Dedicated; smp++){
		TString cinput_ggzz = cinput_qqzz_common + sample_FullSim[smp] + "_Reprocessed.root";
		tggzz->Add(cinput_ggzz);
	}
	TChain* tzx = new TChain(TREE_NAME);
  TString cinput_zx = user_dir_hep;
  cinput_zx.Append("No_SIP/");
  cinput_zx.Append(erg_dir);
  cinput_zx = cinput_zx + "/CR/" + sample_FullSim[kAllSamples - 1] + ".root";
  tzx->Add(cinput_zx);

	int nbinsx = 50;
	double xlow = -1, xhigh = 1;
  double binwidthx = (xhigh-xlow)/nbinsx;
	int nbinsy = 50;
	double ylow=0,yhigh=1;
	TH2F* hqqzz_full = new TH2F("h_template_qqZZ",Form("%i TeV %s q#bar{q}#rightarrow4l Background",erg_tev,user_folder[folder]),nbinsx,xlow,xhigh,nbinsy,ylow,yhigh);
  hqqzz_full->GetXaxis()->SetTitle(Form("tanh(T_{xy} / %.0f mm)", KD_scale));
	hqqzz_full->GetYaxis()->SetTitle("#it{D}_{bkg}");
	TH2F* hggzz_full = (TH2F*) hqqzz_full->Clone();
	TH2F* hzx_full = (TH2F*) hqqzz_full->Clone();
	TH2F* hmc_full = (TH2F*) hqqzz_full->Clone();
	hggzz_full->SetNameTitle("h_template_ggZZ",Form("%i TeV %s gg Background",erg_tev,user_folder[folder]));
	hzx_full->SetNameTitle("h_template_ZX",Form("%i TeV %s Z+X Background",erg_tev,user_folder[folder]));
	hmc_full->SetNameTitle("h_template_totalbkg",Form("%i TeV %s Total Background",erg_tev,user_folder[folder]));
	TH2F* hfull[ntrees+1] = { hqqzz_full, hggzz_full, hzx_full, hmc_full };
	char* ctreefull[ntrees] = { "qqZZ", "ggZZ", "ZX" };
	for (int tt = 0; tt < ntrees; tt++) hfull[tt]->Sumw2();
	TH1F* hproj[ntrees+1];

	TChain* tc[ntrees] = { tqqzz, tggzz, tzx };
	TTree* mytree[ntrees];
	for (int tt = 0; tt < ntrees; tt++){
		cout << "Tree " << tt << endl;
    if (tt != 2){
			tc[tt]->SetBranchAddress("MC_weight_noxsec", &MC_weight_noxsec);
      tc[tt]->SetBranchAddress("MC_weight_Kfactor", &MC_weight_Kfactor);
      tc[tt]->SetBranchAddress("MC_weight_QQBGGProper", MC_weight_QQBGGProper);
      if (tt == 0){
        tc[tt]->SetBranchAddress("MC_weight_QQZZEWK", &MC_weight_QQZZEWK);
      }

      tc[tt]->SetBranchAddress("D_bkg", &D_bkg);
    }
    else{
      tc[tt]->SetBranchAddress("RunNumber", &RunNumber);
      tc[tt]->SetBranchAddress("LumiNumber", &LumiNumber);

      tc[tt]->SetBranchAddress("ZXfake_weight_OS", ZXfake_weight_OS);
      tc[tt]->SetBranchAddress("ZXfake_weight_SS", &ZXfake_weight_SS);

      tc[tt]->SetBranchAddress("CRflag", &CRflag);

      tc[tt]->SetBranchAddress("Lep1combRelIsoPF", &Lep1combRelIsoPF);
      tc[tt]->SetBranchAddress("Lep2combRelIsoPF", &Lep2combRelIsoPF);
      tc[tt]->SetBranchAddress("Lep3combRelIsoPF", &Lep3combRelIsoPF);
      tc[tt]->SetBranchAddress("Lep4combRelIsoPF", &Lep4combRelIsoPF);

      tc[tt]->SetBranchAddress("Lep3isID", &Lep3isID);
      tc[tt]->SetBranchAddress("Lep4isID", &Lep4isID);

      tc[tt]->SetBranchAddress("Z1ids", &Z1ids);
      tc[tt]->SetBranchAddress("Lep3ID", &Lep3ID);
      tc[tt]->SetBranchAddress("Lep4ID", &Lep4ID);

      tc[tt]->SetBranchAddress("p0plus_VAJHU", &p0plus_VAJHU);
      tc[tt]->SetBranchAddress("p0plus_VAMCFM", &p0plus_VAMCFM);
      tc[tt]->SetBranchAddress("bkg_VAMCFM", &bkg_VAMCFM);
      tc[tt]->SetBranchAddress("p0plus_m4l", &p0plus_m4l);
      tc[tt]->SetBranchAddress("bkg_m4l", &bkg_m4l);
      tc[tt]->SetBranchAddress("p0plus_m4l_ScaleUp", &p0plus_m4l_ScaleUp);
      tc[tt]->SetBranchAddress("bkg_m4l_ScaleUp", &bkg_m4l_ScaleUp);
      tc[tt]->SetBranchAddress("p0plus_m4l_ScaleDown", &p0plus_m4l_ScaleDown);
      tc[tt]->SetBranchAddress("bkg_m4l_ScaleDown", &bkg_m4l_ScaleDown);
      tc[tt]->SetBranchAddress("p0plus_m4l_ResUp", &p0plus_m4l_ResUp);
      tc[tt]->SetBranchAddress("bkg_m4l_ResUp", &bkg_m4l_ResUp);
      tc[tt]->SetBranchAddress("p0plus_m4l_ResDown", &p0plus_m4l_ResDown);
      tc[tt]->SetBranchAddress("bkg_m4l_ResDown", &bkg_m4l_ResDown);
    }

		tc[tt]->SetBranchAddress("ZZMass", &ZZMass);
		tc[tt]->SetBranchAddress("ZZPt", &ZZPt);
		tc[tt]->SetBranchAddress("ZZEta", &ZZEta);
		tc[tt]->SetBranchAddress("ZZPhi", &ZZPhi);

		tc[tt]->SetBranchAddress("Lep1SIP", &Lep1SIP);
		tc[tt]->SetBranchAddress("Lep2SIP", &Lep2SIP);
		tc[tt]->SetBranchAddress("Lep3SIP", &Lep3SIP);
		tc[tt]->SetBranchAddress("Lep4SIP", &Lep4SIP);
		tc[tt]->SetBranchAddress("Lep1_Z1SIP", &Lep1_Z1SIP);
		tc[tt]->SetBranchAddress("Lep2_Z1SIP", &Lep2_Z1SIP);
		tc[tt]->SetBranchAddress("Lep3_Z1SIP", &Lep3_Z1SIP);
		tc[tt]->SetBranchAddress("Lep4_Z1SIP", &Lep4_Z1SIP);
		tc[tt]->SetBranchAddress("OfflinePrimaryVtx_x", &OfflinePrimaryVtx_x);
		tc[tt]->SetBranchAddress("OfflinePrimaryVtx_y", &OfflinePrimaryVtx_y);
		tc[tt]->SetBranchAddress("OfflinePrimaryVtx_z", &OfflinePrimaryVtx_z);
		tc[tt]->SetBranchAddress("OfflinePrimaryVtx_cov_xx", &OfflinePrimaryVtx_cov_xx);
		tc[tt]->SetBranchAddress("OfflinePrimaryVtx_cov_xy", &OfflinePrimaryVtx_cov_xy);
		tc[tt]->SetBranchAddress("OfflinePrimaryVtx_cov_xz", &OfflinePrimaryVtx_cov_xz);
		tc[tt]->SetBranchAddress("OfflinePrimaryVtx_cov_yy", &OfflinePrimaryVtx_cov_yy);
		tc[tt]->SetBranchAddress("OfflinePrimaryVtx_cov_yz", &OfflinePrimaryVtx_cov_yz);
		tc[tt]->SetBranchAddress("OfflinePrimaryVtx_cov_zz", &OfflinePrimaryVtx_cov_zz);
		tc[tt]->SetBranchAddress("OfflinePrimaryVtx_ndof", &OfflinePrimaryVtx_ndof);
		tc[tt]->SetBranchAddress("KalmanCandVtx_x", &KalmanCandVtx_x);
		tc[tt]->SetBranchAddress("KalmanCandVtx_y", &KalmanCandVtx_y);
		tc[tt]->SetBranchAddress("KalmanCandVtx_z", &KalmanCandVtx_z);
		tc[tt]->SetBranchAddress("KalmanCandVtx_chi2", &KalmanCandVtx_chi2);
		tc[tt]->SetBranchAddress("KalmanCandVtx_cov_xx", &KalmanCandVtx_cov_xx);
		tc[tt]->SetBranchAddress("KalmanCandVtx_cov_xy", &KalmanCandVtx_cov_xy);
		tc[tt]->SetBranchAddress("KalmanCandVtx_cov_xz", &KalmanCandVtx_cov_xz);
		tc[tt]->SetBranchAddress("KalmanCandVtx_cov_yy", &KalmanCandVtx_cov_yy);
		tc[tt]->SetBranchAddress("KalmanCandVtx_cov_yz", &KalmanCandVtx_cov_yz);
		tc[tt]->SetBranchAddress("KalmanCandVtx_cov_zz", &KalmanCandVtx_cov_zz);
		tc[tt]->SetBranchAddress("Z1CandVtx_x", &Z1CandVtx_x);
		tc[tt]->SetBranchAddress("Z1CandVtx_y", &Z1CandVtx_y);
		tc[tt]->SetBranchAddress("Z1CandVtx_z", &Z1CandVtx_z);

    if (tc[tt]->GetBranchStatus("GenPrimaryVtx_x")){
      tc[tt]->SetBranchAddress("GenPrimaryVtx_x", &GenPrimaryVtx_x);
      tc[tt]->SetBranchAddress("GenPrimaryVtx_y", &GenPrimaryVtx_y);
      tc[tt]->SetBranchAddress("GenPrimaryVtx_z", &GenPrimaryVtx_z);
      tc[tt]->SetBranchAddress("GenIntVtx_x", &GenIntVtx_x);
      tc[tt]->SetBranchAddress("GenIntVtx_y", &GenIntVtx_y);
      tc[tt]->SetBranchAddress("GenIntVtx_z", &GenIntVtx_z);
      tc[tt]->SetBranchAddress("GenHPt", &GenHPt);
      tc[tt]->SetBranchAddress("GenHPhi", &GenHPhi);
      tc[tt]->SetBranchAddress("GenHMass", &GenHMass);
    }

		mytree[tt] = new TTree(Form("T_2D_%s", ctreefull[tt]), "");
		mytree[tt]->Branch("MC_weight",&templateWeight);
		mytree[tt]->Branch("D_bkg",&D_bkg);
    mytree[tt]->Branch("KD", &KD);
    mytree[tt]->Branch("ZZMass", &ZZMass);
//    if (tt==2) mytree[tt]->Branch("CRflag", &CRflag);
  }
	for (int tt = 0; tt < ntrees; tt++){
    int counter_edges[2]={ 0 };
    int smeared_edges[2]={ 1, 1 };
		for (int ev = 0; ev < tc[tt]->GetEntries(); ev++){
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

			tc[tt]->GetEntry(ev);

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
      if (tt==0) wgt *= MC_weight_QQBGGProper[1]*MC_weight_QQZZEWK;
      else if (tt==1) wgt *= MC_weight_Kfactor*MC_weight_QQBGGProper[0];

      if ((tt==2 && Systematics!=0) && (ZZMass >= 220 && ZZMass < 1600)) zx_offshell_sum += wgt;

      if (!(tt==2 && Systematics!=0) && !(ZZMass >= ZZMass_PeakCut[0] && ZZMass < ZZMass_PeakCut[1])) continue;
      else if ((tt==2 && Systematics!=0) && !(ZZMass >= ZZMass_PeakCut[1] && ZZMass < 170)) continue;

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

      float KDSystwgt = 1;
      if (tt!=2) KDSystwgt = get_KD_Syst(smearingVal, smearingErr, Systematics);
      KD = compute_KDAlternatives(Txy_BS, Txy_BS_true, KDSystwgt);
      KD = tanh(KD / KD_scale); // Scaled to due to tanh function
      if (KD!=KD) continue;
      if (KD >= 1.) {
        KD = 1. - 1e-4;
        counter_edges[1]++;
      }
      if (KD <= -1.) {
        KD = -1. + 1e-4;
        counter_edges[0]++;
      }
			if (tt==2) D_bkg = p0plus_VAJHU*p0plus_m4l / (p0plus_VAJHU*p0plus_m4l + bkg_VAMCFM*bkg_m4l);
			if (D_bkg!=D_bkg) continue;

//      if (tt==2) cout << ZZMass << '\t' << Txy_BS << '\t' << MC_weight_noxsec << '\t' << CRflag << endl;

			templateWeight = wgt;
			hfull[tt]->Fill(KD,D_bkg,templateWeight);
		}
    if (tt==2 && Systematics!=0) scale_yield[tt] = targetYield[tt]/zx_offshell_sum;
		else scale_yield[tt] = targetYield[tt]/hfull[tt]->Integral(0, hfull[tt]->GetNbinsX() + 1, 0, hfull[tt]->GetNbinsY() + 1);
    cout << "Scaling " << hfull[tt]->GetName() << " by " << scale_yield[tt] << endl;
		hfull[tt]->Scale(scale_yield[tt]);

    cout << "Counter edges: " << counter_edges[0] << "\t" << counter_edges[1] << endl;

    double sum_filled=0;

		for (int ev = 0; ev < tc[tt]->GetEntries(); ev++){
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

      tc[tt]->GetEntry(ev);

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

      if (!(tt==2 && Systematics!=0) && !(ZZMass >= ZZMass_PeakCut[0] && ZZMass < ZZMass_PeakCut[1])) continue;
      else if ((tt==2 && Systematics!=0) && !(ZZMass >= ZZMass_PeakCut[1] && ZZMass < 170)) continue;

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

      float KDSystwgt = 1;
      if (tt!=2) KDSystwgt = get_KD_Syst(smearingVal, smearingErr, Systematics);
      KD = compute_KDAlternatives(Txy_BS, Txy_BS_true, KDSystwgt);
      KD = tanh(KD / KD_scale); // Scaled to due to tanh function
      if (KD >= 1.) {
        double correctKD = binwidthx*smeared_edges[1]/(counter_edges[1]+1);
//        if (smeared_edges[1] % 2 == 1) correctKD *=-1.;
//        correctKD += binwidthx/2.;
        KD = 1. - correctKD;
        if (KD>=1) KD = 1-0.1*binwidthx;
        smeared_edges[1]++;
      }
      if (KD <= -1.) {
        double correctKD = binwidthx*smeared_edges[0]/(counter_edges[0]+1);
//        if (smeared_edges[0] % 2 == 1) correctKD *=-1.;
//        correctKD += binwidthx/2.;
        KD = -1. + correctKD;
        if (KD<=-1) KD = -1+0.1*binwidthx;
        smeared_edges[0]++;
      }

			if(tt==2) D_bkg = p0plus_VAJHU*p0plus_m4l / (p0plus_VAJHU*p0plus_m4l + bkg_VAMCFM*bkg_m4l);
			if (D_bkg!=D_bkg) continue;

      double wgt = MC_weight_noxsec;
      if (tt==0) wgt *= MC_weight_QQBGGProper[1]*MC_weight_QQZZEWK;
      else if (tt==1) wgt *= MC_weight_Kfactor*MC_weight_QQBGGProper[0];

//      if (tt==2) cout << ZZMass << '\t' << Txy_BS << '\t' << MC_weight_noxsec << '\t' << CRflag << endl;

      templateWeight = wgt*scale_yield[tt];
      sum_filled += templateWeight;

      int fill_Iterator=1;
      if (tt!=2){
        templateWeight *= 0.5;
        fill_Iterator=2;
      }
      mytree[tt]->Fill();
      if (fill_Iterator==2){
        KD = -KD;
        mytree[tt]->Fill();
      }
    }
    cout << mytree[tt]->GetName() << " events: " << (tt==2 ? sum_filled : sum_filled*luminosity[EnergyIndex]) << endl;
	}

	for (int tt = 0; tt < ntrees; tt++) hfull[ntrees]->Add(hfull[tt]);
	for (int tt = 0; tt < ntrees+1; tt++){
		TH1F* hProjKD = (TH1F*) hfull[tt]->ProjectionX();
		TH2F* hfull_conditional = (TH2F*) hfull[tt]->Clone(Form("%s_conditional",hfull[tt]->GetName()));
		hProjKD->SetTitle(Form("%s Projection",hfull[tt]->GetTitle()));
		hfull_conditional->SetTitle(Form("%s Conditional",hfull[tt]->GetTitle()));
		hproj[tt] = hProjKD;

		for (int biny = 0; biny <= hfull_conditional->GetNbinsY()+1; biny++){
			double integralX = hfull_conditional->Integral(0, hfull_conditional->GetNbinsX() + 1,biny,biny);
			if (integralX != 0){
				for (int binx = 0; binx <= hfull_conditional->GetNbinsX()+1; binx++) hfull_conditional->SetBinContent(binx, biny, (hfull_conditional->GetBinContent(binx, biny) / integralX));
			}
		}

		foutput->WriteTObject(hfull[tt]);
		foutput->WriteTObject(hfull_conditional);
		foutput->WriteTObject(hproj[tt]);
		if(tt<ntrees) foutput->WriteTObject(mytree[tt]);
		delete hfull_conditional;
		delete hfull[tt];
		if (tt < ntrees){
			delete mytree[tt];
			delete tc[tt];
		}
	}
  for (int bs=0; bs<2; bs++) delete tBeam[bs];

	foutput->Close();
  delete hRatio;
  fSSinput->Close();
  delete txy_syst;
  finput_syst->Close();
}


void regularizeKDSyst(TH1F* txy_syst[3]){
  for (int ss=1; ss<3; ss++) txy_syst[ss]->Divide(txy_syst[0]);
  int nbins = txy_syst[0]->GetNbinsX();
  int nbins_remainder = nbins % 2;
  int binmiddle = (nbins+1)/2;
  if (nbins_remainder==1) binmiddle = (nbins+1)/2;
  else binmiddle = nbins/2;
  double first_up = txy_syst[1]->GetBinContent(1);
  double first_dn = txy_syst[2]->GetBinContent(1);
  double first_difference = fabs(first_up - first_dn);
  double last_up = txy_syst[1]->GetBinContent(nbins);
  double last_dn = txy_syst[2]->GetBinContent(nbins);
  double last_difference = fabs(last_up - last_dn);
  double middle_up = txy_syst[1]->GetBinContent(binmiddle);
  double middle_dn = txy_syst[2]->GetBinContent(binmiddle);
  if (nbins_remainder==0){
    middle_up += txy_syst[1]->GetBinContent(binmiddle+1);
    middle_dn += txy_syst[2]->GetBinContent(binmiddle+1);
    middle_up /= 2.;
    middle_dn /= 2.;
  }
  double middle_difference = fabs(middle_up - middle_dn);
  double first_maxdiff = max(first_difference, middle_difference);
  double last_maxdiff = max(last_difference, middle_difference);
  for (int binx=2; binx<binmiddle-1; binx++){ // "first" side
    double ratio_up = txy_syst[1]->GetBinContent(binx);
    double ratio_dn = txy_syst[2]->GetBinContent(binx);
    double difference = fabs(ratio_up - ratio_dn);
    if (
      (ratio_up>=1 && ratio_dn>=1) || (ratio_up<=1 && ratio_dn<=1)
      ){ // Someting went statisticlly wrong
      ratio_up = difference * (first_up-1.) / first_difference + 1.;
      ratio_dn = difference * (first_dn-1.) / first_difference + 1.;
    }

    if (difference>(1.5*first_maxdiff)){
      double scale = (1.5*first_maxdiff) / difference;
      ratio_up = (ratio_up-1.)*scale + 1.;
      ratio_dn = (ratio_dn-1.)*scale + 1.;
    }
    txy_syst[1]->SetBinContent(binx, ratio_up);
    txy_syst[2]->SetBinContent(binx, ratio_dn);
  }
  for (int binx=binmiddle+2; binx<nbins; binx++){ // "last" side
    double ratio_up = txy_syst[1]->GetBinContent(binx);
    double ratio_dn = txy_syst[2]->GetBinContent(binx);
    double difference = fabs(ratio_up - ratio_dn);
    if (
      (ratio_up>=1 && ratio_dn>=1) || (ratio_up<=1 && ratio_dn<=1)
      ){ // Someting went statisticlly wrong
      ratio_up = difference * (first_up-1.) / last_difference + 1.;
      ratio_dn = difference * (first_dn-1.) / last_difference + 1.;
    }

    if (difference>(1.5*last_maxdiff)){
      double scale = (1.5*last_maxdiff) / difference;
      ratio_up = (ratio_up-1.)*scale + 1.;
      ratio_dn = (ratio_dn-1.)*scale + 1.;
    }
    txy_syst[1]->SetBinContent(binx, ratio_up);
    txy_syst[2]->SetBinContent(binx, ratio_dn);
  }
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

