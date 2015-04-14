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
#include "RooAbsPdf.h"
#include "RooGaussian.h"
#include "RooBreitWigner.h"
#include "RooGenericPdf.h"
#include "RooRealSumPdf.h"
#include "RooRealVar.h"
#include "RooFormulaVar.h"
#include "RooDataSet.h"
#include "RooCategory.h"
#include "RooPlot.h"
#include "RooMsgService.h"
#include "RooFitResult.h"
#include "RooSimultaneous.h"
#include "RooWorkspace.h"
#include "./Pdfs/RooRealFlooredSumPdf.h"
#include "./data/ZZ4l_125_6_Samples.h"

using namespace RooFit;
using namespace std;
using namespace ROOT::Math;


const int nSystVars=3;
char* cSystVariable[nSystVars] = { "Txy", "Dbkg", "Txy_BeamSpot" };

char* cSystVariable_label[nSystVars] = {
  "T_{xy} (#mum)", "D_{bkg}", "T^{Beam}_{xy} (#mum)"
};
int nbins_SystVariable[nSystVars] = { 113, 20, 121 };
double KD_limits_SystVariable[nSystVars][2] = {
  { -1412.5, 1412.5 },
  { 0, 1 },
  { -2420, 2420 }
};
double visualRange_SystVariable[nSystVars][2] = {
  { -1412.5, 1412.5 },
  { 0, 1 },
  { -2420, 2420 }
};

const int nMZZ_ranges=4;
float systZZMass_range[nMZZ_ranges][2] = {
  { 70, 105.6 }, { 140.6, 170 }, { 170, 800 }, { 105.6, 140.6 }
};


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



float compute_Dxy(
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

	double unit_pt[3] = {
		cos(ZZPhi),sin(ZZPhi),0
	};
	double cand_cov[9] = {
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
	double opv_cov[9] = {
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

float compute_sigmaxy(
	float ZZPhi,
	float cov_xx,
	float cov_xy,
	float cov_xz,
	float cov_yy,
	float cov_yz,
	float cov_zz
	){

	double unit_pt[3] = {
		cos(ZZPhi), sin(ZZPhi), 0
	};
	double cand_cov[9] = {
		cov_xx,
		cov_xy,
		cov_xz,
		cov_xy,
		cov_yy,
		cov_yz,
		cov_xz,
		cov_yz,
		cov_zz
	};
	TMatrixD candCov(3, 3);
	candCov.SetMatrixArray(cand_cov);

	TVectorD unitVector(3, unit_pt);
	TVectorD intVector = unitVector;
	intVector *= candCov;
	double dot_product = intVector*unitVector;
	return sqrt(dot_product);
}



void constructLifetimeSystTrees(int folder, int erg_tev, int iSyst=0){
  if (iSyst>=nSystVars) return;

  cout << "Starting the program..." << endl;

  char TREE_NAME[]="SelectedTree";
  TString erg_dir;
  erg_dir.Form("LHC_%iTeV", erg_tev);
  TString comstring;
  comstring.Form("%iTeV", erg_tev);
  int EnergyIndex=0; if (erg_tev==8) EnergyIndex=1;

  TString OUTPUT_NAME = "LifetimeKD_SystTrees_DataMC_";
  OUTPUT_NAME.Append("newSIP_");
  OUTPUT_NAME.Append(user_folder[folder]);
  OUTPUT_NAME.Append(Form("_%s_", cSystVariable[iSyst]));
  OUTPUT_NAME.Append(comstring);
  OUTPUT_NAME.Append(".root");

  float ggzz_offshell_sum=0, ggzz_80_100_sum=0;
  float zx_offshell_sum=0, zx_80_100_sum=0, zx_signal_sum=0;
  float qqzz_signal_sum=0, ggzz_signal_sum=0;
  float higgs_signal_sum[4]={ 0 };

  // SIP cut scaling
  string processName[8] ={
    "CTau0", "CTau100", "CTau500", "CTau1000",
    "qqZZ", "ggZZ", "CR", "data"
  };

  string cinput_aux = user_dir_hep + "Analysis/Auxiliary/";
  TString SSinputname = "OSoverSS_";
  SSinputname.Append("MCCR_");
  SSinputname += comstring;
  SSinputname.Append(".root");
  SSinputname.Prepend(cinput_aux.c_str());
  TFile* fSSinput = new TFile(SSinputname, "read");
  TH2F* hRatio;
  TString chratio = "hCR_MC_OSSSRatios_";
  chratio.Append("New cut");
  hRatio = (TH2F*)fSSinput->Get(chratio);
  cout << "Obtained ratio histogram " << hRatio->GetName() << endl;

  double ratio_targetsignalyield[8][5][4] ={ { { 0 } } };
  double targetCRcount[2][5][4] ={ { { 0 } } };
  float m4l_lowhigh[4][2]={
    { 70., 105.6 },
    { 140.6, 170. },
    { 170., 3000. },
    { 105.6, 140.6 }
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
        for (int b = 1; b <= 5; b++) ratio_targetsignalyield[p][b-1][ir] = htemp->GetBinContent(b);
        delete htemp;
      }
      else{
        htemp = (TH1F*)finput->Get(Form("h%s_OS", processName[p].c_str()));
        for (int b = 1; b <= 5; b++) ratio_targetsignalyield[p][b-1][ir] = htemp->GetBinContent(b);
        delete htemp;
        htemp = (TH1F*)finput->Get(Form("h%s_SS", processName[p].c_str()));
        for (int b = 1; b <= 5; b++) ratio_targetsignalyield[p][b-1][ir] += htemp->GetBinContent(b);
        delete htemp;

        htemp = (TH1F*)finput->Get(Form("h%s_Unweighted_OS", processName[p].c_str()));
        for (int b = 1; b <= 5; b++) targetCRcount[0][b-1][ir] = htemp->GetBinContent(b);
        delete htemp;
        htemp = (TH1F*)finput->Get(Form("h%s_Unweighted_SS", processName[p].c_str()));
        for (int b = 1; b <= 5; b++) targetCRcount[1][b-1][ir] = htemp->GetBinContent(b);
        delete htemp;
      }
      finput->Close();
    }
  }
  for (int b = 1; b <= 5; b++){
    ratio_targetsignalyield[5][b-1][0] = ratio_targetsignalyield[5][b-1][3];
    ratio_targetsignalyield[6][b-1][0] = ratio_targetsignalyield[6][b-1][3];
    targetCRcount[0][b-1][0] = targetCRcount[0][b-1][3];
    targetCRcount[1][b-1][0] = targetCRcount[1][b-1][3];
  }
  cout << "Received yield scaling for new vertex cut." << endl;


  float MC_weight;
  float ZXfake_weight_SS;
  float ZXfake_weight_OS[6];

  float MC_weight_noxsec;
  float MC_weight_QQZZEWK=1;
  float MC_weight_Kfactor=1;
  float MC_weight_QQBGGProper[4];
  float ZZMass, ZZPt, ZZEta, ZZPhi;
  float KalmanCandVtx_x, KalmanCandVtx_y, KalmanCandVtx_z, KalmanCandVtx_chi2;
  float OfflinePrimaryVtx_x, OfflinePrimaryVtx_y, OfflinePrimaryVtx_z;
  float GenPrimaryVtx_x, GenPrimaryVtx_y, GenPrimaryVtx_z;
  float GenIntVtx_x, GenIntVtx_y, GenIntVtx_z;
  float GenHMass, GenHPt, GenHPhi;
  float KalmanCandVtx_cov_xx;
  float KalmanCandVtx_cov_xy;
  float KalmanCandVtx_cov_xz;
  float KalmanCandVtx_cov_yy;
  float KalmanCandVtx_cov_yz;
  float KalmanCandVtx_cov_zz;
  float OfflinePrimaryVtx_cov_xx;
  float OfflinePrimaryVtx_cov_xy;
  float OfflinePrimaryVtx_cov_xz;
  float OfflinePrimaryVtx_cov_yy;
  float OfflinePrimaryVtx_cov_yz;
  float OfflinePrimaryVtx_cov_zz;

  int RunNumber, LumiNumber, RunNumber_Ref, LumiNumber_Ref;
  float BeamPosX, BeamPosY, BeamPosZ, BeamPosXErr, BeamPosYErr, BeamPosZErr;

  float Lep1_Z1SIP, Lep2_Z1SIP, Lep3_Z1SIP, Lep4_Z1SIP;
  float Lep1SIP, Lep2SIP, Lep3SIP, Lep4SIP;
  int Lep1ID, Lep2ID, Lep3ID, Lep4ID;

  float D_bkg, p0plus_VAJHU, bkg_VAMCFM, p0plus_m4l, bkg_m4l;
  float Txy, Dxy, Txy_true, Dxy_true, delTxy, delDxy, KD;
  float Txy_BS, Dxy_BS;
  float sigmaPV_xy, sigmaInt_xy;
  float CandVtx_x, CandVtx_y, CandVtx_z;
  float CandBSVtx_x, CandBSVtx_y, CandBSVtx_z;
  float TrueCandVtx_x, TrueCandVtx_y, TrueCandVtx_z;

  int CRflag=-1;
  float Lep1combRelIsoPF=0, Lep2combRelIsoPF=0, Lep3combRelIsoPF=0, Lep4combRelIsoPF=0;
  bool Lep3isID=true, Lep4isID=true;
  short Z1ids;

  TString cinput_common = user_dir_hep + erg_dir + "/";
  TString cinput_common_noSIP = user_dir_hep + "No_SIP/" + erg_dir + "/";
  TString coutput_common = user_dir_hep + "Analysis/ShapeSystematics/";
  coutput_common.Append(Form("%s/", cSystVariable[iSyst]));
  TString mkdirCommand = "mkdir -p ";
  mkdirCommand.Append(coutput_common);
  gSystem->Exec(mkdirCommand);
  TString coutput = coutput_common + OUTPUT_NAME;
  TFile* foutput = new TFile(coutput, "recreate");
  cout << "Opened the output file." << endl;

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

  const int ntrees=4;
  TString cinput_qqzz_common = cinput_common + user_folder[folder] + "/";
  TString cinput_ggzz_common = cinput_qqzz_common;
  TString cinput_zx_common = cinput_common;
  TString cinput_qqzz_common_noSIP = cinput_common_noSIP + user_folder[folder] + "/";
  TString cinput_ggzz_common_noSIP = cinput_qqzz_common_noSIP;
  TString cinput_zx_common_noSIP = cinput_common_noSIP;

  TChain* tqqzz = new TChain(TREE_NAME);
  for (int smp = kGGSamples; smp < kQQBZZSamples; smp++){
    TString cinput_qqzz = cinput_qqzz_common_noSIP + sample_FullSim[smp] + "_Reprocessed.root";
    tqqzz->Add(cinput_qqzz);
    cout << "Acquired " << cinput_qqzz << endl;
  }
  TChain* tggzz = new TChain(TREE_NAME);
  for (int smp = kGGOLDSamples; smp < kGGMCFMSamples; smp++){
    TString cinput_ggzz = cinput_ggzz_common_noSIP + sample_FullSim[smp] + "_Reprocessed.root";
    tggzz->Add(cinput_ggzz);
    cout << "Acquired " << cinput_ggzz << endl;
  }
  TChain* tzx = new TChain(TREE_NAME);
  TString cinput_zx = user_dir_hep;
  cinput_zx.Append("No_SIP/");
  cinput_zx.Append(erg_dir);
  cinput_zx = cinput_zx + "/CR/" + sample_FullSim[kAllSamples - 1] + ".root";
  tzx->Add(cinput_zx);
  cout << "Acquired " << cinput_zx << endl;
  
  TString cinput_data = cinput_common_noSIP + user_folder[3] + "/" + data_files[folder] + ".root";
  TChain* tdata = new TChain(TREE_NAME);
  tdata->Add(cinput_data);
  cout << "Acquired " << cinput_data << endl;
  TChain* tc[ntrees] ={ tdata, tqqzz, tggzz, tzx };
  cout << "Obtained tc." << endl;

  const int ntrees_extra=6;
  TString cinput_higgs_common = cinput_qqzz_common;
  TString cinput_higgs_common_noSIP = cinput_qqzz_common_noSIP;
  TChain* tqqzz_extra = new TChain(TREE_NAME);
  for (int smp = kQQBZZSamples; smp < kQQBZZSamples_Dedicated; smp++){
    TString cinput_qqzz = cinput_qqzz_common_noSIP + sample_FullSim[smp] + "_Reprocessed.root";
    tqqzz_extra->Add(cinput_qqzz);
  }
  TChain* tggzz_extra = new TChain(TREE_NAME);
  for (int smp = kGGHSamples; smp < kGGSamples; smp++){
    if (smp >= kGGOLDSamples && smp < kGGMCFMSamples) continue;
    TString cinput_ggzz = cinput_ggzz_common_noSIP + sample_FullSim[smp] + "_Reprocessed.root";
    tggzz_extra->Add(cinput_ggzz);
  }
  TChain* thiggs[kGGHSamples];
  for (int smp = 0; smp < kGGHSamples; smp++){
    TString cinput_higgs;
    cinput_higgs = cinput_higgs_common_noSIP + sample_FullSim[smp] + ".root";
    cout << "Attaching Higgs under " << cinput_higgs << endl;
    thiggs[smp] = new TChain(TREE_NAME);
    thiggs[smp]->Add(cinput_higgs);
  }
  TChain* tc_extra[ntrees_extra] ={ thiggs[0], tqqzz_extra, tggzz_extra, thiggs[1], thiggs[2], thiggs[3] };
  cout << "Obtained tc_extra." << endl;

  int nbins_KD = nbins_SystVariable[iSyst];
  double KD_limits[2] ={ KD_limits_SystVariable[iSyst][0], KD_limits_SystVariable[iSyst][1] };

  TH2F* hdata_full = new TH2F("hdata_full", Form("%i TeV %s Data", erg_tev, user_folder[folder]), nMZZ_ranges, 0, nMZZ_ranges, nbins_KD, KD_limits[0], KD_limits[1]);
  hdata_full->GetYaxis()->SetTitle(cSystVariable_label[iSyst]);
  hdata_full->GetXaxis()->SetTitle("m_{4l} (GeV)");
  hdata_full->GetXaxis()->SetBinLabel(1, "70 - 105.6");
  hdata_full->GetXaxis()->SetBinLabel(2, "140.6 - 170");
  hdata_full->GetXaxis()->SetBinLabel(3, "170 - 800");
  hdata_full->GetXaxis()->SetBinLabel(4, "105.6 - 140.6");
  hdata_full->SetOption("colz");
  TH2F* hqqzz_full = (TH2F*)hdata_full->Clone();
  TH2F* hggzz_full = (TH2F*)hdata_full->Clone();
  TH2F* hzx_full = (TH2F*)hdata_full->Clone();
  TH2F* hmc_full = (TH2F*)hdata_full->Clone();
  hqqzz_full->SetNameTitle("hqqzz_full", Form("%i TeV %s q#bar{q}#rightarrow4l Background", erg_tev, user_folder[folder]));
  hggzz_full->SetNameTitle("hggzz_full", Form("%i TeV %s gg Background", erg_tev, user_folder[folder]));
  hzx_full->SetNameTitle("hzx_full", Form("%i TeV %s Z+X Background", erg_tev, user_folder[folder]));
  hmc_full->SetNameTitle("hmc_full", Form("%i TeV %s Total Background", erg_tev, user_folder[folder]));
  TH2F* hfull[ntrees+1] ={ hdata_full, hqqzz_full, hggzz_full, hzx_full, hmc_full };
  for (int tt = 0; tt < ntrees; tt++) hfull[tt]->Sumw2();
  TH1F* hproj[ntrees+1];

  TH1F* hhiggs[kGGHSamples];
  for (int hh = 0; hh < kGGHSamples; hh++){
    hhiggs[hh] = (TH1F*)hdata_full->ProjectionY();
    hhiggs[hh]->SetNameTitle(Form("hhiggs_CTau%.0f_full", sample_CTau[hh]), Form("%i TeV %s SM Higgs (c#Tau=%.0f #mum)", erg_tev, user_folder[folder], sample_CTau[hh]));
    hhiggs[hh]->Sumw2();
  }

  for (int tt = 0; tt < ntrees; tt++){
    cout << "Tree: " << tt << ", nEvents: " << tc[tt]->GetEntries() << endl;
    if (tt==0 || tt==3){
      tc[tt]->SetBranchAddress("RunNumber", &RunNumber);
      tc[tt]->SetBranchAddress("LumiNumber", &LumiNumber);
    }
    if (tt>0){
      if (tt != 3){
        tc[tt]->SetBranchAddress("MC_weight", &MC_weight);
        tc[tt]->SetBranchAddress("MC_weight_noxsec", &MC_weight_noxsec);
      }
      else{
        if (tc[tt]->GetBranchStatus("ZXfake_weight_SS")){
          tc[tt]->SetBranchAddress("ZXfake_weight_OS", ZXfake_weight_OS);
          tc[tt]->SetBranchAddress("ZXfake_weight_SS", &ZXfake_weight_SS);
        }
        else tc[tt]->SetBranchAddress("ZXfake_weight", &MC_weight);
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
      }
    }
    if (tt == 2 || tt == 1){
      tc[tt]->SetBranchAddress("MC_weight_Kfactor", &MC_weight_Kfactor);
      tc[tt]->SetBranchAddress("MC_weight_QQBGGProper", MC_weight_QQBGGProper);
      if (tt == 1){
        tc[tt]->SetBranchAddress("MC_weight_QQZZEWK", &MC_weight_QQZZEWK);
      }
    }

    tc[tt]->SetBranchAddress("ZZMass", &ZZMass);
    tc[tt]->SetBranchAddress("ZZPt", &ZZPt);
    tc[tt]->SetBranchAddress("ZZEta", &ZZEta);
    tc[tt]->SetBranchAddress("ZZPhi", &ZZPhi);
    tc[tt]->SetBranchAddress("OfflinePrimaryVtx_x", &OfflinePrimaryVtx_x);
    tc[tt]->SetBranchAddress("OfflinePrimaryVtx_y", &OfflinePrimaryVtx_y);
    tc[tt]->SetBranchAddress("OfflinePrimaryVtx_z", &OfflinePrimaryVtx_z);
    tc[tt]->SetBranchAddress("KalmanCandVtx_x", &KalmanCandVtx_x);
    tc[tt]->SetBranchAddress("KalmanCandVtx_y", &KalmanCandVtx_y);
    tc[tt]->SetBranchAddress("KalmanCandVtx_z", &KalmanCandVtx_z);
    tc[tt]->SetBranchAddress("OfflinePrimaryVtx_cov_xx", &OfflinePrimaryVtx_cov_xx);
    tc[tt]->SetBranchAddress("OfflinePrimaryVtx_cov_xy", &OfflinePrimaryVtx_cov_xy);
    tc[tt]->SetBranchAddress("OfflinePrimaryVtx_cov_xz", &OfflinePrimaryVtx_cov_xz);
    tc[tt]->SetBranchAddress("OfflinePrimaryVtx_cov_yy", &OfflinePrimaryVtx_cov_yy);
    tc[tt]->SetBranchAddress("OfflinePrimaryVtx_cov_yz", &OfflinePrimaryVtx_cov_yz);
    tc[tt]->SetBranchAddress("OfflinePrimaryVtx_cov_zz", &OfflinePrimaryVtx_cov_zz);
    tc[tt]->SetBranchAddress("KalmanCandVtx_cov_xx", &KalmanCandVtx_cov_xx);
    tc[tt]->SetBranchAddress("KalmanCandVtx_cov_xy", &KalmanCandVtx_cov_xy);
    tc[tt]->SetBranchAddress("KalmanCandVtx_cov_xz", &KalmanCandVtx_cov_xz);
    tc[tt]->SetBranchAddress("KalmanCandVtx_cov_yy", &KalmanCandVtx_cov_yy);
    tc[tt]->SetBranchAddress("KalmanCandVtx_cov_yz", &KalmanCandVtx_cov_yz);
    tc[tt]->SetBranchAddress("KalmanCandVtx_cov_zz", &KalmanCandVtx_cov_zz);

    if (!tc[tt]->GetBranchStatus("GenPrimaryVtx_x")) cout << "Tree " << tt << " has no gen. vtx" << endl;
    else{
      cout << "Tree " << tt << " setting gen. vtx" << endl;
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
    tc[tt]->SetBranchAddress("KalmanCandVtx_chi2", &KalmanCandVtx_chi2);
    tc[tt]->SetBranchAddress("Lep1_Z1SIP", &Lep1_Z1SIP);
    tc[tt]->SetBranchAddress("Lep2_Z1SIP", &Lep2_Z1SIP);
    tc[tt]->SetBranchAddress("Lep3_Z1SIP", &Lep3_Z1SIP);
    tc[tt]->SetBranchAddress("Lep4_Z1SIP", &Lep4_Z1SIP);
    tc[tt]->SetBranchAddress("Lep1SIP", &Lep1SIP);
    tc[tt]->SetBranchAddress("Lep2SIP", &Lep2SIP);
    tc[tt]->SetBranchAddress("Lep3SIP", &Lep3SIP);
    tc[tt]->SetBranchAddress("Lep4SIP", &Lep4SIP);
    if (tc[tt]->GetBranchStatus("D_bkg")){
      tc[tt]->SetBranchAddress("D_bkg", &D_bkg);
    }
    else if (tc[tt]->GetBranchStatus("p0plus_VAJHU")){
      tc[tt]->SetBranchAddress("p0plus_VAJHU", &p0plus_VAJHU);
      tc[tt]->SetBranchAddress("bkg_VAMCFM", &bkg_VAMCFM);
      tc[tt]->SetBranchAddress("p0plus_m4l", &p0plus_m4l);
      tc[tt]->SetBranchAddress("bkg_m4l", &bkg_m4l);
    }
    else cout << "Could not find p0plus_VAJHU in tc[" << tt << "]!!!" << endl;
  }

  for (int tt = 0; tt < ntrees_extra; tt++){
    cout << "Extra tree: " << tt << ", nEvents: " << tc_extra[tt]->GetEntries() << endl;
    tc_extra[tt]->SetBranchAddress("MC_weight", &MC_weight);
    tc_extra[tt]->SetBranchAddress("MC_weight_noxsec", &MC_weight_noxsec);
    if (tt == 2 || tt == 1){
      tc_extra[tt]->SetBranchAddress("MC_weight_Kfactor", &MC_weight_Kfactor);
      tc_extra[tt]->SetBranchAddress("MC_weight_QQBGGProper", MC_weight_QQBGGProper);
      if (tt == 1){
        tc_extra[tt]->SetBranchAddress("MC_weight_QQZZEWK", &MC_weight_QQZZEWK);
      }
    }
    tc_extra[tt]->SetBranchAddress("ZZMass", &ZZMass);
    tc_extra[tt]->SetBranchAddress("ZZPt", &ZZPt);
    tc_extra[tt]->SetBranchAddress("ZZEta", &ZZEta);
    tc_extra[tt]->SetBranchAddress("ZZPhi", &ZZPhi);
    tc_extra[tt]->SetBranchAddress("OfflinePrimaryVtx_x", &OfflinePrimaryVtx_x);
    tc_extra[tt]->SetBranchAddress("OfflinePrimaryVtx_y", &OfflinePrimaryVtx_y);
    tc_extra[tt]->SetBranchAddress("OfflinePrimaryVtx_z", &OfflinePrimaryVtx_z);
    tc_extra[tt]->SetBranchAddress("KalmanCandVtx_x", &KalmanCandVtx_x);
    tc_extra[tt]->SetBranchAddress("KalmanCandVtx_y", &KalmanCandVtx_y);
    tc_extra[tt]->SetBranchAddress("KalmanCandVtx_z", &KalmanCandVtx_z);
    tc_extra[tt]->SetBranchAddress("OfflinePrimaryVtx_cov_xx", &OfflinePrimaryVtx_cov_xx);
    tc_extra[tt]->SetBranchAddress("OfflinePrimaryVtx_cov_xy", &OfflinePrimaryVtx_cov_xy);
    tc_extra[tt]->SetBranchAddress("OfflinePrimaryVtx_cov_xz", &OfflinePrimaryVtx_cov_xz);
    tc_extra[tt]->SetBranchAddress("OfflinePrimaryVtx_cov_yy", &OfflinePrimaryVtx_cov_yy);
    tc_extra[tt]->SetBranchAddress("OfflinePrimaryVtx_cov_yz", &OfflinePrimaryVtx_cov_yz);
    tc_extra[tt]->SetBranchAddress("OfflinePrimaryVtx_cov_zz", &OfflinePrimaryVtx_cov_zz);
    tc_extra[tt]->SetBranchAddress("KalmanCandVtx_cov_xx", &KalmanCandVtx_cov_xx);
    tc_extra[tt]->SetBranchAddress("KalmanCandVtx_cov_xy", &KalmanCandVtx_cov_xy);
    tc_extra[tt]->SetBranchAddress("KalmanCandVtx_cov_xz", &KalmanCandVtx_cov_xz);
    tc_extra[tt]->SetBranchAddress("KalmanCandVtx_cov_yy", &KalmanCandVtx_cov_yy);
    tc_extra[tt]->SetBranchAddress("KalmanCandVtx_cov_yz", &KalmanCandVtx_cov_yz);
    tc_extra[tt]->SetBranchAddress("KalmanCandVtx_cov_zz", &KalmanCandVtx_cov_zz);

    if (!tc_extra[tt]->GetBranchStatus("GenPrimaryVtx_x")) cout << "Extra tree " << tt << " has no gen. vtx" << endl;
    else{
      cout << "Extra tree " << tt << " setting gen. vtx" << endl;
      tc_extra[tt]->SetBranchAddress("GenPrimaryVtx_x", &GenPrimaryVtx_x);
      tc_extra[tt]->SetBranchAddress("GenPrimaryVtx_y", &GenPrimaryVtx_y);
      tc_extra[tt]->SetBranchAddress("GenPrimaryVtx_z", &GenPrimaryVtx_z);
      tc_extra[tt]->SetBranchAddress("GenIntVtx_x", &GenIntVtx_x);
      tc_extra[tt]->SetBranchAddress("GenIntVtx_y", &GenIntVtx_y);
      tc_extra[tt]->SetBranchAddress("GenIntVtx_z", &GenIntVtx_z);
      tc_extra[tt]->SetBranchAddress("GenHPt", &GenHPt);
      tc_extra[tt]->SetBranchAddress("GenHPhi", &GenHPhi);
      tc_extra[tt]->SetBranchAddress("GenHMass", &GenHMass);
    }
    tc_extra[tt]->SetBranchAddress("Lep1_Z1SIP", &Lep1_Z1SIP);
    tc_extra[tt]->SetBranchAddress("Lep2_Z1SIP", &Lep2_Z1SIP);
    tc_extra[tt]->SetBranchAddress("Lep3_Z1SIP", &Lep3_Z1SIP);
    tc_extra[tt]->SetBranchAddress("Lep4_Z1SIP", &Lep4_Z1SIP);
    tc_extra[tt]->SetBranchAddress("KalmanCandVtx_chi2", &KalmanCandVtx_chi2);
    tc_extra[tt]->SetBranchAddress("Lep1SIP", &Lep1SIP);
    tc_extra[tt]->SetBranchAddress("Lep2SIP", &Lep2SIP);
    tc_extra[tt]->SetBranchAddress("Lep3SIP", &Lep3SIP);
    tc_extra[tt]->SetBranchAddress("Lep4SIP", &Lep4SIP);
    if (tc_extra[tt]->GetBranchStatus("D_bkg")){
      tc_extra[tt]->SetBranchAddress("D_bkg", &D_bkg);
    }
    else if (tc_extra[tt]->GetBranchStatus("p0plus_VAJHU")){
      tc_extra[tt]->SetBranchAddress("p0plus_VAJHU", &p0plus_VAJHU);
      tc_extra[tt]->SetBranchAddress("bkg_VAMCFM", &bkg_VAMCFM);
      tc_extra[tt]->SetBranchAddress("p0plus_m4l", &p0plus_m4l);
      tc_extra[tt]->SetBranchAddress("bkg_m4l", &bkg_m4l);
    }
    else cout << "Could not find p0plus_VAJHU in tc_extra[" << tt << "]!!!" << endl;
  }

  for (int tt = 0; tt < ntrees; tt++){
    for (int ev = 0; ev < tc[tt]->GetEntries(); ev++){
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

      tc[tt]->GetEntry(ev);
      if ((fabs(Lep1_Z1SIP) >= 4 || fabs(Lep2_Z1SIP) >= 4 || fabs(Lep3_Z1SIP) >= 5 || fabs(Lep4_Z1SIP) >= 5 || KalmanCandVtx_chi2 >= 30)) continue;
      if (
        tt==3 &&
        !applyZXloose(CRflag,
        Lep1combRelIsoPF,
        Lep2combRelIsoPF,
        Lep3combRelIsoPF,
        Lep4combRelIsoPF,
        Lep3isID,
        Lep4isID)
        ) continue;
      if (tt==3 && !ZXchannelselection(folder, Z1ids, CRflag)) continue;
      if (tt==3 && (
        (CRflag==6 || CRflag==10) ||
        (CRflag==8 || CRflag==12)
        )){
        int regionIndex = 3;
        for (int ir=0; ir<4; ir++){
          if (ZZMass<m4l_lowhigh[ir][1] && ZZMass>=m4l_lowhigh[ir][0]){
            regionIndex = ir; break;
          }
        }
        if (ZZMass<m4l_lowhigh[0][0]) regionIndex = 0;
        if (ZZMass>=m4l_lowhigh[2][1]) regionIndex = 2;
        if (tc[tt]->GetBranchStatus("ZXfake_weight_OS")) MC_weight = ZXfake_weight_OS[4] * (targetCRcount[0][4][regionIndex] / (targetCRcount[0][4][regionIndex] + targetCRcount[1][4][regionIndex]));
      }
      else if (tt==3 && (
        (CRflag==7 || CRflag==11) ||
        (CRflag==5 || CRflag==9)
        )){
        if (tc[tt]->GetBranchStatus("ZXfake_weight_SS")){
          int OSoverSS_biny = -1;
          if ((CRflag==5 || CRflag==9) && Z1ids==-169) OSoverSS_biny = 1; // 4mu
          if ((CRflag==7 || CRflag==11) && Z1ids==-169) OSoverSS_biny = 3; // 2mu2e
          if ((CRflag==5 || CRflag==9) && Z1ids==-121) OSoverSS_biny = 4; // 2e2mu
          if ((CRflag==7 || CRflag==11) && Z1ids==-121) OSoverSS_biny = 6; // 4e
          float SSwgt_OSoverSS = 1;
          if (OSoverSS_biny>0) SSwgt_OSoverSS = hRatio->GetBinContent(hRatio->GetXaxis()->FindBin(ZZMass), OSoverSS_biny);

          int regionIndex = 3;
          for (int ir=0; ir<4; ir++){
            if (ZZMass<m4l_lowhigh[ir][1] && ZZMass>=m4l_lowhigh[ir][0]){
              regionIndex = ir; break;
            }
          }
          if (ZZMass<m4l_lowhigh[0][0]) regionIndex = 0;
          if (ZZMass>=m4l_lowhigh[2][1]) regionIndex = 2;

          MC_weight = ZXfake_weight_SS*SSwgt_OSoverSS;
          MC_weight *= (targetCRcount[1][4][regionIndex] / (targetCRcount[0][4][regionIndex] + targetCRcount[1][4][regionIndex]));
        }
      }

      if (iSyst>=2 && iSyst<nSystVars){ // BE CAREFUL HERE, HARDCODED BOUNDS AHEAD
        int index_BS = 1;
        bool matchFound = false;
        if (tt==0 || tt==3) index_BS = 0;
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
      }

      CandVtx_x = KalmanCandVtx_x - OfflinePrimaryVtx_x;
      CandVtx_y = KalmanCandVtx_y - OfflinePrimaryVtx_y;
      CandVtx_z = KalmanCandVtx_z - OfflinePrimaryVtx_z;

      TrueCandVtx_x = GenIntVtx_x - GenPrimaryVtx_x;
      TrueCandVtx_y = GenIntVtx_y - GenPrimaryVtx_y;
      TrueCandVtx_z = GenIntVtx_z - GenPrimaryVtx_z;

      Dxy = (CandVtx_x*cos(ZZPhi) + CandVtx_y*sin(ZZPhi));
      Dxy_true = (TrueCandVtx_x*cos(GenHPhi) + TrueCandVtx_y*sin(GenHPhi));
      Txy = Dxy*ZZMass / ZZPt;
      Txy_true = Dxy_true*GenHMass / GenHPt;

      delDxy = compute_Dxy(
        ZZPhi,
        KalmanCandVtx_cov_xx,
        KalmanCandVtx_cov_xy,
        KalmanCandVtx_cov_xz,
        KalmanCandVtx_cov_yy,
        KalmanCandVtx_cov_yz,
        KalmanCandVtx_cov_zz,
        OfflinePrimaryVtx_cov_xx,
        OfflinePrimaryVtx_cov_xy,
        OfflinePrimaryVtx_cov_xz,
        OfflinePrimaryVtx_cov_yy,
        OfflinePrimaryVtx_cov_yz,
        OfflinePrimaryVtx_cov_zz
        );
      sigmaPV_xy = compute_sigmaxy(
        ZZPhi,
        OfflinePrimaryVtx_cov_xx,
        OfflinePrimaryVtx_cov_xy,
        OfflinePrimaryVtx_cov_xz,
        OfflinePrimaryVtx_cov_yy,
        OfflinePrimaryVtx_cov_yz,
        OfflinePrimaryVtx_cov_zz
        );
      sigmaInt_xy = compute_sigmaxy(
        ZZPhi,
        KalmanCandVtx_cov_xx,
        KalmanCandVtx_cov_xy,
        KalmanCandVtx_cov_xz,
        KalmanCandVtx_cov_yy,
        KalmanCandVtx_cov_yz,
        KalmanCandVtx_cov_zz
        );
      delTxy = delDxy*ZZMass / ZZPt;

      if (!tc[tt]->GetBranchStatus("D_bkg")) D_bkg = p0plus_VAJHU*p0plus_m4l / (p0plus_VAJHU*p0plus_m4l + bkg_VAMCFM*bkg_m4l);

      if (iSyst==0) KD = Txy;
      if (iSyst==2) KD = Txy_BS;
// Some more KDs in cm
      KD = KD*10000.; // in 1 um
// Some other unitless KDs
      if (iSyst==1){
        KD = D_bkg; if (KD==1.) KD=0.999999;
      }

      float wgt = 1;

      int biny = hfull[tt]->GetYaxis()->FindBin(KD);
      for (int binx = 0; binx < nMZZ_ranges; binx++){
        if (tt == 3 && binx == 0 && ZZMass >= systZZMass_range[binx][0] + 30 && ZZMass < systZZMass_range[binx][1] + 30){ // ZX special fill
          wgt = MC_weight;
          hfull[tt]->AddBinContent(hfull[tt]->GetBin(binx + 1, biny), wgt);
          hfull[tt]->SetBinError(hfull[tt]->GetBin(binx + 1, biny), sqrt(wgt*wgt + pow(hfull[tt]->GetBinError(hfull[tt]->GetBin(binx + 1, biny)), 2)));
          if (ZZMass<130 && ZZMass>=110){
            zx_80_100_sum += wgt;
          }
          continue;
        }
        else if (tt == 3 && binx == 0) continue;
        if (tt == 2 && binx == 0) continue;
        if (tt == 1 && binx == 0 && ZZMass >= systZZMass_range[binx][0] && ZZMass < systZZMass_range[binx][1]){ // gg bkg special fill
          wgt = MC_weight*MC_weight_Kfactor*MC_weight_QQBGGProper[0];
          hfull[2]->AddBinContent(hfull[2]->GetBin(binx + 1, biny), wgt);
          hfull[2]->SetBinError(hfull[2]->GetBin(binx + 1, biny), sqrt(wgt*wgt + pow(hfull[2]->GetBinError(hfull[2]->GetBin(binx + 1, biny)), 2)));
          if (ZZMass>=80 && ZZMass<100) ggzz_80_100_sum += wgt;
        }
        if (binx < nMZZ_ranges - 1){
          wgt = MC_weight;
          if (tt == 1) wgt *= MC_weight_QQZZEWK;
          if (tt == 2) wgt *= MC_weight_Kfactor;
          if (ZZMass < systZZMass_range[binx][1] && ZZMass >= systZZMass_range[binx][0]){
            hfull[tt]->AddBinContent(hfull[tt]->GetBin(binx + 1, biny), wgt);
            hfull[tt]->SetBinError(hfull[tt]->GetBin(binx + 1, biny), sqrt(wgt*wgt + pow(hfull[tt]->GetBinError(hfull[tt]->GetBin(binx + 1, biny)), 2)));
          }
          if (tt == 2 && ZZMass >= 220 && ZZMass<1600 && binx==2) ggzz_offshell_sum += wgt;
          if (tt == 3 && ZZMass >= 220 && ZZMass<1600 && binx==2) zx_offshell_sum += wgt;
        }
        else{ // Fill signal region
          if (ZZMass >= systZZMass_range[binx][1] || ZZMass < systZZMass_range[binx][0]) continue;
          wgt = ((tt==0 || tt>=3) ? MC_weight : MC_weight_noxsec);
          if (tt == 1){
            wgt = MC_weight_noxsec*MC_weight_Kfactor*MC_weight_QQBGGProper[0];
            hfull[2]->AddBinContent(hfull[2]->GetBin(binx + 1, biny), wgt);
            hfull[2]->SetBinError(hfull[2]->GetBin(binx + 1, biny), sqrt(wgt*wgt + pow(hfull[2]->GetBinError(hfull[2]->GetBin(binx + 1, biny)), 2)));
            ggzz_signal_sum += wgt;
            wgt = MC_weight_noxsec*MC_weight_QQBGGProper[1];
            hfull[1]->AddBinContent(hfull[1]->GetBin(binx + 1, biny), wgt);
            hfull[1]->SetBinError(hfull[1]->GetBin(binx + 1, biny), sqrt(wgt*wgt + pow(hfull[1]->GetBinError(hfull[1]->GetBin(binx + 1, biny)), 2)));
            qqzz_signal_sum += wgt;
          }
          else if (tt == 2){
            wgt *= MC_weight_Kfactor*MC_weight_QQBGGProper[0];
            hfull[tt]->AddBinContent(hfull[tt]->GetBin(binx + 1, biny), wgt);
            hfull[tt]->SetBinError(hfull[tt]->GetBin(binx + 1, biny), sqrt(wgt*wgt + pow(hfull[tt]->GetBinError(hfull[tt]->GetBin(binx + 1, biny)), 2)));
            ggzz_signal_sum += wgt;
          }
          else{
            hfull[tt]->AddBinContent(hfull[tt]->GetBin(binx + 1, biny), wgt);
            hfull[tt]->SetBinError(hfull[tt]->GetBin(binx + 1, biny), sqrt(wgt*wgt + pow(hfull[tt]->GetBinError(hfull[tt]->GetBin(binx + 1, biny)), 2)));
          }
          if (tt == 3) zx_signal_sum += wgt;
        }
      }
    }
  }

  for (int tt = 0; tt < ntrees_extra; tt++){
    for (int ev = 0; ev < tc_extra[tt]->GetEntries(); ev++){
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

      tc_extra[tt]->GetEntry(ev);
      if ((fabs(Lep1_Z1SIP) >= 4 || fabs(Lep2_Z1SIP) >= 4 || fabs(Lep3_Z1SIP) >= 5 || fabs(Lep4_Z1SIP) >= 5 || KalmanCandVtx_chi2 >= 30)) continue;

      int binx = nMZZ_ranges-1;
      if (ZZMass >= systZZMass_range[binx][1] || ZZMass < systZZMass_range[binx][0]) continue;

      if (iSyst>=2 && iSyst<nSystVars){ // BE CAREFUL HERE, HARDCODED BOUNDS AHEAD
        int index_BS = 1; // These are all MC.
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
      }

      CandVtx_x = KalmanCandVtx_x - OfflinePrimaryVtx_x;
      CandVtx_y = KalmanCandVtx_y - OfflinePrimaryVtx_y;
      CandVtx_z = KalmanCandVtx_z - OfflinePrimaryVtx_z;

      TrueCandVtx_x = GenIntVtx_x - GenPrimaryVtx_x;
      TrueCandVtx_y = GenIntVtx_y - GenPrimaryVtx_y;
      TrueCandVtx_z = GenIntVtx_z - GenPrimaryVtx_z;

      Dxy = (CandVtx_x*cos(ZZPhi) + CandVtx_y*sin(ZZPhi));
      Dxy_true = (TrueCandVtx_x*cos(GenHPhi) + TrueCandVtx_y*sin(GenHPhi));
      //			if(Dxy_true<0) cout << "Warning! Dxy_true = " << Dxy_true << " x/y: " << TrueCandVtx_x/TrueCandVtx_y << " cot(phi): " << cos(GenHPhi)/sin(GenHPhi) << endl;
      Txy = Dxy*ZZMass / ZZPt;
      Txy_true = Dxy_true*GenHMass / GenHPt;

      delDxy = compute_Dxy(
        ZZPhi,
        KalmanCandVtx_cov_xx,
        KalmanCandVtx_cov_xy,
        KalmanCandVtx_cov_xz,
        KalmanCandVtx_cov_yy,
        KalmanCandVtx_cov_yz,
        KalmanCandVtx_cov_zz,
        OfflinePrimaryVtx_cov_xx,
        OfflinePrimaryVtx_cov_xy,
        OfflinePrimaryVtx_cov_xz,
        OfflinePrimaryVtx_cov_yy,
        OfflinePrimaryVtx_cov_yz,
        OfflinePrimaryVtx_cov_zz
        );
      sigmaPV_xy = compute_sigmaxy(
        ZZPhi,
        OfflinePrimaryVtx_cov_xx,
        OfflinePrimaryVtx_cov_xy,
        OfflinePrimaryVtx_cov_xz,
        OfflinePrimaryVtx_cov_yy,
        OfflinePrimaryVtx_cov_yz,
        OfflinePrimaryVtx_cov_zz
        );
      sigmaInt_xy = compute_sigmaxy(
        ZZPhi,
        KalmanCandVtx_cov_xx,
        KalmanCandVtx_cov_xy,
        KalmanCandVtx_cov_xz,
        KalmanCandVtx_cov_yy,
        KalmanCandVtx_cov_yz,
        KalmanCandVtx_cov_zz
        );
      delTxy = delDxy*ZZMass / ZZPt;

      if (!tc_extra[tt]->GetBranchStatus("D_bkg")) D_bkg = p0plus_VAJHU*p0plus_m4l / (p0plus_VAJHU*p0plus_m4l + bkg_VAMCFM*bkg_m4l);

      if (iSyst==0) KD = Txy;
      if (iSyst==2) KD = Txy_BS;
      KD = KD*10000.; // in 1 um

      if (iSyst==1){
        KD = D_bkg; if (KD==1.) KD=0.999999;
      }

      float wgt = 1;

      int biny = hhiggs[0]->GetXaxis()->FindBin(KD);
      wgt = ((tt==0 || tt>=3) ? MC_weight : MC_weight_noxsec);
      if (tt == 1){
        wgt = MC_weight_noxsec*MC_weight_Kfactor*MC_weight_QQBGGProper[0];
        hfull[2]->AddBinContent(hfull[2]->GetBin(binx + 1, biny), wgt);
        hfull[2]->SetBinError(hfull[2]->GetBin(binx + 1, biny), sqrt(wgt*wgt + pow(hfull[2]->GetBinError(hfull[2]->GetBin(binx + 1, biny)), 2)));
        ggzz_signal_sum += wgt;
        wgt = MC_weight_noxsec*MC_weight_QQBGGProper[1];
        hfull[1]->AddBinContent(hfull[1]->GetBin(binx + 1, biny), wgt);
        hfull[1]->SetBinError(hfull[1]->GetBin(binx + 1, biny), sqrt(wgt*wgt + pow(hfull[1]->GetBinError(hfull[1]->GetBin(binx + 1, biny)), 2)));
        qqzz_signal_sum += wgt;
      }
      else if (tt == 2){
        wgt *= MC_weight_Kfactor*MC_weight_QQBGGProper[0];
        hfull[tt]->AddBinContent(hfull[tt]->GetBin(binx + 1, biny), wgt);
        hfull[tt]->SetBinError(hfull[tt]->GetBin(binx + 1, biny), sqrt(wgt*wgt + pow(hfull[tt]->GetBinError(hfull[tt]->GetBin(binx + 1, biny)), 2)));
        ggzz_signal_sum += wgt;
      }
      else if (tt == 0){
        hhiggs[tt]->AddBinContent(biny, wgt);
        hhiggs[tt]->SetBinError(biny, sqrt(wgt*wgt + pow(hhiggs[tt]->GetBinError(biny), 2)));
        higgs_signal_sum[tt] += wgt;
      }
      else if (tt >= 3){
        hhiggs[tt-2]->AddBinContent(biny, wgt);
        hhiggs[tt-2]->SetBinError(biny, sqrt(wgt*wgt + pow(hhiggs[tt-2]->GetBinError(biny), 2)));
        higgs_signal_sum[tt-2] += wgt;
      }
    }
  }

  for (int biny = 0; biny <= hggzz_full->GetNbinsY()+1; biny++){
    hggzz_full->SetBinContent(1, biny, (hggzz_full->GetBinContent(1, biny))*(yield_80_100_ggzz[EnergyIndex][folder] / ggzz_80_100_sum / luminosity[EnergyIndex]));
    hggzz_full->SetBinContent(4, biny, (hggzz_full->GetBinContent(4, biny))*(yield_signal_ggzz[EnergyIndex][folder] / ggzz_signal_sum / luminosity[EnergyIndex]));
    hggzz_full->SetBinContent(3, biny, (hggzz_full->GetBinContent(3, biny))*(yield_offshell_ggzz[EnergyIndex][folder] / ggzz_offshell_sum / luminosity[EnergyIndex]));
    hggzz_full->SetBinContent(2, biny, (hggzz_full->GetBinContent(2, biny))*(yield_offshell_ggzz[EnergyIndex][folder] / ggzz_offshell_sum / luminosity[EnergyIndex]));

    hggzz_full->SetBinError(1, biny, (hggzz_full->GetBinError(1, biny))*(yield_80_100_ggzz[EnergyIndex][folder] / ggzz_80_100_sum / luminosity[EnergyIndex]));
    hggzz_full->SetBinError(4, biny, (hggzz_full->GetBinError(4, biny))*(yield_signal_ggzz[EnergyIndex][folder] / ggzz_signal_sum / luminosity[EnergyIndex]));
    hggzz_full->SetBinError(3, biny, (hggzz_full->GetBinError(3, biny))*(yield_offshell_ggzz[EnergyIndex][folder] / ggzz_offshell_sum / luminosity[EnergyIndex]));
    hggzz_full->SetBinError(2, biny, (hggzz_full->GetBinError(2, biny))*(yield_offshell_ggzz[EnergyIndex][folder] / ggzz_offshell_sum / luminosity[EnergyIndex]));
  }
  for (int biny = 0; biny <= hzx_full->GetNbinsY()+1; biny++){
    hzx_full->SetBinContent(1, biny, (hzx_full->GetBinContent(1, biny))*(yield_80_100_zx[EnergyIndex][folder] / zx_80_100_sum / luminosity[EnergyIndex]));
    hzx_full->SetBinContent(4, biny, (hzx_full->GetBinContent(4, biny))*(yield_signal_zx[EnergyIndex][folder] / zx_signal_sum / luminosity[EnergyIndex]));
    hzx_full->SetBinContent(3, biny, (hzx_full->GetBinContent(3, biny))*(yield_offshell_zx[EnergyIndex][folder] / zx_offshell_sum / luminosity[EnergyIndex]));
    hzx_full->SetBinContent(2, biny, (hzx_full->GetBinContent(2, biny))*(yield_offshell_zx[EnergyIndex][folder] / zx_offshell_sum / luminosity[EnergyIndex]));

    hzx_full->SetBinError(1, biny, (hzx_full->GetBinError(1, biny))*(yield_80_100_zx[EnergyIndex][folder] / zx_80_100_sum / luminosity[EnergyIndex]));
    hzx_full->SetBinError(4, biny, (hzx_full->GetBinError(4, biny))*(yield_signal_zx[EnergyIndex][folder] / zx_signal_sum / luminosity[EnergyIndex]));
    hzx_full->SetBinError(3, biny, (hzx_full->GetBinError(3, biny))*(yield_offshell_zx[EnergyIndex][folder] / zx_offshell_sum / luminosity[EnergyIndex]));
    hzx_full->SetBinError(2, biny, (hzx_full->GetBinError(2, biny))*(yield_offshell_zx[EnergyIndex][folder] / zx_offshell_sum / luminosity[EnergyIndex]));
  }
  for (int biny = 0; biny <= hqqzz_full->GetNbinsY()+1; biny++){
    hqqzz_full->SetBinContent(4, biny, (hqqzz_full->GetBinContent(4, biny))*(yield_signal_qqzz[EnergyIndex][folder] / qqzz_signal_sum / luminosity[EnergyIndex]));
    hqqzz_full->SetBinError(4, biny, (hqqzz_full->GetBinError(4, biny))*(yield_signal_qqzz[EnergyIndex][folder] / qqzz_signal_sum / luminosity[EnergyIndex]));
    hqqzz_full->SetBinContent(4, biny, (hqqzz_full->GetBinContent(4, biny))*(ratio_targetsignalyield[4][4][3]/ratio_targetsignalyield[4][0][3]));
    hqqzz_full->SetBinError(4, biny, (hqqzz_full->GetBinError(4, biny))*(ratio_targetsignalyield[4][4][3]/ratio_targetsignalyield[4][0][3]));

    for (int hh=0; hh<kGGHSamples; hh++){
      hhiggs[hh]->SetBinContent(biny, (hhiggs[hh]->GetBinContent(biny))*(yield_signal_higgs[EnergyIndex][folder]*(ratio_targetsignalyield[0][1][3] / ratio_targetsignalyield[0][0][3]) / higgs_signal_sum[hh] / luminosity[EnergyIndex]));
      hhiggs[hh]->SetBinError(biny, (hhiggs[hh]->GetBinError(biny))*(yield_signal_higgs[EnergyIndex][folder]*(ratio_targetsignalyield[0][1][3] / ratio_targetsignalyield[0][0][3]) / higgs_signal_sum[hh] / luminosity[EnergyIndex]));
    }
  }
  for (int hh = 0; hh < kGGHSamples; hh++) hhiggs[hh]->Scale(ratio_targetsignalyield[hh][4][3] / ratio_targetsignalyield[hh][1][3]);
  for (int binx = 1; binx <= hggzz_full->GetNbinsX(); binx++){
    double scale = ratio_targetsignalyield[5][4][binx-1]/ratio_targetsignalyield[5][0][binx-1];
    for (int biny = 0; biny <= hggzz_full->GetNbinsY()+1; biny++){
      double bincontent = hggzz_full->GetBinContent(binx, biny);
      double binerror = hggzz_full->GetBinError(binx, biny);

      bincontent *= scale;
      binerror *= scale;
      hggzz_full->SetBinContent(binx, biny, bincontent);
      hggzz_full->SetBinError(binx, biny, binerror);
    }
  }
  for (int binx = 1; binx <= hzx_full->GetNbinsX(); binx++){
    double scale = ratio_targetsignalyield[6][4][binx-1]/ratio_targetsignalyield[6][0][binx-1];
    for (int biny = 0; biny <= hzx_full->GetNbinsY()+1; biny++){
      double bincontent = hzx_full->GetBinContent(binx, biny);
      double binerror = hzx_full->GetBinError(binx, biny);

      bincontent *= scale;
      binerror *= scale;
      hzx_full->SetBinContent(binx, biny, bincontent);
      hzx_full->SetBinError(binx, biny, binerror);
    }
  }

  cout << "Signal gg norm: " << hggzz_full->Integral(4, 4, 0, hggzz_full->GetNbinsY()+1)*luminosity[EnergyIndex] << "\tqq norm: " << hqqzz_full->Integral(4, 4, 0, hqqzz_full->GetNbinsY()+1)*luminosity[EnergyIndex]
    << "\tZX norm: " << hzx_full->Integral(4, 4, 0, hzx_full->GetNbinsY()+1)*luminosity[EnergyIndex] << "\tHiggs norm: " << hhiggs[0]->Integral(0, hhiggs[0]->GetNbinsX()+1)*luminosity[EnergyIndex] << "\tHiggs (BSM) norm: " << hhiggs[3]->Integral(0, hhiggs[3]->GetNbinsX()+1)*luminosity[EnergyIndex] << endl;

  for (int sb=1; sb<=3; sb++) cout << "SB" << sb << " gg norm: " << hggzz_full->Integral(sb, sb, 0, hggzz_full->GetNbinsY()+1)*luminosity[EnergyIndex]
    << "\tqq norm: " << hqqzz_full->Integral(sb, sb, 0, hqqzz_full->GetNbinsY()+1)*luminosity[EnergyIndex]
    << "\tZX norm: " << hzx_full->Integral(sb, sb, 0, hzx_full->GetNbinsY()+1)*luminosity[EnergyIndex] << endl;

  foutput->cd();
  TTree* tdata_zerosb = new TTree("SB13_Data", "");
  TTree* tbkg_zerosb = new TTree("SB13_MC", "");
  TTree* tbkg_zxsignal = new TTree("SR_ZX", "");
  TTree* tbkg_signal = new TTree("SR_NonZX", "");
  TTree* tbkg_zxsb = new TTree("SB2_ZX", "");
  TTree* tSystContainer[5]={
    tdata_zerosb, tbkg_zerosb, tbkg_signal, tbkg_zxsignal, tbkg_zxsb
  };
  float fullweight=0;
  for (int tt=0; tt<5; tt++){
    tSystContainer[tt]->Branch("weight", &fullweight);

    tSystContainer[tt]->Branch("ZZMass", &ZZMass);
    if (tt==3 || tt==4) tSystContainer[tt]->Branch("CRflag", &CRflag);

    tSystContainer[tt]->Branch(cSystVariable[0], &Txy);
    tSystContainer[tt]->Branch(cSystVariable[1], &D_bkg);
    if(iSyst>=2) tSystContainer[tt]->Branch(cSystVariable[2], &Txy_BS);
  }
  double sum_sb13_mc=0, sum_sr_nonzx=0, sum_sb2_zx=0, sum_signal_zx=0;
  for (int tt = 0; tt < ntrees; tt++){
    for (int ev = 0; ev < tc[tt]->GetEntries(); ev++){
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

      tc[tt]->GetEntry(ev);
      if ((fabs(Lep1_Z1SIP) >= 4 || fabs(Lep2_Z1SIP) >= 4 || fabs(Lep3_Z1SIP) >= 5 || fabs(Lep4_Z1SIP) >= 5 || KalmanCandVtx_chi2 >= 30)) continue;
      if (
        tt==3 &&
        !applyZXloose(CRflag,
        Lep1combRelIsoPF,
        Lep2combRelIsoPF,
        Lep3combRelIsoPF,
        Lep4combRelIsoPF,
        Lep3isID,
        Lep4isID)
        ) continue;
      if (tt==3 && !ZXchannelselection(folder, Z1ids, CRflag)) continue;
      if (tt==3 && (
        (CRflag==6 || CRflag==10) ||
        (CRflag==8 || CRflag==12)
        )){
        int regionIndex = 3;
        for (int ir=0; ir<4; ir++){
          if (ZZMass<m4l_lowhigh[ir][1] && ZZMass>=m4l_lowhigh[ir][0]){
            regionIndex = ir; break;
          }
        }
        if (ZZMass<m4l_lowhigh[0][0]) regionIndex = 0;
        if (ZZMass>=m4l_lowhigh[2][1]) regionIndex = 2;
        if (tc[tt]->GetBranchStatus("ZXfake_weight_OS")) MC_weight = ZXfake_weight_OS[4] * (targetCRcount[0][4][regionIndex] / (targetCRcount[0][4][regionIndex] + targetCRcount[1][4][regionIndex]));
      }
      else if (tt==3 && (
        (CRflag==7 || CRflag==11) ||
        (CRflag==5 || CRflag==9)
        )){
        if (tc[tt]->GetBranchStatus("ZXfake_weight_SS")){
          int OSoverSS_biny = -1;
          if ((CRflag==5 || CRflag==9) && Z1ids==-169) OSoverSS_biny = 1; // 4mu
          if ((CRflag==7 || CRflag==11) && Z1ids==-169) OSoverSS_biny = 3; // 2mu2e
          if ((CRflag==5 || CRflag==9) && Z1ids==-121) OSoverSS_biny = 4; // 2e2mu
          if ((CRflag==7 || CRflag==11) && Z1ids==-121) OSoverSS_biny = 6; // 4e
          float SSwgt_OSoverSS = 1;
          if (OSoverSS_biny>0) SSwgt_OSoverSS = hRatio->GetBinContent(hRatio->GetXaxis()->FindBin(ZZMass), OSoverSS_biny);

          int regionIndex = 3;
          for (int ir=0; ir<4; ir++){
            if (ZZMass<m4l_lowhigh[ir][1] && ZZMass>=m4l_lowhigh[ir][0]){
              regionIndex = ir; break;
            }
          }
          if (ZZMass<m4l_lowhigh[0][0]) regionIndex = 0;
          if (ZZMass>=m4l_lowhigh[2][1]) regionIndex = 2;

          MC_weight = ZXfake_weight_SS*SSwgt_OSoverSS;
          MC_weight *= (targetCRcount[1][4][regionIndex] / (targetCRcount[0][4][regionIndex] + targetCRcount[1][4][regionIndex]));
        }
      }

      if (iSyst>=2 && iSyst<nSystVars){ // BE CAREFUL HERE, HARDCODED BOUNDS AHEAD
        int index_BS = 1;
        bool matchFound = false;
        if (tt==0 || tt==3) index_BS = 0;
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
      }

      CandVtx_x = KalmanCandVtx_x - OfflinePrimaryVtx_x;
      CandVtx_y = KalmanCandVtx_y - OfflinePrimaryVtx_y;
      CandVtx_z = KalmanCandVtx_z - OfflinePrimaryVtx_z;

      TrueCandVtx_x = GenIntVtx_x - GenPrimaryVtx_x;
      TrueCandVtx_y = GenIntVtx_y - GenPrimaryVtx_y;
      TrueCandVtx_z = GenIntVtx_z - GenPrimaryVtx_z;

      Dxy = (CandVtx_x*cos(ZZPhi) + CandVtx_y*sin(ZZPhi));
      Dxy_true = (TrueCandVtx_x*cos(GenHPhi) + TrueCandVtx_y*sin(GenHPhi));
      Txy = Dxy*ZZMass / ZZPt;
      Txy_true = Dxy_true*GenHMass / GenHPt;

      delDxy = compute_Dxy(
        ZZPhi,
        KalmanCandVtx_cov_xx,
        KalmanCandVtx_cov_xy,
        KalmanCandVtx_cov_xz,
        KalmanCandVtx_cov_yy,
        KalmanCandVtx_cov_yz,
        KalmanCandVtx_cov_zz,
        OfflinePrimaryVtx_cov_xx,
        OfflinePrimaryVtx_cov_xy,
        OfflinePrimaryVtx_cov_xz,
        OfflinePrimaryVtx_cov_yy,
        OfflinePrimaryVtx_cov_yz,
        OfflinePrimaryVtx_cov_zz
        );
      sigmaPV_xy = compute_sigmaxy(
        ZZPhi,
        OfflinePrimaryVtx_cov_xx,
        OfflinePrimaryVtx_cov_xy,
        OfflinePrimaryVtx_cov_xz,
        OfflinePrimaryVtx_cov_yy,
        OfflinePrimaryVtx_cov_yz,
        OfflinePrimaryVtx_cov_zz
        );
      sigmaInt_xy = compute_sigmaxy(
        ZZPhi,
        KalmanCandVtx_cov_xx,
        KalmanCandVtx_cov_xy,
        KalmanCandVtx_cov_xz,
        KalmanCandVtx_cov_yy,
        KalmanCandVtx_cov_yz,
        KalmanCandVtx_cov_zz
        );
      delTxy = delDxy*ZZMass / ZZPt;

      if (!tc[tt]->GetBranchStatus("D_bkg")) D_bkg = p0plus_VAJHU*p0plus_m4l / (p0plus_VAJHU*p0plus_m4l + bkg_VAMCFM*bkg_m4l);

      Txy*=10000; // in 1 um
      Txy_BS*=10000; // in 1 um

      if (iSyst==0) KD = Txy;
      if (iSyst==2) KD = Txy_BS;
      if (iSyst==1){
        KD = D_bkg; if (KD==1.) KD=0.999999;
      }

      float wgt = 1;
      for (int binx = 0; binx < nMZZ_ranges-1; binx++){
        if (tt == 2 && binx == 0) continue;
        if (tt != 3 && binx == 1) continue;

        if (tt == 3 && binx == 0 && ZZMass >= systZZMass_range[binx][0] + 30 && ZZMass < systZZMass_range[binx][1] + 30){ // ZX special fill
          wgt = MC_weight;
          fullweight = wgt * (yield_80_100_zx[EnergyIndex][folder] / zx_80_100_sum)*(ratio_targetsignalyield[6][4][binx]/ratio_targetsignalyield[6][0][binx]);
          ZZMass -= 30.;
          tbkg_zerosb->Fill();
          sum_sb13_mc += fullweight;
          ZZMass += 30.;
          continue;
        }
        else if (tt == 3 && binx == 0) continue;

        if (ZZMass >= systZZMass_range[binx][1] || ZZMass < systZZMass_range[binx][0]) continue;

        if (tt == 1 && binx == 0){ // gg bkg special fill from qqZZ
          wgt = MC_weight*MC_weight_Kfactor*MC_weight_QQBGGProper[0];
          fullweight = wgt * (yield_80_100_ggzz[EnergyIndex][folder] / ggzz_80_100_sum) * (ratio_targetsignalyield[5][4][binx]/ratio_targetsignalyield[5][0][binx]);
          wgt = MC_weight*MC_weight_QQZZEWK;
          fullweight += wgt * luminosity[EnergyIndex];
          tbkg_zerosb->Fill();
          sum_sb13_mc += fullweight;
        }
        else{
          wgt = MC_weight;
          if (tt==0){
            fullweight = wgt;
            tdata_zerosb->Fill();
          }
          else if (tt == 1){
            wgt *= MC_weight_QQZZEWK;
            fullweight = wgt * luminosity[EnergyIndex];
            tbkg_zerosb->Fill();
            sum_sb13_mc += fullweight;
          }
          else if (tt == 2){
            wgt *= MC_weight_Kfactor;
            fullweight = wgt * (yield_offshell_ggzz[EnergyIndex][folder] / ggzz_offshell_sum) * (ratio_targetsignalyield[5][4][binx]/ratio_targetsignalyield[5][0][binx]);
            tbkg_zerosb->Fill();
            sum_sb13_mc += fullweight;
          }
          else if (tt == 3){
            if (binx==2){
              fullweight = wgt * (yield_offshell_zx[EnergyIndex][folder] / zx_offshell_sum) * (ratio_targetsignalyield[6][4][binx]/ratio_targetsignalyield[6][0][binx]);
              tbkg_zerosb->Fill();
              sum_sb13_mc += fullweight;
            }
            else if (binx==1){
              fullweight = wgt * (yield_offshell_zx[EnergyIndex][folder] / zx_offshell_sum) * (ratio_targetsignalyield[6][4][binx]/ratio_targetsignalyield[6][0][binx]);
              tbkg_zxsb->Fill();
              sum_sb2_zx += fullweight;
            }
          }
        }
      }
      if (tt==3 && ZZMass >= systZZMass_range[nMZZ_ranges-1][0] && ZZMass < systZZMass_range[nMZZ_ranges-1][1]){ // ZX signal fill
        wgt = MC_weight;
        fullweight = wgt * (yield_signal_zx[EnergyIndex][folder] / zx_signal_sum) * (ratio_targetsignalyield[6][4][3]/ratio_targetsignalyield[6][0][3]);
        tbkg_zxsignal->Fill();
        sum_signal_zx += fullweight;
      }
      else if (tt==2 && ZZMass >= systZZMass_range[nMZZ_ranges-1][0] && ZZMass < systZZMass_range[nMZZ_ranges-1][1]){
        wgt = MC_weight_noxsec*MC_weight_Kfactor*MC_weight_QQBGGProper[0];
        fullweight = wgt * (yield_signal_ggzz[EnergyIndex][folder] / ggzz_signal_sum)*(ratio_targetsignalyield[5][4][3]/ratio_targetsignalyield[5][0][3]);
        sum_sr_nonzx += fullweight;
        tbkg_signal->Fill();
      }
      else if (tt==1 && ZZMass >= systZZMass_range[nMZZ_ranges-1][0] && ZZMass < systZZMass_range[nMZZ_ranges-1][1]){
        wgt = MC_weight_noxsec*MC_weight_Kfactor*MC_weight_QQBGGProper[0]*(yield_signal_ggzz[EnergyIndex][folder] / ggzz_signal_sum)*(ratio_targetsignalyield[5][4][3]/ratio_targetsignalyield[5][0][3]);
        fullweight = wgt;
        wgt = MC_weight_noxsec*MC_weight_QQBGGProper[1]*(yield_signal_qqzz[EnergyIndex][folder] / qqzz_signal_sum)*(ratio_targetsignalyield[4][4][3]/ratio_targetsignalyield[4][0][3]);
        fullweight += wgt;
        sum_sr_nonzx += fullweight;
        tbkg_signal->Fill();
      }
    }
  }

  for (int tt = 1; tt <= 2; tt++){
    for (int ev = 0; ev < tc_extra[tt]->GetEntries(); ev++){
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

      tc_extra[tt]->GetEntry(ev);
      if ((fabs(Lep1_Z1SIP) >= 4 || fabs(Lep2_Z1SIP) >= 4 || fabs(Lep3_Z1SIP) >= 5 || fabs(Lep4_Z1SIP) >= 5 || KalmanCandVtx_chi2 >= 30)) continue;

      int binx = nMZZ_ranges-1;
      if (ZZMass >= systZZMass_range[binx][1] || ZZMass < systZZMass_range[binx][0]) continue;

      if (iSyst>=2 && iSyst<nSystVars){ // BE CAREFUL HERE, HARDCODED BOUNDS AHEAD
        int index_BS = 1; // These are all MC.
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
      }

      CandVtx_x = KalmanCandVtx_x - OfflinePrimaryVtx_x;
      CandVtx_y = KalmanCandVtx_y - OfflinePrimaryVtx_y;
      CandVtx_z = KalmanCandVtx_z - OfflinePrimaryVtx_z;

      TrueCandVtx_x = GenIntVtx_x - GenPrimaryVtx_x;
      TrueCandVtx_y = GenIntVtx_y - GenPrimaryVtx_y;
      TrueCandVtx_z = GenIntVtx_z - GenPrimaryVtx_z;

      Dxy = (CandVtx_x*cos(ZZPhi) + CandVtx_y*sin(ZZPhi));
      Dxy_true = (TrueCandVtx_x*cos(GenHPhi) + TrueCandVtx_y*sin(GenHPhi));
      //			if(Dxy_true<0) cout << "Warning! Dxy_true = " << Dxy_true << " x/y: " << TrueCandVtx_x/TrueCandVtx_y << " cot(phi): " << cos(GenHPhi)/sin(GenHPhi) << endl;
      Txy = Dxy*ZZMass / ZZPt;
      Txy_true = Dxy_true*GenHMass / GenHPt;

      delDxy = compute_Dxy(
        ZZPhi,
        KalmanCandVtx_cov_xx,
        KalmanCandVtx_cov_xy,
        KalmanCandVtx_cov_xz,
        KalmanCandVtx_cov_yy,
        KalmanCandVtx_cov_yz,
        KalmanCandVtx_cov_zz,
        OfflinePrimaryVtx_cov_xx,
        OfflinePrimaryVtx_cov_xy,
        OfflinePrimaryVtx_cov_xz,
        OfflinePrimaryVtx_cov_yy,
        OfflinePrimaryVtx_cov_yz,
        OfflinePrimaryVtx_cov_zz
        );
      sigmaPV_xy = compute_sigmaxy(
        ZZPhi,
        OfflinePrimaryVtx_cov_xx,
        OfflinePrimaryVtx_cov_xy,
        OfflinePrimaryVtx_cov_xz,
        OfflinePrimaryVtx_cov_yy,
        OfflinePrimaryVtx_cov_yz,
        OfflinePrimaryVtx_cov_zz
        );
      sigmaInt_xy = compute_sigmaxy(
        ZZPhi,
        KalmanCandVtx_cov_xx,
        KalmanCandVtx_cov_xy,
        KalmanCandVtx_cov_xz,
        KalmanCandVtx_cov_yy,
        KalmanCandVtx_cov_yz,
        KalmanCandVtx_cov_zz
        );
      delTxy = delDxy*ZZMass / ZZPt;

      if (!tc_extra[tt]->GetBranchStatus("D_bkg")) D_bkg = p0plus_VAJHU*p0plus_m4l / (p0plus_VAJHU*p0plus_m4l + bkg_VAMCFM*bkg_m4l);

      Txy*=10000; // in 1 um
      Txy_BS*=10000; // in 1 um

      if (iSyst==0) KD = Txy;
      if (iSyst==2) KD = Txy_BS;
      if (iSyst==1){
        KD = D_bkg; if (KD==1.) KD=0.999999;
      }

      float wgt = 1;
      if (tt == 1){
        wgt = MC_weight_noxsec*MC_weight_Kfactor*MC_weight_QQBGGProper[0]*(yield_signal_ggzz[EnergyIndex][folder] / ggzz_signal_sum)*(ratio_targetsignalyield[5][4][3]/ratio_targetsignalyield[5][0][3]);
        fullweight = wgt;
        wgt = MC_weight_noxsec*MC_weight_QQBGGProper[1]*(yield_signal_qqzz[EnergyIndex][folder] / qqzz_signal_sum)*(ratio_targetsignalyield[4][4][3]/ratio_targetsignalyield[4][0][3]);
        fullweight += wgt;
        sum_sr_nonzx += fullweight;
        tbkg_signal->Fill();
      }
      else if (tt == 2){
        wgt = MC_weight_noxsec*MC_weight_Kfactor*MC_weight_QQBGGProper[0];
        fullweight = wgt * (yield_signal_ggzz[EnergyIndex][folder] / ggzz_signal_sum)*(ratio_targetsignalyield[5][4][3]/ratio_targetsignalyield[5][0][3]);
        sum_sr_nonzx += fullweight;
        tbkg_signal->Fill();
      }
    }
  }

  cout << "Signal ZX: " << sum_signal_zx << '\t' << "Bkg ZX: " << sum_sb2_zx << '\t' << "SB13 MC: " << sum_sb13_mc << '\t' << "SR Non-ZX: " << sum_sr_nonzx << endl;
  for (int tt=0; tt<5; tt++){
    foutput->WriteTObject(tSystContainer[tt]);
    delete tSystContainer[tt];
  }

  for (int tt = 1; tt < ntrees; tt++){
    if (tt==3 && iSyst>=7 && iSyst<=11) continue;
    hfull[ntrees]->Add(hfull[tt]);
  }
  for (int tt = 0; tt < ntrees+1; tt++){
    if (tt == ntrees){
      for (int hh = 0; hh < kGGHSamples; hh++) foutput->WriteTObject(hhiggs[hh]);
    }
    foutput->WriteTObject(hfull[tt]);

    TH1F* hProjKD = (TH1F*)hfull[tt]->ProjectionY(Form("%s_py", hfull[tt]->GetName()), 0, -1, "e");
    hProjKD->SetTitle(Form("%s Projection", hfull[tt]->GetTitle()));
    hproj[tt] = hProjKD;

    foutput->WriteTObject(hproj[tt]);

    for (int binx = 1; binx <= hfull[tt]->GetNbinsX(); binx++){
      TH1F* hSlice = (TH1F*)hProjKD->Clone((binx<hfull[tt]->GetNbinsX() ? Form("%s_SidebandRegion%i", hproj[tt]->GetName(), binx) : Form("%s_SignalRegion", hproj[tt]->GetName())));
      for (int biny = 1; biny <= hfull[tt]->GetNbinsY(); biny++){
        hSlice->SetBinContent(biny, hfull[tt]->GetBinContent(binx, biny)); hSlice->SetBinError(biny, hfull[tt]->GetBinError(binx, biny));
      }
      hSlice->SetTitle(Form("%s %s GeV Projection", hfull[tt]->GetTitle(), hfull[tt]->GetXaxis()->GetBinLabel(binx)));
      foutput->WriteTObject(hSlice);
      delete hSlice;
    }
  }
  for (int tt = 0; tt < ntrees_extra; tt++) delete tc_extra[tt];
  for (int tt = 0; tt < ntrees+1; tt++){
    delete hproj[tt];
    delete hfull[tt];
    if (tt<ntrees) delete tc[tt];
  }
  for (int bs=0; bs<2; bs++) delete tBeam[bs];

  foutput->Close();
  delete hRatio;
  fSSinput->Close();
}

void fitLifetimeSystTrees(int folder, int iSyst=0){
  if (iSyst>=nSystVars) return;

  RooFit::Verbose(false);
  RooFit::PrintLevel(-1000);
  RooMsgService::instance().setStreamStatus(1, false);
  RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);

  string folder_id;
  if (folder==0) folder_id = "4#mu";
  if (folder==1) folder_id = "4e";
  if (folder==2) folder_id = "2e+2#mu";

  TString OUTPUT_NAME = "LifetimeKD_SystFits_DataMC_";
  OUTPUT_NAME.Append("newSIP_");
  OUTPUT_NAME.Append(user_folder[folder]);
  OUTPUT_NAME.Append(Form("_%s_", cSystVariable[iSyst]));
  OUTPUT_NAME.Append("AllTeV");
  OUTPUT_NAME.Append(".root");
  TString coutput_common = user_dir_hep + "Analysis/ShapeSystematics/";
  coutput_common.Append(Form("%s/", cSystVariable[iSyst]));
  TString coutput = coutput_common + OUTPUT_NAME;
  TString plotDir = coutput_common + "Plots/";
  TString mkdirCommand = "mkdir -p ";
  mkdirCommand.Append(plotDir);
  gSystem->Exec(mkdirCommand);
  TFile* foutput = new TFile(coutput, "recreate");

  TFile* finput[2];
  TTree* tSystContainer[5][2];
  TH1F* hzx=0;
  TH1F* hzx_sb=0;
  TH1F* hmc=0;
  TH1F* hdata=0;
  TH1F* hmc_sb=0;
  cout << "Starting to loop over input files... ";
  for (int erg_tev=7; erg_tev<9; erg_tev++){
    TString erg_dir;
    erg_dir.Form("LHC_%iTeV", erg_tev);
    TString comstring;
    comstring.Form("%iTeV", erg_tev);
    cout << erg_tev << " TeV ";
    int EnergyIndex=0; if (erg_tev==8) EnergyIndex=1;
    TString INPUT_NAME = "LifetimeKD_SystTrees_DataMC_";
    INPUT_NAME.Append("newSIP_");
    INPUT_NAME.Append(user_folder[folder]);
    INPUT_NAME.Append(Form("_%s_", cSystVariable[iSyst]));
    INPUT_NAME.Append(comstring);
    INPUT_NAME.Append(".root");
    TString cinput_common = user_dir_hep + "Analysis/ShapeSystematics/";
    cinput_common.Append(Form("%s/", cSystVariable[iSyst]));
    TString cinput = cinput_common + INPUT_NAME;
    finput[EnergyIndex] = new TFile(cinput, "read");

    if(finput[EnergyIndex]!=0 && !finput[EnergyIndex]->IsZombie()) cout << "... input file open... ";
    else cout << "... file not open...";

    TH1F* htemp = (TH1F*)finput[EnergyIndex]->Get("hzx_full_py_SignalRegion");
    htemp->Scale(luminosity[EnergyIndex]);
    if (hzx!=0) hzx->Add(htemp);
    else hzx = (TH1F*)htemp->Clone(Form("%s_clone", htemp->GetName()));
    delete htemp;
    htemp = (TH1F*)finput[EnergyIndex]->Get("hzx_full_py_SidebandRegion2");
    htemp->Scale(luminosity[EnergyIndex]);
    if (hzx_sb!=0) hzx_sb->Add(htemp);
    else hzx_sb = (TH1F*)htemp->Clone(Form("%s_clone", htemp->GetName()));
    delete htemp;

    cout << "... Z+X... ";

    htemp = (TH1F*)finput[EnergyIndex]->Get("hdata_full_py_SidebandRegion1");
    if (hdata!=0) hdata->Add(htemp);
    else hdata = (TH1F*)htemp->Clone(Form("%s_clone", htemp->GetName()));
    delete htemp;
    htemp = (TH1F*)finput[EnergyIndex]->Get("hdata_full_py_SidebandRegion3");
    if (hdata!=0) hdata->Add(htemp);
    else hdata = (TH1F*)htemp->Clone(Form("%s_clone", htemp->GetName()));
    delete htemp;

    cout << "... SB data... ";

    htemp = (TH1F*)finput[EnergyIndex]->Get("hmc_full_py_SidebandRegion1");
    htemp->Scale(luminosity[EnergyIndex]);
    if (hmc_sb!=0) hmc_sb->Add(htemp);
    else hmc_sb = (TH1F*)htemp->Clone(Form("%s_clone", htemp->GetName()));
    delete htemp;
    htemp = (TH1F*)finput[EnergyIndex]->Get("hmc_full_py_SidebandRegion3");
    htemp->Scale(luminosity[EnergyIndex]);
    if (hmc_sb!=0) hmc_sb->Add(htemp);
    else hmc_sb = (TH1F*)htemp->Clone(Form("%s_clone", htemp->GetName()));
    delete htemp;
    cout << "... SB MC... ";
    htemp = (TH1F*)finput[EnergyIndex]->Get("hmc_full_py_SignalRegion");
    htemp->Scale(luminosity[EnergyIndex]);
    if (hmc!=0) hmc->Add(htemp);
    else hmc = (TH1F*)htemp->Clone(Form("%s_clone", htemp->GetName()));
    delete htemp;
    cout << "... Signal MC... ";

    TTree* tdata_zerosb = (TTree*)finput[EnergyIndex]->Get("SB13_Data");
    TTree* tbkg_zerosb = (TTree*)finput[EnergyIndex]->Get("SB13_MC");
    TTree* tbkg_signal = (TTree*)finput[EnergyIndex]->Get("SR_NonZX");
    TTree* tbkg_zxsignal = (TTree*)finput[EnergyIndex]->Get("SR_ZX");
    TTree* tbkg_zxsb = (TTree*)finput[EnergyIndex]->Get("SB2_ZX");
    tSystContainer[0][EnergyIndex]=tbkg_zerosb;
    tSystContainer[1][EnergyIndex]=tdata_zerosb;
    tSystContainer[2][EnergyIndex]=tbkg_signal;
    tSystContainer[3][EnergyIndex]=tbkg_zxsignal;
    tSystContainer[4][EnergyIndex]=tbkg_zxsb;
    for (int tt=0; tt<5; tt++) tSystContainer[tt][EnergyIndex]->SetName(Form("%s_%i",tSystContainer[tt][EnergyIndex]->GetName(), EnergyIndex));
    cout << "done...";
  }
  hzx_sb->Scale(hzx->Integral()/hzx_sb->Integral());
  cout << endl;

  hdata->SetMarkerColor(kBlack);
  hdata->SetMarkerStyle(20);
  hdata->SetLineStyle(1);
  hdata->SetLineColor(kBlack);
  hdata->SetLineWidth(1);
  hdata->Sumw2();
  hmc->SetLineStyle(1);
  hmc->SetLineColor(kBlack);
  hmc->SetLineWidth(2);
  hmc_sb->SetLineStyle(1);
  hmc_sb->SetLineColor(kBlue);
  hmc_sb->SetLineWidth(2);

  TString OUTPUT_W_NAME = "LifetimeKD_SystFitsRecord_DataMC_";
  OUTPUT_W_NAME.Append("newSIP_");
  OUTPUT_W_NAME.Append(user_folder[folder]);
  OUTPUT_W_NAME.Append(Form("_%s_", cSystVariable[iSyst]));
  OUTPUT_W_NAME.Append("AllTeV");
  OUTPUT_W_NAME.Append(".root");
  TString coutput_w = coutput_common + OUTPUT_W_NAME;
  foutput->cd();
  RooWorkspace* w = new RooWorkspace("w", kTRUE);
  w->importClassCode(RooRealFlooredSumPdf::Class(), kTRUE);

  RooRealVar* KD = new RooRealVar(cSystVariable[iSyst], cSystVariable_label[iSyst], 0, -1.0e6, 1.0e6);
  RooRealVar* weight = new RooRealVar("weight", "weight", 0, 0, 1.0e10);

  const int nGaussians=6;
  RooRealVar* varA[nGaussians][2];
  RooFormulaVar* varF[nGaussians][2];
  RooRealVar* varB[nGaussians][2];
  RooFormulaVar* varS[nGaussians][2];
  RooRealVar* extra_smear = new RooRealVar("extra_smear", "extra_smear", 1, 0.0001, 10);
  RooRealVar* center = new RooRealVar("C", "C", 0);
  for (int gg=0; gg<nGaussians; gg++){
    for (int e=0; e<2; e++){
      if (gg==0){
        varA[gg][e] = new RooRealVar(Form("A%i_cat%i", gg, e), Form("A%i_cat%i", gg, e), 1);
        varF[gg][e] = new RooFormulaVar(Form("Norm%i_cat%i", gg, e), Form("Norm%i_cat%i", gg, e), "@0", RooArgList(*(varA[gg][e])));
        varB[gg][e] = new RooRealVar(Form("B%i_cat%i", gg, e), Form("B%i_cat%i", gg, e), 10, 5, 1000);
        varS[gg][e] = new RooFormulaVar(Form("Smear%i_cat%i", gg, e), Form("Smear%i_cat%i", gg, e), "@0*@1", RooArgList(*(varB[gg][e]),*(extra_smear)));
      }
      else{
        varA[gg][e] = new RooRealVar(Form("A%i_cat%i", gg, e), Form("A%i_cat%i", gg, e), 0.1, 0, 1);
        varF[gg][e] = new RooFormulaVar(Form("Norm%i_cat%i", gg, e), Form("Norm%i_cat%i", gg, e), "@0*@1", RooArgList(*(varF[gg-1][e]), *(varA[gg][e])));
        varB[gg][e] = new RooRealVar(Form("B%i_cat%i", gg, e), Form("B%i_cat%i", gg, e), 4,1,10);
        varS[gg][e] = new RooFormulaVar(Form("Smear%i_cat%i", gg, e), Form("Smear%i_cat%i", gg, e), "@0*@1", RooArgList(*(varS[gg-1][e]), *(varB[gg][e])));
      }
    }
  }

  RooDataSet* fitdata;
  TH1F* hSmearingResults = new TH1F("hCommonSmearing", "", 2, 0, 2);

  TString strtout = Form("%slogFits_%s_%s.txt", coutput_common.Data(), cSystVariable[iSyst], user_folder[folder]);
  cout << "Opening " << strtout << endl;
  ofstream tout(strtout.Data(), ios::out);
  for (int tt=0; tt<5; tt++){
    if (tt==0){
      extra_smear->setVal(1);
      extra_smear->setConstant(true);
    }
    else if (tt>=1 && tt<=2){
      extra_smear->setConstant(false);
      extra_smear->setVal(1);
      for (int gg=0; gg<nGaussians; gg++){
        for (int e=0; e<2; e++){
          varA[gg][e]->setConstant(true);
          varB[gg][e]->setConstant(true);
        }
      }
    }

    RooArgSet* fitArgs = new RooArgSet(*KD, *weight);
    RooDataSet* tmpdata[2]={
      new RooDataSet(Form("%s_fitdata", tSystContainer[tt][0]->GetName()), Form("%s_0_fitdata", tSystContainer[tt][0]->GetName()), tSystContainer[tt][0], *fitArgs, 0, "weight"),
      new RooDataSet(Form("%s_fitdata", tSystContainer[tt][1]->GetName()), Form("%s_1_fitdata", tSystContainer[tt][1]->GetName()), tSystContainer[tt][1], *fitArgs, 0, "weight")
    };
    fitdata = tmpdata[0];
    fitdata->append(*tmpdata[1]);
    cout << "Constructed datasets" << endl;
    if (tt==2) fitdata->Print("v");

    RooBreitWigner* breit[2];
    RooFormulaVar* breit_width[2];
    RooGaussian* gaussians[nGaussians-2];
    RooRealFlooredSumPdf* tmp_pdf;
    breit_width[0] = new RooFormulaVar(Form("%s_bwWidth0", tSystContainer[tt][0]->GetName()), "", "@0*2", RooArgList(*(varS[0][0])));
    breit_width[1] = new RooFormulaVar(Form("%s_bwWidth1", tSystContainer[tt][0]->GetName()), "", "@0*2", RooArgList(*(varS[1][0])));
    breit[0] = new RooBreitWigner(Form("%s_bw_0", tSystContainer[tt][0]->GetName()), "", *KD, *center, *breit_width[0]);
    breit[1] = new RooBreitWigner(Form("%s_bw_1", tSystContainer[tt][0]->GetName()), "", *KD, *center, *breit_width[1]);

    for (int gg=0; gg<nGaussians; gg++){
      RooGaussian* my_g = new RooGaussian(Form("%s_gauss%i", tSystContainer[tt][0]->GetName(), gg), "", *KD, *center, *varS[gg][0]);
      gaussians[gg]=my_g;
    }
    cout << "Set the Gaussians" << endl;

    int nShapes = (folder!=0 ? 5 : 4);
    RooArgList funcList;
    RooArgList coeffList;
    funcList.add(*breit[0]);
    for (int gg=1; gg<nShapes; gg++){
      funcList.add(*(gaussians[gg]));
    }
    for (int gg=0; gg<nShapes; gg++){
        coeffList.add(*(varF[gg][0]));
    }
    tmp_pdf = new RooRealFlooredSumPdf(Form("%s_pdf", tSystContainer[tt][0]->GetName()), "", funcList, coeffList);

    RooRealFlooredSumPdf* sumpdf = tmp_pdf;
//    if (tt!=1 && tt<=2) sumpdf->fitTo(*fitdata, SumW2Error(kFALSE));
//    else if (tt<=2) sumpdf->fitTo(*fitdata, SumW2Error(kTRUE));
    if (tt<=2) sumpdf->fitTo(*fitdata, SumW2Error(kTRUE));
    w->import(*sumpdf);

    if (tt==1 || tt==2){
      hSmearingResults->SetBinContent(tt, extra_smear->getVal());
      hSmearingResults->SetBinError(tt, extra_smear->getError());
    }
    if (tt==0){
      cout << "\n\n";
      for (int ss=0; ss<nShapes; ss++){
        cout << "| "
          << "A[" << ss << "] = " << varA[ss][0]->getVal() << " +" << varA[ss][0]->getErrorHi() << "/" << varA[ss][0]->getErrorLo() << " | "
          << "B[" << ss << "] = " << varB[ss][0]->getVal() << " +" << varB[ss][0]->getErrorHi() << "/" << varB[ss][0]->getErrorLo() << " |" << endl;
        tout << "| "
          << "A[" << ss << "] = " << varA[ss][0]->getVal() << " +" << varA[ss][0]->getErrorHi() << "/" << varA[ss][0]->getErrorLo() << " | "
          << "B[" << ss << "] = " << varB[ss][0]->getVal() << " +" << varB[ss][0]->getErrorHi() << "/" << varB[ss][0]->getErrorLo() << " |" << endl;
      }
      cout << endl;
      tout << endl;
    }
    else if (tt==1 || tt==2){
      cout << "\n\nIteration " << tt << " Smear: "
        << extra_smear->getVal() << " +" << extra_smear->getErrorHi() << "/" << extra_smear->getErrorLo()
        << endl;
      tout << "Iteration " << tt << " Smear: "
        << extra_smear->getVal() << " +" << extra_smear->getErrorHi() << "/" << extra_smear->getErrorLo()
        << endl;
    }

    TString canvasname_2D = Form("cCanvas_FitComparisons_%s_%s", tSystContainer[tt][0]->GetName(), user_folder[folder]);
    canvasname_2D.Append("_AllTeV");
    TString yTitle = Form("Events / %.0f GeV", (visualRange_SystVariable[iSyst][1] - visualRange_SystVariable[iSyst][0]) / nbins_SystVariable[iSyst]);
    TString xTitle = cSystVariable_label[iSyst];

    if (tt==1){
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
      //      if (erg_tev==7) text = pt->AddText(0.837, 0.45, "#font[42]{ 5.1 fb^{-1} (7 TeV)}");
      //      if (erg_tev==8) text = pt->AddText(0.837, 0.45, "#font[42]{19.7 fb^{-1} (8 TeV)}");
      text = pt->AddText(0.537, 0.45, "#font[42]{19.7 fb^{-1} (8 TeV) + 5.1 fb^{-1} (7 TeV)}");
      text->SetTextSize(0.0315);

      TCanvas* c2D = new TCanvas(canvasname_2D, "", 8, 30, 800, 800);
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

      RooPlot* rplot = new RooPlot(Form("rplot_%s", tSystContainer[tt][0]->GetName()), "", *KD, visualRange_SystVariable[iSyst][0], visualRange_SystVariable[iSyst][1], nbins_SystVariable[iSyst]);
      fitdata->plotOn(rplot, MarkerColor(kBlack), MarkerStyle(20), MarkerSize(1.0), LineWidth(1), XErrorSize(0), DataError(RooAbsData::Poisson));
      sumpdf->plotOn(rplot, LineColor(kRed), LineWidth(2), LineStyle(7));
      TH1F* hpdf = (TH1F*)sumpdf->createHistogram(Form("%s_histo", sumpdf->GetName()), *KD, Binning(nbins_SystVariable[iSyst], visualRange_SystVariable[iSyst][0], visualRange_SystVariable[iSyst][1]));
      hpdf->Scale(hdata->Integral()/hpdf->Integral());

      TH1F* hmc_clone = (TH1F*)hmc->Clone(Form("%s_clone", hmc->GetName()));
      hmc_clone->Scale(hdata->Integral()/hmc_clone->Integral());
      TH1F* hmc_sb_clone = (TH1F*)hmc_sb->Clone(Form("%s_clone", hmc_sb->GetName()));
      hmc_sb_clone->Scale(hdata->Integral()/hmc_sb_clone->Integral());

      TLegend *l2D = new TLegend(0.20, 0.65, 0.58, 0.88);
      l2D->SetBorderSize(0);
      l2D->SetTextFont(42);
      l2D->SetTextSize(0.03);
      l2D->SetLineColor(1);
      l2D->SetLineStyle(1);
      l2D->SetLineWidth(1);
      l2D->SetFillColor(0);
      l2D->SetFillStyle(0);

      rplot->GetXaxis()->SetTitle(xTitle);
      rplot->GetYaxis()->SetTitle(yTitle);
      rplot->GetXaxis()->SetNdivisions(505);
      rplot->GetXaxis()->SetLabelFont(42);
      rplot->GetXaxis()->SetLabelOffset(0.007);
      rplot->GetXaxis()->SetLabelSize(0.04);
      rplot->GetXaxis()->SetTitleSize(0.06);
      rplot->GetXaxis()->SetTitleOffset(0.9);
      rplot->GetXaxis()->SetTitleFont(42);
      rplot->GetYaxis()->SetNdivisions(505);
      rplot->GetYaxis()->SetLabelFont(42);
      rplot->GetYaxis()->SetLabelOffset(0.007);
      rplot->GetYaxis()->SetLabelSize(0.04);
      rplot->GetYaxis()->SetTitleSize(0.06);
      rplot->GetYaxis()->SetTitleOffset(1.1);
      rplot->GetYaxis()->SetTitleFont(42);
      hpdf->SetLineStyle(1);
      hpdf->SetLineColor(kRed);
      hpdf->SetLineWidth(2);
      TH1F* hpdf_clone = (TH1F*)hpdf->Clone("hpdf_temp");
      hpdf_clone->SetLineStyle(7);

      double max_plot = max(hmc->GetMaximum(), max(hpdf->GetMaximum(), max(hmc_sb->GetMaximum(), hdata->GetMaximum())));
      rplot->GetYaxis()->SetRangeUser(0, max_plot*1.6);
      rplot->GetXaxis()->SetRangeUser((folder!=0 ? -762.5 : -312.5), (folder!=0 ? 762.5 : 312.5));


      l2D->AddEntry(hdata, "Observed sideband", "ep");
      l2D->AddEntry(hpdf_clone, "Fit to sideband data", "l");
//      l2D->AddEntry(hpdf, "Fit (histogram)", "l");
      l2D->AddEntry(hmc_sb_clone, "Bkg. in sideband", "l");

      rplot->Draw();
//      hpdf->Draw("histsame");
      hmc_sb_clone->Draw("histsame");
//      hmc_clone->Draw("histsame");
      l2D->Draw();
      pt->Draw();

      TString strSidebandRegion = "Sidebands 1 and 3";
      TPaveText *pt10 = new TPaveText(0.68, 0.84, 0.90, 0.92, "brNDC");
      pt10->SetBorderSize(0);
      pt10->SetTextAlign(12);
      pt10->SetTextSize(0.03);
      pt10->SetFillStyle(0);
      pt10->SetTextFont(42);
      TText* text10 = pt10->AddText(0.01, 0.01, strSidebandRegion);
      pt10->Draw();

      TString strSIPcut = "New cut";
      strSIPcut.Prepend(", ");
      strSIPcut.Prepend(folder_id.c_str());
      TPaveText* pt20;
      pt20 = new TPaveText(0.68, 0.80, 0.90, 0.84, "brNDC");
      pt20->SetBorderSize(0);
      pt20->SetTextAlign(12);
      pt20->SetTextSize(0.03);
      pt20->SetFillStyle(0);
      pt20->SetTextFont(42);
      TText* text20 = pt20->AddText(0.01, 0.01, strSIPcut);
      pt20->Draw();

      TString canvasname = canvasname_2D;
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
      c2D->SaveAs(canvasname_pdf);
      c2D->SaveAs(canvasname_eps);
      c2D->SaveAs(canvasname_png);
      c2D->SaveAs(canvasname_root);
      c2D->SaveAs(canvasname_c);

      foutput->WriteTObject(c2D);
      delete pt20;
      delete pt10;
      delete l2D;
      delete hmc_clone;
      delete hmc_sb_clone;
      delete hpdf_clone;
      delete hpdf;
      delete rplot;
      c2D->Close();
      delete pt;
    }
    else if (tt==0){
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
      //      if (erg_tev==7) text = pt->AddText(0.837, 0.45, "#font[42]{ 5.1 fb^{-1} (7 TeV)}");
      //      if (erg_tev==8) text = pt->AddText(0.837, 0.45, "#font[42]{19.7 fb^{-1} (8 TeV)}");
      text = pt->AddText(0.537, 0.45, "#font[42]{19.7 fb^{-1} (8 TeV) + 5.1 fb^{-1} (7 TeV)}");
      text->SetTextSize(0.0315);

      TCanvas* c2D = new TCanvas(canvasname_2D, "", 8, 30, 800, 800);
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

      RooPlot* rplot = new RooPlot(Form("rplot_%s", tSystContainer[tt][0]->GetName()), "", *KD, visualRange_SystVariable[iSyst][0], visualRange_SystVariable[iSyst][1], nbins_SystVariable[iSyst]);
      fitdata->plotOn(rplot, MarkerColor(kBlue), MarkerStyle(1), MarkerSize(1.0), LineWidth(2), XErrorSize(0), DataError(RooAbsData::None));
      sumpdf->plotOn(rplot, LineColor(kRed), LineWidth(2), LineStyle(7));
      TH1F* hpdf = (TH1F*)sumpdf->createHistogram(Form("%s_histo", sumpdf->GetName()), *KD, Binning(nbins_SystVariable[iSyst], visualRange_SystVariable[iSyst][0], visualRange_SystVariable[iSyst][1]));

      TH1F* hmc_sb_clone = (TH1F*)hmc_sb->Clone(Form("%s_clone", hmc_sb->GetName()));
      hpdf->Scale(hmc_sb_clone->Integral()/hpdf->Integral());
      TH1F* hmc_clone = (TH1F*)hmc->Clone(Form("%s_clone", hmc->GetName()));
      hmc_clone->Scale(hmc_sb_clone->Integral()/hmc_clone->Integral());

      TLegend *l2D = new TLegend(0.20, 0.65, 0.58, 0.88);
      l2D->SetBorderSize(0);
      l2D->SetTextFont(42);
      l2D->SetTextSize(0.03);
      l2D->SetLineColor(1);
      l2D->SetLineStyle(1);
      l2D->SetLineWidth(1);
      l2D->SetFillColor(0);
      l2D->SetFillStyle(0);

      rplot->GetXaxis()->SetTitle(xTitle);
      rplot->GetYaxis()->SetTitle(yTitle);
      rplot->GetXaxis()->SetNdivisions(505);
      rplot->GetXaxis()->SetLabelFont(42);
      rplot->GetXaxis()->SetLabelOffset(0.007);
      rplot->GetXaxis()->SetLabelSize(0.04);
      rplot->GetXaxis()->SetTitleSize(0.06);
      rplot->GetXaxis()->SetTitleOffset(0.9);
      rplot->GetXaxis()->SetTitleFont(42);
      rplot->GetYaxis()->SetNdivisions(505);
      rplot->GetYaxis()->SetLabelFont(42);
      rplot->GetYaxis()->SetLabelOffset(0.007);
      rplot->GetYaxis()->SetLabelSize(0.04);
      rplot->GetYaxis()->SetTitleSize(0.06);
      rplot->GetYaxis()->SetTitleOffset(1.1);
      rplot->GetYaxis()->SetTitleFont(42);
      hpdf->SetLineStyle(1);
      hpdf->SetLineColor(kRed);
      hpdf->SetLineWidth(2);
      TH1F* hpdf_clone = (TH1F*)hpdf->Clone("hpdf_temp");
      hpdf_clone->SetLineStyle(7);

//      hmc_sb->Scale(hpdf->Integral()/hmc_sb->Integral());

      double max_plot = max(max(hpdf->GetMaximum(), hmc_sb_clone->GetMaximum()), hmc_clone->GetMaximum());
      rplot->GetYaxis()->SetRangeUser(0, max_plot*1.6);
      rplot->GetXaxis()->SetRangeUser((folder!=0 ? -762.5 : -312.5), (folder!=0 ? 762.5 : 312.5));

      l2D->AddEntry(hmc_sb_clone, "Bkg. in sideband", "l");
      l2D->AddEntry(hpdf_clone, "Fit to sideband bkg.", "l");
//      l2D->AddEntry(hpdf, "Fit (histogram)", "l");
      l2D->AddEntry(hmc_clone, "Bkg. (shape) in signal", "l");

      rplot->Draw();
//      hpdf->Draw("histsame");
      hmc_clone->Draw("histsame");
      hmc_sb_clone->Draw("histsame");
      l2D->Draw();
      pt->Draw();

      TString strSidebandRegion = "Sidebands 1 and 3";
      TPaveText *pt10 = new TPaveText(0.68, 0.84, 0.90, 0.92, "brNDC");
      pt10->SetBorderSize(0);
      pt10->SetTextAlign(12);
      pt10->SetTextSize(0.03);
      pt10->SetFillStyle(0);
      pt10->SetTextFont(42);
      TText* text10 = pt10->AddText(0.01, 0.01, strSidebandRegion);
      pt10->Draw();

      TString strSIPcut = "New cut";
      strSIPcut.Prepend(", ");
      strSIPcut.Prepend(folder_id.c_str());
      TPaveText* pt20;
      pt20 = new TPaveText(0.68, 0.80, 0.90, 0.84, "brNDC");
      pt20->SetBorderSize(0);
      pt20->SetTextAlign(12);
      pt20->SetTextSize(0.03);
      pt20->SetFillStyle(0);
      pt20->SetTextFont(42);
      TText* text20 = pt20->AddText(0.01, 0.01, strSIPcut);
      pt20->Draw();

      TString canvasname = canvasname_2D;
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
      c2D->SaveAs(canvasname_pdf);
      c2D->SaveAs(canvasname_eps);
      c2D->SaveAs(canvasname_png);
      c2D->SaveAs(canvasname_root);
      c2D->SaveAs(canvasname_c);

      foutput->WriteTObject(c2D);
      delete pt20;
      delete pt10;
      delete l2D;
      delete hmc_clone;
      delete hmc_sb_clone;
      delete hpdf_clone;
      delete hpdf;
      delete rplot;
      c2D->Close();
      delete pt;
    }
    else if (tt==2){
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
      //      if (erg_tev==7) text = pt->AddText(0.837, 0.45, "#font[42]{ 5.1 fb^{-1} (7 TeV)}");
      //      if (erg_tev==8) text = pt->AddText(0.837, 0.45, "#font[42]{19.7 fb^{-1} (8 TeV)}");
      text = pt->AddText(0.537, 0.45, "#font[42]{19.7 fb^{-1} (8 TeV) + 5.1 fb^{-1} (7 TeV)}");
      text->SetTextSize(0.0315);

      TCanvas* c2D = new TCanvas(canvasname_2D, "", 8, 30, 800, 800);
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

      RooPlot* rplot = new RooPlot(Form("rplot_%s", tSystContainer[tt][0]->GetName()), "", *KD, visualRange_SystVariable[iSyst][0], visualRange_SystVariable[iSyst][1], nbins_SystVariable[iSyst]);
      fitdata->plotOn(rplot, MarkerColor(kBlack), MarkerStyle(1), MarkerSize(1.0), LineWidth(2), XErrorSize(0), DataError(RooAbsData::None));
      sumpdf->plotOn(rplot, LineColor(kRed), LineWidth(2), LineStyle(7));
      TH1F* hpdf = (TH1F*)sumpdf->createHistogram(Form("%s_histo", sumpdf->GetName()), *KD, Binning(nbins_SystVariable[iSyst], visualRange_SystVariable[iSyst][0], visualRange_SystVariable[iSyst][1]));

      TH1F* hmc_clone = (TH1F*)hmc->Clone(Form("%s_clone", hmc->GetName()));
      hmc_clone->Add(hzx, -1);
      TH1F* hmc_sb_clone = (TH1F*)hmc_sb->Clone(Form("%s_clone", hmc_sb->GetName()));
      hmc_sb_clone->Scale(hmc_clone->Integral()/hmc_sb_clone->Integral());
      hpdf->Scale(hmc_clone->Integral()/hpdf->Integral());

      TLegend *l2D = new TLegend(0.20, 0.65, 0.58, 0.88);
      l2D->SetBorderSize(0);
      l2D->SetTextFont(42);
      l2D->SetTextSize(0.03);
      l2D->SetLineColor(1);
      l2D->SetLineStyle(1);
      l2D->SetLineWidth(1);
      l2D->SetFillColor(0);
      l2D->SetFillStyle(0);

      rplot->GetXaxis()->SetTitle(xTitle);
      rplot->GetYaxis()->SetTitle(yTitle);
      rplot->GetXaxis()->SetNdivisions(505);
      rplot->GetXaxis()->SetLabelFont(42);
      rplot->GetXaxis()->SetLabelOffset(0.007);
      rplot->GetXaxis()->SetLabelSize(0.04);
      rplot->GetXaxis()->SetTitleSize(0.06);
      rplot->GetXaxis()->SetTitleOffset(0.9);
      rplot->GetXaxis()->SetTitleFont(42);
      rplot->GetYaxis()->SetNdivisions(505);
      rplot->GetYaxis()->SetLabelFont(42);
      rplot->GetYaxis()->SetLabelOffset(0.007);
      rplot->GetYaxis()->SetLabelSize(0.04);
      rplot->GetYaxis()->SetTitleSize(0.06);
      rplot->GetYaxis()->SetTitleOffset(1.1);
      rplot->GetYaxis()->SetTitleFont(42);
      hpdf->SetLineStyle(1);
      hpdf->SetLineColor(kRed);
      hpdf->SetLineWidth(2);
      TH1F* hpdf_clone = (TH1F*)hpdf->Clone("hpdf_temp");
      hpdf_clone->SetLineStyle(7);

      hmc_sb->Scale(hpdf->Integral()/hmc_sb->Integral());

      double max_plot = max(hpdf->GetMaximum(), hmc_sb->GetMaximum());
      rplot->GetYaxis()->SetRangeUser(0, max_plot*1.6);
      rplot->GetXaxis()->SetRangeUser((folder!=0 ? -762.5 : -312.5), (folder!=0 ? 762.5 : 312.5));

      l2D->AddEntry(hmc_clone, "Prompt bkg. in signal", "l");
      l2D->AddEntry(hpdf_clone, "Fit to bkg. in signal", "l");
//      l2D->AddEntry(hpdf, "Fit (histogram)", "l");
      l2D->AddEntry(hmc_sb_clone, "Bkg. (shape) in sideband", "l");

      rplot->Draw();
//      hpdf->Draw("histsame");
      hmc_clone->Draw("histsame");
      hmc_sb_clone->Draw("histsame");
      l2D->Draw();
      pt->Draw();

      TString strSidebandRegion = "Signal region";
      TPaveText *pt10 = new TPaveText(0.68, 0.84, 0.90, 0.92, "brNDC");
      pt10->SetBorderSize(0);
      pt10->SetTextAlign(12);
      pt10->SetTextSize(0.03);
      pt10->SetFillStyle(0);
      pt10->SetTextFont(42);
      TText* text10 = pt10->AddText(0.01, 0.01, strSidebandRegion);
      pt10->Draw();

      TString strSIPcut = "New cut";
      strSIPcut.Prepend(", ");
      strSIPcut.Prepend(folder_id.c_str());
      TPaveText* pt20;
      pt20 = new TPaveText(0.68, 0.80, 0.90, 0.84, "brNDC");
      pt20->SetBorderSize(0);
      pt20->SetTextAlign(12);
      pt20->SetTextSize(0.03);
      pt20->SetFillStyle(0);
      pt20->SetTextFont(42);
      TText* text20 = pt20->AddText(0.01, 0.01, strSIPcut);
      pt20->Draw();

      TString canvasname = canvasname_2D;
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
      c2D->SaveAs(canvasname_pdf);
      c2D->SaveAs(canvasname_eps);
      c2D->SaveAs(canvasname_png);
      c2D->SaveAs(canvasname_root);
      c2D->SaveAs(canvasname_c);

      foutput->WriteTObject(c2D);
      delete pt20;
      delete pt10;
      delete l2D;
      delete hmc_clone;
      delete hmc_sb_clone;
      delete hpdf_clone;
      delete hpdf;
      delete rplot;
      c2D->Close();
      delete pt;
    }
    else if (tt==4){
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
      //      if (erg_tev==7) text = pt->AddText(0.837, 0.45, "#font[42]{ 5.1 fb^{-1} (7 TeV)}");
      //      if (erg_tev==8) text = pt->AddText(0.837, 0.45, "#font[42]{19.7 fb^{-1} (8 TeV)}");
      text = pt->AddText(0.537, 0.45, "#font[42]{19.7 fb^{-1} (8 TeV) + 5.1 fb^{-1} (7 TeV)}");
      text->SetTextSize(0.0315);

      TCanvas* c2D = new TCanvas(canvasname_2D, "", 8, 30, 800, 800);
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

      TLegend *l2D = new TLegend(0.20, 0.65, 0.58, 0.88);
      l2D->SetBorderSize(0);
      l2D->SetTextFont(42);
      l2D->SetTextSize(0.03);
      l2D->SetLineColor(1);
      l2D->SetLineStyle(1);
      l2D->SetLineWidth(1);
      l2D->SetFillColor(0);
      l2D->SetFillStyle(0);

//      TString yTitle = Form("Events / %.1f", (visualRange_SystVariable[iSyst][1] - visualRange_SystVariable[iSyst][0]) / nbins_SystVariable[iSyst]);
//      TString xTitle = cSystVariable_label[iSyst];
      hzx->SetTitle("");
      hzx->GetXaxis()->SetTitle(xTitle);
      hzx->GetYaxis()->SetTitle(yTitle);
      hzx->GetXaxis()->SetNdivisions(505);
      hzx->GetXaxis()->SetLabelFont(42);
      hzx->GetXaxis()->SetLabelOffset(0.007);
      hzx->GetXaxis()->SetLabelSize(0.04);
      hzx->GetXaxis()->SetTitleSize(0.06);
      hzx->GetXaxis()->SetTitleOffset(0.9);
      hzx->GetXaxis()->SetTitleFont(42);
      hzx->GetYaxis()->SetNdivisions(505);
      hzx->GetYaxis()->SetLabelFont(42);
      hzx->GetYaxis()->SetLabelOffset(0.007);
      hzx->GetYaxis()->SetLabelSize(0.04);
      hzx->GetYaxis()->SetTitleSize(0.06);
      hzx->GetYaxis()->SetTitleOffset(1.1);
      hzx->GetYaxis()->SetTitleFont(42);
      hzx_sb->SetLineStyle(1);
      hzx_sb->SetLineColor(kRed);
      hzx_sb->SetLineWidth(2);
      hzx->SetLineStyle(1);
      hzx->SetLineColor(kBlack);
      hzx->SetLineWidth(2);
      TH1F* hzx_sb_clone = (TH1F*)hzx_sb->Clone(Form("%s_clone", hzx_sb->GetName()));
      TH1F* hzx_clone = (TH1F*)hzx->Clone(Form("%s_clone", hzx->GetName()));
      hzx_clone->SetFillStyle(3002);

      double max_plot = max(hzx->GetMaximum(), hzx_sb->GetMaximum());
      hzx->GetYaxis()->SetRangeUser(0, max_plot*1.5);
      hzx->GetXaxis()->SetRangeUser((folder!=0 ? -762.5 : -312.5), (folder!=0 ? 762.5 : 312.5));

      l2D->AddEntry(hzx, "Z+X in signal", "l");
      l2D->AddEntry(hzx_sb, "Z+X in sideband", "l");

      hzx->Draw("hist");
      hzx_clone->Draw("e2same");
      hzx_sb->Draw("histsame");
      l2D->Draw();
      pt->Draw();

      TString strSidebandRegion = "Z+X Sideband 2";
      TPaveText *pt10 = new TPaveText(0.68, 0.84, 0.90, 0.92, "brNDC");
      pt10->SetBorderSize(0);
      pt10->SetTextAlign(12);
      pt10->SetTextSize(0.03);
      pt10->SetFillStyle(0);
      pt10->SetTextFont(42);
      TText* text10 = pt10->AddText(0.01, 0.01, strSidebandRegion);
      pt10->Draw();

      TString strSIPcut = "New cut";
      strSIPcut.Prepend(", ");
      strSIPcut.Prepend(folder_id.c_str());
      TPaveText* pt20;
      pt20 = new TPaveText(0.68, 0.80, 0.90, 0.84, "brNDC");
      pt20->SetBorderSize(0);
      pt20->SetTextAlign(12);
      pt20->SetTextSize(0.03);
      pt20->SetFillStyle(0);
      pt20->SetTextFont(42);
      TText* text20 = pt20->AddText(0.01, 0.01, strSIPcut);
      pt20->Draw();

      TString canvasname = canvasname_2D;
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
      c2D->SaveAs(canvasname_pdf);
      c2D->SaveAs(canvasname_eps);
      c2D->SaveAs(canvasname_png);
      c2D->SaveAs(canvasname_root);
      c2D->SaveAs(canvasname_c);

      foutput->WriteTObject(c2D);
      delete pt20;
      delete pt10;
      delete l2D;
      delete hzx_clone;
      delete hzx_sb_clone;
      c2D->Close();
      delete pt;
    }

    delete tmp_pdf;
    delete breit[0];
    delete breit[1];
    for (int gg=0; gg<nGaussians; gg++) gaussians[gg];
    for (int e=0; e<2; e++) delete tmpdata[e];
  }
  w->writeToFile(coutput_w);
  foutput->WriteTObject(hSmearingResults);
  tout.close();
  delete hSmearingResults;
  delete center;
  delete extra_smear;
  for (int gg=0; gg<nGaussians; gg++){
    for (int e=0; e<2; e++){
      delete varA[gg][e];
      delete varF[gg][e];
      delete varB[gg][e];
      delete varS[gg][e];
    }
  }
  delete KD;
  delete weight;
  delete w;

  for (int e=0; e<2; e++) finput[e]->Close();
  foutput->Close();
}

/*
void fitLifetimeSystTrees(int folder, int iSyst=0){
  if (iSyst>=nSystVars) return;

  RooFit::Verbose(false);
  RooFit::PrintLevel(-1000);
  RooMsgService::instance().setStreamStatus(1, false);
  RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);

  string folder_id;
  if (folder==0) folder_id = "4#mu";
  if (folder==1) folder_id = "4e";
  if (folder==2) folder_id = "2e2#mu";

  TString OUTPUT_NAME = "LifetimeKD_SystFits_DataMC_";
  OUTPUT_NAME.Append("newSIP_");
  OUTPUT_NAME.Append(user_folder[folder]);
  OUTPUT_NAME.Append(Form("_%s_", cSystVariable[iSyst]));
  OUTPUT_NAME.Append("AllTeV");
  OUTPUT_NAME.Append(".root");
  TString coutput_common = user_dir_hep + "Analysis/ShapeSystematics/";
  coutput_common.Append(Form("%s/", cSystVariable[iSyst]));
  TString coutput = coutput_common + OUTPUT_NAME;
  TFile* foutput = new TFile(coutput, "recreate");

  TFile* finput[2];
  TTree* tSystContainer[5][2];
  TH1F* hzx=0;
  TH1F* hzx_sb=0;
  TH1F* hmc=0;
  TH1F* hdata=0;
  TH1F* hmc_sb=0;
  cout << "Starting to loop over input files... ";
  for (int erg_tev=7; erg_tev<9; erg_tev++){
    TString erg_dir;
    erg_dir.Form("LHC_%iTeV", erg_tev);
    TString comstring;
    comstring.Form("%iTeV", erg_tev);
    cout << erg_tev << " TeV ";
    int EnergyIndex=0; if (erg_tev==8) EnergyIndex=1;
    TString INPUT_NAME = "LifetimeKD_SystTrees_DataMC_";
    INPUT_NAME.Append("newSIP_");
    INPUT_NAME.Append(user_folder[folder]);
    INPUT_NAME.Append(Form("_%s_", cSystVariable[iSyst]));
    INPUT_NAME.Append(comstring);
    INPUT_NAME.Append(".root");
    TString cinput_common = user_dir_hep + "Analysis/ShapeSystematics/";
    cinput_common.Append(Form("%s/", cSystVariable[iSyst]));
    TString cinput = cinput_common + INPUT_NAME;
    finput[EnergyIndex] = new TFile(cinput, "read");

    if (finput[EnergyIndex]!=0 && !finput[EnergyIndex]->IsZombie()) cout << "... input file open... ";
    else cout << "... file not open...";

    TH1F* htemp = (TH1F*)finput[EnergyIndex]->Get("hzx_full_py_SignalRegion");
    if (hzx!=0) hzx->Add(htemp);
    else hzx = (TH1F*)htemp->Clone(Form("%s_clone", htemp->GetName()));
    delete htemp;
    htemp = (TH1F*)finput[EnergyIndex]->Get("hzx_full_py_SidebandRegion2");
    if (hzx_sb!=0) hzx_sb->Add(htemp);
    else hzx_sb = (TH1F*)htemp->Clone(Form("%s_clone", htemp->GetName()));
    delete htemp;

    cout << "... Z+X... ";

    htemp = (TH1F*)finput[EnergyIndex]->Get("hdata_full_py_SidebandRegion1");
    if (hdata!=0) hdata->Add(htemp);
    else hdata = (TH1F*)htemp->Clone(Form("%s_clone", htemp->GetName()));
    delete htemp;
    htemp = (TH1F*)finput[EnergyIndex]->Get("hdata_full_py_SidebandRegion3");
    if (hdata!=0) hdata->Add(htemp);
    else hdata = (TH1F*)htemp->Clone(Form("%s_clone", htemp->GetName()));
    delete htemp;

    cout << "... SB data... ";

    htemp = (TH1F*)finput[EnergyIndex]->Get("hmc_full_py_SidebandRegion1");
    htemp->Scale(luminosity[EnergyIndex]);
    if (hmc_sb!=0) hmc_sb->Add(htemp);
    else hmc_sb = (TH1F*)htemp->Clone(Form("%s_clone", htemp->GetName()));
    delete htemp;
    htemp = (TH1F*)finput[EnergyIndex]->Get("hmc_full_py_SidebandRegion3");
    htemp->Scale(luminosity[EnergyIndex]);
    if (hmc_sb!=0) hmc_sb->Add(htemp);
    else hmc_sb = (TH1F*)htemp->Clone(Form("%s_clone", htemp->GetName()));
    delete htemp;
    cout << "... SB MC... ";
    htemp = (TH1F*)finput[EnergyIndex]->Get("hmc_full_py_SignalRegion");
    htemp->Scale(luminosity[EnergyIndex]);
    if (hmc!=0) hmc->Add(htemp);
    else hmc = (TH1F*)htemp->Clone(Form("%s_clone", htemp->GetName()));
    delete htemp;
    cout << "... Signal MC... ";

    TTree* tdata_zerosb = (TTree*)finput[EnergyIndex]->Get("SB13_Data");
    TTree* tbkg_zerosb = (TTree*)finput[EnergyIndex]->Get("SB13_MC");
    TTree* tbkg_signal = (TTree*)finput[EnergyIndex]->Get("SR_NonZX");
    TTree* tbkg_zxsignal = (TTree*)finput[EnergyIndex]->Get("SR_ZX");
    TTree* tbkg_zxsb = (TTree*)finput[EnergyIndex]->Get("SB2_ZX");
    tSystContainer[0][EnergyIndex]=tbkg_zerosb;
    tSystContainer[1][EnergyIndex]=tdata_zerosb;
    tSystContainer[2][EnergyIndex]=tbkg_signal;
    tSystContainer[3][EnergyIndex]=tbkg_zxsignal;
    tSystContainer[4][EnergyIndex]=tbkg_zxsb;
    cout << "done...";
  }
  hzx_sb->Scale(hzx->Integral()/hzx_sb->Integral());
  cout << endl;

  hdata->SetMarkerColor(kBlack);
  hdata->SetMarkerStyle(20);
  hdata->SetLineStyle(1);
  hdata->SetLineColor(kBlack);
  hdata->SetLineWidth(1);
  hdata->Sumw2();
  hmc->SetLineStyle(1);
  hmc->SetLineColor(kBlack);
  hmc->SetLineWidth(2);
  hmc_sb->SetLineStyle(1);
  hmc_sb->SetLineColor(kBlue);
  hmc_sb->SetLineWidth(2);

  TString OUTPUT_W_NAME = "LifetimeKD_SystFitsRecord_DataMC_";
  OUTPUT_W_NAME.Append("newSIP_");
  OUTPUT_W_NAME.Append(user_folder[folder]);
  OUTPUT_W_NAME.Append(Form("_%s_", cSystVariable[iSyst]));
  OUTPUT_W_NAME.Append("AllTeV");
  OUTPUT_W_NAME.Append(".root");
  TString coutput_w = coutput_common + OUTPUT_W_NAME;
  foutput->cd();
  RooWorkspace* w = new RooWorkspace("w", kTRUE);
  w->importClassCode(RooRealFlooredSumPdf::Class(), kTRUE);

  RooCategory cat("channel", "channel");
  for (int erg_tev=7; erg_tev<9; erg_tev++){
    TString comstring;
    comstring.Form("%iTeV", erg_tev);
    cat.defineType(comstring, erg_tev);
    cat.setLabel(comstring);
  };


  RooRealVar* KD = new RooRealVar(cSystVariable[iSyst], cSystVariable_label[iSyst], 0, -1.0e6, 1.0e6);
  RooRealVar* weight = new RooRealVar("weight", "weight", 0, 0, 1.0e10);

  const int nGaussians=6;
  RooRealVar* varA[nGaussians][2];
  RooFormulaVar* varF[nGaussians][2];
  RooRealVar* varB[nGaussians][2];
  RooFormulaVar* varS[nGaussians][2];
  RooRealVar* extra_smear = new RooRealVar("extra_smear", "extra_smear", 1, 0.0001, 10);
  RooRealVar* center = new RooRealVar("C", "C", 0);
  for (int gg=0; gg<nGaussians; gg++){
    for (int e=0; e<2; e++){
      if (gg==0){
        varA[gg][e] = new RooRealVar(Form("A%i_cat%i", gg, e), Form("A%i_cat%i", gg, e), 1);
        varF[gg][e] = new RooFormulaVar(Form("Norm%i_cat%i", gg, e), Form("Norm%i_cat%i", gg, e), "@0", RooArgList(*(varA[gg][e])));
        varB[gg][e] = new RooRealVar(Form("B%i_cat%i", gg, e), Form("B%i_cat%i", gg, e), 80, 10, 1000);
        varS[gg][e] = new RooFormulaVar(Form("Smear%i_cat%i", gg, e), Form("Smear%i_cat%i", gg, e), "2*@0*@1", RooArgList(*(varB[gg][e]), *(extra_smear)));
      }
      else{
        varA[gg][e] = new RooRealVar(Form("A%i_cat%i", gg, e), Form("A%i_cat%i", gg, e), 0.1, 0, 1);
        varF[gg][e] = new RooFormulaVar(Form("Norm%i_cat%i", gg, e), Form("Norm%i_cat%i", gg, e), "@0*@1", RooArgList(*(varF[gg-1][e]), *(varA[gg][e])));
        varB[gg][e] = new RooRealVar(Form("B%i_cat%i", gg, e), Form("B%i_cat%i", gg, e), 3, 1, 20);
        varS[gg][e] = new RooFormulaVar(Form("Smear%i_cat%i", gg, e), Form("Smear%i_cat%i", gg, e), "2*@0*@1", RooArgList(*(varS[gg-1][e]), *(varB[gg][e])));
      }
    }
  }

  RooDataSet* fitdata;
  for (int tt=0; tt<5; tt++){
    if (tt==0){
      extra_smear->setVal(1);
      extra_smear->setConstant(true);
    }
    else if (tt==1 || tt==2){
      extra_smear->setConstant(false);
      extra_smear->setVal(1);
      for (int gg=0; gg<nGaussians; gg++){
        for (int e=0; e<2; e++){
          varA[gg][e]->setConstant(true);
          varB[gg][e]->setConstant(true);
        }
      }
    }

    RooArgSet* fitArgs = new RooArgSet(*KD, *weight);
    RooDataSet* tmpdata[2]={
      new RooDataSet(Form("%s_fitdata", tSystContainer[tt][0]->GetName()), Form("%s_0_fitdata", tSystContainer[tt][0]->GetName()), tSystContainer[tt][0], *fitArgs, 0, "weight"),
      new RooDataSet(Form("%s_fitdata", tSystContainer[tt][1]->GetName()), Form("%s_1_fitdata", tSystContainer[tt][1]->GetName()), tSystContainer[tt][1], *fitArgs, 0, "weight")
    };
    fitArgs->add(cat);
    fitdata = new RooDataSet(Form("%s_combdata", tSystContainer[tt][0]->GetName()), Form("%s_combdata", tSystContainer[tt][0]->GetName()), *fitArgs, Index(cat), WeightVar("weight"), Import("7TeV", *tmpdata[0]), Import("8TeV", *tmpdata[1]));
    cout << "Constructed datasets" << endl;

    RooBreitWigner* breit[2];
    RooGaussian* gaussians[nGaussians-1][2];
    RooRealFlooredSumPdf* tmp_pdf[2];
    for (int e=0; e<2; e++){
      RooBreitWigner* bw = new RooBreitWigner(Form("%s_bw", tSystContainer[tt][e]->GetName()), "", *KD, *center, *varS[0][e]);
      breit[e]=bw;
      for (int gg=1; gg<nGaussians; gg++){
        RooGaussian* my_g = new RooGaussian(Form("%s_gauss%i", tSystContainer[tt][e]->GetName(), gg), "", *KD, *center, *varS[gg][e]);
        gaussians[gg-1][e]=my_g;
      }
      cout << "Set the Gaussians" << endl;

      RooArgList funcList(*bw);
      RooArgList coeffList(*varF[0][e]);
      for (int gg=1; gg<4; gg++){
        funcList.add(*(gaussians[gg-1][e]));
        coeffList.add(*(varF[gg][e]));
      }
      tmp_pdf[e] = new RooRealFlooredSumPdf(Form("%s_%i_pdf", tSystContainer[tt][e]->GetName(), e), "", funcList, coeffList);
    }

    RooSimultaneous* sumpdf = new RooSimultaneous(Form("%s_pdf", tSystContainer[tt][0]->GetName()), "", cat);
    sumpdf->addPdf(*tmp_pdf[0], "7TeV");
    sumpdf->addPdf(*tmp_pdf[1], "8TeV");
    if (tt!=1) sumpdf->fitTo(*fitdata, SumW2Error(kFALSE));
    else sumpdf->fitTo(*fitdata, SumW2Error(kTRUE));
    w->import(*sumpdf);

    cout << "\n\nSmear: "
      << extra_smear->getVal()
      << endl;

    TString canvasname_2D = Form("cCanvas_FitComparisons_%s_%s", tSystContainer[tt][0]->GetName(), user_folder[folder]);
    canvasname_2D.Append("_AllTeV");
    TString yTitle = Form("Events / %.1f", (visualRange_SystVariable[iSyst][1] - visualRange_SystVariable[iSyst][0]) / nbins_SystVariable[iSyst]);
    TString xTitle = cSystVariable_label[iSyst];

    if (tt==1){
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
      //      if (erg_tev==7) text = pt->AddText(0.837, 0.45, "#font[42]{ 5.1 fb^{-1} (7 TeV)}");
      //      if (erg_tev==8) text = pt->AddText(0.837, 0.45, "#font[42]{19.7 fb^{-1} (8 TeV)}");
      text = pt->AddText(0.537, 0.45, "#font[42]{19.7 fb^{-1} (8 TeV) + 5.1 fb^{-1} (7 TeV)}");
      text->SetTextSize(0.0315);

      TCanvas* c2D = new TCanvas(canvasname_2D, "", 8, 30, 800, 800);
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

      RooPlot* rplot = new RooPlot(Form("rplot_%s", tSystContainer[tt][0]->GetName()), "", *KD, visualRange_SystVariable[iSyst][0], visualRange_SystVariable[iSyst][1], nbins_SystVariable[iSyst]);
      fitdata->plotOn(rplot, MarkerColor(kBlack), MarkerStyle(20), MarkerSize(1.0), LineWidth(1), XErrorSize(0), DataError(RooAbsData::Poisson));
      sumpdf->plotOn(rplot, LineColor(kRed), LineWidth(2), LineStyle(7));
      TH1F* hpdf = (TH1F*)sumpdf->createHistogram(Form("%s_histo", sumpdf->GetName()), *KD, Binning(nbins_SystVariable[iSyst], visualRange_SystVariable[iSyst][0], visualRange_SystVariable[iSyst][1]));
      hpdf->Scale(hdata->Integral()/hpdf->Integral());

      TH1F* hmc_clone = (TH1F*)hmc->Clone(Form("%s_clone", hmc->GetName()));
      hmc_clone->Scale(hdata->Integral()/hmc_clone->Integral());
      TH1F* hmc_sb_clone = (TH1F*)hmc_sb->Clone(Form("%s_clone", hmc_sb->GetName()));
      hmc_sb_clone->Scale(hdata->Integral()/hmc_sb_clone->Integral());

      TLegend *l2D = new TLegend(0.20, 0.65, 0.58, 0.88);
      l2D->SetBorderSize(0);
      l2D->SetTextFont(42);
      l2D->SetTextSize(0.03);
      l2D->SetLineColor(1);
      l2D->SetLineStyle(1);
      l2D->SetLineWidth(1);
      l2D->SetFillColor(0);
      l2D->SetFillStyle(0);

      rplot->GetXaxis()->SetTitle(xTitle);
      rplot->GetYaxis()->SetTitle(yTitle);
      rplot->GetXaxis()->SetNdivisions(505);
      rplot->GetXaxis()->SetLabelFont(42);
      rplot->GetXaxis()->SetLabelOffset(0.007);
      rplot->GetXaxis()->SetLabelSize(0.04);
      rplot->GetXaxis()->SetTitleSize(0.06);
      rplot->GetXaxis()->SetTitleOffset(0.9);
      rplot->GetXaxis()->SetTitleFont(42);
      rplot->GetYaxis()->SetNdivisions(505);
      rplot->GetYaxis()->SetLabelFont(42);
      rplot->GetYaxis()->SetLabelOffset(0.007);
      rplot->GetYaxis()->SetLabelSize(0.04);
      rplot->GetYaxis()->SetTitleSize(0.06);
      rplot->GetYaxis()->SetTitleOffset(1.1);
      rplot->GetYaxis()->SetTitleFont(42);
      hpdf->SetLineStyle(1);
      hpdf->SetLineColor(kRed);
      hpdf->SetLineWidth(2);
      TH1F* hpdf_clone = (TH1F*)hpdf->Clone("hpdf_temp");
      hpdf_clone->SetLineStyle(7);

      double max_plot = max(hmc->GetMaximum(), max(hpdf->GetMaximum(), max(hmc_sb->GetMaximum(), hdata->GetMaximum())));
      rplot->GetYaxis()->SetRangeUser(0, max_plot*1.5);
      rplot->GetXaxis()->SetRangeUser(-750, 750);


      l2D->AddEntry(hpdf_clone, "Fit to sideband data", "l");
      l2D->AddEntry(hpdf, "Fit projection", "l");
      l2D->AddEntry(hdata, "Observed sideband", "ep");
      //        l2D->AddEntry(hmc_sb_clone, "Bkg. (sideband regions)", "l");
      l2D->AddEntry(hmc_clone, "Bkg. (signal region)", "l");

      rplot->Draw();
      hpdf->Draw("histsame");
      //        hmc_sb_clone->Draw("histsame");
      hmc_clone->Draw("histsame");
      l2D->Draw();
      pt->Draw();

      TString strSidebandRegion = "Sidebands 1 and 3";
      TPaveText *pt10 = new TPaveText(0.68, 0.84, 0.90, 0.92, "brNDC");
      pt10->SetBorderSize(0);
      pt10->SetTextAlign(12);
      pt10->SetTextSize(0.03);
      pt10->SetFillStyle(0);
      pt10->SetTextFont(42);
      TText* text10 = pt10->AddText(0.01, 0.01, strSidebandRegion);
      pt10->Draw();

      TString strSIPcut = "New cut";
      strSIPcut.Prepend(", ");
      strSIPcut.Prepend(folder_id.c_str());
      TPaveText* pt20;
      pt20 = new TPaveText(0.68, 0.80, 0.90, 0.84, "brNDC");
      pt20->SetBorderSize(0);
      pt20->SetTextAlign(12);
      pt20->SetTextSize(0.03);
      pt20->SetFillStyle(0);
      pt20->SetTextFont(42);
      TText* text20 = pt20->AddText(0.01, 0.01, strSIPcut);
      pt20->Draw();

      foutput->WriteTObject(c2D);
      delete pt20;
      delete pt10;
      delete l2D;
      delete hmc_clone;
      delete hmc_sb_clone;
      delete hpdf_clone;
      delete hpdf;
      delete rplot;
      c2D->Close();
      delete pt;
    }
    else if (tt==0){
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
      //      if (erg_tev==7) text = pt->AddText(0.837, 0.45, "#font[42]{ 5.1 fb^{-1} (7 TeV)}");
      //      if (erg_tev==8) text = pt->AddText(0.837, 0.45, "#font[42]{19.7 fb^{-1} (8 TeV)}");
      text = pt->AddText(0.537, 0.45, "#font[42]{19.7 fb^{-1} (8 TeV) + 5.1 fb^{-1} (7 TeV)}");
      text->SetTextSize(0.0315);

      TCanvas* c2D = new TCanvas(canvasname_2D, "", 8, 30, 800, 800);
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

      RooPlot* rplot = new RooPlot(Form("rplot_%s", tSystContainer[tt][0]->GetName()), "", *KD, visualRange_SystVariable[iSyst][0], visualRange_SystVariable[iSyst][1], nbins_SystVariable[iSyst]);
      fitdata->plotOn(rplot, MarkerColor(kBlack), MarkerStyle(20), MarkerSize(1.0), LineWidth(1), XErrorSize(0), DataError(RooAbsData::SumW2));
      sumpdf->plotOn(rplot, LineColor(kRed), LineWidth(2), LineStyle(7));
      TH1F* hpdf = (TH1F*)sumpdf->createHistogram(Form("%s_histo", sumpdf->GetName()), *KD, Binning(nbins_SystVariable[iSyst], visualRange_SystVariable[iSyst][0], visualRange_SystVariable[iSyst][1]));

      TH1F* hmc_sb_clone = (TH1F*)hmc_sb->Clone(Form("%s_clone", hmc_sb->GetName()));
      hpdf->Scale(hmc_sb_clone->Integral()/hpdf->Integral());
      TH1F* hmc_clone = (TH1F*)hmc->Clone(Form("%s_clone", hmc->GetName()));
      hmc_clone->Scale(hmc_sb_clone->Integral()/hmc_clone->Integral());

      TLegend *l2D = new TLegend(0.20, 0.65, 0.58, 0.88);
      l2D->SetBorderSize(0);
      l2D->SetTextFont(42);
      l2D->SetTextSize(0.03);
      l2D->SetLineColor(1);
      l2D->SetLineStyle(1);
      l2D->SetLineWidth(1);
      l2D->SetFillColor(0);
      l2D->SetFillStyle(0);

      rplot->GetXaxis()->SetTitle(xTitle);
      rplot->GetYaxis()->SetTitle(yTitle);
      rplot->GetXaxis()->SetNdivisions(505);
      rplot->GetXaxis()->SetLabelFont(42);
      rplot->GetXaxis()->SetLabelOffset(0.007);
      rplot->GetXaxis()->SetLabelSize(0.04);
      rplot->GetXaxis()->SetTitleSize(0.06);
      rplot->GetXaxis()->SetTitleOffset(0.9);
      rplot->GetXaxis()->SetTitleFont(42);
      rplot->GetYaxis()->SetNdivisions(505);
      rplot->GetYaxis()->SetLabelFont(42);
      rplot->GetYaxis()->SetLabelOffset(0.007);
      rplot->GetYaxis()->SetLabelSize(0.04);
      rplot->GetYaxis()->SetTitleSize(0.06);
      rplot->GetYaxis()->SetTitleOffset(1.1);
      rplot->GetYaxis()->SetTitleFont(42);
      hpdf->SetLineStyle(1);
      hpdf->SetLineColor(kRed);
      hpdf->SetLineWidth(2);
      TH1F* hpdf_clone = (TH1F*)hpdf->Clone("hpdf_temp");
      hpdf_clone->SetLineStyle(7);

      hmc_sb->Scale(hpdf->Integral()/hmc_sb->Integral());

      double max_plot = max(hpdf->GetMaximum(), hmc_sb->GetMaximum());
      rplot->GetYaxis()->SetRangeUser(0, max_plot*1.5);

      l2D->AddEntry(hpdf_clone, "Fit to sideband data", "l");
      l2D->AddEntry(hpdf, "Fit to data (histogram)", "l");
      l2D->AddEntry(hmc_sb_clone, "Bkg. sideband", "l");

      rplot->Draw();
      hpdf->Draw("histsame");
      hmc_clone->Draw("histsame");
      hmc_sb_clone->Draw("histsame");
      l2D->Draw();
      pt->Draw();

      TString strSidebandRegion = "Sidebands 1 and 3";
      TPaveText *pt10 = new TPaveText(0.68, 0.84, 0.90, 0.92, "brNDC");
      pt10->SetBorderSize(0);
      pt10->SetTextAlign(12);
      pt10->SetTextSize(0.03);
      pt10->SetFillStyle(0);
      pt10->SetTextFont(42);
      TText* text10 = pt10->AddText(0.01, 0.01, strSidebandRegion);
      pt10->Draw();

      TString strSIPcut = "New cut";
      strSIPcut.Prepend(", ");
      strSIPcut.Prepend(folder_id.c_str());
      TPaveText* pt20;
      pt20 = new TPaveText(0.68, 0.80, 0.90, 0.84, "brNDC");
      pt20->SetBorderSize(0);
      pt20->SetTextAlign(12);
      pt20->SetTextSize(0.03);
      pt20->SetFillStyle(0);
      pt20->SetTextFont(42);
      TText* text20 = pt20->AddText(0.01, 0.01, strSIPcut);
      pt20->Draw();

      foutput->WriteTObject(c2D);
      delete pt20;
      delete pt10;
      delete l2D;
      delete hmc_clone;
      delete hmc_sb_clone;
      delete hpdf_clone;
      delete hpdf;
      delete rplot;
      c2D->Close();
      delete pt;
    }
    else if (tt==4){
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
      //      if (erg_tev==7) text = pt->AddText(0.837, 0.45, "#font[42]{ 5.1 fb^{-1} (7 TeV)}");
      //      if (erg_tev==8) text = pt->AddText(0.837, 0.45, "#font[42]{19.7 fb^{-1} (8 TeV)}");
      text = pt->AddText(0.537, 0.45, "#font[42]{19.7 fb^{-1} (8 TeV) + 5.1 fb^{-1} (7 TeV)}");
      text->SetTextSize(0.0315);

      TCanvas* c2D = new TCanvas(canvasname_2D, "", 8, 30, 800, 800);
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

      TLegend *l2D = new TLegend(0.20, 0.65, 0.58, 0.88);
      l2D->SetBorderSize(0);
      l2D->SetTextFont(42);
      l2D->SetTextSize(0.03);
      l2D->SetLineColor(1);
      l2D->SetLineStyle(1);
      l2D->SetLineWidth(1);
      l2D->SetFillColor(0);
      l2D->SetFillStyle(0);

      TString yTitle = Form("Events / %.1f", (visualRange_SystVariable[iSyst][1] - visualRange_SystVariable[iSyst][0]) / nbins_SystVariable[iSyst]);
      TString xTitle = cSystVariable_label[iSyst];
      hzx->SetTitle("");
      hzx->GetXaxis()->SetTitle(xTitle);
      hzx->GetYaxis()->SetTitle(yTitle);
      hzx->GetXaxis()->SetNdivisions(505);
      hzx->GetXaxis()->SetLabelFont(42);
      hzx->GetXaxis()->SetLabelOffset(0.007);
      hzx->GetXaxis()->SetLabelSize(0.04);
      hzx->GetXaxis()->SetTitleSize(0.06);
      hzx->GetXaxis()->SetTitleOffset(0.9);
      hzx->GetXaxis()->SetTitleFont(42);
      hzx->GetYaxis()->SetNdivisions(505);
      hzx->GetYaxis()->SetLabelFont(42);
      hzx->GetYaxis()->SetLabelOffset(0.007);
      hzx->GetYaxis()->SetLabelSize(0.04);
      hzx->GetYaxis()->SetTitleSize(0.06);
      hzx->GetYaxis()->SetTitleOffset(1.1);
      hzx->GetYaxis()->SetTitleFont(42);
      hzx_sb->SetLineStyle(1);
      hzx_sb->SetLineColor(kRed);
      hzx_sb->SetLineWidth(2);
      hzx->SetLineStyle(1);
      hzx->SetLineColor(kBlack);
      hzx->SetLineWidth(2);
      TH1F* hzx_sb_clone = (TH1F*)hzx_sb->Clone(Form("%s_clone", hzx_sb->GetName()));
      TH1F* hzx_clone = (TH1F*)hzx->Clone(Form("%s_clone", hzx->GetName()));
      hzx_clone->SetFillStyle(3002);

      double max_plot = max(hzx->GetMaximum(), hzx_sb->GetMaximum());
      hzx->GetYaxis()->SetRangeUser(0, max_plot*1.5);
      hzx->GetXaxis()->SetRangeUser(-750, 750);

      l2D->AddEntry(hzx, "Z+X in signal region", "l");
      l2D->AddEntry(hzx_sb, "Z+X in sideband region", "l");

      hzx->Draw("hist");
      hzx_clone->Draw("e2same");
      hzx_sb->Draw("histsame");
      l2D->Draw();
      pt->Draw();

      TString strSidebandRegion = "Z+X Sideband 2";
      TPaveText *pt10 = new TPaveText(0.68, 0.84, 0.90, 0.92, "brNDC");
      pt10->SetBorderSize(0);
      pt10->SetTextAlign(12);
      pt10->SetTextSize(0.03);
      pt10->SetFillStyle(0);
      pt10->SetTextFont(42);
      TText* text10 = pt10->AddText(0.01, 0.01, strSidebandRegion);
      pt10->Draw();

      TString strSIPcut = "New cut";
      strSIPcut.Prepend(", ");
      strSIPcut.Prepend(folder_id.c_str());
      TPaveText* pt20;
      pt20 = new TPaveText(0.68, 0.80, 0.90, 0.84, "brNDC");
      pt20->SetBorderSize(0);
      pt20->SetTextAlign(12);
      pt20->SetTextSize(0.03);
      pt20->SetFillStyle(0);
      pt20->SetTextFont(42);
      TText* text20 = pt20->AddText(0.01, 0.01, strSIPcut);
      pt20->Draw();

      foutput->WriteTObject(c2D);
      delete pt20;
      delete pt10;
      delete l2D;
      delete hzx_clone;
      delete hzx_sb_clone;
      c2D->Close();
      delete pt;
    }

    delete sumpdf;
    for (int e=0; e<2; e++){
      delete breit[e];
      for (int gg=0; gg<nGaussians-1; gg++) gaussians[gg][e];
    }
    delete fitdata;
  }
  w->writeToFile(coutput_w);
  foutput->Close();
}

*/

/*
void fitLifetimeSystTrees(int folder, int erg_tev, int iSyst=0){
  if (iSyst>=nSystVars) return;

  RooFit::Verbose(false);
  RooFit::PrintLevel(-1000);
  RooMsgService::instance().setStreamStatus(1, false);
  RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);

  TString erg_dir;
  erg_dir.Form("LHC_%iTeV", erg_tev);
  TString comstring;
  comstring.Form("%iTeV", erg_tev);
  int EnergyIndex=0; if (erg_tev==8) EnergyIndex=1;

  string folder_id;
  if (folder==0) folder_id = "4#mu";
  if (folder==1) folder_id = "4e";
  if (folder==2) folder_id = "2e2#mu";

  TString INPUT_NAME = "LifetimeKD_SystTrees_DataMC_";
  INPUT_NAME.Append("newSIP_");
  INPUT_NAME.Append(user_folder[folder]);
  INPUT_NAME.Append(Form("_%s_", cSystVariable[iSyst]));
  INPUT_NAME.Append(comstring);
  INPUT_NAME.Append(".root");
  TString cinput_common = user_dir_hep + "Analysis/ShapeSystematics/";
  cinput_common.Append(Form("%s/", cSystVariable[iSyst]));
  TString cinput = cinput_common + INPUT_NAME;
  TFile* finput = new TFile(cinput, "read");

  TH1F* hzx = (TH1F*)finput->Get("hzx_full_py_SignalRegion");
  TH1F* hzx_sb = (TH1F*)finput->Get("hzx_full_py_SidebandRegion2");
  hzx_sb->Scale(hzx->Integral()/hzx_sb->Integral());

  TH1F* hmc = (TH1F*)finput->Get("hmc_full_py_SignalRegion");
  TH1F* hdata = (TH1F*)finput->Get("hdata_full_py_SidebandRegion1");
  hdata->Add((TH1F*)finput->Get("hdata_full_py_SidebandRegion3"));
  TH1F* hmc_sb = (TH1F*)finput->Get("hmc_full_py_SidebandRegion1");
  hmc_sb->Add((TH1F*)finput->Get("hmc_full_py_SidebandRegion3"));
  hdata->SetMarkerColor(kBlack);
  hdata->SetMarkerStyle(20);
  hdata->SetLineStyle(1);
  hdata->SetLineColor(kBlack);
  hdata->SetLineWidth(1);
  hdata->Sumw2();
  hmc->SetLineStyle(1);
  hmc->SetLineColor(kBlack);
  hmc->SetLineWidth(2);
  hmc_sb->SetLineStyle(1);
  hmc_sb->SetLineColor(kBlue);
  hmc_sb->SetLineWidth(2);

  hmc->Scale(luminosity[EnergyIndex]);
  hmc_sb->Scale(luminosity[EnergyIndex]);

  TTree* tdata_zerosb = (TTree*)finput->Get("SB13_Data");
  TTree* tbkg_zerosb = (TTree*)finput->Get("SB13_MC");
  TTree* tbkg_signal = (TTree*)finput->Get("SR_NonZX");
  TTree* tbkg_zxsignal = (TTree*)finput->Get("SR_ZX");
  TTree* tbkg_zxsb = (TTree*)finput->Get("SB2_ZX");
  TTree* tSystContainer[5]={
    tbkg_zerosb, tdata_zerosb, tbkg_signal, tbkg_zxsignal, tbkg_zxsb
  };

  TString OUTPUT_NAME = "LifetimeKD_SystFits_DataMC_";
  OUTPUT_NAME.Append("newSIP_");
  OUTPUT_NAME.Append(user_folder[folder]);
  OUTPUT_NAME.Append(Form("_%s_", cSystVariable[iSyst]));
  OUTPUT_NAME.Append(comstring);
  OUTPUT_NAME.Append(".root");
  TString coutput_common = user_dir_hep + "Analysis/ShapeSystematics/";
  coutput_common.Append(Form("%s/", cSystVariable[iSyst]));
  TString coutput = coutput_common + OUTPUT_NAME;

  TString OUTPUT_W_NAME = "LifetimeKD_SystFitsRecord_DataMC_";
  OUTPUT_W_NAME.Append("newSIP_");
  OUTPUT_W_NAME.Append(user_folder[folder]);
  OUTPUT_W_NAME.Append(Form("_%s_", cSystVariable[iSyst]));
  OUTPUT_W_NAME.Append(comstring);
  OUTPUT_W_NAME.Append(".root");
  TString coutput_w = coutput_common + OUTPUT_W_NAME;

  TFile* foutput = new TFile(coutput, "recreate");

  foutput->cd();
  RooWorkspace* w = new RooWorkspace("w", kTRUE);
  w->importClassCode(RooRealFlooredSumPdf::Class(), kTRUE);
  RooRealVar* KD = new RooRealVar(cSystVariable[iSyst], cSystVariable_label[iSyst], 0, -1.0e6, 1.0e6);
  RooRealVar* weight = new RooRealVar("weight", "weight", 0, 0, 1.0e10);

  const int nGaussians=6;
  RooRealVar* A1 = new RooRealVar("A1", "A1", 1);
  A1->setConstant(true);
  RooRealVar* A2 = new RooRealVar("A2", "A2", 0.1, 1e-10, 1);
  RooRealVar* A3 = new RooRealVar("A3", "A3", 0.1, 1e-10, 1);
  RooRealVar* A4 = new RooRealVar("A4", "A4", 0.1, 1e-10, 1);
  RooRealVar* A5 = new RooRealVar("A5", "A5", 0.1, 1e-10, 1);
  RooRealVar* A6 = new RooRealVar("A6", "A6", 0.1, 1e-10, 1);
  RooRealVar* varA[nGaussians]={ A1, A2, A3, A4, A5, A6 };

  RooFormulaVar* F3 = new RooFormulaVar("F3", "F3", "@0*@1", RooArgList(*A2, *A3));
  RooFormulaVar* F4 = new RooFormulaVar("F4", "F4", "@0*@1", RooArgList(*F3, *A4));
  RooFormulaVar* F5 = new RooFormulaVar("F5", "F5", "@0*@1", RooArgList(*F4, *A5));
  RooFormulaVar* F6 = new RooFormulaVar("F6", "F6", "@0*@1", RooArgList(*F5, *A6));

  RooRealVar* extra_smear = new RooRealVar("extra_smear", "extra_smear", 1, 0.0001, 10);
  RooRealVar* B1 = new RooRealVar("B1", "B1", 80, 10, 1000);
  RooRealVar* B2 = new RooRealVar("B2", "B2", 4, 1, 15);
  RooRealVar* B3 = new RooRealVar("B3", "B3", 2, 1, 15);
  RooRealVar* B4 = new RooRealVar("B4", "B4", 3, 1, 20);
  RooRealVar* B5 = new RooRealVar("B5", "B5", 3, 1, 20);
  RooRealVar* B6 = new RooRealVar("B6", "B6", 3, 1, 20);
  RooRealVar* varB[nGaussians]={ B1, B2, B3, B4, B5, B6 };

  RooFormulaVar* S0 = new RooFormulaVar("S0", "S0", "2*@0", RooArgList(*B1));
  RooFormulaVar* S1 = new RooFormulaVar("S1", "S1", "@0*@1", RooArgList(*S0, *extra_smear));
  RooFormulaVar* S2 = new RooFormulaVar("S2", "S2", "@0*@1", RooArgList(*S1, *B2));
  RooFormulaVar* S3 = new RooFormulaVar("S3", "S3", "@0*@1", RooArgList(*S2, *B3));
  RooFormulaVar* S4 = new RooFormulaVar("S4", "S4", "@0*@1", RooArgList(*S3, *B4));
  RooFormulaVar* S5 = new RooFormulaVar("S5", "S5", "@0*@1", RooArgList(*S4, *B5));
  RooFormulaVar* S6 = new RooFormulaVar("S6", "S6", "@0*@1", RooArgList(*S5, *B6));

  RooRealVar* C1 = new RooRealVar("C1", "C1", 0);
  RooRealVar* C2 = new RooRealVar("C2", "C2", 0);
  RooRealVar* C3 = new RooRealVar("C3", "C3", 0);
  RooRealVar* C4 = new RooRealVar("C4", "C4", 0);
  RooRealVar* C5 = new RooRealVar("C5", "C5", 0);
  RooRealVar* C6 = new RooRealVar("C6", "C6", 0);

  RooDataSet* fitdata;
  for (int tt=0; tt<5; tt++){
    if (tt==0){
      extra_smear->setVal(1);
      extra_smear->setConstant(true);
    }
    else if (tt==1 || tt==2){
      extra_smear->setConstant(false);
      extra_smear->setVal(1);
      for (int vv=0; vv<nGaussians; vv++){
        varA[vv]->setConstant(true);
        varB[vv]->setConstant(true);
      }
    }

    fitdata = new RooDataSet(Form("%s_fitdata", tSystContainer[tt]->GetName()), Form("%s_fitdata", tSystContainer[tt]->GetName()), tSystContainer[tt], RooArgSet(*KD, *weight), 0, "weight");
    RooBreitWigner* bw = new RooBreitWigner(Form("%s_bw", tSystContainer[tt]->GetName()), "", *KD, *C1, *S1);
    RooGaussian* g1 = new RooGaussian(Form("%s_gauss1", tSystContainer[tt]->GetName()), "", *KD, *C1, *S1);
    RooGaussian* g2 = new RooGaussian(Form("%s_gauss2", tSystContainer[tt]->GetName()), "", *KD, *C2, *S2);
    RooGaussian* g3 = new RooGaussian(Form("%s_gauss3", tSystContainer[tt]->GetName()), "", *KD, *C3, *S3);
    RooGaussian* g4 = new RooGaussian(Form("%s_gauss4", tSystContainer[tt]->GetName()), "", *KD, *C4, *S4);
    RooGaussian* g5 = new RooGaussian(Form("%s_gauss5", tSystContainer[tt]->GetName()), "", *KD, *C5, *S5);
    RooGaussian* g6 = new RooGaussian(Form("%s_gauss6", tSystContainer[tt]->GetName()), "", *KD, *C6, *S6);

    cout << "Set the Gaussians" << endl;


    RooArgList* funcList;
    RooArgList* coeffList;
    if (folder==0){
      funcList = new RooArgList(*bw, *g2, *g3, *g4);
      coeffList = new RooArgList(*A1, *A2, *F3, *F4);
    }
    else{
      funcList = new RooArgList(*bw, *g2, *g3, *g4);
      coeffList = new RooArgList(*A1, *A2, *F3, *F4);
    }
    RooRealFlooredSumPdf* sumpdf = new RooRealFlooredSumPdf(Form("%s_sumpdf", tSystContainer[tt]->GetName()), "", *funcList, *coeffList);
    sumpdf->fitTo(*fitdata, SumW2Error(kFALSE));
    w->import(*sumpdf);

    cout << "\n\n"
      << extra_smear->getVal() << '\t'
      << A1->getVal() << '\t'
      << A2->getVal() << '\t'
      << A3->getVal() << '\t'
      << A4->getVal() << '\t'
      << A5->getVal() << '\t'
      << A6->getVal() << '\t'
      << B1->getVal() << '\t'
      << B2->getVal() << '\t'
      << B3->getVal() << '\t'
      << B4->getVal() << '\t'
      << B5->getVal() << '\t'
      << B6->getVal()
      << endl;

    TString canvasname_2D = Form("cCanvas_FitComparisons_%s_%s", tSystContainer[tt]->GetName(), user_folder[folder]);
    canvasname_2D.Append("_");
    canvasname_2D.Append(comstring);
    TString yTitle = Form("Events / %.1f", (visualRange_SystVariable[iSyst][1] - visualRange_SystVariable[iSyst][0]) / nbins_SystVariable[iSyst]);
    TString xTitle = cSystVariable_label[iSyst];

    if (tt==1){
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
      if (erg_tev==7) text = pt->AddText(0.837, 0.45, "#font[42]{ 5.1 fb^{-1} (7 TeV)}");
      if (erg_tev==8) text = pt->AddText(0.837, 0.45, "#font[42]{19.7 fb^{-1} (8 TeV)}");
      text->SetTextSize(0.0315);

      TCanvas* c2D = new TCanvas(canvasname_2D, "", 8, 30, 800, 800);
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

      RooPlot* rplot = new RooPlot(Form("rplot_%s", tSystContainer[tt]->GetName()), "", *KD, visualRange_SystVariable[iSyst][0], visualRange_SystVariable[iSyst][1], nbins_SystVariable[iSyst]);
      fitdata->plotOn(rplot, MarkerColor(kBlack), MarkerStyle(20), MarkerSize(1.0), LineWidth(1), XErrorSize(0), DataError(RooAbsData::Poisson));
      sumpdf->plotOn(rplot, LineColor(kRed), LineWidth(2), LineStyle(7));
      TH1F* hpdf = (TH1F*)sumpdf->createHistogram(Form("%s_histo", sumpdf->GetName()), *KD, Binning(nbins_SystVariable[iSyst], visualRange_SystVariable[iSyst][0], visualRange_SystVariable[iSyst][1]));
      hpdf->Scale(hdata->Integral()/hpdf->Integral());

      TH1F* hmc_clone = (TH1F*)hmc->Clone(Form("%s_clone", hmc->GetName()));
      hmc_clone->Scale(hdata->Integral()/hmc_clone->Integral());
      TH1F* hmc_sb_clone = (TH1F*)hmc_sb->Clone(Form("%s_clone", hmc_sb->GetName()));
      hmc_sb_clone->Scale(hdata->Integral()/hmc_sb_clone->Integral());

      TLegend *l2D = new TLegend(0.20, 0.65, 0.58, 0.88);
      l2D->SetBorderSize(0);
      l2D->SetTextFont(42);
      l2D->SetTextSize(0.03);
      l2D->SetLineColor(1);
      l2D->SetLineStyle(1);
      l2D->SetLineWidth(1);
      l2D->SetFillColor(0);
      l2D->SetFillStyle(0);

      rplot->GetXaxis()->SetTitle(xTitle);
      rplot->GetYaxis()->SetTitle(yTitle);
      rplot->GetXaxis()->SetNdivisions(505);
      rplot->GetXaxis()->SetLabelFont(42);
      rplot->GetXaxis()->SetLabelOffset(0.007);
      rplot->GetXaxis()->SetLabelSize(0.04);
      rplot->GetXaxis()->SetTitleSize(0.06);
      rplot->GetXaxis()->SetTitleOffset(0.9);
      rplot->GetXaxis()->SetTitleFont(42);
      rplot->GetYaxis()->SetNdivisions(505);
      rplot->GetYaxis()->SetLabelFont(42);
      rplot->GetYaxis()->SetLabelOffset(0.007);
      rplot->GetYaxis()->SetLabelSize(0.04);
      rplot->GetYaxis()->SetTitleSize(0.06);
      rplot->GetYaxis()->SetTitleOffset(1.1);
      rplot->GetYaxis()->SetTitleFont(42);
      hpdf->SetLineStyle(1);
      hpdf->SetLineColor(kRed);
      hpdf->SetLineWidth(2);
      TH1F* hpdf_clone = (TH1F*)hpdf->Clone("hpdf_temp");
      hpdf_clone->SetLineStyle(7);

      double max_plot = max(hmc->GetMaximum(), max(hpdf->GetMaximum(), max(hmc_sb->GetMaximum(), hdata->GetMaximum())));
      rplot->GetYaxis()->SetRangeUser(0, max_plot*1.5);
      rplot->GetXaxis()->SetRangeUser(-750, 750);


      l2D->AddEntry(hpdf_clone, "Fit to sideband data", "l");
      l2D->AddEntry(hpdf, "Fit projection", "l");
      l2D->AddEntry(hdata, "Observed sideband", "ep");
      //        l2D->AddEntry(hmc_sb_clone, "Bkg. (sideband regions)", "l");
      l2D->AddEntry(hmc_clone, "Bkg. (signal region)", "l");

      rplot->Draw();
      hpdf->Draw("histsame");
      //        hmc_sb_clone->Draw("histsame");
      hmc_clone->Draw("histsame");
      l2D->Draw();
      pt->Draw();

      TString strSidebandRegion = "Sidebands 1 and 3";
      TPaveText *pt10 = new TPaveText(0.68, 0.84, 0.90, 0.92, "brNDC");
      pt10->SetBorderSize(0);
      pt10->SetTextAlign(12);
      pt10->SetTextSize(0.03);
      pt10->SetFillStyle(0);
      pt10->SetTextFont(42);
      TText* text10 = pt10->AddText(0.01, 0.01, strSidebandRegion);
      pt10->Draw();

      TString strSIPcut = "New cut";
      strSIPcut.Prepend(", ");
      strSIPcut.Prepend(folder_id.c_str());
      TPaveText* pt20;
      pt20 = new TPaveText(0.68, 0.80, 0.90, 0.84, "brNDC");
      pt20->SetBorderSize(0);
      pt20->SetTextAlign(12);
      pt20->SetTextSize(0.03);
      pt20->SetFillStyle(0);
      pt20->SetTextFont(42);
      TText* text20 = pt20->AddText(0.01, 0.01, strSIPcut);
      pt20->Draw();

      foutput->WriteTObject(c2D);
      delete pt20;
      delete pt10;
      delete l2D;
      delete hmc_clone;
      delete hmc_sb_clone;
      delete hpdf_clone;
      delete hpdf;
      delete rplot;
      c2D->Close();
      delete pt;
    }
    else if (tt==0){
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
      if (erg_tev==7) text = pt->AddText(0.837, 0.45, "#font[42]{ 5.1 fb^{-1} (7 TeV)}");
      if (erg_tev==8) text = pt->AddText(0.837, 0.45, "#font[42]{19.7 fb^{-1} (8 TeV)}");
      text->SetTextSize(0.0315);

      TCanvas* c2D = new TCanvas(canvasname_2D, "", 8, 30, 800, 800);
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

      RooPlot* rplot = new RooPlot(Form("rplot_%s", tSystContainer[tt]->GetName()), "", *KD, visualRange_SystVariable[iSyst][0], visualRange_SystVariable[iSyst][1], nbins_SystVariable[iSyst]);
      fitdata->plotOn(rplot, MarkerColor(kBlack), MarkerStyle(20), MarkerSize(1.0), LineWidth(1), XErrorSize(0), DataError(RooAbsData::SumW2));
      sumpdf->plotOn(rplot, LineColor(kRed), LineWidth(2), LineStyle(7));
      TH1F* hpdf = (TH1F*)sumpdf->createHistogram(Form("%s_histo", sumpdf->GetName()), *KD, Binning(nbins_SystVariable[iSyst], visualRange_SystVariable[iSyst][0], visualRange_SystVariable[iSyst][1]));

      TH1F* hmc_sb_clone = (TH1F*)hmc_sb->Clone(Form("%s_clone", hmc_sb->GetName()));
      hpdf->Scale(hmc_sb_clone->Integral()/hpdf->Integral());
      TH1F* hmc_clone = (TH1F*)hmc->Clone(Form("%s_clone", hmc->GetName()));
      hmc_clone->Scale(hmc_sb_clone->Integral()/hmc_clone->Integral());

      TLegend *l2D = new TLegend(0.20, 0.65, 0.58, 0.88);
      l2D->SetBorderSize(0);
      l2D->SetTextFont(42);
      l2D->SetTextSize(0.03);
      l2D->SetLineColor(1);
      l2D->SetLineStyle(1);
      l2D->SetLineWidth(1);
      l2D->SetFillColor(0);
      l2D->SetFillStyle(0);

      rplot->GetXaxis()->SetTitle(xTitle);
      rplot->GetYaxis()->SetTitle(yTitle);
      rplot->GetXaxis()->SetNdivisions(505);
      rplot->GetXaxis()->SetLabelFont(42);
      rplot->GetXaxis()->SetLabelOffset(0.007);
      rplot->GetXaxis()->SetLabelSize(0.04);
      rplot->GetXaxis()->SetTitleSize(0.06);
      rplot->GetXaxis()->SetTitleOffset(0.9);
      rplot->GetXaxis()->SetTitleFont(42);
      rplot->GetYaxis()->SetNdivisions(505);
      rplot->GetYaxis()->SetLabelFont(42);
      rplot->GetYaxis()->SetLabelOffset(0.007);
      rplot->GetYaxis()->SetLabelSize(0.04);
      rplot->GetYaxis()->SetTitleSize(0.06);
      rplot->GetYaxis()->SetTitleOffset(1.1);
      rplot->GetYaxis()->SetTitleFont(42);
      hpdf->SetLineStyle(1);
      hpdf->SetLineColor(kRed);
      hpdf->SetLineWidth(2);
      TH1F* hpdf_clone = (TH1F*)hpdf->Clone("hpdf_temp");
      hpdf_clone->SetLineStyle(7);

      hmc_sb->Scale(hpdf->Integral()/hmc_sb->Integral());

      double max_plot = max(hpdf->GetMaximum(), hmc_sb->GetMaximum());
      rplot->GetYaxis()->SetRangeUser(0, max_plot*1.5);

      l2D->AddEntry(hpdf_clone, "Fit to sideband data", "l");
      l2D->AddEntry(hpdf, "Fit to data (histogram)", "l");
      l2D->AddEntry(hmc_sb_clone, "Bkg. sideband", "l");

      rplot->Draw();
      hpdf->Draw("histsame");
      hmc_clone->Draw("histsame");
      hmc_sb_clone->Draw("histsame");
      l2D->Draw();
      pt->Draw();

      TString strSidebandRegion = "Sidebands 1 and 3";
      TPaveText *pt10 = new TPaveText(0.68, 0.84, 0.90, 0.92, "brNDC");
      pt10->SetBorderSize(0);
      pt10->SetTextAlign(12);
      pt10->SetTextSize(0.03);
      pt10->SetFillStyle(0);
      pt10->SetTextFont(42);
      TText* text10 = pt10->AddText(0.01, 0.01, strSidebandRegion);
      pt10->Draw();

      TString strSIPcut = "New cut";
      strSIPcut.Prepend(", ");
      strSIPcut.Prepend(folder_id.c_str());
      TPaveText* pt20;
      pt20 = new TPaveText(0.68, 0.80, 0.90, 0.84, "brNDC");
      pt20->SetBorderSize(0);
      pt20->SetTextAlign(12);
      pt20->SetTextSize(0.03);
      pt20->SetFillStyle(0);
      pt20->SetTextFont(42);
      TText* text20 = pt20->AddText(0.01, 0.01, strSIPcut);
      pt20->Draw();

      foutput->WriteTObject(c2D);
      delete pt20;
      delete pt10;
      delete l2D;
      delete hmc_clone;
      delete hmc_sb_clone;
      delete hpdf_clone;
      delete hpdf;
      delete rplot;
      c2D->Close();
      delete pt;
    }
    else if (tt==4){
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
      if (erg_tev==7) text = pt->AddText(0.837, 0.45, "#font[42]{ 5.1 fb^{-1} (7 TeV)}");
      if (erg_tev==8) text = pt->AddText(0.837, 0.45, "#font[42]{19.7 fb^{-1} (8 TeV)}");
      text->SetTextSize(0.0315);

      TCanvas* c2D = new TCanvas(canvasname_2D, "", 8, 30, 800, 800);
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

      TLegend *l2D = new TLegend(0.20, 0.65, 0.58, 0.88);
      l2D->SetBorderSize(0);
      l2D->SetTextFont(42);
      l2D->SetTextSize(0.03);
      l2D->SetLineColor(1);
      l2D->SetLineStyle(1);
      l2D->SetLineWidth(1);
      l2D->SetFillColor(0);
      l2D->SetFillStyle(0);

      TString yTitle = Form("Events / %.1f", (visualRange_SystVariable[iSyst][1] - visualRange_SystVariable[iSyst][0]) / nbins_SystVariable[iSyst]);
      TString xTitle = cSystVariable_label[iSyst];
      hzx->SetTitle("");
      hzx->GetXaxis()->SetTitle(xTitle);
      hzx->GetYaxis()->SetTitle(yTitle);
      hzx->GetXaxis()->SetNdivisions(505);
      hzx->GetXaxis()->SetLabelFont(42);
      hzx->GetXaxis()->SetLabelOffset(0.007);
      hzx->GetXaxis()->SetLabelSize(0.04);
      hzx->GetXaxis()->SetTitleSize(0.06);
      hzx->GetXaxis()->SetTitleOffset(0.9);
      hzx->GetXaxis()->SetTitleFont(42);
      hzx->GetYaxis()->SetNdivisions(505);
      hzx->GetYaxis()->SetLabelFont(42);
      hzx->GetYaxis()->SetLabelOffset(0.007);
      hzx->GetYaxis()->SetLabelSize(0.04);
      hzx->GetYaxis()->SetTitleSize(0.06);
      hzx->GetYaxis()->SetTitleOffset(1.1);
      hzx->GetYaxis()->SetTitleFont(42);
      hzx_sb->SetLineStyle(1);
      hzx_sb->SetLineColor(kRed);
      hzx_sb->SetLineWidth(2);
      hzx->SetLineStyle(1);
      hzx->SetLineColor(kBlack);
      hzx->SetLineWidth(2);
      TH1F* hzx_sb_clone = (TH1F*)hzx_sb->Clone(Form("%s_clone", hzx_sb->GetName()));
      TH1F* hzx_clone = (TH1F*)hzx->Clone(Form("%s_clone", hzx->GetName()));
      hzx_clone->SetFillStyle(3002);

      double max_plot = max(hzx->GetMaximum(), hzx_sb->GetMaximum());
      hzx->GetYaxis()->SetRangeUser(0, max_plot*1.5);
      hzx->GetXaxis()->SetRangeUser(-750, 750);

      l2D->AddEntry(hzx, "Z+X in signal region", "l");
      l2D->AddEntry(hzx_sb, "Z+X in sideband region", "l");

      hzx->Draw("hist");
      hzx_clone->Draw("e2same");
      hzx_sb->Draw("histsame");
      l2D->Draw();
      pt->Draw();

      TString strSidebandRegion = "Z+X Sideband 2";
      TPaveText *pt10 = new TPaveText(0.68, 0.84, 0.90, 0.92, "brNDC");
      pt10->SetBorderSize(0);
      pt10->SetTextAlign(12);
      pt10->SetTextSize(0.03);
      pt10->SetFillStyle(0);
      pt10->SetTextFont(42);
      TText* text10 = pt10->AddText(0.01, 0.01, strSidebandRegion);
      pt10->Draw();

      TString strSIPcut = "New cut";
      strSIPcut.Prepend(", ");
      strSIPcut.Prepend(folder_id.c_str());
      TPaveText* pt20;
      pt20 = new TPaveText(0.68, 0.80, 0.90, 0.84, "brNDC");
      pt20->SetBorderSize(0);
      pt20->SetTextAlign(12);
      pt20->SetTextSize(0.03);
      pt20->SetFillStyle(0);
      pt20->SetTextFont(42);
      TText* text20 = pt20->AddText(0.01, 0.01, strSIPcut);
      pt20->Draw();

      foutput->WriteTObject(c2D);
      delete pt20;
      delete pt10;
      delete l2D;
      delete hzx_clone;
      delete hzx_sb_clone;
      c2D->Close();
      delete pt;
    }

    delete sumpdf;
    delete bw;
    delete g1;
    delete g2;
    delete g3;
    delete g4;
    delete g5;
    delete g6;
    delete fitdata;
  }
  w->writeToFile(coutput_w);
  foutput->Close();
}

*/