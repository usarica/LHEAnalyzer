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
#include "./data/ZZ4l_125_6_Samples.h"

using namespace std;
using namespace ROOT::Math;


const int nSystVars=32;
char* cSystVariable[nSystVars] = {
"Txy", "Dxy", "delTxy", "delDxy", "pToverMZZ",
"pullTxy","pullDxy",

"resoPV",
"resoIntV",
"pull4lVtx_x",
"pull4lVtx_y",
"pull4lVtx_z",

"diffRecoTruePV",
"diffRecoTrueIntV",

"reso4lVtx_x",
"reso4lVtx_y",
"reso4lVtx_z",

"resoOPVtx_x",
"resoOPVtx_y",
"resoOPVtx_z",

"tanhTxy",
"Dbkg",

"diffRecoOPV_BeamSpot_x",
"diffRecoOPV_BeamSpot_y",
"diffRecoOPV_BeamSpot_z",

"diffTrueOPV_BeamSpot_x",
"diffTrueOPV_BeamSpot_y",
"diffTrueOPV_BeamSpot_z",

"Txy_BeamSpot",
"diffTxy_BeamSpot_Txy",
"diffTxy_BeamSpot_TxyTrue",

"nOPVTracks"
};

char* cSystVariable_label[nSystVars] = {
  "T_{xy} (#mum)", "D_{xy} (#mum)", "#deltaT_{xy} (#mum)", "#deltaD_{xy} (#mum)", "p_{#lower[-0.2]{T}} / m_{4l}",
"(T_{xy} - c#Deltat_{true}) / #deltaT_{xy}",
"(D_{xy} - D^{true}_{xy}) / #deltaD_{xy}",

"(r^{reco}_{PV} - r^{true}_{PV}) / #sigma^{PV}_{T} . #hat{p}^{true}_{#lower[-0.2]{T}}",
"(r^{reco}_{4l} - r^{true}_{4l}) / #sigma^{4l}_{T} . #hat{p}^{true}_{#lower[-0.2]{T}}",
"(x^{reco}_{4l} - x^{true}_{4l}) / #sigma^{4l}_{x}",
"(y^{reco}_{4l} - y^{true}_{4l}) / #sigma^{4l}_{y}",
"(z^{reco}_{4l} - z^{true}_{4l}) / #sigma^{4l}_{z}",

"(r^{reco}_{#lower[-0.3]{PV}} - r^{#lower[-0.15]{true}}_{PV}) . #hat{p}^{#lower[0.2]{true}}_{#lower[-0.2]{T}} (#mum)",
"(r^{#lower[0.1]{reco}}_{#lower[-0.3]{4l}} - r^{#lower[-0.15]{true}}_{4l}) . #hat{p}^{#lower[0.2]{true}}_{#lower[-0.2]{T}} (#mum)",

"x^{reco}_{4l} - x^{true}_{4l} (#mum)",
"y^{reco}_{4l} - y^{true}_{4l} (#mum)",
"z^{reco}_{4l} - z^{true}_{4l} (#mum)",

"x^{reco}_{PV} - x^{true}_{PV} (#mum)",
"y^{reco}_{PV} - y^{true}_{PV} (#mum)",
"z^{reco}_{PV} - z^{true}_{PV} (#mum)",

"tanh(T_{xy} / 800 #mum)",
"D_{bkg}",

"x^{reco}_{PV} - x^{reco}_{Beam} (#mum)",
"y^{reco}_{PV} - y^{reco}_{Beam} (#mum)",
"z^{reco}_{PV} - z^{reco}_{Beam} (mm)",

"x^{true}_{PV} - x^{reco}_{Beam} (#mum)",
"y^{true}_{PV} - y^{reco}_{Beam} (#mum)",
"z^{true}_{PV} - z^{reco}_{Beam} (mm)",

"c#Deltat (#mum)",

"c#Deltat - T_{xy} (#mum)", "c#Deltat - c#Deltat_{true} (#mum)",

"N^{PV}_{tracks}"
};

/*
"T_{xy} (#mum)", // 0
"D_{xy} (#mum)", // 1
"#deltaT_{xy} (#mum)", // 2
"#deltaD_{xy} (#mum)", //3
"p_{T} / m_{4l}", //4
"(T_{xy} - c#Deltat_{true}) / #deltaT_{xy}", //5
"(D_{xy} - D^{true}_{xy}) / #deltaD_{xy}", //6
"(r^{reco}_{PV} - r^{true}_{PV}) / #sigma^{PV}_{T} . #hat{p}^{true}_{T}", //7
"(r^{reco}_{4l} - r^{true}_{4l}) / #sigma^{4l}_{T} . #hat{p}^{true}_{T}", //8
"(x^{reco}_{4l} - x^{true}_{4l}) / #sigma^{4l}_{x}", //9
"(y^{reco}_{4l} - y^{true}_{4l}) / #sigma^{4l}_{y}", //10
"(z^{reco}_{4l} - z^{true}_{4l}) / #sigma^{4l}_{z}", //11
"(r^{reco}_{PV} - r^{true}_{PV}) . #hat{p}^{true}_{T} (#mum)", //12
"(r^{reco}_{4l} - r^{true}_{4l}) . #hat{p}^{true}_{T} (#mum)", //13
"x^{reco}_{4l} - x^{true}_{4l} (#mum)", //14
"y^{reco}_{4l} - y^{true}_{4l} (#mum)", //15
"z^{reco}_{4l} - z^{true}_{4l} (#mum)", //16
"x^{reco}_{PV} - x^{true}_{PV} (#mum)", //17
"y^{reco}_{PV} - y^{true}_{PV} (#mum)", //18
"z^{reco}_{PV} - z^{true}_{PV} (#mum)", //19
"tanh(T_{xy} / 800 #mum)", //20
"D_{bkg}", //21
"x^{reco}_{PV} - x^{reco}_{Beam} (#mum)", //22
"y^{reco}_{PV} - y^{reco}_{Beam} (#mum)", //23
"z^{reco}_{PV} - z^{reco}_{Beam} (mm)", //24
"x^{true}_{PV} - x^{reco}_{Beam} (#mum)", //25
"y^{true}_{PV} - y^{reco}_{Beam} (#mum)", //26
"z^{true}_{PV} - z^{reco}_{Beam} (mm)", //27
"c#Deltat (#mum)", //28
"c#Deltat - T_{xy} (#mum)", //29
"c#Deltat - c#Deltat_{true} (#mum)", //30
"N^{tracks}_{PV}" //31
*/

int nbins_SystVariable[nSystVars] = {
60, 40, 150, 200, 100,
60, 60,

60,
60,
48,
48,
30,

62,
62,

60,
60,
60,

60,
60,
60,

50,
20,

60,
60,
30,

60,
60,
30,

120,

120,60,

50
};
double KD_limits_SystVariable[nSystVars][2] = {
	{ -1400, 1400 },
	{ -400, 400 },
	{ 0, 1280 },
	{ 0, 200 },
	{ 0, 4 },

	{ -15, 15 },
	{ -15, 15 },

  { -15, 15 },
	{ -15, 15 },
	{ -12, 12 },
	{ -12, 12 },
	{ -12, 12 },

	{ -310, 310 },
	{ -310, 310 },

  { -90, 90 },
  { -90, 90 },
  { -90, 90 },

  { -90, 90 },
  { -90, 90 },
  { -90, 90 },

  { -1, 1 },
  { 0, 1 },

  { -200, 200 },
  { -200, 200 },
  { -180, 180 }, // mm

  { -200, 200 },
  { -200, 200 },
  { -180, 180 }, // mm

  { -2800, 2800 },

  { -2400, 2400 },
  { -1200, 1200 },

  { 0, 200 }
};
double visualRange_SystVariable[nSystVars][2] = {
	{ -1400, 1400 },
	{ -84, 84 },
	{ 0, 640 },
	{ 2, 20 },
	{ 0, 1.65 },
	{ -15, 15 },
	{ -15, 15 },
	{ -15, 15 },
	{ -15, 15 },
	{ -12, 12 },
	{ -12, 12 },
	{ -12, 12 },

  { -305, 305 },
  { -305, 305 },

  { -90, 90 },
  { -90, 90 },
  { -90, 90 },
  { -90, 90 },
  { -90, 90 },
  { -90, 90 },

  { -1, 1 },
  { 0, 1 },

  { -200, 200 },
  { -200, 200 },
  { -180, 180 }, // mm

  { -200, 200 },
  { -200, 200 },
  { -180, 180 }, // mm

  { -2800, 2800 },

  { -2400, 2400 },
  { -1200, 1200 },

  { 0, 150 }
};

const int nMZZ_ranges=4;
float systZZMass_range[nMZZ_ranges][2] = {
	{ 70, 105.6 }, { 140.6, 170 }, { 170, 800 }, { 105.6, 140.6 }
};

TString cutLabel[6] ={
  "Old cut",
  "No cut",
  "Z1-SIP",
  "chi**2",
  "New cut",
  "New + old"
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

void symmetrize_PromptTemplates(TH1F* hrepair);


void produce_KDSystVariables(int folder, int erg_tev, int iSyst=0, int removeSIP=0){
	char TREE_NAME[]="SelectedTree";
	TString erg_dir;
	erg_dir.Form("LHC_%iTeV",erg_tev);
	TString comstring;
	comstring.Form("%iTeV",erg_tev);
	int EnergyIndex=0;if(erg_tev==8) EnergyIndex=1;

	TString OUTPUT_NAME = "LifetimeKD_SystVariables_DataQQBZZ_";
	if(removeSIP==1) OUTPUT_NAME.Append("noSIP_");
	else if(removeSIP==2) OUTPUT_NAME.Append("newSIP_");
	OUTPUT_NAME.Append(user_folder[folder]);
	OUTPUT_NAME.Append(Form("_%s_",cSystVariable[iSyst]));
	OUTPUT_NAME.Append(comstring);
	OUTPUT_NAME.Append(".root");

  float ggzz_offshell_sum=0,ggzz_80_100_sum=0;
	float zx_offshell_sum=0,zx_80_100_sum=0,zx_signal_sum=0;
  float qqzz_signal_sum=0, ggzz_signal_sum=0;
  float higgs_signal_sum[4]={ 0 };

// SIP cut scaling
	string processName[8] = {
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
  if (removeSIP==0) chratio.Append("Old cut");
  else if (removeSIP==1) chratio.Append("No cut");
  else if (removeSIP==2) chratio.Append("New cut");
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

	float MC_weight;
  float ZXfake_weight_SS;
  float ZXfake_weight_OS[6];

  float MC_weight_noxsec;
	float MC_weight_QQZZEWK=1;
	float MC_weight_Kfactor=1;
	float MC_weight_QQBGGProper[4];
	float ZZMass,ZZPt,ZZEta,ZZPhi;
	float KalmanCandVtx_x, KalmanCandVtx_y, KalmanCandVtx_z, KalmanCandVtx_chi2;
  float OfflinePrimaryVtx_x, OfflinePrimaryVtx_y, OfflinePrimaryVtx_z, OfflinePrimaryVtx_ndof;
	float GenPrimaryVtx_x,GenPrimaryVtx_y,GenPrimaryVtx_z;
	float GenIntVtx_x,GenIntVtx_y,GenIntVtx_z;
	float GenHMass,GenHPt,GenHPhi;
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
	float sigmaPV_xy,sigmaInt_xy;
  float CandVtx_x, CandVtx_y, CandVtx_z;
  float CandBSVtx_x, CandBSVtx_y, CandBSVtx_z;
  float TrueCandVtx_x, TrueCandVtx_y, TrueCandVtx_z;


  int CRflag=-1;
  float Lep1combRelIsoPF=0, Lep2combRelIsoPF=0, Lep3combRelIsoPF=0, Lep4combRelIsoPF=0;
  bool Lep3isID=true, Lep4isID=true;
  short Z1ids;

	TString cinput_common = user_dir_hep + erg_dir + "/";
	TString cinput_common_noSIP = user_dir_hep + "No_SIP/" + erg_dir + "/";
	TString coutput_common = user_dir_hep + "Analysis/Plots/";
	coutput_common.Append(Form("%s/",cSystVariable[iSyst]));
	TString mkdirCommand = "mkdir -p ";
	mkdirCommand.Append(coutput_common);
	gSystem->Exec(mkdirCommand);
	TString coutput = coutput_common + OUTPUT_NAME;
	TFile* foutput = new TFile(coutput,"recreate");

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
	}
	TChain* tggzz = new TChain(TREE_NAME);
	for (int smp = kGGOLDSamples; smp < kGGMCFMSamples; smp++){
		TString cinput_ggzz = cinput_ggzz_common_noSIP + sample_FullSim[smp] + "_Reprocessed.root";
		tggzz->Add(cinput_ggzz);
	}
	TChain* tzx = new TChain(TREE_NAME);
  TString cinput_zx = cinput_zx_common_noSIP;
	cinput_zx = cinput_zx + "CR/" + sample_FullSim[kAllSamples - 1] + ".root";
	tzx->Add(cinput_zx);

  TString cinput_data = cinput_common_noSIP + user_folder[3] + "/" + data_files[folder] + ".root";
	TChain* tdata = new TChain(TREE_NAME);
	tdata->Add(cinput_data);
	TChain* tc[ntrees] = { tdata, tqqzz, tggzz, tzx };

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
		if(smp >= kGGOLDSamples && smp < kGGMCFMSamples) continue;
		TString cinput_ggzz = cinput_ggzz_common_noSIP + sample_FullSim[smp] + "_Reprocessed.root";
		tggzz_extra->Add(cinput_ggzz);
	}
	TChain* thiggs[kGGHSamples];
	for (int smp = 0; smp < kGGHSamples; smp++){
		TString cinput_higgs = cinput_higgs_common_noSIP + sample_FullSim[smp] + ".root";
		cout << "Attaching Higgs under " << cinput_higgs << endl;
		thiggs[smp] = new TChain(TREE_NAME);
		thiggs[smp]->Add(cinput_higgs);
	}
	TChain* tc_extra[ntrees_extra] = { thiggs[0], tqqzz_extra, tggzz_extra, thiggs[1], thiggs[2], thiggs[3] };
	
	int nbins_KD = nbins_SystVariable[iSyst];
	double KD_limits[2] = { KD_limits_SystVariable[iSyst][0], KD_limits_SystVariable[iSyst][1] };

	TH2F* hdata_full = new TH2F("hdata_full", Form("%i TeV %s Data", erg_tev, user_folder[folder]), nMZZ_ranges, 0, nMZZ_ranges, nbins_KD, KD_limits[0], KD_limits[1]);
	hdata_full->GetYaxis()->SetTitle(cSystVariable_label[iSyst]);
	hdata_full->GetXaxis()->SetTitle("m_{4l} (GeV)");
	hdata_full->GetXaxis()->SetBinLabel(1,"70 - 105.6");
	hdata_full->GetXaxis()->SetBinLabel(2,"140.6 - 170");
	hdata_full->GetXaxis()->SetBinLabel(3,"170 - 800");
	hdata_full->GetXaxis()->SetBinLabel(4,"105.6 - 140.6");
	hdata_full->SetOption("colz");
	TH2F* hqqzz_full = (TH2F*) hdata_full->Clone();
	TH2F* hggzz_full = (TH2F*) hdata_full->Clone();
	TH2F* hzx_full = (TH2F*) hdata_full->Clone();
	TH2F* hmc_full = (TH2F*) hdata_full->Clone();
	hqqzz_full->SetNameTitle("hqqzz_full",Form("%i TeV %s q#bar{q}#rightarrow4l Background",erg_tev,user_folder[folder]));
	hggzz_full->SetNameTitle("hggzz_full",Form("%i TeV %s gg Background",erg_tev,user_folder[folder]));
	hzx_full->SetNameTitle("hzx_full",Form("%i TeV %s Z+X Background",erg_tev,user_folder[folder]));
	hmc_full->SetNameTitle("hmc_full",Form("%i TeV %s Total Background",erg_tev,user_folder[folder]));
	TH2F* hfull[ntrees+1] = { hdata_full, hqqzz_full, hggzz_full, hzx_full, hmc_full };
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
//        if (tc[tt]->GetBranchStatus("ZXfake_weightProper")) tc[tt]->SetBranchAddress("ZXfake_weightProper", &MC_weight);
//        else{
//          if (tc[tt]->GetBranchStatus("ZXfake_weight_SS")){
            tc[tt]->SetBranchAddress("ZXfake_weight_OS", ZXfake_weight_OS);
            tc[tt]->SetBranchAddress("ZXfake_weight_SS", &ZXfake_weight_SS);
//          }
//          else tc[tt]->SetBranchAddress("ZXfake_weight", &MC_weight);
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
//        }
      }
      if (tt == 2 || tt == 1){
        tc[tt]->SetBranchAddress("MC_weight_Kfactor", &MC_weight_Kfactor);
        tc[tt]->SetBranchAddress("MC_weight_QQBGGProper", MC_weight_QQBGGProper);
        if (tt == 1){
          tc[tt]->SetBranchAddress("MC_weight_QQZZEWK", &MC_weight_QQZZEWK);
        }
      }
    }
    tc[tt]->SetBranchAddress("ZZMass", &ZZMass);
    tc[tt]->SetBranchAddress("ZZPt", &ZZPt);
    tc[tt]->SetBranchAddress("ZZEta", &ZZEta);
    tc[tt]->SetBranchAddress("ZZPhi", &ZZPhi);
    tc[tt]->SetBranchAddress("OfflinePrimaryVtx_x", &OfflinePrimaryVtx_x);
    tc[tt]->SetBranchAddress("OfflinePrimaryVtx_y", &OfflinePrimaryVtx_y);
    tc[tt]->SetBranchAddress("OfflinePrimaryVtx_z", &OfflinePrimaryVtx_z);
    tc[tt]->SetBranchAddress("OfflinePrimaryVtx_ndof", &OfflinePrimaryVtx_ndof);
    tc[tt]->SetBranchAddress("KalmanCandVtx_x", &KalmanCandVtx_x);
    tc[tt]->SetBranchAddress("KalmanCandVtx_y", &KalmanCandVtx_y);
    tc[tt]->SetBranchAddress("KalmanCandVtx_z", &KalmanCandVtx_z);
    tc[tt]->SetBranchAddress("KalmanCandVtx_chi2", &KalmanCandVtx_chi2);
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
    tc_extra[tt]->SetBranchAddress("OfflinePrimaryVtx_ndof", &OfflinePrimaryVtx_ndof);
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
      if (removeSIP==2 && (fabs(Lep1_Z1SIP) >= 4 || fabs(Lep2_Z1SIP) >= 4 || fabs(Lep3_Z1SIP) >= 5 || fabs(Lep4_Z1SIP) >= 5 || KalmanCandVtx_chi2 >= 30)) continue;
      if (removeSIP==0 && (fabs(Lep1SIP) >= 4 || fabs(Lep2SIP) >= 4 || fabs(Lep3SIP) >= 4 || fabs(Lep4SIP) >= 4)) continue;
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
        if (removeSIP==0 && tc[tt]->GetBranchStatus("ZXfake_weight_OS")) MC_weight = ZXfake_weight_OS[0] * (targetCRcount[0][0][regionIndex] / (targetCRcount[0][0][regionIndex] + targetCRcount[1][0][regionIndex]));
        else if (removeSIP==1 && tc[tt]->GetBranchStatus("ZXfake_weight_OS")) MC_weight = ZXfake_weight_OS[1] * (targetCRcount[0][1][regionIndex] / (targetCRcount[0][1][regionIndex] + targetCRcount[1][1][regionIndex]));
        else if (removeSIP==2 && tc[tt]->GetBranchStatus("ZXfake_weight_OS")) MC_weight = ZXfake_weight_OS[4] * (targetCRcount[0][4][regionIndex] / (targetCRcount[0][4][regionIndex] + targetCRcount[1][4][regionIndex]));
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

//          if (removeSIP==2 && tt==3 && folder==0 && ZZMass>=110 && ZZMass<130) cout << ZZMass << '\t' << SSwgt_OSoverSS<< '\t' << ZXfake_weight_SS << '\t' << (targetCRcount[1][4][regionIndex] / (targetCRcount[0][4][regionIndex] + targetCRcount[1][4][regionIndex])) << endl;

          MC_weight = ZXfake_weight_SS*SSwgt_OSoverSS;
          if (removeSIP==0) MC_weight *= (targetCRcount[1][0][regionIndex] / (targetCRcount[0][0][regionIndex] + targetCRcount[1][0][regionIndex]));
          else if (removeSIP==1) MC_weight *= (targetCRcount[1][1][regionIndex] / (targetCRcount[0][1][regionIndex] + targetCRcount[1][1][regionIndex]));
          else if (removeSIP==2) MC_weight *= (targetCRcount[1][4][regionIndex] / (targetCRcount[0][4][regionIndex] + targetCRcount[1][4][regionIndex]));
        }
      }

      if (iSyst>=22 && iSyst<nSystVars-1){ // BE CAREFUL HERE, HARDCODED BOUNDS AHEAD
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

			if(iSyst==0 || iSyst==20) KD = Txy;
			if(iSyst==1) KD = Dxy;
			if(iSyst==2) KD = delTxy;
			if(iSyst==3) KD = delDxy;
			if(iSyst==12) KD = ( OfflinePrimaryVtx_x - GenPrimaryVtx_x )*cos(GenHPhi) + ( OfflinePrimaryVtx_y - GenPrimaryVtx_y )*sin(GenHPhi);
			if(iSyst==13) KD = ( KalmanCandVtx_x - GenIntVtx_x )*cos(GenHPhi) + ( KalmanCandVtx_y - GenIntVtx_y )*sin(GenHPhi);
      if (iSyst==14) KD = (KalmanCandVtx_x-GenIntVtx_x);
      if (iSyst==15) KD = (KalmanCandVtx_y-GenIntVtx_y);
      if (iSyst==16) KD = (KalmanCandVtx_z-GenIntVtx_z);
      if (iSyst==17) KD = (OfflinePrimaryVtx_x - GenPrimaryVtx_x);
      if (iSyst==18) KD = (OfflinePrimaryVtx_y - GenPrimaryVtx_y);
      if (iSyst==19) KD = (OfflinePrimaryVtx_z - GenPrimaryVtx_z);
      if (iSyst==22) KD = (OfflinePrimaryVtx_x - BeamPosX);
      if (iSyst==23) KD = (OfflinePrimaryVtx_y - BeamPosY);
      if (iSyst==24) KD = (OfflinePrimaryVtx_z - BeamPosZ)/1000.;
      if (iSyst==25) KD = (GenPrimaryVtx_x - BeamPosX);
      if (iSyst==26) KD = (GenPrimaryVtx_y - BeamPosY);
      if (iSyst==27) KD = (GenPrimaryVtx_z - BeamPosZ)/1000.;
      if (iSyst==28) KD = Txy_BS;
      if (iSyst==29) KD = Txy_BS - Txy;
      if (iSyst==30) KD = Txy_BS - Txy_true;

      KD = KD*10000.; // in 1 um

      if (iSyst==20){
        KD = tanh(KD/800.); if (KD==1.) KD=0.999999; if (KD==-1.) KD=-0.999999;
      }
      if (iSyst==21){
        KD = D_bkg; if (KD==1.) KD=0.999999;
      }

			if(iSyst==4) KD = ZZPt/ZZMass;
			if(iSyst==5) KD = (Txy-Txy_true)/delTxy;
			if(iSyst==6) KD = (Dxy-Dxy_true)/delDxy;
			if(iSyst==7) KD = ( ( OfflinePrimaryVtx_x - GenPrimaryVtx_x )*cos(GenHPhi) + ( OfflinePrimaryVtx_y - GenPrimaryVtx_y )*sin(GenHPhi) )/sigmaPV_xy;
			if(iSyst==8) KD = ( ( KalmanCandVtx_x - GenIntVtx_x )*cos(GenHPhi) + ( KalmanCandVtx_y - GenIntVtx_y )*sin(GenHPhi) )/sigmaInt_xy;
			if(iSyst==9) KD = (KalmanCandVtx_x-GenIntVtx_x)/sqrt(KalmanCandVtx_cov_xx);
			if(iSyst==10) KD = (KalmanCandVtx_y-GenIntVtx_y)/sqrt(KalmanCandVtx_cov_yy);
			if(iSyst==11) KD = (KalmanCandVtx_z-GenIntVtx_z)/sqrt(KalmanCandVtx_cov_zz);
      if(iSyst==31) KD = (OfflinePrimaryVtx_ndof+3.)/2.;

			float wgt = 1;

			int biny = hfull[tt]->GetYaxis()->FindBin(KD);
			for (int binx = 0; binx < nMZZ_ranges; binx++){
				if (tt == 3 && binx == 0 && ZZMass >= systZZMass_range[binx][0] + 30 && ZZMass < systZZMass_range[binx][1] + 30){ // ZX special fill
					wgt = MC_weight;
					hfull[tt]->AddBinContent(hfull[tt]->GetBin(binx + 1, biny), wgt);
          hfull[tt]->SetBinError(hfull[tt]->GetBin(binx + 1, biny), sqrt(wgt*wgt + pow(hfull[tt]->GetBinError(hfull[tt]->GetBin(binx + 1, biny)), 2)));
          if (ZZMass<130 && ZZMass>=110){
            zx_80_100_sum += wgt;
//            cout << ZZMass << '\t' << KalmanCandVtx_chi2<< '\t' << CRflag << '\t' << wgt << endl;
          }
//          cout << ZZMass << '\t' << KalmanCandVtx_chi2<< '\t' << CRflag << '\t' << wgt << endl;
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

			Lep1_Z1SIP=0;Lep2_Z1SIP=0;Lep3_Z1SIP=0;Lep4_Z1SIP=0;KalmanCandVtx_chi2=0;

			tc_extra[tt]->GetEntry(ev);
			if (removeSIP==0 && (fabs(Lep1SIP) >= 4 || fabs(Lep2SIP) >= 4 || fabs(Lep3SIP) >= 4 || fabs(Lep4SIP) >= 4)) continue;
      if (removeSIP==2 && (fabs(Lep1_Z1SIP) >= 4 || fabs(Lep2_Z1SIP) >= 4 || fabs(Lep3_Z1SIP) >= 5 || fabs(Lep4_Z1SIP) >= 5 || KalmanCandVtx_chi2 >= 30)) continue;

			int binx = nMZZ_ranges-1;
			if (ZZMass >= systZZMass_range[binx][1] || ZZMass < systZZMass_range[binx][0]) continue;

      if (iSyst>=22 && iSyst<nSystVars-1){ // BE CAREFUL HERE, HARDCODED BOUNDS AHEAD
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

      if (iSyst==0 || iSyst==20) KD = Txy;
      if (iSyst==1) KD = Dxy;
      if (iSyst==2) KD = delTxy;
      if (iSyst==3) KD = delDxy;
      if (iSyst==12) KD = (OfflinePrimaryVtx_x - GenPrimaryVtx_x)*cos(GenHPhi) + (OfflinePrimaryVtx_y - GenPrimaryVtx_y)*sin(GenHPhi);
      if (iSyst==13) KD = (KalmanCandVtx_x - GenIntVtx_x)*cos(GenHPhi) + (KalmanCandVtx_y - GenIntVtx_y)*sin(GenHPhi);
      if (iSyst==14) KD = (KalmanCandVtx_x-GenIntVtx_x);
      if (iSyst==15) KD = (KalmanCandVtx_y-GenIntVtx_y);
      if (iSyst==16) KD = (KalmanCandVtx_z-GenIntVtx_z);
      if (iSyst==17) KD = (OfflinePrimaryVtx_x - GenPrimaryVtx_x);
      if (iSyst==18) KD = (OfflinePrimaryVtx_y - GenPrimaryVtx_y);
      if (iSyst==19) KD = (OfflinePrimaryVtx_z - GenPrimaryVtx_z);

      if (iSyst==22) KD = (OfflinePrimaryVtx_x - BeamPosX);
      if (iSyst==23) KD = (OfflinePrimaryVtx_y - BeamPosY);
      if (iSyst==24) KD = (OfflinePrimaryVtx_z - BeamPosZ)/1000.;
      if (iSyst==25) KD = (GenPrimaryVtx_x - BeamPosX);
      if (iSyst==26) KD = (GenPrimaryVtx_y - BeamPosY);
      if (iSyst==27) KD = (GenPrimaryVtx_z - BeamPosZ)/1000.;
      if (iSyst==28) KD = Txy_BS;
      if (iSyst==29) KD = Txy_BS - Txy;
      if (iSyst==30) KD = Txy_BS - Txy_true;

      KD = KD*10000.; // in 1 um

      if (iSyst==20){
        KD = tanh(KD/800.); if (KD==1.) KD=0.999999; if (KD==-1.) KD=-0.999999;
      }
      if (iSyst==21){
        KD = D_bkg; if (KD==1.) KD=0.999999;;
      }

      if (iSyst==4) KD = ZZPt/ZZMass;
			if(iSyst==5) KD = (Txy-Txy_true)/delTxy;
			if(iSyst==6) KD = (Dxy-Dxy_true)/delDxy;
			if(iSyst==7) KD = ( ( OfflinePrimaryVtx_x - GenPrimaryVtx_x )*cos(GenHPhi) + ( OfflinePrimaryVtx_y - GenPrimaryVtx_y )*sin(GenHPhi) )/sigmaPV_xy;
			if(iSyst==8) KD = ( ( KalmanCandVtx_x - GenIntVtx_x )*cos(GenHPhi) + ( KalmanCandVtx_y - GenIntVtx_y )*sin(GenHPhi) )/sigmaInt_xy;
			if(iSyst==9) KD = (KalmanCandVtx_x-GenIntVtx_x)/sqrt(KalmanCandVtx_cov_xx);
			if(iSyst==10) KD = (KalmanCandVtx_y-GenIntVtx_y)/sqrt(KalmanCandVtx_cov_yy);
			if(iSyst==11) KD = (KalmanCandVtx_z-GenIntVtx_z)/sqrt(KalmanCandVtx_cov_zz);
      if(iSyst==31) KD = (OfflinePrimaryVtx_ndof+3.)/2.;

      float wgt = 1;

			int biny = hhiggs[0]->GetXaxis()->FindBin(KD);
			wgt = ( (tt==0 || tt>=3) ? MC_weight : MC_weight_noxsec);
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

  if (removeSIP==2) cout << zx_80_100_sum << '\t' << hzx_full->Integral(1, 1, 0, hzx_full->GetNbinsY()+1) << '\t' << yield_80_100_zx[EnergyIndex][folder] << '\t' << ratio_targetsignalyield[6][4][0]/ratio_targetsignalyield[6][0][0] << endl;

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
    if (removeSIP >= 1){
      hqqzz_full->SetBinContent(4, biny, (hqqzz_full->GetBinContent(4, biny))*(ratio_targetsignalyield[4][1][3]/ratio_targetsignalyield[4][0][3]));
      hqqzz_full->SetBinError(4, biny, (hqqzz_full->GetBinError(4, biny))*(ratio_targetsignalyield[4][1][3]/ratio_targetsignalyield[4][0][3]));
      if (removeSIP == 2){
        hqqzz_full->SetBinContent(4, biny, (hqqzz_full->GetBinContent(4, biny))*(ratio_targetsignalyield[4][4][3]/ratio_targetsignalyield[4][1][3]));
        hqqzz_full->SetBinError(4, biny, (hqqzz_full->GetBinError(4, biny))*(ratio_targetsignalyield[4][4][3]/ratio_targetsignalyield[4][1][3]));
      }
    }
    for (int hh=0; hh<kGGHSamples; hh++){
      hhiggs[hh]->SetBinContent(biny, (hhiggs[hh]->GetBinContent(biny))*(yield_signal_higgs[EnergyIndex][folder]*(ratio_targetsignalyield[0][1][3] / ratio_targetsignalyield[0][0][3]) / higgs_signal_sum[hh] / luminosity[EnergyIndex]));
      hhiggs[hh]->SetBinError(biny, (hhiggs[hh]->GetBinError(biny))*(yield_signal_higgs[EnergyIndex][folder]*(ratio_targetsignalyield[0][1][3] / ratio_targetsignalyield[0][0][3]) / higgs_signal_sum[hh] / luminosity[EnergyIndex]));
    }
  }
  for (int hh = 0; hh < kGGHSamples; hh++){
    if (removeSIP == 0) hhiggs[hh]->Scale(ratio_targetsignalyield[hh][0][3] / ratio_targetsignalyield[hh][1][3]);
    else if (removeSIP == 2) hhiggs[hh]->Scale(ratio_targetsignalyield[hh][4][3] / ratio_targetsignalyield[hh][1][3]);
  }
  if (removeSIP >= 1){
    for (int binx = 1; binx <= hggzz_full->GetNbinsX(); binx++){
      double scale = ratio_targetsignalyield[5][1][binx-1]/ratio_targetsignalyield[5][0][binx-1];
      if (removeSIP==2) scale = ratio_targetsignalyield[5][4][binx-1]/ratio_targetsignalyield[5][0][binx-1];

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
      double scale = ratio_targetsignalyield[6][1][binx-1]/ratio_targetsignalyield[6][0][binx-1];
      if (removeSIP==2) scale = ratio_targetsignalyield[6][4][binx-1]/ratio_targetsignalyield[6][0][binx-1];

      for (int biny = 0; biny <= hzx_full->GetNbinsY()+1; biny++){
        double bincontent = hzx_full->GetBinContent(binx, biny);
        double binerror = hzx_full->GetBinError(binx, biny);

        bincontent *= scale;
        binerror *= scale;
        hzx_full->SetBinContent(binx, biny, bincontent);
        hzx_full->SetBinError(binx, biny, binerror);
      }
    }
	}

  cout << "Signal gg norm: " << hggzz_full->Integral(4, 4, 0, hggzz_full->GetNbinsY()+1)*luminosity[EnergyIndex] << "\tqq norm: " << hqqzz_full->Integral(4, 4, 0, hqqzz_full->GetNbinsY()+1)*luminosity[EnergyIndex]
    << "\tZX norm: " << hzx_full->Integral(4, 4, 0, hzx_full->GetNbinsY()+1)*luminosity[EnergyIndex] << "\tHiggs norm: " << hhiggs[0]->Integral(0, hhiggs[0]->GetNbinsX()+1)*luminosity[EnergyIndex] << "\tHiggs (BSM) norm: " << hhiggs[3]->Integral(0, hhiggs[3]->GetNbinsX()+1)*luminosity[EnergyIndex] << endl;

  for(int sb=1;sb<=3;sb++) cout << "SB" << sb << " gg norm: " << hggzz_full->Integral(sb, sb, 0, hggzz_full->GetNbinsY()+1)*luminosity[EnergyIndex]
    << "\tqq norm: " << hqqzz_full->Integral(sb, sb, 0, hqqzz_full->GetNbinsY()+1)*luminosity[EnergyIndex]
    << "\tZX norm: " << hzx_full->Integral(sb, sb, 0, hzx_full->GetNbinsY()+1)*luminosity[EnergyIndex] << endl;

	foutput->cd();
	for (int tt = 1; tt < ntrees; tt++){
    if (tt==3 && (iSyst>=7 && iSyst<20)) continue;
    if (tt==3 && iSyst==27) continue;
    hfull[ntrees]->Add(hfull[tt]);
	}
	for (int tt = 0; tt < ntrees+1; tt++){
		if (tt == ntrees){
			for (int hh = 0; hh < kGGHSamples; hh++) foutput->WriteTObject(hhiggs[hh]);
		}
		foutput->WriteTObject(hfull[tt]);

    TH1F* hProjKD = (TH1F*)hfull[tt]->ProjectionY(Form("%s_py",hfull[tt]->GetName()),0,-1,"e");
		hProjKD->SetTitle(Form("%s Projection",hfull[tt]->GetTitle()));
		hproj[tt] = hProjKD;

		foutput->WriteTObject(hproj[tt]);

		for (int binx = 1; binx <= hfull[tt]->GetNbinsX(); binx++){
			TH1F* hSlice = (TH1F*) hProjKD->Clone( (binx<hfull[tt]->GetNbinsX() ? Form("%s_SidebandRegion%i",hproj[tt]->GetName(),binx) : Form("%s_SignalRegion",hproj[tt]->GetName()) ) );
      for (int biny = 1; biny <= hfull[tt]->GetNbinsY(); biny++){
        hSlice->SetBinContent(biny, hfull[tt]->GetBinContent(binx, biny)); hSlice->SetBinError(biny, hfull[tt]->GetBinError(binx, biny));
      }
			hSlice->SetTitle(Form("%s %s GeV Projection",hfull[tt]->GetTitle(),hfull[tt]->GetXaxis()->GetBinLabel(binx)));
			foutput->WriteTObject(hSlice);
			delete hSlice;
		}
	}
	for (int tt = 0; tt < ntrees_extra; tt++) delete tc_extra[tt];
	for (int tt = 0; tt < ntrees+1; tt++){
		delete hproj[tt];
		delete hfull[tt];
		if(tt<ntrees) delete tc[tt];
	}
  for (int bs=0; bs<2; bs++) delete tBeam[bs];

	foutput->Close();
  delete hRatio;
  fSSinput->Close();
}


void compare_KDShapeVariation_Untransformed(int iSyst=0, int normOption=0, int SR_noBkg=0, int removeSIP=0){
	TString OUTPUT_NAME = "LifetimeKD_ShapeSystVars_";
	if(SR_noBkg==1) OUTPUT_NAME.Append("NoBkgSROnly_");
	if(removeSIP==1) OUTPUT_NAME.Append("NoSIP_");
	else if(removeSIP==2) OUTPUT_NAME.Append("NewSIP_");
	if(normOption==1) OUTPUT_NAME.Append("AreaNormalized_");
	OUTPUT_NAME.Append(Form("%s_Comparison.root",cSystVariable[iSyst]));
	TString coutput_common = user_dir_hep + "Analysis/Plots/";
	coutput_common.Append(Form("%s/",cSystVariable[iSyst]));
	TString coutput = coutput_common + OUTPUT_NAME;
	TFile* foutput = new TFile(coutput, "recreate");
	const int ntrees = 2;
	char* chProj[ntrees]={"hdata_full_py", "hmc_full_py"};
	TH1F* hProj[nMZZ_ranges*ntrees][3] = { { 0 } };
	TH1F* hhiggs[kGGHSamples][3] = { { 0 } };
	TGraphAsymmErrors* tgdata[nMZZ_ranges*3];
	TFile* finput[3][2] = { { 0 } };
	double max_plot[nMZZ_ranges] = { 0 };
  bool received_ergtev[2] ={ 0 };

	for (int folder = 0; folder < 3; folder++){
		for (int erg_tev = 7; erg_tev < 9; erg_tev++){
			TString erg_dir;
			erg_dir.Form("LHC_%iTeV", erg_tev);
			TString comstring;
			comstring.Form("%iTeV", erg_tev);
			int EnergyIndex = 0; if (erg_tev == 8) EnergyIndex = 1;
			TString cinput_common = user_dir_hep + "Analysis/Plots/";
			cinput_common.Append(Form("%s/",cSystVariable[iSyst]));

			TString INPUT_NAME = "LifetimeKD_SystVariables_DataQQBZZ_";
			if(removeSIP==1) INPUT_NAME.Append("noSIP_");
			if(removeSIP==2) INPUT_NAME.Append("newSIP_");
			INPUT_NAME.Append(user_folder[folder]);
			INPUT_NAME.Append(Form("_%s_", cSystVariable[iSyst]));
			INPUT_NAME.Append(comstring);
			INPUT_NAME.Append(".root");
			cout << INPUT_NAME << endl;
			TString cinput = cinput_common + INPUT_NAME;
			cout << "Attempting to open " << cinput << endl;
			finput[folder][EnergyIndex] = new TFile(cinput, "read");

			if (finput[folder][EnergyIndex] == 0 || finput[folder][EnergyIndex]->IsZombie()) continue;
      else received_ergtev[EnergyIndex] = true;
			for (int sb = 0; sb < nMZZ_ranges; sb++){
				for (int hh = 0; hh < ntrees; hh++){
					TString hname = Form("%s_SidebandRegion%i", chProj[hh], sb + 1);
					if(sb==nMZZ_ranges-1) hname = Form("%s_SignalRegion", chProj[hh]);
					cout << hname << endl;
					TH1F* htemp = (TH1F*)finput[folder][EnergyIndex]->Get(hname);
					if(hh==1) htemp->Scale(luminosity[EnergyIndex]);
					if (hProj[nMZZ_ranges * hh + sb][folder] == 0){
						hProj[nMZZ_ranges*hh+sb][folder] = (TH1F*)htemp->Clone(Form("acc_%s", htemp->GetName()));
					}
					else hProj[nMZZ_ranges*hh+sb][folder]->Add(htemp);
					delete htemp;
				}
			}
			for (int hh = 0; hh < kGGHSamples; hh++){
				TString hname = Form("hhiggs_CTau%.0f_full", sample_CTau[hh]);
				TH1F* htemp = (TH1F*)finput[folder][EnergyIndex]->Get(hname);
				htemp->Scale(luminosity[EnergyIndex]);
				htemp->SetLineColor(kBlack);
				htemp->SetLineWidth(2);
				htemp->SetTitle("");
				if(hh==0) htemp->SetLineStyle(7);
				else if(hh==1) htemp->SetLineStyle(6);
				else if(hh==2) htemp->SetLineStyle(2);
				else if(hh==3) htemp->SetLineStyle(3);
				if (hhiggs[hh][folder] == 0 && htemp!=0){
					hhiggs[hh][folder] = (TH1F*)htemp->Clone(Form("acc_%s", htemp->GetName()));
					hhiggs[hh][folder]->SetLineColor(kBlack);
				}
				else if(htemp!=0) hhiggs[hh][folder]->Add(htemp);
				if(htemp!=0) delete htemp;
			}
		}

		for (int sb = 0; sb < nMZZ_ranges; sb++){
			if (normOption == 1){
				cout << "Scaling MC by 1./" << hProj[nMZZ_ranges + sb][folder]->Integral() << endl;
				hProj[nMZZ_ranges + sb][folder]->Scale(1. / hProj[nMZZ_ranges + sb][folder]->Integral());
				if (sb == nMZZ_ranges - 1){
					for (int hh = 0; hh < kGGHSamples; hh++) hhiggs[hh][folder]->Scale(1. / hhiggs[hh][folder]->Integral());
				}
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
			if(normOption==1) integral_data = hProj[sb][folder]->Integral();

			for (int bin = 1; bin <= hProj[sb][folder]->GetNbinsX(); bin++){
				double bincenter = hProj[sb][folder]->GetBinCenter(bin);
				double bincontent = hProj[sb][folder]->GetBinContent(bin);

				if (bincontent > 0){
					xx_data[ndata] = bincenter;
					yy_data[ndata] = bincontent / integral_data;
					xu_data[ndata] = 0;
					xd_data[ndata] = 0;
					yu_data[ndata] = (ROOT::Math::chisquared_quantile_c(quant, 2 * (bincontent + 1)) / 2. - bincontent) / integral_data;
					yd_data[ndata] = ((bincontent == 0) ? 0 : (bincontent - ROOT::Math::chisquared_quantile_c(1 - quant, 2 * bincontent) / 2.)) / integral_data;

					double high_data = yy_data[ndata] + yu_data[ndata];
					if (high_data > max_plot[sb] && sb<(nMZZ_ranges - 1)) max_plot[sb] = high_data; // TO BE REVISED LATER
					ndata++;
				}
			}
			tgdata[nMZZ_ranges * folder + sb] = new TGraphAsymmErrors(ndata, xx_data, yy_data, xd_data, xu_data, yd_data, yu_data);
			tgdata[nMZZ_ranges * folder + sb]->SetName(Form("tgdata_%s_SB%i", user_folder[folder], sb + 1));
			if(folder==2) tgdata[nMZZ_ranges * folder + sb]->SetMarkerSize(1.2);
			if(folder==1) tgdata[nMZZ_ranges * folder + sb]->SetMarkerSize(1.4);
			if(folder==0) tgdata[nMZZ_ranges * folder + sb]->SetMarkerSize(1.5);
			if(folder==2) tgdata[nMZZ_ranges * folder + sb]->SetMarkerStyle(20);
			if(folder==1) tgdata[nMZZ_ranges * folder + sb]->SetMarkerStyle(33);
			if(folder==0) tgdata[nMZZ_ranges * folder + sb]->SetMarkerStyle(34);
			tgdata[nMZZ_ranges * folder + sb]->SetMarkerColor(kBlack);
			tgdata[nMZZ_ranges * folder + sb]->SetLineColor(kBlack);
			tgdata[nMZZ_ranges * folder + sb]->SetLineWidth(1);

			if(normOption==1) hProj[sb][folder]->Scale(1. / hProj[sb][folder]->Integral());
		}
	}

	gStyle->SetTitleFont(62, "t");
	gROOT->SetStyle(gStyle->GetName());
	gROOT->ForceStyle();

	for (int sb = 0; sb < nMZZ_ranges; sb++){
		foutput->cd();
		if (sb < nMZZ_ranges - 1 && iSyst >= 7 && iSyst <= 11) continue;
		if (SR_noBkg == 1 && sb < nMZZ_ranges - 1) continue;
//		if (removeSIP >= 1 && sb < nMZZ_ranges - 1) continue;
		for (int folder = 0; folder < 3; folder++){
			max_plot[sb] = TMath::Max(max_plot[sb], hProj[nMZZ_ranges + sb][folder]->GetMaximum());
			if (sb == nMZZ_ranges - 1){
				for (int hh = 0; hh < kGGHSamples; hh++) max_plot[sb] = TMath::Max(max_plot[sb], hhiggs[hh][folder]->GetMaximum());
			}
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
    if (received_ergtev[0] && received_ergtev[1]) text = pt->AddText(0.537, 0.45, "#font[42]{19.7 fb^{-1} (8 TeV) + 5.1 fb^{-1} (7 TeV)}");
    else if (received_ergtev[1]) text = pt->AddText(0.837, 0.45, "#font[42]{19.7 fb^{-1} (8 TeV)}");
    else if (received_ergtev[0]) text = pt->AddText(0.837, 0.45, "#font[42]{ 5.1 fb^{-1} (7 TeV)}");
		text->SetTextSize(0.0315);

		tgdata[nMZZ_ranges * 2 + sb]->SetLineColor(kOrange + 10);
		tgdata[nMZZ_ranges * 2 + sb]->SetMarkerColor(kOrange + 10);
		hProj[nMZZ_ranges * 1 + sb][2]->SetLineColor(kOrange + 10);
		hProj[nMZZ_ranges * 1 + sb][2]->SetLineWidth(2);
		hProj[nMZZ_ranges * 1 + sb][2]->SetLineStyle(1);
		tgdata[nMZZ_ranges * 1 + sb]->SetLineColor(kBlue);
		tgdata[nMZZ_ranges * 1 + sb]->SetMarkerColor(kBlue);
		hProj[nMZZ_ranges * 1 + sb][1]->SetLineColor(kBlue);
		hProj[nMZZ_ranges * 1 + sb][1]->SetLineWidth(2);
		hProj[nMZZ_ranges * 1 + sb][1]->SetLineStyle(1);
		tgdata[nMZZ_ranges * 0 + sb]->SetLineColor(TColor::GetColor("#669966"));
		tgdata[nMZZ_ranges * 0 + sb]->SetMarkerColor(TColor::GetColor("#669966"));
		hProj[nMZZ_ranges * 1 + sb][0]->SetLineColor(TColor::GetColor("#669966"));
		hProj[nMZZ_ranges * 1 + sb][0]->SetLineWidth(2);
		hProj[nMZZ_ranges * 1 + sb][0]->SetLineStyle(1);
		for (int folder = 0; folder < 3; folder++){
			for (int hh = 0; hh < ntrees; hh++){
				hProj[nMZZ_ranges * hh + sb][folder]->SetTitle("");
				tgdata[nMZZ_ranges * folder + sb]->SetTitle("");
			}
		}

		foutput->cd();
		TString canvasname_2D = "cCompare_DataMC_AllChannels_VtxSystVar_";
		if (removeSIP == 1) canvasname_2D.Append("NoSIP_");
		if (removeSIP == 2) canvasname_2D.Append("NewSIP_");
		if (SR_noBkg == 1) canvasname_2D.Append("NoBkg_");
		if (sb < nMZZ_ranges - 1) canvasname_2D.Append(Form("%s_SB%i", cSystVariable[iSyst], sb + 1));
		else canvasname_2D.Append(Form("%s_SR", cSystVariable[iSyst]));
		if (normOption == 1) canvasname_2D.Append("_AreaNormalized");
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

		TLegend *l2D = new TLegend(0.20, 0.57, 0.58, 0.90);
		l2D->SetBorderSize(0);
		l2D->SetTextFont(42);
		l2D->SetTextSize(0.03);
		l2D->SetLineColor(1);
		l2D->SetLineStyle(1);
		l2D->SetLineWidth(1);
		l2D->SetFillColor(0);
		l2D->SetFillStyle(0);

		TString yTitle = (normOption == 0 ? Form("Events / %.3f", hProj[nMZZ_ranges * 1 + sb][0]->GetXaxis()->GetBinWidth(1)) : Form("Rate / %.3f", hProj[nMZZ_ranges * 1 + sb][0]->GetXaxis()->GetBinWidth(1)));
		if (iSyst == 3) yTitle = (normOption == 0 ? Form("Events / %.2e", hProj[nMZZ_ranges * 1 + sb][0]->GetXaxis()->GetBinWidth(1)) : Form("Rate / %.2e", hProj[nMZZ_ranges * 1 + sb][0]->GetXaxis()->GetBinWidth(1)));
		TString xTitle = cSystVariable_label[iSyst];

		hProj[nMZZ_ranges * 1 + sb][0]->GetXaxis()->SetTitle(xTitle);
		hProj[nMZZ_ranges * 1 + sb][0]->GetYaxis()->SetTitle(yTitle);
		hProj[nMZZ_ranges * 1 + sb][0]->GetYaxis()->SetRangeUser(0, max_plot[sb] * 1.25);
		if (iSyst == 4 && sb < 2) hProj[nMZZ_ranges * 1 + sb][0]->GetYaxis()->SetRangeUser(0, max_plot[sb] * 1.75);
		if (iSyst == 3) hProj[nMZZ_ranges * 1 + sb][0]->GetYaxis()->SetRangeUser(0, max_plot[sb] * 1.85);
		if (iSyst == 2) hProj[nMZZ_ranges * 1 + sb][0]->GetYaxis()->SetRangeUser(0, max_plot[sb] * 1.75);
		hProj[nMZZ_ranges * 1 + sb][0]->GetXaxis()->SetNdivisions(505);
		hProj[nMZZ_ranges * 1 + sb][0]->GetXaxis()->SetLabelFont(42);
		hProj[nMZZ_ranges * 1 + sb][0]->GetXaxis()->SetLabelOffset(0.007);
		hProj[nMZZ_ranges * 1 + sb][0]->GetXaxis()->SetLabelSize(0.04);
		hProj[nMZZ_ranges * 1 + sb][0]->GetXaxis()->SetTitleSize(0.06);
		hProj[nMZZ_ranges * 1 + sb][0]->GetXaxis()->SetTitleOffset(0.9);
		hProj[nMZZ_ranges * 1 + sb][0]->GetXaxis()->SetTitleFont(42);
		hProj[nMZZ_ranges * 1 + sb][0]->GetYaxis()->SetNdivisions(505);
		hProj[nMZZ_ranges * 1 + sb][0]->GetYaxis()->SetLabelFont(42);
		hProj[nMZZ_ranges * 1 + sb][0]->GetYaxis()->SetLabelOffset(0.007);
		hProj[nMZZ_ranges * 1 + sb][0]->GetYaxis()->SetLabelSize(0.04);
		hProj[nMZZ_ranges * 1 + sb][0]->GetYaxis()->SetTitleSize(0.06);
		hProj[nMZZ_ranges * 1 + sb][0]->GetYaxis()->SetTitleOffset(1.1);
		hProj[nMZZ_ranges * 1 + sb][0]->GetYaxis()->SetTitleFont(42);
		hProj[nMZZ_ranges * 1 + sb][0]->GetXaxis()->SetRangeUser(visualRange_SystVariable[iSyst][0], visualRange_SystVariable[iSyst][1]);
		if (sb == nMZZ_ranges - 1 && SR_noBkg == 1){
			hhiggs[0][0]->GetXaxis()->SetTitle(xTitle);
			hhiggs[0][0]->GetYaxis()->SetTitle(yTitle);
			hhiggs[0][0]->GetYaxis()->SetRangeUser(0, max_plot[sb] * 1.25);
			if (iSyst == 3) hhiggs[0][0]->GetYaxis()->SetRangeUser(0, max_plot[sb] * 1.85);
			if (iSyst == 2) hhiggs[0][0]->GetYaxis()->SetRangeUser(0, max_plot[sb] * 1.75);
			hhiggs[0][0]->GetXaxis()->SetNdivisions(505);
			hhiggs[0][0]->GetXaxis()->SetLabelFont(42);
			hhiggs[0][0]->GetXaxis()->SetLabelOffset(0.007);
			hhiggs[0][0]->GetXaxis()->SetLabelSize(0.04);
			hhiggs[0][0]->GetXaxis()->SetTitleSize(0.06);
			hhiggs[0][0]->GetXaxis()->SetTitleOffset(0.9);
			hhiggs[0][0]->GetXaxis()->SetTitleFont(42);
			hhiggs[0][0]->GetYaxis()->SetNdivisions(505);
			hhiggs[0][0]->GetYaxis()->SetLabelFont(42);
			hhiggs[0][0]->GetYaxis()->SetLabelOffset(0.007);
			hhiggs[0][0]->GetYaxis()->SetLabelSize(0.04);
			hhiggs[0][0]->GetYaxis()->SetTitleSize(0.06);
			hhiggs[0][0]->GetYaxis()->SetTitleOffset(1.1);
			hhiggs[0][0]->GetYaxis()->SetTitleFont(42);
			hhiggs[0][0]->GetXaxis()->SetRangeUser(visualRange_SystVariable[iSyst][0], visualRange_SystVariable[iSyst][1]);
		}

		TH1F* higgsclone[kGGHSamples];
		if (sb < nMZZ_ranges - 1){
			l2D->AddEntry(tgdata[nMZZ_ranges * 2 + sb], "Observed 2e2#mu", "ep");
			l2D->AddEntry(tgdata[nMZZ_ranges * 1 + sb], "Observed 4e", "ep");
			l2D->AddEntry(tgdata[nMZZ_ranges * 0 + sb], "Observed 4#mu", "ep");
		}
		else{
			for (int hh = 0; hh < kGGHSamples; hh++){
				if (hhiggs[hh][0] != 0) higgsclone[hh] = (TH1F*)hhiggs[hh][0]->Clone(Form("%s_CLONED", hhiggs[hh][0]->GetName()));
				higgsclone[hh]->SetLineColor(kBlack);
			}
			if(SR_noBkg == 1) higgsclone[0]->SetLineStyle(1);
      l2D->AddEntry(higgsclone[0], "Higgs c#tau_{H}=0#mum", "l");
			if (SR_noBkg == 1){
					l2D->AddEntry(higgsclone[1], "Higgs c#tau_{H}=100 #mum", "l");
					l2D->AddEntry(higgsclone[2], "Higgs c#tau_{H}=500 #mum", "l");
					l2D->AddEntry(higgsclone[3], "Higgs c#tau_{H}=1000 #mum", "l");
			}
      else if (normOption==0){
        l2D->AddEntry(higgsclone[1], "Higgs c#tau_{H}=100 #mum", "l");
      }
		}
		if (SR_noBkg == 0){
			l2D->AddEntry(hProj[nMZZ_ranges * 1 + sb][2], "Bkg. 2e2#mu", "l");
			l2D->AddEntry(hProj[nMZZ_ranges * 1 + sb][1], "Bkg. 4e", "l");
			l2D->AddEntry(hProj[nMZZ_ranges * 1 + sb][0], "Bkg. 4#mu", "l");
			hProj[nMZZ_ranges * 1 + sb][0]->Draw("hist");
			hProj[nMZZ_ranges * 1 + sb][1]->Draw("histsame");
			hProj[nMZZ_ranges * 1 + sb][2]->Draw("histsame");
			l2D->Draw("same");
			if (sb < nMZZ_ranges - 1){
				tgdata[nMZZ_ranges * 2 + sb]->Draw("e1psame");
				tgdata[nMZZ_ranges * 0 + sb]->Draw("e1psame");
				tgdata[nMZZ_ranges * 1 + sb]->Draw("e1psame");
			}
			else{
				cout << "Drawing Higgs MC now" << endl;
        for (int folder = 0; folder < 3; folder++){
          hhiggs[0][folder]->SetLineColor(hProj[nMZZ_ranges * 1 + sb][folder]->GetLineColor());
          hhiggs[0][folder]->Draw("histsame");
          if (normOption==0){
            hhiggs[1][folder]->SetLineColor(hProj[nMZZ_ranges * 1 + sb][folder]->GetLineColor());
            hhiggs[1][folder]->Draw("histsame");
          }
        }
			}
		}
		else if (sb == nMZZ_ranges - 1){
			cout << "Drawing Higgs MC now" << endl;
			for (int hh = 0; hh < kGGHSamples; hh++){
				for (int folder = 0; folder < 3; folder++){
					if (hhiggs[hh][folder] == 0){
						cout << "Skipping draw for " << hh << "\t" << folder << endl;
						continue;
					}
					hhiggs[hh][folder]->SetLineColor(hProj[nMZZ_ranges * 1 + sb][folder]->GetLineColor());
					if(hh==0) hhiggs[hh][folder]->SetLineStyle(1);
					if(hh==0 && folder==0) hhiggs[hh][folder]->Draw("hist");
					else  hhiggs[hh][folder]->Draw("histsame");
				}
			}
			l2D->Draw("same");
		}
	

		TF1* guideFunction;
		if (normOption == 1 && (iSyst == 5 || iSyst == 6)){
			guideFunction = new TF1("guideFunction","exp(-pow(x,2)/(2*[0]))/sqrt(2*TMath::Pi()*[0])*[1]",visualRange_SystVariable[iSyst][0],visualRange_SystVariable[iSyst][1]);
			guideFunction->SetParameters(1,hProj[nMZZ_ranges * 1 + sb][0]->GetXaxis()->GetBinWidth(1));
			guideFunction->SetLineColor(kBlack);
			guideFunction->SetMarkerColor(kBlack);
			guideFunction->SetLineWidth(2);
			guideFunction->SetLineStyle(7);
			guideFunction->Draw("csame");
		}
		pt->Draw();

		double strSRlow = systZZMass_range[sb][0];
		double strSRhigh = systZZMass_range[sb][1];
		int strIntSRlow = (int) systZZMass_range[sb][0];
		int strIntSRhigh = (int) systZZMass_range[sb][1];
		TString strSidebandRegion;
		if(strSRlow>strIntSRlow) strSidebandRegion = Form("%.1f < m_{4l}",systZZMass_range[sb][0]);
		else strSidebandRegion = Form("%.0f < m_{4l}",systZZMass_range[sb][0]);
		if(strSRhigh>strIntSRhigh) strSidebandRegion.Append(Form(" < %.1f GeV",systZZMass_range[sb][1]));
		else strSidebandRegion.Append(Form(" < %.0f GeV",systZZMass_range[sb][1]));

		TPaveText *pt10 = new TPaveText(0.61,0.84,0.90,0.92,"brNDC");
		pt10->SetBorderSize(0);
		pt10->SetTextAlign(12);
		pt10->SetTextSize(0.03);
		pt10->SetFillStyle(0);
		pt10->SetTextFont(42);
		TText* text10 = pt10->AddText(0.01,0.01,strSidebandRegion);
		pt10->Draw();

    TString strSIPcut = "";
    if (removeSIP == 1) strSIPcut = "No cut";
    else if (removeSIP == 2) strSIPcut = "New cut";
    TPaveText* pt20;
    if (removeSIP<2) pt20 = new TPaveText(0.82, 0.80, 0.90, 0.84, "brNDC");
    if (removeSIP==2) pt20 = new TPaveText(0.81, 0.80, 0.90, 0.84, "brNDC");
    pt20->SetBorderSize(0);
    pt20->SetTextAlign(12);
    pt20->SetTextSize(0.03);
    pt20->SetFillStyle(0);
    pt20->SetTextFont(42);
    TText* text20 = pt20->AddText(0.01, 0.01, strSIPcut);
    if(removeSIP>0) pt20->Draw();


		c2D->RedrawAxis();
		c2D->Modified();
		c2D->Update();
		foutput->WriteTObject(c2D);

		TString canvasDir = coutput_common;
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

		if (normOption == 1 && (iSyst == 5 || iSyst == 6)) delete guideFunction;
		delete l2D;
		c2D->Close();
		delete pt;

		for (int folder = 0; folder < 3; folder++){
			if (normOption == 1){
				hProj[nMZZ_ranges * 1 + sb][folder]->SetName(Form("%s_%s",hProj[nMZZ_ranges * 1 + sb][folder]->GetName(),user_folder[folder]));
				tgdata[nMZZ_ranges * folder + sb]->SetName(Form("%s_%s",tgdata[nMZZ_ranges * folder + sb]->GetName(),user_folder[folder]));
				foutput->WriteTObject(hProj[nMZZ_ranges * 1 + sb][folder]);
				foutput->WriteTObject(tgdata[nMZZ_ranges * folder + sb]);
			}
			for (int hh = 0; hh < ntrees; hh++) delete hProj[nMZZ_ranges * hh + sb][folder];
			delete tgdata[nMZZ_ranges * folder + sb];
		}
	}

	for (int folder = 0; folder < 3; folder++){
		for (int e = 0; e < 2; e++){
			if (finput[folder][e] == 0 || finput[folder][e]->IsZombie()) continue;
			else finput[folder][e]->Close();
		}
	}
	foutput->Close();
}

void compare_KDShapeVariation_Untransformed_perChannel(int folder, int iSyst=0, int normOption=0, int SR_noBkg=0, int removeSIP=0){
  TString OUTPUT_NAME = "LifetimeKD_ShapeSystVars_";
  OUTPUT_NAME.Append(user_folder[folder]);
  if (SR_noBkg==1) OUTPUT_NAME.Append("NoBkgSROnly_");
  if (removeSIP==1) OUTPUT_NAME.Append("NoSIP_");
  else if (removeSIP==2) OUTPUT_NAME.Append("NewSIP_");
  if (normOption==1) OUTPUT_NAME.Append("AreaNormalized_");
  OUTPUT_NAME.Append(Form("%s_Comparison.root", cSystVariable[iSyst]));
  TString coutput_common = user_dir_hep + "Analysis/Plots/";
  coutput_common.Append(Form("%s/", cSystVariable[iSyst]));
  TString coutput = coutput_common + OUTPUT_NAME;
  TFile* foutput = new TFile(coutput, "recreate");
  const int ntrees = 2;
  char* chProj[ntrees]={ "hdata_full_py", "hmc_full_py" };
  TH1F* hProj[nMZZ_ranges*ntrees][3] ={ { 0 } };
  TH1F* hhiggs[kGGHSamples][3] ={ { 0 } };
  TGraphAsymmErrors* tgdata[nMZZ_ranges*3];
  TFile* finput[3][2] ={ { 0 } };
  double max_plot[nMZZ_ranges] ={ 0 };
  bool received_ergtev[2] = { 0 };

  for (int erg_tev = 7; erg_tev < 9; erg_tev++){
    TString erg_dir;
    erg_dir.Form("LHC_%iTeV", erg_tev);
    TString comstring;
    comstring.Form("%iTeV", erg_tev);
    int EnergyIndex = 0; if (erg_tev == 8) EnergyIndex = 1;
    TString cinput_common = user_dir_hep + "Analysis/Plots/";
    cinput_common.Append(Form("%s/", cSystVariable[iSyst]));

    TString INPUT_NAME = "LifetimeKD_SystVariables_DataQQBZZ_";
    if (removeSIP==1) INPUT_NAME.Append("noSIP_");
    if (removeSIP==2) INPUT_NAME.Append("newSIP_");
    INPUT_NAME.Append(user_folder[folder]);
    INPUT_NAME.Append(Form("_%s_", cSystVariable[iSyst]));
    INPUT_NAME.Append(comstring);
    INPUT_NAME.Append(".root");
    cout << INPUT_NAME << endl;
    TString cinput = cinput_common + INPUT_NAME;
    cout << "Attempting to open " << cinput << endl;
    finput[folder][EnergyIndex] = new TFile(cinput, "read");

    if (finput[folder][EnergyIndex] == 0 || finput[folder][EnergyIndex]->IsZombie()) continue;
    else received_ergtev[EnergyIndex]=true;
    for (int sb = 0; sb < nMZZ_ranges; sb++){
      for (int hh = 0; hh < ntrees; hh++){
        TString hname = Form("%s_SidebandRegion%i", chProj[hh], sb + 1);
        if (sb==nMZZ_ranges-1) hname = Form("%s_SignalRegion", chProj[hh]);
        cout << hname << endl;
        TH1F* htemp = (TH1F*)finput[folder][EnergyIndex]->Get(hname);
        if (hh==1) htemp->Scale(luminosity[EnergyIndex]);
        if (hProj[nMZZ_ranges * hh + sb][folder] == 0){
          hProj[nMZZ_ranges*hh+sb][folder] = (TH1F*)htemp->Clone(Form("acc_%s", htemp->GetName()));
        }
        else hProj[nMZZ_ranges*hh+sb][folder]->Add(htemp);
        delete htemp;
      }
    }
    for (int hh = 0; hh < kGGHSamples; hh++){
      TString hname = Form("hhiggs_CTau%.0f_full", sample_CTau[hh]);
      TH1F* htemp = (TH1F*)finput[folder][EnergyIndex]->Get(hname);
      htemp->Scale(luminosity[EnergyIndex]);
      htemp->SetLineColor(kBlack);
      htemp->SetLineWidth(2);
      htemp->SetTitle("");
      if (hh==0) htemp->SetLineStyle(7);
      else if (hh==1) htemp->SetLineStyle(6);
      else if (hh==2) htemp->SetLineStyle(2);
      else if (hh==3) htemp->SetLineStyle(3);
      if (hhiggs[hh][folder] == 0 && htemp!=0){
        hhiggs[hh][folder] = (TH1F*)htemp->Clone(Form("acc_%s", htemp->GetName()));
        hhiggs[hh][folder]->SetLineColor(kBlack);
      }
      else if (htemp!=0) hhiggs[hh][folder]->Add(htemp);
      if (htemp!=0) delete htemp;
    }


    for (int sb = 0; sb < nMZZ_ranges; sb++){
      if (normOption == 1){
        cout << "Scaling MC by 1./" << hProj[nMZZ_ranges + sb][folder]->Integral() << endl;
        hProj[nMZZ_ranges + sb][folder]->Scale(1. / hProj[nMZZ_ranges + sb][folder]->Integral());
        if (sb == nMZZ_ranges - 1){
          for (int hh = 0; hh < kGGHSamples; hh++) hhiggs[hh][folder]->Scale(1. / hhiggs[hh][folder]->Integral());
        }
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
      if (normOption==1) integral_data = hProj[sb][folder]->Integral();

      for (int bin = 1; bin <= hProj[sb][folder]->GetNbinsX(); bin++){
        double bincenter = hProj[sb][folder]->GetBinCenter(bin);
        double bincontent = hProj[sb][folder]->GetBinContent(bin);

        if (bincontent > 0){
          xx_data[ndata] = bincenter;
          yy_data[ndata] = bincontent / integral_data;
          xu_data[ndata] = 0;
          xd_data[ndata] = 0;
          yu_data[ndata] = (ROOT::Math::chisquared_quantile_c(quant, 2 * (bincontent + 1)) / 2. - bincontent) / integral_data;
          yd_data[ndata] = ((bincontent == 0) ? 0 : (bincontent - ROOT::Math::chisquared_quantile_c(1 - quant, 2 * bincontent) / 2.)) / integral_data;

          double high_data = yy_data[ndata] + yu_data[ndata];
          if (high_data > max_plot[sb] && sb<(nMZZ_ranges - 1)) max_plot[sb] = high_data; // TO BE REVISED LATER
          ndata++;
        }
      }
      tgdata[nMZZ_ranges * folder + sb] = new TGraphAsymmErrors(ndata, xx_data, yy_data, xd_data, xu_data, yd_data, yu_data);
      tgdata[nMZZ_ranges * folder + sb]->SetName(Form("tgdata_%s_SB%i", user_folder[folder], sb + 1));
      tgdata[nMZZ_ranges * folder + sb]->SetMarkerSize(1.2);
      tgdata[nMZZ_ranges * folder + sb]->SetMarkerStyle(20);
      tgdata[nMZZ_ranges * folder + sb]->SetMarkerColor(kBlack);
      tgdata[nMZZ_ranges * folder + sb]->SetLineColor(kBlack);
      tgdata[nMZZ_ranges * folder + sb]->SetLineWidth(1);

      if (normOption==1) hProj[sb][folder]->Scale(1. / hProj[sb][folder]->Integral());
    }
  }

  gStyle->SetTitleFont(62, "t");
  gROOT->SetStyle(gStyle->GetName());
  gROOT->ForceStyle();

  for (int sb = 0; sb < nMZZ_ranges; sb++){
    foutput->cd();
    if (sb < nMZZ_ranges - 1 && iSyst >= 7 && iSyst <= 11) continue;
    if (SR_noBkg == 1 && sb < nMZZ_ranges - 1) continue;
    //		if (removeSIP >= 1 && sb < nMZZ_ranges - 1) continue;
    max_plot[sb] = TMath::Max(max_plot[sb], hProj[nMZZ_ranges + sb][folder]->GetMaximum());
    if (sb == nMZZ_ranges - 1){
      for (int hh = 0; hh < kGGHSamples; hh++) max_plot[sb] = TMath::Max(max_plot[sb], hhiggs[hh][folder]->GetMaximum());
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
    if (received_ergtev[0] && received_ergtev[1]) text = pt->AddText(0.537, 0.45, "#font[42]{19.7 fb^{-1} (8 TeV) + 5.1 fb^{-1} (7 TeV)}");
    else if (received_ergtev[1]) text = pt->AddText(0.837, 0.45, "#font[42]{19.7 fb^{-1} (8 TeV)}");
    else if (received_ergtev[0]) text = pt->AddText(0.837, 0.45, "#font[42]{ 5.1 fb^{-1} (7 TeV)}");
    text->SetTextSize(0.0315);

    tgdata[nMZZ_ranges * folder + sb]->SetLineColor(kOrange + 10);
    tgdata[nMZZ_ranges * folder + sb]->SetMarkerColor(kOrange + 10);
    hProj[nMZZ_ranges * 1 + sb][folder]->SetLineColor(kOrange + 10);
    hProj[nMZZ_ranges * 1 + sb][folder]->SetLineWidth(2);
    hProj[nMZZ_ranges * 1 + sb][folder]->SetLineStyle(1);
    for (int hh = 0; hh < ntrees; hh++){
      hProj[nMZZ_ranges * hh + sb][folder]->SetTitle("");
      tgdata[nMZZ_ranges * folder + sb]->SetTitle("");
    }


    foutput->cd();
    TString canvasname_2D = "cCompare_DataMC_AllChannels_VtxSystVar_";
    canvasname_2D.Append(user_folder[folder]);
    canvasname_2D.Append("_");
    if (removeSIP == 1) canvasname_2D.Append("NoSIP_");
    if (removeSIP == 2) canvasname_2D.Append("NewSIP_");
    if (SR_noBkg == 1) canvasname_2D.Append("NoBkg_");
    if (sb < nMZZ_ranges - 1) canvasname_2D.Append(Form("%s_SB%i", cSystVariable[iSyst], sb + 1));
    else canvasname_2D.Append(Form("%s_SR", cSystVariable[iSyst]));
    if (normOption == 1) canvasname_2D.Append("_AreaNormalized");
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

    TLegend *l2D = new TLegend(0.20, 0.57, 0.58, 0.90);
    l2D->SetBorderSize(0);
    l2D->SetTextFont(42);
    l2D->SetTextSize(0.03);
    l2D->SetLineColor(1);
    l2D->SetLineStyle(1);
    l2D->SetLineWidth(1);
    l2D->SetFillColor(0);
    l2D->SetFillStyle(0);

    TString yTitle = (normOption == 0 ? Form("Events / %.3f", hProj[nMZZ_ranges * 1 + sb][folder]->GetXaxis()->GetBinWidth(1)) : Form("Rate / %.3f", hProj[nMZZ_ranges * 1 + sb][folder]->GetXaxis()->GetBinWidth(1)));
    if (iSyst == 3) yTitle = (normOption == 0 ? Form("Events / %.2e", hProj[nMZZ_ranges * 1 + sb][folder]->GetXaxis()->GetBinWidth(1)) : Form("Rate / %.2e", hProj[nMZZ_ranges * 1 + sb][folder]->GetXaxis()->GetBinWidth(1)));
    TString xTitle = cSystVariable_label[iSyst];

    hProj[nMZZ_ranges * 1 + sb][folder]->GetXaxis()->SetTitle(xTitle);
    hProj[nMZZ_ranges * 1 + sb][folder]->GetYaxis()->SetTitle(yTitle);
    hProj[nMZZ_ranges * 1 + sb][folder]->GetYaxis()->SetRangeUser(0, max_plot[sb] * 1.25);
    if (iSyst == 4 && sb < 2) hProj[nMZZ_ranges * 1 + sb][folder]->GetYaxis()->SetRangeUser(0, max_plot[sb] * 1.75);
    if (iSyst == 3) hProj[nMZZ_ranges * 1 + sb][folder]->GetYaxis()->SetRangeUser(0, max_plot[sb] * 1.85);
    if (iSyst == 2) hProj[nMZZ_ranges * 1 + sb][folder]->GetYaxis()->SetRangeUser(0, max_plot[sb] * 1.75);
    hProj[nMZZ_ranges * 1 + sb][folder]->GetXaxis()->SetNdivisions(505);
    hProj[nMZZ_ranges * 1 + sb][folder]->GetXaxis()->SetLabelFont(42);
    hProj[nMZZ_ranges * 1 + sb][folder]->GetXaxis()->SetLabelOffset(0.007);
    hProj[nMZZ_ranges * 1 + sb][folder]->GetXaxis()->SetLabelSize(0.04);
    hProj[nMZZ_ranges * 1 + sb][folder]->GetXaxis()->SetTitleSize(0.06);
    hProj[nMZZ_ranges * 1 + sb][folder]->GetXaxis()->SetTitleOffset(0.9);
    hProj[nMZZ_ranges * 1 + sb][folder]->GetXaxis()->SetTitleFont(42);
    hProj[nMZZ_ranges * 1 + sb][folder]->GetYaxis()->SetNdivisions(505);
    hProj[nMZZ_ranges * 1 + sb][folder]->GetYaxis()->SetLabelFont(42);
    hProj[nMZZ_ranges * 1 + sb][folder]->GetYaxis()->SetLabelOffset(0.007);
    hProj[nMZZ_ranges * 1 + sb][folder]->GetYaxis()->SetLabelSize(0.04);
    hProj[nMZZ_ranges * 1 + sb][folder]->GetYaxis()->SetTitleSize(0.06);
    hProj[nMZZ_ranges * 1 + sb][folder]->GetYaxis()->SetTitleOffset(1.1);
    hProj[nMZZ_ranges * 1 + sb][folder]->GetYaxis()->SetTitleFont(42);
    hProj[nMZZ_ranges * 1 + sb][folder]->GetXaxis()->SetRangeUser(visualRange_SystVariable[iSyst][0], visualRange_SystVariable[iSyst][1]);

    if (sb == nMZZ_ranges - 1 && SR_noBkg == 1){
      hhiggs[0][folder]->GetXaxis()->SetTitle(xTitle);
      hhiggs[0][folder]->GetYaxis()->SetTitle(yTitle);
      hhiggs[0][folder]->GetYaxis()->SetRangeUser(0, max_plot[sb] * 1.25);
      if (iSyst == 3) hhiggs[0][folder]->GetYaxis()->SetRangeUser(0, max_plot[sb] * 1.85);
      if (iSyst == 2) hhiggs[0][folder]->GetYaxis()->SetRangeUser(0, max_plot[sb] * 1.75);
      hhiggs[0][folder]->GetXaxis()->SetNdivisions(505);
      hhiggs[0][folder]->GetXaxis()->SetLabelFont(42);
      hhiggs[0][folder]->GetXaxis()->SetLabelOffset(0.007);
      hhiggs[0][folder]->GetXaxis()->SetLabelSize(0.04);
      hhiggs[0][folder]->GetXaxis()->SetTitleSize(0.06);
      hhiggs[0][folder]->GetXaxis()->SetTitleOffset(0.9);
      hhiggs[0][folder]->GetXaxis()->SetTitleFont(42);
      hhiggs[0][folder]->GetYaxis()->SetNdivisions(505);
      hhiggs[0][folder]->GetYaxis()->SetLabelFont(42);
      hhiggs[0][folder]->GetYaxis()->SetLabelOffset(0.007);
      hhiggs[0][folder]->GetYaxis()->SetLabelSize(0.04);
      hhiggs[0][folder]->GetYaxis()->SetTitleSize(0.06);
      hhiggs[0][folder]->GetYaxis()->SetTitleOffset(1.1);
      hhiggs[0][folder]->GetYaxis()->SetTitleFont(42);
      hhiggs[0][folder]->GetXaxis()->SetRangeUser(visualRange_SystVariable[iSyst][0], visualRange_SystVariable[iSyst][1]);
    }

    string folder_id;
    if (folder==0) folder_id = "4#mu";
    if (folder==1) folder_id = "4e";
    if (folder==2) folder_id = "2e2#mu";
    TH1F* higgsclone[kGGHSamples];
    if (sb < nMZZ_ranges - 1){
      l2D->AddEntry(tgdata[nMZZ_ranges * folder + sb], Form("Observed %s", folder_id.c_str()), "ep");
    }
    else{
      for (int hh = 0; hh < kGGHSamples; hh++){
        if (hhiggs[hh][folder] != 0) higgsclone[hh] = (TH1F*)hhiggs[hh][folder]->Clone(Form("%s_CLONED", hhiggs[hh][folder]->GetName()));
        higgsclone[hh]->SetLineColor(kBlack);
      }
      if (SR_noBkg == 1) higgsclone[0]->SetLineStyle(1);
      l2D->AddEntry(higgsclone[0], "Higgs c#tau_{H}=0 #mum", "l");
      if (SR_noBkg == 1){
        higgsclone[1]->SetLineStyle(1);
        higgsclone[1]->SetLineColor(kBlue);
        l2D->AddEntry(higgsclone[1], "Higgs c#tau_{H}=100 #mum", "l");
        higgsclone[2]->SetLineStyle(1);
        higgsclone[2]->SetLineColor(kGreen+2);
//        l2D->AddEntry(higgsclone[2], "Higgs c#tau_{H}=500 #mum", "l");
        higgsclone[3]->SetLineStyle(1);
        higgsclone[3]->SetLineColor(kRed);
        l2D->AddEntry(higgsclone[3], "Higgs c#tau_{H}=1000 #mum", "l");
      }
      else if (normOption==0){
        l2D->AddEntry(higgsclone[1], "Higgs c#tau_{H}=100 #mum", "l");
      }
    }
    if (SR_noBkg == 0){
      l2D->AddEntry(hProj[nMZZ_ranges * 1 + sb][folder], Form("Bkg. %s", folder_id.c_str()), "l");
      hProj[nMZZ_ranges * 1 + sb][folder]->Draw("hist");
      l2D->Draw("same");
      if (sb < nMZZ_ranges - 1){
        tgdata[nMZZ_ranges * folder + sb]->Draw("e1psame");
      }
      else{
        cout << "Drawing Higgs MC now" << endl;
        hhiggs[0][folder]->SetLineColor(hProj[nMZZ_ranges * 1 + sb][folder]->GetLineColor());
        hhiggs[0][folder]->Draw("histsame");
        if (normOption==0){
          hhiggs[1][folder]->SetLineColor(hProj[nMZZ_ranges * 1 + sb][folder]->GetLineColor());
          hhiggs[1][folder]->Draw("histsame");
        }
        cout << "Passed draw" << endl;
      }
    }
    else if (sb == nMZZ_ranges - 1){
      cout << "Drawing Higgs MC now" << endl;
      for (int hh = 0; hh < kGGHSamples; hh++){
        if (hh==2)continue;
        if (hhiggs[hh][folder] == 0){
          cout << "Skipping draw for " << hh << "\t" << folder << endl;
          continue;
        }
        if (hh==0) hhiggs[hh][folder]->SetLineColor(kBlack);
        if (hh==1) hhiggs[hh][folder]->SetLineColor(kBlue);
        if (hh==2) hhiggs[hh][folder]->SetLineColor(kGreen+2);
        if (hh==3) hhiggs[hh][folder]->SetLineColor(kRed);
        hhiggs[hh][folder]->SetLineStyle(1);
        if (hh==0) hhiggs[hh][folder]->Draw("hist");
        else  hhiggs[hh][folder]->Draw("histsame");
      }
      l2D->Draw("same");
    }

    pt->Draw();

    double strSRlow = systZZMass_range[sb][0];
    double strSRhigh = systZZMass_range[sb][1];
    int strIntSRlow = (int)systZZMass_range[sb][0];
    int strIntSRhigh = (int)systZZMass_range[sb][1];
    TString strSidebandRegion;
    if (strSRlow>strIntSRlow) strSidebandRegion = Form("%.1f < m_{4l}", systZZMass_range[sb][0]);
    else strSidebandRegion = Form("%.0f < m_{4l}", systZZMass_range[sb][0]);
    if (strSRhigh>strIntSRhigh) strSidebandRegion.Append(Form(" < %.1f GeV", systZZMass_range[sb][1]));
    else strSidebandRegion.Append(Form(" < %.0f GeV", systZZMass_range[sb][1]));

    TPaveText *pt10 = new TPaveText(0.61, 0.84, 0.90, 0.92, "brNDC");
    pt10->SetBorderSize(0);
    pt10->SetTextAlign(12);
    pt10->SetTextSize(0.03);
    pt10->SetFillStyle(0);
    pt10->SetTextFont(42);
    TText* text10 = pt10->AddText(0.01, 0.01, strSidebandRegion);
    pt10->Draw();

    TString strSIPcut = "";
    if (removeSIP == 1) strSIPcut = "No SIP";
    else if (removeSIP == 2) strSIPcut = "New cut";
    strSIPcut.Prepend(", ");
    strSIPcut.Prepend(folder_id.c_str());
    TPaveText* pt20;
    if (removeSIP<2) pt20 = new TPaveText(0.62, 0.80, 0.90, 0.84, "brNDC");
    if (removeSIP==2) pt20 = new TPaveText(0.61, 0.80, 0.90, 0.84, "brNDC");
    pt20->SetBorderSize(0);
    pt20->SetTextAlign(12);
    pt20->SetTextSize(0.03);
    pt20->SetFillStyle(0);
    pt20->SetTextFont(42);
    TText* text20 = pt20->AddText(0.01, 0.01, strSIPcut);
    if (removeSIP>0) pt20->Draw();


    c2D->RedrawAxis();
    c2D->Modified();
    c2D->Update();
    foutput->WriteTObject(c2D);

    TString canvasDir = coutput_common;
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

    delete l2D;
    c2D->Close();
    delete pt;

    if (normOption == 1){
      hProj[nMZZ_ranges * 1 + sb][folder]->SetName(Form("%s_%s", hProj[nMZZ_ranges * 1 + sb][folder]->GetName(), user_folder[folder]));
      tgdata[nMZZ_ranges * folder + sb]->SetName(Form("%s_%s", tgdata[nMZZ_ranges * folder + sb]->GetName(), user_folder[folder]));
      foutput->WriteTObject(hProj[nMZZ_ranges * 1 + sb][folder]);
      foutput->WriteTObject(tgdata[nMZZ_ranges * folder + sb]);
    }
    for (int hh = 0; hh < ntrees; hh++) delete hProj[nMZZ_ranges * hh + sb][folder];
    delete tgdata[nMZZ_ranges * folder + sb];

  }

  for (int e = 0; e < 2; e++){
    if (finput[folder][e] == 0 || finput[folder][e]->IsZombie()) continue;
    else finput[folder][e]->Close();
  }
  foutput->Close();
}

void compare_KDShapeVariation_SIPVariation(int process = 3, int iSyst = 0, int normOption = 0){
	const int ntrees = 7;
	const int nSIPCuts = 3;
	TString chprocess[ntrees] = { "hqqzz_full_py_SignalRegion", "hggzz_full_py_SignalRegion", "hzx_full_py_SignalRegion", "hhiggs_CTau0_full", "hhiggs_CTau100_full", "hhiggs_CTau500_full", "hhiggs_CTau1000_full" };
	char* cnameprocess[ntrees] = { "qqzz", "ggzz", "zx", "CTau0", "CTau100", "CTau500", "CTau1000" };
	TString ctitleprocess[ntrees] = { "q#bar{q}#rightarrow4l bkg.", "gg#rightarrow4l bkg.", "Z+X bkg.", "Higgs c#tau_{H}=0 #mum", "Higgs c#tau_{H}=100 #mum", "Higgs c#tau_{H}=500 #mum", "Higgs c#tau_{H}=1000 #mum" };
  char* ctitle_SIPCut[nSIPCuts] = { "SIP_{3D}<4", "No cut", "New cut" };

	TString OUTPUT_NAME = "LifetimeKD_ShapeSystVars_SIPEffects_";
	if (normOption == 1) OUTPUT_NAME.Append("AreaNormalized_");
	OUTPUT_NAME.Append(Form("%s_%s_Comparison.root", cSystVariable[iSyst],cnameprocess[process]));
	TString coutput_common = user_dir_hep + "Analysis/Plots/";
	coutput_common.Append(Form("%s/", cSystVariable[iSyst]));
	TString coutput = coutput_common + OUTPUT_NAME;
	TFile* foutput = new TFile(coutput, "recreate");
	TH1F* hprocess[nSIPCuts] = { 0 };
	TFile* finput[nSIPCuts*3][2] = { { 0 } };
  bool received_ergtev[2] = { 0 };

	double max_plot = 0;
  for (int sipcut = 0; sipcut < nSIPCuts; sipcut++){
		for (int folder = 0; folder < 3; folder++){
			for (int erg_tev = 7; erg_tev < 9; erg_tev++){
				TString erg_dir;
				erg_dir.Form("LHC_%iTeV", erg_tev);
				TString comstring;
				comstring.Form("%iTeV", erg_tev);
				int EnergyIndex = 0; if (erg_tev == 8) EnergyIndex = 1;
				TString cinput_common = user_dir_hep + "Analysis/Plots/";
				cinput_common.Append(Form("%s/", cSystVariable[iSyst]));

				TString INPUT_NAME = "LifetimeKD_SystVariables_DataQQBZZ_";
				if (sipcut == 1) INPUT_NAME.Append("noSIP_");
				if (sipcut == 2) INPUT_NAME.Append("newSIP_");
				INPUT_NAME.Append(user_folder[folder]);
				INPUT_NAME.Append(Form("_%s_", cSystVariable[iSyst]));
				INPUT_NAME.Append(comstring);
				INPUT_NAME.Append(".root");
				cout << INPUT_NAME << endl;
				TString cinput = cinput_common + INPUT_NAME;
				cout << "Attempting to open " << cinput << endl;
				finput[3*sipcut+folder][EnergyIndex] = new TFile(cinput, "read");

				if (finput[3*sipcut+folder][EnergyIndex] == 0 || finput[3*sipcut+folder][EnergyIndex]->IsZombie()) continue;
        else received_ergtev[EnergyIndex]=true;
				TString hname = chprocess[process];
				cout << hname << endl;
				TH1F* htemp = (TH1F*)finput[3*sipcut+folder][EnergyIndex]->Get(hname);
				htemp->Scale(luminosity[EnergyIndex]);
				if (hprocess[sipcut] == 0){
					if(sipcut==0) hprocess[sipcut] = (TH1F*)htemp->Clone(Form("acc_%s_sip4", htemp->GetName()));
					else if(sipcut==1) hprocess[sipcut] = (TH1F*)htemp->Clone(Form("acc_%s_nosip", htemp->GetName()));
					else if(sipcut==2) hprocess[sipcut] = (TH1F*)htemp->Clone(Form("acc_%s_newcut", htemp->GetName()));
				}
				else hprocess[sipcut]->Add(htemp);
				delete htemp;
			}
		}
		cout << "Integral: " << hprocess[sipcut]->Integral(0, hprocess[sipcut]->GetNbinsX() + 1) << endl;
		if (normOption == 1){
			cout << "Scaling MC by 1./" << hprocess[sipcut]->Integral(0, hprocess[sipcut]->GetNbinsX() + 1) << endl;
			hprocess[sipcut]->Scale(1. / hprocess[sipcut]->Integral(0, hprocess[sipcut]->GetNbinsX() + 1) );
		}
		double maxY = hprocess[sipcut]->GetMaximum();
		if(maxY>max_plot) max_plot=maxY;
		hprocess[sipcut]->SetTitle("");
	}

	gStyle->SetTitleFont(62, "t");
	gROOT->SetStyle(gStyle->GetName());
	gROOT->ForceStyle();

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
  if (received_ergtev[0] && received_ergtev[1]) text = pt->AddText(0.537, 0.45, "#font[42]{19.7 fb^{-1} (8 TeV) + 5.1 fb^{-1} (7 TeV)}");
  else if (received_ergtev[1]) text = pt->AddText(0.837, 0.45, "#font[42]{19.7 fb^{-1} (8 TeV)}");
  else if (received_ergtev[0]) text = pt->AddText(0.837, 0.45, "#font[42]{ 5.1 fb^{-1} (7 TeV)}");
  text->SetTextSize(0.0315);

  hprocess[0]->SetLineColor(kBlue);
  hprocess[0]->SetMarkerColor(kBlue);
  hprocess[0]->SetLineWidth(2);
	hprocess[0]->SetLineStyle(1);
	hprocess[1]->SetLineColor(kBlack);
  hprocess[1]->SetMarkerColor(kBlack);
  hprocess[1]->SetLineWidth(2);
	hprocess[1]->SetLineStyle(1);
	hprocess[2]->SetLineColor(kRed);
  hprocess[2]->SetMarkerColor(kRed);
  hprocess[2]->SetLineWidth(2);
	hprocess[2]->SetLineStyle(1);

	foutput->cd();
	TString canvasname_2D = "cCompare_DataMC_AllChannels_VtxSystVar_";
	canvasname_2D.Append(Form("%s_SignalRegion_SIPEffect_", cSystVariable[iSyst]));
	canvasname_2D.Append(cnameprocess[process]);
	if(normOption==1) canvasname_2D.Append("_AreaNormalized");

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

	TLegend *l2D = new TLegend(0.20, 0.70, 0.58, 0.90);
	l2D->SetBorderSize(0);
	l2D->SetTextFont(42);
	l2D->SetTextSize(0.03);
	l2D->SetLineColor(1);
	l2D->SetLineStyle(1);
	l2D->SetLineWidth(1);
	l2D->SetFillColor(0);
	l2D->SetFillStyle(0);

	TString yTitle = (normOption == 0 ? Form("Events / %.3f", hprocess[0]->GetXaxis()->GetBinWidth(1)) : Form("Rate / %.3f", hprocess[0]->GetXaxis()->GetBinWidth(1)));
	if (iSyst == 3) yTitle = (normOption == 0 ? Form("Events / %.2e", hprocess[0]->GetXaxis()->GetBinWidth(1)) : Form("Rate / %.2e", hprocess[0]->GetXaxis()->GetBinWidth(1)));
	TString xTitle = cSystVariable_label[iSyst];

	hprocess[0]->GetXaxis()->SetTitle(xTitle);
	hprocess[0]->GetYaxis()->SetTitle(yTitle);
	hprocess[0]->GetYaxis()->SetRangeUser(0, max_plot * 1.5);
	if (iSyst == 3) hprocess[0]->GetYaxis()->SetRangeUser(0, max_plot * 1.85);
	if (iSyst == 2) hprocess[0]->GetYaxis()->SetRangeUser(0, max_plot * 1.75);
	hprocess[0]->GetXaxis()->SetNdivisions(505);
	hprocess[0]->GetXaxis()->SetLabelFont(42);
	hprocess[0]->GetXaxis()->SetLabelOffset(0.007);
	hprocess[0]->GetXaxis()->SetLabelSize(0.04);
	hprocess[0]->GetXaxis()->SetTitleSize(0.06);
	hprocess[0]->GetXaxis()->SetTitleOffset(0.9);
	hprocess[0]->GetXaxis()->SetTitleFont(42);
	hprocess[0]->GetYaxis()->SetNdivisions(505);
	hprocess[0]->GetYaxis()->SetLabelFont(42);
	hprocess[0]->GetYaxis()->SetLabelOffset(0.007);
	hprocess[0]->GetYaxis()->SetLabelSize(0.04);
	hprocess[0]->GetYaxis()->SetTitleSize(0.06);
	hprocess[0]->GetYaxis()->SetTitleOffset(1.1);
	hprocess[0]->GetYaxis()->SetTitleFont(42);
	hprocess[0]->GetXaxis()->SetRangeUser(visualRange_SystVariable[iSyst][0], visualRange_SystVariable[iSyst][1]);

  TH1F* htempfill;
	l2D->AddEntry(hprocess[1], Form("%s", ctitle_SIPCut[1]), "l");
	l2D->AddEntry(hprocess[0], Form("%s", ctitle_SIPCut[0]), "l");
	l2D->AddEntry(hprocess[2], Form("%s", ctitle_SIPCut[2]), "l");
	for (int sipcut = 0; sipcut < nSIPCuts; sipcut++){
		if(sipcut==0) hprocess[sipcut]->Draw("hist");
		else hprocess[sipcut]->Draw("histsame");
    if (sipcut==2){
      htempfill = (TH1F*) hprocess[sipcut]->Clone("htempfillcloned");
      htempfill->SetFillColor(kRed);
      htempfill->SetFillStyle(3002);
      htempfill->Draw("e2same");
    }
	}
	l2D->Draw("same");
	pt->Draw();

	TPaveText *pt10 = new TPaveText(0.61, 0.84, 0.90, 0.92, "brNDC");
	pt10->SetBorderSize(0);
	pt10->SetTextAlign(12);
	pt10->SetTextSize(0.03);
	pt10->SetFillStyle(0);
	pt10->SetTextFont(42);
	TText* text10 = pt10->AddText(0.01, 0.01, ctitleprocess[process]);
	pt10->Draw();

	c2D->RedrawAxis();
	c2D->Modified();
	c2D->Update();
	foutput->WriteTObject(c2D);

	TString canvasDir = coutput_common;
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

  delete htempfill;
	delete l2D;
	c2D->Close();
	delete pt;

	for (int sipcut = 0; sipcut < nSIPCuts; sipcut++){
		delete hprocess[sipcut];
		for (int folder = 0; folder < 3; folder++){
			for (int e = 0; e < 2; e++){
				if (finput[3 * sipcut + folder][e] == 0 || finput[3 * sipcut + folder][e]->IsZombie()) continue;
				else finput[3 * sipcut + folder][e]->Close();
			}
		}
	}
	foutput->Close();
}

void compare_HiggsProductions(int iSyst=0, int forANPAS=0){
  char TREE_NAME[]="SelectedTree";
  const int ntrees = 7; // ggH, VBF, WH, ZH, ttH, ggH minlo, ggH powheg+jhugen125
  TString chprocess[ntrees] ={ "hhiggs_CTau0_full", "hhiggs_CTau0_full_vbf", "hhiggs_CTau0_full_wh", "hhiggs_CTau0_full_zh", "hhiggs_CTau0_full_tth", "hhiggs_CTau0_full_ggh_minlo", "hhiggs_CTau0_full_ggh_powheg" };
  char* cnameprocess[ntrees] ={ "powheg15jhuGenV3-0PMH125.6", "VBFH125", "WH125", "ZH125", "ttH125", "minloH125", "powheg15jhuGenV3H125" };
  TString ctitleprocess[ntrees] ={ "ggH", "VBF H", "WH", "ZH", "t#bar{t}H", "ggH", "ggH" };
  TString ctitlemassprocess[ntrees] ={ " (125.6 GeV)", " (125 GeV)", " (125 GeV)", " (125 GeV)", " (125 GeV)", " (125 GeV, MINLO)", " (125 GeV)" };
  char ctitle_SIPCut[] = "New cut";

  TString strBeamSpot[2]={ "MC_START44_V13", "MC_START53_V23" };

  TString OUTPUT_NAME = "LifetimeKD_ShapeSystVars_HiggsProductions_";
  OUTPUT_NAME.Append(Form("%s_%s_Comparison.root", cSystVariable[iSyst], ctitle_SIPCut));
  TString coutput_common = user_dir_hep + "Analysis/Plots/HiggsProductions/";
  coutput_common.Append(Form("%s/", cSystVariable[iSyst]));
  TString mkdirCommand = "mkdir -p ";
  mkdirCommand.Append(coutput_common);
  gSystem->Exec(mkdirCommand);
  TString coutput = coutput_common + OUTPUT_NAME;
  TFile* foutput = new TFile(coutput, "recreate");
  TH1F* hprocess[ntrees] ={ 0 };
  TH1F* hprocess_rewgt[4] ={ 0 }; // VBF - ttH from ggH
  bool received_ergtev[2] ={ 0 };

  float MC_weight;
  float MC_weight_noxsec;
  float MC_weight_QQZZEWK=1;
  float MC_weight_Kfactor=1;
  float MC_weight_QQBGGProper[4];
  float ZZMass, ZZPt, ZZEta, ZZPhi;
  float KalmanCandVtx_x, KalmanCandVtx_y, KalmanCandVtx_z, KalmanCandVtx_chi2;
  float OfflinePrimaryVtx_x, OfflinePrimaryVtx_y, OfflinePrimaryVtx_z, OfflinePrimaryVtx_ndof;
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

  float Lep1combRelIsoPF=0, Lep2combRelIsoPF=0, Lep3combRelIsoPF=0, Lep4combRelIsoPF=0;
  bool Lep3isID=true, Lep4isID=true;
  short Z1ids;

  double max_plot = 0;
  TChain* thiggs;
  for (int erg_tev = 7; erg_tev < 9; erg_tev++){
    TString erg_dir;
    erg_dir.Form("LHC_%iTeV", erg_tev);
    TString comstring;
    comstring.Form("%iTeV", erg_tev);
    int EnergyIndex = 0; if (erg_tev == 8) EnergyIndex = 1;

    // BeamSpot info
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

    TFile* finput_pTrewgt = new TFile(Form("./data/compareProductionModes_pT_%iTeV%s",erg_tev,".root"), "read");
    TF1* fit_pT[4];
    for (int prod=0; prod<4; prod++){
      TString fitname = cnameprocess[prod+1];
      fitname.Append("_pToverMZZ_ratio_fit");
      fit_pT[prod] = (TF1*)finput_pTrewgt->Get(fitname);
    }

    for (int folder = 0; folder < 3; folder++){
      for (int prod = 0; prod < ntrees; prod++){
        foutput->cd();
        TString cinput_common;
        if (prod==0){
          cinput_common = user_dir_hep + "Analysis/Plots/";
          cinput_common.Append(Form("%s/", cSystVariable[iSyst]));

          TString INPUT_NAME = "LifetimeKD_SystVariables_DataQQBZZ_";
          INPUT_NAME.Append("newSIP_");
          INPUT_NAME.Append(user_folder[folder]);
          INPUT_NAME.Append(Form("_%s_", cSystVariable[iSyst]));
          INPUT_NAME.Append(comstring);
          INPUT_NAME.Append(".root");
          cout << INPUT_NAME << endl;
          TString cinput = cinput_common + INPUT_NAME;
          cout << "Attempting to open " << cinput << endl;
          TFile* finput = new TFile(cinput, "read");

          if (finput == 0 || finput->IsZombie()){
            if (finput->IsZombie()) delete finput;
            continue;
          }
          else received_ergtev[EnergyIndex]=true;
          TString hname = chprocess[prod];
          cout << hname << endl;
          TH1F* htemp = (TH1F*)finput->Get(hname);
          htemp->Scale(luminosity[EnergyIndex]);
          if (hprocess[prod] == 0){
            foutput->cd();
            hprocess[prod] = (TH1F*)htemp->Clone(Form("acc_%s_newcut", htemp->GetName()));
          }
          else hprocess[prod]->Add(htemp);
          delete htemp;
          finput->Close();
        }
        else{
          cout << endl;
          cout << "Starting the acquisition of " << cnameprocess[prod] << endl;
          cinput_common = user_dir_hep + "No_SIP/" + erg_dir + "/" + user_folder[folder] + "/";
          TString cinput = cinput_common + "HZZ4lTree_" + cnameprocess[prod] + ".root";

          thiggs = new TChain(TREE_NAME);
          thiggs->Add(cinput);

          cout << "nEntries: " << thiggs->GetEntries() << endl;

          thiggs->SetBranchAddress("MC_weight", &MC_weight);
          thiggs->SetBranchAddress("MC_weight_noxsec", &MC_weight_noxsec);
          thiggs->SetBranchAddress("ZZMass", &ZZMass);
          thiggs->SetBranchAddress("ZZPt", &ZZPt);
          thiggs->SetBranchAddress("ZZEta", &ZZEta);
          thiggs->SetBranchAddress("ZZPhi", &ZZPhi);
          thiggs->SetBranchAddress("OfflinePrimaryVtx_x", &OfflinePrimaryVtx_x);
          thiggs->SetBranchAddress("OfflinePrimaryVtx_y", &OfflinePrimaryVtx_y);
          thiggs->SetBranchAddress("OfflinePrimaryVtx_z", &OfflinePrimaryVtx_z);
          thiggs->SetBranchAddress("KalmanCandVtx_x", &KalmanCandVtx_x);
          thiggs->SetBranchAddress("KalmanCandVtx_y", &KalmanCandVtx_y);
          thiggs->SetBranchAddress("KalmanCandVtx_z", &KalmanCandVtx_z);
          thiggs->SetBranchAddress("OfflinePrimaryVtx_cov_xx", &OfflinePrimaryVtx_cov_xx);
          thiggs->SetBranchAddress("OfflinePrimaryVtx_cov_xy", &OfflinePrimaryVtx_cov_xy);
          thiggs->SetBranchAddress("OfflinePrimaryVtx_cov_xz", &OfflinePrimaryVtx_cov_xz);
          thiggs->SetBranchAddress("OfflinePrimaryVtx_cov_yy", &OfflinePrimaryVtx_cov_yy);
          thiggs->SetBranchAddress("OfflinePrimaryVtx_cov_yz", &OfflinePrimaryVtx_cov_yz);
          thiggs->SetBranchAddress("OfflinePrimaryVtx_cov_zz", &OfflinePrimaryVtx_cov_zz);
          thiggs->SetBranchAddress("KalmanCandVtx_cov_xx", &KalmanCandVtx_cov_xx);
          thiggs->SetBranchAddress("KalmanCandVtx_cov_xy", &KalmanCandVtx_cov_xy);
          thiggs->SetBranchAddress("KalmanCandVtx_cov_xz", &KalmanCandVtx_cov_xz);
          thiggs->SetBranchAddress("KalmanCandVtx_cov_yy", &KalmanCandVtx_cov_yy);
          thiggs->SetBranchAddress("KalmanCandVtx_cov_yz", &KalmanCandVtx_cov_yz);
          thiggs->SetBranchAddress("KalmanCandVtx_cov_zz", &KalmanCandVtx_cov_zz);

          thiggs->SetBranchAddress("OfflinePrimaryVtx_ndof", &OfflinePrimaryVtx_ndof);


          if (!thiggs->GetBranchStatus("GenPrimaryVtx_x")) cout << "Tree has no gen. vtx" << endl;
          else{
            thiggs->SetBranchAddress("GenPrimaryVtx_x", &GenPrimaryVtx_x);
            thiggs->SetBranchAddress("GenPrimaryVtx_y", &GenPrimaryVtx_y);
            thiggs->SetBranchAddress("GenPrimaryVtx_z", &GenPrimaryVtx_z);
            thiggs->SetBranchAddress("GenIntVtx_x", &GenIntVtx_x);
            thiggs->SetBranchAddress("GenIntVtx_y", &GenIntVtx_y);
            thiggs->SetBranchAddress("GenIntVtx_z", &GenIntVtx_z);
            thiggs->SetBranchAddress("GenHPt", &GenHPt);
            thiggs->SetBranchAddress("GenHPhi", &GenHPhi);
            thiggs->SetBranchAddress("GenHMass", &GenHMass);
          }
          thiggs->SetBranchAddress("Lep1_Z1SIP", &Lep1_Z1SIP);
          thiggs->SetBranchAddress("Lep2_Z1SIP", &Lep2_Z1SIP);
          thiggs->SetBranchAddress("Lep3_Z1SIP", &Lep3_Z1SIP);
          thiggs->SetBranchAddress("Lep4_Z1SIP", &Lep4_Z1SIP);
          thiggs->SetBranchAddress("KalmanCandVtx_chi2", &KalmanCandVtx_chi2);
          thiggs->SetBranchAddress("Lep1SIP", &Lep1SIP);
          thiggs->SetBranchAddress("Lep2SIP", &Lep2SIP);
          thiggs->SetBranchAddress("Lep3SIP", &Lep3SIP);
          thiggs->SetBranchAddress("Lep4SIP", &Lep4SIP);
          if (thiggs->GetBranchStatus("D_bkg")){
            thiggs->SetBranchAddress("D_bkg", &D_bkg);
          }
          else if (thiggs->GetBranchStatus("p0plus_VAJHU")){
            thiggs->SetBranchAddress("p0plus_VAJHU", &p0plus_VAJHU);
            thiggs->SetBranchAddress("bkg_VAMCFM", &bkg_VAMCFM);
            thiggs->SetBranchAddress("p0plus_m4l", &p0plus_m4l);
            thiggs->SetBranchAddress("bkg_m4l", &bkg_m4l);
          }
          else cout << "Could not find p0plus_VAJHU." << endl;

          TH1F* hhiggs = (TH1F*)hprocess[0]->Clone(chprocess[prod]);
          hhiggs->Reset("ICESM");
          hhiggs->Sumw2();

          double sum_sipcut[3] ={ 0 };
          for (int ev = 0; ev < thiggs->GetEntries(); ev++){
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

            thiggs->GetEntry(ev);
            if (!(ZZMass >= 105 && ZZMass < 140)) continue;
            float wgt = MC_weight*luminosity[EnergyIndex];

            sum_sipcut[0] += wgt;
            if (!(fabs(Lep1SIP) >= 4 || fabs(Lep2SIP) >= 4 || fabs(Lep3SIP) >= 4 || fabs(Lep4SIP) >= 4)){
              sum_sipcut[1] += wgt;
            }
            if (fabs(Lep1_Z1SIP) >= 4 || fabs(Lep2_Z1SIP) >= 4 || fabs(Lep3_Z1SIP) >= 5 || fabs(Lep4_Z1SIP) >= 5 || KalmanCandVtx_chi2 >= 30) continue;
            sum_sipcut[2] += wgt;

            if (iSyst>=22 && iSyst<nSystVars-1){ // BE CAREFUL HERE, HARDCODED BOUNDS AHEAD
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

            if (!thiggs->GetBranchStatus("D_bkg")) D_bkg = p0plus_VAJHU*p0plus_m4l / (p0plus_VAJHU*p0plus_m4l + bkg_VAMCFM*bkg_m4l);

            if (iSyst==0 || iSyst==20) KD = Txy;
            if (iSyst==1) KD = Dxy;
            if (iSyst==2) KD = delTxy;
            if (iSyst==3) KD = delDxy;
            if (iSyst==12) KD = (OfflinePrimaryVtx_x - GenPrimaryVtx_x)*cos(GenHPhi) + (OfflinePrimaryVtx_y - GenPrimaryVtx_y)*sin(GenHPhi);
            if (iSyst==13) KD = (KalmanCandVtx_x - GenIntVtx_x)*cos(GenHPhi) + (KalmanCandVtx_y - GenIntVtx_y)*sin(GenHPhi);
            if (iSyst==14) KD = (KalmanCandVtx_x-GenIntVtx_x);
            if (iSyst==15) KD = (KalmanCandVtx_y-GenIntVtx_y);
            if (iSyst==16) KD = (KalmanCandVtx_z-GenIntVtx_z);
            if (iSyst==17) KD = (OfflinePrimaryVtx_x - GenPrimaryVtx_x);
            if (iSyst==18) KD = (OfflinePrimaryVtx_y - GenPrimaryVtx_y);
            if (iSyst==19) KD = (OfflinePrimaryVtx_z - GenPrimaryVtx_z);

            if (iSyst==22) KD = (OfflinePrimaryVtx_x - BeamPosX);
            if (iSyst==23) KD = (OfflinePrimaryVtx_y - BeamPosY);
            if (iSyst==24) KD = (OfflinePrimaryVtx_z - BeamPosZ)/1000.;
            if (iSyst==25) KD = (GenPrimaryVtx_x - BeamPosX);
            if (iSyst==26) KD = (GenPrimaryVtx_y - BeamPosY);
            if (iSyst==27) KD = (GenPrimaryVtx_z - BeamPosZ)/1000.;
            if (iSyst==28) KD = Txy_BS;
            if (iSyst==29) KD = Txy_BS - Txy;
            if (iSyst==30) KD = Txy_BS - Txy_true;

            KD = KD*10000.; // in 1 um

            if (iSyst==20){
              KD = tanh(KD/800.); if (KD==1.) KD=0.999999; if (KD==-1.) KD=-0.999999;
            }
            if (iSyst==21){
              KD = D_bkg; if (KD==1.) KD=0.999999;;
            }

            if (iSyst==4) KD = ZZPt/ZZMass;
            if (iSyst==5) KD = (Txy-Txy_true)/delTxy;
            if (iSyst==6) KD = (Dxy-Dxy_true)/delDxy;
            if (iSyst==7) KD = ((OfflinePrimaryVtx_x - GenPrimaryVtx_x)*cos(GenHPhi) + (OfflinePrimaryVtx_y - GenPrimaryVtx_y)*sin(GenHPhi))/sigmaPV_xy;
            if (iSyst==8) KD = ((KalmanCandVtx_x - GenIntVtx_x)*cos(GenHPhi) + (KalmanCandVtx_y - GenIntVtx_y)*sin(GenHPhi))/sigmaInt_xy;
            if (iSyst==9) KD = (KalmanCandVtx_x-GenIntVtx_x)/sqrt(KalmanCandVtx_cov_xx);
            if (iSyst==10) KD = (KalmanCandVtx_y-GenIntVtx_y)/sqrt(KalmanCandVtx_cov_yy);
            if (iSyst==11) KD = (KalmanCandVtx_z-GenIntVtx_z)/sqrt(KalmanCandVtx_cov_zz);
            if (iSyst==31) KD = (OfflinePrimaryVtx_ndof+3.)/2.;

            hhiggs->Fill(KD, wgt);
          }

          if (hprocess[prod] == 0){
            foutput->cd();
            hprocess[prod] = (TH1F*)hhiggs->Clone(Form("acc_%s_newcut", hhiggs->GetName()));
          }
          else hprocess[prod]->Add(hhiggs);

          cout << user_folder[folder] << " " << comstring << " " << hhiggs->GetName() << endl;
          cout << "New/All: " << sum_sipcut[2]/sum_sipcut[0] << " | Old/All: " << sum_sipcut[1]/sum_sipcut[0] <<  " | New/Old: " << sum_sipcut[2]/sum_sipcut[1] << endl;

          delete hhiggs;
          delete thiggs;

          cout << "Finished " << cnameprocess[prod] << " counting." << endl;

          // Start reweighting from ggH
          if (prod<ntrees-2){
            foutput->cd();

            cout << endl;
            cout << "Reweighting ggH to " << cnameprocess[prod] << " at prod = " << prod << endl;
            cinput = cinput_common + "HZZ4lTree_" + cnameprocess[0] + ".root";

            thiggs = new TChain(TREE_NAME);
            thiggs->Add(cinput);
            cout << "Acquired " << cinput << endl;

            thiggs->SetBranchAddress("MC_weight", &MC_weight);
            thiggs->SetBranchAddress("MC_weight_noxsec", &MC_weight_noxsec);
            thiggs->SetBranchAddress("ZZMass", &ZZMass);
            thiggs->SetBranchAddress("ZZPt", &ZZPt);
            thiggs->SetBranchAddress("ZZEta", &ZZEta);
            thiggs->SetBranchAddress("ZZPhi", &ZZPhi);
            thiggs->SetBranchAddress("OfflinePrimaryVtx_x", &OfflinePrimaryVtx_x);
            thiggs->SetBranchAddress("OfflinePrimaryVtx_y", &OfflinePrimaryVtx_y);
            thiggs->SetBranchAddress("OfflinePrimaryVtx_z", &OfflinePrimaryVtx_z);
            thiggs->SetBranchAddress("KalmanCandVtx_x", &KalmanCandVtx_x);
            thiggs->SetBranchAddress("KalmanCandVtx_y", &KalmanCandVtx_y);
            thiggs->SetBranchAddress("KalmanCandVtx_z", &KalmanCandVtx_z);
            thiggs->SetBranchAddress("OfflinePrimaryVtx_cov_xx", &OfflinePrimaryVtx_cov_xx);
            thiggs->SetBranchAddress("OfflinePrimaryVtx_cov_xy", &OfflinePrimaryVtx_cov_xy);
            thiggs->SetBranchAddress("OfflinePrimaryVtx_cov_xz", &OfflinePrimaryVtx_cov_xz);
            thiggs->SetBranchAddress("OfflinePrimaryVtx_cov_yy", &OfflinePrimaryVtx_cov_yy);
            thiggs->SetBranchAddress("OfflinePrimaryVtx_cov_yz", &OfflinePrimaryVtx_cov_yz);
            thiggs->SetBranchAddress("OfflinePrimaryVtx_cov_zz", &OfflinePrimaryVtx_cov_zz);
            thiggs->SetBranchAddress("KalmanCandVtx_cov_xx", &KalmanCandVtx_cov_xx);
            thiggs->SetBranchAddress("KalmanCandVtx_cov_xy", &KalmanCandVtx_cov_xy);
            thiggs->SetBranchAddress("KalmanCandVtx_cov_xz", &KalmanCandVtx_cov_xz);
            thiggs->SetBranchAddress("KalmanCandVtx_cov_yy", &KalmanCandVtx_cov_yy);
            thiggs->SetBranchAddress("KalmanCandVtx_cov_yz", &KalmanCandVtx_cov_yz);
            thiggs->SetBranchAddress("KalmanCandVtx_cov_zz", &KalmanCandVtx_cov_zz);

            thiggs->SetBranchAddress("OfflinePrimaryVtx_ndof", &OfflinePrimaryVtx_ndof);


            if (!thiggs->GetBranchStatus("GenPrimaryVtx_x")) cout << "Tree has no gen. vtx" << endl;
            else{
              thiggs->SetBranchAddress("GenPrimaryVtx_x", &GenPrimaryVtx_x);
              thiggs->SetBranchAddress("GenPrimaryVtx_y", &GenPrimaryVtx_y);
              thiggs->SetBranchAddress("GenPrimaryVtx_z", &GenPrimaryVtx_z);
              thiggs->SetBranchAddress("GenIntVtx_x", &GenIntVtx_x);
              thiggs->SetBranchAddress("GenIntVtx_y", &GenIntVtx_y);
              thiggs->SetBranchAddress("GenIntVtx_z", &GenIntVtx_z);
              thiggs->SetBranchAddress("GenHPt", &GenHPt);
              thiggs->SetBranchAddress("GenHPhi", &GenHPhi);
              thiggs->SetBranchAddress("GenHMass", &GenHMass);
            }
            thiggs->SetBranchAddress("Lep1_Z1SIP", &Lep1_Z1SIP);
            thiggs->SetBranchAddress("Lep2_Z1SIP", &Lep2_Z1SIP);
            thiggs->SetBranchAddress("Lep3_Z1SIP", &Lep3_Z1SIP);
            thiggs->SetBranchAddress("Lep4_Z1SIP", &Lep4_Z1SIP);
            thiggs->SetBranchAddress("KalmanCandVtx_chi2", &KalmanCandVtx_chi2);
            thiggs->SetBranchAddress("Lep1SIP", &Lep1SIP);
            thiggs->SetBranchAddress("Lep2SIP", &Lep2SIP);
            thiggs->SetBranchAddress("Lep3SIP", &Lep3SIP);
            thiggs->SetBranchAddress("Lep4SIP", &Lep4SIP);
            if (thiggs->GetBranchStatus("D_bkg")){
              thiggs->SetBranchAddress("D_bkg", &D_bkg);
            }
            else if (thiggs->GetBranchStatus("p0plus_VAJHU")){
              thiggs->SetBranchAddress("p0plus_VAJHU", &p0plus_VAJHU);
              thiggs->SetBranchAddress("bkg_VAMCFM", &bkg_VAMCFM);
              thiggs->SetBranchAddress("p0plus_m4l", &p0plus_m4l);
              thiggs->SetBranchAddress("bkg_m4l", &bkg_m4l);
            }
            else cout << "Could not find p0plus_VAJHU." << endl;

            hhiggs = (TH1F*)hprocess[0]->Clone(Form("%s_rewgt",chprocess[prod].Data()));
            hhiggs->Reset("ICESM");
            hhiggs->Sumw2();

            double sum_sipcut_rewgt[3] ={ 0 };
            cout << "nEntries:  " << thiggs->GetEntries() << endl;
            for (int ev = 0; ev < thiggs->GetEntries(); ev++){
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

              thiggs->GetEntry(ev);
              if (!(ZZMass >= 105 && ZZMass < 140)) continue;
              float wgt = MC_weight*luminosity[EnergyIndex];
              if (fit_pT[prod-1]!=0) wgt *= fit_pT[prod-1]->Eval(GenHPt/GenHMass);
              else cerr << "Warning! Fit could not be evaluated." << endl;
              if (fit_pT[prod-1]!=0 && ev == 0) cout << "Address of fit is " << fit_pT[prod-1] << endl;

              sum_sipcut_rewgt[0] += wgt;
              if (!(fabs(Lep1SIP) >= 4 || fabs(Lep2SIP) >= 4 || fabs(Lep3SIP) >= 4 || fabs(Lep4SIP) >= 4)){
                sum_sipcut_rewgt[1] += wgt;
              }
              if (fabs(Lep1_Z1SIP) >= 4 || fabs(Lep2_Z1SIP) >= 4 || fabs(Lep3_Z1SIP) >= 5 || fabs(Lep4_Z1SIP) >= 5 || KalmanCandVtx_chi2 >= 30) continue;
              sum_sipcut_rewgt[2] += wgt;

              if (iSyst>=22 && iSyst<nSystVars-1){ // BE CAREFUL HERE, HARDCODED BOUNDS AHEAD
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

              if (!thiggs->GetBranchStatus("D_bkg")) D_bkg = p0plus_VAJHU*p0plus_m4l / (p0plus_VAJHU*p0plus_m4l + bkg_VAMCFM*bkg_m4l);

              if (iSyst==0 || iSyst==20) KD = Txy;
              if (iSyst==1) KD = Dxy;
              if (iSyst==2) KD = delTxy;
              if (iSyst==3) KD = delDxy;
              if (iSyst==12) KD = (OfflinePrimaryVtx_x - GenPrimaryVtx_x)*cos(GenHPhi) + (OfflinePrimaryVtx_y - GenPrimaryVtx_y)*sin(GenHPhi);
              if (iSyst==13) KD = (KalmanCandVtx_x - GenIntVtx_x)*cos(GenHPhi) + (KalmanCandVtx_y - GenIntVtx_y)*sin(GenHPhi);
              if (iSyst==14) KD = (KalmanCandVtx_x-GenIntVtx_x);
              if (iSyst==15) KD = (KalmanCandVtx_y-GenIntVtx_y);
              if (iSyst==16) KD = (KalmanCandVtx_z-GenIntVtx_z);
              if (iSyst==17) KD = (OfflinePrimaryVtx_x - GenPrimaryVtx_x);
              if (iSyst==18) KD = (OfflinePrimaryVtx_y - GenPrimaryVtx_y);
              if (iSyst==19) KD = (OfflinePrimaryVtx_z - GenPrimaryVtx_z);

              if (iSyst==22) KD = (OfflinePrimaryVtx_x - BeamPosX);
              if (iSyst==23) KD = (OfflinePrimaryVtx_y - BeamPosY);
              if (iSyst==24) KD = (OfflinePrimaryVtx_z - BeamPosZ)/1000.;
              if (iSyst==25) KD = (GenPrimaryVtx_x - BeamPosX);
              if (iSyst==26) KD = (GenPrimaryVtx_y - BeamPosY);
              if (iSyst==27) KD = (GenPrimaryVtx_z - BeamPosZ)/1000.;
              if (iSyst==28) KD = Txy_BS;
              if (iSyst==29) KD = Txy_BS - Txy;
              if (iSyst==30) KD = Txy_BS - Txy_true;

              KD = KD*10000.; // in 1 um

              if (iSyst==20){
                KD = tanh(KD/800.); if (KD==1.) KD=0.999999; if (KD==-1.) KD=-0.999999;
              }
              if (iSyst==21){
                KD = D_bkg; if (KD==1.) KD=0.999999;;
              }

              if (iSyst==4) KD = ZZPt/ZZMass;
              if (iSyst==5) KD = (Txy-Txy_true)/delTxy;
              if (iSyst==6) KD = (Dxy-Dxy_true)/delDxy;
              if (iSyst==7) KD = ((OfflinePrimaryVtx_x - GenPrimaryVtx_x)*cos(GenHPhi) + (OfflinePrimaryVtx_y - GenPrimaryVtx_y)*sin(GenHPhi))/sigmaPV_xy;
              if (iSyst==8) KD = ((KalmanCandVtx_x - GenIntVtx_x)*cos(GenHPhi) + (KalmanCandVtx_y - GenIntVtx_y)*sin(GenHPhi))/sigmaInt_xy;
              if (iSyst==9) KD = (KalmanCandVtx_x-GenIntVtx_x)/sqrt(KalmanCandVtx_cov_xx);
              if (iSyst==10) KD = (KalmanCandVtx_y-GenIntVtx_y)/sqrt(KalmanCandVtx_cov_yy);
              if (iSyst==11) KD = (KalmanCandVtx_z-GenIntVtx_z)/sqrt(KalmanCandVtx_cov_zz);
              if (iSyst==31) KD = (OfflinePrimaryVtx_ndof+3.)/2.;

              hhiggs->Fill(KD, wgt);
            }
            cout << "Accumulated weights: " << hhiggs->Integral(0, hhiggs->GetNbinsX()) << endl;

            double scale_rewgt = sum_sipcut[2] / sum_sipcut_rewgt[2];
            cout << "Scale: " << scale_rewgt << endl;
            hhiggs->Scale(scale_rewgt);

            if (hprocess_rewgt[prod-1] == 0){
              foutput->cd();
              cout << "Writing new rewgt histo" << endl;
              hprocess_rewgt[prod-1] = (TH1F*)hhiggs->Clone(Form("acc_%s_newcut", hhiggs->GetName()));
              cout << "Wrote new rewgt histo" << endl;
            }
            else{
              cout << "Adding to rewgt histo" << endl;
              hprocess_rewgt[prod-1]->Add(hhiggs);
              cout << "Added to rewgt histo" << endl;
            }
            cout << "Cleaning up reweighting" << endl;

            delete hhiggs;
            cout << "Cleaning up reweighted tree" << endl;
            delete thiggs;
            cout << "Clean up complete" << endl;
          }
        }
        cout << "Closed prod " << prod << endl;

      }

    }

    cout << "Cleaning " << comstring << endl;
    for (int prod=0; prod<4; prod++) delete fit_pT[prod];
    finput_pTrewgt->Close();
    delete tBeam;
  }

  foutput->cd();
  max_plot = 0;
  for (int prod=0; prod<ntrees; prod++){
    if (forANPAS==0){
      if (iSyst<2 || iSyst>=12 &&iSyst<=20 || iSyst>=22 && iSyst<=28) symmetrize_PromptTemplates(hprocess[prod]);
    }
    hprocess[prod]->SetTitle("");
    hprocess[prod]->Scale(1./hprocess[prod]->Integral(0, hprocess[prod]->GetNbinsX()+1));
    max_plot = max(max_plot, hprocess[prod]->GetMaximum());
    foutput->WriteTObject(hprocess[prod]);
  }
  for (int prod=0; prod<4; prod++){
    hprocess_rewgt[prod]->SetTitle("");
    hprocess_rewgt[prod]->Scale(1./hprocess_rewgt[prod]->Integral(0, hprocess_rewgt[prod]->GetNbinsX()+1));
    max_plot = max(max_plot, hprocess_rewgt[prod]->GetMaximum());
    foutput->WriteTObject(hprocess_rewgt[prod]);
  }

  gStyle->SetTitleFont(62, "t");
  gROOT->SetStyle(gStyle->GetName());
  gROOT->ForceStyle();

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
  if (received_ergtev[0] && received_ergtev[1]) text = pt->AddText(0.537, 0.45, "#font[42]{19.7 fb^{-1} (8 TeV) + 5.1 fb^{-1} (7 TeV)}");
  else if (received_ergtev[1]) text = pt->AddText(0.837, 0.45, "#font[42]{19.7 fb^{-1} (8 TeV)}");
  else if (received_ergtev[0]) text = pt->AddText(0.837, 0.45, "#font[42]{ 5.1 fb^{-1} (7 TeV)}");
  text->SetTextSize(0.0315);

  hprocess[0]->SetLineColor(kBlack);
  hprocess[0]->SetMarkerColor(kBlack);
  hprocess[0]->SetLineWidth(2);
  hprocess[0]->SetLineStyle(1);
  hprocess[1]->SetLineColor(kViolet);
  hprocess[1]->SetMarkerColor(kViolet);
  hprocess[1]->SetLineWidth(2);
  hprocess[1]->SetLineStyle(1);
  hprocess[2]->SetLineColor(kBlue);
  hprocess[2]->SetMarkerColor(kBlue);
  hprocess[2]->SetLineWidth(2);
  hprocess[2]->SetLineStyle(1);
  hprocess[3]->SetLineColor(kRed);
  hprocess[3]->SetMarkerColor(kRed);
  hprocess[3]->SetLineWidth(2);
  hprocess[3]->SetLineStyle(1);
  hprocess[4]->SetLineColor(kGreen+2);
  hprocess[4]->SetMarkerColor(kGreen+2);
  hprocess[4]->SetLineWidth(2);
  hprocess[4]->SetLineStyle(1);
  hprocess[5]->SetLineColor(kBlack);
  hprocess[5]->SetMarkerColor(kBlack);
  hprocess[5]->SetLineWidth(2);
  hprocess[5]->SetLineStyle(7);
  hprocess[6]->SetLineColor(kBlack);
  hprocess[6]->SetMarkerColor(kBlack);
  hprocess[6]->SetLineWidth(2);
  if (forANPAS==1) hprocess[6]->SetLineStyle(3);
  else hprocess[6]->SetLineStyle(1);

  hprocess_rewgt[0]->SetLineColor(kViolet);
  hprocess_rewgt[0]->SetMarkerColor(kViolet);
  hprocess_rewgt[0]->SetLineWidth(2);
  hprocess_rewgt[0]->SetLineStyle(7);
  hprocess_rewgt[1]->SetLineColor(kBlue);
  hprocess_rewgt[1]->SetMarkerColor(kBlue);
  hprocess_rewgt[1]->SetLineWidth(2);
  hprocess_rewgt[1]->SetLineStyle(7);
  hprocess_rewgt[2]->SetLineColor(kRed);
  hprocess_rewgt[2]->SetMarkerColor(kRed);
  hprocess_rewgt[2]->SetLineWidth(2);
  hprocess_rewgt[2]->SetLineStyle(7);
  hprocess_rewgt[3]->SetLineColor(kGreen+2);
  hprocess_rewgt[3]->SetMarkerColor(kGreen+2);
  hprocess_rewgt[3]->SetLineWidth(2);
  hprocess_rewgt[3]->SetLineStyle(7);

  foutput->cd();
  TString canvasname_2D = "cCompare_HiggsProductions_VtxSystVar_";
  canvasname_2D.Append(Form("%s_SignalRegion", cSystVariable[iSyst]));
  if (forANPAS==0) canvasname_2D.Append("_Twiki");

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


  float lxmin = 0.20;
  float lxwidth = 0.38;
  float lymax = 0.9;
  float lywidth = 0.04;
  float lxmax = lxmin + lxwidth;
  float lymin = lymax;
  lymin -= lywidth*3.;

  float lxmin2 = lxmin+0.41;
  if (forANPAS==0) lxmin2 += 0.1;
  float lymax2 = lymax;
  float lxmax2 = lxmin2 + lxwidth;
  float lymin2 = lymax2;
  lymin2 -= lywidth*4.;
  if (forANPAS==0) lymin2 -= lywidth;

  TLegend *l2D = new TLegend(lxmin, lymin, lxmax, lymax);
  TLegend *l2D_2 = new TLegend(lxmin2, lymin2, lxmax2, lymax2);
  l2D->SetBorderSize(0);
  l2D->SetTextFont(42);
  l2D->SetTextSize(0.03);
  l2D->SetLineColor(1);
  l2D->SetLineStyle(1);
  l2D->SetLineWidth(1);
  l2D->SetFillColor(0);
  l2D->SetFillStyle(0);
  l2D_2->SetBorderSize(0);
  l2D_2->SetTextFont(42);
  l2D_2->SetTextSize(0.03);
  l2D_2->SetLineColor(1);
  l2D_2->SetLineStyle(1);
  l2D_2->SetLineWidth(1);
  l2D_2->SetFillColor(0);
  l2D_2->SetFillStyle(0);

  double xbinwidth = hprocess[0]->GetXaxis()->GetBinWidth(1);
  int xbinwidth_int = (int)xbinwidth;
  double xbinwidth_intdouble = (double)xbinwidth_int;
  TString yTitle = Form("Rate / %.2f", xbinwidth);
  if (xbinwidth_intdouble==xbinwidth) yTitle = Form("Rate / %.0f", xbinwidth);
  else if (xbinwidth>=0.1) yTitle = Form("Rate / %.1f", xbinwidth);
  if (iSyst==31) yTitle.Append(" tracks");
  else if (iSyst<4 || iSyst>=12 && iSyst<=19 || iSyst>=22 && iSyst<27 || iSyst>27 && iSyst<31) yTitle.Append(" #mum");
  else if (iSyst==27) yTitle.Append(" mm");
  TString xTitle = cSystVariable_label[iSyst];

  for (int prod=0; prod<ntrees; prod++){
    hprocess[prod]->GetXaxis()->SetTitle(xTitle);
    hprocess[prod]->GetYaxis()->SetTitle(yTitle);
    hprocess[prod]->GetYaxis()->SetRangeUser(0, max_plot * 1.65);
    if (forANPAS==0) hprocess[prod]->GetYaxis()->SetRangeUser(0, max_plot * 1.3);
    hprocess[prod]->GetXaxis()->SetNdivisions(505);
    hprocess[prod]->GetXaxis()->SetLabelFont(42);
    hprocess[prod]->GetXaxis()->SetLabelOffset(0.007);
    hprocess[prod]->GetXaxis()->SetLabelSize(0.04);
    hprocess[prod]->GetXaxis()->SetTitleSize(0.06);
    hprocess[prod]->GetXaxis()->SetTitleOffset(0.9);
    hprocess[prod]->GetXaxis()->SetTitleFont(42);
    hprocess[prod]->GetYaxis()->SetNdivisions(505);
    hprocess[prod]->GetYaxis()->SetLabelFont(42);
    hprocess[prod]->GetYaxis()->SetLabelOffset(0.007);
    hprocess[prod]->GetYaxis()->SetLabelSize(0.04);
    hprocess[prod]->GetYaxis()->SetTitleSize(0.06);
    hprocess[prod]->GetYaxis()->SetTitleOffset(1.1);
    if (max_plot<0.2) hprocess[prod]->GetYaxis()->SetTitleOffset(1.2);
    hprocess[prod]->GetYaxis()->SetTitleFont(42);
    hprocess[prod]->GetXaxis()->SetRangeUser(visualRange_SystVariable[iSyst][0], visualRange_SystVariable[iSyst][1]);
    if (iSyst==28) hprocess[prod]->GetXaxis()->SetRangeUser(-800, 800);
    if (iSyst==12 || iSyst==13) hprocess[prod]->GetXaxis()->SetRangeUser(-105, 105);
  }

  TH1F* htempfill[ntrees]={ 0 };
  int ctr=0;
  TString ltitle = ctitleprocess[6];
  if (forANPAS==0) l2D_2->AddEntry(hprocess[6], ltitle, "l");
  for (int prod=0; prod<ntrees; prod++){
    if (forANPAS==0 && (prod==0 || prod==5)) continue;
    if (forANPAS==1){
      ltitle = ctitleprocess[prod];
      ltitle.Append(ctitlemassprocess[prod]);
      if (prod==0 || prod>4) l2D->AddEntry(hprocess[prod], ltitle, "l");
      else l2D_2->AddEntry(hprocess[prod], ltitle, "l");
    }
    else{
      ltitle = ctitleprocess[prod];
      if (prod!=6) l2D_2->AddEntry(hprocess[prod], ltitle, "l");
    }
    if (ctr==0) hprocess[prod]->Draw("hist");
    else hprocess[prod]->Draw("histsame");
    htempfill[prod] = (TH1F*)hprocess[prod]->Clone(Form("%s_clone", hprocess[prod]->GetName()));
    htempfill[prod]->SetFillColor(hprocess[prod]->GetLineColor());
    htempfill[prod]->SetFillStyle(3002);
    if(forANPAS==1) htempfill[prod]->Draw("e2same");
    ctr++;
  }
  if (forANPAS==1){
    for (int prod=0; prod<4; prod++) hprocess_rewgt[prod]->Draw("histsame");
  }
  if (forANPAS==1) l2D->Draw("same");
  l2D_2->Draw("same");
  pt->Draw();

  c2D->RedrawAxis();
  c2D->Modified();
  c2D->Update();
  foutput->WriteTObject(c2D);

  TString canvasDir = coutput_common;
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

  for (int prod=0; prod<ntrees; prod++){
    if (htempfill[prod]!=0) delete htempfill[prod];
  }
  delete l2D_2;
  delete l2D;
  c2D->Close();
  delete pt;

  for (int prod = 0; prod < 4; prod++) delete hprocess_rewgt[prod];
  for (int prod = 0; prod < ntrees; prod++) delete hprocess[prod];
  foutput->Close();
}


void plotSystematics(int initialSyst=0, int finalSyst = nSystVars){
  for (int s=initialSyst; s<finalSyst; s++){
    for (int r=0; r<=2; r++){
      for (int e=7; e<=8; e++){
        for (int f=0; f<3; f++){
          produce_KDSystVariables(f, e, s, r);
        }
      }
    }
  }
  for (int s=initialSyst; s<finalSyst; s++){
    for (int r=0; r<=2; r++){
      for (int n=0; n<2; n++){
        for (int bb=0; bb<2; bb++){
          compare_KDShapeVariation_Untransformed(s, n, bb, r);
          for (int f=0; f<3; f++){
            compare_KDShapeVariation_Untransformed_perChannel(f, s, n, bb, r);
          }
        }
      }
    }
  }
  for (int s=initialSyst; s<finalSyst; s++){
    for (int p=0; p<=6; p++){
      for (int n=0; n<2; n++){
        compare_KDShapeVariation_SIPVariation(p, s, n);
      }
    }
  }
}

void generic_Histo2DPlotter(TFile* foutput, TString canvasDir, TH2* hSlice, int erg_tev, float lumi, int channel, int removeSIP, TGraph* tgdata = 0, TPaveText* ptextra = 0){
	gStyle->SetOptStat(0);
	gStyle->SetPalette(1);
	gStyle->SetPadRightMargin(0.20);
	gROOT->ForceStyle();
	foutput->cd();
	gStyle->SetOptStat(0);
	TString canvasname_2D = Form("c%s_%s_", hSlice->GetName(),user_folder[channel]);
  if (removeSIP==1) canvasname_2D = Form("c%s_NoSIP_%s_", hSlice->GetName(), user_folder[channel]);
  if (removeSIP==2) canvasname_2D = Form("c%s_NewSIP_%s_", hSlice->GetName(), user_folder[channel]);
  if (erg_tev>0){
    canvasname_2D.Append(Form("%iTeV", erg_tev));
  }
  else if (erg_tev<0){
    canvasname_2D.Append("AllTeV");
  }

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
	hSlice->GetXaxis()->SetTitleOffset(1);
	hSlice->GetYaxis()->SetTitleOffset(1.1);
  if ((hSlice->GetXaxis()->GetXmax() - hSlice->GetXaxis()->GetXmin())>3000) hSlice->GetXaxis()->SetNdivisions(505);
	hSlice->Draw("colz");
  if (tgdata!=0) tgdata->Draw("psame");
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
	TString cErgTev = Form("#font[42]{%.1f fb^{-1} (%i TeV)}",lumi,erg_tev);
  if (erg_tev<0) cErgTev = "#font[42]{19.7 fb^{-1} (8 TeV) + 5.1 fb^{-1} (7 TeV)}";
	if(channel==0) cErgTev.Prepend("4#mu ");
	if(channel==1) cErgTev.Prepend("4e ");
	if(channel==2) cErgTev.Prepend("2e+2#mu ");
  float posTextX = 0.575;
  if (channel<0 && erg_tev>0) posTextX = 0.6;
  else if (channel>=0 && channel<2 && erg_tev>0) posTextX = 0.575;
  else if (channel==2 && erg_tev>0) posTextX = 0.545;
  else if (channel>=0 && channel<2 && erg_tev<0) posTextX = 0.35;
  else if (channel==2 && erg_tev<0) posTextX = 0.32;
  else if (channel<0 && erg_tev<0) posTextX = 0.41;
  text = pt->AddText(posTextX, 0.45, cErgTev);
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
  if (ptextra!=0) ptextra->Draw("same");
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

void floorHistogram(TH2F* h){
  for (int binx=0; binx<=h->GetNbinsX()+1; binx++){
    for (int biny=0; biny<=h->GetNbinsY()+1; biny++){
      if (h->GetBinContent(binx, biny)==0) h->SetBinContent(binx, biny,1e-15);
    }
  }
}

void produce_KDSystVariables_2D(int folder, int erg_tev, int iSyst_1=0, int iSyst_2=1, int normOption=0, int removeSIP=2){
	gROOT->ProcessLine(".x tdrstyle.cc");
	gStyle->SetOptStat(0);
	gStyle->SetPalette(1);
	gStyle->SetPadRightMargin(0.20);
	gROOT->ForceStyle();

  // SIP cut scaling
  string processName[8] ={
    "CTau0", "CTau100", "CTau500", "CTau1000",
    "qqZZ", "ggZZ", "CR", "data"
  };

	char TREE_NAME[]="SelectedTree";
	TString erg_dir;
	erg_dir.Form("LHC_%iTeV",erg_tev);
	TString comstring;
	comstring.Form("%iTeV",erg_tev);
	int EnergyIndex=0;if(erg_tev==8) EnergyIndex=1;

  TString OUTPUT_NAME = "LifetimeKD_SystVariables_DataQQBZZ_2D_";
  if (removeSIP==1) OUTPUT_NAME.Append("noSIP_");
  else if (removeSIP==2) OUTPUT_NAME.Append("newSIP_");
	OUTPUT_NAME.Append(user_folder[folder]);
	OUTPUT_NAME.Append(Form("_%s_",cSystVariable[iSyst_1]));
	OUTPUT_NAME.Append(Form("_%s_",cSystVariable[iSyst_2]));
	OUTPUT_NAME.Append(comstring);
	if(normOption==1) OUTPUT_NAME.Append("_Conditional");
	OUTPUT_NAME.Append(".root");

	float ggzz_offshell_sum=0,ggzz_80_100_sum=0;
	float zx_offshell_sum=0,zx_80_100_sum=0,zx_signal_sum=0;
  float qqzz_signal_sum=0, ggzz_signal_sum=0;
  float higgs_signal_sum[4]={ 0 };

  string cinput_aux = user_dir_hep + "Analysis/Auxiliary/";
  TString SSinputname = "OSoverSS_";
  SSinputname.Append("MCCR_");
  SSinputname += comstring;
  SSinputname.Append(".root");
  SSinputname.Prepend(cinput_aux.c_str());
  TFile* fSSinput = new TFile(SSinputname, "read");
  TH2F* hRatio;
  TString chratio = "hCR_MC_OSSSRatios_";
  if (removeSIP==0) chratio.Append("Old cut");
  else if (removeSIP==1) chratio.Append("No cut");
  else if (removeSIP==2) chratio.Append("New cut");
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


  float ZXfake_weight_SS;
  float ZXfake_weight_OS[6];

  float MC_weight;
  float MC_weight_noxsec;
  float MC_weight_QQZZEWK=1;
  float MC_weight_Kfactor=1;
  float MC_weight_QQBGGProper[4];
  float ZZMass, ZZPt, ZZEta, ZZPhi;
  float KalmanCandVtx_x, KalmanCandVtx_y, KalmanCandVtx_z, KalmanCandVtx_chi2;
  float OfflinePrimaryVtx_x, OfflinePrimaryVtx_y, OfflinePrimaryVtx_z, OfflinePrimaryVtx_ndof;
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
  float Txy, Dxy, Txy_true, Dxy_true, delTxy, delDxy, KD1, KD2;
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
	TString coutput_common = user_dir_hep + "Analysis/Plots/";
	coutput_common.Append(Form("%s_vs_%s/",cSystVariable[iSyst_1],cSystVariable[iSyst_2]));
	TString mkdirCommand = "mkdir -p ";
	mkdirCommand.Append(coutput_common);
	gSystem->Exec(mkdirCommand);
	TString coutput = coutput_common + OUTPUT_NAME;
	TFile* foutput = new TFile(coutput,"recreate");

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
  }
  TChain* tggzz = new TChain(TREE_NAME);
  for (int smp = kGGOLDSamples; smp < kGGMCFMSamples; smp++){
    TString cinput_ggzz = cinput_ggzz_common_noSIP + sample_FullSim[smp] + "_Reprocessed.root";
    tggzz->Add(cinput_ggzz);
  }
  TChain* tzx = new TChain(TREE_NAME);
  TString cinput_zx = cinput_zx_common_noSIP;
  cinput_zx = cinput_zx + "CR/" + sample_FullSim[kAllSamples - 1] + ".root";
  tzx->Add(cinput_zx);

  TString cinput_data = cinput_common_noSIP + user_folder[3] + "/" + data_files[folder] + ".root";
  TChain* tdata = new TChain(TREE_NAME);
  tdata->Add(cinput_data);
  TChain* tc[ntrees] ={ tdata, tqqzz, tggzz, tzx };

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
    TString cinput_higgs = cinput_higgs_common_noSIP + sample_FullSim[smp] + ".root";
    cout << "Attaching Higgs under " << cinput_higgs << endl;
    thiggs[smp] = new TChain(TREE_NAME);
    thiggs[smp]->Add(cinput_higgs);
  }
  TChain* tc_extra[ntrees_extra] ={ thiggs[0], tqqzz_extra, tggzz_extra, thiggs[1], thiggs[2], thiggs[3] };
	
	int nbins_KD1 = nbins_SystVariable[iSyst_1];
	double KD1_limits[2] = { KD_limits_SystVariable[iSyst_1][0], KD_limits_SystVariable[iSyst_1][1] };
	int nbins_KD2 = nbins_SystVariable[iSyst_2]/2;
  if (nbins_SystVariable[iSyst_2] % 2 == 1) nbins_KD2++;
	double KD2_limits[2] = { KD_limits_SystVariable[iSyst_2][0], KD_limits_SystVariable[iSyst_2][1] };

  TGraph* tgdata[nMZZ_ranges];
  std::vector<float> dataX[nMZZ_ranges];
  std::vector<float> dataY[nMZZ_ranges];
  float* dataArrayX[nMZZ_ranges];
  float* dataArrayY[nMZZ_ranges];

	TH3F* hdata_full = new TH3F("hdata_full", Form("%i TeV %s Data", erg_tev, user_folder[folder]), nMZZ_ranges, 0, nMZZ_ranges, nbins_KD1, KD1_limits[0], KD1_limits[1], nbins_KD2, KD2_limits[0], KD2_limits[1]);
	hdata_full->GetYaxis()->SetTitle(cSystVariable_label[iSyst_1]);
	hdata_full->GetZaxis()->SetTitle(cSystVariable_label[iSyst_2]);
	hdata_full->GetXaxis()->SetTitle("m_{4l} (GeV)");
	hdata_full->GetXaxis()->SetBinLabel(1,"70 - 105.6");
	hdata_full->GetXaxis()->SetBinLabel(2,"140.6 - 170");
	hdata_full->GetXaxis()->SetBinLabel(3,"170 - 800");
	hdata_full->GetXaxis()->SetBinLabel(4,"105.6 - 140.6");
	hdata_full->SetOption("colz");
	TH3F* hqqzz_full = (TH3F*) hdata_full->Clone();
	TH3F* hggzz_full = (TH3F*) hdata_full->Clone();
	TH3F* hzx_full = (TH3F*) hdata_full->Clone();
	TH3F* hmc_full = (TH3F*) hdata_full->Clone();
	hqqzz_full->SetNameTitle("hqqzz_full",Form("%i TeV %s q#bar{q}#rightarrow4l Background",erg_tev,user_folder[folder]));
	hggzz_full->SetNameTitle("hggzz_full",Form("%i TeV %s gg Background",erg_tev,user_folder[folder]));
	hzx_full->SetNameTitle("hzx_full",Form("%i TeV %s Z+X Background",erg_tev,user_folder[folder]));
	hmc_full->SetNameTitle("hmc_full",Form("%i TeV %s Total Background",erg_tev,user_folder[folder]));
	TH3F* hfull[ntrees+1] = { hdata_full, hqqzz_full, hggzz_full, hzx_full, hmc_full };
	for (int tt = 0; tt < ntrees; tt++) hfull[tt]->Sumw2();
	TH2F* hproj[ntrees+1];

  TH2F* hhiggs[kGGHSamples];
  for (int hh = 0; hh < kGGHSamples; hh++){
    hhiggs[hh] = (TH2F*)hdata_full->Project3D("zy");
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
        //        }
      }
      if (tt == 2 || tt == 1){
        tc[tt]->SetBranchAddress("MC_weight_Kfactor", &MC_weight_Kfactor);
        tc[tt]->SetBranchAddress("MC_weight_QQBGGProper", MC_weight_QQBGGProper);
        if (tt == 1){
          tc[tt]->SetBranchAddress("MC_weight_QQZZEWK", &MC_weight_QQZZEWK);
        }
      }
    }

    tc[tt]->SetBranchAddress("ZZMass", &ZZMass);
    tc[tt]->SetBranchAddress("ZZPt", &ZZPt);
    tc[tt]->SetBranchAddress("ZZEta", &ZZEta);
    tc[tt]->SetBranchAddress("ZZPhi", &ZZPhi);
    tc[tt]->SetBranchAddress("OfflinePrimaryVtx_x", &OfflinePrimaryVtx_x);
    tc[tt]->SetBranchAddress("OfflinePrimaryVtx_y", &OfflinePrimaryVtx_y);
    tc[tt]->SetBranchAddress("OfflinePrimaryVtx_z", &OfflinePrimaryVtx_z);
    tc[tt]->SetBranchAddress("OfflinePrimaryVtx_ndof", &OfflinePrimaryVtx_ndof);
    tc[tt]->SetBranchAddress("KalmanCandVtx_x", &KalmanCandVtx_x);
    tc[tt]->SetBranchAddress("KalmanCandVtx_y", &KalmanCandVtx_y);
    tc[tt]->SetBranchAddress("KalmanCandVtx_z", &KalmanCandVtx_z);
    tc[tt]->SetBranchAddress("KalmanCandVtx_chi2", &KalmanCandVtx_chi2);
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
    }
    if (tt == 1){
      tc_extra[tt]->SetBranchAddress("MC_weight_QQZZEWK", &MC_weight_QQZZEWK);
    }
    tc_extra[tt]->SetBranchAddress("ZZMass", &ZZMass);
    tc_extra[tt]->SetBranchAddress("ZZPt", &ZZPt);
    tc_extra[tt]->SetBranchAddress("ZZEta", &ZZEta);
    tc_extra[tt]->SetBranchAddress("ZZPhi", &ZZPhi);
    tc_extra[tt]->SetBranchAddress("OfflinePrimaryVtx_x", &OfflinePrimaryVtx_x);
    tc_extra[tt]->SetBranchAddress("OfflinePrimaryVtx_y", &OfflinePrimaryVtx_y);
    tc_extra[tt]->SetBranchAddress("OfflinePrimaryVtx_z", &OfflinePrimaryVtx_z);
    tc_extra[tt]->SetBranchAddress("OfflinePrimaryVtx_ndof", &OfflinePrimaryVtx_ndof);
    tc_extra[tt]->SetBranchAddress("KalmanCandVtx_x", &KalmanCandVtx_x);
    tc_extra[tt]->SetBranchAddress("KalmanCandVtx_y", &KalmanCandVtx_y);
    tc_extra[tt]->SetBranchAddress("KalmanCandVtx_z", &KalmanCandVtx_z);
    tc_extra[tt]->SetBranchAddress("KalmanCandVtx_chi2", &KalmanCandVtx_chi2);
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
      if (removeSIP==2 && (fabs(Lep1_Z1SIP) >= 4 || fabs(Lep2_Z1SIP) >= 4 || fabs(Lep3_Z1SIP) >= 5 || fabs(Lep4_Z1SIP) >= 5 || KalmanCandVtx_chi2 >= 30)) continue;
      if (removeSIP==0 && (fabs(Lep1SIP) >= 4 || fabs(Lep2SIP) >= 4 || fabs(Lep3SIP) >= 4 || fabs(Lep4SIP) >= 4)) continue;
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
        if (removeSIP==0 && tc[tt]->GetBranchStatus("ZXfake_weight_OS")) MC_weight = ZXfake_weight_OS[0] * (targetCRcount[0][0][regionIndex] / (targetCRcount[0][0][regionIndex] + targetCRcount[1][0][regionIndex]));
        else if (removeSIP==1 && tc[tt]->GetBranchStatus("ZXfake_weight_OS")) MC_weight = ZXfake_weight_OS[1] * (targetCRcount[0][1][regionIndex] / (targetCRcount[0][1][regionIndex] + targetCRcount[1][1][regionIndex]));
        else if (removeSIP==2 && tc[tt]->GetBranchStatus("ZXfake_weight_OS")) MC_weight = ZXfake_weight_OS[4] * (targetCRcount[0][4][regionIndex] / (targetCRcount[0][4][regionIndex] + targetCRcount[1][4][regionIndex]));
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

          //          if (removeSIP==2 && tt==3 && folder==0 && ZZMass>=110 && ZZMass<130) cout << ZZMass << '\t' << SSwgt_OSoverSS<< '\t' << ZXfake_weight_SS << '\t' << (targetCRcount[1][4][regionIndex] / (targetCRcount[0][4][regionIndex] + targetCRcount[1][4][regionIndex])) << endl;

          MC_weight = ZXfake_weight_SS*SSwgt_OSoverSS;
          if (removeSIP==0) MC_weight *= (targetCRcount[1][0][regionIndex] / (targetCRcount[0][0][regionIndex] + targetCRcount[1][0][regionIndex]));
          else if (removeSIP==1) MC_weight *= (targetCRcount[1][1][regionIndex] / (targetCRcount[0][1][regionIndex] + targetCRcount[1][1][regionIndex]));
          else if (removeSIP==2) MC_weight *= (targetCRcount[1][4][regionIndex] / (targetCRcount[0][4][regionIndex] + targetCRcount[1][4][regionIndex]));
        }
      }

      if ((iSyst_1>=22 && iSyst_1<nSystVars)||(iSyst_2>=22 && iSyst_2<nSystVars)){ // BE CAREFUL HERE, HARDCODED BOUNDS AHEAD
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

      if (iSyst_1==0 || iSyst_1==20) KD1 = Txy;
      if (iSyst_1==1) KD1 = Dxy;
      if (iSyst_1==2) KD1 = delTxy;
      if (iSyst_1==3) KD1 = delDxy;
      if (iSyst_1==12) KD1 = (OfflinePrimaryVtx_x - GenPrimaryVtx_x)*cos(GenHPhi) + (OfflinePrimaryVtx_y - GenPrimaryVtx_y)*sin(GenHPhi);
      if (iSyst_1==13) KD1 = (KalmanCandVtx_x - GenIntVtx_x)*cos(GenHPhi) + (KalmanCandVtx_y - GenIntVtx_y)*sin(GenHPhi);
      if (iSyst_1==14) KD1 = (KalmanCandVtx_x-GenIntVtx_x);
      if (iSyst_1==15) KD1 = (KalmanCandVtx_y-GenIntVtx_y);
      if (iSyst_1==16) KD1 = (KalmanCandVtx_z-GenIntVtx_z);
      if (iSyst_1==17) KD1 = (OfflinePrimaryVtx_x - GenPrimaryVtx_x);
      if (iSyst_1==18) KD1 = (OfflinePrimaryVtx_y - GenPrimaryVtx_y);
      if (iSyst_1==19) KD1 = (OfflinePrimaryVtx_z - GenPrimaryVtx_z);
      if (iSyst_1==22) KD1 = (OfflinePrimaryVtx_x - BeamPosX);
      if (iSyst_1==23) KD1 = (OfflinePrimaryVtx_y - BeamPosY);
      if (iSyst_1==24) KD1 = (OfflinePrimaryVtx_z - BeamPosZ)/1000.;
      if (iSyst_1==25) KD1 = (GenPrimaryVtx_x - BeamPosX);
      if (iSyst_1==26) KD1 = (GenPrimaryVtx_y - BeamPosY);
      if (iSyst_1==27) KD1 = (GenPrimaryVtx_z - BeamPosZ)/1000.;
      if (iSyst_1==28) KD1 = Txy_BS;
      if (iSyst_1==29) KD1 = Txy_BS - Txy;
      if (iSyst_1==30) KD1 = Txy_BS - Txy_true;
      KD1 = KD1*10000.; // in 1 um
      if (iSyst_1==20){
        KD1 = tanh(KD1/800.); if (KD1==1.) KD1=0.999999; if (KD1==-1.) KD1=-0.999999;
      }
      if (iSyst_1==21){
        KD1 = D_bkg; if (KD1==1.) KD1=0.999999;
      }
      if (iSyst_1==4) KD1 = ZZPt/ZZMass;
      if (iSyst_1==5) KD1 = (Txy-Txy_true)/delTxy;
      if (iSyst_1==6) KD1 = (Dxy-Dxy_true)/delDxy;
      if (iSyst_1==7) KD1 = ((OfflinePrimaryVtx_x - GenPrimaryVtx_x)*cos(GenHPhi) + (OfflinePrimaryVtx_y - GenPrimaryVtx_y)*sin(GenHPhi))/sigmaPV_xy;
      if (iSyst_1==8) KD1 = ((KalmanCandVtx_x - GenIntVtx_x)*cos(GenHPhi) + (KalmanCandVtx_y - GenIntVtx_y)*sin(GenHPhi))/sigmaInt_xy;
      if (iSyst_1==9) KD1 = (KalmanCandVtx_x-GenIntVtx_x)/sqrt(KalmanCandVtx_cov_xx);
      if (iSyst_1==10) KD1 = (KalmanCandVtx_y-GenIntVtx_y)/sqrt(KalmanCandVtx_cov_yy);
      if (iSyst_1==11) KD1 = (KalmanCandVtx_z-GenIntVtx_z)/sqrt(KalmanCandVtx_cov_zz);
      if (iSyst_1==31) KD1 = (OfflinePrimaryVtx_ndof+3.)/2.;

      if (iSyst_2==0 || iSyst_2==20) KD2 = Txy;
      if (iSyst_2==1) KD2 = Dxy;
      if (iSyst_2==2) KD2 = delTxy;
      if (iSyst_2==3) KD2 = delDxy;
      if (iSyst_2==12) KD2 = (OfflinePrimaryVtx_x - GenPrimaryVtx_x)*cos(GenHPhi) + (OfflinePrimaryVtx_y - GenPrimaryVtx_y)*sin(GenHPhi);
      if (iSyst_2==13) KD2 = (KalmanCandVtx_x - GenIntVtx_x)*cos(GenHPhi) + (KalmanCandVtx_y - GenIntVtx_y)*sin(GenHPhi);
      if (iSyst_2==14) KD2 = (KalmanCandVtx_x-GenIntVtx_x);
      if (iSyst_2==15) KD2 = (KalmanCandVtx_y-GenIntVtx_y);
      if (iSyst_2==16) KD2 = (KalmanCandVtx_z-GenIntVtx_z);
      if (iSyst_2==17) KD2 = (OfflinePrimaryVtx_x - GenPrimaryVtx_x);
      if (iSyst_2==18) KD2 = (OfflinePrimaryVtx_y - GenPrimaryVtx_y);
      if (iSyst_2==19) KD2 = (OfflinePrimaryVtx_z - GenPrimaryVtx_z);
      if (iSyst_2==22) KD2 = (OfflinePrimaryVtx_x - BeamPosX);
      if (iSyst_2==23) KD2 = (OfflinePrimaryVtx_y - BeamPosY);
      if (iSyst_2==24) KD2 = (OfflinePrimaryVtx_z - BeamPosZ)/1000.;
      if (iSyst_2==25) KD2 = (GenPrimaryVtx_x - BeamPosX);
      if (iSyst_2==26) KD2 = (GenPrimaryVtx_y - BeamPosY);
      if (iSyst_2==27) KD2 = (GenPrimaryVtx_z - BeamPosZ)/1000.;
      if (iSyst_2==28) KD2 = Txy_BS;
      if (iSyst_2==29) KD2 = Txy_BS - Txy;
      if (iSyst_2==30) KD2 = Txy_BS - Txy_true;
      KD2 = KD2*10000.; // in 1 um
      if (iSyst_2==20){
        KD2 = tanh(KD2/800.); if (KD2==1.) KD2=0.999999; if (KD2==-1.) KD2=-0.999999;
      }
      if (iSyst_2==21){
        KD2 = D_bkg; if (KD2==1.) KD2=0.999999;
      }
      if (iSyst_2==4) KD2 = ZZPt/ZZMass;
      if (iSyst_2==5) KD2 = (Txy-Txy_true)/delTxy;
      if (iSyst_2==6) KD2 = (Dxy-Dxy_true)/delDxy;
      if (iSyst_2==7) KD2 = ((OfflinePrimaryVtx_x - GenPrimaryVtx_x)*cos(GenHPhi) + (OfflinePrimaryVtx_y - GenPrimaryVtx_y)*sin(GenHPhi))/sigmaPV_xy;
      if (iSyst_2==8) KD2 = ((KalmanCandVtx_x - GenIntVtx_x)*cos(GenHPhi) + (KalmanCandVtx_y - GenIntVtx_y)*sin(GenHPhi))/sigmaInt_xy;
      if (iSyst_2==9) KD2 = (KalmanCandVtx_x-GenIntVtx_x)/sqrt(KalmanCandVtx_cov_xx);
      if (iSyst_2==10) KD2 = (KalmanCandVtx_y-GenIntVtx_y)/sqrt(KalmanCandVtx_cov_yy);
      if (iSyst_2==11) KD2 = (KalmanCandVtx_z-GenIntVtx_z)/sqrt(KalmanCandVtx_cov_zz);
      if (iSyst_2==31) KD2 = (OfflinePrimaryVtx_ndof+3.)/2.;

			float wgt = 1;

			int biny = hfull[tt]->GetYaxis()->FindBin(KD1);
			int binz = hfull[tt]->GetZaxis()->FindBin(KD2);
			for (int binx = 0; binx < nMZZ_ranges; binx++){
        if (tt == 3 && binx == 0 && ZZMass >= systZZMass_range[binx][0] + 30 && ZZMass < systZZMass_range[binx][1] + 30){ // ZX special fill
					wgt = MC_weight;
					hfull[tt]->AddBinContent(hfull[tt]->GetBin(binx + 1, biny,binz), wgt);
          hfull[tt]->SetBinError(hfull[tt]->GetBin(binx + 1, biny,binz), sqrt(wgt*wgt + pow(hfull[tt]->GetBinError(hfull[tt]->GetBin(binx + 1, biny,binz)), 2)));
          if (ZZMass<130 && ZZMass>=110){
            zx_80_100_sum += wgt;
          }
          continue;
        }
        else if (tt == 3 && binx == 0) continue;
        if (tt == 2 && binx == 0) continue;
        if (tt == 1 && binx == 0 && ZZMass >= systZZMass_range[binx][0] && ZZMass < systZZMass_range[binx][1]){ // gg bkg special fill
          wgt = MC_weight*MC_weight_Kfactor*MC_weight_QQBGGProper[0];
          hfull[2]->AddBinContent(hfull[2]->GetBin(binx + 1, biny, binz), wgt);
          hfull[2]->SetBinError(hfull[2]->GetBin(binx + 1, biny, binz), sqrt(wgt*wgt + pow(hfull[2]->GetBinError(hfull[2]->GetBin(binx + 1, biny, binz)), 2)));
          if (ZZMass>=80 && ZZMass<100) ggzz_80_100_sum += wgt;
        }
				if (binx < nMZZ_ranges - 1){
					wgt = MC_weight;
          if (tt==0) wgt=1;
          if (tt == 1) wgt *= MC_weight_QQZZEWK;
					if (tt == 2) wgt *= MC_weight_Kfactor;
          if (ZZMass < systZZMass_range[binx][1] && ZZMass >= systZZMass_range[binx][0]){
            hfull[tt]->AddBinContent(hfull[tt]->GetBin(binx + 1, biny, binz), wgt);
            hfull[tt]->SetBinError(hfull[tt]->GetBin(binx + 1, biny, binz), sqrt(wgt*wgt + pow(hfull[tt]->GetBinError(hfull[tt]->GetBin(binx + 1, biny, binz)), 2)));
            if (tt==0){
              dataX[binx].push_back(KD1);
              dataY[binx].push_back(KD2);
            }
          }
          if (tt == 2 && ZZMass >= 220 && ZZMass<1600 && binx==2) ggzz_offshell_sum += wgt;
					if (tt == 3 && ZZMass >= 220 && ZZMass<1600 && binx==2) zx_offshell_sum += wgt;
				}
				else{ // Fill signal region
          if (ZZMass >= systZZMass_range[binx][1] || ZZMass < systZZMass_range[binx][0]) continue;
          wgt = (tt==3 ? MC_weight : MC_weight_noxsec);
          if (tt==0) wgt=1;
					if (tt == 1){
						wgt = MC_weight_noxsec*MC_weight_Kfactor*MC_weight_QQBGGProper[0];
						hfull[2]->AddBinContent(hfull[2]->GetBin(binx + 1, biny,binz), wgt);
            hfull[2]->SetBinError(hfull[2]->GetBin(binx + 1, biny, binz), sqrt(wgt*wgt + pow(hfull[2]->GetBinError(hfull[2]->GetBin(binx + 1, biny, binz)), 2)));
            ggzz_signal_sum += wgt;
						wgt = MC_weight_noxsec*MC_weight_QQBGGProper[1];
						hfull[1]->AddBinContent(hfull[1]->GetBin(binx + 1, biny,binz), wgt);
            hfull[1]->SetBinError(hfull[1]->GetBin(binx + 1, biny, binz), sqrt(wgt*wgt + pow(hfull[1]->GetBinError(hfull[1]->GetBin(binx + 1, biny, binz)), 2)));
            qqzz_signal_sum += wgt;
					}
					else if (tt == 2){
						wgt *= MC_weight_Kfactor*MC_weight_QQBGGProper[0];
						hfull[tt]->AddBinContent(hfull[tt]->GetBin(binx + 1, biny,binz), wgt);
            hfull[tt]->SetBinError(hfull[tt]->GetBin(binx + 1, biny, binz), sqrt(wgt*wgt + pow(hfull[tt]->GetBinError(hfull[tt]->GetBin(binx + 1, biny, binz)), 2)));
            ggzz_signal_sum += wgt;
					}
          else{
            hfull[tt]->AddBinContent(hfull[tt]->GetBin(binx + 1, biny, binz), wgt);
            hfull[tt]->SetBinError(hfull[tt]->GetBin(binx + 1, biny, binz), sqrt(wgt*wgt + pow(hfull[tt]->GetBinError(hfull[tt]->GetBin(binx + 1, biny, binz)), 2)));
            if (tt==0){
              dataX[binx].push_back(KD1);
              dataY[binx].push_back(KD2);
            }
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
      if (removeSIP==0 && (fabs(Lep1SIP) >= 4 || fabs(Lep2SIP) >= 4 || fabs(Lep3SIP) >= 4 || fabs(Lep4SIP) >= 4)) continue;
      if (removeSIP==2 && (fabs(Lep1_Z1SIP) >= 4 || fabs(Lep2_Z1SIP) >= 4 || fabs(Lep3_Z1SIP) >= 5 || fabs(Lep4_Z1SIP) >= 5 || KalmanCandVtx_chi2 >= 30)) continue;

      int binx = nMZZ_ranges-1;
      if (ZZMass >= systZZMass_range[binx][1] || ZZMass < systZZMass_range[binx][0]) continue;

      if ((iSyst_1>=22 && iSyst_1<nSystVars)||(iSyst_2>=22 && iSyst_2<nSystVars)){ // BE CAREFUL HERE, HARDCODED BOUNDS AHEAD
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
      
      if (!tc_extra[tt]->GetBranchStatus("D_bkg")) D_bkg = p0plus_VAJHU*p0plus_m4l / (p0plus_VAJHU*p0plus_m4l + bkg_VAMCFM*bkg_m4l);

      if (iSyst_1==0 || iSyst_1==20) KD1 = Txy;
      if (iSyst_1==1) KD1 = Dxy;
      if (iSyst_1==2) KD1 = delTxy;
      if (iSyst_1==3) KD1 = delDxy;
      if (iSyst_1==12) KD1 = (OfflinePrimaryVtx_x - GenPrimaryVtx_x)*cos(GenHPhi) + (OfflinePrimaryVtx_y - GenPrimaryVtx_y)*sin(GenHPhi);
      if (iSyst_1==13) KD1 = (KalmanCandVtx_x - GenIntVtx_x)*cos(GenHPhi) + (KalmanCandVtx_y - GenIntVtx_y)*sin(GenHPhi);
      if (iSyst_1==14) KD1 = (KalmanCandVtx_x-GenIntVtx_x);
      if (iSyst_1==15) KD1 = (KalmanCandVtx_y-GenIntVtx_y);
      if (iSyst_1==16) KD1 = (KalmanCandVtx_z-GenIntVtx_z);
      if (iSyst_1==17) KD1 = (OfflinePrimaryVtx_x - GenPrimaryVtx_x);
      if (iSyst_1==18) KD1 = (OfflinePrimaryVtx_y - GenPrimaryVtx_y);
      if (iSyst_1==19) KD1 = (OfflinePrimaryVtx_z - GenPrimaryVtx_z);
      if (iSyst_1==22) KD1 = (OfflinePrimaryVtx_x - BeamPosX);
      if (iSyst_1==23) KD1 = (OfflinePrimaryVtx_y - BeamPosY);
      if (iSyst_1==24) KD1 = (OfflinePrimaryVtx_z - BeamPosZ)/1000.;
      if (iSyst_1==25) KD1 = (GenPrimaryVtx_x - BeamPosX);
      if (iSyst_1==26) KD1 = (GenPrimaryVtx_y - BeamPosY);
      if (iSyst_1==27) KD1 = (GenPrimaryVtx_z - BeamPosZ)/1000.;
      if (iSyst_1==28) KD1 = Txy_BS;
      if (iSyst_1==29) KD1 = Txy_BS - Txy;
      if (iSyst_1==30) KD1 = Txy_BS - Txy_true;
      KD1 = KD1*10000.; // in 1 um
      if (iSyst_1==20){
        KD1 = tanh(KD1/800.); if (KD1==1.) KD1=0.999999; if (KD1==-1.) KD1=-0.999999;
      }
      if (iSyst_1==21){
        KD1 = D_bkg; if (KD1==1.) KD1=0.999999;
      }
      if (iSyst_1==4) KD1 = ZZPt/ZZMass;
      if (iSyst_1==5) KD1 = (Txy-Txy_true)/delTxy;
      if (iSyst_1==6) KD1 = (Dxy-Dxy_true)/delDxy;
      if (iSyst_1==7) KD1 = ((OfflinePrimaryVtx_x - GenPrimaryVtx_x)*cos(GenHPhi) + (OfflinePrimaryVtx_y - GenPrimaryVtx_y)*sin(GenHPhi))/sigmaPV_xy;
      if (iSyst_1==8) KD1 = ((KalmanCandVtx_x - GenIntVtx_x)*cos(GenHPhi) + (KalmanCandVtx_y - GenIntVtx_y)*sin(GenHPhi))/sigmaInt_xy;
      if (iSyst_1==9) KD1 = (KalmanCandVtx_x-GenIntVtx_x)/sqrt(KalmanCandVtx_cov_xx);
      if (iSyst_1==10) KD1 = (KalmanCandVtx_y-GenIntVtx_y)/sqrt(KalmanCandVtx_cov_yy);
      if (iSyst_1==11) KD1 = (KalmanCandVtx_z-GenIntVtx_z)/sqrt(KalmanCandVtx_cov_zz);
      if (iSyst_1==31) KD1 = (OfflinePrimaryVtx_ndof+3.)/2.;

      if (iSyst_2==0 || iSyst_2==20) KD2 = Txy;
      if (iSyst_2==1) KD2 = Dxy;
      if (iSyst_2==2) KD2 = delTxy;
      if (iSyst_2==3) KD2 = delDxy;
      if (iSyst_2==12) KD2 = (OfflinePrimaryVtx_x - GenPrimaryVtx_x)*cos(GenHPhi) + (OfflinePrimaryVtx_y - GenPrimaryVtx_y)*sin(GenHPhi);
      if (iSyst_2==13) KD2 = (KalmanCandVtx_x - GenIntVtx_x)*cos(GenHPhi) + (KalmanCandVtx_y - GenIntVtx_y)*sin(GenHPhi);
      if (iSyst_2==14) KD2 = (KalmanCandVtx_x-GenIntVtx_x);
      if (iSyst_2==15) KD2 = (KalmanCandVtx_y-GenIntVtx_y);
      if (iSyst_2==16) KD2 = (KalmanCandVtx_z-GenIntVtx_z);
      if (iSyst_2==17) KD2 = (OfflinePrimaryVtx_x - GenPrimaryVtx_x);
      if (iSyst_2==18) KD2 = (OfflinePrimaryVtx_y - GenPrimaryVtx_y);
      if (iSyst_2==19) KD2 = (OfflinePrimaryVtx_z - GenPrimaryVtx_z);
      if (iSyst_2==22) KD2 = (OfflinePrimaryVtx_x - BeamPosX);
      if (iSyst_2==23) KD2 = (OfflinePrimaryVtx_y - BeamPosY);
      if (iSyst_2==24) KD2 = (OfflinePrimaryVtx_z - BeamPosZ)/1000.;
      if (iSyst_2==25) KD2 = (GenPrimaryVtx_x - BeamPosX);
      if (iSyst_2==26) KD2 = (GenPrimaryVtx_y - BeamPosY);
      if (iSyst_2==27) KD2 = (GenPrimaryVtx_z - BeamPosZ)/1000.;
      if (iSyst_2==28) KD2 = Txy_BS;
      if (iSyst_2==29) KD2 = Txy_BS - Txy;
      if (iSyst_2==30) KD2 = Txy_BS - Txy_true;
      KD2 = KD2*10000.; // in 1 um
      if (iSyst_2==20){
        KD2 = tanh(KD2/800.); if (KD2==1.) KD2=0.999999; if (KD2==-1.) KD2=-0.999999;
      }
      if (iSyst_2==21){
        KD2 = D_bkg; if (KD2==1.) KD2=0.999999;
      }
      if (iSyst_2==4) KD2 = ZZPt/ZZMass;
      if (iSyst_2==5) KD2 = (Txy-Txy_true)/delTxy;
      if (iSyst_2==6) KD2 = (Dxy-Dxy_true)/delDxy;
      if (iSyst_2==7) KD2 = ((OfflinePrimaryVtx_x - GenPrimaryVtx_x)*cos(GenHPhi) + (OfflinePrimaryVtx_y - GenPrimaryVtx_y)*sin(GenHPhi))/sigmaPV_xy;
      if (iSyst_2==8) KD2 = ((KalmanCandVtx_x - GenIntVtx_x)*cos(GenHPhi) + (KalmanCandVtx_y - GenIntVtx_y)*sin(GenHPhi))/sigmaInt_xy;
      if (iSyst_2==9) KD2 = (KalmanCandVtx_x-GenIntVtx_x)/sqrt(KalmanCandVtx_cov_xx);
      if (iSyst_2==10) KD2 = (KalmanCandVtx_y-GenIntVtx_y)/sqrt(KalmanCandVtx_cov_yy);
      if (iSyst_2==11) KD2 = (KalmanCandVtx_z-GenIntVtx_z)/sqrt(KalmanCandVtx_cov_zz);
      if (iSyst_2==31) KD2 = (OfflinePrimaryVtx_ndof+3.)/2.;

			float wgt = 1;

			int biny = hhiggs[0]->GetXaxis()->FindBin(KD1);
      int binz = hhiggs[0]->GetYaxis()->FindBin(KD2);
      wgt = ((tt==0 || tt>=3) ? MC_weight : MC_weight_noxsec);
      if (tt == 1){
				wgt = MC_weight_noxsec*MC_weight_Kfactor*MC_weight_QQBGGProper[0];
				hfull[2]->AddBinContent(hfull[2]->GetBin(binx + 1, biny,binz), wgt);
        hfull[2]->SetBinError(hfull[2]->GetBin(binx + 1, biny, binz), sqrt(wgt*wgt + pow(hfull[2]->GetBinError(hfull[2]->GetBin(binx + 1, biny, binz)), 2)));
        ggzz_signal_sum += wgt;
				wgt = MC_weight_noxsec*MC_weight_QQBGGProper[1];
				hfull[1]->AddBinContent(hfull[1]->GetBin(binx + 1, biny,binz), wgt);
        hfull[1]->SetBinError(hfull[1]->GetBin(binx + 1, biny, binz), sqrt(wgt*wgt + pow(hfull[1]->GetBinError(hfull[1]->GetBin(binx + 1, biny, binz)), 2)));
        qqzz_signal_sum += wgt;
			}
			else if (tt == 2){
				wgt *= MC_weight_Kfactor*MC_weight_QQBGGProper[0];
				hfull[tt]->AddBinContent(hfull[tt]->GetBin(binx + 1, biny,binz), wgt);
        hfull[tt]->SetBinError(hfull[tt]->GetBin(binx + 1, biny, binz), sqrt(wgt*wgt + pow(hfull[tt]->GetBinError(hfull[tt]->GetBin(binx + 1, biny, binz)), 2)));
        ggzz_signal_sum += wgt;
			}
			else if (tt == 0){
        hhiggs[tt]->AddBinContent(hhiggs[tt]->GetBin(biny, binz), wgt);
        hhiggs[tt]->SetBinError(hhiggs[tt]->GetBin(biny, binz), sqrt(wgt*wgt + pow(hhiggs[tt]->GetBinError(hhiggs[tt]->GetBin(biny, binz)), 2)));
        higgs_signal_sum[tt] += wgt;
			}
      else if (tt >= 3){
        hhiggs[tt-2]->AddBinContent(hhiggs[tt-2]->GetBin(biny, binz), wgt);
        hhiggs[tt-2]->SetBinError(hhiggs[tt-2]->GetBin(biny, binz), sqrt(wgt*wgt + pow(hhiggs[tt-2]->GetBinError(hhiggs[tt-2]->GetBin(biny, binz)), 2)));
        higgs_signal_sum[tt-2] += wgt;
      }
    }
	}

	for (int biny = 0; biny <= hggzz_full->GetNbinsY()+1; biny++){
		for (int binz = 0; binz <= hggzz_full->GetNbinsZ() + 1; binz++){
			hggzz_full->SetBinContent(1, biny, binz, (hggzz_full->GetBinContent(1, biny, binz))*(yield_80_100_ggzz[EnergyIndex][folder] / ggzz_80_100_sum / luminosity[EnergyIndex]));
			hggzz_full->SetBinContent(4, biny, binz, (hggzz_full->GetBinContent(4, biny, binz))*(yield_signal_ggzz[EnergyIndex][folder] / ggzz_signal_sum / luminosity[EnergyIndex]));
			hggzz_full->SetBinContent(3, biny, binz, (hggzz_full->GetBinContent(3, biny, binz))*(yield_offshell_ggzz[EnergyIndex][folder] / ggzz_offshell_sum / luminosity[EnergyIndex]));
			hggzz_full->SetBinContent(2, biny, binz, (hggzz_full->GetBinContent(2, biny, binz))*(yield_offshell_ggzz[EnergyIndex][folder] / ggzz_offshell_sum / luminosity[EnergyIndex]));

      hggzz_full->SetBinError(1, biny, binz, (hggzz_full->GetBinError(1, biny, binz))*(yield_80_100_ggzz[EnergyIndex][folder] / ggzz_80_100_sum / luminosity[EnergyIndex]));
      hggzz_full->SetBinError(4, biny, binz, (hggzz_full->GetBinError(4, biny, binz))*(yield_signal_ggzz[EnergyIndex][folder] / ggzz_signal_sum / luminosity[EnergyIndex]));
      hggzz_full->SetBinError(3, biny, binz, (hggzz_full->GetBinError(3, biny, binz))*(yield_offshell_ggzz[EnergyIndex][folder] / ggzz_offshell_sum / luminosity[EnergyIndex]));
      hggzz_full->SetBinError(2, biny, binz, (hggzz_full->GetBinError(2, biny, binz))*(yield_offshell_ggzz[EnergyIndex][folder] / ggzz_offshell_sum / luminosity[EnergyIndex]));
    }
	}
  if (removeSIP==2) cout << zx_80_100_sum << '\t' << hzx_full->Integral(1, 1, 0, hzx_full->GetNbinsY()+1, 0, hzx_full->GetNbinsZ()+1) << '\t' << yield_80_100_zx[EnergyIndex][folder] << '\t' << ratio_targetsignalyield[6][4][0]/ratio_targetsignalyield[6][0][0] << endl;
  for (int biny = 0; biny <= hzx_full->GetNbinsY() + 1; biny++){
		for (int binz = 0; binz <= hzx_full->GetNbinsZ() + 1; binz++){
			hzx_full->SetBinContent(1, biny, binz, (hzx_full->GetBinContent(1, biny, binz))*(yield_80_100_zx[EnergyIndex][folder] / zx_80_100_sum / luminosity[EnergyIndex]));
			hzx_full->SetBinContent(4, biny, binz, (hzx_full->GetBinContent(4, biny, binz))*(yield_signal_zx[EnergyIndex][folder] / zx_signal_sum / luminosity[EnergyIndex]));
			hzx_full->SetBinContent(3, biny, binz, (hzx_full->GetBinContent(3, biny, binz))*(yield_offshell_zx[EnergyIndex][folder] / zx_offshell_sum / luminosity[EnergyIndex]));
			hzx_full->SetBinContent(2, biny, binz, (hzx_full->GetBinContent(2, biny, binz))*(yield_offshell_zx[EnergyIndex][folder] / zx_offshell_sum / luminosity[EnergyIndex]));

      hzx_full->SetBinError(1, biny, binz, (hzx_full->GetBinError(1, biny, binz))*(yield_80_100_zx[EnergyIndex][folder] / zx_80_100_sum / luminosity[EnergyIndex]));
      hzx_full->SetBinError(4, biny, binz, (hzx_full->GetBinError(4, biny, binz))*(yield_signal_zx[EnergyIndex][folder] / zx_signal_sum / luminosity[EnergyIndex]));
      hzx_full->SetBinError(3, biny, binz, (hzx_full->GetBinError(3, biny, binz))*(yield_offshell_zx[EnergyIndex][folder] / zx_offshell_sum / luminosity[EnergyIndex]));
      hzx_full->SetBinError(2, biny, binz, (hzx_full->GetBinError(2, biny, binz))*(yield_offshell_zx[EnergyIndex][folder] / zx_offshell_sum / luminosity[EnergyIndex]));
    }
	}
  for (int biny = 0; biny <= hqqzz_full->GetNbinsY()+1; biny++){
		for (int binz = 0; binz <= hqqzz_full->GetNbinsZ() + 1; binz++){
      hqqzz_full->SetBinContent(4, biny, binz, (hqqzz_full->GetBinContent(4, biny, binz))*(yield_signal_qqzz[EnergyIndex][folder] / qqzz_signal_sum / luminosity[EnergyIndex]));
      hqqzz_full->SetBinError(4, biny, binz, (hqqzz_full->GetBinError(4, biny, binz))*(yield_signal_qqzz[EnergyIndex][folder] / qqzz_signal_sum / luminosity[EnergyIndex]));
      if (removeSIP >= 1){
        hqqzz_full->SetBinContent(4, biny, binz, (hqqzz_full->GetBinContent(4, biny, binz))*(ratio_targetsignalyield[4][1][3]/ratio_targetsignalyield[4][0][3]));
        hqqzz_full->SetBinError(4, biny, binz, (hqqzz_full->GetBinError(4, biny, binz))*(ratio_targetsignalyield[4][1][3]/ratio_targetsignalyield[4][0][3]));
        if (removeSIP == 2){
          hqqzz_full->SetBinContent(4, biny, binz, (hqqzz_full->GetBinContent(4, biny, binz))*(ratio_targetsignalyield[4][4][3]/ratio_targetsignalyield[4][1][3]));
          hqqzz_full->SetBinError(4, biny, binz, (hqqzz_full->GetBinError(4, biny, binz))*(ratio_targetsignalyield[4][4][3]/ratio_targetsignalyield[4][1][3]));
        }
      }
      for (int hh=0; hh<kGGHSamples; hh++){
        hhiggs[hh]->SetBinContent(biny, binz, (hhiggs[hh]->GetBinContent(biny, binz))*(yield_signal_higgs[EnergyIndex][folder]*(ratio_targetsignalyield[0][1][3] / ratio_targetsignalyield[0][0][3]) / higgs_signal_sum[hh] / luminosity[EnergyIndex]));
        hhiggs[hh]->SetBinError(biny, binz, (hhiggs[hh]->GetBinError(biny, binz))*(yield_signal_higgs[EnergyIndex][folder]*(ratio_targetsignalyield[0][1][3] / ratio_targetsignalyield[0][0][3]) / higgs_signal_sum[hh] / luminosity[EnergyIndex]));
      }
    }
	}
  for (int hh = 0; hh < kGGHSamples; hh++){
    if (removeSIP == 0) hhiggs[hh]->Scale(ratio_targetsignalyield[hh][0][3] / ratio_targetsignalyield[hh][1][3]);
    else if (removeSIP == 2) hhiggs[hh]->Scale(ratio_targetsignalyield[hh][4][3] / ratio_targetsignalyield[hh][1][3]);
  }
  if (removeSIP >= 1){
    for (int binx = 1; binx <= hggzz_full->GetNbinsX(); binx++){
      double scale = ratio_targetsignalyield[5][1][binx-1]/ratio_targetsignalyield[5][0][binx-1];
      if (removeSIP==2) scale = ratio_targetsignalyield[5][4][binx-1]/ratio_targetsignalyield[5][0][binx-1];

      for (int biny = 0; biny <= hggzz_full->GetNbinsY()+1; biny++){
        for (int binz = 0; binz <= hggzz_full->GetNbinsZ()+1; binz++){
          double bincontent = hggzz_full->GetBinContent(binx, biny, binz);
          double binerror = hggzz_full->GetBinError(binx, biny, binz);

          bincontent *= scale;
          binerror *= scale;
          hggzz_full->SetBinContent(binx, biny, binz, bincontent);
          hggzz_full->SetBinError(binx, biny, binz, binerror);
        }
      }
    }
    for (int binx = 1; binx <= hzx_full->GetNbinsX(); binx++){
      double scale = ratio_targetsignalyield[6][1][binx-1]/ratio_targetsignalyield[6][0][binx-1];
      if (removeSIP==2) scale = ratio_targetsignalyield[6][4][binx-1]/ratio_targetsignalyield[6][0][binx-1];

      for (int biny = 0; biny <= hzx_full->GetNbinsY()+1; biny++){
        for (int binz = 0; binz <= hzx_full->GetNbinsZ()+1; binz++){
          double bincontent = hzx_full->GetBinContent(binx, biny, binz);
          double binerror = hzx_full->GetBinError(binx, biny, binz);

          bincontent *= scale;
          binerror *= scale;
          hzx_full->SetBinContent(binx, biny, binz, bincontent);
          hzx_full->SetBinError(binx, biny, binz, binerror);
        }
      }
    }
  }

  cout << "Signal gg norm: " << hggzz_full->Integral(4, 4, 0, hggzz_full->GetNbinsY()+1, 0, hggzz_full->GetNbinsZ()+1)*luminosity[EnergyIndex] << "\tqq norm: " << hqqzz_full->Integral(4, 4, 0, hqqzz_full->GetNbinsY()+1, 0, hqqzz_full->GetNbinsZ()+1)*luminosity[EnergyIndex]
    << "\tZX norm: " << hzx_full->Integral(4, 4, 0, hzx_full->GetNbinsY()+1, 0, hzx_full->GetNbinsZ()+1)*luminosity[EnergyIndex] << "\tHiggs norm: " << hhiggs[0]->Integral(0, hhiggs[0]->GetNbinsX()+1, 0, hhiggs[0]->GetNbinsY()+1)*luminosity[EnergyIndex] << "\tHiggs (BSM) norm: " << hhiggs[3]->Integral(0, hhiggs[3]->GetNbinsX()+1, 0, hhiggs[3]->GetNbinsY()+1)*luminosity[EnergyIndex] << endl;

  for (int sb=1; sb<=3; sb++) cout << "SB" << sb << " gg norm: " << hggzz_full->Integral(sb, sb, 0, hggzz_full->GetNbinsY()+1, 0, hggzz_full->GetNbinsZ()+1)*luminosity[EnergyIndex]
    << "\tqq norm: " << hqqzz_full->Integral(sb, sb, 0, hqqzz_full->GetNbinsY()+1, 0, hqqzz_full->GetNbinsZ()+1)*luminosity[EnergyIndex]
    << "\tZX norm: " << hzx_full->Integral(sb, sb, 0, hzx_full->GetNbinsY()+1, 0, hzx_full->GetNbinsZ()+1)*luminosity[EnergyIndex] << endl;

  for (int hh = 0; hh < kGGHSamples; hh++){
    hhiggs[hh]->Scale(luminosity[EnergyIndex]);
  }
  hqqzz_full->Scale(luminosity[EnergyIndex]);
	hggzz_full->Scale(luminosity[EnergyIndex]);
	hzx_full->Scale(luminosity[EnergyIndex]);

	foutput->cd();

  for (int binx=0; binx<nMZZ_ranges; binx++){
    int ndatabins = dataX[binx].size();
    cout << "Observed: " << ndatabins << " for region " << binx << endl;
    dataArrayX[binx] = new float[ndatabins];
    dataArrayY[binx] = new float[ndatabins];
    for (int db=0; db<ndatabins; db++){
      dataArrayX[binx][db] = dataX[binx].at(db);
      dataArrayY[binx][db] = dataY[binx].at(db);
    }
    tgdata[binx] = new TGraph(ndatabins, dataArrayX[binx], dataArrayY[binx]);
    tgdata[binx]->SetName(binx<nMZZ_ranges-1 ? Form("%s_SidebandRegion%i", "tgdata", binx+1) : Form("%s_SignalRegion", "tgdata"));
    foutput->WriteTObject(tgdata[binx]);
  }

  for (int tt = 1; tt < ntrees; tt++){
    if (tt==3 && ((iSyst_1>=7 && iSyst_1<20)||(iSyst_2>=7 && iSyst_2<20))) continue;
    if (tt==3 && (iSyst_1==27 || iSyst_2==27)) continue;
    hfull[ntrees]->Add(hfull[tt]);
  }
	for (int tt = 0; tt < ntrees+1; tt++){
		if (tt == ntrees){
      for (int hh = 0; hh < kGGHSamples; hh++){
        if (normOption == 1){
          hhiggs[hh]->SetName(Form("%s_Conditional", hhiggs[hh]->GetName()));
          hhiggs[hh]->SetTitle(Form("%s (Conditional in %s)", hhiggs[hh]->GetTitle(), hhiggs[hh]->GetYaxis()->GetTitle()));
          for (int biny = 1; biny <= hhiggs[hh]->GetNbinsY(); biny++){
            double integralY = hhiggs[hh]->Integral(0, hhiggs[hh]->GetNbinsX()+1, biny, biny);
            int nConditionalized=0;
            for (int binx = 1; binx <= hhiggs[hh]->GetNbinsX(); binx++){
              double bincontent = hhiggs[hh]->GetBinContent(binx, biny);
              if (normOption==1 && integralY!=0){
                bincontent/=integralY;
                if (bincontent>0) nConditionalized++;
              }
              hhiggs[hh]->SetBinContent(binx, biny, bincontent);
            }
            if (normOption==1 && integralY!=0 && nConditionalized<=10){
              for (int binx = 1; binx <= hhiggs[hh]->GetNbinsX(); binx++) hhiggs[hh]->SetBinContent(binx, biny, 0);
            }
          }
        }
        foutput->WriteTObject(hhiggs[hh]);
        generic_Histo2DPlotter(foutput, coutput_common, hhiggs[hh], erg_tev, luminosity[EnergyIndex], folder, removeSIP);
      }
		}
		foutput->WriteTObject(hfull[tt]);

		TH2F* hProjKD = (TH2F*) hfull[tt]->Project3D("zy");
		hProjKD->SetTitle(Form("%s Projection",hfull[tt]->GetTitle()));
		hproj[tt] = hProjKD;

		foutput->WriteTObject(hproj[tt]);

		for (int binx = 1; binx <= hfull[tt]->GetNbinsX(); binx++){
			TH2F* hSlice = (TH2F*) hProjKD->Clone( (binx<hfull[tt]->GetNbinsX() ? Form("%s_SidebandRegion%i",hproj[tt]->GetName(),binx) : Form("%s_SignalRegion",hproj[tt]->GetName()) ) );
			hSlice->SetTitle(Form("%s %s GeV Projection",hfull[tt]->GetTitle(),hfull[tt]->GetXaxis()->GetBinLabel(binx)));
			if (normOption == 1){
				hSlice->SetName(Form("%s_Conditional", hSlice->GetName()));
				hSlice->SetTitle(Form("%s (Conditional in %s)", hSlice->GetTitle(), hSlice->GetYaxis()->GetTitle()));
			}
			for (int binz = 0; binz <= hfull[tt]->GetNbinsZ()+1; binz++){
				double integralZ = hfull[tt]->Integral(binx,binx,0,hfull[tt]->GetNbinsY()+1,binz,binz);
        int nConditionalized=0;
        for (int biny = 0; biny <= hfull[tt]->GetNbinsY()+1; biny++){
          double bincontent = hfull[tt]->GetBinContent(binx, biny, binz);
          double binerror = hfull[tt]->GetBinError(binx, biny, binz);
          if (normOption==1 && integralZ!=0){
            bincontent/=integralZ; binerror/=integralZ;
            if (bincontent>0) nConditionalized++;
          }
          hSlice->SetBinContent(biny, binz, bincontent);
          hSlice->SetBinError(biny, binz, 0);
        }
        if (normOption==1 && integralZ!=0 && nConditionalized<=10){
          for (int biny = 0; biny <= hfull[tt]->GetNbinsY()+1; biny++){
            hSlice->SetBinContent(biny, binz, 0);
          }
        }
      }

      if (tt==0) cout << "Integral of data slice: " << hSlice->Integral() << endl;
			foutput->WriteTObject(hSlice);

      if (tt!=ntrees) generic_Histo2DPlotter(foutput, coutput_common, hSlice, erg_tev, luminosity[EnergyIndex], folder, removeSIP);
      else generic_Histo2DPlotter(foutput, coutput_common, hSlice, erg_tev, luminosity[EnergyIndex], folder, removeSIP, tgdata[binx-1]);

			delete hSlice;
		}
	}
  for (int binx=0; binx<nMZZ_ranges; binx++){
    delete tgdata[binx];
    delete[] dataArrayX[binx];
    delete[] dataArrayY[binx];
  }
  for (int tt = 0; tt < ntrees_extra; tt++) delete tc_extra[tt];
	for (int tt = 0; tt < ntrees+1; tt++){
		delete hproj[tt];
		delete hfull[tt];
		if(tt<ntrees) delete tc[tt];
	}

  foutput->Close();
  delete hRatio;
  fSSinput->Close();
}

void compare_KDShapeVariation_Untransformed_perChannel_2D(int folder, int iSyst_1=0, int iSyst_2=1, int SR_noBkg=0, int normOption=0, float normSlice1_ValFirst = 0, float normSlice1_ValLast = 0, int removeSIP=2){
  gROOT->ProcessLine(".x tdrstyle.cc");
  gStyle->SetOptStat(0);
  gStyle->SetPalette(1);
  gStyle->SetPadRightMargin(0.20);
  gROOT->ForceStyle();

  TString OUTPUT_NAME = "LifetimeKD_ShapeSystVars_2D_";
  OUTPUT_NAME.Append(user_folder[folder]);
  if (SR_noBkg==1) OUTPUT_NAME.Append("_NoBkgSROnly_");
  if (removeSIP==1) OUTPUT_NAME.Append("_NoSIP_");
  else if (removeSIP==2) OUTPUT_NAME.Append("_NewSIP_");
  if (normOption==1) OUTPUT_NAME.Append("Sliced_");
  OUTPUT_NAME.Append(Form("%s_", cSystVariable[iSyst_1]));
  OUTPUT_NAME.Append(Form("%s_Comparison.root", cSystVariable[iSyst_2]));
  TString coutput_common = user_dir_hep + "Analysis/Plots/";
  coutput_common.Append(Form("%s_vs_%s/", cSystVariable[iSyst_1], cSystVariable[iSyst_2]));
  TString coutput = coutput_common + OUTPUT_NAME;
  TFile* foutput = new TFile(coutput, "recreate");
  const int ntrees = 2;
  char* chProj[ntrees]={ "hdata_full_zy", "hmc_full_zy" };
  TH2F* hProj[nMZZ_ranges*ntrees][3] ={ { 0 } };
  TH2F* hhiggs[kGGHSamples][3] ={ { 0 } };
  TGraph* tgdata[nMZZ_ranges*3][2];
  TFile* finput[3][2] ={ { 0 } };
  double max_plot[nMZZ_ranges] ={ 0 };
  bool received_ergtev[2] ={ 0 };
  int normSlice1Bin[2] ={ 0 };

  for (int erg_tev = 7; erg_tev < 9; erg_tev++){
    TString erg_dir;
    erg_dir.Form("LHC_%iTeV", erg_tev);
    TString comstring;
    comstring.Form("%iTeV", erg_tev);
    int EnergyIndex = 0; if (erg_tev == 8) EnergyIndex = 1;
    TString cinput_common = user_dir_hep + "Analysis/Plots/";
    cinput_common.Append(Form("%s_vs_%s/", cSystVariable[iSyst_1], cSystVariable[iSyst_2]));

    TString INPUT_NAME = "LifetimeKD_SystVariables_DataQQBZZ_2D_";
    if (removeSIP==1) INPUT_NAME.Append("noSIP_");
    if (removeSIP==2) INPUT_NAME.Append("newSIP_");
    INPUT_NAME.Append(user_folder[folder]);
    INPUT_NAME.Append(Form("_%s_", cSystVariable[iSyst_1]));
    INPUT_NAME.Append(Form("_%s_", cSystVariable[iSyst_2]));
    INPUT_NAME.Append(comstring);
    INPUT_NAME.Append(".root");
    cout << INPUT_NAME << endl;
    TString cinput = cinput_common + INPUT_NAME;
    cout << "Attempting to open " << cinput << endl;
    finput[folder][EnergyIndex] = new TFile(cinput, "read");

    if (finput[folder][EnergyIndex] == 0 || finput[folder][EnergyIndex]->IsZombie()) continue;
    else received_ergtev[EnergyIndex]=true;
    for (int sb = 0; sb < nMZZ_ranges; sb++){
      for (int hh = 0; hh < ntrees; hh++){
        TString hname = Form("%s_SidebandRegion%i", chProj[hh], sb + 1);
        if (sb==nMZZ_ranges-1) hname = Form("%s_SignalRegion", chProj[hh]);
        cout << hname << endl;
        TH2F* htemp = (TH2F*)finput[folder][EnergyIndex]->Get(hname);
        if (htemp==0) cout << "Warning! " << hname << " DNE!" << endl;
        if (hProj[nMZZ_ranges * hh + sb][folder] == 0){
          hProj[nMZZ_ranges*hh+sb][folder] = (TH2F*)htemp->Clone(Form("Summed_%s", htemp->GetName()));
        }
        else hProj[nMZZ_ranges*hh+sb][folder]->Add(htemp);
        cout << htemp->Integral() << endl;
        delete htemp;
      }
      TString tgname = Form("%s_SidebandRegion%i", "tgdata", sb + 1);
      if (sb==nMZZ_ranges-1) tgname = Form("%s_SignalRegion", "tgdata");
      cout << tgname << endl;
      TGraph* tgtemp = (TGraph*)finput[folder][EnergyIndex]->Get(tgname);
      if(tgtemp!=0) tgdata[nMZZ_ranges * folder + sb][EnergyIndex] = tgtemp;
      tgname.Append("_");
      tgname.Append(comstring);
      if (tgdata[nMZZ_ranges * folder + sb][EnergyIndex]!=0) tgdata[nMZZ_ranges * folder + sb][EnergyIndex]->SetName(tgname);
    }
    for (int hh = 0; hh < kGGHSamples; hh++){
      TString hname = Form("hhiggs_CTau%.0f_full", sample_CTau[hh]);
      TH2F* htemp = (TH2F*)finput[folder][EnergyIndex]->Get(hname);
      htemp->SetLineColor(kBlack);
      htemp->SetLineWidth(2);
      htemp->SetTitle("");
      if (hh==0) htemp->SetLineStyle(7);
      else if (hh==1) htemp->SetLineStyle(6);
      else if (hh==2) htemp->SetLineStyle(2);
      else if (hh==3) htemp->SetLineStyle(3);
      if (hhiggs[hh][folder] == 0 && htemp!=0){
        hhiggs[hh][folder] = (TH2F*)htemp->Clone(Form("Summed_%s", htemp->GetName()));
        hhiggs[hh][folder]->SetLineColor(kBlack);
      }
      else if (htemp!=0) hhiggs[hh][folder]->Add(htemp);
      if (htemp!=0) delete htemp;
    }
  }
  cout << "Accumulated the histos" << endl;

  if (normOption == 1){
    normSlice1Bin[0] = hProj[0][folder]->GetYaxis()->FindBin(normSlice1_ValFirst);
    normSlice1Bin[1] = hProj[0][folder]->GetYaxis()->FindBin(normSlice1_ValLast);
  }
  cout << "Slice boundaries: " << normSlice1Bin[0] << " " << normSlice1Bin[1] << endl;

  TString canvasDir = coutput_common;
  foutput->cd();

  for (int sb = 0; sb < nMZZ_ranges; sb++){
    cout << "Rebuilding region " << sb  << " tgdata" << endl;
    double* dataX[2];
    double* dataY[2];
    int dataN[2] ={ 0 };
    TString tgname;
    for (int ee=0; ee<2; ee++){
      if (tgdata[nMZZ_ranges * folder + sb][ee]!=0){
        dataX[ee] = tgdata[nMZZ_ranges * folder + sb][ee]->GetX();
        dataY[ee] = tgdata[nMZZ_ranges * folder + sb][ee]->GetY();
        dataN[ee] = tgdata[nMZZ_ranges * folder + sb][ee]->GetN();
        cout << "Ndata[" << ee << "]: " << dataN[ee] << endl;
        tgname = tgdata[nMZZ_ranges * folder + sb][ee]->GetName();
        cout << tgname << endl;
      }
    }
    int dataN_new = dataN[0] + dataN[1];
    if (dataN_new!=0){
      cout << "New TGraph nbins: " << dataN_new <<  endl;
      cout << "Compare to hdata: " << hProj[nMZZ_ranges * 0 + sb][folder]->Integral() << endl;
      double* dataX_new = new double[dataN_new];
      double* dataY_new = new double[dataN_new];
      for (int ee=0; ee<2; ee++){
        for (int bin=0; bin<dataN[ee]; bin++){
          int newbin = bin;
          if (ee==1) newbin += dataN[0];
          dataX_new[newbin] = dataX[ee][bin];
          dataY_new[newbin] = dataY[ee][bin];
        }
      }
      for (int ee=0; ee<2; ee++){
        dataX[ee] = 0;
        dataY[ee] = 0;
        delete tgdata[nMZZ_ranges * folder + sb][ee];
        tgdata[nMZZ_ranges * folder + sb][ee] = 0;
      }
      tgdata[nMZZ_ranges * folder + sb][0] = new TGraph(dataN_new, dataX_new, dataY_new);
      tgdata[nMZZ_ranges * folder + sb][0]->SetName(tgname);
      tgdata[nMZZ_ranges * folder + sb][0]->SetTitle("");
      delete[] dataX_new;
      delete[] dataY_new;
    }
  }

  if (normOption==0){
    for (int sb = 0; sb < nMZZ_ranges+1; sb++){
      foutput->cd();
      gStyle->SetOptStat(0);
      if (SR_noBkg == 1 && sb < nMZZ_ranges - 1) continue;

      TString strSidebandRegion = Form("Sideband %i", sb+1);
      if (sb==nMZZ_ranges-1) strSidebandRegion = "Signal Region";
      else if (sb==nMZZ_ranges) strSidebandRegion = "Sidebands 1&3";
      float ptRegion_xhigh = 0.75;
      if (sb==nMZZ_ranges-1) ptRegion_xhigh = 0.76;
      if (sb==nMZZ_ranges) ptRegion_xhigh = 0.77;
      TPaveText* ptRegion = new TPaveText(0.61, 0.86, ptRegion_xhigh, 0.90, "brNDC");
      ptRegion->SetBorderSize(0);
      ptRegion->SetTextAlign(12);
      ptRegion->SetTextSize(0.03);
      ptRegion->SetFillStyle(1001);
      ptRegion->SetFillColor(0);
      ptRegion->SetTextFont(42);
      TText* textRegion = ptRegion->AddText(0.0, 0.0, strSidebandRegion);

      int sbIndex = sb;
      if (sb==nMZZ_ranges){
        sbIndex=0;
        int sb3Index = nMZZ_ranges-2;
        for (int hh = 0; hh < ntrees; hh++){
          hProj[nMZZ_ranges * hh + sbIndex][folder]->Add(hProj[nMZZ_ranges * hh + sb3Index][folder], 1);
          hProj[nMZZ_ranges * hh + sbIndex][folder]->SetName(Form("%s%s", hProj[nMZZ_ranges * hh + sbIndex][folder]->GetName(),"3"));
        }
        double* dataX[2];
        double* dataY[2];
        int dataN[2] ={ 0 };
        int tgsbInd[2] = {sbIndex, sb3Index};
        TString tgname;
        for (int ee=0; ee<2; ee++){
          if (tgdata[nMZZ_ranges * folder + tgsbInd[ee]][0]!=0){
            dataX[ee] = tgdata[nMZZ_ranges * folder + tgsbInd[ee]][0]->GetX();
            dataY[ee] = tgdata[nMZZ_ranges * folder + tgsbInd[ee]][0]->GetY();
            dataN[ee] = tgdata[nMZZ_ranges * folder + tgsbInd[ee]][0]->GetN();
            cout << "Ndata[" << ee << "]: " << dataN[ee] << endl;
            tgname = tgdata[nMZZ_ranges * folder + tgsbInd[ee]][0]->GetName();
            cout << tgname << endl;
          }
        }
        int dataN_new = dataN[0] + dataN[1];
        if (dataN_new!=0){
          double* dataX_new = new double[dataN_new];
          double* dataY_new = new double[dataN_new];
          for (int ee=0; ee<2; ee++){
            for (int bin=0; bin<dataN[ee]; bin++){
              int newbin = bin;
              if (ee==1) newbin += dataN[0];
              dataX_new[newbin] = dataX[ee][bin];
              dataY_new[newbin] = dataY[ee][bin];
            }
          }
          for (int ee=0; ee<2; ee++){
            dataX[ee] = 0;
            dataY[ee] = 0;
            delete tgdata[nMZZ_ranges * folder + tgsbInd[ee]][0];
            tgdata[nMZZ_ranges * folder + tgsbInd[ee]][0] = 0;
          }
          tgdata[nMZZ_ranges * folder + tgsbInd[0]][0] = new TGraph(dataN_new, dataX_new, dataY_new);
          tgdata[nMZZ_ranges * folder + tgsbInd[0]][0]->SetName(tgname);
          tgdata[nMZZ_ranges * folder + tgsbInd[0]][0]->SetTitle("");
          delete[] dataX_new;
          delete[] dataY_new;
        }
      }

      for (int hh = 0; hh < ntrees; hh++){
        floorHistogram(hProj[nMZZ_ranges * hh + sbIndex][folder]);
        generic_Histo2DPlotter(foutput, canvasDir, hProj[nMZZ_ranges * hh + sbIndex][folder], -1, 0, folder, removeSIP, tgdata[nMZZ_ranges * folder + sbIndex][0], ptRegion);
      }
      delete ptRegion;
      /*
            if (sb == nMZZ_ranges - 1 && SR_noBkg == 1){
            hhiggs[0][folder]->GetXaxis()->SetTitle(xTitle);
            hhiggs[0][folder]->GetYaxis()->SetTitle(yTitle);
            hhiggs[0][folder]->GetYaxis()->SetRangeUser(0, max_plot[sb] * 1.25);
            if (iSyst == 3) hhiggs[0][folder]->GetYaxis()->SetRangeUser(0, max_plot[sb] * 1.85);
            if (iSyst == 2) hhiggs[0][folder]->GetYaxis()->SetRangeUser(0, max_plot[sb] * 1.75);
            hhiggs[0][folder]->GetXaxis()->SetNdivisions(505);
            hhiggs[0][folder]->GetXaxis()->SetLabelFont(42);
            hhiggs[0][folder]->GetXaxis()->SetLabelOffset(0.007);
            hhiggs[0][folder]->GetXaxis()->SetLabelSize(0.04);
            hhiggs[0][folder]->GetXaxis()->SetTitleSize(0.06);
            hhiggs[0][folder]->GetXaxis()->SetTitleOffset(0.9);
            hhiggs[0][folder]->GetXaxis()->SetTitleFont(42);
            hhiggs[0][folder]->GetYaxis()->SetNdivisions(505);
            hhiggs[0][folder]->GetYaxis()->SetLabelFont(42);
            hhiggs[0][folder]->GetYaxis()->SetLabelOffset(0.007);
            hhiggs[0][folder]->GetYaxis()->SetLabelSize(0.04);
            hhiggs[0][folder]->GetYaxis()->SetTitleSize(0.06);
            hhiggs[0][folder]->GetYaxis()->SetTitleOffset(1.1);
            hhiggs[0][folder]->GetYaxis()->SetTitleFont(42);
            hhiggs[0][folder]->GetXaxis()->SetRangeUser(visualRange_SystVariable[iSyst][0], visualRange_SystVariable[iSyst][1]);
            }

            string folder_id;
            if (folder==0) folder_id = "4#mu";
            if (folder==1) folder_id = "4e";
            if (folder==2) folder_id = "2e2#mu";
            TH1F* higgsclone[kGGHSamples];
            if (sb < nMZZ_ranges - 1){
            if (tgdata[nMZZ_ranges * folder + sb]!=0) l2D->AddEntry(tgdata[nMZZ_ranges * folder + sb], Form("Observed %s", folder_id.c_str()), "ep");
            }
            else{
            for (int hh = 0; hh < kGGHSamples; hh++){
            if (hhiggs[hh][folder] != 0) higgsclone[hh] = (TH1F*)hhiggs[hh][folder]->Clone(Form("%s_CLONED", hhiggs[hh][folder]->GetName()));
            higgsclone[hh]->SetLineColor(kBlack);
            }
            if (SR_noBkg == 1) higgsclone[0]->SetLineStyle(1);
            l2D->AddEntry(higgsclone[0], "Higgs c#tau_{H}=0 #mum", "l");
            if (SR_noBkg == 1){
            higgsclone[1]->SetLineStyle(1);
            higgsclone[1]->SetLineColor(kBlue);
            l2D->AddEntry(higgsclone[1], "Higgs c#tau_{H}=100 #mum", "l");
            higgsclone[2]->SetLineStyle(1);
            higgsclone[2]->SetLineColor(kGreen+2);
            //        l2D->AddEntry(higgsclone[2], "Higgs c#tau_{H}=500 #mum", "l");
            higgsclone[3]->SetLineStyle(1);
            higgsclone[3]->SetLineColor(kRed);
            l2D->AddEntry(higgsclone[3], "Higgs c#tau_{H}=1000 #mum", "l");
            }
            else if (normOption==0){
            l2D->AddEntry(higgsclone[1], "Higgs c#tau_{H}=100 #mum", "l");
            }
            }
            if (SR_noBkg == 0){
            l2D->AddEntry(hProj[nMZZ_ranges * 1 + sb][folder], Form("Bkg. %s", folder_id.c_str()), "l");
            hProj[nMZZ_ranges * 1 + sb][folder]->Draw("hist");
            l2D->Draw("same");
            if (sb < nMZZ_ranges - 1){
            if (tgdata[nMZZ_ranges * folder + sb]!=0) tgdata[nMZZ_ranges * folder + sb]->Draw("e1psame");
            }
            else{
            cout << "Drawing Higgs MC now" << endl;
            hhiggs[0][folder]->SetLineColor(hProj[nMZZ_ranges * 1 + sb][folder]->GetLineColor());
            hhiggs[0][folder]->Draw("histsame");
            if (normOption==0){
            hhiggs[1][folder]->SetLineColor(hProj[nMZZ_ranges * 1 + sb][folder]->GetLineColor());
            hhiggs[1][folder]->Draw("histsame");
            }
            cout << "Passed draw" << endl;
            }
            }
            else if (sb == nMZZ_ranges - 1){
            cout << "Drawing Higgs MC now" << endl;
            for (int hh = 0; hh < kGGHSamples; hh++){
            if (hh==2)continue;
            if (hhiggs[hh][folder] == 0){
            cout << "Skipping draw for " << hh << "\t" << folder << endl;
            continue;
            }
            if (hh==0) hhiggs[hh][folder]->SetLineColor(kBlack);
            if (hh==1) hhiggs[hh][folder]->SetLineColor(kBlue);
            if (hh==2) hhiggs[hh][folder]->SetLineColor(kGreen+2);
            if (hh==3) hhiggs[hh][folder]->SetLineColor(kRed);
            hhiggs[hh][folder]->SetLineStyle(1);
            if (hh==0) hhiggs[hh][folder]->Draw("hist");
            else  hhiggs[hh][folder]->Draw("histsame");
            }
            l2D->Draw("same");
            }

            pt->Draw();

            */
    }
  }
  else{
    cout << "Recognized normOpt=0" << endl;
    TH1F* hSlice[ntrees][2] ={ { 0 } };
    TH1F* hSlice_higgs[kGGHSamples][2] ={ { 0 } };
    for (int sb = 0; sb < nMZZ_ranges+1; sb++){
      int sbIndex = sb;
      double max_plot_1D = 0;
      if (sb==nMZZ_ranges){
        sbIndex=0;
        int sb3Index = nMZZ_ranges-2;
        for (int hh = 0; hh < ntrees; hh++){
          if (hProj[nMZZ_ranges * hh + sbIndex][folder]==0 || hProj[nMZZ_ranges * hh + sb3Index][folder]==0) cerr << "Sideband 1 or 3 histo for index " << hh << " is 0." << endl;
          else cout << "Processing proj " << hh << endl;
          hProj[nMZZ_ranges * hh + sbIndex][folder]->Add(hProj[nMZZ_ranges * hh + sb3Index][folder], 1);
          hProj[nMZZ_ranges * hh + sbIndex][folder]->SetName(Form("%s%s", hProj[nMZZ_ranges * hh + sbIndex][folder]->GetName(), "3"));
        }
      }
      if (sb==nMZZ_ranges-1){
        for (int hh=0; hh<kGGHSamples; hh++){
//          hhiggs[hh][folder]->Add(hProj[nMZZ_ranges * 1 + sbIndex][folder]);
          if (hhiggs[hh][folder]==0) cerr << "Higgs " << hh << " is 0." << endl;
          else cout << "Processing Higgs " << hh << endl;
          hSlice_higgs[hh][0] = (TH1F*)hhiggs[hh][folder]->ProjectionX(Form("%s_Slice0", hhiggs[hh][folder]->GetName()), normSlice1Bin[0], normSlice1Bin[1]);
          hSlice_higgs[hh][1] = (TH1F*)hhiggs[hh][folder]->ProjectionX(Form("%s_Slice1", hhiggs[hh][folder]->GetName()));
          hSlice_higgs[hh][1]->SetLineStyle(7);
          for (int ss=0; ss<2; ss++){
            hSlice_higgs[hh][ss]->Scale(1./hSlice_higgs[hh][ss]->Integral(0, hSlice_higgs[hh][ss]->GetNbinsX()+1));
            hSlice_higgs[hh][ss]->SetLineWidth(2);
          }
        }
        for (int ss=0; ss<2; ss++){
          hSlice_higgs[0][ss]->SetLineColor(kBlack);
          hSlice_higgs[1][ss]->SetLineColor(kViolet);
          hSlice_higgs[2][ss]->SetLineColor(kViolet);
          hSlice_higgs[3][ss]->SetLineColor(kViolet);
          max_plot_1D = max(max_plot_1D, max(hSlice_higgs[0][ss]->GetMaximum(), hSlice_higgs[3][ss]->GetMaximum()));
        }
      }
      for (int hh = 0; hh < ntrees; hh++){
        if (hProj[nMZZ_ranges * hh + sbIndex][folder]==0) cerr << "Proj histo " << hh << " is 0." << endl;
        hSlice[hh][0] = (TH1F*)hProj[nMZZ_ranges * hh + sbIndex][folder]->ProjectionX(Form("%s_Slice0", hProj[nMZZ_ranges * hh + sbIndex][folder]->GetName()), normSlice1Bin[0], normSlice1Bin[1]);
        hSlice[hh][1] = (TH1F*)hProj[nMZZ_ranges * hh + sbIndex][folder]->ProjectionX(Form("%s_Slice1", hProj[nMZZ_ranges * hh + sbIndex][folder]->GetName()));
        hSlice[hh][1]->Add(hSlice[hh][0], -1.);
        for (int ss=0; ss<2; ss++) hSlice[hh][ss]->Scale(1./hSlice[hh][ss]->Integral(0, hSlice[hh][ss]->GetNbinsX()+1));

        hSlice[hh][0]->SetLineColor(kRed);
        hSlice[hh][1]->SetLineColor(kRed);
        hSlice[hh][0]->SetLineWidth(2);
        hSlice[hh][1]->SetLineWidth(2);
        hSlice[hh][1]->SetLineStyle(7);

        max_plot_1D = max(max_plot_1D, max(hSlice[hh][0]->GetMaximum(), hSlice[hh][1]->GetMaximum()));
      }

      cout << "Max plot: " << max_plot_1D << endl;

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
      if (received_ergtev[0] && received_ergtev[1]) text = pt->AddText(0.537, 0.45, "#font[42]{19.7 fb^{-1} (8 TeV) + 5.1 fb^{-1} (7 TeV)}");
      else if (received_ergtev[1]) text = pt->AddText(0.837, 0.45, "#font[42]{19.7 fb^{-1} (8 TeV)}");
      else if (received_ergtev[0]) text = pt->AddText(0.837, 0.45, "#font[42]{ 5.1 fb^{-1} (7 TeV)}");
      text->SetTextSize(0.0315);

      TString canvasname_2D = "cCompare_VtxSystVar_Slices_";
      canvasname_2D.Append(user_folder[folder]);
      canvasname_2D.Append("_");
      if (removeSIP == 1) canvasname_2D.Append("NoSIP_");
      if (removeSIP == 2) canvasname_2D.Append("NewSIP_");
      if (SR_noBkg == 1) canvasname_2D.Append("NoBkg_");
      if (sb < nMZZ_ranges - 1) canvasname_2D.Append(Form("%s_SB%i", cSystVariable[iSyst_1], sb + 1));
      else if (sb == nMZZ_ranges - 1) canvasname_2D.Append(Form("%s_SR", cSystVariable[iSyst_1]));
      else if (sb == nMZZ_ranges) canvasname_2D.Append(Form("%s_SB13", cSystVariable[iSyst_1]));
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

      TLegend *l2D = new TLegend(0.20, 0.57, 0.58, 0.90);
      l2D->SetBorderSize(0);
      l2D->SetTextFont(42);
      l2D->SetTextSize(0.03);
      l2D->SetLineColor(1);
      l2D->SetLineStyle(1);
      l2D->SetLineWidth(1);
      l2D->SetFillColor(0);
      l2D->SetFillStyle(0);

      TString yTitle = Form("Rate / %.2f", hSlice[1][0]->GetXaxis()->GetBinWidth(1));
      TString xTitle = cSystVariable_label[iSyst_1];


      for (int hh=0; hh<ntrees; hh++){
        for (int ss=0; ss<2; ss++){
          hSlice[hh][ss]->GetXaxis()->SetTitle(xTitle);
          hSlice[hh][ss]->GetYaxis()->SetTitle(yTitle);
          hSlice[hh][ss]->GetYaxis()->SetRangeUser(0, max_plot_1D * 1.25);
          hSlice[hh][ss]->GetXaxis()->SetNdivisions(505);
          hSlice[hh][ss]->GetXaxis()->SetLabelFont(42);
          hSlice[hh][ss]->GetXaxis()->SetLabelOffset(0.007);
          hSlice[hh][ss]->GetXaxis()->SetLabelSize(0.04);
          hSlice[hh][ss]->GetXaxis()->SetTitleSize(0.06);
          hSlice[hh][ss]->GetXaxis()->SetTitleOffset(0.9);
          hSlice[hh][ss]->GetXaxis()->SetTitleFont(42);
          hSlice[hh][ss]->GetYaxis()->SetNdivisions(505);
          hSlice[hh][ss]->GetYaxis()->SetLabelFont(42);
          hSlice[hh][ss]->GetYaxis()->SetLabelOffset(0.007);
          hSlice[hh][ss]->GetYaxis()->SetLabelSize(0.04);
          hSlice[hh][ss]->GetYaxis()->SetTitleSize(0.06);
          hSlice[hh][ss]->GetYaxis()->SetTitleOffset(1.1);
          hSlice[hh][ss]->GetYaxis()->SetTitleFont(42);
          hSlice[hh][ss]->GetXaxis()->SetRangeUser(visualRange_SystVariable[iSyst_1][0], visualRange_SystVariable[iSyst_1][1]);
        }
      }

      string folder_id;
      if (folder==0) folder_id = "4#mu";
      if (folder==1) folder_id = "4e";
      if (folder==2) folder_id = "2e2#mu";
      float sliceLow = hProj[nMZZ_ranges * 0 + sbIndex][folder]->GetYaxis()->GetBinLowEdge(normSlice1Bin[0]);
      float sliceHigh = hProj[nMZZ_ranges * 0 + sbIndex][folder]->GetYaxis()->GetBinUpEdge(normSlice1Bin[1]);
      TString strLabel_Slice1 = Form(" ([%.1f, %.1f)", sliceLow, sliceHigh);
      if (normSlice1Bin[0]==0) strLabel_Slice1 = Form(" (<%.1f", sliceHigh);
      if (normSlice1Bin[1]==hSlice[1][0]->GetNbinsX()+1) strLabel_Slice1 = Form(" (>%.1f", sliceLow);
      if (strstr(cSystVariable_label[iSyst_2], "#mum")!=0) strLabel_Slice1.Append(" #mum");
      TString strLabel_Slice2 = " (else";
      if (normSlice1Bin[0]==0) strLabel_Slice2 = Form(" (>%.1f", sliceHigh);
      if (normSlice1Bin[1]==hSlice[1][0]->GetNbinsX()+1) strLabel_Slice2 = Form(" (<%.1f", sliceLow);
      if ((normSlice1Bin[0]==0 || normSlice1Bin[1]==hSlice[1][0]->GetNbinsX()+1) && strstr(cSystVariable_label[iSyst_2], "#mum")!=0) strLabel_Slice2.Append(" #mum");
      strLabel_Slice1.Append(")");
      strLabel_Slice2.Append(")");

      l2D->AddEntry(hSlice[1][0], Form("Bkg. %s", strLabel_Slice1.Data()), "l");
      l2D->AddEntry(hSlice[1][1], Form("Bkg. %s", strLabel_Slice2.Data()), "l");
      hSlice[1][0]->Draw("hist");
      hSlice[1][1]->Draw("histsame");
      if (sb==nMZZ_ranges-1){
        l2D->AddEntry(hSlice_higgs[0][0], Form("c#tau_{H}=0 #mum %s", strLabel_Slice1.Data()), "l");
        l2D->AddEntry(hSlice_higgs[0][1], Form("c#tau_{H}=0 #mum %s", strLabel_Slice2.Data()), "l");
        l2D->AddEntry(hSlice_higgs[3][0], Form("c#tau_{H}=1000 #mum %s", strLabel_Slice1.Data()), "l");
        l2D->AddEntry(hSlice_higgs[3][1], Form("c#tau_{H}=1000 #mum %s", strLabel_Slice2.Data()), "l");
        for (int ss=0; ss<2; ss++){
          hSlice_higgs[0][ss]->Draw("histsame");
          hSlice_higgs[3][ss]->Draw("histsame");
        }
      }
      l2D->Draw("same");

      pt->Draw();

      TString strSidebandRegion = Form("Sideband %i", sb+1);
      if (sb==nMZZ_ranges-1) strSidebandRegion = "Signal Region";
      else if (sb==nMZZ_ranges) strSidebandRegion = "Sidebands 1&3";

      TPaveText *pt10 = new TPaveText(0.67, 0.84, 0.90, 0.92, "brNDC");
      pt10->SetBorderSize(0);
      pt10->SetTextAlign(12);
      pt10->SetTextSize(0.03);
      pt10->SetFillStyle(0);
      pt10->SetTextFont(42);
      TText* text10 = pt10->AddText(0.01, 0.01, strSidebandRegion);
      pt10->Draw();

      TString strSIPcut = "";
      if (removeSIP == 1) strSIPcut = "No SIP";
      else if (removeSIP == 2) strSIPcut = "New cut";
      strSIPcut.Prepend(", ");
      strSIPcut.Prepend(folder_id.c_str());
      TPaveText* pt20;
      pt20 = new TPaveText(0.67, 0.80, 0.90, 0.84, "brNDC");
      pt20->SetBorderSize(0);
      pt20->SetTextAlign(12);
      pt20->SetTextSize(0.03);
      pt20->SetFillStyle(0);
      pt20->SetTextFont(42);
      TText* text20 = pt20->AddText(0.01, 0.01, strSIPcut);
      if (removeSIP>0) pt20->Draw();


      c2D->RedrawAxis();
      c2D->Modified();
      c2D->Update();
      foutput->WriteTObject(c2D);

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


      for (int hh = 0; hh < ntrees; hh++){
        for (int ss=0; ss<2; ss++) delete hSlice[hh][ss];
      }
      if (sb==nMZZ_ranges-1){
        for (int hh=0; hh<kGGHSamples; hh++){
          for (int ss=0; ss<2; ss++) delete hSlice_higgs[hh][ss];
        }
      }
    }
  }
  for (int sb = 0; sb < nMZZ_ranges; sb++){
    for (int hh = 0; hh < ntrees; hh++) delete hProj[nMZZ_ranges * hh + sb][folder];
    for (int ee = 0; ee < 2; ee++){
      if (tgdata[nMZZ_ranges * folder + sb][ee]!=0) delete tgdata[nMZZ_ranges * folder + sb][ee];
    }
  }

  for (int e = 0; e < 2; e++){
    if (finput[folder][e] == 0 || finput[folder][e]->IsZombie()) continue;
    else finput[folder][e]->Close();
  }
  foutput->Close();
}


void plotSystematics_2D(int iSyst_1=0, int iSyst_2=1, int normOption=0, int removeSIP=2){
  for (int e=7; e<=8; e++){
    for (int f=0; f<3; f++) produce_KDSystVariables_2D(f, e, iSyst_1, iSyst_2, normOption, removeSIP);
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

