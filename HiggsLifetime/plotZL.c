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


const int nSystVars=25;
char* cSystVariable[nSystVars] = { "Txy", "Dxy", "delTxy", "delDxy", "pToverMZZ", "pullTxy","pullDxy",

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

"KalmanCandVtx_chi2",
"Lep12_Z1SIP",
"Lep3_Z1SIP",
"PVSIP"
};

char* cSystVariable_label[nSystVars] = { "T_{xy} (#mum)", "D_{xy} (#mum)", "#deltaT_{xy} (#mum)", "#deltaD_{xy} (#mum)", "p^{T}_{3l} / m_{3l}",
"(T^{reco}_{xy} - T^{true}_{xy}) / #deltaT_{xy}",
"(D^{reco}_{xy} - D^{true}_{xy}) / #deltaD_{xy}",

"(r^{reco}_{PV,T} - r^{true}_{PV,T}) / #sigma_{PV,T} . #hat{p}^{true}_{T}",
"(r^{reco}_{3l,T} - r^{true}_{3l,T}) / #sigma_{3l,T} . #hat{p}^{true}_{T}",
"(x^{reco}_{3l} - x^{true}_{3l}) / #sigma_{3l,x}",
"(y^{reco}_{3l} - y^{true}_{3l}) / #sigma_{3l,y}",
"(z^{reco}_{3l} - z^{true}_{3l}) / #sigma_{3l,z}",

"(r^{reco}_{PV} - r^{true}_{PV}) . #hat{p}^{true}_{T} (#mum)",
"(r^{reco}_{3l} - r^{true}_{3l}) . #hat{p}^{true}_{T} (#mum)",

"x^{reco}_{3l} - x^{true}_{3l} (#mum)",
"y^{reco}_{3l} - y^{true}_{3l} (#mum)",
"z^{reco}_{3l} - z^{true}_{3l} (#mum)",

"x^{reco}_{OPV} - x^{true}_{OPV} (#mum)",
"y^{reco}_{OPV} - y^{true}_{OPV} (#mum)",
"z^{reco}_{OPV} - z^{true}_{OPV} (#mum)",

"tanh(T_{xy} / 800 #mum)",

"#chi^{2}_{3l}",
"#it{SIP}^{Z1}_{1,2} (Signed)",
"#it{SIP}^{Z1}_{3} (Signed)",
"#it{SIP}_{3D}"
};
int nbins_SystVariable[nSystVars] = {
60, 40, 150, 200, 100,
60, 60,

60,
60,
48,
48,
30,

60,
60,

60,
60,
60,

60,
60,
60,

50,

41,
52,
52,
31
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

	{ -1200, 1200 },
	{ -1200, 1200 },

  { -90, 90 },
  { -90, 90 },
  { -90, 90 },
  { -90, 90 },
  { -90, 90 },
  { -90, 90 },

  { -1, 1 },

  { 0, 41 },
  { -5.2, 5.2 },
  { -5.2, 5.2 },
  { 0, 31 }
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

	{ -1200, 1200 },
	{ -1200, 1200 },

  { -90, 90 },
  { -90, 90 },
  { -90, 90 },
  { -90, 90 },
  { -90, 90 },
  { -90, 90 },

  { -1, 1 },

  { 0, 41 },
  { -5.2, 5.2 },
  { -5.2, 5.2 },
  { 0, 31 }
};

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
  "#it{SIP}_{3D}<4",
  "No cut",
  "#it{SIP}^{Z1} cut",
  "#chi^{2}<30",
  "New cut",
  "New + old"
};
TString strPFMET_label[2]={
  "E^{miss}_{T}<25 GeV",
  "E^{miss}_{T}#geq25 GeV"
};
TString strZ1Category_label[2]={
  "Z->#mu#mu",
  "Z->ee"
};
TString strLep3Category_label[2]={
  " + #mu",
  " + e"
};
TString strLepIsoCut_label[3] ={
  "Tight",
  "2p2f",
  "3p1f"
};
TString Z1masswidth_label[3]={
  "40<m_{1}<120",
  "|m_{1}-m_{Z}|<10",
  "|m_{3l}-m_{Z}|<5"
};
const int nProcesses = 7;
char* processName[nProcesses] ={
  "qqZZ", "DY_NoB", "DY_B", "ttbar", "WZ", "WW", "CRZL"
};



double applyCRselection(
  int option[6],

  float Lep1combRelIsoPF, float Lep2combRelIsoPF, float Lep3combRelIsoPF,
  bool Lep1isID, bool Lep2isID, bool Lep3isID,
  short Z1ids,
  float PFMET,

  float Lep_Z1SIP[3],
  float KalmanCandVtx_chi2,
  float LepSIP[3],
  int LepID[3],

  float Z1Mass, float ZZMass
  ){

  int iCR = option[0]; // PFMet<25; <=25 
  int icut = option[1]; // SIP selection
  int catZ1 = option[2]; // Z->mumu; Z->ee
  int catZ2 = option[3]; // +mu, +e
  int isocut = option[4]; // Z tight l loose, Z+l tight
  int iM1cut = option[5]; // No M1 cut, |M1-mZ|<10; |m3l-mZ|<5

  int Lep1ID = LepID[0];
  int Lep2ID = LepID[1];
  int Lep3ID = LepID[2];

  float Lep1_Z1SIP = Lep_Z1SIP[0];
  float Lep2_Z1SIP = Lep_Z1SIP[1];
  float Lep3_Z1SIP = Lep_Z1SIP[2];

  float Lep1SIP = LepSIP[0];
  float Lep2SIP = LepSIP[1];
  float Lep3SIP = LepSIP[2];

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
  bool strLep3Category[2]={
    fabs(Lep3ID) == 13,
    fabs(Lep3ID) == 11
  };
  bool strLepIsoCut = Lep3combRelIsoPF<0.4 && Lep3isID;  // Lep1,2 already have the cut; what is seen on ntuples is the relIso before FSR

  bool strcut = strnewSIP_CRZL[icut];
  strcut = strcut && strZ1Category[catZ1];
  strcut = strcut && strLep3Category[catZ2];

  if (isocut==0){
    strcut = strcut;
  }
  else if (isocut==1){
    strcut = strcut && strLepIsoCut;
  }
  if (iM1cut==1){
    strcut = strcut && fabs(Z1Mass-Z1pole)<Z1masswidth[iM1cut-1];
  }
  else if (iM1cut==2){
    strcut = strcut && fabs(ZZMass-Z1pole)<Z1masswidth[iM1cut-1];
  }
  if (iCR==1){
    strcut = strcut && PFMET>=25;
  }
  else{
    strcut = strcut && PFMET<25;
  }

  double result = (strcut ? 1 : 0);
  return result;
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

void plot1DNewOld(int iSyst, int iCR, int isocut, int iM1cut, int catZ1, int catZ2, int icut, int erg_tev, TH1F** hfull, TFile* foutput, bool isNorm);

TString produceCRname(int iCR, int isocut, int iM1cut, int catZ1, int catZ2, int icut/*, char* cappend*/){
  TString name = "ZL_";
  if (iCR==0) name.Append("METNonPrompt");
  else name.Append("METPrompt");

  if (catZ1==0) name.Append("_Zmumu");
  else name.Append("_Zee");

  if (catZ2==0) name.Append("_mu");
  if (catZ2==1) name.Append("_e");

  name.Append("_");
  name += cutNames[icut];

  if (isocut==0) name.Append("_AllLoose");
  else{ name.Append("_"); name += strLepIsoCut_label[isocut-1]; }

  if (iM1cut==1) name.Append("_ZPeak");
  else if (iM1cut==2) name.Append("_ZGPeak");
/*
  name.Append("_");
  name.Append(cappend);
*/
  return name;
}


void produce_KDSystVariables(
  int erg_tev,
  int iSyst=0
){

	char TREE_NAME[]="SelectedTree";
	TString erg_dir;
	erg_dir.Form("LHC_%iTeV",erg_tev);
	TString comstring;
	comstring.Form("%iTeV",erg_tev);
	int EnergyIndex=0;if(erg_tev==8) EnergyIndex=1;

	TString OUTPUT_NAME = "LifetimeKD_SystVariables_DataMC_ZL";
	OUTPUT_NAME.Append(Form("_%s_",cSystVariable[iSyst]));
	OUTPUT_NAME.Append(comstring);
	OUTPUT_NAME.Append(".root");

  float MC_weight;
  float MC_weight_noxsec;
	float PFMET,Z1Mass,ZZMass,ZZPt,ZZEta,ZZPhi;
	float KalmanCandVtx_x, KalmanCandVtx_y, KalmanCandVtx_z, KalmanCandVtx_chi2;
	float OfflinePrimaryVtx_x,OfflinePrimaryVtx_y,OfflinePrimaryVtx_z;
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

	float Lep1_Z1SIP, Lep2_Z1SIP, Lep3_Z1SIP;
  float Lep1SIP, Lep2SIP, Lep3SIP;
  int Lep1ID, Lep2ID, Lep3ID;

	float Txy,Dxy,Txy_true,Dxy_true,delTxy,delDxy,KD;
	float sigmaPV_xy,sigmaInt_xy;
	float CandVtx_x,CandVtx_y,CandVtx_z;
	float TrueCandVtx_x,TrueCandVtx_y,TrueCandVtx_z;

  float Lep1combRelIsoPF=0, Lep2combRelIsoPF=0, Lep3combRelIsoPF=0;
  bool Lep1isID=true, Lep2isID=true, Lep3isID=true;
  short Z1ids;

	TString cinput_common_noSIP = user_dir_hep + "No_SIP/" + erg_dir + "/CR/";
	TString coutput_common = user_dir_hep + "Analysis/Plots/ZL_MC/";
	coutput_common.Append(Form("%s/",cSystVariable[iSyst]));
	TString mkdirCommand = "mkdir -p ";
	mkdirCommand.Append(coutput_common);
	gSystem->Exec(mkdirCommand);
	TString coutput = coutput_common + OUTPUT_NAME;
	TFile* foutput = new TFile(coutput,"recreate");

  const int ntrees=nProcesses;

  TString cinput;
  TChain* tc[ntrees];
  for(int tt=0;tt<ntrees;tt++) tc[tt] = new TChain(TREE_NAME);
	for (int smp = kGGSamples; smp < kQQBZZSamples; smp++){
		cinput = cinput_common_noSIP + sample_FullSim[smp] + "_CRZLTree.root";
		tc[0]->Add(cinput);
	}
  for (int smp = 0; smp < 2; smp++){
    cinput = cinput_common_noSIP + sample_CR_MC[smp] + "_CRZLTree.root";
    tc[1]->Add(cinput);
  }
  for (int smp = 2; smp < 4; smp++){
    cinput = cinput_common_noSIP + sample_CR_MC[smp] + "_CRZLTree.root";
    tc[2]->Add(cinput);
  }
  for (int smp = 4; smp < nCRZLLMC; smp++){
    cinput = cinput_common_noSIP + sample_CR_MC[smp] + "_CRZLTree.root";
    tc[smp-1]->Add(cinput);
  }
  cinput = cinput_common_noSIP + data_files[2] + "_CRZLTree.root";
  tc[6]->Add(cinput);

  for (int tt = 0; tt < ntrees; tt++){
    cout << "Tree: " << tt << ", nEvents: " << tc[tt]->GetEntries() << endl;
    if (tt != ntrees-1){
      tc[tt]->SetBranchAddress("MC_weight", &MC_weight);
      tc[tt]->SetBranchAddress("MC_weight_noxsec", &MC_weight_noxsec);
    }
    tc[tt]->SetBranchAddress("Lep1combRelIsoPF", &Lep1combRelIsoPF);
    tc[tt]->SetBranchAddress("Lep2combRelIsoPF", &Lep2combRelIsoPF);
    tc[tt]->SetBranchAddress("Lep3combRelIsoPF", &Lep3combRelIsoPF);

    tc[tt]->SetBranchAddress("Z1ids", &Z1ids);
    tc[tt]->SetBranchAddress("PFMET", &PFMET);
    tc[tt]->SetBranchAddress("Z1Mass", &Z1Mass);
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
    if (tc[tt]->GetBranchStatus("KalmanCandVtx_chi2")){
      if (tc[tt]->GetBranchStatus("Lep1_Z1SIP")){
        tc[tt]->SetBranchAddress("Lep1_Z1SIP", &Lep1_Z1SIP);
        tc[tt]->SetBranchAddress("Lep2_Z1SIP", &Lep2_Z1SIP);
        tc[tt]->SetBranchAddress("Lep3_Z1SIP", &Lep3_Z1SIP);
      }
      tc[tt]->SetBranchAddress("KalmanCandVtx_chi2", &KalmanCandVtx_chi2);
    }
    if (tc[tt]->GetBranchStatus("Lep1SIP")){
      tc[tt]->SetBranchAddress("Lep1SIP", &Lep1SIP);
      tc[tt]->SetBranchAddress("Lep2SIP", &Lep2SIP);
      tc[tt]->SetBranchAddress("Lep3SIP", &Lep3SIP);
    }
    tc[tt]->SetBranchAddress("Lep1ID", &Lep1ID);
    tc[tt]->SetBranchAddress("Lep2ID", &Lep2ID);
    tc[tt]->SetBranchAddress("Lep3ID", &Lep3ID);
    tc[tt]->SetBranchAddress("Lep1isID", &Lep1isID);
    tc[tt]->SetBranchAddress("Lep2isID", &Lep2isID);
    tc[tt]->SetBranchAddress("Lep3isID", &Lep3isID);
  }

  int nbins_KD = nbins_SystVariable[iSyst];
  double KD_limits[2] ={ KD_limits_SystVariable[iSyst][0], KD_limits_SystVariable[iSyst][1] };

  for (int iCR=0; iCR<2; iCR++){
    for (int icut=0; icut<nSIPCuts; icut++){
      for (int catZ1=0; catZ1<2; catZ1++){
        for (int catZ2=0; catZ2<2; catZ2++){
          for (int isocut=0; isocut<2; isocut++){
            for (int iM1cut=0; iM1cut<2; iM1cut++){
//              if (icut==2 || icut==3 || icut==5) continue;
              if (icut!=4 && icut!=0 && icut!=1 && iSyst!=21 && iSyst!=22 && iSyst!=23 && iSyst!=24) continue;
              if (icut>2 && iSyst==21) continue;
              if (icut!=0 && icut!=1 && icut!=3 && (iSyst==22 || iSyst==23)) continue;
              if (icut==0 && iSyst==24) continue;
              int option[6]={
                iCR,
                icut,
                catZ1,
                catZ2,
                isocut,
                iM1cut
              };

              TH1F* hfull[ntrees+1];
              TString appendName = produceCRname(iCR, isocut, iM1cut, catZ1, catZ2, icut);
              cout << "Now producing family " << appendName << endl;
              TString hname = appendName;
              hname.Append("_data");
              TH1F* hdata_full = new TH1F(hname, "", nbins_KD, KD_limits[0], KD_limits[1]);
              hdata_full->GetXaxis()->SetTitle(cSystVariable_label[iSyst]);
              hdata_full->GetYaxis()->SetTitle("Events / bin");
              for (int hh=0; hh<ntrees; hh++){
                hfull[hh] = (TH1F*)hdata_full->Clone();
                hname = appendName;
                if (hh<ntrees-1) hname.Append(Form("_%s", processName[hh]));
                else hname.Append("_MCTotal");
                hfull[hh]->SetName(hname);
                hfull[hh]->Sumw2();
              }
              hfull[ntrees] = hdata_full;
              hfull[ntrees]->Sumw2();

              for (int tt = 0; tt < ntrees; tt++){
                for (int ev = 0; ev < tc[tt]->GetEntries(); ev++){
                  MC_weight = 1;
                  MC_weight_noxsec = 1;

                  GenPrimaryVtx_x=0;
                  GenPrimaryVtx_y=0;
                  GenPrimaryVtx_z=0;
                  GenIntVtx_x=0;
                  GenIntVtx_y=0;
                  GenIntVtx_z=0;

                  GenHMass=125.6;
                  GenHPt=20;
                  GenHPhi=0;
                  Lep1_Z1SIP=0; Lep2_Z1SIP=0; Lep3_Z1SIP=0; KalmanCandVtx_chi2=0;
                  Lep1SIP=0; Lep2SIP=0; Lep3SIP=0;
                  Lep1ID=0; Lep2ID=0; Lep3ID=0;

                  tc[tt]->GetEntry(ev);

                  int LepID[3] ={ Lep1ID, Lep2ID, Lep3ID };
                  float Lep_Z1SIP[3] ={ Lep1_Z1SIP, Lep2_Z1SIP, Lep3_Z1SIP };
                  float LepSIP[3] ={ Lep1SIP, Lep2SIP, Lep3SIP };
                  double selection_weight = applyCRselection(
                    option,

                    Lep1combRelIsoPF, Lep2combRelIsoPF, Lep3combRelIsoPF,
                    Lep1isID, Lep2isID, Lep3isID,
                    Z1ids,
                    PFMET,

                    Lep_Z1SIP,
                    KalmanCandVtx_chi2,
                    LepSIP,
                    LepID,

                    Z1Mass, ZZMass
                    );
                  if (selection_weight==0)continue;

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
                  KD = KD*10000.; // in 1 um

                  if (iSyst==20){
                    KD = tanh(KD/800.); if (KD==1.) KD=0.999999; if (KD==-1.) KD=-0.999999;
                  }
                  if (iSyst==21){
                    KD = KalmanCandVtx_chi2; if (KD>=KD_limits_SystVariable[iSyst][1]) KD=KD_limits_SystVariable[iSyst][1]-0.0001;
                  }
                  if (iSyst==22){
                    KD = Lep1_Z1SIP; if (KD>=KD_limits_SystVariable[iSyst][1]) KD=KD_limits_SystVariable[iSyst][1]-0.0001; if (KD<=KD_limits_SystVariable[iSyst][0]) KD=KD_limits_SystVariable[iSyst][0]+0.0001;
                  }
                  if (iSyst==23){
                    KD = Lep3_Z1SIP; if (KD>=KD_limits_SystVariable[iSyst][1]) KD=KD_limits_SystVariable[iSyst][1]-0.0001; if (KD<=KD_limits_SystVariable[iSyst][0]) KD=KD_limits_SystVariable[iSyst][0]+0.0001;
                  }
                  if (iSyst==24){
                    KD = max(Lep1SIP,max(Lep2SIP,Lep3SIP)); if (KD>=KD_limits_SystVariable[iSyst][1]) KD=KD_limits_SystVariable[iSyst][1]-0.0001;
                  }

                  if (iSyst==4) KD = ZZPt/ZZMass;
                  if (iSyst==5) KD = (Txy-Txy_true)/delTxy;
                  if (iSyst==6) KD = (Dxy-Dxy_true)/delDxy;
                  if (iSyst==7) KD = ((OfflinePrimaryVtx_x - GenPrimaryVtx_x)*cos(GenHPhi) + (OfflinePrimaryVtx_y - GenPrimaryVtx_y)*sin(GenHPhi))/sigmaPV_xy;
                  if (iSyst==8) KD = ((KalmanCandVtx_x - GenIntVtx_x)*cos(GenHPhi) + (KalmanCandVtx_y - GenIntVtx_y)*sin(GenHPhi))/sigmaInt_xy;
                  if (iSyst==9) KD = (KalmanCandVtx_x-GenIntVtx_x)/sqrt(KalmanCandVtx_cov_xx);
                  if (iSyst==10) KD = (KalmanCandVtx_y-GenIntVtx_y)/sqrt(KalmanCandVtx_cov_yy);
                  if (iSyst==11) KD = (KalmanCandVtx_z-GenIntVtx_z)/sqrt(KalmanCandVtx_cov_zz);

                  float wgt = selection_weight*MC_weight;
                  if (tt<ntrees-1) wgt *= luminosity[EnergyIndex];
                  if (iSyst==22) wgt /= 2.;
                  if (tt<ntrees-1) hfull[tt]->Fill(KD, wgt);
                  else hfull[ntrees]->Fill(KD, wgt);

                  if (iSyst==22){
                    KD = Lep2_Z1SIP; if (KD>=KD_limits_SystVariable[iSyst][1]) KD=KD_limits_SystVariable[iSyst][1]-0.0001; if (KD<=KD_limits_SystVariable[iSyst][0]) KD=KD_limits_SystVariable[iSyst][0]+0.0001;
                    if (tt<ntrees-1) hfull[tt]->Fill(KD, wgt);
                    else hfull[ntrees]->Fill(KD, wgt);
                  }

                }
              }

              foutput->cd();
              for (int tt = 0; tt < ntrees-1; tt++) hfull[ntrees-1]->Add(hfull[tt]);
              plot1DNewOld(iSyst, iCR, isocut, iM1cut, catZ1, catZ2, icut, erg_tev, hfull, foutput, false);
              plot1DNewOld(iSyst, iCR, isocut, iM1cut, catZ1, catZ2, icut, erg_tev, hfull, foutput, true);
              for (int tt = 0; tt < ntrees+1; tt++){
                foutput->WriteTObject(hfull[tt]);
                delete hfull[tt];
              }
            }
          }
        }
      }
    }
  }
	for (int tt = 0; tt < ntrees; tt++){
		delete tc[tt];
	}

	foutput->Close();
}

void plot1DNewOld(int iSyst, int iCR, int isocut, int iM1cut, int catZ1, int catZ2, int icut, int erg_tev, TH1F** hfull, TFile* foutput, bool isNorm){
  TString comstring;
  comstring.Form("%iTeV", erg_tev);
  TString coutput_common = user_dir_hep + "Analysis/Plots/ZL_MC/";
  coutput_common.Append(Form("%s/", cSystVariable[iSyst]));
  coutput_common.Append(comstring);
  coutput_common.Append("/");
  if (isocut==0) coutput_common.Append("AllLoose");
  else coutput_common.Append(strLepIsoCut_label[isocut-1]);
  coutput_common.Append("/");
  if (iM1cut==0) coutput_common.Append("m1_All/");
  else if (iM1cut==1) coutput_common.Append("m1_ZPeak/");
  else if (iM1cut==2) coutput_common.Append("m1_ZGPeak/");

  TString canvasdir = coutput_common;
  canvasdir.Append(cutNames[icut]);
  canvasdir.Append("/");

  TString mkdirCommand = "mkdir -p ";
  mkdirCommand.Append(canvasdir);
  gSystem->Exec(mkdirCommand);

  TString cutlabel = cutLabel[icut];
  TH1F* htemp[4] ={
    (TH1F*)hfull[nProcesses]->Clone("datatemp"),
    (TH1F*)hfull[nProcesses-1]->Clone("mctemp"), // +DY+ttb
    (TH1F*)hfull[4]->Clone("mctemp2"), // +WZ+WW
    (TH1F*)hfull[0]->Clone("mctemp3") // qqb
  };

  double data_integral = htemp[0]->Integral(0, htemp[0]->GetNbinsX()+1);
  double mc_integral = htemp[1]->Integral(0, htemp[1]->GetNbinsX()+1);
  TGraphAsymmErrors* tgdata;
  int ndata = 0;
  double xx_data[1000];
  double xu_data[1000];
  double xd_data[1000];
  double yy_data[1000];
  double yu_data[1000];
  double yd_data[1000];
  const double quant = (1.0-0.6827)/2.0;
  double max_plot=htemp[1]->GetMaximum() + htemp[1]->GetBinError(htemp[1]->GetMaximumBin());
  if (isNorm) max_plot/=mc_integral;
  else max_plot/=1000.;

  for (int bin=1; bin<=htemp[0]->GetNbinsX(); bin++){
    double bincenter = htemp[0]->GetBinCenter(bin);
    double bincontent = htemp[0]->GetBinContent(bin);

    if (bincontent>=0){
      xx_data[ndata] = bincenter;
      yy_data[ndata] = bincontent;
      xu_data[ndata] = 0;
      xd_data[ndata] = 0;
      yu_data[ndata] = ROOT::Math::chisquared_quantile_c(quant, 2*(bincontent+1))/2.-bincontent;
      yd_data[ndata] = (bincontent==0) ? 0 : (bincontent-ROOT::Math::chisquared_quantile_c(1-quant, 2*bincontent)/2.);
      if (isNorm){
        yy_data[ndata] /= data_integral;
        yu_data[ndata] /= data_integral;
        yd_data[ndata] /= data_integral;
      }
      else{
        yy_data[ndata] /= 1000.;
        yu_data[ndata] /= 1000.;
        yd_data[ndata] /= 1000.;
      }
      double high_data = yy_data[ndata]+yu_data[ndata];
      if (high_data > max_plot) max_plot = high_data;
      ndata++;
    }
  }
  tgdata = new TGraphAsymmErrors(ndata, xx_data, yy_data, xd_data, xu_data, yd_data, yu_data);
  tgdata->SetMarkerSize(0.8);
  tgdata->SetMarkerStyle(20);
  tgdata->SetMarkerColor(kBlack);
  tgdata->SetLineColor(kBlack);
  tgdata->SetLineWidth(2);

  htemp[2]->Add(hfull[5]);
  htemp[2]->Add(htemp[3]);

  if (isNorm){
    htemp[0]->Scale(1./data_integral);
    htemp[0]->GetYaxis()->SetTitle("Rate / bin");
  }
  else{
    htemp[0]->Scale(1./1000.);
    htemp[0]->GetYaxis()->SetTitle("Events / bin / 1000");
  }
  htemp[0]->SetLineWidth(2);
  for (int hh=3; hh>0; hh--){
    if (isNorm){
      htemp[hh]->Scale(1./mc_integral);
      htemp[hh]->GetYaxis()->SetTitle("Rate / bin");
    }
    else{
      htemp[hh]->Scale(1./1000.);
      htemp[hh]->GetYaxis()->SetTitle("Events / bin / 1000");
    }
    htemp[hh]->SetLineWidth(2);
    if(hh!=1) htemp[hh]->SetFillStyle(1001);
  }
  TString identifier_label = "#frac{";
  identifier_label = identifier_label + strZ1Category_label[catZ1];
  identifier_label = identifier_label + strLep3Category_label[catZ2];
  identifier_label = identifier_label + ", " + cutlabel;
  identifier_label = identifier_label + "}{";
  if (isocut>1) identifier_label = identifier_label + strLepIsoCut_label[isocut-1];
  else if (isocut==1) identifier_label = identifier_label + "3p";
  if (isocut!=0) identifier_label = identifier_label + ", ";
  identifier_label = identifier_label + Z1masswidth_label[iM1cut] + ", " + strPFMET_label[iCR] + "}";

  TString canvasname = produceCRname(iCR, isocut, iM1cut, catZ1, catZ2, icut);
  canvasname = canvasname + "_" + comstring;
  canvasname.Prepend("cCompare_DataMC_");
  if (isNorm) canvasname.Append("_AreaNormalized");

  foutput->cd();

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

  TLegend *l2D = new TLegend(0.66, 0.69, 0.98, 0.90);
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

  htemp[1]->GetYaxis()->SetRangeUser(0, max_plot * 1.5);
  htemp[1]->GetXaxis()->SetNdivisions(505);
  htemp[1]->GetXaxis()->SetLabelFont(42);
  htemp[1]->GetXaxis()->SetLabelOffset(0.007);
  htemp[1]->GetXaxis()->SetLabelSize(0.04);
  htemp[1]->GetXaxis()->SetTitleSize(0.06);
  htemp[1]->GetXaxis()->SetTitleOffset(0.9);
  htemp[1]->GetXaxis()->SetTitleFont(42);
  htemp[1]->GetYaxis()->SetNdivisions(505);
  htemp[1]->GetYaxis()->SetLabelFont(42);
  htemp[1]->GetYaxis()->SetLabelOffset(0.007);
  htemp[1]->GetYaxis()->SetLabelSize(0.04);
  htemp[1]->GetYaxis()->SetTitleSize(0.06);
  htemp[1]->GetYaxis()->SetTitleOffset(1.1);
  htemp[1]->GetYaxis()->SetTitleFont(42);

  htemp[0]->SetLineColor(kBlack);
  htemp[1]->SetLineColor(kOrange + 10);
//  htemp[2]->SetLineColor(kBlack);
//  htemp[3]->SetLineColor(kBlack);
  htemp[2]->SetLineColor(kAzure-2);
  htemp[3]->SetLineColor(TColor::GetColor("#669966"));
  htemp[2]->SetFillColor(kAzure-2);
  htemp[3]->SetFillColor(TColor::GetColor("#669966"));

  TH1F* hclone_mc = (TH1F*)htemp[1]->Clone("clone_mc");
  hclone_mc->SetFillColor(kOrange + 10);
  hclone_mc->SetFillStyle(3005);

  l2D->AddEntry(tgdata, "Observed", "ep");
  l2D->AddEntry(htemp[3], "q#bar{q}", "f");
  l2D->AddEntry(htemp[2], "+WW/Z+jets", "f");
  l2D->AddEntry(hclone_mc, "+DY+t#bar{t}", "fl");

  htemp[1]->Draw("hist");
  htemp[2]->Draw("histsame");
  htemp[3]->Draw("histsame");
  hclone_mc->Draw("e2same");
  tgdata->Draw("e1psame");

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

  foutput->WriteTObject(c2D);

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

  delete tgdata;
  delete hclone_mc;
  for(int hh=0;hh<4;hh++) htemp[hh];
  c2D->Close();
  foutput->cd();
}
