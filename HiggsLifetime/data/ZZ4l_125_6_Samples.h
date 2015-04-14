#include <iostream>
#include <string>
#include "TMath.h"
#include "TH1D.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TDirectory.h"
#include "TCanvas.h"

const float PI_VAL = TMath::Pi();

enum sample {
	k0um,
	k100um,
	k500um,
	k1000um,

	kNumSamples
};
string sampleName[kNumSamples] = {
	"0microns",
	"100microns",
	"500microns",
	"1000microns"
};
char* sample_suffix[kNumSamples]={
	"0microns",
	"100microns",
	"500microns",
	"1000microns"
};
double sample_CTau[kNumSamples] = {
	0,
	100,
	500,
	1000
};
char* user_folder[4]={
	"4mu",
	"4e",
	"2mu2e",
	"data"
};
const int kAllSamples=23;
const int kGGHSamples=4;
const int kGGMCFMSamples=9;
const int kGGOLDSamples=6;
const int kGGSamples=11;
const int kQQBZZSamples=17;
const int kQQBZZSamples_Dedicated=22;
char* sample_FullSim[kAllSamples+1]={
	"HZZ4lTree_powheg15jhuGenV3-0PMH125.6",
	"HZZ4lTree_powheg15jhuGenV3-CTau100-0PMH125.6",
	"HZZ4lTree_powheg15jhuGenV3-CTau500-0PMH125.6",
	"HZZ4lTree_powheg15jhuGenV3-CTau1000-0PMH125.6",

	"HZZ4lTree_ggZZ4l",
	"HZZ4lTree_ggZZ2l2l",

	"HZZ4lTree_ggTo4mu_Contin-MCFM67",
	"HZZ4lTree_ggTo4e_Contin-MCFM67",
	"HZZ4lTree_ggTo2e2mu_Contin-MCFM67",

	"HZZ4lTree_ggTo2l2l_Continuum",
	"HZZ4lTree_ggTo4l_Continuum",

	"HZZ4lTree_ZZTo4mu",
	"HZZ4lTree_ZZTo4e",
	"HZZ4lTree_ZZTo2e2mu",
	"HZZ4lTree_ZZTo2mu2tau",
	"HZZ4lTree_ZZTo2e2tau",
	"HZZ4lTree_ZZTo4tau",
	"HZZ4lTree_ZZ95-160To2e2mu",
	"HZZ4lTree_ZZ95-160To2mu2tau",
	"HZZ4lTree_ZZ95-160To4e",
	"HZZ4lTree_ZZ95-160To4mu",
	"HZZ4lTree_ZZ95-160To4tau",

	"HZZ4lTree_DoubleOr_CRZLLTree",
	"HZZ4lTree_DoubleOr_CRZLTree"
};
char* data_files[3]={
	"HZZ4lTree_DoubleMu",
	"HZZ4lTree_DoubleEle",
	"HZZ4lTree_DoubleOr"
};
const int nCRZLLMC = 7;
char* sample_CR_MC[nCRZLLMC]={
  "HZZ4lTree_DYJetsToLLTuneZ2M10-NoB",
  "HZZ4lTree_DYJetsToLLTuneZ2M50-NoB",
  "HZZ4lTree_DYJetsToLLTuneZ2M10-B",
  "HZZ4lTree_DYJetsToLLTuneZ2M50-B",
  "HZZ4lTree_TTTo2L2Nu2B",
  "HZZ4lTree_WZ",
  "HZZ4lTree_WWJets"
};


string user_dir = "/scratch0/hep/usarical/HiggsLifetime/";
string user_dir_hep = "/scratch0/hep/usarical/HiggsLifetime/";
string user_dir_hhpc = "/scratch0/hhpc-hn1/usarical/ggtoH-PWGSamples-125_6/";

double yield_offshell_zx[2][3] ={
  { 0.1000, 0.4000, 0.34 },
  { 0.5500, 1.7800, 1.38 }
};
double yield_offshell_ggzz[2][3] ={
  { 1.7132, 1.1836, 2.8654 },
  { 8.9407, 6.2543, 14.8431 }
};
double yield_80_100_zx[2][3] ={
  { 0.0103, 0.0522, 0.1402 },
  { 0.0549, 0.2320, 0.5663 }
};
double yield_80_100_ggzz[2][3] ={
  { 0.0222, 0.0106, 0.0113 },
  { 0.2441, 0.0918, 0.0649 }
};
double yield_signal_zx[2][3] ={
  { 0.2230, 0.3412, 1.0249 },
  { 1.1878, 1.5167, 4.1398 }
};
double yield_signal_ggzz[2][3] ={
  { 0.0625, 0.0341, 0.0741 },
  { 0.4131, 0.2041, 0.5005 }
};
double yield_signal_qqzz[2][3] ={
  { 1.7971, 0.8386, 2.2456 },
  { 7.6478, 2.9364, 8.8585 }
};

double yield_signal_higgs[2][3] ={
  { 1.243885272, 0.695146432, 1.666241752 },
  { 5.947054319, 3.089788181, 7.680661487 }
};
double yield_signal_ggh[2][3]={
  { 1.0902689, 0.6087736, 1.4452079 },
  { 5.1963998, 2.6944639, 6.6562629 }
};
double yield_signal_vbfh[2][3]={
  { 0.092458836, 0.051755897, 0.128619212 },
  { 0.467988070, 0.247885532, 0.617816891 }
};
double yield_signal_wh[2][3]={
  { 0.032020833, 0.018214085, 0.048193985 },
  { 0.143887487, 0.07605782, 0.208553654 }
};
double yield_signal_zh[2][3]={
  { 0.024986956, 0.013932042, 0.037813666 },
  { 0.114986662, 0.058084659, 0.162819182 }
};
double yield_signal_tth[2][3]={
  { 0.004149795, 0.002470837, 0.006407011 },
  { 0.023792303, 0.01329626, 0.035208888 }
};

double luminosity[2] ={ 5.051, 19.712 };
double XSEC_ggH[2] ={ 14.99, 19.09 };
double XSEC_VBFH[2] ={ 1.214, 1.572 };
double XSEC_WH[2] ={ 0.5688, 0.6931 };
double XSEC_ZH[2] ={ 0.3299, 0.4091 };
double XSEC_ttH[2] ={ 0.08508, 0.1274 };
double BR_flavor[4] ={ 3.45e-05, 3.45e-05, 6.27e-05, 0.0279 }; // 4mu, 4e, 2e2mu, ZZ

const int kNumFiles_GGHVVBSM=24;
char* sample_FullSim_GGHVVBSM[kNumFiles_GGHVVBSM]={
  "HZZ4lTree_powheg15jhuGenV3-0PMH125.6",
  "HZZ4lTree_powheg15jhuGenV3-0PHH125.6",
  "HZZ4lTree_powheg15jhuGenV3-0MH125.6",
  "HZZ4lTree_powheg15jhuGenV3-0L1H125.6",
  "HZZ4lTree_powheg15jhuGenV3-0PHf05ph0H125.6",
  "HZZ4lTree_powheg15jhuGenV3-0Mf05ph0H125.6",
  "HZZ4lTree_powheg15jhuGenV3-0PHf05ph0Mf05ph0H125.6",
  "HZZ4lTree_powheg15jhuGenV3-0L1f05ph180H125.6",
  "HZZ4lTree_powheg15jhuGenV3-0PHf033ph0Mf033ph0H125.6",
  "HZZ4lTree_powheg15jhuGenV3-0PHf01ph0H125.6",
  "HZZ4lTree_powheg15jhuGenV3-0Mf01ph0H125.6",
  "HZZ4lTree_powheg15jhuGenV3-0PHf01ph0Mf01ph0H125.6",
  "HZZ4lTree_powheg15jhuGenV3-0L1f05ph0H125.6",
  "HZZ4lTree_powheg15jhuGenV3-0L1f01ph0H125.6",
  "HZZ4lTree_powheg15jhuGenV3-0PHf05ph90H125.6",
  "HZZ4lTree_powheg15jhuGenV3-0Mf05ph90H125.6",
  "HZZ4lTree_powheg15jhuGenV3-0PHf05ph0Mf05ph90H125.6",
  "HZZ4lTree_powheg15jhuGenV3-0PHf033ph0Mf033ph90H125.6",
  "HZZ4lTree_powheg15jhuGenV3-0PHf01ph90H125.6",
  "HZZ4lTree_powheg15jhuGenV3-0Mf01ph90H125.6",
  "HZZ4lTree_powheg15jhuGenV3-0PHf01ph0Mf01ph90H125.6",
  "HZZ4lTree_powheg15jhuGenV3-0PHf05ph180H125.6",
  "HZZ4lTree_powheg15jhuGenV3-0Mf05ph180H125.6",
  "HZZ4lTree_powheg15jhuGenV3-0PHf05ph180Mf05ph0H125.6"
};
const int kNumFiles_VBFHVVBSM=3;
char* sample_FullSim_VBFHVVBSM[kNumFiles_VBFHVVBSM]={
  "HZZ4lTree_VBF0P_H125.6",
  "HZZ4lTree_VBF0PH_H125.6",
  "HZZ4lTree_VBF0M_H125.6"
};
const int kNumFiles_WHVVBSM=2;
char* sample_FullSim_WHVVBSM[kNumFiles_WHVVBSM]={
  "HZZ4lTree_WHiggs0P_H125.6",
  "HZZ4lTree_WHiggs0PH_H125.6"
};
const int kNumFiles_ZHVVBSM=3;
char* sample_FullSim_ZHVVBSM[kNumFiles_ZHVVBSM]={
  "HZZ4lTree_ZHiggs0P_H125.6",
  "",
  "HZZ4lTree_ZHiggs0M_H125.6"
};

enum sample_GGHVVBSM {
  kfg2_0_fg4_0,
  kfg2_1_fg4_0,
  kfg2_0_fg4_1,
  kfLambda1_1,

  kfg2_05_fg4_0,
  kfg2_0_fg4_05,
  kfg2_05_fg4_05,
  kfLambda1_m05,

  kfg2_33_fg4_33,
  kfg2_01_fg4_0,
  kfg2_0_fg4_01,
  kfg2_01_fg4_01,

  kfLambda1_05,
  kfLambda1_03, // No sample here
  kfLambda1_01,

  kfg2_05_fg4_0_p290,
  kfg2_0_fg4_05_p390,
  kfg2_05_fg4_05_p390,

  kfLambda1_05_pL190,

  kfg2_33_fg4_33_p390,
  kfg2_01_fg4_0_p290,
  kfg2_0_fg4_01_p390,
  kfg2_01_fg4_01_p390,

  kfg2_05_fg4_0_p2Pi,
  kfg2_0_fg4_05_p3Pi,
  kfg2_05_fg4_05_p2Pi,

  kfg2_05_fLambda1_m05,
  kfg2_05_fLambda1_m05_p2270,
  kfg2_05_fLambda1_05,
  kfg2_33_fLambda1_33_p2Pi,

  kfg4_05_fLambda1_m05,
  kfg4_05_fLambda1_m05_p3270,
  kfg4_05_fLambda1_05,
  kfg4_33_fLambda1_33_p3Pi,

  // ZZ LambdaQ and ZZ imaginary g1
  kfLambdaQ_1,
  kfLambdaQ_1_pLQ90,
  kfg1_p190,

  kfg1_05_fLambdaQ_05,
  kfg1_05_fLambdaQ_05_pLQ90,

  // ZG and GG
  kfZG_1_fGG_0,
  kfZG_0_fGG_1,
  kfZG_05_fGG_0,
  kfZG_0_fGG_05,

  kfMZG_1_fMGG_0,
  kfMZG_0_fMGG_1,
  kfMZG_05_fMGG_0,
  kfMZG_0_fMGG_05,

  kfZG_05_fMZG_05,
  kfGG_05_fMGG_05,

  kfZG_SM_fGG_SM,

  kfZG_Lambda1_1,
  kfZG_Lambda1_05,

  kNumGGHVVBSM
};

TString sample_GGHVVBSM_label[kNumGGHVVBSM]={
  "SM",
  "f_{a2}=1",
  "f_{a3}=1",
  "f_{#Lambda1}=1",

  "f_{a2}=0.5 #phi_{a2}=0",
  "f_{a3}=0.5 #phi_{a3}=0",
  "f_{a2}=f_{a3}=0.5 #phi_{a2}=#phi_{a3}=0",
  "f_{#Lambda1}=-0.5",

  "f_{a2}=f_{a3}=0.33 #phi_{a2}=#phi_{a3}=0",
  "f_{a2}=0.1 #phi_{a2}=0",
  "f_{a3}=0.1 #phi_{a3}=0",
  "f_{a2}=f_{a3}=0.1 #phi_{a2}=#phi_{a3}=0",

  "f_{#Lambda1}=0.5",
  "f_{#Lambda1}=0.3",
  "f_{#Lambda1}=0.1",

  "f_{a2}=0.5 #phi_{a2}=#pi/2",
  "f_{a3}=0.5 #phi_{a3}=#pi/2",
  "f_{a2}=f_{a3}=0.5 #phi_{a3}=#pi/2",

  "f_{#Lambda1}=0.5 #phi_{#Lambda1}=#pi/2",

  "f_{a2}=f_{a3}=0.33 #phi_{a3}=#pi/2",
  "f_{a2}=0.1 #phi_{a2}=#pi/2",
  "f_{a3}=0.1 #phi_{a3}=#pi/2",
  "f_{a2}=f_{a3}=0.1 #phi_{a3}=#pi/2",

  "f_{a2}=0.5 #phi_{a2}=#pi",
  "f_{a3}=0.5 #phi_{a3}=#pi",
  "f_{a2}=f_{a3}=0.5 #phi_{a2}=#pi",

  "f_{a2}=f_{#Lambda1}=0.5 #phi_{#Lambda1}=#pi",
  "f_{a2}=f_{#Lambda1}=0.5 #phi_{a2}=3#pi/2 #phi_{#Lambda1}=#pi",
  "f_{a2}=f_{#Lambda1}=0.5 #phi_{#Lambda1}=0",
  "f_{a2}=f_{#Lambda1}=0.33 #phi_{a2}=#pi",

  "f_{a3}=f_{#Lambda1}=0.5 #phi_{#Lambda1}=#pi",
  "f_{a3}=f_{#Lambda1}=0.5 #phi_{a3}=3#pi/2 #phi_{#Lambda1}=#pi",
  "f_{a3}=f_{#Lambda1}=0.5 #phi_{#Lambda1}=0",
  "f_{a3}=f_{#Lambda1}=0.33 #phi_{a3}=#pi",

  // ZZ LambdaQ and ZZ imaginary g1
  "f_{#Lambda Q}=1",
  "f_{#Lambda Q}=1 #phi_{#Lambda Q}=#pi/2",
  "f_{a1}=1 #phi_{a1}=#pi/2",

  "f_{#Lambda Q}=0.5",
  "f_{#Lambda Q}=0.5 #phi_{#Lambda Q}=#pi/2",

  // ZG and GG
  "f_{Z#gamma}=1, f_{#gamma#gamma}=0",
  "f_{Z#gamma}=0, f_{#gamma#gamma}=1",
  "f_{Z#gamma}=0.5, f_{#gamma#gamma}=0",
  "f_{Z#gamma}=0, f_{#gamma#gamma}=0.5",

  "f_{0M-Z#gamma}=1, f_{0M-#gamma#gamma}=0",
  "f_{0M-Z#gamma}=0, f_{0M-#gamma#gamma}=1",
  "f_{0M-Z#gamma}=0.5, f_{0M-#gamma#gamma}=0",
  "f_{0M-Z#gamma}=0, f_{0M-#gamma#gamma}=0.5",

  "f_{Z#gamma}=0.5, f_{0M-Z#gamma}=0.5",
  "f_{#gamma#gamma}=0.5, f_{0M-#gamma#gamma}=0.5",

  "f_{Z#gamma}, f_{#gamma#gamma} Full SM",

  "f_{Z#gamma, #Lambda1}=1",
  "f_{Z#gamma, #Lambda1}=0.5"
};

