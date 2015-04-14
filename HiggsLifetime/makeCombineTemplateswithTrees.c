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


void makeCombineSignalTemplateswithTrees_GenSmeared_single(int folder, int erg_tev);
bool testLepId_GenLevel(int channel, int Lep1ID, int Lep2ID, int Lep3ID, int Lep4ID);

double formula_TxySyst(double* x, double* par);

//void makeCombineSignalTemplateswithTrees_FullSim_single(int folder, int erg_tev, int Systematics);
void produce_KDShapeVariation(int folder, int erg_tev);

void produce_KDShapeVariation_Untransformed(int folder, int erg_tev);


void makeCombineSignalTemplateswithTrees(int isRecoLevel=0){
	const int kNumSyst=7;
	int systematics[kNumSyst]={0,1,-1,2,-2,3,-3};
//	for(int CoM=7;CoM<9;++CoM){
	for(int CoM=8;CoM<9;++CoM){
		for(int channel=0;channel<3;channel++){
			if(isRecoLevel==0) makeCombineSignalTemplateswithTrees_GenSmeared_single(channel,CoM);	
		}
	}
}

void makeCombineSignalTemplateswithTrees_GenSmeared_single(int folder, int erg_tev){
	char TREE_NAME[]="SelectedTree";
	TString INPUT_NAME = "HZZ4lTree_powheg15jhuGenV3-CTau";
	TString OUTPUT_NAME = "_templates_GenLevel.root";
	OUTPUT_NAME.Prepend(user_folder[folder]);
	TString comstring;
	comstring.Form("%iTeV",erg_tev);
	TString erg_dir;
	erg_dir.Form("LHC_%iTeV",erg_tev);

	const int nSyst=7;
	int kSyst=nSyst;
	const int kSystNominal=3;
	const int iSystRes=1;
	const int iSystScale=2;
	const int iSystTxy=3;
	char* cSyst[nSyst] = {
		"_TxyDown","_ScaleDown","_ResDown",
		"",
		"_ResUp","_ScaleUp","_TxyUp"
	};

	int EnergyIndex=1;
	if(erg_tev==7) EnergyIndex=0;
	float ZZMass_PeakCut[2]={105.6,140.6};
	//use ggH yields directly
 	double nSM_ScaledPeak[2][3]={
 		{1.0902689,0.6087736,1.4452079},
 		{5.1963998,2.6944639,6.6562629}
	};
	double luminosity[2] = { 5.051, 19.712 };
	for (int e = 0; e < 2; e++){ for (int ff = 0; ff < 3; ff++) nSM_ScaledPeak[e][ff] /= luminosity[e]; }

	float templateWeight[nSyst] = { 1 };
	double ctauWeight=1;
	float GenCTau;
	float MC_weight;
	float ZZMass,ZZPt,ZZEta,ZZPhi;
	float CandVtx_x,CandVtx_y,CandVtx_z;

	int Lep1ID;
	int Lep2ID;
	int Lep3ID;
	int Lep4ID;

	int FileID = -1;

	float p0plus_VAJHU;
	float p0plus_VAMCFM;
	float bkg_VAMCFM;
	float p0plus_m4l, p0plus_m4l_ScaleUp, p0plus_m4l_ScaleDown, p0plus_m4l_ResUp, p0plus_m4l_ResDown;
	float bkg_m4l, bkg_m4l_ScaleUp, bkg_m4l_ScaleDown, bkg_m4l_ResUp, bkg_m4l_ResDown;

	float Txy=-99;
	float KD=-99;
	float D_bkg=-99;
	float D_bkg_ScaleUp=-99;
	float D_bkg_ScaleDown=-99;
	float D_bkg_ResUp=-99;
	float D_bkg_ResDown=-99;

	double ctau_limits[2] = { 0, 1000 };
	const int nctaus = 40;
	double ctau_gridsize = (ctau_limits[1]-ctau_limits[0])/nctaus;
//	int nbinsx = 400;
//	double xlow = -10000, xhigh = 10000;
	int nbinsx = 80;
	double xlow = -1, xhigh = 1;
	int nbinsy = 50;
	double ylow=0,yhigh=1;

	TString INPUT_SYST_NAME = "LifetimeKD_ShapeVariation_DataQQBZZ_Untransformed_Comparison.root";
	TString cinput_syst = user_dir_hep + "Analysis/Plots/" + INPUT_SYST_NAME;
	TFile* finput_syst = new TFile(cinput_syst, "read");
	TF1* txy_syst = (TF1*) finput_syst->Get(Form("tgratio_%s_fit",user_folder[folder])); // [0]-[1]*exp(-pow((x-[2])/[3],2))
	const double Txy_shift_base = txy_syst->GetParameter(2); // in 2 mm
	const double Txy_smear_base = txy_syst->GetParameter(3); // in 2 mm

	TString cinput_common = user_dir_hep + erg_dir + "/";
	TString coutput_common = user_dir_hep + "Analysis/Templates/" + comstring + "/";
	TString coutput = coutput_common + OUTPUT_NAME;
	TFile* foutput = new TFile(coutput,"recreate");
	cout << "Starting to create " << coutput << endl;

	TChain* tree[kNumSamples] = { 0 };
	for (int f = 0; f < kNumSamples; f++){
		TString cinput = cinput_common + sampleName[f] + "/" + INPUT_NAME + sample_suffix[f] + "-0PMH125.6_GenLevel_Reprocessed.root";
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
		tree[f]->SetBranchAddress("CandVtx_x", &CandVtx_x);
		tree[f]->SetBranchAddress("CandVtx_y", &CandVtx_y);
		tree[f]->SetBranchAddress("CandVtx_z", &CandVtx_z);
		tree[f]->SetBranchAddress("Lep1ID", &Lep1ID);
		tree[f]->SetBranchAddress("Lep2ID", &Lep2ID);
		tree[f]->SetBranchAddress("Lep3ID", &Lep3ID);
		tree[f]->SetBranchAddress("Lep4ID", &Lep4ID);
		tree[f]->SetBranchAddress("GenCTau", &GenCTau);
	}

	for (int genctau = 0; genctau <= nctaus; genctau++){
		double target_ctau = ctau_limits[0] + ctau_gridsize*genctau;
		double Txy_smear = sqrt(pow(Txy_smear_base,2)+2.*pow(target_ctau,2)); // 2* appears since smear**2 = 2*RMS**2
		double Txy_shift = Txy_smear/Txy_smear_base*Txy_shift_base + target_ctau;
		txy_syst->SetParameter(2,Txy_shift);txy_syst->SetParameter(3,Txy_smear);
		cout << "Template for cTau = " << target_ctau << " of " << user_folder[folder] << " @ " << comstring << endl;
		cout << "Txy shift = " << Txy_shift << ", Txy smear " << Txy_smear << endl;

		TH2F* hfill[nSyst];
		for (int ss = 0; ss < nSyst; ss++){
			hfill[ss] = new TH2F(Form("H_2D_CTau%.0f%s", target_ctau, cSyst[ss]), "", nbinsx, xlow, xhigh, nbinsy, ylow, yhigh);
			hfill[ss]->SetTitle( Form("%.0f #mum %s at %i TeV", target_ctau, user_folder[folder], erg_tev) );
			hfill[ss]->GetXaxis()->SetTitle("tanh(T_{xy} / 2 mm)");
			hfill[ss]->GetYaxis()->SetTitle("#it{D}_{bkg}");
		}
		TH2F* hfill_conditional;
		TH1F* hfill_projX;
		TTree* mytree = new TTree(Form("T_2D_CTau%.0f", target_ctau), "");
		mytree->Branch("kSystematics",&kSyst);
		mytree->Branch("MC_weight",templateWeight,"MC_weight[kSystematics]/F");
		mytree->Branch("D_bkg",&D_bkg);
		mytree->Branch("D_bkg_ScaleUp",&D_bkg_ScaleUp);
		mytree->Branch("D_bkg_ScaleDown",&D_bkg_ScaleDown);
		mytree->Branch("D_bkg_ResUp",&D_bkg_ResUp);
		mytree->Branch("D_bkg_ResDown",&D_bkg_ResDown);
		mytree->Branch("KD",&KD);
		mytree->Branch("FileID",&FileID);

		TH2F* hCounters = new TH2F(Form("hCounters_2D_CTau%.0f", target_ctau), "", kNumSamples, 0, kNumSamples, 2, 0, 2);
		hCounters->SetOption("colz");
		for (int f = 0; f < kNumSamples; f++){
			if (f == 0 && genctau>0) continue;
			if (f > 0 && genctau==0) break;
			double initial_ctau = sample_CTau[f];
			double effectiveInvCTau = (f>0 ? 1. / target_ctau - 1. / initial_ctau : 0);

			for (int ev = 0; ev < tree[f]->GetEntries(); ev++){
				tree[f]->GetEntry(ev);

				templateWeight[kSystNominal] = 1;
				ctauWeight = exp(-GenCTau*effectiveInvCTau);
				templateWeight[kSystNominal] *= ctauWeight;

				hCounters->AddBinContent(hCounters->GetBin(f+1,1),1.);
				hCounters->AddBinContent(hCounters->GetBin(f+1,2),templateWeight[kSystNominal]);
				templateWeight[kSystNominal]=1;
			}
		}
		foutput->WriteTObject(hCounters);

		for (int f = 0; f < kNumSamples; f++){
			if (f == 0 && genctau>0) continue;
			if (f > 0 && genctau==0) continue;
			double initial_ctau = sample_CTau[f];
			if(initial_ctau<target_ctau) continue;
			double effectiveInvCTau = (f>0 ? 1. / target_ctau - 1. / initial_ctau : 0);
			double xsec_ctau_wgt = hCounters->GetBinContent(f+1, 1) / hCounters->GetBinContent(f+1, 2);
			cout << "CTau xsec reweighting from " << initial_ctau << " to " << target_ctau << " = " << xsec_ctau_wgt << endl;

			for (int ev = 0; ev < tree[f]->GetEntries(); ev++){
				tree[f]->GetEntry(ev);
				if (!(ZZMass >= ZZMass_PeakCut[0] && ZZMass < ZZMass_PeakCut[1])) continue;
				if (!(testLepId_GenLevel(folder, Lep1ID, Lep2ID, Lep3ID, Lep4ID))) continue;

				Txy = (CandVtx_x*cos(ZZPhi) + CandVtx_y*sin(ZZPhi))*ZZMass / ZZPt;

				D_bkg = p0plus_VAJHU*p0plus_m4l / (p0plus_VAJHU*p0plus_m4l + bkg_VAMCFM*bkg_m4l);
				D_bkg_ScaleUp = (p0plus_VAJHU*p0plus_m4l_ScaleUp) / (p0plus_VAJHU*p0plus_m4l_ScaleUp + bkg_VAMCFM*bkg_m4l_ScaleUp);
				D_bkg_ScaleDown = (p0plus_VAJHU*p0plus_m4l_ScaleDown) / (p0plus_VAJHU*p0plus_m4l_ScaleDown + bkg_VAMCFM*bkg_m4l_ScaleDown);
				D_bkg_ResUp = (p0plus_VAJHU*p0plus_m4l_ResUp) / (p0plus_VAJHU*p0plus_m4l_ResUp + bkg_VAMCFM*bkg_m4l_ResUp);
				D_bkg_ResDown = (p0plus_VAJHU*p0plus_m4l_ResDown) / (p0plus_VAJHU*p0plus_m4l_ResDown + bkg_VAMCFM*bkg_m4l_ResDown);

				float smd[nSyst] = {
					D_bkg,D_bkg_ScaleDown,D_bkg_ResDown,
					D_bkg,
					D_bkg_ResUp,D_bkg_ScaleUp,D_bkg
				};
				for (int ss = 0; ss < nSyst; ss++){
					KD = Txy / 2000.; // in 2 mm
					double TxySystwgt = 1;
					if (ss == (kSystNominal + iSystTxy)) TxySystwgt = txy_syst->Eval(KD);
					if (ss == (kSystNominal - iSystTxy)) TxySystwgt = 1. / txy_syst->Eval(KD);
					KD = tanh(KD);
					if (KD >= 1.) KD = 1. - 1.0e-10;
					if (KD <= -1.) KD = -1. + 1.0e-10;

					templateWeight[ss] = TxySystwgt;
					ctauWeight = exp(-GenCTau*effectiveInvCTau);
					templateWeight[ss] *= ctauWeight*xsec_ctau_wgt;
					hfill[ss]->Fill(KD, smd[ss], templateWeight[ss]);
				}
			}
		}
		double integral_rescale[nSyst] = { 1 };
		for (int ss = 0; ss < nSyst; ss++){
			hfill[ss]->SetOption("colz");
			double integral_rewgt = hfill[ss]->Integral();
			double integral_target = nSM_ScaledPeak[EnergyIndex][folder];
			integral_rescale[ss] = integral_target / integral_rewgt;
			cout << "Template for cTau = " << target_ctau << " of " << user_folder[folder] << " @ " << comstring << " is rescaled with " << integral_rescale[ss] << " for systematic " << ss << "..." << endl;
			hfill[ss]->Scale(integral_rescale[ss]);
		}

		for (int f = 0; f < kNumSamples; f++){
			if (f == 0 && genctau>0) continue;
			if (f > 0 && genctau==0) continue;
			double initial_ctau = sample_CTau[f];
			if(initial_ctau<target_ctau) continue;
			double effectiveInvCTau = (f>0 ? 1. / target_ctau - 1. / initial_ctau : 0);
			double xsec_ctau_wgt = hCounters->GetBinContent(f+1, 1) / hCounters->GetBinContent(f+1, 2);

			for (int ev = 0; ev < tree[f]->GetEntries(); ev++){
				tree[f]->GetEntry(ev);

				if( !(ZZMass>=ZZMass_PeakCut[0] && ZZMass<ZZMass_PeakCut[1]) ) continue;
				if( !( testLepId_GenLevel(folder,Lep1ID,Lep2ID,Lep3ID,Lep4ID) ) ) continue;

				FileID = f;

				Txy = ( CandVtx_x*cos(ZZPhi) + CandVtx_y*sin(ZZPhi) )*ZZMass/ZZPt;
//				Tz = CandVtx_z*ZZMass/( ZZPt*sinh(ZZEta) );

				D_bkg = p0plus_VAJHU*p0plus_m4l / (p0plus_VAJHU*p0plus_m4l + bkg_VAMCFM*bkg_m4l);
				D_bkg_ScaleUp = (p0plus_VAJHU*p0plus_m4l_ScaleUp) / (p0plus_VAJHU*p0plus_m4l_ScaleUp + bkg_VAMCFM*bkg_m4l_ScaleUp);
				D_bkg_ScaleDown = (p0plus_VAJHU*p0plus_m4l_ScaleDown) / (p0plus_VAJHU*p0plus_m4l_ScaleDown + bkg_VAMCFM*bkg_m4l_ScaleDown);
				D_bkg_ResUp = (p0plus_VAJHU*p0plus_m4l_ResUp) / (p0plus_VAJHU*p0plus_m4l_ResUp + bkg_VAMCFM*bkg_m4l_ResUp);
				D_bkg_ResDown = (p0plus_VAJHU*p0plus_m4l_ResDown) / (p0plus_VAJHU*p0plus_m4l_ResDown + bkg_VAMCFM*bkg_m4l_ResDown);

				float smd[nSyst] = {
					D_bkg,D_bkg_ScaleDown,D_bkg_ResDown,
					D_bkg,
					D_bkg_ResUp,D_bkg_ScaleUp,D_bkg
				};
				for (int ss = 0; ss < nSyst; ss++){
					KD = Txy / 2000.; // in 2 mm
					double TxySystwgt = 1;
					if (ss == (kSystNominal + iSystTxy)) TxySystwgt = txy_syst->Eval(KD);
					if (ss == (kSystNominal - iSystTxy)) TxySystwgt = 1. / txy_syst->Eval(KD);
					KD = tanh(KD);
					if (KD >= 1.) KD = 1. - 1.0e-10;
					if (KD <= -1.) KD = -1. + 1.0e-10;

					templateWeight[ss] = TxySystwgt;
					ctauWeight = exp(-GenCTau*effectiveInvCTau);
					templateWeight[ss] *= ctauWeight*xsec_ctau_wgt;
					templateWeight[ss] *= integral_rescale[ss];
				}
				mytree->Fill();
			}
		}
		cout << "Filled the tree" << endl;
		for (int ss = 0; ss < nSyst; ss++){
			hfill_conditional = (TH2F*)hfill[ss]->Clone(Form("%s_Conditional", hfill[ss]->GetName()));
			hfill_conditional->SetTitle(Form("%s Conditional on tanh(T_{xy})", hfill[ss]->GetTitle()));
			for (int biny = 1; biny <= hfill_conditional->GetNbinsY(); biny++){
				double integralX = hfill_conditional->Integral(0, hfill_conditional->GetNbinsX() + 1, biny, biny);
				if (integralX != 0){
					for (int binx = 1; binx <= hfill_conditional->GetNbinsX(); binx++) hfill_conditional->SetBinContent(binx, biny, (hfill_conditional->GetBinContent(binx, biny) / integralX));
				}
			}

			hfill_projX = (TH1F*) hfill[ss]->ProjectionX();
			hfill_projX->SetTitle(Form("%s Projection on tanh(T_{xy})", hfill[ss]->GetTitle()));

			foutput->WriteTObject(hfill[ss]);
			foutput->WriteTObject(hfill_conditional);
			foutput->WriteTObject(hfill_projX);
			delete hfill_projX;
			delete hfill_conditional;
			delete hfill[ss];
		}
		foutput->WriteTObject(mytree);
		delete mytree;
	}
	for(int f=0;f<kNumSamples;f++){if(tree[f]!=0) delete tree[f];}

	foutput->Close();
	cout << "Closed " << coutput << endl;
	finput_syst->Close();
}


bool testLepId_GenLevel(int channel, int Lep1ID, int Lep2ID, int Lep3ID, int Lep4ID){
	bool result = true;
	if( channel==0
		&& !(
				abs(Lep1ID)==abs(Lep2ID)
				&& abs(Lep1ID)==abs(Lep3ID)
				&& abs(Lep1ID)==abs(Lep4ID)
				&& abs(Lep1ID)==13
			)
		) result = false;
	else if( channel==1
		&& !(
				abs(Lep1ID)==abs(Lep2ID)
				&& abs(Lep1ID)==abs(Lep3ID)
				&& abs(Lep1ID)==abs(Lep4ID)
				&& abs(Lep1ID)==11
			)
		) result = false;
	else if( channel==2
		&& !(
				(abs(Lep1ID)==11 && abs(Lep3ID)==13)
				|| (abs(Lep1ID)==13 && abs(Lep3ID)==11)
			)
		) result = false;
	return result;
}

void produce_KDShapeVariation(int folder, int erg_tev){
	char TREE_NAME[]="SelectedTree";
	TString OUTPUT_NAME = "LifetimeKD_ShapeVariation_DataQQBZZ_";
	OUTPUT_NAME.Append(user_folder[folder]);
	OUTPUT_NAME.Append(".root");
	TString erg_dir;
	erg_dir.Form("LHC_%iTeV",erg_tev);
	TString comstring;
	comstring.Form("%iTeV",erg_tev);
	int EnergyIndex=0;if(erg_tev==8) EnergyIndex=1;

	float yield_offshell_zx[2][3] = {
		{ 0.1000, 0.4000, 0.34 },
		{ 0.5500, 1.7800, 1.38 }
	};
	float yield_offshell_ggzz[2][3]={
		{ 1.7132, 1.1836, 2.8654 },
		{ 8.9407, 6.2543, 14.8431 }
	};
	float yield_80_100_zx[2][3] = {
		{ 0.0103, 0.1174, 0.1465 },
		{ 0.0549, 0.5219, 0.5916 }
	};
	float yield_80_100_ggzz[2][3]={
		{ 0.0222, 0.0106, 0.0113 },
		{ 0.2441, 0.0918, 0.0649 }
	};
	float ggzz_offshell_sum=0,ggzz_80_100_sum=0;
	float zx_offshell_sum=0,zx_80_100_sum=0;
	float luminosity[2] = { 5.051, 19.712 };

	float MC_weight;
	float MC_weight_QQZZEWK=1;
	float MC_weight_Kfactor=1;
	float MC_weight_QQBGGProper[4];
	float ZZMass,ZZPt,ZZEta,ZZPhi;
	float KalmanCandVtx_x,KalmanCandVtx_y,KalmanCandVtx_z;
	float OfflinePrimaryVtx_x,OfflinePrimaryVtx_y,OfflinePrimaryVtx_z;

	float Txy,Tz,KD;
	float CandVtx_x,CandVtx_y,CandVtx_z;

	TString cinput_common = user_dir_hep + erg_dir + "/";
	TString coutput_common = user_dir_hep + "Analysis/Templates/" + comstring + "/";
	TString coutput = coutput_common + OUTPUT_NAME;
	TFile* foutput = new TFile(coutput,"recreate");

	const int ntrees=4;
	TString cinput_qqzz_common = cinput_common + user_folder[folder] + "/";
	TString cinput_ggzz_common = cinput_qqzz_common;
	TString cinput_zx_common = cinput_common;

	TChain* tqqzz = new TChain(TREE_NAME);
	for (int smp = kGGSamples; smp < kQQBZZSamples; smp++){
		TString cinput_qqzz = cinput_qqzz_common + sample_FullSim[smp] + "_Reprocessed.root";
		tqqzz->Add(cinput_qqzz);
	}
	TChain* tggzz = new TChain(TREE_NAME);
	for (int smp = kGGOLDSamples; smp < kGGMCFMSamples; smp++){
		TString cinput_ggzz = cinput_qqzz_common + sample_FullSim[smp] + "_Reprocessed.root";
		tggzz->Add(cinput_ggzz);
	}
	TChain* tzx = new TChain(TREE_NAME);
	for (int ee = 7; ee < 9; ee++){
		for (int f = 0; f < 3; f++){
			TString cinput_zx = user_dir_hep + "LHC_";
			cinput_zx.Append(Form("%iTeV/", ee));
			cinput_zx = cinput_zx + user_folder[folder] + "/" + sample_FullSim[kAllSamples - 1] + "_Reprocessed.root";
			tzx->Add(cinput_zx);
		}
	}

	TString cinput_data = cinput_common + user_folder[3] + "/" + data_files[folder] + "_Reprocessed.root";
	TChain* tdata = new TChain(TREE_NAME);
	tdata->Add(cinput_data);
	
	float systZZMass_range[3][2] = {
		{ 70, 105.6 }, { 140.6, 170 }, { 170, 800 }
	};
	const int nbins_KD = 12;
	TH2F* hdata_full = new TH2F("hdata_full",Form("%i TeV %s Data",erg_tev,user_folder[folder]),3,0,3,nbins_KD,-1,1);
	hdata_full->GetYaxis()->SetTitle("tanh(T_{xy})");
	hdata_full->GetXaxis()->SetTitle("m_{4l} (GeV)");
	hdata_full->GetXaxis()->SetBinLabel(1,"70 - 105.6");
	hdata_full->GetXaxis()->SetBinLabel(2,"140.6 - 170");
	hdata_full->GetXaxis()->SetBinLabel(3,"170 - 800");
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

	TChain* tc[ntrees] = { tdata, tqqzz, tggzz, tzx };
	for (int tt = 0; tt < ntrees; tt++){
		cout << "Tree " << tt << endl;
		if (tt>0){
			if (tt != 3) tc[tt]->SetBranchAddress("MC_weight", &MC_weight);
			else tc[tt]->SetBranchAddress("ZXfake_weightProper", &MC_weight);
			if (tt == 2 || tt == 1) tc[tt]->SetBranchAddress("MC_weight_Kfactor", &MC_weight_Kfactor);
			if (tt == 1){
				tc[tt]->SetBranchAddress("MC_weight_QQZZEWK", &MC_weight_QQZZEWK);
				tc[tt]->SetBranchAddress("MC_weight_QQBGGProper", &MC_weight_QQBGGProper);
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
	}
	for (int tt = 0; tt < ntrees; tt++){
		for (int ev = 0; ev < tc[tt]->GetEntries(); ev++){
			MC_weight = 1;
			MC_weight_QQZZEWK = 1;
			MC_weight_Kfactor = 1;

			tc[tt]->GetEntry(ev);

			CandVtx_x = KalmanCandVtx_x - OfflinePrimaryVtx_x;
			CandVtx_y = KalmanCandVtx_y - OfflinePrimaryVtx_y;
			CandVtx_z = KalmanCandVtx_z - OfflinePrimaryVtx_z;

			Txy = (CandVtx_x*cos(ZZPhi) + CandVtx_y*sin(ZZPhi))*ZZMass / ZZPt;
//			Tz = CandVtx_z*ZZMass/( ZZPt*sinh(ZZEta) );

			KD = Txy;
			KD = tanh(KD*5.); // Different from gen. level since full sim is in cm while gen. level is in microns
			if (KD >= 1.) KD = 1. - 1.0e-10;
			if (KD <= -1.) KD = -1. + 1.0e-10;

			float wgt = 1;

			int biny = hfull[tt]->GetYaxis()->FindBin(KD);
			for (int binx = 0; binx < 3; binx++){
				if (tt == 3 && binx == 0 && ZZMass >= systZZMass_range[binx][0] + 30 && ZZMass < systZZMass_range[binx][1] + 30){ // ZX special fill
					wgt = MC_weight;
					hfull[tt]->AddBinContent(hfull[tt]->GetBin(binx + 1, biny), wgt);
					if(ZZMass<systZZMass_range[binx][0] + 50) zx_80_100_sum += wgt;
				}
				else if (tt == 3 && binx == 0) continue;
				if (tt == 2 && binx == 0) continue;
				if (tt == 1 && binx == 0 && ZZMass >= systZZMass_range[binx][0] && ZZMass < systZZMass_range[binx][1]){ // gg bkg special fill
					wgt = MC_weight*MC_weight_Kfactor*MC_weight_QQBGGProper[0];
					hfull[2]->AddBinContent(hfull[tt]->GetBin(binx + 1, biny), wgt);
					if(ZZMass>=80 && ZZMass<100) ggzz_80_100_sum += wgt;
				}
				if (ZZMass >= systZZMass_range[binx][1] || ZZMass < systZZMass_range[binx][0]) continue;
				wgt = MC_weight;
				if (tt == 1) wgt *= MC_weight_QQZZEWK;
				if (tt == 2) wgt *= MC_weight_Kfactor;
				hfull[tt]->AddBinContent(hfull[tt]->GetBin(binx + 1, biny), wgt);
				if(tt==2 && ZZMass>=220) ggzz_offshell_sum += wgt;
				if(tt==3 && ZZMass>=220) zx_offshell_sum += wgt;
			}
		}
	}

	for (int biny = 1; biny <= hggzz_full->GetNbinsY(); biny++){
		hggzz_full->SetBinContent(1, biny, (hggzz_full->GetBinContent(1, biny))*(yield_80_100_ggzz[EnergyIndex][folder] / ggzz_80_100_sum / luminosity[EnergyIndex]));
		hggzz_full->SetBinContent(3, biny, (hggzz_full->GetBinContent(3, biny))*(yield_offshell_ggzz[EnergyIndex][folder] / ggzz_offshell_sum / luminosity[EnergyIndex]));
		hggzz_full->SetBinContent(2, biny, (hggzz_full->GetBinContent(2, biny))*(yield_offshell_ggzz[EnergyIndex][folder] / ggzz_offshell_sum / luminosity[EnergyIndex]));
	}
	for (int biny = 1; biny <= hzx_full->GetNbinsY(); biny++){
		hzx_full->SetBinContent(1, biny, (hzx_full->GetBinContent(1, biny))*(yield_80_100_zx[EnergyIndex][folder] / zx_80_100_sum / luminosity[EnergyIndex]));
		hzx_full->SetBinContent(3, biny, (hzx_full->GetBinContent(3, biny))*(yield_offshell_zx[EnergyIndex][folder] / zx_offshell_sum / luminosity[EnergyIndex]));
		hzx_full->SetBinContent(2, biny, (hzx_full->GetBinContent(2, biny))*(yield_offshell_zx[EnergyIndex][folder] / zx_offshell_sum / luminosity[EnergyIndex]));
	}

	for (int tt = 1; tt < ntrees; tt++) hfull[ntrees]->Add(hfull[tt]);
	for (int tt = 0; tt < ntrees+1; tt++){
		TH1F* hProjKD = (TH1F*) hfull[tt]->ProjectionY();
		hProjKD->SetTitle(Form("%s Projection",hfull[tt]->GetName()));
//		hProjKD->Scale(1./hProjKD->Integral(0,hProjKD->GetNbinsX()));
		hproj[tt] = hProjKD;

		for (int binx = 1; binx <= hfull[tt]->GetNbinsX(); binx++){
			double integralY = hfull[tt]->Integral(binx, binx, 0, hfull[tt]->GetNbinsY() + 1);
			if (integralY != 0){
//				for (int biny = 1; biny <= hfull[tt]->GetNbinsY(); biny++) hfull[tt]->SetBinContent(binx, biny, (hfull[tt]->GetBinContent(binx, biny) / integralY));
			}
		}

		foutput->WriteTObject(hfull[tt]);
		foutput->WriteTObject(hproj[tt]);
		delete hfull[tt];
		if(tt<ntrees) delete tc[tt];
	}

	foutput->Close();
}

void compare_KDShapeVariation(){
	TString OUTPUT_NAME = "LifetimeKD_ShapeVariation_DataQQBZZ_Comparison.root";
	TString coutput_common = user_dir_hep + "Analysis/Plots/";
	TString coutput = coutput_common + OUTPUT_NAME;
	TFile* foutput = new TFile(coutput, "recreate");
	const int ntrees = 2;
	char* chProj[ntrees]={"hdata_full_py", "hmc_full_py"};
	TH1F* hProj[ntrees][3] = { { 0 } };
	TGraphAsymmErrors* tgdata[3];
	TGraphAsymmErrors* tgratio[3];
	double max_plot=0;
	TFile* finput[3][2] = { { 0 } };

	for (int folder = 0; folder < 3; folder++){
		for (int erg_tev = 7; erg_tev < 9; erg_tev++){
			TString erg_dir;
			erg_dir.Form("LHC_%iTeV", erg_tev);
			TString comstring;
			comstring.Form("%iTeV", erg_tev);
			int EnergyIndex = 0; if (erg_tev == 8) EnergyIndex = 1;
			TString cinput_common = user_dir_hep + "Analysis/Templates/" + comstring + "/";

			TString INPUT_NAME = "LifetimeKD_ShapeVariation_DataQQBZZ_";
			INPUT_NAME.Append(user_folder[folder]);
			INPUT_NAME.Append(".root");
			cout << INPUT_NAME << endl;
			TString cinput = cinput_common + INPUT_NAME;
			finput[folder][EnergyIndex] = new TFile(cinput, "read");

			if (finput[folder][EnergyIndex]==0 || finput[folder][EnergyIndex]->IsZombie()) continue;
			for (int hh = 0; hh < ntrees; hh++){
				TH1F* htemp = (TH1F*) finput[folder][EnergyIndex]->Get(chProj[hh]);
				if (hProj[hh][folder] == 0){
					hProj[hh][folder] = (TH1F*)htemp->Clone(Form("acc_%s", htemp->GetName()));
				}
				else hProj[hh][folder]->Add(htemp);
				delete htemp;
			}
		}

		hProj[1][folder]->Scale(1./hProj[1][folder]->Integral());

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
		const double quant = (1.0-0.6827)/2.0;
		double integral_data=hProj[0][folder]->Integral();

		for(int bin=1;bin<=hProj[0][folder]->GetNbinsX();bin++){
			double bincenter = hProj[0][folder]->GetBinCenter(bin);
			double bincontent = hProj[0][folder]->GetBinContent(bin);

			if(bincontent>0){
				xx_data[ndata] = bincenter;
				yy_data[ndata] = bincontent/integral_data;
				xu_data[ndata] = 0;
				xd_data[ndata] = 0;
				yu_data[ndata] = (ROOT::Math::chisquared_quantile_c(quant,2*(bincontent+1))/2.-bincontent)/integral_data;
				yd_data[ndata] = ((bincontent==0)?0:(bincontent-ROOT::Math::chisquared_quantile_c(1-quant,2*bincontent)/2.))/integral_data;

				rr_data[ndata] = yy_data[ndata] / hProj[1][folder]->GetBinContent(bin);
				ru_data[ndata] = yu_data[ndata] / hProj[1][folder]->GetBinContent(bin);
				rd_data[ndata] = yd_data[ndata] / hProj[1][folder]->GetBinContent(bin);

				double high_data = yy_data[ndata]+yu_data[ndata];
				if (high_data > max_plot) max_plot = high_data;
				ndata++;
			}
		}
		tgdata[folder] = new TGraphAsymmErrors(ndata,xx_data,yy_data,xd_data,xu_data,yd_data,yu_data);
		tgratio[folder] = new TGraphAsymmErrors(ndata,xx_data,rr_data,xd_data,xu_data,rd_data,ru_data);
		tgdata[folder]->SetName(Form("tgdata_%s",user_folder[folder]));
		tgratio[folder]->SetName(Form("tgratio_%s",user_folder[folder]));
		tgdata[folder]->SetMarkerSize(1.2);
		tgdata[folder]->SetMarkerStyle(20);
		tgdata[folder]->SetMarkerColor(kBlack);
		tgdata[folder]->SetLineColor(kBlack);
		tgdata[folder]->SetLineWidth(1);

		hProj[0][folder]->Scale(1./hProj[0][folder]->Integral());
	}

	gStyle->SetTitleFont(62, "t");
	gROOT->SetStyle(gStyle->GetName());
	gROOT->ForceStyle();

	TPaveText* pt = new TPaveText(0.15,0.93,0.85,1,"brNDC");
	pt->SetBorderSize(0);
	pt->SetFillStyle(0);
	pt->SetTextAlign(12);
	pt->SetTextFont(42);
	pt->SetTextSize(0.045);
	TText* text = pt->AddText(0.025,0.45,"#font[61]{CMS}");
	text->SetTextSize(0.044);
	text = pt->AddText(0.165, 0.42, "#font[52]{Preliminary}");
	text->SetTextSize(0.0315);
	text = pt->AddText(0.837,0.45,"#font[42]{19.7 fb^{-1} (8 TeV)}");
//	text = pt->AddText(0.537,0.45,"#font[42]{19.7 fb^{-1} (8 TeV) + 5.1 fb^{-1} (7 TeV)}");
	text->SetTextSize(0.0315);

	tgdata[2]->SetLineColor(kOrange + 10);
	tgdata[2]->SetMarkerColor(kOrange + 10);
	hProj[1][2]->SetLineColor(kOrange + 10);
	hProj[1][2]->SetLineWidth(2);
	hProj[1][2]->SetLineStyle(1);
	tgdata[1]->SetLineColor(kBlue);
	tgdata[1]->SetMarkerColor(kBlue);
	hProj[1][1]->SetLineColor(kBlue);
	hProj[1][1]->SetLineWidth(2);
	hProj[1][1]->SetLineStyle(1);
	tgdata[0]->SetLineColor(TColor::GetColor("#669966"));
	tgdata[0]->SetMarkerColor(TColor::GetColor("#669966"));
	hProj[1][0]->SetLineColor(TColor::GetColor("#669966"));
	hProj[1][0]->SetLineWidth(2);
	hProj[1][0]->SetLineStyle(1);
	for (int folder = 0; folder < 3; folder++){
		for (int hh = 0; hh < ntrees; hh++){
			hProj[hh][folder]->SetTitle("");
			tgdata[folder]->SetTitle("");
		}
	}

	foutput->cd();
	TString canvasname_2D = "cCompare_DataMC_AllChannels_VtxSyst_tanh";
	TCanvas* c2D = new TCanvas(canvasname_2D,"",8,30,800,800);
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

	TLegend *l2D = new TLegend(0.20,0.57,0.58,0.90);
	l2D->SetBorderSize(0);
	l2D->SetTextFont(42);
	l2D->SetTextSize(0.04);
	l2D->SetLineColor(1);
	l2D->SetLineStyle(1);
	l2D->SetLineWidth(1);
	l2D->SetFillColor(0);
	l2D->SetFillStyle(0);

	char yTitle[] = "Sideband Rate";
	char xTitle[] = "tanh(T_{xy})";

	hProj[1][0]->GetXaxis()->SetTitle(xTitle);
	hProj[1][0]->GetYaxis()->SetTitle(yTitle);
	hProj[1][0]->GetYaxis()->SetRangeUser(0,max_plot*1.25);
	hProj[1][0]->GetXaxis()->SetNdivisions(505);
	hProj[1][0]->GetXaxis()->SetLabelFont(42);
	hProj[1][0]->GetXaxis()->SetLabelOffset(0.007);
	hProj[1][0]->GetXaxis()->SetLabelSize(0.04);
	hProj[1][0]->GetXaxis()->SetTitleSize(0.06);
	hProj[1][0]->GetXaxis()->SetTitleOffset(0.9);
	hProj[1][0]->GetXaxis()->SetTitleFont(42);
	hProj[1][0]->GetYaxis()->SetNdivisions(505);
	hProj[1][0]->GetYaxis()->SetLabelFont(42);
	hProj[1][0]->GetYaxis()->SetLabelOffset(0.007);
	hProj[1][0]->GetYaxis()->SetLabelSize(0.04);
	hProj[1][0]->GetYaxis()->SetTitleSize(0.06);
	hProj[1][0]->GetYaxis()->SetTitleOffset(1.1);
	hProj[1][0]->GetYaxis()->SetTitleFont(42);

	l2D->AddEntry(tgdata[2],"Observed 2e2#mu","ep");
	l2D->AddEntry(tgdata[1],"Observed 4e","ep");
	l2D->AddEntry(tgdata[0],"Observed 4#mu","ep");
	l2D->AddEntry(hProj[1][2],"Bkg. 2e2#mu","l");
	l2D->AddEntry(hProj[1][1],"Bkg. 4e","l");
	l2D->AddEntry(hProj[1][0],"Bkg. 4#mu","l");

	hProj[1][0]->Draw("hist");
	hProj[1][1]->Draw("same");
	hProj[1][2]->Draw("same");
	tgdata[0]->Draw("e1psame");
	tgdata[1]->Draw("e1psame");
	tgdata[2]->Draw("e1psame");
	l2D->Draw("same");
	pt->Draw();

	c2D->RedrawAxis();
	c2D->Modified();
	c2D->Update();
	foutput->WriteTObject(c2D);

	TString canvasDir = coutput_common;
	canvasname_2D.Prepend(canvasDir);
	TString canvasname_2D_pdf=canvasname_2D;
	TString canvasname_2D_eps=canvasname_2D;
	TString canvasname_2D_png=canvasname_2D;
	TString canvasname_2D_root=canvasname_2D;
	TString canvasname_2D_c=canvasname_2D;
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
//	delete text;
	delete pt;

	for (int folder = 0; folder < 3; folder++){
		foutput->WriteTObject(tgdata[folder]);
		foutput->WriteTObject(tgratio[folder]);
		for (int hh = 0; hh < ntrees; hh++) delete hProj[hh][folder];
		delete tgratio[folder];
		delete tgdata[folder];
	}
	for (int folder = 0; folder < 3; folder++){
		for (int e = 0; e < 2; e++){
			if (finput[folder][e] == 0 || finput[folder][e]->IsZombie()) continue;
			else finput[folder][e]->Close();
		}
	}
	foutput->Close();
}


void produce_KDShapeVariation_Untransformed(int folder, int erg_tev){
	char TREE_NAME[]="SelectedTree";
	TString OUTPUT_NAME = "LifetimeKD_ShapeVariation_DataQQBZZ_Untransformed_";
	OUTPUT_NAME.Append(user_folder[folder]);
	OUTPUT_NAME.Append(".root");
	TString erg_dir;
	erg_dir.Form("LHC_%iTeV",erg_tev);
	TString comstring;
	comstring.Form("%iTeV",erg_tev);
	int EnergyIndex=0;if(erg_tev==8) EnergyIndex=1;

	float yield_offshell_zx[2][3] = {
		{ 0.1000, 0.4000, 0.34 },
		{ 0.5500, 1.7800, 1.38 }
	};
	float yield_offshell_ggzz[2][3]={
		{ 1.7132, 1.1836, 2.8654 },
		{ 8.9407, 6.2543, 14.8431 }
	};
	float yield_80_100_zx[2][3] = {
		{ 0.0103, 0.1174, 0.1465 },
		{ 0.0549, 0.5219, 0.5916 }
	};
	float yield_80_100_ggzz[2][3]={
		{ 0.0222, 0.0106, 0.0113 },
		{ 0.2441, 0.0918, 0.0649 }
	};
	float ggzz_offshell_sum=0,ggzz_80_100_sum=0;
	float zx_offshell_sum=0,zx_80_100_sum=0;
	float luminosity[2] = { 5.051, 19.712 };

	float MC_weight;
	float MC_weight_QQZZEWK=1;
	float MC_weight_Kfactor=1;
	float MC_weight_QQBGGProper[4];
	float ZZMass,ZZPt,ZZEta,ZZPhi;
	float KalmanCandVtx_x,KalmanCandVtx_y,KalmanCandVtx_z;
	float OfflinePrimaryVtx_x,OfflinePrimaryVtx_y,OfflinePrimaryVtx_z;

	float Txy,Tz,KD;
	float CandVtx_x,CandVtx_y,CandVtx_z;

	TString cinput_common = user_dir_hep + erg_dir + "/";
	TString coutput_common = user_dir_hep + "Analysis/Templates/" + comstring + "/";
	TString coutput = coutput_common + OUTPUT_NAME;
	TFile* foutput = new TFile(coutput,"recreate");

	const int ntrees=4;
	TString cinput_qqzz_common = cinput_common + user_folder[folder] + "/";
	TString cinput_ggzz_common = cinput_qqzz_common;
	TString cinput_zx_common = cinput_common;

	TChain* tqqzz = new TChain(TREE_NAME);
	for (int smp = kGGSamples; smp < kQQBZZSamples; smp++){
		TString cinput_qqzz = cinput_qqzz_common + sample_FullSim[smp] + "_Reprocessed.root";
		tqqzz->Add(cinput_qqzz);
	}
	TChain* tggzz = new TChain(TREE_NAME);
	for (int smp = kGGOLDSamples; smp < kGGMCFMSamples; smp++){
		TString cinput_ggzz = cinput_qqzz_common + sample_FullSim[smp] + "_Reprocessed.root";
		tggzz->Add(cinput_ggzz);
	}
	TChain* tzx = new TChain(TREE_NAME);
	for (int ee = 7; ee < 9; ee++){
		for (int f = 0; f < 3; f++){
			TString cinput_zx = user_dir_hep + "LHC_";
			cinput_zx.Append(Form("%iTeV/", ee));
			cinput_zx = cinput_zx + user_folder[folder] + "/" + sample_FullSim[kAllSamples - 1] + "_Reprocessed.root";
			tzx->Add(cinput_zx);
		}
	}

	TString cinput_data = cinput_common + user_folder[3] + "/" + data_files[folder] + "_Reprocessed.root";
	TChain* tdata = new TChain(TREE_NAME);
	tdata->Add(cinput_data);
	
	float systZZMass_range[3][2] = {
		{ 70, 105.6 }, { 140.6, 170 }, { 170, 800 }
	};
//	const int nbins_KD = 200;
	const int nbins_KD = 8;
	double KD_limits[2] = { -20., 20. };
	double KD_bins[nbins_KD+1] = { -20, -1, -0.4, -0.15, 0, 0.15, 1., 3., 20 };
	double xbins[4] = { 0, 1, 2, 3 };
//	TH2F* hdata_full = new TH2F("hdata_full",Form("%i TeV %s Data",erg_tev,user_folder[folder]),3,0,3,nbins_KD,KD_limits[0],KD_limits[1]);
	TH2F* hdata_full = new TH2F("hdata_full",Form("%i TeV %s Data",erg_tev,user_folder[folder]),3,xbins,nbins_KD,KD_bins);
	hdata_full->GetYaxis()->SetTitle("T_{xy} (2 mm)");
	hdata_full->GetXaxis()->SetTitle("m_{4l} (GeV)");
	hdata_full->GetXaxis()->SetBinLabel(1,"70 - 105.6");
	hdata_full->GetXaxis()->SetBinLabel(2,"140.6 - 170");
	hdata_full->GetXaxis()->SetBinLabel(3,"170 - 800");
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

	TChain* tc[ntrees] = { tdata, tqqzz, tggzz, tzx };
	for (int tt = 0; tt < ntrees; tt++){
		cout << "Tree " << tt << endl;
		if (tt>0){
			if (tt != 3) tc[tt]->SetBranchAddress("MC_weight", &MC_weight);
			else tc[tt]->SetBranchAddress("ZXfake_weightProper", &MC_weight);
			if (tt == 2 || tt == 1) tc[tt]->SetBranchAddress("MC_weight_Kfactor", &MC_weight_Kfactor);
			if (tt == 1){
				tc[tt]->SetBranchAddress("MC_weight_QQZZEWK", &MC_weight_QQZZEWK);
				tc[tt]->SetBranchAddress("MC_weight_QQBGGProper", &MC_weight_QQBGGProper);
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
	}
	for (int tt = 0; tt < ntrees; tt++){
		for (int ev = 0; ev < tc[tt]->GetEntries(); ev++){
			MC_weight = 1;
			MC_weight_QQZZEWK = 1;
			MC_weight_Kfactor = 1;

			tc[tt]->GetEntry(ev);

			CandVtx_x = KalmanCandVtx_x - OfflinePrimaryVtx_x;
			CandVtx_y = KalmanCandVtx_y - OfflinePrimaryVtx_y;
			CandVtx_z = KalmanCandVtx_z - OfflinePrimaryVtx_z;

			Txy = (CandVtx_x*cos(ZZPhi) + CandVtx_y*sin(ZZPhi))*ZZMass / ZZPt;

			KD = Txy;
			KD = KD*5.; // Different from gen. level since full sim is in cm while gen. level is in microns

			float wgt = 1;

			int biny = hfull[tt]->GetYaxis()->FindBin(KD);
			for (int binx = 0; binx < 3; binx++){
				if (tt == 3 && binx == 0 && ZZMass >= systZZMass_range[binx][0] + 30 && ZZMass < systZZMass_range[binx][1] + 30){ // ZX special fill
					wgt = MC_weight;
					hfull[tt]->AddBinContent(hfull[tt]->GetBin(binx + 1, biny), wgt);
					if(ZZMass<systZZMass_range[binx][0] + 50) zx_80_100_sum += wgt;
				}
				else if (tt == 3 && binx == 0) continue;
				if (tt == 2 && binx == 0) continue;
				if (tt == 1 && binx == 0 && ZZMass >= systZZMass_range[binx][0] && ZZMass < systZZMass_range[binx][1]){ // gg bkg special fill
					wgt = MC_weight*MC_weight_Kfactor*MC_weight_QQBGGProper[0];
					hfull[2]->AddBinContent(hfull[tt]->GetBin(binx + 1, biny), wgt);
					if(ZZMass>=80 && ZZMass<100) ggzz_80_100_sum += wgt;
				}
				if (ZZMass >= systZZMass_range[binx][1] || ZZMass < systZZMass_range[binx][0]) continue;
				wgt = MC_weight;
				if (tt == 1) wgt *= MC_weight_QQZZEWK;
				if (tt == 2) wgt *= MC_weight_Kfactor;
				hfull[tt]->AddBinContent(hfull[tt]->GetBin(binx + 1, biny), wgt);
				if(tt==2 && ZZMass>=220) ggzz_offshell_sum += wgt;
				if(tt==3 && ZZMass>=220) zx_offshell_sum += wgt;
			}
		}
	}

	for (int biny = 1; biny <= hggzz_full->GetNbinsY(); biny++){
		hggzz_full->SetBinContent(1, biny, (hggzz_full->GetBinContent(1, biny))*(yield_80_100_ggzz[EnergyIndex][folder] / ggzz_80_100_sum / luminosity[EnergyIndex]));
		hggzz_full->SetBinContent(3, biny, (hggzz_full->GetBinContent(3, biny))*(yield_offshell_ggzz[EnergyIndex][folder] / ggzz_offshell_sum / luminosity[EnergyIndex]));
		hggzz_full->SetBinContent(2, biny, (hggzz_full->GetBinContent(2, biny))*(yield_offshell_ggzz[EnergyIndex][folder] / ggzz_offshell_sum / luminosity[EnergyIndex]));
	}
	for (int biny = 1; biny <= hzx_full->GetNbinsY(); biny++){
		hzx_full->SetBinContent(1, biny, (hzx_full->GetBinContent(1, biny))*(yield_80_100_zx[EnergyIndex][folder] / zx_80_100_sum / luminosity[EnergyIndex]));
		hzx_full->SetBinContent(3, biny, (hzx_full->GetBinContent(3, biny))*(yield_offshell_zx[EnergyIndex][folder] / zx_offshell_sum / luminosity[EnergyIndex]));
		hzx_full->SetBinContent(2, biny, (hzx_full->GetBinContent(2, biny))*(yield_offshell_zx[EnergyIndex][folder] / zx_offshell_sum / luminosity[EnergyIndex]));
	}

	for (int tt = 1; tt < ntrees; tt++) hfull[ntrees]->Add(hfull[tt]);
	for (int tt = 0; tt < ntrees+1; tt++){
		TH1F* hProjKD = (TH1F*) hfull[tt]->ProjectionY();
		hProjKD->SetTitle(Form("%s Projection",hfull[tt]->GetName()));
		hproj[tt] = hProjKD;

		for (int binx = 1; binx <= hfull[tt]->GetNbinsX(); binx++){
			double integralY = hfull[tt]->Integral(binx, binx, 0, hfull[tt]->GetNbinsY() + 1);
			if (integralY != 0){
//				for (int biny = 1; biny <= hfull[tt]->GetNbinsY(); biny++) hfull[tt]->SetBinContent(binx, biny, (hfull[tt]->GetBinContent(binx, biny) / integralY));
			}
		}

		foutput->WriteTObject(hfull[tt]);
		foutput->WriteTObject(hproj[tt]);
		delete hfull[tt];
		if(tt<ntrees) delete tc[tt];
	}

	foutput->Close();
}


void compare_KDShapeVariation_Untransformed(){
	TString OUTPUT_NAME = "LifetimeKD_ShapeVariation_DataQQBZZ_Untransformed_Comparison.root";
	TString coutput_common = user_dir_hep + "Analysis/Plots/";
	TString coutput = coutput_common + OUTPUT_NAME;
	TFile* foutput = new TFile(coutput, "recreate");
	const int ntrees = 2;
	char* chProj[ntrees]={"hdata_full_py", "hmc_full_py"};
	TH1F* hProj[ntrees][3] = { { 0 } };
	TGraphAsymmErrors* tgdata[3];
	TGraphAsymmErrors* tgratio[3];
	double max_plot=0;
	TFile* finput[3][2] = { { 0 } };

	for (int folder = 0; folder < 3; folder++){
		for (int erg_tev = 7; erg_tev < 9; erg_tev++){
			TString erg_dir;
			erg_dir.Form("LHC_%iTeV", erg_tev);
			TString comstring;
			comstring.Form("%iTeV", erg_tev);
			int EnergyIndex = 0; if (erg_tev == 8) EnergyIndex = 1;
			TString cinput_common = user_dir_hep + "Analysis/Templates/" + comstring + "/";

			TString INPUT_NAME = "LifetimeKD_ShapeVariation_DataQQBZZ_Untransformed_";
			INPUT_NAME.Append(user_folder[folder]);
			INPUT_NAME.Append(".root");
			cout << INPUT_NAME << endl;
			TString cinput = cinput_common + INPUT_NAME;
			finput[folder][EnergyIndex] = new TFile(cinput, "read");

			if (finput[folder][EnergyIndex]==0 || finput[folder][EnergyIndex]->IsZombie()) continue;
			for (int hh = 0; hh < ntrees; hh++){
				TH1F* htemp = (TH1F*) finput[folder][EnergyIndex]->Get(chProj[hh]);
				if (hProj[hh][folder] == 0){
					hProj[hh][folder] = (TH1F*)htemp->Clone(Form("acc_%s", htemp->GetName()));
				}
				else hProj[hh][folder]->Add(htemp);
				delete htemp;
			}
		}

		hProj[1][folder]->Scale(1./hProj[1][folder]->Integral());

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
		const double quant = (1.0-0.6827)/2.0;
		double integral_data=hProj[0][folder]->Integral();

		for(int bin=1;bin<=hProj[0][folder]->GetNbinsX();bin++){
			double bincenter = hProj[0][folder]->GetBinCenter(bin);
			double bincontent = hProj[0][folder]->GetBinContent(bin);

			if(bincontent>0){
				xx_data[ndata] = bincenter;
				yy_data[ndata] = bincontent/integral_data;
				xu_data[ndata] = 0;
				xd_data[ndata] = 0;
				yu_data[ndata] = (ROOT::Math::chisquared_quantile_c(quant,2*(bincontent+1))/2.-bincontent)/integral_data;
				yd_data[ndata] = ((bincontent==0)?0:(bincontent-ROOT::Math::chisquared_quantile_c(1-quant,2*bincontent)/2.))/integral_data;

				rr_data[ndata] = yy_data[ndata] / hProj[1][folder]->GetBinContent(bin);
				ru_data[ndata] = yu_data[ndata] / hProj[1][folder]->GetBinContent(bin);
				rd_data[ndata] = yd_data[ndata] / hProj[1][folder]->GetBinContent(bin);

				double high_data = yy_data[ndata]+yu_data[ndata];
				if (high_data > max_plot) max_plot = high_data;
				ndata++;
			}
		}
		tgdata[folder] = new TGraphAsymmErrors(ndata,xx_data,yy_data,xd_data,xu_data,yd_data,yu_data);
		tgratio[folder] = new TGraphAsymmErrors(ndata,xx_data,rr_data,xd_data,xu_data,rd_data,ru_data);
		tgdata[folder]->SetName(Form("tgdata_%s",user_folder[folder]));
		tgratio[folder]->SetName(Form("tgratio_%s",user_folder[folder]));
		tgdata[folder]->SetMarkerSize(1.2);
		tgdata[folder]->SetMarkerStyle(20);
		tgdata[folder]->SetMarkerColor(kBlack);
		tgdata[folder]->SetLineColor(kBlack);
		tgdata[folder]->SetLineWidth(1);

		hProj[0][folder]->Scale(1./hProj[0][folder]->Integral());
	}

	gStyle->SetTitleFont(62, "t");
	gROOT->SetStyle(gStyle->GetName());
	gROOT->ForceStyle();

	TPaveText* pt = new TPaveText(0.15,0.93,0.85,1,"brNDC");
	pt->SetBorderSize(0);
	pt->SetFillStyle(0);
	pt->SetTextAlign(12);
	pt->SetTextFont(42);
	pt->SetTextSize(0.045);
	TText* text = pt->AddText(0.025,0.45,"#font[61]{CMS}");
	text->SetTextSize(0.044);
	text = pt->AddText(0.165, 0.42, "#font[52]{Preliminary}");
	text->SetTextSize(0.0315);
	text = pt->AddText(0.837,0.45,"#font[42]{19.7 fb^{-1} (8 TeV)}");
//	text = pt->AddText(0.537,0.45,"#font[42]{19.7 fb^{-1} (8 TeV) + 5.1 fb^{-1} (7 TeV)}");
	text->SetTextSize(0.0315);

	tgdata[2]->SetLineColor(kOrange + 10);
	tgdata[2]->SetMarkerColor(kOrange + 10);
	hProj[1][2]->SetLineColor(kOrange + 10);
	hProj[1][2]->SetLineWidth(2);
	hProj[1][2]->SetLineStyle(1);
	tgdata[1]->SetLineColor(kBlue);
	tgdata[1]->SetMarkerColor(kBlue);
	hProj[1][1]->SetLineColor(kBlue);
	hProj[1][1]->SetLineWidth(2);
	hProj[1][1]->SetLineStyle(1);
	tgdata[0]->SetLineColor(TColor::GetColor("#669966"));
	tgdata[0]->SetMarkerColor(TColor::GetColor("#669966"));
	hProj[1][0]->SetLineColor(TColor::GetColor("#669966"));
	hProj[1][0]->SetLineWidth(2);
	hProj[1][0]->SetLineStyle(1);
	for (int folder = 0; folder < 3; folder++){
		for (int hh = 0; hh < ntrees; hh++){
			hProj[hh][folder]->SetTitle("");
			tgdata[folder]->SetTitle("");
		}
	}

	foutput->cd();
	TString canvasname_2D = "cCompare_DataMC_AllChannels_VtxSyst_Untransformed";
	TCanvas* c2D = new TCanvas(canvasname_2D,"",8,30,800,800);
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

	TLegend *l2D = new TLegend(0.20,0.57,0.58,0.90);
	l2D->SetBorderSize(0);
	l2D->SetTextFont(42);
	l2D->SetTextSize(0.04);
	l2D->SetLineColor(1);
	l2D->SetLineStyle(1);
	l2D->SetLineWidth(1);
	l2D->SetFillColor(0);
	l2D->SetFillStyle(0);

	char yTitle[] = "Sideband Rate";
	char xTitle[] = "T_{xy} (2 mm)";

	hProj[1][0]->GetXaxis()->SetTitle(xTitle);
	hProj[1][0]->GetYaxis()->SetTitle(yTitle);
	hProj[1][0]->GetYaxis()->SetRangeUser(0,max_plot*1.25);
	hProj[1][0]->GetXaxis()->SetNdivisions(505);
	hProj[1][0]->GetXaxis()->SetLabelFont(42);
	hProj[1][0]->GetXaxis()->SetLabelOffset(0.007);
	hProj[1][0]->GetXaxis()->SetLabelSize(0.04);
	hProj[1][0]->GetXaxis()->SetTitleSize(0.06);
	hProj[1][0]->GetXaxis()->SetTitleOffset(0.9);
	hProj[1][0]->GetXaxis()->SetTitleFont(42);
	hProj[1][0]->GetYaxis()->SetNdivisions(505);
	hProj[1][0]->GetYaxis()->SetLabelFont(42);
	hProj[1][0]->GetYaxis()->SetLabelOffset(0.007);
	hProj[1][0]->GetYaxis()->SetLabelSize(0.04);
	hProj[1][0]->GetYaxis()->SetTitleSize(0.06);
	hProj[1][0]->GetYaxis()->SetTitleOffset(1.1);
	hProj[1][0]->GetYaxis()->SetTitleFont(42);

	l2D->AddEntry(tgdata[2],"Observed 2e2#mu","ep");
	l2D->AddEntry(tgdata[1],"Observed 4e","ep");
	l2D->AddEntry(tgdata[0],"Observed 4#mu","ep");
	l2D->AddEntry(hProj[1][2],"Bkg. 2e2#mu","l");
	l2D->AddEntry(hProj[1][1],"Bkg. 4e","l");
	l2D->AddEntry(hProj[1][0],"Bkg. 4#mu","l");

	hProj[1][0]->GetXaxis()->SetRangeUser(-3,3);
	hProj[1][0]->Draw("hist");
	hProj[1][1]->Draw("same");
	hProj[1][2]->Draw("same");
	tgdata[0]->Draw("e1psame");
	tgdata[1]->Draw("e1psame");
	tgdata[2]->Draw("e1psame");
	l2D->Draw("same");
	pt->Draw();

	c2D->RedrawAxis();
	c2D->Modified();
	c2D->Update();
	foutput->WriteTObject(c2D);

	TString canvasDir = coutput_common;
	canvasname_2D.Prepend(canvasDir);
	TString canvasname_2D_pdf=canvasname_2D;
	TString canvasname_2D_eps=canvasname_2D;
	TString canvasname_2D_png=canvasname_2D;
	TString canvasname_2D_root=canvasname_2D;
	TString canvasname_2D_c=canvasname_2D;
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
//	delete text;
	delete pt;
	cout << "Writing the canvas now..." << endl;

	for (int folder = 0; folder < 3; folder++){
		double xmin,ymin,xmax,ymax;
//		if (folder != 2){ tgratio[folder]->RemovePoint(0); tgratio[folder]->RemovePoint(0); }
//		else tgratio[folder]->RemovePoint(tgratio[folder]->GetN()-1);
		tgratio[folder]->GetPoint(0,xmin,ymin);
		tgratio[folder]->GetPoint(tgratio[folder]->GetN()-1,xmax,ymax);
		int parn=20;
		TF1* fratio = new TF1(Form("%s_fit",tgratio[folder]->GetName()),"[0]-[1]*exp(-pow((x-[2])/[3],2))",-100,100);
//		TF1* fratio = new TF1(Form("%s_fit",tgratio[folder]->GetName()),formula_TxySyst,-100,100,7);
/*		fratio->SetParameters(1.,-0.06,0.31,0.29,(double) parn,xmin,xmax);
		if(folder==1 || folder==0) fratio->SetParameters(1,0,0.3,0.29,(double) parn,xmin,xmax);
		fratio->SetParLimits(0,0.8,1.2);
		fratio->SetParLimits(1,-0.1,0.1);
		fratio->SetParLimits(2,0.26,0.34);
		fratio->SetParLimits(3,0.26,0.34);
		fratio->FixParameter(4,(double) parn);
		fratio->FixParameter(5,xmin);
		fratio->FixParameter(6,xmax);
*/
		if(folder==2) fratio->SetParameters(3.4,2.5,0.2,1.2);
		if(folder==1) fratio->SetParameters(1.3,1,1,0.6);
		if(folder==0) fratio->SetParameters(15,9,0.1,1.0);
		if(folder==0) fratio->FixParameter(0,15);
		fratio->SetNpx(8000);
		tgratio[folder]->Fit(fratio,"ME");
/*
		double par0 = fratio->GetParameter(0);
		double par1 = fratio->GetParameter(1);
		double par2 = fratio->GetParameter(2);
		double par3 = fratio->GetParameter(3);
		ymin = fratio->Eval(xmin);
		ymax = fratio->Eval(xmax);
		double smin = -2.*ymin*( (xmin-par1)/pow(par2,2)-xmin/pow(par3,2) );
		double smax = -2.*ymax*( (xmax-par1)/pow(par2,2)-xmax/pow(par3,2) );
		double parbmin = -smin*pow(xmin,parn+1)/parn;
		double parbmax = -smax*pow(xmax,parn+1)/parn;
		double paramin = ymin + smin*xmin/parn;
		double paramax = ymax + smax*xmax/parn;
		cout << "Parameter min b: " << parbmin << "; ";
		cout << "Parameter min a: " << paramin << endl;
		cout << "Parameter max b: " << parbmax << "; ";
		cout << "Parameter max a: " << paramax << endl;
*/
//		foutput->WriteTObject(tgdata[folder]);
		foutput->WriteTObject(tgratio[folder]);
		foutput->WriteTObject(fratio);
		for (int hh = 0; hh < ntrees; hh++) delete hProj[hh][folder];
		delete fratio;
		delete tgratio[folder];
		delete tgdata[folder];
	}
	for (int folder = 0; folder < 3; folder++){
		for (int e = 0; e < 2; e++){
			if (finput[folder][e] == 0 || finput[folder][e]->IsZombie()) continue;
			else finput[folder][e]->Close();
		}
	}
	foutput->Close();
}


double formula_TxySyst(double* x, double* par){
	double xx = x[0];
	double par0 = par[0];
	double par1 = par[1];
	double par2 = par[2];
	double par3 = par[3];
	int parn = (int) par[4];
	double xmin = par[5];
	double xmax = par[6];
	double ymin = par0*exp( -pow((xmin-par1)/par2,2)+pow(xmin/par3,2) );
	double ymax = par0*exp( -pow((xmax-par1)/par2,2)+pow(xmax/par3,2) );
	double smin = -2.*ymin*( (xmin-par1)/pow(par2,2)-xmin/pow(par3,2) );
	double smax = -2.*ymax*( (xmax-par1)/pow(par2,2)-xmax/pow(par3,2) );
	double parbmin = -smin*pow(xmin,parn+1)/parn;
	double parbmax = -smax*pow(xmax,parn+1)/parn;
	double paramin = ymin + smin*xmin/parn;
	double paramax = ymax + smax*xmax/parn;

//	if(paramin>15 || paramax>15) return 0;

	if (xx <= xmax && xx >= xmin)
		return par0*exp( -pow((xx-par1)/par2,2)+pow(xx/par3,2) );
	else if(xx<xmin)
		return paramin + parbmin/pow(xx,parn);
	else
		return paramax + parbmax/pow(xx,parn);
}


void makeCombineSignalTemplateswithTrees_Bkg_single(int folder, int erg_tev, int Systematics){
	char TREE_NAME[]="SelectedTree";
	TString OUTPUT_NAME = "_templates_background_";
	OUTPUT_NAME.Prepend(user_folder[folder]);
	TString comstring;
	comstring.Form("%iTeV",erg_tev);
	TString erg_dir;
	erg_dir.Form("LHC_%iTeV",erg_tev);
	int EnergyIndex=0;if(erg_tev==8) EnergyIndex=1;

	if (Systematics == 0) OUTPUT_NAME += "Nominal.root";
	if (Systematics == 1) OUTPUT_NAME += "TxyUp.root";
	if (Systematics == -1) OUTPUT_NAME += "TxyDown.root";

	TString INPUT_SYST_NAME = "LifetimeKD_ShapeVariation_DataQQBZZ_Untransformed_Comparison.root";
	TString cinput_syst = user_dir_hep + "Analysis/Plots/" + INPUT_SYST_NAME;
	TFile* finput_syst = new TFile(cinput_syst, "read");
	TF1* txy_syst = (TF1*) finput_syst->Get(Form("tgratio_%s_fit",user_folder[folder]));

	float ZZMass_PeakCut[2]={105.6,140.6};
	double yield_zx[2][3] = {
		{ 0.2230,	0.6226,	1.0628 },
		{ 1.1878,	2.7676,	4.2929 }
	};
	double yield_ggzz[2][3]={
		{ 0.0625,	0.0341,	0.0741 },
		{ 0.4131,	0.2041,	0.5005 }
	};
	double yield_qqzz[2][3]={
		{ 1.7971,	0.8386,	2.2456 },
		{ 7.6478,	2.9364,	8.8585 }
	};
	double yield_signal[2][3]={
 		{1.0902689,	0.6087736,	1.4452079},
 		{5.1963998,	2.6944639,	6.6562629}
	};
	double luminosity[2] = { 5.051, 19.712 };

	float templateWeight=1;
	float MC_weight;
	float MC_weight_noxsec;
	float MC_weight_Kfactor=1;
	float MC_weight_QQBGGProper[4];
	float ZZMass,ZZPt,ZZEta,ZZPhi;
	float D_bkg;
	float KalmanCandVtx_x,KalmanCandVtx_y,KalmanCandVtx_z;
	float OfflinePrimaryVtx_x,OfflinePrimaryVtx_y,OfflinePrimaryVtx_z;

	float Txy,Tz,KD;
	float CandVtx_x,CandVtx_y,CandVtx_z;

	TString cinput_common = user_dir_hep + erg_dir + "/";
	TString coutput_common = user_dir_hep + "Analysis/Templates/" + comstring + "/";
	TString coutput = coutput_common + OUTPUT_NAME;
	TFile* foutput = new TFile(coutput,"recreate");

	const int ntrees=3;
	TString cinput_qqzz_common = cinput_common + user_folder[folder] + "/";
	TString cinput_ggzz_common = cinput_qqzz_common;
	TString cinput_zx_common = cinput_common;

	double expected_yield[ntrees] = {
		yield_qqzz[EnergyIndex][folder]/luminosity[EnergyIndex],
		yield_ggzz[EnergyIndex][folder]/luminosity[EnergyIndex],
		yield_zx[EnergyIndex][folder]/luminosity[EnergyIndex]
	};
	double scale_yield[ntrees] = {1,1,1};

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
	for (int ee = 7; ee < 9; ee++){
		for (int f = 0; f < 3; f++){
			TString cinput_zx = user_dir_hep + "LHC_";
			cinput_zx.Append(Form("%iTeV/", ee));
			cinput_zx = cinput_zx + user_folder[folder] + "/" + sample_FullSim[kAllSamples - 1] + "_Reprocessed.root";
			tzx->Add(cinput_zx);
		}
	}
	
	int nbinsx = 80;
	double xlow = -1, xhigh = 1;
	int nbinsy = 50;
	double ylow=0,yhigh=1;
	TH2F* hqqzz_full = new TH2F("h_template_qqZZ",Form("%i TeV %s q#bar{q}#rightarrow4l Background",erg_tev,user_folder[folder]),nbinsx,xlow,xhigh,nbinsy,ylow,yhigh);
	hqqzz_full->GetXaxis()->SetTitle("tanh(T_{xy} / 2 mm)");
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
			tc[tt]->SetBranchAddress("MC_weight_QQBGGProper", &MC_weight_QQBGGProper);
		}
		else tc[tt]->SetBranchAddress("ZXfake_weightProper", &MC_weight_noxsec);
		if (tt == 1) tc[tt]->SetBranchAddress("MC_weight_Kfactor", &MC_weight_Kfactor);

		tc[tt]->SetBranchAddress("D_bkg", &D_bkg);
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
		
		mytree[tt] = new TTree(Form("T_2D_%s", ctreefull[tt]), "");
		mytree[tt]->Branch("MC_weight",&templateWeight);
		mytree[tt]->Branch("D_bkg",&D_bkg);
		mytree[tt]->Branch("KD",&KD);
	}
	for (int tt = 0; tt < ntrees; tt++){
		for (int ev = 0; ev < tc[tt]->GetEntries(); ev++){
			MC_weight_noxsec = 1;
			MC_weight_Kfactor = 1;
			MC_weight_QQBGGProper[0]=1;MC_weight_QQBGGProper[1]=1;

			tc[tt]->GetEntry(ev);

			CandVtx_x = KalmanCandVtx_x - OfflinePrimaryVtx_x;
			CandVtx_y = KalmanCandVtx_y - OfflinePrimaryVtx_y;
			CandVtx_z = KalmanCandVtx_z - OfflinePrimaryVtx_z;

			Txy = (CandVtx_x*cos(ZZPhi) + CandVtx_y*sin(ZZPhi))*ZZMass / ZZPt;

			KD = Txy*5.; // per 2 mm
			double Txy_syst_wgt = txy_syst->Eval(KD);			
			KD = tanh(KD);
			if (KD >= 1.) KD = 1. - 1.0e-10;
			if (KD <= -1.) KD = -1. + 1.0e-10;

			if (ZZMass >= ZZMass_PeakCut[1] || ZZMass < ZZMass_PeakCut[0]) continue;
			double wgt = MC_weight_noxsec;
			if(tt==0) wgt *= MC_weight_QQBGGProper[1];
			if(tt==1) wgt *= MC_weight_Kfactor*MC_weight_QQBGGProper[0];

			if(Systematics==1) wgt *= Txy_syst_wgt;
			if(Systematics==-1) wgt *= 1./Txy_syst_wgt;

			templateWeight = wgt;
			hfull[tt]->Fill(KD,D_bkg,templateWeight);
		}
		scale_yield[tt] = expected_yield[tt]/hfull[tt]->Integral(0, hfull[tt]->GetNbinsX() + 1, 0, hfull[tt]->GetNbinsY() + 1);
		hfull[tt]->Scale(scale_yield[tt]);

		for (int ev = 0; ev < tc[tt]->GetEntries(); ev++){
			MC_weight_noxsec = 1;
			MC_weight_Kfactor = 1;
			MC_weight_QQBGGProper[0]=1;MC_weight_QQBGGProper[1]=1;

			tc[tt]->GetEntry(ev);

			CandVtx_x = KalmanCandVtx_x - OfflinePrimaryVtx_x;
			CandVtx_y = KalmanCandVtx_y - OfflinePrimaryVtx_y;
			CandVtx_z = KalmanCandVtx_z - OfflinePrimaryVtx_z;

			Txy = (CandVtx_x*cos(ZZPhi) + CandVtx_y*sin(ZZPhi))*ZZMass / ZZPt;

			KD = Txy*5.; // per 2 mm
			double Txy_syst_wgt = txy_syst->Eval(KD);			
			KD = tanh(KD);
			if (KD >= 1.) KD = 1. - 1.0e-10;
			if (KD <= -1.) KD = -1. + 1.0e-10;

			if (ZZMass >= ZZMass_PeakCut[1] || ZZMass < ZZMass_PeakCut[0]) continue;
			double wgt = MC_weight_noxsec;
			if(tt==0) wgt *= MC_weight_QQBGGProper[1];
			if(tt==1) wgt *= MC_weight_Kfactor*MC_weight_QQBGGProper[0];

			if(Systematics==1) wgt *= Txy_syst_wgt;
			if(Systematics==-1) wgt *= 1./Txy_syst_wgt;

			templateWeight = wgt*scale_yield[tt];
			mytree[tt]->Fill();
		}
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

	foutput->Close();
	finput_syst->Close();
}

