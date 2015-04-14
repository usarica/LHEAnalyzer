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
TString cutNames[nSIPCuts] = {
	"PVSIPcut",
	"NoSIPcut",
	"Z1SIPcut",
	"4lchi2cut",
	"Z1SIP_4lchi2cut",
	"CombinedCut"
};
TString strnewSIP[nSIPCuts] ={
	" max(Lep1SIP,max(Lep2SIP,max(Lep3SIP,Lep4SIP)))<4 ",
	" ",
	" abs(Lep1_Z1SIP)<4 && abs(Lep2_Z1SIP)<4 && abs(Lep3_Z1SIP)<5 && abs(Lep4_Z1SIP)<5 ",
	" KalmanCandVtx_chi2<30 ",
	" abs(Lep1_Z1SIP)<4 && abs(Lep2_Z1SIP)<4 && abs(Lep3_Z1SIP)<5 && abs(Lep4_Z1SIP)<5 && KalmanCandVtx_chi2<30 ",
	" abs(Lep1_Z1SIP)<4 && abs(Lep2_Z1SIP)<4 && abs(Lep3_Z1SIP)<5 && abs(Lep4_Z1SIP)<5 && KalmanCandVtx_chi2<30 && max(Lep1SIP,max(Lep2SIP,max(Lep3SIP,Lep4SIP)))<4 "
};
TString strnewSIP_CRZL[nSIPCuts] ={
	" max(Lep1SIP,max(Lep2SIP,Lep3SIP))<4 ",
	" ",
	" abs(Lep1_Z1SIP)<4 && abs(Lep2_Z1SIP)<4 && abs(Lep3_Z1SIP)<5 ",
	" KalmanCandVtx_chi2<30 ",
	" abs(Lep1_Z1SIP)<4 && abs(Lep2_Z1SIP)<4 && abs(Lep3_Z1SIP)<5 && KalmanCandVtx_chi2<30 ",
	" abs(Lep1_Z1SIP)<4 && abs(Lep2_Z1SIP)<4 && abs(Lep3_Z1SIP)<5 && KalmanCandVtx_chi2<30 && max(Lep1SIP,max(Lep2SIP,Lep3SIP))<4 "
};
TString strZ1Category[2]={
	"Z1ids == -169",
	"Z1ids == -121"
};
TString strZ2Category[6]={
	"(CRflag==6 || CRflag==10)", // DO require these, otherwise you will need mass cut, iso cut, ID as well, and has to be the best Z
	"Z2ids == -143",
	"(CRflag==8 || CRflag==12)",
	"(CRflag==7 || CRflag==11)",
	"Z2ids == 143",
	"(CRflag==5 || CRflag==9)"
};
TString strLep3Category[2]={
  "abs(Lep3ID) == 13",
  "abs(Lep3ID) == 11"
};
const int nCR=2;
TString strLepIsoCut[nCR]={
  " max(Lep1combRelIsoPF,max(Lep2combRelIsoPF,max(Lep3combRelIsoPF,Lep4combRelIsoPF)))<0.4 && Lep3isID && Lep4isID ",
//  " max(Lep1combRelIsoPF,max(Lep2combRelIsoPF,Lep3combRelIsoPF))<0.4 "
  " Lep3combRelIsoPF<0.4 && Lep3isID " // Lep1,2 already have the cut; what is seen on ntuples is the relIso before FSR
};
TString strLepIsoCut_pf[2]={
  " Lep3combRelIsoPF>=0.4 && Lep4combRelIsoPF>=0.4 && !Lep3isID && !Lep4isID ", // 2p2f
  " ( !(Lep3combRelIsoPF<0.4 && Lep4combRelIsoPF<0.4 && Lep3isID && Lep4isID) && !(Lep3combRelIsoPF>=0.4 && Lep4combRelIsoPF>=0.4 && !Lep3isID && !Lep4isID) ) " // 3p1f
};
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


void countSIPYields_CR(int erg_tev, int iCR=0, double m4l_low = 105.6, double m4l_high = 140.6);

TString produceCRname(int iCR, int iwgt, int isocut, int iM1cut, int catZ1, int catZ2, int icut, char* cappend);

void plot1DNewOld(int iCR, int iwgt, int isocut, int iM1cut, int catZ1, int catZ2, char* cappend, int erg_tev, TFile* frecord);
void plot2DNewOld(int iCR, int iwgt, int isocut, int iM1cut, int catZ1, int catZ2, char* cappend, int erg_tev, TFile* frecord);
void plot2DNewOldDoubleRatio(int iCR, int iwgt, int isocut, int iM1cut, int catZ1, int catZ2, char* cappend, int erg_tev, TFile* frecord);

void countSIPYields(int erg_tev=8, double m4l_low = 105.6, double m4l_high = 140.6){
	for (int iCR=0; iCR<2; iCR++) countSIPYields_CR(erg_tev, iCR, m4l_low, m4l_high);
}

void countSIPYields_CR(int erg_tev, int iCR, double m4l_low, double m4l_high){
	if (iCR>1) return;
	char TREE_NAME[]="SelectedTree";
	TString strWgt[2] ={ "", "ZXfake_weight" };

	TString erg_dir;
	erg_dir.Form("LHC_%iTeV", erg_tev);
	TString comstring;
	comstring.Form("%iTeV", erg_tev);
//	int EnergyIndex = 0; if (erg_tev == 8) EnergyIndex = 1;
	int sampleindex = kAllSamples-1 + iCR;

	TChain* tc = new TChain(TREE_NAME);

  TString cinput_common_noSIP = user_dir_hep + "No_SIP/" + erg_dir + "/CR/";
  TString cinput = cinput_common_noSIP;
	cinput = cinput + sample_FullSim[sampleindex] + ".root";
	tc->Add(cinput);
	if (tc->GetEntries() == 0){
		cout << "Could not find any files, aborting..." << endl; return;
	}

  TString coutput_common = user_dir_hep + "Analysis/Auxiliary/";

  TFile* frecord;
  TString fname = "CRdistributions_";
  if (iCR==0) fname.Append("ZLL_");
  if (iCR==1) fname.Append("ZL_");
  fname += comstring;
  fname.Append(Form("_m4l_%.1f_%.1f%s", m4l_low, m4l_high, ".root"));
  fname.Prepend(coutput_common);
  frecord = new TFile(fname, "recreate");

  TFile* foutput[3];
  TH1F* hYield[3]={ 0 };
  ofstream tout("logCR.log",ios::out);

	if (iCR==0 && m4l_low==105.6 && m4l_high==140.6){
		for (int ff=0; ff<3; ff++){
			TString OUTPUT_NAME = "LifetimeKD_RelativeSIPYields_";
			OUTPUT_NAME.Append(Form("%s", "CR"));
			OUTPUT_NAME.Append(Form("_%s_", user_folder[ff]));
			OUTPUT_NAME.Append(comstring);
			OUTPUT_NAME.Append(".root");
			TString coutput = coutput_common + OUTPUT_NAME;
			foutput[ff] = new TFile(coutput, "recreate");
			hYield[ff] = new TH1F("hCR", "Yield ratio to No SIP", nSIPCuts - 1, 0, nSIPCuts - 1);
		}
	}

  frecord->cd();

	TString hname = "htemp";
	TString strdraw = "Z1Mass>>"; strdraw.Append(hname);

  TString strdraw_Txy = "((KalmanCandVtx_x - OfflinePrimaryVtx_x)*cos(ZZPhi) + (KalmanCandVtx_x - OfflinePrimaryVtx_x)*sin(ZZPhi))*10000.0*ZZMass / ZZPt>>";
  TString strdraw_mZZ = "ZZMass>>";
  TString strdraw_mZ1 = "Z1Mass>>";
  TString strdraw_mZ2 = "Z2Mass>>";
  TString strdraw_pteta_ld[2] ={
    "(abs(Lep3Eta)*(Lep3Pt>Lep4Pt)+abs(Lep4Eta)*(Lep3Pt<=Lep4Pt)):(Lep3Pt*(Lep3Pt>Lep4Pt)+Lep4Pt*(Lep3Pt<=Lep4Pt))>>",
    "abs(Lep3Eta):Lep3Pt>>"
  };
  TString strdraw_pteta_subld = "(abs(Lep3Eta)*(Lep3Pt<=Lep4Pt)+abs(Lep4Eta)*(Lep3Pt>Lep4Pt)):(Lep3Pt*(Lep3Pt<=Lep4Pt)+Lep4Pt*(Lep3Pt>Lep4Pt))>>";

  for (int iwgt=0; iwgt<2; iwgt++){
    if (iCR==1 && iwgt==1) continue;
    for (int isocut=0; isocut<4; isocut++){
      if (iCR==1 && isocut>=2) continue;
      for (int iM1cut=0; iM1cut<3; iM1cut++){
        if (iM1cut==2 && iCR==0) continue;
        if (iwgt==0) tout << "Un-weighted";
        else tout << "Weighted";
        if (isocut>=1) tout << ", " << strLepIsoCut_label[isocut-1];
        tout << ", " << Z1masswidth_label[iM1cut] << endl;
        tout << "\t"; for (int icut=0; icut<nSIPCuts; icut++) tout << cutLabel[icut] << '\t';
        tout << endl;
        for (int catZ1=0; catZ1<2; catZ1++){
          for (int catZ2=0; catZ2<(iCR==0 ? 6 : 2); catZ2++){
            tout << strZ1Category_label[catZ1] << " ";
            if (iCR!=1) tout << strZ2Category_label[catZ2];
            else tout << strLep3Category_label[catZ2];
            for (int icut=0; icut<nSIPCuts; icut++){
              gStyle->SetTitleFont(62, "t");
              gROOT->SetStyle(gStyle->GetName());
              gROOT->ForceStyle();

              TH1F* hsip = new TH1F(hname, "", 1, 0, 1000);
              hsip->Sumw2();
              TString strcut = strnewSIP[icut];
              if (iCR==1) strcut = strnewSIP_CRZL[icut];
              if (iCR==0 && icut!=1) strcut.Append(Form("&& ZZMass>=%.1f && ZZMass<%.1f ", m4l_low, m4l_high));
              else if (iCR==0) strcut.Append(Form("ZZMass>=%.1f && ZZMass<%.1f ", m4l_low, m4l_high));
              if ((iCR==1 && icut!=1) || iCR==0) strcut.Append(" && ");
              strcut.Append(strZ1Category[catZ1]);
              strcut.Append(" && ");
              if (iCR==0) strcut.Append(strZ2Category[catZ2]);
              else strcut.Append(strLep3Category[catZ2]);
              if (isocut==0 && iCR==1){
                strcut.Append(" && ");
                strcut.Append("max(Lep1combRelIsoPF,Lep2combRelIsoPF)<0.4");
              }
              else if (isocut==1){
                strcut.Append(" && ");
                strcut.Append(strLepIsoCut[iCR]);
              }
              else if (isocut>1){
                strcut.Append(" && ");
                strcut.Append(strLepIsoCut_pf[isocut-2]);
              }
              if (iM1cut==1){
                strcut.Append(" && ");
                strcut.Append(Form("abs(Z1Mass-%.4f)<%.1f", Z1pole, Z1masswidth[iM1cut-1]));
              }
              else if (iM1cut==2){
                strcut.Append(" && ");
                strcut.Append(Form("abs(ZZMass-%.4f)<%.1f", Z1pole, Z1masswidth[iM1cut-1]));
              }
              if (iCR==1){
                strcut.Append(" && PFMET<25");
              }


              strcut.Append(")");
              strcut.Prepend("(");
              if (iwgt>0){
                strcut.Prepend("*");
                strcut.Prepend(strWgt[iwgt]);
              }
              tc->Draw(strdraw, strcut);

              double ptbins[9]={ 0, 5, 7, 10, 20, 30, 40, 50, 80 };
              double etabins[3]={ 0, 1.2, 2.5 };

              if ((iCR==1 && catZ2==1) || (iCR==0 && (catZ2==2 || catZ2==3))) etabins[1]=1.45;

              TH1F* hTxy = new TH1F(produceCRname(iCR, iwgt, isocut, iM1cut, catZ1, catZ2, icut, "Txy"), "", 60, -1400,1400);
              hTxy->Sumw2();
              TString newdraw = strdraw_Txy;
              newdraw.Append(hTxy->GetName());
              tc->Draw(newdraw, strcut);
              hTxy->SetXTitle("T_{xy} (#mum)");
              hTxy->SetYTitle("Events / 46.7 #mum");
              frecord->WriteTObject(hTxy);
              delete hTxy;
              TH1F* hmZZ = new TH1F(produceCRname(iCR, iwgt, isocut, iM1cut, catZ1, catZ2, icut, "mZZ"), "", (iCR!=1 ? 5 : 10), m4l_low, m4l_high);
              hmZZ->Sumw2();
              newdraw = strdraw_mZZ;
              newdraw.Append(hmZZ->GetName());
              tc->Draw(newdraw, strcut);
              if (iCR==0){
                hmZZ->SetXTitle("m_{4l} (GeV)");
                hmZZ->SetYTitle("Events / 7 GeV");
              }
              else if(iCR==1){
                hmZZ->SetXTitle("m_{3l} (GeV)");
                hmZZ->SetYTitle("Events / 3.5 GeV");
              }
              frecord->WriteTObject(hmZZ);
              delete hmZZ;
              if (iM1cut==0){
                TH1F* hmZ1 = new TH1F(produceCRname(iCR, iwgt, isocut, iM1cut, catZ1, catZ2, icut, "mZ1"), "", (iCR!=1 ? 10 : 20), 40, 120);
                hmZ1->Sumw2();
                newdraw = strdraw_mZ1;
                newdraw.Append(hmZ1->GetName());
                tc->Draw(newdraw, strcut);
                hmZ1->SetXTitle("m_{1} (GeV)");
                if (iCR==0) hmZ1->SetYTitle("Events / 8 GeV");
                else if (iCR==1) hmZ1->SetYTitle("Events / 4 GeV");
                frecord->WriteTObject(hmZ1);
                delete hmZ1;
              }
              if (iCR==0){
                TH1F* hmZ2 = new TH1F(produceCRname(iCR, iwgt, isocut, iM1cut, catZ1, catZ2, icut, "mZ2"), "", 17, 0, 136);
                hmZ2->Sumw2();
                newdraw = strdraw_mZ2;
                newdraw.Append(hmZ2->GetName());
                tc->Draw(newdraw, strcut);
                hmZ2->SetXTitle("m_{2} (GeV)");
                hmZ2->SetYTitle("Events / 8 GeV");
                frecord->WriteTObject(hmZ2);
                delete hmZ2;
              }
              TH2F* hpteta = new TH2F(produceCRname(iCR, iwgt, isocut, iM1cut, catZ1, catZ2, icut, "pTeta_leading"), "", 8, ptbins, 2, etabins);
              hpteta->Sumw2();
              newdraw = strdraw_pteta_ld[iCR];
              newdraw.Append(hpteta->GetName());
              tc->Draw(newdraw, strcut);
              if (iCR==0) hpteta->SetXTitle("Leading p^{3,4}_{T} (GeV)");
              else if (iCR==1) hpteta->SetXTitle("p^{3}_{T} (GeV)");
              hpteta->SetYTitle("|#eta|");
              hpteta->SetZTitle("Events / bin");
              frecord->WriteTObject(hpteta);
              delete hpteta;
              if (iCR==0){
                hpteta = new TH2F(produceCRname(iCR, iwgt, isocut, iM1cut, catZ1, catZ2, icut, "pTeta_subleading"), "", 8, ptbins, 2, etabins);
                hpteta->Sumw2();
                newdraw = strdraw_pteta_subld;
                newdraw.Append(hpteta->GetName());
                tc->Draw(newdraw, strcut);
                hpteta->SetXTitle("Sub-leading p^{3,4}_{T} (GeV)");
                hpteta->SetYTitle("|#eta|");
                hpteta->SetZTitle("Events / bin");
                frecord->WriteTObject(hpteta);
                delete hpteta;
              }

              tout << '\t' << hsip->Integral();

							if (iCR==0 && m4l_low==105.6 && m4l_high==140.6){
								if (icut<nSIPCuts-1 && !(catZ2==1 || catZ2==4) && isocut==0 && iwgt==1 && iM1cut==0){
									int channeltype=2; // 2e2mu + 2mu2e
									if (catZ1==0 && (catZ2==0 || catZ2==5)) channeltype=0; // 4mu
									else if (catZ1==1 && (catZ2==2 || catZ2==3)) channeltype=1; // 4e
									else if ( (catZ1==1 && (catZ2==0 || catZ2==5))
											|| (catZ1==0 && (catZ2==2 || catZ2==3)) ) channeltype=2; // 2e2mu
									else{
										cout << "Invalid cut" << endl; channeltype=-1;
									}
									if(channeltype>=0) hYield[channeltype]->AddBinContent(icut+1, hsip->Integral());
								}
							}

              delete hsip;
            }
            tout << endl;
            plot1DNewOld(iCR, iwgt, isocut, iM1cut, catZ1, catZ2, "Txy", erg_tev, frecord);
            plot1DNewOld(iCR, iwgt, isocut, iM1cut, catZ1, catZ2, "mZZ", erg_tev, frecord);
            if (iM1cut==0) plot1DNewOld(iCR, iwgt, isocut, iM1cut, catZ1, catZ2, "mZ1", erg_tev, frecord);
            if (iCR==0) plot1DNewOld(iCR, iwgt, isocut, iM1cut, catZ1, catZ2, "mZ2", erg_tev, frecord);
            plot2DNewOld(iCR, iwgt, isocut, iM1cut, catZ1, catZ2, "pTeta_leading", erg_tev, frecord);
            if(iCR==0) plot2DNewOld(iCR, iwgt, isocut, iM1cut, catZ1, catZ2, "pTeta_subleading", erg_tev, frecord);
          }
        }
        tout << endl;
      }
      tout << endl;
    }
    tout << endl;

    for (int isocut=1; isocut<4; isocut++){
      if (iCR==1 && isocut>=2) continue;
      for (int catZ1=0; catZ1<2; catZ1++){
        for (int catZ2=0; catZ2<(iCR==0 ? 6 : 2); catZ2++){
          for (int iM1cut=0; iM1cut<3; iM1cut++){
            if (iM1cut==2 && iCR==0) continue;
            plot2DNewOldDoubleRatio(iCR, iwgt, isocut, iM1cut, catZ1, catZ2, "pTeta_leading", erg_tev, frecord);
            if (iCR==0) plot2DNewOldDoubleRatio(iCR, iwgt, isocut, iM1cut, catZ1, catZ2, "pTeta_subleading", erg_tev, frecord);
          }
        }
      }
    }

  }
  delete tc;

	if (iCR==0 && m4l_low==105.6 && m4l_high==140.6){
		for (int ff=0; ff<3; ff++){
			double nocutabs = hYield[ff]->GetBinContent(2);
			for (int bin=1; bin<=nSIPCuts-1; bin++) hYield[ff]->SetBinContent(bin, hYield[ff]->GetBinContent(bin) / nocutabs);
			tout << "Yield scaling ratio for channel " << ff << ": " << hYield[ff]->GetBinContent(nSIPCuts-1)/hYield[ff]->GetBinContent(1) << endl;
      foutput[ff]->cd();
			foutput[ff]->WriteTObject(hYield[ff]);
			delete hYield[ff];
			foutput[ff]->Close();
		}
	}
  tout.close();
  frecord->Close();
}


TString produceCRname(int iCR, int iwgt, int isocut, int iM1cut, int catZ1, int catZ2, int icut, char* cappend){
  TString name;
  if (iCR==0) name = "ZLL";
  else name = "ZL";

  if (catZ1==0) name.Append("_Zmumu");
  else name.Append("_Zee");
  if (iCR==0){
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

void plot1DNewOld(int iCR, int iwgt, int isocut, int iM1cut, int catZ1, int catZ2, char* cappend, int erg_tev, TFile* frecord){
  TString comstring;
  comstring.Form("%iTeV", erg_tev);
  TString coutput_common = user_dir_hep + "Analysis/Plots/CR_FakeRates/Distributions/";
  coutput_common.Append(comstring);
  coutput_common.Append("/");
  if (iCR==0) coutput_common = coutput_common + "ZLL/";
  else coutput_common = coutput_common +  "ZL/";
  if (isocut==0) coutput_common.Append("AllLoose");
  else coutput_common.Append(strLepIsoCut_label[isocut-1]);
  coutput_common.Append("/");
  if (iM1cut==0) coutput_common.Append("m1_All/");
  else if (iM1cut==1) coutput_common.Append("m1_ZPeak/");
  else if (iM1cut==2) coutput_common.Append("m1_ZGPeak/");

  TString canvasdir = coutput_common;
  canvasdir.Append(Form("%s/", cappend));

  TString mkdirCommand = "mkdir -p ";
  mkdirCommand.Append(canvasdir);
  gSystem->Exec(mkdirCommand);

  TString hname[2] = {
    produceCRname(iCR, iwgt, isocut, iM1cut, catZ1, catZ2, 4, cappend),
    produceCRname(iCR, iwgt, isocut, iM1cut, catZ1, catZ2, 0, cappend)
  };
  TString hlabel[2] ={
    "New cut", "SIP_{3D}<4"
  };
  TH1F* htemp[2] ={
    (TH1F*)frecord->Get(hname[0]),
    (TH1F*)frecord->Get(hname[1])
  };

  if (iCR==1){
    htemp[0]->Scale(0.001);
    htemp[1]->Scale(0.001);
    htemp[0]->SetYTitle(Form("%s / 1000", htemp[0]->GetYaxis()->GetTitle()));
    htemp[1]->SetYTitle(Form("%s / 1000", htemp[1]->GetYaxis()->GetTitle()));
  }

  TString identifier_label = "#frac{";
  identifier_label = identifier_label + strZ1Category_label[catZ1];
  if (iCR!=1) identifier_label = identifier_label + strZ2Category_label[catZ2];
  else identifier_label = identifier_label + strLep3Category_label[catZ2];
  identifier_label = identifier_label + "}{";
  if (isocut>1) identifier_label = identifier_label + strLepIsoCut_label[isocut-1];
  else if (isocut==1 && iCR==0) identifier_label = identifier_label + "4p";
  else if (isocut==1 && iCR==1) identifier_label = identifier_label + "3p";
  if(isocut!=0) identifier_label = identifier_label + ", ";
  identifier_label = identifier_label + Z1masswidth_label[iM1cut] + "}";

  TString canvasname = produceCRname(iCR, iwgt, isocut, iM1cut, catZ1, catZ2, 5, cappend);
  canvasname = canvasname + "_" + comstring;

  frecord->cd();

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

  TLegend *l2D = new TLegend(0.76, 0.74, 0.98, 0.90);
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

  htemp[0]->GetYaxis()->SetRangeUser(0, (TMath::Max(htemp[0]->GetMaximum() + htemp[0]->GetBinError(htemp[0]->GetMaximumBin()), htemp[1]->GetMaximum() + htemp[1]->GetBinError(htemp[1]->GetMaximumBin()))) * 1.5);
  htemp[0]->GetXaxis()->SetNdivisions(505);
  htemp[0]->GetXaxis()->SetLabelFont(42);
  htemp[0]->GetXaxis()->SetLabelOffset(0.007);
  htemp[0]->GetXaxis()->SetLabelSize(0.04);
  htemp[0]->GetXaxis()->SetTitleSize(0.06);
  htemp[0]->GetXaxis()->SetTitleOffset(0.9);
  htemp[0]->GetXaxis()->SetTitleFont(42);
  htemp[0]->GetYaxis()->SetNdivisions(505);
  htemp[0]->GetYaxis()->SetLabelFont(42);
  htemp[0]->GetYaxis()->SetLabelOffset(0.007);
  htemp[0]->GetYaxis()->SetLabelSize(0.04);
  htemp[0]->GetYaxis()->SetTitleSize(0.06);
  htemp[0]->GetYaxis()->SetTitleOffset(1.1);
  htemp[0]->GetYaxis()->SetTitleFont(42);

  htemp[0]->SetLineWidth(2);
  htemp[1]->SetLineWidth(2);
  htemp[0]->SetLineColor(kBlue);
  htemp[1]->SetLineColor(kRed);

  l2D->AddEntry(htemp[0], hlabel[0], "lf");
  l2D->AddEntry(htemp[1], hlabel[1], "lf");

  htemp[0]->Draw("hist");
  htemp[1]->Draw("histsame");

  TH1F* hclone[2]={
    (TH1F*)htemp[0]->Clone("h1"),
    (TH1F*)htemp[1]->Clone("h2")
  };
    
  hclone[0]->SetFillColor(kBlue);
  hclone[1]->SetFillColor(kRed);
  hclone[0]->SetFillStyle(3004);
  hclone[1]->SetFillStyle(3005);
  hclone[0]->Draw("e2same");
  hclone[1]->Draw("e2same");
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

  frecord->WriteTObject(c2D);

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

  delete hclone[0];
  delete hclone[1];
  c2D->Close();
}

void plot2DNewOld(int iCR, int iwgt, int isocut, int iM1cut, int catZ1, int catZ2, char* cappend, int erg_tev, TFile* frecord){
  TString comstring;
  comstring.Form("%iTeV", erg_tev);
  TString coutput_common = user_dir_hep + "Analysis/Plots/CR_FakeRates/Distributions/";
  coutput_common.Append(comstring);
  coutput_common.Append("/");
  if (iCR==0) coutput_common = coutput_common + "ZLL/";
  else coutput_common = coutput_common +  "ZL/";
  if (isocut==0) coutput_common.Append("AllLoose");
  else coutput_common.Append(strLepIsoCut_label[isocut-1]);
  coutput_common.Append("/");
  if (iM1cut==0) coutput_common.Append("m1_All/");
  else if (iM1cut==1) coutput_common.Append("m1_ZPeak/");
  else if (iM1cut==2) coutput_common.Append("m1_ZGPeak/");

  TString canvasdir = coutput_common;
  canvasdir.Append(Form("%s/", cappend));

  TString mkdirCommand = "mkdir -p ";
  mkdirCommand.Append(canvasdir);
  gSystem->Exec(mkdirCommand);

  TString hname[2] ={
    produceCRname(iCR, iwgt, isocut, iM1cut, catZ1, catZ2, 4, cappend),
    produceCRname(iCR, iwgt, isocut, iM1cut, catZ1, catZ2, 0, cappend)
  };
  TString hlabel[4] ={
    "New cut (All)", "SIP_{3D}<4 (All)",
    "New cut (Endcap)", "SIP_{3D}<4 (Endcap)"
  };
  TH2F* htemp[2] ={
    (TH2F*)frecord->Get(hname[0]),
    (TH2F*)frecord->Get(hname[1])
  };

  if (iCR==1){
    htemp[0]->Scale(0.001);
    htemp[1]->Scale(0.001);
    htemp[0]->GetZaxis()->SetTitle(Form("%s / 1000", htemp[0]->GetZaxis()->GetTitle()));
    htemp[1]->GetZaxis()->SetTitle(Form("%s / 1000", htemp[1]->GetZaxis()->GetTitle()));
  }

  TH1F* htemp_1D[4] ={
    (TH1F*)htemp[0]->ProjectionX("hnew0_12", 0, -1, "e"),
    (TH1F*)htemp[1]->ProjectionX("hnew1_12", 0, -1, "e"),
    (TH1F*)htemp[0]->ProjectionX("hnew0_2", 0, -1, "e"),
    (TH1F*)htemp[1]->ProjectionX("hnew1_2", 0, -1, "e")
  };
  TH1F* hclone[4];
  for (int hh=0; hh<2; hh++){
    for (int binx=0; binx<=htemp_1D[hh]->GetNbinsX()+1; binx++){
      double inclusive = htemp_1D[hh]->GetBinContent(binx);
      double inclusive_err = htemp_1D[hh]->GetBinError(binx);
      double exclusive = htemp[hh]->GetBinContent(binx,1);
      double exclusive_err = htemp[hh]->GetBinError(binx,1);
      double bincontent = inclusive - exclusive;
      double binerror = sqrt(pow(inclusive_err, 2) - pow(exclusive_err, 2));
      htemp_1D[2+hh]->SetBinContent(binx, bincontent);
      htemp_1D[2+hh]->SetBinError(binx, binerror);
    }
  }

  TString identifier_label = "#frac{";
  identifier_label = identifier_label + strZ1Category_label[catZ1];
  if (iCR!=1) identifier_label = identifier_label + strZ2Category_label[catZ2];
  else identifier_label = identifier_label + strLep3Category_label[catZ2];
  identifier_label = identifier_label + "}{";
  if (isocut>1) identifier_label = identifier_label + strLepIsoCut_label[isocut-1];
  else if (isocut==1 && iCR==0) identifier_label = identifier_label + "4p";
  else if (isocut==1 && iCR==1) identifier_label = identifier_label + "3p";
  if (isocut!=0) identifier_label = identifier_label + ", ";
  identifier_label = identifier_label + Z1masswidth_label[iM1cut] + "}";

  TString canvasname = produceCRname(iCR, iwgt, isocut, iM1cut, catZ1, catZ2, 5, cappend);
  canvasname = canvasname + "_" + comstring;

  frecord->cd();

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

  TLegend *l2D = new TLegend(0.61, 0.74, 0.98, 0.90);
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

  for (int hh=0; hh<2; hh++){
    htemp_1D[2*hh]->GetYaxis()->SetRangeUser(0, (TMath::Max(htemp_1D[2*hh]->GetMaximum() + htemp_1D[2*hh]->GetBinError(htemp_1D[2*hh]->GetMaximumBin()), htemp_1D[2*hh+1]->GetMaximum() + htemp_1D[2*hh+1]->GetBinError(htemp_1D[2*hh+1]->GetMaximumBin()))) * 1.5);
    htemp_1D[2*hh]->GetXaxis()->SetNdivisions(505);
    htemp_1D[2*hh]->GetXaxis()->SetLabelFont(42);
    htemp_1D[2*hh]->GetXaxis()->SetLabelOffset(0.007);
    htemp_1D[2*hh]->GetXaxis()->SetLabelSize(0.04);
    htemp_1D[2*hh]->GetXaxis()->SetTitleSize(0.06);
    htemp_1D[2*hh]->GetXaxis()->SetTitleOffset(0.9);
    htemp_1D[2*hh]->GetXaxis()->SetTitleFont(42);
    htemp_1D[2*hh]->GetYaxis()->SetNdivisions(505);
    htemp_1D[2*hh]->GetYaxis()->SetLabelFont(42);
    htemp_1D[2*hh]->GetYaxis()->SetLabelOffset(0.007);
    htemp_1D[2*hh]->GetYaxis()->SetLabelSize(0.04);
    htemp_1D[2*hh]->GetYaxis()->SetTitleSize(0.06);
    htemp_1D[2*hh]->GetYaxis()->SetTitleOffset(1.1);
    htemp_1D[2*hh]->GetYaxis()->SetTitleFont(42);

    htemp_1D[2*hh]->SetYTitle(htemp[0]->GetZaxis()->GetTitle());
    htemp_1D[2*hh+1]->SetYTitle(htemp[0]->GetZaxis()->GetTitle());

    htemp_1D[2*hh]->SetLineWidth(2);
    htemp_1D[2*hh+1]->SetLineWidth(2);
    htemp_1D[2*hh]->SetLineColor(kBlue);
    htemp_1D[2*hh+1]->SetLineColor(kRed);

    if (hh==1){
      htemp_1D[2*hh]->SetLineStyle(7);
      htemp_1D[2*hh+1]->SetLineStyle(7);
    }

    l2D->AddEntry(htemp_1D[2*hh], hlabel[2*hh], "lf");
    l2D->AddEntry(htemp_1D[2*hh+1], hlabel[2*hh+1], "lf");

    if (hh==0) htemp_1D[2*hh]->Draw("hist");
    else htemp_1D[2*hh]->Draw("histsame");
    htemp_1D[2*hh+1]->Draw("histsame");

    hclone[2*hh]=(TH1F*) htemp_1D[2*hh]->Clone(Form("h%i", 2*hh));
    hclone[2*hh+1]=(TH1F*) htemp_1D[2*hh+1]->Clone(Form("h%i", 2*hh+1));

    hclone[2*hh]->SetFillColor(kBlue);
    hclone[2*hh+1]->SetFillColor(kRed);
    hclone[2*hh]->SetFillStyle(3004);
    hclone[2*hh+1]->SetFillStyle(3005);
    hclone[2*hh]->Draw("e2same");
    hclone[2*hh+1]->Draw("e2same");
  }
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

  frecord->WriteTObject(c2D);

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

  for (int hh=0; hh<4; hh++){
    delete hclone[hh];
    delete htemp_1D[hh];
  }
  c2D->Close();
}

void plot2DNewOldDoubleRatio(int iCR, int iwgt, int isocut, int iM1cut, int catZ1, int catZ2, char* cappend, int erg_tev, TFile* frecord){
  TString comstring;
  comstring.Form("%iTeV", erg_tev);
  TString coutput_common = user_dir_hep + "Analysis/Plots/CR_FakeRates/DoubleRatios/";
  coutput_common.Append(comstring);
  coutput_common.Append("/");
  if (iCR==0) coutput_common = coutput_common + "ZLL/";
  else coutput_common = coutput_common +  "ZL/";
  if (isocut==0) return;
  else coutput_common.Append(strLepIsoCut_label[isocut-1]);
  coutput_common.Append("/");
  if (iM1cut==0) coutput_common.Append("m1_All/");
  else if (iM1cut==1) coutput_common.Append("m1_ZPeak/");
  else if (iM1cut==2) coutput_common.Append("m1_ZGPeak/");

  TString canvasdir = coutput_common;
  canvasdir.Append(Form("%s/", cappend));

  TString mkdirCommand = "mkdir -p ";
  mkdirCommand.Append(canvasdir);
  gSystem->Exec(mkdirCommand);

  TString hname[6] ={
    produceCRname(iCR, iwgt, isocut, iM1cut, catZ1, catZ2, 4, cappend),
    produceCRname(iCR, iwgt, isocut, iM1cut, catZ1, catZ2, 0, cappend),
    produceCRname(iCR, iwgt, isocut, iM1cut, catZ1, catZ2, 5, cappend),
    produceCRname(iCR, iwgt, 0, iM1cut, catZ1, catZ2, 4, cappend),
    produceCRname(iCR, iwgt, 0, iM1cut, catZ1, catZ2, 0, cappend),
    produceCRname(iCR, iwgt, 0, iM1cut, catZ1, catZ2, 5, cappend)
  };
  TString hlabel[2][2] ={
    { "New cut (Barrel)", "New cut (Endcap)" },
    { "SIP_{3D}<4 (Barrel)", "SIP_{3D}<4 (Endcap)" }
  };
  TH2F* htemp[6];
  for(int hh=0;hh<6;hh++){
    htemp[hh] = (TH2F*)frecord->Get(hname[hh]);
    if (iCR==1){
      htemp[hh]->Scale(0.001);
      htemp[hh]->GetZaxis()->SetTitle(Form("%s / 1000", htemp[hh]->GetZaxis()->GetTitle()));
    }
  }

  double maxplot=0;
  TH1F* htemp_1D[2][2];
  for(int hh=0;hh<2;hh++){
    for (int it=0; it<2; it++){
      TString hrationame = Form("%s_%s_ratio", htemp[hh]->GetName(), "Barrel");
      if(it==1) hrationame = Form("%s_%s_ratio", htemp[hh]->GetName(), "Endcap");
      htemp_1D[hh][it] = (TH1F*)htemp[hh]->ProjectionX(hrationame);
      htemp_1D[hh][it]->SetXTitle(htemp[hh]->GetXaxis()->GetTitle());
      htemp_1D[hh][it]->SetYTitle("Fake rate");
      if(iCR==0) htemp_1D[hh][it]->SetYTitle("Ratio");
      for (int bin=0; bin<=htemp[hh]->GetNbinsX()+1; bin++){
        double tightval[2] ={ htemp[hh]->GetBinContent(bin, it+1), htemp[hh]->GetBinError(bin, it+1) };
        double looseval[2] ={ htemp[hh+3]->GetBinContent(bin, it+1), htemp[hh+3]->GetBinError(bin, it+1) };

        double commontightval[2] ={ htemp[2]->GetBinContent(bin, it+1), htemp[2]->GetBinError(bin, it+1) };
        double commonlooseval[2] ={ htemp[5]->GetBinContent(bin, it+1), htemp[5]->GetBinError(bin, it+1) };

        double bincontent = 0;
        double binerror = 0;
        if (looseval[0]>0){
          bincontent = tightval[0] / looseval[0];
          double binerror_num = sqrt(pow(tightval[1]*(looseval[0] - tightval[0]), 2) + pow(tightval[0]*sqrt(pow(looseval[1], 2) - pow(tightval[1], 2)), 2));
          binerror = binerror_num / pow(looseval[0], 2);
        }
        htemp_1D[hh][it]->SetBinContent(bin, bincontent);
        htemp_1D[hh][it]->SetBinError(bin, binerror);

        if (maxplot<(bincontent + binerror))maxplot=(bincontent + binerror);
      }
    }
  }
//  cout << maxplot << endl;
  TH1F* hclone[2][2];

  TString identifier_label = "#frac{";
  identifier_label = identifier_label + strZ1Category_label[catZ1];
  if (iCR!=1) identifier_label = identifier_label + strZ2Category_label[catZ2];
  else identifier_label = identifier_label + strLep3Category_label[catZ2];
  identifier_label = identifier_label + "}{";
  if (isocut>1) identifier_label = identifier_label + strLepIsoCut_label[isocut-1];
  else if (isocut==1 && iCR==0) identifier_label = identifier_label + "4p";
  else if (isocut==1 && iCR==1) identifier_label = identifier_label + "3p";
  if (isocut!=0) identifier_label = identifier_label + ", ";
  identifier_label = identifier_label + Z1masswidth_label[iM1cut] + "}";

  TString canvasname = produceCRname(iCR, iwgt, isocut, iM1cut, catZ1, catZ2, 5, cappend);
  canvasname = canvasname + "_" + comstring + "_RatioComparison";

  frecord->cd();

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

  TLegend *l2D = new TLegend(0.61, 0.74, 0.98, 0.90);
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

  for (int hh=0; hh<2; hh++){
    htemp_1D[hh][0]->GetYaxis()->SetRangeUser(0,maxplot * 1.5);
    htemp_1D[hh][0]->GetXaxis()->SetNdivisions(505);
    htemp_1D[hh][0]->GetXaxis()->SetLabelFont(42);
    htemp_1D[hh][0]->GetXaxis()->SetLabelOffset(0.007);
    htemp_1D[hh][0]->GetXaxis()->SetLabelSize(0.04);
    htemp_1D[hh][0]->GetXaxis()->SetTitleSize(0.06);
    htemp_1D[hh][0]->GetXaxis()->SetTitleOffset(0.9);
    htemp_1D[hh][0]->GetXaxis()->SetTitleFont(42);
    htemp_1D[hh][0]->GetYaxis()->SetNdivisions(505);
    htemp_1D[hh][0]->GetYaxis()->SetLabelFont(42);
    htemp_1D[hh][0]->GetYaxis()->SetLabelOffset(0.007);
    htemp_1D[hh][0]->GetYaxis()->SetLabelSize(0.04);
    htemp_1D[hh][0]->GetYaxis()->SetTitleSize(0.06);
    htemp_1D[hh][0]->GetYaxis()->SetTitleOffset(1.1);
    htemp_1D[hh][0]->GetYaxis()->SetTitleFont(42);


    htemp_1D[hh][0]->SetLineWidth(2);
    htemp_1D[hh][1]->SetLineWidth(2);

    if (hh==0){
      htemp_1D[hh][0]->SetLineColor(kBlue);
      htemp_1D[hh][1]->SetLineColor(kBlue);
    }
    else if (hh==1){
      htemp_1D[hh][0]->SetLineColor(kRed);
      htemp_1D[hh][1]->SetLineColor(kRed);
    }
    htemp_1D[hh][1]->SetLineStyle(7);

    frecord->WriteTObject(htemp_1D[hh][0]);
    frecord->WriteTObject(htemp_1D[hh][1]);

    l2D->AddEntry(htemp_1D[hh][0], hlabel[hh][0], "lf");
    l2D->AddEntry(htemp_1D[hh][1], hlabel[hh][1], "lf");

    if (hh==0) htemp_1D[hh][0]->Draw("hist");
    else htemp_1D[hh][0]->Draw("histsame");
    htemp_1D[hh][1]->Draw("histsame");

    hclone[hh][0]=(TH1F*)htemp_1D[hh][0]->Clone(Form("h%i_0", hh));
    hclone[hh][1]=(TH1F*)htemp_1D[hh][1]->Clone(Form("h%i_1", hh+1));

    if (hh==0){
      hclone[hh][0]->SetFillColor(kBlue);
      hclone[hh][1]->SetFillColor(kBlue);
      hclone[hh][0]->SetFillStyle(3004);
      hclone[hh][1]->SetFillStyle(3006);
    }
    else if (hh==1){
      hclone[hh][0]->SetFillColor(kRed);
      hclone[hh][1]->SetFillColor(kRed);
      hclone[hh][0]->SetFillStyle(3005);
      hclone[hh][1]->SetFillStyle(3007);
    }
    hclone[hh][0]->Draw("e2same");
    hclone[hh][1]->Draw("e2same");
  }
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

  frecord->WriteTObject(c2D);

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

  for (int hh=0; hh<2; hh++){
    for (int it=0; it<2; it++){
      delete hclone[hh][it];
      delete htemp_1D[hh][it];
    }
  }
  c2D->Close();
}
