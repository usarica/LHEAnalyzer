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


double ZZMass_Range[2]={ 105.6, 140.6 };

TString cutLabel[6] ={
  "Old cut",
  "No cut",
  "Z1-SIP",
  "chi**2",
  "New cut",
  "New + old"
};

void listDataEvents_single(int folder, int erg_tev, int removeSIP=0);

void listDataEvents(){
  for (int ee=7; ee<9; ee++){
    for (int f=0; f<3; f++){
      for (int r=0; r<3; r++) listDataEvents_single(f, ee, r);
    }
  }
}

void listDataEvents_single(int folder, int erg_tev, int removeSIP){
	char TREE_NAME[]="SelectedTree";
	TString erg_dir;
	erg_dir.Form("LHC_%iTeV",erg_tev);
	TString comstring;
	comstring.Form("%iTeV",erg_tev);
	int EnergyIndex=0;if(erg_tev==8) EnergyIndex=1;

	float ZZMass,ZZPt,ZZEta,ZZPhi;
	float KalmanCandVtx_x, KalmanCandVtx_y, KalmanCandVtx_z, KalmanCandVtx_chi2;
  float OfflinePrimaryVtx_x, OfflinePrimaryVtx_y, OfflinePrimaryVtx_z, OfflinePrimaryVtx_ndof;
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
  Long64_t EventNumber;
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
  TString coutput_common = user_dir_hep + "/Analysis/Auxiliary/";
  TString OUTPUT_NAME = "logData_";
  OUTPUT_NAME.Append(comstring);
  OUTPUT_NAME.Append("_");
  OUTPUT_NAME.Append(user_folder[folder]);
  OUTPUT_NAME.Append("_");
  if (removeSIP==0) OUTPUT_NAME.Append("OldSIP");
  else if (removeSIP==1) OUTPUT_NAME.Append("NoSIP");
  else if (removeSIP==2) OUTPUT_NAME.Append("NewSIP");
  OUTPUT_NAME.Append(".log");
  OUTPUT_NAME.Prepend(coutput_common);
  ofstream tout(OUTPUT_NAME.Data(), ios::out);

  // BeamSpot info
  TChain* tBeam;
  TString strBeamSpot[2]={
    "data_GR_R_44_V15C",
    "data_FT_53_V21_AN4"
  };
  tBeam = new TChain("BeamSpotRecord");
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

  const int ntrees=1;
  TString cinput_data = cinput_common_noSIP + user_folder[3] + "/" + data_files[folder] + ".root";
	TChain* tdata = new TChain(TREE_NAME);
	tdata->Add(cinput_data);
	TChain* tc[ntrees] = { tdata };

  for (int tt = 0; tt < ntrees; tt++){
    cout << "Tree: " << tt << ", nEvents: " << tc[tt]->GetEntries() << endl;

    tc[tt]->SetBranchAddress("RunNumber", &RunNumber);
    tc[tt]->SetBranchAddress("LumiNumber", &LumiNumber);
    tc[tt]->SetBranchAddress("EventNumber", &EventNumber);
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

  tout
    << "RunNumber" << "\t"
    << "LumiNumber" << "\t"
    << "EventNumber" << "\t"

    << "ZZMass" << "\t"
    << "ZZPt" << "\t"

    << "Lep1SIP" << "\t"
    << "Lep2SIP" << "\t"
    << "Lep3SIP" << "\t"
    << "Lep4SIP" << "\t"

    << "Lep1_Z1SIP" << "\t"
    << "Lep2_Z1SIP" << "\t"
    << "Lep3_Z1SIP" << "\t"
    << "Lep4_Z1SIP" << "\t"
    << "KalmanCandVtx_chi2" << "\t"


    << "BeamPosX" << "\t"
    << "BeamPosY" << "\t"
    << "BeamPosZ" << "\t"
    << "KalmanCandVtx_x" << "\t"
    << "KalmanCandVtx_y" << "\t"
    << "KalmanCandVtx_y" << "\t"
    << "OfflinePrimaryVtx_x" << "\t"
    << "OfflinePrimaryVtx_y" << "\t"
    << "OfflinePrimaryVtx_z" << "\t"
    << "(OfflinePrimaryVtx_ndof+3.)/2." << "\t"

    << "Txy" << "\t"
    << "Dxy" << "\t"
    << "Txy_BS" << "\t"
    << "D_bkg" << "\t"

    << endl;


  for (int tt = 0; tt < ntrees; tt++){
    int ctrEvents = 0;
    for (int ev = 0; ev < tc[tt]->GetEntries(); ev++){
      RunNumber = 1;
      LumiNumber = 1;

      Lep1_Z1SIP=0; Lep2_Z1SIP=0; Lep3_Z1SIP=0; Lep4_Z1SIP=0; KalmanCandVtx_chi2=0;
      Lep1SIP=0; Lep2SIP=0; Lep3SIP=0; Lep4SIP=0;

      tc[tt]->GetEntry(ev);

      if (ZZMass>=ZZMass_Range[1] || ZZMass<ZZMass_Range[0]) continue;
      if (removeSIP==2 && (fabs(Lep1_Z1SIP) >= 4 || fabs(Lep2_Z1SIP) >= 4 || fabs(Lep3_Z1SIP) >= 5 || fabs(Lep4_Z1SIP) >= 5 || KalmanCandVtx_chi2 >= 30)) continue;
      if (removeSIP==0 && (fabs(Lep1SIP) >= 4 || fabs(Lep2SIP) >= 4 || fabs(Lep3SIP) >= 4 || fabs(Lep4SIP) >= 4)) continue;

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

      CandVtx_x = KalmanCandVtx_x - OfflinePrimaryVtx_x;
      CandVtx_y = KalmanCandVtx_y - OfflinePrimaryVtx_y;
      CandVtx_z = KalmanCandVtx_z - OfflinePrimaryVtx_z;

      Dxy = (CandVtx_x*cos(ZZPhi) + CandVtx_y*sin(ZZPhi));
      Txy = Dxy*ZZMass / ZZPt;

      if (!tc[tt]->GetBranchStatus("D_bkg")) D_bkg = p0plus_VAJHU*p0plus_m4l / (p0plus_VAJHU*p0plus_m4l + bkg_VAMCFM*bkg_m4l);

      tout
        << RunNumber << "\t"
        << LumiNumber << "\t"
        << EventNumber << "\t"

        << ZZMass << "\t"
        << ZZPt << "\t"

        << Lep1SIP << "\t"
        << Lep2SIP << "\t"
        << Lep3SIP << "\t"
        << Lep4SIP << "\t"

        << Lep1_Z1SIP << "\t"
        << Lep2_Z1SIP << "\t"
        << Lep3_Z1SIP << "\t"
        << Lep4_Z1SIP << "\t"
        << KalmanCandVtx_chi2 << "\t"


        << BeamPosX*10000. << "\t"
        << BeamPosY*10000. << "\t"
        << BeamPosZ*10000. << "\t"
        << KalmanCandVtx_x*10000. << "\t"
        << KalmanCandVtx_y*10000. << "\t"
        << KalmanCandVtx_y*10000. << "\t"
        << OfflinePrimaryVtx_x*10000. << "\t"
        << OfflinePrimaryVtx_y*10000. << "\t"
        << OfflinePrimaryVtx_z*10000. << "\t"
        << (OfflinePrimaryVtx_ndof+3.)/2. << "\t"

        << Txy*10000. << "\t"
        << Dxy*10000. << "\t"
        << Txy_BS*10000. << "\t"
        << D_bkg << "\t"

        << endl;

      ctrEvents++;
    }
    cout << "Number of events that passed the cut: " << ctrEvents << endl;
    tout << "\nTotal: " << ctrEvents << endl;
  }

  tout.close();

	for (int tt = 0; tt < ntrees; tt++){
		delete tc[tt];
	}
  delete tBeam;
}

void listEfficiencies_single(int folder, int erg_tev){
  TString erg_dir;
  erg_dir.Form("LHC_%iTeV", erg_tev);
  TString comstring;
  comstring.Form("%iTeV", erg_tev);
  int EnergyIndex=0; if (erg_tev==8) EnergyIndex=1;

  TString processName[8] ={
    "CTau0", "CTau100", "CTau500", "CTau1000",
    "qqZZ", "ggZZ", "CR", "data"
  };

  TString cinput_common = user_dir_hep + "Analysis/Auxiliary/";

  cout << endl;
  cout << comstring << "\t" << user_folder[folder] << endl;
  cout << "Process\t";
  for (int c=0; c<5; c++){
    if (c==1) continue;
    cout << cutLabel[c] << '\t';
  }
  cout << endl;

  for (int p=0; p<6; p++){
    TString INPUT_NAME = "LifetimeKD_RelativeSIPYields_";
    INPUT_NAME.Append(processName[p]);
    INPUT_NAME.Append("_");
    INPUT_NAME.Append(user_folder[folder]);
    INPUT_NAME.Append("_");
    INPUT_NAME.Append(comstring);
    INPUT_NAME.Append("_m4l105.6_140.6.root");
    INPUT_NAME.Prepend(cinput_common);
    TFile* finput = new TFile(INPUT_NAME, "read");
    TH1F* hyield = (TH1F*)finput->Get(Form("h%s", processName[p].Data()));
    cout << processName[p] << '\t';
    for (int c=1; c<6; c++){
      if (c==2) continue;
      else cout << hyield->GetBinContent(c) << '\t';
    }
    cout << endl;
    delete hyield;
    finput->Close();
  }

}