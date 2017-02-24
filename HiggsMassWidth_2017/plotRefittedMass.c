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
#include "TLorentzVector.h"
#include "TPaveText.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TColor.h"
#include "TDirectory.h"
#include "./data/Samples.h"

using namespace std;
using namespace ROOT::Math;

enum{
  vTrueM4l,
  vObsM4l,
  vRefitKinZM4l,
  vRefitHMCM4l,

  vObsM4lErr,
  vPrefitKinZM4lErr,
  vPrefitHMCM4lErr,
  vRefitKinZM4lErr,
  vRefitHMCM4lErr,

  vObsTrueM4lDiffOverErr,
  vPrefitKinZTrueM4lDiffOverErr,
  vPrefitHMCTrueM4lDiffOverErr,
  vRefitKinZTrueM4lDiffOverErr,
  vRefitHMCTrueM4lDiffOverErr,

  vRefitKinZObsM4lDiffOverPrefitErr,
  vRefitHMCObsM4lDiffOverPrefitErr,

  nVariables
};
TString varNames[nVariables] ={
  "True_m4l",
  "Obs_m4l",
  "KinZ_m4l",
  "HMC_m4l",

  "Obs_m4lErr",
  "KinZ_m4lErrPre",
  "HMC_m4lErrPre",
  "KinZ_m4lErr",
  "HMC_m4lErr",

  "ObsTrue_m4lDiffOverErr",
  "KinZTrue_m4lDiffOverErrPre",
  "HMCTrue_m4lDiffOverErrPre",
  "KinZTrue_m4lDiffOverErr",
  "HMCTrue_m4lDiffOverErr",

  "KinZObs_m4lDiffOverPrefitErr",
  "HMCObs_m4lDiffOverPrefitErr"
};
TString varTitles[nVariables] ={
  "True",
  "Observed",
  "KinZFitter",
  "HMC",
  "Observed",
  "KinZFitter (Prefit)",
  "HMC (Prefit)",
  "KinZFitter",
  "HMC",
  "Observed",
  "KinZFitter (Prefit)",
  "HMC (Prefit)",
  "KinZFitter",
  "HMC",
  "KinZFitter",
  "HMC"
};

const int nM4lBins=50;
float m4lRange[nSamples][2] ={
  //{ 100, 200 },
  { 95, 130 },
  { 105, 140 },
  { 130, 170 },
  { 500, 1000 }
};
const int nM4lErrBins=50;
float m4lErrRange[nSamples][2] ={
  //{ 0, 6 },
  { 0, 6 },
  { 0, 6 },
  { 0, 6 },
  { 0, 20 }
};
const int nM4lTrueDiffOverErrBins=50;
float m4lTrueDiffOverErrRange[nSamples][2] ={
  //{ -7, 7 },
  { -7, 7 },
  { -7, 7 },
  { -7, 7 },
  { -7, 7 }
};
const int nM4lObsDiffOverPrefitErrBins=50;
float m4lObsDiffOverPrefitErrRange[nSamples][2] ={
  //{ -2, 2 },
  { -2, 2 },
  { -2, 2 },
  { -2, 2 },
  { -2, 2 }
};

void plotRefittedMass_single(int iSample){
  gROOT->ProcessLine(".x tdrstyle.cc");

  TString TREE_NAME = "ZZTree/candTree";
  if (strSamples[iSample].second=="DY") TREE_NAME = "CRZLLTree/candTree";

  Float_t overallEventWeight=1;
  Float_t KFactor_QCD_ggZZ_Nominal=1;
  Float_t KFactor_EW_qqZZ=1;
  Float_t KFactor_QCD_qqZZ_M=1;
  Float_t xsec=1;
  Short_t ZZsel;
  Float_t GenHMass=0;
  Float_t ZZMass;
  Float_t ZZMassRefit;
  Float_t ZZMassHMCRefit;
  Float_t ZZMassErrCorr;
  Float_t ZZMassUnrefitErr;
  Float_t ZZMassHMCUnrefitErr;
  Float_t ZZMassRefitErr;
  Float_t ZZMassHMCRefitErr;
  Float_t* ZZMassRef[vRefitHMCM4l+1]={
    &GenHMass,
    &ZZMass,
    &ZZMassRefit,
    &ZZMassHMCRefit
  };
  Float_t* ZZMassErrRef[vRefitHMCM4lErr-vObsM4lErr+1]={
    &ZZMassErrCorr,
    &ZZMassUnrefitErr,
    &ZZMassHMCUnrefitErr,
    &ZZMassRefitErr,
    &ZZMassHMCRefitErr
  };
  Short_t Z1Flav;
  Short_t Z2Flav;

  TString OUTPUT_NAME = "M4L_Comparison_";
  OUTPUT_NAME.Append(Form("%s", strSamples[iSample].second.Data()));
  //OUTPUT_NAME.Append(Form("_%s_", channame[ichan]));
  OUTPUT_NAME.Append("_13TeV");
  OUTPUT_NAME.Append(".root");

  TString cinput_common = user_dir;
  TString coutput_common = "./";

  TString coutput = coutput_common; coutput += OUTPUT_NAME;
  TFile* foutput = TFile::Open(coutput, "recreate");

  TString cinput = cinput_common; cinput += strSamples[iSample].first;
  TFile* fin = TFile::Open(cinput, "read");
  TTree* tin = (TTree*)fin->Get(TREE_NAME);

  if (tin->GetBranchStatus("KFactor_QCD_ggZZ_Nominal")) tin->SetBranchAddress("KFactor_QCD_ggZZ_Nominal", &KFactor_QCD_ggZZ_Nominal);
  if (tin->GetBranchStatus("KFactor_EW_qqZZ")){
    tin->SetBranchAddress("KFactor_EW_qqZZ", &KFactor_EW_qqZZ);
    tin->SetBranchAddress("KFactor_QCD_qqZZ_M", &KFactor_QCD_qqZZ_M);
  }
  tin->SetBranchAddress("Z1Flav", &Z1Flav);
  tin->SetBranchAddress("Z2Flav", &Z2Flav);
  tin->SetBranchAddress("ZZsel", &ZZsel);
  tin->SetBranchAddress("ZZMass", &ZZMass);
  tin->SetBranchAddress("ZZMassRefit", &ZZMassRefit);
  tin->SetBranchAddress("ZZMassHMCRefit", &ZZMassHMCRefit);
  tin->SetBranchAddress("ZZMassErrCorr", &ZZMassErrCorr);
  tin->SetBranchAddress("ZZMassUnrefitErr", &ZZMassUnrefitErr);
  tin->SetBranchAddress("ZZMassHMCUnrefitErr", &ZZMassHMCUnrefitErr);
  tin->SetBranchAddress("ZZMassRefitErr", &ZZMassRefitErr);
  tin->SetBranchAddress("ZZMassHMCRefitErr", &ZZMassHMCRefitErr);
  if (tin->GetBranchStatus("GenHMass")){
    tin->SetBranchAddress("GenHMass", &GenHMass);
    tin->SetBranchAddress("overallEventWeight", &overallEventWeight);
    tin->SetBranchAddress("xsec", &xsec);
  }

  TH1F* hvars[nChannels][nVariables] ={ { 0 } };
  for (int ichan = 0; ichan < nChannels; ichan++){
    TString hname_core = "h_";
    for (int v = 0; v < nVariables; v++){
      int nbins;
      float xmin, xmax;
      TString xtitle;

      if (v<=vRefitHMCM4l){
        nbins = nM4lBins;
        xmin = m4lRange[iSample][0];
        xmax = m4lRange[iSample][1];
        xtitle = "m_{4l} (GeV)";
      }
      else if (v<=vRefitHMCM4lErr){
        nbins = nM4lErrBins;
        if (strSamples[iSample].second=="ggH750") nbins *= 2;
        xmin = m4lErrRange[iSample][0];
        xmax = m4lErrRange[iSample][1];
        xtitle = "#sigma_{m4l} (GeV)";
      }
      else if (v<=vRefitHMCTrueM4lDiffOverErr){
        nbins = nM4lTrueDiffOverErrBins;
        xmin = m4lTrueDiffOverErrRange[iSample][0];
        xmax = m4lTrueDiffOverErrRange[iSample][1];
        if (strSamples[iSample].second=="ggH750"){ xmax *= 1.5; xmin *= 1.5; }
        xtitle = "(m_{4l} - m^{true}_{4l})/#sigma_{m4l}";
      }
      else if (v<=vRefitHMCObsM4lDiffOverPrefitErr){
        nbins = nM4lObsDiffOverPrefitErrBins;
        xmin = m4lObsDiffOverPrefitErrRange[iSample][0];
        xmax = m4lObsDiffOverPrefitErrRange[iSample][1];
        if (strSamples[iSample].second=="ggH750"){ xmax *= 1.5; xmin *= 1.5; }
        xtitle = "(m_{4l} - m^{obs}_{4l})/#sigma^{prefit}_{m4l}";
      }
      float binwidth = (xmax-xmin)/((float)nbins);

      TString hname = hname_core;
      hname.Append(varNames[v]);
      hname.Append("_");
      hname.Append(channame[ichan]);
      hvars[ichan][v] = new TH1F(hname, "", nbins, xmin, xmax);
      hvars[ichan][v]->Sumw2();
      hvars[ichan][v]->SetXTitle(xtitle);
      hvars[ichan][v]->SetYTitle("Rate / bin");
      hvars[ichan][v]->GetXaxis()->SetNdivisions(505);
      hvars[ichan][v]->GetXaxis()->SetLabelFont(42);
      hvars[ichan][v]->GetXaxis()->SetLabelOffset(0.007);
      hvars[ichan][v]->GetXaxis()->SetLabelSize(0.04);
      hvars[ichan][v]->GetXaxis()->SetTitleSize(0.06);
      hvars[ichan][v]->GetXaxis()->SetTitleOffset(0.9);
      hvars[ichan][v]->GetXaxis()->SetTitleFont(42);
      hvars[ichan][v]->GetYaxis()->SetNdivisions(505);
      hvars[ichan][v]->GetYaxis()->SetLabelFont(42);
      hvars[ichan][v]->GetYaxis()->SetLabelOffset(0.007);
      hvars[ichan][v]->GetYaxis()->SetLabelSize(0.04);
      hvars[ichan][v]->GetYaxis()->SetTitleSize(0.06);
      hvars[ichan][v]->GetYaxis()->SetTitleOffset(1.1);
      hvars[ichan][v]->GetYaxis()->SetTitleFont(42);

      hvars[ichan][v]->SetLineWidth(2);
      if (v==vTrueM4l){
        hvars[ichan][v]->SetLineStyle(1);
        hvars[ichan][v]->SetLineColor(kBlack);
        hvars[ichan][v]->SetMarkerColor(kBlack);
      }
      else if (v==vObsM4l || v==vObsM4lErr || v==vObsTrueM4lDiffOverErr){
        hvars[ichan][v]->SetLineStyle(1);
        hvars[ichan][v]->SetLineColor(kBlue);
        hvars[ichan][v]->SetMarkerColor(kBlue);
      }
      else if (v==vPrefitKinZM4lErr || v==vPrefitKinZTrueM4lDiffOverErr){
        hvars[ichan][v]->SetLineStyle(1);
        hvars[ichan][v]->SetLineColor(kCyan-7);
        hvars[ichan][v]->SetMarkerColor(kCyan-7);
      }
      else if (v==vPrefitHMCM4lErr || v==vPrefitHMCTrueM4lDiffOverErr){
        hvars[ichan][v]->SetLineStyle(1);
        hvars[ichan][v]->SetLineColor(kRed);
        hvars[ichan][v]->SetMarkerColor(kRed);
      }
      else if (v==vRefitKinZM4l || v==vRefitKinZM4lErr || v==vRefitKinZTrueM4lDiffOverErr || v==vRefitKinZObsM4lDiffOverPrefitErr){
        hvars[ichan][v]->SetLineStyle(2);
        hvars[ichan][v]->SetLineColor(kGreen+2);
        hvars[ichan][v]->SetMarkerColor(kGreen+2);
      }
      else if (v==vRefitHMCM4l || v==vRefitHMCM4lErr || v==vRefitHMCTrueM4lDiffOverErr || v==vRefitHMCObsM4lDiffOverPrefitErr){
        hvars[ichan][v]->SetLineStyle(2);
        hvars[ichan][v]->SetLineColor(kViolet);
        hvars[ichan][v]->SetMarkerColor(kViolet);
      }
    }
  }

  double nEntries[nChannels][nVariables]={ { 0 } };
  double sumWgt[nChannels][nVariables]={ { 0 } };
  double mean[nChannels][nVariables]={ { 0 } };
  double rms[nChannels][nVariables]={ { 0 } };

  for (int ev=0; ev<tin->GetEntries(); ev++){
    tin->GetEntry(ev);

    int ichan=-1;
    if (abs(Z1Flav)==169 && abs(Z2Flav)==169) ichan=0;
    else if (abs(Z1Flav)==121 && abs(Z2Flav)==121) ichan=1;
    else if (abs(Z1Flav)*abs(Z2Flav)==121*169) ichan=2;
    else continue;

    double weight = overallEventWeight*xsec*KFactor_QCD_ggZZ_Nominal*KFactor_EW_qqZZ*KFactor_QCD_qqZZ_M;
    for (int v = 0; v < nVariables; v++){
      Float_t varVal=0;

      if (v<=vRefitHMCM4l){
        int vMass = v;
        varVal = *(ZZMassRef[vMass]);
        //if (varVal<m4lRange[iSample][0] || varVal>=m4lRange[iSample][1]) continue;
      }
      else if (v<=vRefitHMCM4lErr){
        int vMass = v-vRefitHMCM4l-1;
        varVal = *(ZZMassErrRef[vMass]);
        //if (varVal<m4lErrRange[iSample][0] || varVal>=m4lErrRange[iSample][1]) continue;
      }
      else if (v<=vRefitHMCTrueM4lDiffOverErr){
        int vMass = 1;
        if (v>=vRefitKinZTrueM4lDiffOverErr) vMass = v - vRefitKinZTrueM4lDiffOverErr + vRefitKinZM4l;
        int vErr = v-vRefitHMCM4lErr-1;
        varVal = (*(ZZMassRef[vMass]) - *(ZZMassRef[vTrueM4l]))/(*(ZZMassErrRef[vErr]));
        //if (varVal<m4lTrueDiffOverErrRange[iSample][0] || varVal>=m4lTrueDiffOverErrRange[iSample][1]) continue;
      }
      else if (v<=vRefitHMCObsM4lDiffOverPrefitErr){
        int vMass = v - vRefitKinZObsM4lDiffOverPrefitErr + vRefitKinZM4l;
        int vErr = v-vRefitKinZObsM4lDiffOverPrefitErr+1;
        varVal = (*(ZZMassRef[vMass]) - *(ZZMassRef[vObsM4l]))/(*(ZZMassErrRef[vErr]));
        //if (varVal<m4lTrueDiffOverErrRange[iSample][0] || varVal>=m4lTrueDiffOverErrRange[iSample][1]) continue;
      }
      if (varVal<hvars[ichan][v]->GetXaxis()->GetBinLowEdge(1) || varVal>=hvars[ichan][v]->GetXaxis()->GetBinUpEdge(hvars[ichan][v]->GetNbinsX())) continue;

      nEntries[ichan][v]++;
      sumWgt[ichan][v] += weight;
      mean[ichan][v] += weight*varVal;
      rms[ichan][v] += weight*pow(varVal, 2);

      hvars[ichan][v]->Fill(varVal, weight);
    }
  }

  for (int ichan = 0; ichan < nChannels; ichan++){
    for (int ivt=0; ivt<4; ivt++){
      TString canvasname;
      int beginVar, endVar;

      if (ivt==0){
        canvasname = "cM4L_Comparison_";
        beginVar=(strSamples[iSample].second=="ggH750" ? 0 : 1);
        endVar=vRefitHMCM4l;
      }
      else if (ivt==1){
        canvasname = "cM4LErr_Comparison_";
        beginVar=vObsM4lErr;
        endVar=vRefitHMCM4lErr;
      }
      else if (ivt==2){
        canvasname = "cM4LTrueDiffOverErr_Comparison_";
        beginVar=vObsTrueM4lDiffOverErr;
        endVar=vRefitHMCTrueM4lDiffOverErr;
      }
      else if (ivt==3){
        canvasname = "cM4LObsDiffOverPrefitErr_Comparison_";
        beginVar=vRefitKinZObsM4lDiffOverPrefitErr;
        endVar=vRefitHMCObsM4lDiffOverPrefitErr;
      }
      canvasname.Append(channame[ichan]);
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

      TLegend* ll;
      float lxmin = 0.18, lxwidth = 0.30;
      float lymax = 0.9, lywidth = 0.2;
      float lxmax = lxmin + lxwidth;
      float lymin = lymax - lywidth;
      ll = new TLegend(lxmin, lymin, lxmax, lymax);
      ll->SetBorderSize(0);
      ll->SetTextFont(42);
      ll->SetTextSize(0.03);
      ll->SetLineColor(1);
      ll->SetLineStyle(1);
      ll->SetLineWidth(1);
      ll->SetFillColor(0);
      ll->SetFillStyle(0);

      double maxplot=0;
      for (int v=beginVar; v<=endVar; v++){
        for (int bin = 1; bin <= hvars[ichan][v]->GetNbinsX(); bin++){
          double bincontent = hvars[ichan][v]->GetBinContent(bin);
          double binerror = hvars[ichan][v]->GetBinError(bin);
          bincontent += binerror;
          if (bincontent>maxplot) maxplot = bincontent;
        }
      }
      maxplot *= 1.3;

      for (int v=beginVar; v<=endVar; v++){
        hvars[ichan][v]->GetYaxis()->SetRangeUser(0, maxplot);
        TString strlabel = varTitles[v];

        cout << "NSample " << strSamples[iSample].second << " variable " << varNames[v] << " for channel " << channame[ichan] << endl;
        cout << "N: " << nEntries[ichan][v] << endl;
        cout << "sumWgt: " << sumWgt[ichan][v] << endl;
        float xbar = mean[ichan][v]/sumWgt[ichan][v];
        float sxbar = sqrt(sumWgt[ichan][v]*rms[ichan][v]-pow(mean[ichan][v], 2))/sumWgt[ichan][v];
        cout << "Mean: " << xbar << endl;
        cout << "RMS: " << sxbar << endl;

        if (
          v<=vRefitHMCM4l
          ||
          (v>=vObsTrueM4lDiffOverErr && v<=vRefitHMCTrueM4lDiffOverErr)
          ||
          (v>=vRefitKinZObsM4lDiffOverPrefitErr && v<=vRefitHMCObsM4lDiffOverPrefitErr)
          ){
          TString strMeanRMS = Form(" %.2f %.2f %.0f", xbar, sxbar, nEntries[ichan][v]);
          strlabel += strMeanRMS;
        }
        ll->AddEntry(hvars[ichan][v], strlabel, "l");

        if (v==beginVar) hvars[ichan][v]->Draw("hist");
        else hvars[ichan][v]->Draw("histsame");
      }
      ll->Draw("same");

      float pt_xmin = 0.71, pt_xwidth = 0.21;
      float pt_ymax = 0.92, pt_ywidth = 0.08;
      float pt_xmax = pt_xmin + pt_xwidth;
      float pt_ymin = pt_ymax - pt_ywidth;
      TPaveText *pt10 = new TPaveText(pt_xmin, pt_ymin, pt_xmax, pt_ymax, "brNDC");
      pt10->SetBorderSize(0);
      pt10->SetTextAlign(12);
      pt10->SetTextSize(0.04);
      pt10->SetFillStyle(0);
      pt10->SetTextFont(42);
      TText* text10 = pt10->AddText(0.01, 0.01, chanlabel[ichan]);
      pt10->Draw();

      TString plotDir = coutput_common;
      plotDir.Append("plots/");
      plotDir.Append(Form("%s/", strSamples[iSample].second.Data()));
      TString plotDir_1D=plotDir;
      plotDir_1D.Append("1D/");
      TString mkdirCommand = "mkdir -p ";
      TString mkdirCommand_1D = mkdirCommand;
      mkdirCommand_1D.Append(plotDir_1D);
      gSystem->Exec(mkdirCommand_1D);

      cc->RedrawAxis();
      cc->Modified();
      cc->Update();
      foutput->WriteTObject(cc);
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
      //cc->SaveAs(canvasname_eps);
      //cc->SaveAs(canvasname_png);
      //cc->SaveAs(canvasname_root);
      //cc->SaveAs(canvasname_c);

      delete pt10;
      delete ll;
      cc->Close();
    }

    for (int v = 0; v < nVariables; v++){
      foutput->WriteTObject(hvars[ichan][v]);
      delete hvars[ichan][v];
    }
  }
  foutput->Close();
}
