#include <iostream>
#include <cmath>
#include <string>
#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include "TF1.h"
#include "TMatrixD.h"
#include "TVectorD.h"
#include "TString.h"
#include "TSpline.h"
#include "TGraph.h"
#include "TGraphAsymmErrors.h"
#include "TProfile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TPaveText.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TColor.h"
#include "TDirectory.h"

using namespace std;


char location_primaryTrees[] = "/afs/cern.ch/work/l/lfeng/public_write/usarica/Production/PrimaryTrees/HZZ_out/";
const int nProdModes=5;
char* sampleName[nProdModes] = {
	"powheg15jhuGenV3-0PMH125.6",
  "VBFH125",
  "WH125",
  "ZH125",
  "ttH125" 
};
char* prodName[nProdModes] = {
	"ggH",
  "VBFH125",
  "WH125",
  "ZH125",
  "ttH125"
};
char* prodTitle[nProdModes]={
  "ggH (125.6 GeV)", "VBF H (125 GeV)", "WH (125 GeV)", "ZH (125 GeV)", "t#bar{t}H (125 GeV)"
};


void makeGenLevel_pTRatios(int erg_tev){
  char TREE_NAME[]="candTree";
  TString OUTPUT_NAME = Form("compareProductionModes_pT_%iTeV.root", erg_tev);
  char csuffix[] = "_pToverMZZ";
  char cxtitle[] = "p_{T} / m_{4l}";

  TFile* foutput = new TFile(OUTPUT_NAME, "recreate");
  TH1F* h_pt[nProdModes];
  TH1F* hratio_pt[nProdModes-1];
  TProfile* hpr_pt[nProdModes-1];
  TGraphAsymmErrors* tgratio_pt[nProdModes-1];
  TF1* fratio_pt[nProdModes];
  double xmin = 0;
  double xmax = 6;
  const int nbins = 23;
  double xbins[nbins+1];
  for (int bin=0; bin<=nbins; bin++){
    if (bin==0) xbins[bin] = xmin;
    else if (xbins[bin-1]<1.2) xbins[bin] = xbins[bin-1] + 0.1; // 12
    else if (xbins[bin-1]<2.4) xbins[bin] = xbins[bin-1] + 0.2; // 6
    else if (xbins[bin-1]<4) xbins[bin] = xbins[bin-1] + 0.4; // 4
    else if (xbins[bin-1]<xmax) xbins[bin] = xbins[bin-1] + 2.0; // 1
  }
  TString cytitle = "Rate / bin";
  for (int p = 0; p < nProdModes; p++){
    TString hname = Form("%s%s", prodName[p], csuffix);
//    TString htitle = Form("%s %i TeV", prodName[p], erg_tev);
    h_pt[p] = new TH1F(hname, "", nbins, xbins);
    h_pt[p]->Sumw2();
    h_pt[p]->GetXaxis()->SetTitle(cxtitle);
    h_pt[p]->GetYaxis()->SetTitle(cytitle);

    if (p>0){
      hname = Form("%s%s_profile", prodName[p], csuffix);
      hpr_pt[p-1] = new TProfile(hname, "", nbins, xbins);
      hpr_pt[p-1]->Sumw2();
    }
  }
  int genFinalState = -1;
  float GenHMass = 0;
  float GenHPt = 0;
  TString cinput_common = location_primaryTrees;
  if (erg_tev==7) cinput_common.Append("141217/");
  else if (erg_tev==8) cinput_common.Append("141021/");
  for (int p = 0; p < nProdModes; p++){
    TString cinput = cinput_common;
//    cinput.Append(Form("%s/ZZ4lAnalysis_%s%s", sampleName[p], sampleName[p], ".root"));
    cinput.Append(Form("ZZ4lAnalysis_%s%s", sampleName[p], ".root"));
    TFile* finput = new TFile(cinput, "read");
    cout << "Production " << prodName[p] << endl;
    TTree* tin[3] ={
      (TTree*)finput->Get(Form("ZZ4muTree/%s", TREE_NAME)),
      (TTree*)finput->Get(Form("ZZ4eTree/%s", TREE_NAME)),
      (TTree*)finput->Get(Form("ZZ2e2muTree/%s", TREE_NAME))
    };
    for (int f = 0; f < 3; f++){
      cout << "Channel " << f << endl;
      cout << "Nentries: " << tin[f]->GetEntries() << endl;
      tin[f]->SetBranchAddress("genFinalState", &genFinalState);
      tin[f]->SetBranchAddress("GenHMass", &GenHMass);
      tin[f]->SetBranchAddress("GenHPt", &GenHPt);
      for (int ev = 0; ev < tin[f]->GetEntries(); ev++){
        //			for (int ev = 0; ev < 1000; ev++){
        genFinalState = -1;
        GenHMass = 0;
        GenHPt = 0;
        tin[f]->GetEntry(ev);
        if (f != genFinalState) continue;
        double kd = GenHPt / GenHMass;
        h_pt[p]->Fill(kd);
        if (p==0){
          for (int pp=0; pp<nProdModes-1; pp++) hpr_pt[pp]->Fill(kd, kd);
        }
        else hpr_pt[p-1]->Fill(kd, kd);
      }
    }
    finput->Close();
    cout << "File closed" << endl;
  }
  foutput->cd();
/*  double zhxsec, whxsec;
  if (erg_tev == 7){ whxsec = 0.5688; zhxsec = 0.3299; }
  else if (erg_tev == 8){ whxsec = 0.6931; zhxsec = 0.4091; }
  h_pt[nProdModes]->Add(h_pt[nProdModes-1], whxsec);
  h_pt[nProdModes]->Add(h_pt[nProdModes-2], zhxsec);
*/  

  double max_plot = 0;
  for (int p = 0; p < nProdModes; p++){
    h_pt[p]->Scale(1. / h_pt[p]->Integral(0, h_pt[p]->GetNbinsX() + 1));

    h_pt[p]->GetYaxis()->SetRangeUser(0, 1);
    h_pt[p]->SetLineWidth(2);
    h_pt[p]->SetLineStyle(1);
    if (p==0) h_pt[p]->SetLineColor(kBlack);
    else if (p==3) h_pt[p]->SetLineColor(kRed);
    else if (p==4) h_pt[p]->SetLineColor(kGreen+2);
    else if (p==2) h_pt[p]->SetLineColor(kBlue);
    else if (p==1) h_pt[p]->SetLineColor(kViolet);
    if (p==0) h_pt[p]->SetMarkerColor(kBlack);
    else if (p==3) h_pt[p]->SetMarkerColor(kRed);
    else if (p==4) h_pt[p]->SetMarkerColor(kGreen+2);
    else if (p==2) h_pt[p]->SetMarkerColor(kBlue);
    else if (p==1) h_pt[p]->SetMarkerColor(kViolet);
    cout << "Writing hpt at " << p << endl;
    foutput->WriteTObject(h_pt[p]);
    if (p > 0){
      cout << "Creating ratio at " << p << endl;
      hratio_pt[p - 1] = (TH1F*)h_pt[p]->Clone(Form("%s_ratio", h_pt[p]->GetName()));
      hratio_pt[p - 1]->Divide(h_pt[0]);
      hratio_pt[p - 1]->GetYaxis()->SetTitle("Ratio to ggH");
      hratio_pt[p - 1]->SetTitle(Form("%s / %s", prodName[p], h_pt[0]->GetTitle()));

      double xx[nbins];
      double yy[nbins];
      double xx_up[nbins];
      double yy_up[nbins];
      double xx_dn[nbins];
      double yy_dn[nbins];
      int nAcc=0;
      for (int bin=1; bin<=nbins; bin++){
        double bincontent = hratio_pt[p - 1]->GetBinContent(bin);
        double binerror = hratio_pt[p - 1]->GetBinError(bin);
        if (bincontent == 0) continue;
        else{
          xx[nAcc] = hpr_pt[p-1]->GetBinContent(bin);
          xx_up[nAcc] = hpr_pt[p-1]->GetBinError(bin);
          xx_dn[nAcc] = hpr_pt[p-1]->GetBinError(bin);
          yy[nAcc] = bincontent;
          yy_up[nAcc] = binerror;
          yy_dn[nAcc] = binerror;
          nAcc++;
          max_plot = max(max_plot, (bincontent+binerror));
        }
      }
      tgratio_pt[p-1] = new TGraphAsymmErrors(nAcc, xx, yy, xx_dn, xx_up, yy_dn, yy_up);
      tgratio_pt[p-1]->SetName(Form("tg_%s_ratio", h_pt[p]->GetName()));
      tgratio_pt[p-1]->GetYaxis()->SetTitle("Ratio to ggH");
      tgratio_pt[p-1]->GetXaxis()->SetTitle(cxtitle);
      tgratio_pt[p-1]->GetYaxis()->SetRangeUser(0, 30);
      tgratio_pt[p-1]->SetLineWidth(2);
      tgratio_pt[p-1]->SetLineStyle(1);
      if (p==3) tgratio_pt[p-1]->SetMarkerColor(kRed);
      else if (p==4) tgratio_pt[p-1]->SetMarkerColor(kGreen+2);
      else if (p==2) tgratio_pt[p-1]->SetMarkerColor(kBlue);
      else if (p==1) tgratio_pt[p-1]->SetMarkerColor(kViolet);
      if (p==3) tgratio_pt[p-1]->SetLineColor(kRed);
      else if (p==4) tgratio_pt[p-1]->SetLineColor(kGreen+2);
      else if (p==2) tgratio_pt[p-1]->SetLineColor(kBlue);
      else if (p==1) tgratio_pt[p-1]->SetLineColor(kViolet);

      TF1* fitratio;
      if (p != 4){
        fitratio = new TF1(Form("%s_fit", hratio_pt[p - 1]->GetName()), "([0]-[1]*exp(-pow(x/[2],2)))*exp(-x/[3])", hratio_pt[p - 1]->GetXaxis()->GetXmin(), hratio_pt[p - 1]->GetXaxis()->GetXmax());
        if (p==1) fitratio->SetParameters(8.76, 8.69, 0.98, 4.3);
        else fitratio->SetParameters(8.76, 8.69, 0.98, -24);
      }
      else{
        fitratio = new TF1(Form("%s_fit", hratio_pt[p - 1]->GetName()), "([0]-[1]*exp(-pow(x/[2],2)))", hratio_pt[p - 1]->GetXaxis()->GetXmin(), hratio_pt[p - 1]->GetXaxis()->GetXmax());
        fitratio->SetParameters(4.2, 4.0, 0.7);
      }
      fitratio->SetParameters(6, 5.77, 1, 5);
      fitratio->SetTitle("");
      fitratio->SetLineColor(hratio_pt[p - 1]->GetLineColor());
      fitratio->SetLineStyle(7);
      fitratio->SetLineWidth(3);
      hratio_pt[p - 1]->GetYaxis()->SetRangeUser(1e-3, 1000);
      tgratio_pt[p - 1]->Fit(fitratio, "N");
      cout << "Writing ratio fit at " << p << endl;
      foutput->WriteTObject(fitratio);
      cout << "Writing ratio at " << p << endl;
      foutput->WriteTObject(hratio_pt[p - 1]);
      foutput->WriteTObject(tgratio_pt[p - 1]);
      fratio_pt[p - 1] = fitratio;
    }
  }
  foutput->cd();
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
  if (erg_tev==8) text = pt->AddText(0.837, 0.45, "#font[42]{               8 TeV}");
  else if (erg_tev==7) text = pt->AddText(0.837, 0.45, "#font[42]{               7 TeV}");
//  if (erg_tev==7) text = pt->AddText(0.837, 0.45, "#font[42]{5.1 fb^{-1} (7 TeV)}");
  text->SetTextSize(0.0315);


  TString canvasname_2D = Form("cCompare_SignalProductionMC_AllChannels_%iTeV", erg_tev);
  canvasname_2D.Append(Form("_%s_ratio", csuffix));
  TCanvas* c2D = new TCanvas(canvasname_2D, "", 8, 30, 800, 800);
  c2D->cd();
  gStyle->SetOptStat(0);
  c2D->SetFillColor(0);
  c2D->SetBorderMode(0);
  c2D->SetBorderSize(2);
  c2D->SetTickx(1);
  c2D->SetTicky(1);
  //		c2D->SetLogy();
  c2D->SetLeftMargin(0.17);
  c2D->SetRightMargin(0.05);
  c2D->SetTopMargin(0.07);
  c2D->SetBottomMargin(0.13);
  c2D->SetFrameFillStyle(0);
  c2D->SetFrameBorderMode(0);
  c2D->SetFrameFillStyle(0);
  c2D->SetFrameBorderMode(0);

  TLegend *l2D = new TLegend(0.20, 0.67, 0.58, 0.90);
  l2D->SetBorderSize(0);
  l2D->SetTextFont(42);
  l2D->SetTextSize(0.03);
  l2D->SetLineColor(1);
  l2D->SetLineStyle(1);
  l2D->SetLineWidth(1);
  l2D->SetFillColor(0);
  l2D->SetFillStyle(0);
/*
  for (int p = 0; p < nProdModes-1; p++){
    hratio_pt[p]->SetTitle("");
    hratio_pt[p]->GetXaxis()->SetNdivisions(505);
    hratio_pt[p]->GetXaxis()->SetLabelFont(42);
    hratio_pt[p]->GetXaxis()->SetLabelOffset(0.007);
    hratio_pt[p]->GetXaxis()->SetLabelSize(0.04);
    hratio_pt[p]->GetXaxis()->SetTitleSize(0.06);
    hratio_pt[p]->GetXaxis()->SetTitleOffset(0.9);
    hratio_pt[p]->GetXaxis()->SetTitleFont(42);
    hratio_pt[p]->GetYaxis()->SetNdivisions(505);
    hratio_pt[p]->GetYaxis()->SetLabelFont(42);
    hratio_pt[p]->GetYaxis()->SetLabelOffset(0.007);
    hratio_pt[p]->GetYaxis()->SetLabelSize(0.04);
    hratio_pt[p]->GetYaxis()->SetTitleSize(0.06);
    hratio_pt[p]->GetYaxis()->SetTitleOffset(1.1);
    hratio_pt[p]->GetYaxis()->SetTitleFont(42);
    hratio_pt[p]->GetYaxis()->SetRangeUser(0, 30);
    hratio_pt[p]->GetXaxis()->SetRangeUser(0, 5);
    l2D->AddEntry(hratio_pt[p], prodTitle[p], "l");
    if (p == 0) hratio_pt[p]->Draw("e1p");
    else hratio_pt[p]->Draw("e1psame");
    fratio_pt[p]->Draw("csame");
  }
*/
  for (int p = 0; p < nProdModes-1; p++){
    tgratio_pt[p]->SetTitle("");
    tgratio_pt[p]->GetXaxis()->SetNdivisions(505);
    tgratio_pt[p]->GetXaxis()->SetLabelFont(42);
    tgratio_pt[p]->GetXaxis()->SetLabelOffset(0.007);
    tgratio_pt[p]->GetXaxis()->SetLabelSize(0.04);
    tgratio_pt[p]->GetXaxis()->SetTitleSize(0.06);
    tgratio_pt[p]->GetXaxis()->SetTitleOffset(0.9);
    tgratio_pt[p]->GetXaxis()->SetTitleFont(42);
    tgratio_pt[p]->GetYaxis()->SetNdivisions(505);
    tgratio_pt[p]->GetYaxis()->SetLabelFont(42);
    tgratio_pt[p]->GetYaxis()->SetLabelOffset(0.007);
    tgratio_pt[p]->GetYaxis()->SetLabelSize(0.04);
    tgratio_pt[p]->GetYaxis()->SetTitleSize(0.06);
    tgratio_pt[p]->GetYaxis()->SetTitleOffset(1.1);
    tgratio_pt[p]->GetYaxis()->SetTitleFont(42);
    tgratio_pt[p]->GetYaxis()->SetRangeUser(0, max_plot*1.25);
    tgratio_pt[p]->GetXaxis()->SetRangeUser(0, 6);
    l2D->AddEntry(tgratio_pt[p], prodTitle[p+1], "l");
    if (p == 0) tgratio_pt[p]->Draw("ae1p");
    else tgratio_pt[p]->Draw("e1psame");
    fratio_pt[p]->Draw("csame");
  }


  l2D->Draw("same");
  pt->Draw();

  c2D->RedrawAxis();
  c2D->Modified();
  c2D->Update();
  foutput->WriteTObject(c2D);

  c2D->Close();
  foutput->cd();

  for (int p = 0; p < nProdModes; p++){
    delete h_pt[p];
    if (p>0){
      delete tgratio_pt[p-1];
      delete hpr_pt[p-1];
      delete hratio_pt[p-1];
    }
  }
  foutput->Close();
}
