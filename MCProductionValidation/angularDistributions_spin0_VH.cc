#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include <iostream>
#include <cmath>
#include <string>
#include <algorithm>
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussModel.h"
#include "RooDecay.h"
#include "RooBMixDecay.h"
#include "RooCategory.h"
#include "RooBinning.h"
#include "RooPlot.h"
#include "TSystem.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TH1.h"
#include "TGaxis.h"
#include "TString.h"
#include "TChain.h"
#include "src/ScalarPdfFactory_VH.cc"


using namespace RooFit;
using namespace std;

void angularDistributions_spin0_VH(TString INPUT_NAME, double sqrts = 13, double g1Re=0, double g2Re=0, double g4Re=0, double g2Im=0, double g4Im=0, int isLeptonic=0, int do3D=1, int nbins=80){
  bool bdo3D = (do3D==0 ? false : true);
  const int nVars=6;
  float kd_vars[nVars];
  TString strKDs[nVars]={
    "GenHMass",
    "GenhelcosthetaV1", "GenhelcosthetaV2", "GenphistarV1",
    "Gencosthetastar", "Genhelphi"
  };
  for (int v=1; v<nVars; v++){
    if (isLeptonic==0) strKDs[v].Append("_VHhadronic");
    else if (isLeptonic==1) strKDs[v].Append("_VHleptonic");
  }

  RooRealVar* mzz = new RooRealVar(strKDs[0], "m_{ZZ} (GeV)", 125, 125.-0.02, 125.+0.02);
  RooRealVar* h1 = new RooRealVar(strKDs[1], "cos#theta_{V*}", -1, 1);
  RooRealVar* h2 = new RooRealVar(strKDs[2], "cos#theta_{V}", -1, 1);
  RooRealVar* Phi1 = new RooRealVar(strKDs[3], "#Phi_{V*}", -TMath::Pi(), TMath::Pi());
  RooRealVar* hs = new RooRealVar(strKDs[4], "cos#theta^{*}", -1, 1);
  RooRealVar* Phi = new RooRealVar(strKDs[5], "#Phi", -TMath::Pi(), TMath::Pi());
  RooRealVar* Sqrts;
  RooRealVar* mV;
  if (INPUT_NAME.Contains("ZH")) mV = new RooRealVar("mV", "mV", 91.1876);
  else mV = new RooRealVar("mV", "mV", 80.399);
  if (isLeptonic==1) Sqrts = new RooRealVar("GenDileptonVVMass", "GenDileptonVVMass", 460, 0, 14000);
  else Sqrts = new RooRealVar("GenDijetVVMass", "GenDijetVVMass", 460, 0, 14000);

  RooArgSet treeargs(*mzz, *h1, *h2, *Phi, *hs, *Phi1);
  RooRealVar* measurables[nVars]={ h1, h2, Phi, hs, Phi1, mzz };

  ScalarPdfFactory_VH* someHiggs = new ScalarPdfFactory_VH(h1, h2, hs, Phi, Phi1, mzz, mV, Sqrts, kRealImag_Gs, bdo3D, false);
  someHiggs->g1Val->setVal(g1Re);
  someHiggs->g2Val->setVal(g2Re);
  someHiggs->g3Val->setVal(0);
  someHiggs->g4Val->setVal(g4Re);
  someHiggs->g2ValIm->setVal(g2Im);
  someHiggs->g4ValIm->setVal(g4Im);
  someHiggs->makeParamsConst(true);

  string coutput_common = "/scratch0/hep/usarical/SpinWidthPaper_2015/CMSSW_6_1_1/src/Analysis/Validation/Plots/ProductionValidation/";
  string cinput = INPUT_NAME.Data();
  size_t lastSlash = cinput.find_last_of("/\\");
  string finName = cinput.substr(lastSlash+1);
  finName = finName.substr(0,finName.find(".root"));
  finName = finName + "/";
  string coutput = coutput_common + finName;
  string strCmd = "mkdir -p ";
  strCmd.append(coutput);
  gSystem->Exec(strCmd.c_str());

  TChain* tree = new TChain("SelectedTree");
  tree->Add(cinput.c_str());
  TTree* reducedTree = new TTree("ReducedTree", "");

  float GenVMass=85.6;
  float GenTripleVMass=125+85.6;
  TString strGenVMass;
  if (isLeptonic==1){
    strGenVMass = "GenDileptonMass";
    tree->SetBranchAddress(strGenVMass, &GenVMass);
  }
  else if (isLeptonic==0){
    strGenVMass = "GenDijetMass";
    tree->SetBranchAddress(strGenVMass, &GenVMass);
  }
  tree->SetBranchAddress(Sqrts->GetName(), &GenTripleVMass);
//  reducedTree->Branch(Sqrts->GetName(), &GenTripleVMass);
  for (int v=0; v<nVars; v++){
    tree->SetBranchAddress(strKDs[v], (kd_vars+v));
    reducedTree->Branch(strKDs[v], (kd_vars+v));
  }
  double avgGenTripleVMass=0;
  int nfilled=0;
  for (int ev=0; ev<tree->GetEntries(); ev++){
    tree->GetEntry(ev);
    if (GenVMass>0){
      reducedTree->Fill();
      nfilled++;
      avgGenTripleVMass+=GenTripleVMass;
    }
  }
  avgGenTripleVMass /= nfilled;
  cout << "Average mVH: " << avgGenTripleVMass << endl;

  RooDataSet* dataSM = new RooDataSet("data", "data", reducedTree, treeargs);
  for (int plotIndex=0; plotIndex<(bdo3D ? 3 : nVars-1); plotIndex++){
    cout << plotIndex << endl;

    RooPlot* plot = measurables[plotIndex]->frame(nbins);
    plot->GetXaxis()->CenterTitle();
    plot->GetYaxis()->SetTitleOffset(1.2);
    plot->GetYaxis()->CenterTitle();
    plot->GetYaxis()->SetTitle("Number of Events");
    plot->GetXaxis()->SetNdivisions(-505);

    string m_name = measurables[plotIndex]->GetName();
    plot->SetTitle(m_name.c_str());

    Sqrts->setConstant(kFALSE);
    dataSM->plotOn(plot, MarkerColor(kRed), MarkerStyle(3), MarkerSize(1.2), LineWidth(0), XErrorSize(0), DataError(RooAbsData::Poisson));
    Sqrts->setVal(avgGenTripleVMass);
    Sqrts->setConstant(kTRUE);
    if (!bdo3D) someHiggs->PDF->plotOn(plot, LineColor(kRed), LineWidth(2));
    else someHiggs->PDF_3D->plotOn(plot, LineColor(kRed), LineWidth(2));

    TGaxis::SetMaxDigits(3);

    TCanvas* can = new TCanvas("can", "can", 600, 600);
    plot->Draw();

    char cname[200];
    sprintf(cname, "%s", m_name.c_str());
    string cname_pdf=cname;
    string cname_eps=cname;
    string cname_png=cname;
    cname_pdf = (coutput + cname_pdf) + ".pdf";
    cname_eps = (coutput + cname_eps) + ".eps";
    cname_png = (coutput + cname_png) + ".png";

    can->SaveAs(cname_pdf.c_str());
    can->SaveAs(cname_eps.c_str());
    can->SaveAs(cname_png.c_str());
    can->Close();
  }

  delete dataSM;
  delete reducedTree;
  delete tree;
  delete someHiggs;
  delete mzz;
  delete Sqrts;
  delete mV;
  delete hs;
  delete h1;
  delete h2;
  delete Phi;
  delete Phi1;
}
