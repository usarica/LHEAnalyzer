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
#include "src/ScalarPdfFactory_ggH.cc"


using namespace RooFit;
using namespace std;

void angularDistributions_spin0_ggH(TString INPUT_NAME, double g1Re=0, double g2Re=0, double g4Re=0, double vfL1=0, double g2Im=0, double g4Im=0, int nbins=80){
  RooRealVar* mzz = new RooRealVar("GenHMass", "M_{ZZ} (GeV)", 125, 125.-0.02, 125.+0.02);
  RooRealVar* z1mass = new RooRealVar("GenZ1Mass", "m_{Z1} (GeV)", 0.0, 120);
  RooRealVar* z2mass = new RooRealVar("GenZ2Mass", "m_{Z2} (GeV)", 0.0, 70);
  RooRealVar* hs = new RooRealVar("Gencosthetastar", "cos#theta^{*}", -1, 1);
  RooRealVar* h1 = new RooRealVar("GenhelcosthetaZ1", "cos#theta_{Z1}", -1, 1);
  RooRealVar* h2 = new RooRealVar("GenhelcosthetaZ2", "cos#theta_{Z2}", -1, 1);
  RooRealVar* Phi = new RooRealVar("Genhelphi", "#Phi", -TMath::Pi(), TMath::Pi());
  RooRealVar* Phi1 = new RooRealVar("GenphistarZ1", "#Phi_{Z1}", -TMath::Pi(), TMath::Pi());

  RooArgSet treeargs(*mzz, *z1mass, *z2mass, *hs, *h1, *h2, *Phi, *Phi1);
  RooRealVar* measurables[8]={ z1mass, z2mass, h1, h2, hs, Phi, Phi1, mzz };
  float kd_vars[8];
  TString strKDs[8]={
    "GenHMass", "GenZ1Mass", "GenZ2Mass",
    "GenhelcosthetaZ1", "GenhelcosthetaZ2", "GenphistarZ1",
    "Gencosthetastar", "Genhelphi"
  };
  TString str_genFinalState = "genFinalState";
  TString strGenLepId[2]={ "GenLep1Id", "GenLep3Id" };
  int genFinalState=2;
  int GenLepId[2]={ 0 };

  ScalarPdfFactory_ggH* someHiggs = new ScalarPdfFactory_ggH(z1mass, z2mass, hs, h1, h2, Phi, Phi1, mzz, 1, false, true);
  someHiggs->_modelParams.fL1->setVal(vfL1);
  someHiggs->_modelParams.g1Val->setVal(g1Re);
  someHiggs->_modelParams.g2Val->setVal(g2Re);
  someHiggs->_modelParams.g3Val->setVal(0);
  someHiggs->_modelParams.g4Val->setVal(g4Re);
  someHiggs->_modelParams.g2ValIm->setVal(g2Im);
  someHiggs->_modelParams.g4ValIm->setVal(g4Im);
  someHiggs->makeParamsConst(true);

  string coutput_common = "/scratch0/hep/usarical/SpinWidthPaper_2015/CMSSW_6_1_1/src/Analysis/Validation/Plots/";
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
  for (int v=0; v<8; v++){
    tree->SetBranchAddress(strKDs[v], (kd_vars+v));
    reducedTree->Branch(strKDs[v], (kd_vars+v));
  }
  for (int v=0; v<2; v++) tree->SetBranchAddress(strGenLepId[v], (GenLepId+v));
  tree->SetBranchAddress(str_genFinalState, &genFinalState);

  for (int ev=0; ev<tree->GetEntries(); ev++){
    tree->GetEntry(ev);
    // Do not select events with lepton interference
    if (genFinalState!=2 && genFinalState!=4 && genFinalState!=5) continue;
    // Only select 4l events
    if (
      (GenLepId[0]==11 || GenLepId[0]==13 || GenLepId[0]==15)
      &&
      (GenLepId[1]==11 || GenLepId[1]==13 || GenLepId[1]==15)
      ) reducedTree->Fill();
  }

  RooDataSet* dataSM = new RooDataSet("data", "data", reducedTree, treeargs);
  for (int plotIndex=0; plotIndex<7; plotIndex++){
    cout << plotIndex << endl;

    RooPlot* plot = measurables[plotIndex]->frame(nbins);
    plot->GetXaxis()->CenterTitle();
    plot->GetYaxis()->SetTitleOffset(1.2);
    plot->GetYaxis()->CenterTitle();
    plot->GetYaxis()->SetTitle("Number of Events");
    plot->GetXaxis()->SetNdivisions(-505);

    string m_name = measurables[plotIndex]->GetName();
    plot->SetTitle(m_name.c_str());

    dataSM->plotOn(plot, MarkerColor(kRed), MarkerStyle(3), MarkerSize(1.2), LineWidth(0), XErrorSize(0), DataError(RooAbsData::Poisson));
    someHiggs->PDF->plotOn(plot, LineColor(kRed), LineWidth(2));

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
  delete z1mass;
  delete z2mass;
  delete hs;
  delete h1;
  delete h2;
  delete Phi;
  delete Phi1;
}
