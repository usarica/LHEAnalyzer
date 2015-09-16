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
#include "include/ScalarPdfFactory_ggH.h"

using namespace RooFit;
using namespace std;

void angularDistributions_spin0_ggH(string cinput, string coutdir, double g1Re=1, double g2Re=0, double g4Re=0, double g1L1Re=0, double g2Im=0, double g4Im=0, double g1L1Im=0, int nbins=80){
  RooRealVar* mzz = new RooRealVar("GenHMass", "M_{ZZ} (GeV)", 125, 125.-0.02, 125.+0.02);
  RooRealVar* z1mass = new RooRealVar("GenZ1Mass", "m_{Z1} (GeV)", 0.0, 120);
  RooRealVar* z2mass = new RooRealVar("GenZ2Mass", "m_{Z2} (GeV)", 0.0, 70);
  RooRealVar* hs = new RooRealVar("Gencosthetastar", "cos#theta^{*}", -1, 1);
  RooRealVar* h1 = new RooRealVar("GenhelcosthetaZ1", "cos#theta_{Z1}", -1, 1);
  RooRealVar* h2 = new RooRealVar("GenhelcosthetaZ2", "cos#theta_{Z2}", -1, 1);
  RooRealVar* Phi = new RooRealVar("Genhelphi", "#Phi", -TMath::Pi(), TMath::Pi());
  RooRealVar* Phi1 = new RooRealVar("GenphistarZ1", "#Phi_{Z1}", -TMath::Pi(), TMath::Pi());
  RooRealVar* Y = new RooRealVar("GenY", "Y", 0);

  RooSpinZero::modelMeasurables measurables_;
  measurables_.h1 = h1;
  measurables_.h2 = h2;
  measurables_.Phi = Phi;
  measurables_.m1 = z1mass;
  measurables_.m2 = z2mass;
  measurables_.m12 = mzz;
  measurables_.hs = hs;
  measurables_.Phi1 = Phi1;
  measurables_.Y = Y;

  const int nVars = 8;
  RooArgSet treeargs(*z1mass, *z2mass, *hs, *h1, *h2, *Phi, *Phi1);
  RooRealVar* measurables[nVars-1]={ z1mass, z2mass, h1, h2, hs, Phi, Phi1 };
  float kd_vars[nVars];
  TString strKDs[nVars]={
    "GenHMass", "GenZ1Mass", "GenZ2Mass",
    "GenhelcosthetaZ1", "GenhelcosthetaZ2", "GenphistarZ1",
    "Gencosthetastar", "Genhelphi"
  };
  TString str_genFinalState = "genFinalState";
  TString strGenLepId[2]={ "GenLep1Id", "GenLep3Id" };
  int genFinalState=2;
  int GenLepId[2]={ 0 };

  ScalarPdfFactory_ggH* someHiggs = new ScalarPdfFactory_ggH(measurables_);
  someHiggs->makeParamsConst(false);
  RooRealVar* g1List[8][2];
  RooRealVar* g2List[8][2];
  //RooRealVar* g3List[8][2];
  RooRealVar* g4List[8][2];
  for (int gg=0; gg<8; gg++){
    for (int im=0; im<2; im++){
      g1List[gg][im] = (RooRealVar*)someHiggs->parameters.g1List[gg][im];
      g2List[gg][im] = (RooRealVar*)someHiggs->parameters.g2List[gg][im];
      //g3List[gg][im] = (RooRealVar*)someHiggs->parameters.g3List[gg][im];
      g4List[gg][im] = (RooRealVar*)someHiggs->parameters.g4List[gg][im];
    }
  }
  g1List[0][0]->setVal(g1Re);
  g1List[2][0]->setVal(g1L1Re);
  g2List[0][0]->setVal(g2Re);
  g4List[0][0]->setVal(g4Re);
  g1List[0][1]->setVal(0);
  g1List[2][1]->setVal(g1L1Im);
  g2List[0][1]->setVal(g2Im);
  g4List[0][1]->setVal(g4Im);
  someHiggs->makeParamsConst(true);

  size_t lastSlash = cinput.find_last_of("/\\");
  string finName = cinput.substr(lastSlash+1);
  finName = finName.substr(0, finName.find(".root"));
  finName = finName + "/";
  string coutput = coutdir + finName;
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
    if ((kd_vars[0]-125.)>=0.02) continue;
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
    RooSpinZero_7DComplex_withAccep_ggH* pdf = (RooSpinZero_7DComplex_withAccep_ggH*)someHiggs->getPDF();
    pdf->plotOn(plot, LineColor(kRed), LineWidth(2));

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
  delete Y;
}
