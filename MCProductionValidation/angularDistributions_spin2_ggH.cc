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
#include "include/TensorPdfFactory_HVV.h"

using namespace RooFit;
using namespace std;

const double my_bList[10][2]={ // For production
  { 1, 0 }, // b1
  { 0, 0 }, // b2
  { 0, 0 }, // b3
  { 0, 0 }, // b4
  { 0, 0 }, // b5
  { 0, 0 }, // b6
  { 0, 0 }, // b7
  { 0, 0 }, // b8
  { 0, 0 }, // b9
  { 0, 0 }  // b10
};

void angularDistributions_spin2_ggH(string cinput, string coutdir, double fz1, double fz2, int nbins=80, bool isHM=false){
  RooRealVar* mzz = new RooRealVar("GenHMass", "M_{ZZ} (GeV)", 125, 125.-0.02, 125.+0.02);
  RooRealVar* z1mass = new RooRealVar("GenZ1Mass", "m_{Z1} (GeV)", 0.0, 120);
  RooRealVar* z2mass = new RooRealVar("GenZ2Mass", "m_{Z2} (GeV)", 0.0, (isHM ? 160. : 70.));
  RooRealVar* hs = new RooRealVar("Gencosthetastar", "cos#theta^{*}", -1, 1);
  RooRealVar* h1 = new RooRealVar("GenhelcosthetaZ1", "cos#theta_{Z1}", -1, 1);
  RooRealVar* h2 = new RooRealVar("GenhelcosthetaZ2", "cos#theta_{Z2}", -1, 1);
  RooRealVar* Phi = new RooRealVar("Genhelphi", "#Phi", -TMath::Pi(), TMath::Pi());
  RooRealVar* Phi1 = new RooRealVar("GenphistarZ1", "#Phi_{Z1}", -TMath::Pi(), TMath::Pi());
  RooRealVar* Y = new RooRealVar("GenY", "Y", 0);

  RooSpinTwo::modelMeasurables measurables_;
  measurables_.h1 = h1;
  measurables_.h2 = h2;
  measurables_.Phi = Phi;
  measurables_.m1 = z1mass;
  measurables_.m2 = z2mass;
  measurables_.m12 = mzz;
  measurables_.hs = hs;
  measurables_.Phi1 = Phi1;
  //measurables_.Y = Y;

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
  int Vdecay=(cinput.find("WW")!=string::npos ? -1 : 1);

  TensorPdfFactory_HVV* someHiggs = new TensorPdfFactory_HVV(measurables_, Vdecay, Vdecay);
  someHiggs->makeParamsConst(false);
  for (int gg=0; gg<10; gg++){
    for (int im=0; im<2; im++){ ((RooRealVar*)someHiggs->parameters.bList[gg][im])->setVal(my_bList[gg][im]); }
  }
  ((RooRealVar*)someHiggs->parameters.f_spinz1)->setVal(fz1);
  ((RooRealVar*)someHiggs->parameters.f_spinz2)->setVal(fz2);
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
      (
      (GenLepId[0]==11 || GenLepId[0]==13 || GenLepId[0]==15)
      &&
      (GenLepId[1]==11 || GenLepId[1]==13 || GenLepId[1]==15)
      && Vdecay==1
      ) || Vdecay==-1
      ) reducedTree->Fill();
  }

  RooDataSet* data = new RooDataSet("data", "data", reducedTree, treeargs);
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

    data->plotOn(plot, MarkerColor(kRed), MarkerStyle(3), MarkerSize(1.2), LineWidth(0), XErrorSize(0), DataError(RooAbsData::Poisson));
    RooSpinTwo_7DComplex_HVV* pdf = (RooSpinTwo_7DComplex_HVV*)someHiggs->getPDF();
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

  delete data;
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
