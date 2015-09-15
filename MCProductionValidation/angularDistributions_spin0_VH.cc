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
#include "TFile.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TH1.h"
#include "TGaxis.h"
#include "TString.h"
#include "TChain.h"
#include "include/ScalarPdfFactory_VH.h"

using namespace RooFit;
using namespace std;

void angularDistributions_spin0_VH(string cinput, string coutdir, double sqrts = 13, double g1Re=1, double g2Re=0, double g4Re=0, double g1L1Re=0, double g2Im=0, double g4Im=0, double g1L1Im=0, int isLeptonic=0, int nbins=80, int PDFType=1){
  sqrts *= 1e3;
  bool bdo3D = (PDFType==1 ? true : false);
  int VHmode=-1;
  if (cinput.find("ZH")!=string::npos) VHmode=1;
  const int nVars=8;
  float kd_vars[nVars];
  TString strKDs[nVars]={
    "GenHMass",
    "GenhelcosthetaV1", "GenhelcosthetaV2", "GenphistarV1",
    "Gencosthetastar", "Genhelphi",
    "", "GenAssociatedVMass"
  };
  for (int v=1; v<nVars-2; v++){
    if (isLeptonic==0) strKDs[v].Append("_VHhadronic");
    else if (isLeptonic==1) strKDs[v].Append("_VHleptonic");
  }
  if (isLeptonic==0) strKDs[6] = "GenDijetVVMass";
  else if (isLeptonic==1) strKDs[6] = "GenDileptonVVMass";

  RooRealVar* m12 = new RooRealVar(strKDs[0], "m_{H} (GeV)", 125, 125.-0.02, 125.+0.02);
  RooRealVar* m1 = new RooRealVar(strKDs[6], (VHmode==1 ? "m_{ZH} (GeV)" : "m_{WH} (GeV)"), 460, (VHmode==1 ? 125.+79. : 125.+70.), 1000);
  RooRealVar* m2 = new RooRealVar(strKDs[7], (VHmode==1 ? "m_{Z} (GeV)" : "m_{W} (GeV)"), (VHmode==1 ? 91.1876 : 80.399), 4, 120);
  RooRealVar* h1 = new RooRealVar(strKDs[1], "cos#theta_{V*}", -1, 1);
  RooRealVar* h2 = new RooRealVar(strKDs[2], "cos#theta_{V}", -1, 1);
  RooRealVar* Phi1 = new RooRealVar(strKDs[3], "#Phi_{V*}", -TMath::Pi(), TMath::Pi());
  RooRealVar* hs = new RooRealVar(strKDs[4], "cos#theta^{*}", -1, 1);
  RooRealVar* Phi = new RooRealVar(strKDs[5], "#Phi", -TMath::Pi(), TMath::Pi());
  RooRealVar* Y = new RooRealVar("GenY", "Y", 0, 0, sqrts);

  RooSpinZero::modelMeasurables measurables_;
  measurables_.h1 = h1;
  measurables_.h2 = h2;
  measurables_.Phi = Phi;
  measurables_.m1 = m1;
  measurables_.m2 = m2;
  measurables_.m12 = m12;
  measurables_.hs = hs;
  measurables_.Phi1 = Phi1;
  measurables_.Y = Y;

  RooArgSet treeargs(*h1, *h2, *Phi, *hs, *Phi1, *m12, *m1, *m2);
  RooRealVar* measurables[nVars]={ h1, h2, Phi, hs, Phi1, m12, m1, m2 };

  cout << "Initiating vh" << endl;

  ScalarPdfFactory_VH* someHiggs = new ScalarPdfFactory_VH(measurables_, sqrts, PDFType, VHmode);
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

  cout << "params set" << endl;

  size_t lastSlash = cinput.find_last_of("/\\");
  string finName = cinput.substr(lastSlash+1);
  finName = finName.substr(0,finName.find(".root"));
  finName = finName + "/";
  string coutput = coutdir + finName;
  string strCmd = "mkdir -p ";
  strCmd.append(coutput);
  gSystem->Exec(strCmd.c_str());

  TChain* tree = new TChain("SelectedTree");
  tree->Add(cinput.c_str());
  TFile* foutput = new TFile((coutput + "/plots.root").c_str(), "recreate");
  TTree* reducedTree = new TTree("ReducedTree", "");

  cout << tree->GetEntries() << endl;

  for (int v=0; v<nVars-1; v++){
    tree->SetBranchAddress(strKDs[v], (kd_vars+v));
    reducedTree->Branch(strKDs[v], (kd_vars+v));
  }

  vector<double> * GenAssociatedVMass = 0;
  vector<double> * GenMotherPz = 0;
  double GenVMass = 0;
  double GenY = 0;
  tree->SetBranchAddress(strKDs[nVars-1], &GenAssociatedVMass);
  tree->SetBranchAddress("GenMotherPz", &GenMotherPz);
  reducedTree->Branch(strKDs[nVars-1], &GenVMass);
  reducedTree->Branch(Y->GetName(), &GenY);

  int nfilled=0;
  for (int ev=0; ev<tree->GetEntries(); ev++){
    tree->GetEntry(ev);
    GenVMass=0;
    GenY=0;
    if (GenAssociatedVMass->size()>0) GenVMass = GenAssociatedVMass->at(0);
    if (GenVMass>0){
      if (GenMotherPz->size()>1 && kd_vars[6]>0) GenY = (GenMotherPz->at(0)+GenMotherPz->at(1))/kd_vars[6];
      else continue;
      reducedTree->Fill();
      nfilled++;
    }
  }
  foutput->WriteTObject(reducedTree);

  cout << "datasm" << endl;


  RooDataSet* dataSM = new RooDataSet("data", "data", reducedTree, treeargs);
  cout << "datasmdone" << endl;
  for (int plotIndex=0; plotIndex<nVars; plotIndex++){
    if (plotIndex==nVars-3) continue;
    if (plotIndex==nVars-1) continue;

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
    someHiggs->getPDF()->plotOn(plot, LineColor(kRed), LineWidth(2));

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

    foutput->WriteTObject(can);
    can->SaveAs(cname_pdf.c_str());
    can->SaveAs(cname_eps.c_str());
    can->SaveAs(cname_png.c_str());
    can->Close();
  }

  delete dataSM;
  delete reducedTree;
  foutput->Close();
  delete tree;
  delete someHiggs;
  delete m12;
  delete m1;
  delete m2;
  delete hs;
  delete h1;
  delete h2;
  delete Phi;
  delete Phi1;
  delete Y;
}
