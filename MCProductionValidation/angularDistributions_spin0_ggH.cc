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

void angularDistributions_spin0_ggH(string cinput, double g1Re=1, double g2Re=0, double g4Re=0, double g1L1Re=0, double g2Im=0, double g4Im=0, double g1L1Im=0, int nbins=80, double mPOLE = 125., int useTaus=1){
  RooRealVar* mzz = new RooRealVar("GenHMass", "M_{ZZ} (GeV)", mPOLE, mPOLE-0.02, mPOLE+0.02);
  RooRealVar* z1mass = new RooRealVar("GenZ1Mass", "m_{Z1} (GeV)", 0.0, min(120., mPOLE));
  RooRealVar* z2mass = new RooRealVar("GenZ2Mass", "m_{Z2} (GeV)", 0.0, min(120., (mPOLE-90.)*20./35.+55.));
  RooRealVar* hs = new RooRealVar("Gencosthetastar", "cos#theta^{*}", -1, 1);
  RooRealVar* h1 = new RooRealVar("GenhelcosthetaZ1", "cos#theta_{Z1}", -1, 1);
  RooRealVar* h2 = new RooRealVar("GenhelcosthetaZ2", "cos#theta_{Z2}", -1, 1);
  RooRealVar* Phi = new RooRealVar("Genhelphi", "#Phi", -TMath::Pi(), TMath::Pi());
  RooRealVar* Phi1 = new RooRealVar("GenphistarZ1", "#Phi_{Z1}", -TMath::Pi(), TMath::Pi());
  RooRealVar* Y = new RooRealVar("GenY", "Y", 0);

  RooSpin::modelMeasurables measurables_;
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
  const int nPlots=nVars-1;
  RooArgSet treeargs;
  RooRealVar* measurables[nPlots]={ hs, Phi1, z1mass, h1, z2mass, h2, Phi };

  float kd_vars[nVars];
  TString strKDs[nVars]={
    "GenHMass", "GenZ1Mass", "GenZ2Mass",
    "GenhelcosthetaZ1", "GenhelcosthetaZ2", "GenphistarZ1",
    "Gencosthetastar", "Genhelphi"
  };
  TString strGenLepId[4]={ "GenLep1Id", "GenLep2Id", "GenLep3Id", "GenLep4Id" };
  int GenLepId[4]={ 0 };
  RooSpin::VdecayType Vdecay1=(cinput.find("WW")!=string::npos ? RooSpin::kVdecayType_Wany : RooSpin::kVdecayType_Zll);
  RooSpin::VdecayType Vdecay2=(cinput.find("WW")!=string::npos ? RooSpin::kVdecayType_Wany : RooSpin::kVdecayType_Zll);
  int nToPlot=nPlots;
  if (cinput.find("ZG")!=string::npos){
    nToPlot-=3;
    Vdecay2=RooSpin::kVdecayType_GammaOnshell;
    if (cinput.find("2l")!=string::npos) Vdecay1=RooSpin::kVdecayType_Zll;
    else if (cinput.find("2nu")!=string::npos) Vdecay1=RooSpin::kVdecayType_Znn;
    else if (cinput.find("2q")!=string::npos) Vdecay1=RooSpin::kVdecayType_Zud;
    else{
      cerr << "Could not find the Z decays! Exiting." << endl;
      return;
    }
  }
  else if (cinput.find("GG")!=string::npos && cinput.find("GGto")==string::npos){ nToPlot-=5; Vdecay1=RooSpin::kVdecayType_GammaOnshell; Vdecay2=RooSpin::kVdecayType_GammaOnshell; }
  else if (cinput.find("ZZ")!=string::npos){
    if (cinput.find("4l")!=string::npos){ Vdecay1=RooSpin::kVdecayType_Zll; Vdecay2=RooSpin::kVdecayType_Zll; }
    else if (cinput.find("4nu")!=string::npos){ Vdecay1=RooSpin::kVdecayType_Znn; Vdecay2=RooSpin::kVdecayType_Znn; }
    else if (cinput.find("4q")!=string::npos){ Vdecay1=RooSpin::kVdecayType_Zud; Vdecay2=RooSpin::kVdecayType_Zud; }
    else if (cinput.find("2l2q")!=string::npos || cinput.find("2q2l")!=string::npos){ Vdecay1=RooSpin::kVdecayType_Zll; Vdecay2=RooSpin::kVdecayType_Zud; }
    else if (cinput.find("2l2nu")!=string::npos || cinput.find("2nu2l")!=string::npos){ Vdecay1=RooSpin::kVdecayType_Zll; Vdecay2=RooSpin::kVdecayType_Znn; }
    else if (cinput.find("2q2nu")!=string::npos || cinput.find("2nu2q")!=string::npos){ Vdecay1=RooSpin::kVdecayType_Zll; Vdecay2=RooSpin::kVdecayType_Zud; }
    else{
      cerr << "Could not find the Z decays! Exiting." << endl;
      return;
    }
  }
  else if (cinput.find("WW")==string::npos){
    cerr << "Could not find the V decays! Exiting." << endl;
    return;
  }
  for (int r=0; r<nToPlot; r++) treeargs.add(*(measurables[r]));
  if (Vdecay1==RooSpin::kVdecayType_GammaOnshell){ z1mass->setVal(0); z1mass->setConstant(true); }
  if (Vdecay2==RooSpin::kVdecayType_GammaOnshell){ z2mass->setVal(0); z2mass->setConstant(true); }

  cout << "Decay modes found are " << Vdecay1 << '\t' << Vdecay2 << endl;
  ScalarPdfFactory_ggH* someHiggs = new ScalarPdfFactory_ggH(measurables_, false, Vdecay1, Vdecay2);
  someHiggs->makeParamsConst(false);

  if ((Vdecay1==RooSpin::kVdecayType_GammaOnshell || Vdecay2==RooSpin::kVdecayType_GammaOnshell) && !(Vdecay1==RooSpin::kVdecayType_GammaOnshell && Vdecay2==RooSpin::kVdecayType_GammaOnshell)){
    ((RooRealVar*)someHiggs->couplings.gzgs1List[0][0])->setVal(g1L1Re);
    ((RooRealVar*)someHiggs->couplings.gzgs2List[0][0])->setVal(g2Re);
    ((RooRealVar*)someHiggs->couplings.gzgs4List[0][0])->setVal(g4Re);
    ((RooRealVar*)someHiggs->couplings.gzgs1List[0][1])->setVal(g1L1Im);
    ((RooRealVar*)someHiggs->couplings.gzgs2List[0][1])->setVal(g2Im);
    ((RooRealVar*)someHiggs->couplings.gzgs4List[0][1])->setVal(g4Im);
  }
  else if (!(Vdecay1==RooSpin::kVdecayType_GammaOnshell && Vdecay2==RooSpin::kVdecayType_GammaOnshell)){
    ((RooRealVar*)someHiggs->couplings.g1List[0][0])->setVal(g1Re);
    ((RooRealVar*)someHiggs->couplings.g1List[2][0])->setVal(g1L1Re);
    ((RooRealVar*)someHiggs->couplings.g2List[0][0])->setVal(g2Re);
    ((RooRealVar*)someHiggs->couplings.g4List[0][0])->setVal(g4Re);
    ((RooRealVar*)someHiggs->couplings.g1List[0][1])->setVal(0);
    ((RooRealVar*)someHiggs->couplings.g1List[2][1])->setVal(g1L1Im);
    ((RooRealVar*)someHiggs->couplings.g2List[0][1])->setVal(g2Im);
    ((RooRealVar*)someHiggs->couplings.g4List[0][1])->setVal(g4Im);
  }
  else{
    ((RooRealVar*)someHiggs->couplings.ggsgs2List[0][0])->setVal(g2Re);
    ((RooRealVar*)someHiggs->couplings.ggsgs4List[0][0])->setVal(g4Re);
    ((RooRealVar*)someHiggs->couplings.ggsgs2List[0][1])->setVal(g2Im);
    ((RooRealVar*)someHiggs->couplings.ggsgs4List[0][1])->setVal(g4Im);
  }
  someHiggs->makeParamsConst(true);

  size_t lastSlash = cinput.find_last_of("/\\");
  string finName;
  finName=cinput;
  finName.resize(lastSlash+1);
  string coutput = finName + "Validation/";
  cout << "Output folder is " << coutput << endl;
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
  for (int v=0; v<4; v++) tree->SetBranchAddress(strGenLepId[v], (GenLepId+v));

  int nRecorded=0;
  for (int ev=0; ev<tree->GetEntries(); ev++){
    tree->GetEntry(ev);
    if ((kd_vars[0]-mPOLE)>=0.02) continue;
    if (
      ((GenLepId[0]==11 || GenLepId[0]==13 || GenLepId[0]==15) && (GenLepId[2]==11 || GenLepId[2]==13 || GenLepId[2]==15) && (Vdecay1==RooSpin::kVdecayType_Zll && Vdecay2==RooSpin::kVdecayType_Zll))
      ||
      ((GenLepId[0]==12 || GenLepId[0]==14 || GenLepId[0]==16) && (GenLepId[2]==12 || GenLepId[2]==14 || GenLepId[2]==16) && (Vdecay1==RooSpin::kVdecayType_Znn && Vdecay2==RooSpin::kVdecayType_Znn))
      ||
      ((GenLepId[0]>0 && GenLepId[0]<6) && (GenLepId[2]>0 && GenLepId[2]<6) && (Vdecay1==RooSpin::kVdecayType_Zud && Vdecay2==RooSpin::kVdecayType_Zud))
      ||
      ((((GenLepId[0]==11 || GenLepId[0]==13 || GenLepId[0]==15) && (GenLepId[2]==12 || GenLepId[2]==14 || GenLepId[2]==16)) || ((GenLepId[2]==11 || GenLepId[2]==13 || GenLepId[2]==15) && (GenLepId[0]==12 || GenLepId[0]==14 || GenLepId[0]==16))) && (Vdecay1==RooSpin::kVdecayType_Zll && Vdecay2==RooSpin::kVdecayType_Znn))
      ||
      ((((GenLepId[0]==11 || GenLepId[0]==13 || GenLepId[0]==15) && (GenLepId[2]>0 && GenLepId[2]<6)) || ((GenLepId[2]==11 || GenLepId[2]==13 || GenLepId[2]==15) && (GenLepId[0]>0 && GenLepId[0]<6))) && (Vdecay1==RooSpin::kVdecayType_Zll && Vdecay2==RooSpin::kVdecayType_Zud))
      ||
      ((((GenLepId[0]==12 || GenLepId[0]==14 || GenLepId[0]==16) && (GenLepId[2]>0 && GenLepId[2]<6)) || ((GenLepId[2]==12 || GenLepId[2]==14 || GenLepId[2]==16) && (GenLepId[0]>0 && GenLepId[0]<6))) && (Vdecay1==RooSpin::kVdecayType_Znn && Vdecay2==RooSpin::kVdecayType_Zud))
      ||
      (Vdecay1==RooSpin::kVdecayType_Wany && Vdecay2==RooSpin::kVdecayType_Wany)
      ||
      (Vdecay1==RooSpin::kVdecayType_GammaOnshell || Vdecay2==RooSpin::kVdecayType_GammaOnshell)
      ){
      if (
        ((GenLepId[2]==11 || GenLepId[2]==13 || GenLepId[2]==15) && (GenLepId[0]==12 || GenLepId[0]==14 || GenLepId[0]==16) && (Vdecay1==RooSpin::kVdecayType_Zll && Vdecay2==RooSpin::kVdecayType_Znn))
        ||
        ((GenLepId[2]==11 || GenLepId[2]==13 || GenLepId[2]==15) && (GenLepId[0]>0 && GenLepId[0]<6) && (Vdecay1==RooSpin::kVdecayType_Zll && Vdecay2==RooSpin::kVdecayType_Zud))
        ||
        ((GenLepId[2]==12 || GenLepId[2]==14 || GenLepId[2]==16) && (GenLepId[0]>0 && GenLepId[0]<6) && (Vdecay1==RooSpin::kVdecayType_Znn && Vdecay2==RooSpin::kVdecayType_Zud))
        ){
        //swap(strKDs[2], strKDs[2]);
        //swap(strKDs[3], strKDs[4]);
        continue; // Don't bother to recalculate angles etc.
      }
      if (GenLepId[0]==GenLepId[2] && GenLepId[1]==GenLepId[3] && Vdecay1!=RooSpin::kVdecayType_GammaOnshell && Vdecay2!=RooSpin::kVdecayType_GammaOnshell) continue;
      if (useTaus==0 && (GenLepId[0]==15 || GenLepId[1]==-15 || GenLepId[2]==15 || GenLepId[3]==-15)) continue;
      reducedTree->Fill();
      nRecorded++;
    }
  }
  cout << "Number of entries recorded: " << nRecorded << '/' << tree->GetEntries() << endl;

  RooDataSet* dataSM = new RooDataSet("data", "data", reducedTree, treeargs);
  for (int plotIndex=0; plotIndex<(int)treeargs.getSize(); plotIndex++){
    cout << plotIndex << endl;

    string m_name = measurables[plotIndex]->GetName();
    int nbins_plot=nbins;
    if (m_name.find("Mass")!=string::npos) nbins_plot*=2;
    RooPlot* plot = measurables[plotIndex]->frame(nbins_plot);
    plot->GetXaxis()->CenterTitle();
    plot->GetYaxis()->SetTitleOffset(1.2);
    plot->GetYaxis()->CenterTitle();
    plot->GetYaxis()->SetTitle("Number of Events");
    plot->GetXaxis()->SetNdivisions(-505);
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

/*
// Example commands
angularDistributions_spin0_ggH("/scratch0/hep/usarical/ggH-JHUGEN125/LHC_13TeV_HZZ4l/0+m/GGtoH125toZZ4l_13TeV.root");
angularDistributions_spin0_ggH("/scratch0/hep/usarical/ggH-JHUGEN125/LHC_13TeV_HZZ4l/0+h/GGtoH125toZZ4l_13TeV.root",0,1);
angularDistributions_spin0_ggH("/scratch0/hep/usarical/ggH-JHUGEN125/LHC_13TeV_HZZ4l/0+L1/GGtoH125toZZ4l_13TeV.root",0,0,0,10000);
angularDistributions_spin0_ggH("/scratch0/hep/usarical/ggH-JHUGEN125/LHC_13TeV_HZZ4l/0-/GGtoH125toZZ4l_13TeV.root"0,0,1,0);
*/

