#ifndef SCALAR_PDF_FACTORY_GGH
#define SCALAR_PDF_FACTORY_GGH

#include "../include/RooSpinZero_7DComplex_withAccep_ggH.h"
#include "RooPolynomial.h"
#include "RooProdPdf.h"
#include "TF1.h"
#include "TMath.h"

class ScalarPdfFactory_ggH{

public:

  RooSpinZero_7DComplex_withAccep_ggH::measurables _measurables;
  RooSpinZero_7DComplex_withAccep_ggH::modelParameters _modelParams;
  RooSpinZero_7DComplex_withAccep_ggH::accepParameters _accepParams;

  /*
  RooRealVar* aM1;
  RooRealVar* bM1;
  RooRealVar* cM1;
  RooRealVar* dM1;

  RooRealVar* aM2;
  RooRealVar* bM2;
  RooRealVar* cM2;
  RooRealVar* dM2;

  RooPolynomial* m1Accep;
  RooPolynomial* m2Accep;
  */

  //RooSpinZero_7DComplex_withAccep_ggH* PDFangleAccep;
  RooSpinZero_7DComplex_withAccep_ggH* PDF;

  bool acceptance_;

  //RooProdPdf* PDF;

  ScalarPdfFactory_ggH(){};

  ScalarPdfFactory_ggH(RooRealVar* m1, RooRealVar* m2, RooRealVar* hs, RooRealVar* h1, RooRealVar* h2, RooRealVar* Phi, RooRealVar* Phi1, RooRealVar* mZZ, int para, bool acceptance, bool pmf_applied=false){

    _measurables.m1   = m1;
    _measurables.m2   = m2;
    _measurables.h1   = h1;
    _measurables.h2   = h2;
    _measurables.hs   = hs;
    _measurables.Phi  = Phi;
    _measurables.Phi1 = Phi1;
    _modelParams.mX   = mZZ;

    // model parameters
    _modelParams.mZ = new RooRealVar("mZ", "mZ", 91.1876);
    _modelParams.gamZ= new RooRealVar("gamZ", "gamZ", 2.49);

    _modelParams.a1Val = new RooRealVar("a1Val", "a1Val", 0., -1e15, 1e15);
    _modelParams.phi1Val = new RooRealVar("phi1Val", "phi1Val", 0., -2.*TMath::Pi(), 2*TMath::Pi());
    _modelParams.a2Val = new RooRealVar("a2Val", "a2Val", 0., -1e15, 1e15);
    _modelParams.phi2Val = new RooRealVar("phi2Val", "phi2Val", 0., -2.*TMath::Pi(), 2*TMath::Pi());
    _modelParams.a3Val = new RooRealVar("a3Val", "a3Val", 0., -1e15, 1e15);
    _modelParams.phi3Val = new RooRealVar("phi3Val", "phi3Val", 0., -2.*TMath::Pi(), 2*TMath::Pi());

    _modelParams.parameterization = para;

    _modelParams.g1Val = new RooRealVar("g1Val", "g1Val", 1., -1e15, 1e15);
    _modelParams.g2Val = new RooRealVar("g2Val", "g2Val", 0., -1e15, 1e15);
    _modelParams.g3Val = new RooRealVar("g3Val", "g3Val", 0., -1e15, 1e15);
    _modelParams.g4Val = new RooRealVar("g4Val", "g4Val", 0., -1e15, 1e15);

    _modelParams.g1ValIm = new RooRealVar("g1ValIm", "g1ValIm", 0., -1e15, 1e15);
    _modelParams.g2ValIm = new RooRealVar("g2ValIm", "g2ValIm", 0., -1e15, 1e15);
    _modelParams.g3ValIm = new RooRealVar("g3ValIm", "g3ValIm", 0., -1e15, 1e15);
    _modelParams.g4ValIm = new RooRealVar("g4ValIm", "g4ValIm", 0., -1e15, 1e15);

    if (!pmf_applied){
      _modelParams.fL1 = new RooRealVar("fL1", "f_{#Lambda1}", 0., 0., 1.0);
      _modelParams.fa2 = new RooRealVar("fa2", "f_{g2}", 0., 0., 1.0);
      _modelParams.fa3 = new RooRealVar("fa3", "f_{g4}", 0., 0., 1.0);
      _modelParams.phia2 = new RooRealVar("phia2", "#phi_{g2}", 0., -2.*TMath::Pi(), 2*TMath::Pi());
      _modelParams.phia3 = new RooRealVar("phia3", "#phi_{g4}", 0., -2.*TMath::Pi(), 2*TMath::Pi());
    }
    else{
      _modelParams.fL1 = new RooRealVar("fL1", "f_{#Lambda1}", 0., -1.0, 1.0);
      _modelParams.fa2 = new RooRealVar("fa2", "f_{g2}", 0., -1.0, 1.0);
      _modelParams.fa3 = new RooRealVar("fa3", "f_{g4}", 0., -1.0, 1.0);
      _modelParams.phia2 = new RooRealVar("phia2", "#phi_{g2}", 0.);
      _modelParams.phia3 = new RooRealVar("phia3", "#phi_{g4}", 0.);
    }
    _modelParams.R1Val = new RooRealVar("R1Val", "R1Val", 0.15);
    _modelParams.R2Val = new RooRealVar("R2Val", "R2Val", 0.15);

    // acceptance parameters
    acceptance_ = acceptance;

    if (acceptance_){
      _accepParams.aPhi = new RooRealVar("aPhi", "aPhi", 1.);
      _accepParams.bPhi = new RooRealVar("bPhi", "bPhi", 4.88199e-03);
      _accepParams.cPhi = new RooRealVar("cPhi", "cPhi", 3.69579e-02);
      _accepParams.dPhi = new RooRealVar("dPhi", "dPhi", 0.);
      _accepParams.ePhi = new RooRealVar("ePhi", "ePhi", 0.);

      _accepParams.aPhi1 = new RooRealVar("aPhi1", "aPhi1", 1.);
      _accepParams.bPhi1 = new RooRealVar("bPhi1", "bPhi1", -1.27958e-02);
      _accepParams.cPhi1 = new RooRealVar("cPhi1", "cPhi1", -1.64892e-01);
      _accepParams.dPhi1 = new RooRealVar("dPhi1", "dPhi1", 0.);
      _accepParams.ePhi1 = new RooRealVar("ePhi1", "ePhi1", 0.);

      _accepParams.aH1 = new RooRealVar("aH1", "aH1", 1.);
      _accepParams.bH1 = new RooRealVar("bH1", "bH1", 2.64540e-02);
      _accepParams.cH1 = new RooRealVar("cH1", "cH1", 0.);
      _accepParams.dH1 = new RooRealVar("dH1", "dH1", 0.);
      _accepParams.eH1 = new RooRealVar("eH1", "eH1", 0.);

      _accepParams.aH2 = new RooRealVar("aH2", "aH2", 1.);
      _accepParams.bH2 = new RooRealVar("bH2", "bH2", -3.73167e-01);
      _accepParams.cH2 = new RooRealVar("cH2", "cH2", 0.);
      _accepParams.dH2 = new RooRealVar("dH2", "dH2", 0.);
      _accepParams.eH2 = new RooRealVar("eH2", "eH2", 0.);

      _accepParams.aHs = new RooRealVar("aHs", "aHs", 1.);
      _accepParams.bHs = new RooRealVar("bHs", "bHs", -1.55528e-01);
      _accepParams.cHs = new RooRealVar("cHs", "cHs", 0.);
      _accepParams.dHs = new RooRealVar("dHs", "dHs", 0.);
      _accepParams.eHs = new RooRealVar("eHs", "eHs", 0.);

      _accepParams.aM1 = new RooRealVar("aM1", "aM1", 1.);
      _accepParams.bM1 = new RooRealVar("bM1", "bM1", -1.26554e-02);
      _accepParams.cM1 = new RooRealVar("cM1", "cM1", 3.13526e-05);
      _accepParams.dM1 = new RooRealVar("dM1", "dM1", 0.);

      _accepParams.aM2 = new RooRealVar("aM2", "aM2", 1.);
      _accepParams.bM2 = new RooRealVar("bM2", "bM2", 5.75519e-04);
      _accepParams.cM2 = new RooRealVar("cM2", "cM2", -7.74696e-05);
      _accepParams.dM2 = new RooRealVar("dM2", "dM2", 0.);

    }

    else {

      _accepParams.aPhi = new RooRealVar("aPhi", "aPhi", 1.);
      _accepParams.bPhi = new RooRealVar("bPhi", "bPhi", 0.);
      _accepParams.cPhi = new RooRealVar("cPhi", "cPhi", 0.);
      _accepParams.dPhi = new RooRealVar("dPhi", "dPhi", 0.);
      _accepParams.ePhi = new RooRealVar("ePhi", "ePhi", 0.);

      _accepParams.aPhi1 = new RooRealVar("aPhi1", "aPhi1", 1.);
      _accepParams.bPhi1 = new RooRealVar("bPhi1", "bPhi1", 0.);
      _accepParams.cPhi1 = new RooRealVar("cPhi1", "cPhi1", 0.);
      _accepParams.dPhi1 = new RooRealVar("dPhi1", "dPhi1", 0.);
      _accepParams.ePhi1 = new RooRealVar("ePhi1", "ePhi1", 0.);

      _accepParams.aH1 = new RooRealVar("aH1", "aH1", 1.);
      _accepParams.bH1 = new RooRealVar("bH1", "bH1", 0.);
      _accepParams.cH1 = new RooRealVar("cH1", "cH1", 0.);
      _accepParams.dH1 = new RooRealVar("dH1", "dH1", 0.);
      _accepParams.eH1 = new RooRealVar("eH1", "eH1", 0.);

      _accepParams.aH2 = new RooRealVar("aH2", "aH2", 1.);
      _accepParams.bH2 = new RooRealVar("bH2", "bH2", 0.);
      _accepParams.cH2 = new RooRealVar("cH2", "cH2", 0.);
      _accepParams.dH2 = new RooRealVar("dH2", "dH2", 0.);
      _accepParams.eH2 = new RooRealVar("eH2", "eH2", 0.);

      _accepParams.aHs = new RooRealVar("aHs", "aHs", 1.);
      _accepParams.bHs = new RooRealVar("bHs", "bHs", 0.);
      _accepParams.cHs = new RooRealVar("cHs", "cHs", 0.);
      _accepParams.dHs = new RooRealVar("dHs", "dHs", 0.);
      _accepParams.eHs = new RooRealVar("eHs", "eHs", 0.);

      _accepParams.aM1 = new RooRealVar("aM1", "aM1", 0.);
      _accepParams.bM1 = new RooRealVar("bM1", "bM1", 0.);
      _accepParams.cM1 = new RooRealVar("cM1", "cM1", 0.);
      _accepParams.dM1 = new RooRealVar("dM1", "dM1", 0.);

      _accepParams.aM2 = new RooRealVar("aM2", "aM2", 0.);
      _accepParams.bM2 = new RooRealVar("bM2", "bM2", 0.);
      _accepParams.cM2 = new RooRealVar("cM2", "cM2", 0.);
      _accepParams.dM2 = new RooRealVar("dM2", "dM2", 0.);

      _accepParams.aPhi->setConstant(kTRUE);
      _accepParams.bPhi->setConstant(kTRUE);
      _accepParams.cPhi->setConstant(kTRUE);
      _accepParams.dPhi->setConstant(kTRUE);
      _accepParams.ePhi->setConstant(kTRUE);

      _accepParams.aPhi1->setConstant(kTRUE);
      _accepParams.bPhi1->setConstant(kTRUE);
      _accepParams.cPhi1->setConstant(kTRUE);
      _accepParams.dPhi1->setConstant(kTRUE);
      _accepParams.ePhi1->setConstant(kTRUE);

      _accepParams.aH1->setConstant(kTRUE);
      _accepParams.bH1->setConstant(kTRUE);
      _accepParams.cH1->setConstant(kTRUE);
      _accepParams.dH1->setConstant(kTRUE);
      _accepParams.eH1->setConstant(kTRUE);

      _accepParams.aH2->setConstant(kTRUE);
      _accepParams.bH2->setConstant(kTRUE);
      _accepParams.cH2->setConstant(kTRUE);
      _accepParams.dH2->setConstant(kTRUE);
      _accepParams.eH2->setConstant(kTRUE);

      _accepParams.aHs->setConstant(kTRUE);
      _accepParams.bHs->setConstant(kTRUE);
      _accepParams.cHs->setConstant(kTRUE);
      _accepParams.dHs->setConstant(kTRUE);
      _accepParams.eHs->setConstant(kTRUE);

      _accepParams.aM1->setConstant(kTRUE);
      _accepParams.bM1->setConstant(kTRUE);
      _accepParams.cM1->setConstant(kTRUE);
      _accepParams.dM1->setConstant(kTRUE);

      _accepParams.aM2->setConstant(kTRUE);
      _accepParams.bM2->setConstant(kTRUE);
      _accepParams.cM2->setConstant(kTRUE);
      _accepParams.dM2->setConstant(kTRUE);
    }
    //m1Accep = new RooPolynomial("m1Accep","m1Accep",*_measurables.m1,RooArgList(*aM1,*bM1,*cM1,*dM1),4);
    //m2Accep = new RooPolynomial("m2Accep","m2Accep",*_measurables.m2,RooArgList(*aM2,*bM2,*cM2,*dM2),4);

    PDF = new RooSpinZero_7DComplex_withAccep_ggH("PDF", "PDF",
      _measurables,
      _modelParams,
      _accepParams);

    //PDF = new RooProdPdf("PDF","PDF",RooArgList(*m1Accep,*m2Accep,*PDFangleAccep));


    //std::cout << "ScalarPdfFactory_ggH::ScalarPdfFactory_ggH - done " << std::endl;

  }

  ~ScalarPdfFactory_ggH(){

    //std::cout << "~ScalarPdfFactory_ggH" << std::endl;

    delete PDF;
    //delete PDFangleAccep;
    //delete m1Accep;
    //delete m2Accep;

    delete _accepParams.aM1;
    delete _accepParams.bM1;
    delete _accepParams.cM1;
    delete _accepParams.dM1;

    delete _accepParams.aM2;
    delete _accepParams.bM2;
    delete _accepParams.cM2;
    delete _accepParams.dM2;

    delete _modelParams.a1Val;
    delete _modelParams.phi1Val;
    delete _modelParams.a2Val;
    delete _modelParams.phi2Val;
    delete _modelParams.a3Val;
    delete _modelParams.phi3Val;
    delete _modelParams.g1Val;
    delete _modelParams.g2Val;
    delete _modelParams.g3Val;
    delete _modelParams.g4Val;
    delete _modelParams.g1ValIm;
    delete _modelParams.g2ValIm;
    delete _modelParams.g3ValIm;
    delete _modelParams.g4ValIm;
    delete _modelParams.fL1;
    delete _modelParams.fa2;
    delete _modelParams.fa3;
    delete _modelParams.phia2;
    delete _modelParams.phia3;

    delete _modelParams.R1Val;
    delete _modelParams.R2Val;

    delete _modelParams.mZ;
    delete _modelParams.gamZ;

    delete _accepParams.aPhi;
    delete _accepParams.bPhi;
    delete _accepParams.cPhi;
    delete _accepParams.dPhi;
    delete _accepParams.ePhi;
    delete _accepParams.aPhi1;
    delete _accepParams.bPhi1;
    delete _accepParams.cPhi1;
    delete _accepParams.dPhi1;
    delete _accepParams.ePhi1;
    delete _accepParams.aH1;
    delete _accepParams.bH1;
    delete _accepParams.cH1;
    delete _accepParams.dH1;
    delete _accepParams.eH1;
    delete _accepParams.aH2;
    delete _accepParams.bH2;
    delete _accepParams.cH2;
    delete _accepParams.dH2;
    delete _accepParams.eH2;
    delete _accepParams.aHs;
    delete _accepParams.bHs;
    delete _accepParams.cHs;
    delete _accepParams.dHs;
    delete _accepParams.eHs;

    //std::cout << "~ScalarPdfFactory_ggH - end " << std::endl;

  }

  void makeSMHiggs(){

    //std::cout << "ScalarPdfFactory_ggH::makeSMHiggs()" << std::endl;

    _modelParams.parameterization=1;

    _modelParams.fL1->setVal(0.0);

    _modelParams.g1Val->setVal(1.0);
    _modelParams.g2Val->setVal(0.0);
    _modelParams.g3Val->setVal(0.0);
    _modelParams.g4Val->setVal(0.0);

    _modelParams.g1ValIm->setVal(0.0);
    _modelParams.g2ValIm->setVal(0.0);
    _modelParams.g3ValIm->setVal(0.0);
    _modelParams.g4ValIm->setVal(0.0);

  }

  void makeLGHiggs(){

    //std::cout << "ScalarPdfFactory_ggH::makeLGHiggs()" << std::endl;

    _modelParams.parameterization=1;

    _modelParams.fL1->setVal(0.0);

    _modelParams.g1Val->setVal(0.0);        // need to calculate the proper normalizations
    _modelParams.g2Val->setVal(1.0);
    _modelParams.g3Val->setVal(0.0);
    _modelParams.g4Val->setVal(0.0);

    _modelParams.g1ValIm->setVal(0.0);        // need to calculate the proper normalizations
    _modelParams.g2ValIm->setVal(0.0);
    _modelParams.g3ValIm->setVal(0.0);
    _modelParams.g4ValIm->setVal(0.0);

  }


  void makePSHiggs(){

    //std::cout << "ScalarPdfFactory_ggH::makePSHiggs()" << std::endl;

    _modelParams.parameterization=1;

    _modelParams.fL1->setVal(0.0);

    _modelParams.g1Val->setVal(0.0);
    _modelParams.g2Val->setVal(0.0);
    _modelParams.g3Val->setVal(0.0);
    _modelParams.g4Val->setVal(1.0);

    _modelParams.g1ValIm->setVal(0.0);
    _modelParams.g2ValIm->setVal(0.0);
    _modelParams.g3ValIm->setVal(0.0);
    _modelParams.g4ValIm->setVal(0.0);


  }

  void makeCustom(double a1, double a2, double a3,
    double phi1, double phi2, double phi3){

    //std::cout << "ScalarPdfFactory_ggH::makeCustom()" << std::endl;

    _modelParams.parameterization=0;

    _modelParams.fL1->setVal(0.0);

    _modelParams.a1Val->setVal(a1);
    _modelParams.phi1Val->setVal(phi1);
    _modelParams.a2Val->setVal(a2);
    _modelParams.phi2Val->setVal(phi2);
    _modelParams.a3Val->setVal(a3);
    _modelParams.phi3Val->setVal(phi3);

  }

  void makeParamsConst(bool yesNo=true){

    //std::cout << "ScalarPdfFactory_ggH::makeParamsConst()" << std::endl;

    if (yesNo){
      _modelParams.a1Val->setConstant(kTRUE);
      _modelParams.phi1Val->setConstant(kTRUE);
      _modelParams.a2Val->setConstant(kTRUE);
      _modelParams.phi2Val->setConstant(kTRUE);
      _modelParams.a3Val->setConstant(kTRUE);
      _modelParams.phi3Val->setConstant(kTRUE);
      _modelParams.g1Val->setConstant(kTRUE);
      _modelParams.g2Val->setConstant(kTRUE);
      _modelParams.g3Val->setConstant(kTRUE);
      _modelParams.g4Val->setConstant(kTRUE);
      _modelParams.g1ValIm->setConstant(kTRUE);
      _modelParams.g2ValIm->setConstant(kTRUE);
      _modelParams.g3ValIm->setConstant(kTRUE);
      _modelParams.g4ValIm->setConstant(kTRUE);
      _modelParams.fL1->setConstant(kTRUE);
      _modelParams.fa2->setConstant(kTRUE);
      _modelParams.fa3->setConstant(kTRUE);
      _modelParams.phia2->setConstant(kTRUE);
      _modelParams.phia3->setConstant(kTRUE);
      _modelParams.gamZ->setConstant(kTRUE);
      _modelParams.mZ->setConstant(kTRUE);
      _modelParams.R1Val->setConstant(kTRUE);
      _modelParams.R2Val->setConstant(kTRUE);

      _accepParams.aPhi->setConstant(kTRUE);
      _accepParams.bPhi->setConstant(kTRUE);
      _accepParams.cPhi->setConstant(kTRUE);
      _accepParams.dPhi->setConstant(kTRUE);
      _accepParams.ePhi->setConstant(kTRUE);

      _accepParams.aPhi1->setConstant(kTRUE);
      _accepParams.bPhi1->setConstant(kTRUE);
      _accepParams.cPhi1->setConstant(kTRUE);
      _accepParams.dPhi1->setConstant(kTRUE);
      _accepParams.ePhi1->setConstant(kTRUE);

      _accepParams.aH1->setConstant(kTRUE);
      _accepParams.bH1->setConstant(kTRUE);
      _accepParams.cH1->setConstant(kTRUE);
      _accepParams.dH1->setConstant(kTRUE);
      _accepParams.eH1->setConstant(kTRUE);

      _accepParams.aH2->setConstant(kTRUE);
      _accepParams.bH2->setConstant(kTRUE);
      _accepParams.cH2->setConstant(kTRUE);
      _accepParams.dH2->setConstant(kTRUE);
      _accepParams.eH2->setConstant(kTRUE);

      _accepParams.aHs->setConstant(kTRUE);
      _accepParams.bHs->setConstant(kTRUE);
      _accepParams.cHs->setConstant(kTRUE);
      _accepParams.dHs->setConstant(kTRUE);
      _accepParams.eHs->setConstant(kTRUE);

      _accepParams.aM1->setConstant(kTRUE);
      _accepParams.bM1->setConstant(kTRUE);
      _accepParams.cM1->setConstant(kTRUE);
      _accepParams.dM1->setConstant(kTRUE);

      _accepParams.aM2->setConstant(kTRUE);
      _accepParams.bM2->setConstant(kTRUE);
      _accepParams.cM2->setConstant(kTRUE);
      _accepParams.dM2->setConstant(kTRUE);

    }
    else{
      _modelParams.a1Val->setConstant(kFALSE);
      _modelParams.phi1Val->setConstant(kFALSE);
      _modelParams.a2Val->setConstant(kFALSE);
      _modelParams.phi2Val->setConstant(kFALSE);
      _modelParams.a3Val->setConstant(kFALSE);
      _modelParams.phi3Val->setConstant(kFALSE);
      _modelParams.g1Val->setConstant(kFALSE);
      _modelParams.g2Val->setConstant(kFALSE);
      _modelParams.g3Val->setConstant(kFALSE);
      _modelParams.g4Val->setConstant(kFALSE);
      _modelParams.g1ValIm->setConstant(kFALSE);
      _modelParams.g2ValIm->setConstant(kFALSE);
      _modelParams.g3ValIm->setConstant(kFALSE);
      _modelParams.g4ValIm->setConstant(kFALSE);
      _modelParams.fa2->setConstant(kFALSE);
      _modelParams.fa3->setConstant(kFALSE);
      _modelParams.phia2->setConstant(kFALSE);
      _modelParams.phia3->setConstant(kFALSE);
      _modelParams.gamZ->setConstant(kFALSE);
      _modelParams.mZ->setConstant(kFALSE);
      _modelParams.R1Val->setConstant(kFALSE);
      _modelParams.R2Val->setConstant(kFALSE);

      _accepParams.aPhi->setConstant(kFALSE);
      _accepParams.bPhi->setConstant(kFALSE);
      _accepParams.cPhi->setConstant(kFALSE);
      _accepParams.dPhi->setConstant(kFALSE);
      _accepParams.ePhi->setConstant(kFALSE);

      _accepParams.aPhi1->setConstant(kFALSE);
      _accepParams.bPhi1->setConstant(kFALSE);
      _accepParams.cPhi1->setConstant(kFALSE);
      _accepParams.dPhi1->setConstant(kFALSE);
      _accepParams.ePhi1->setConstant(kFALSE);

      _accepParams.aH1->setConstant(kFALSE);
      _accepParams.bH1->setConstant(kFALSE);
      _accepParams.cH1->setConstant(kFALSE);
      _accepParams.dH1->setConstant(kFALSE);
      _accepParams.eH1->setConstant(kFALSE);

      _accepParams.aH2->setConstant(kFALSE);
      _accepParams.bH2->setConstant(kFALSE);
      _accepParams.cH2->setConstant(kFALSE);
      _accepParams.dH2->setConstant(kFALSE);
      _accepParams.eH2->setConstant(kFALSE);

      _accepParams.aHs->setConstant(kFALSE);
      _accepParams.bHs->setConstant(kFALSE);
      _accepParams.cHs->setConstant(kFALSE);
      _accepParams.dHs->setConstant(kFALSE);
      _accepParams.eHs->setConstant(kFALSE);

      _accepParams.aM1->setConstant(kFALSE);
      _accepParams.bM1->setConstant(kFALSE);
      _accepParams.cM1->setConstant(kFALSE);
      _accepParams.dM1->setConstant(kFALSE);

      _accepParams.aM2->setConstant(kFALSE);
      _accepParams.bM2->setConstant(kFALSE);
      _accepParams.cM2->setConstant(kFALSE);
      _accepParams.dM2->setConstant(kFALSE);

    }
    if (!acceptance_){
      _accepParams.aPhi->setConstant(kTRUE);
      _accepParams.bPhi->setConstant(kTRUE);
      _accepParams.cPhi->setConstant(kTRUE);
      _accepParams.dPhi->setConstant(kTRUE);
      _accepParams.ePhi->setConstant(kTRUE);

      _accepParams.aPhi1->setConstant(kTRUE);
      _accepParams.bPhi1->setConstant(kTRUE);
      _accepParams.cPhi1->setConstant(kTRUE);
      _accepParams.dPhi1->setConstant(kTRUE);
      _accepParams.ePhi1->setConstant(kTRUE);

      _accepParams.aH1->setConstant(kTRUE);
      _accepParams.bH1->setConstant(kTRUE);
      _accepParams.cH1->setConstant(kTRUE);
      _accepParams.dH1->setConstant(kTRUE);
      _accepParams.eH1->setConstant(kTRUE);

      _accepParams.aH2->setConstant(kTRUE);
      _accepParams.bH2->setConstant(kTRUE);
      _accepParams.cH2->setConstant(kTRUE);
      _accepParams.dH2->setConstant(kTRUE);
      _accepParams.eH2->setConstant(kTRUE);

      _accepParams.aHs->setConstant(kTRUE);
      _accepParams.bHs->setConstant(kTRUE);
      _accepParams.cHs->setConstant(kTRUE);
      _accepParams.dHs->setConstant(kTRUE);
      _accepParams.eHs->setConstant(kTRUE);

      _accepParams.aM1->setConstant(kTRUE);
      _accepParams.bM1->setConstant(kTRUE);
      _accepParams.cM1->setConstant(kTRUE);
      _accepParams.dM1->setConstant(kTRUE);

      _accepParams.aM2->setConstant(kTRUE);
      _accepParams.bM2->setConstant(kTRUE);
      _accepParams.cM2->setConstant(kTRUE);
      _accepParams.dM2->setConstant(kTRUE);
    }
  }

};




#endif



