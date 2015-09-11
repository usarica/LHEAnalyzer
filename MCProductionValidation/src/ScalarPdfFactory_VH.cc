#include "../include/RooSpinZero_5D_VH.h"
#include "../include/RooSpinZero_3D_withAccep_VH.h"
#include "RooRealVar.h"
#include "TF1.h"
#include "TMath.h"

ScalarPdfFactory_VH::ScalarPdfFactory_VH(RooSpinZero::modelMeasurables measurables_, int PDFType_, bool acceptance_, int V1decay_, int V2decay_) :
ScalarPdfFactory(measurables_, acceptance_, V1decay_, V2decay_),
PDFType(PDFType_)
{
  initAcceptanceParams();
  makeParamsConst(true);
  initPDF();
}
ScalarPdfFactory_VH::ScalarPdfFactory_VH(RooSpinZero::modelMeasurables measurables_, double gRatio_[4][8], int PDFType_, bool pmf_applied_, bool acceptance_, int V1decay_, int V2decay_) :
ScalarPdfFactory(measurables_, gRatio_, pmf_applied_, acceptance_, V1decay_, V2decay_),
PDFType(PDFType_)
{
  initAcceptanceParams();
  makeParamsConst(true);
  initPDF();
}
ScalarPdfFactory_VH::~ScalarPdfFactory_VH(){
  destroyPDF();
  destroyAcceptanceParams();
}

void ScalarPdfFactory_VH::initAcceptanceParams(){
  if (acceptance){
    accepParams.b2 = new RooRealVar("b2", "b2", -7.41520e-02, -10, 10);
    accepParams.cgaus = new RooRealVar("cgaus", "cgaus", -3.04678e-01, -10, 10);
    accepParams.sgaus = new RooRealVar("sgaus", "sgaus", 4.97867e-02, 0, 1.);
  }
  else{
    accepParams.b2 = new RooRealVar("b2", "b2", 0);
    accepParams.cgaus = new RooRealVar("cgaus", "cgaus", 0);
    accepParams.sgaus = new RooRealVar("sgaus", "sgaus", 1.);
  }
}
void ScalarPdfFactory_VH::destroyAcceptanceParams(){
  delete accepParams.b2;
  delete accepParams.cgaus;
  delete accepParams.sgaus;
}

void ScalarPdfFactory_VH::makeParamsConst(bool yesNo){
  parameters.Lambda->setConstant(true);
  parameters.Lambda_z1->setConstant(true);
  parameters.Lambda_z2->setConstant(true);
  parameters.Lambda_z3->setConstant(true);
  parameters.Lambda_z4->setConstant(true);
  parameters.Lambda_Q->setConstant(true);

  parameters.mX->setConstant(yesNo);
  parameters.gamX->setConstant(yesNo);
  parameters.mV->setConstant(yesNo);
  parameters.gamV->setConstant(yesNo);
  parameters.R1Val->setConstant(yesNo);
  parameters.R2Val->setConstant(yesNo);

  if (acceptance && !yesNo){
    accepParams.b2->setConstant(kFALSE);
    accepParams.cgaus->setConstant(kFALSE);
    accepParams.sgaus->setConstant(kFALSE);
  }
  else{
    accepParams.b2->setConstant(kTRUE);
    accepParams.cgaus->setConstant(kTRUE);
    accepParams.sgaus->setConstant(kTRUE);
  }
}

void ScalarPdfFactory_VH::initPDF(){
  PDF_ILC_5D=0;
  PDF_ILC_3D=0;
  if (PDFType==1) PDF_ILC_5D = new RooSpinZero_5D_VH(
    "PDF", "PDF",
    measurables,
    parameters
    );
  else if (PDFType==2) PDF_ILC_3D = new RooSpinZero_3D_withAccep_VH(
    "PDF", "PDF",
    measurables,
    parameters,
    accepParams,
    acceptance
    );
}
void ScalarPdfFactory_VH::destroyPDF(){
  if (PDF_ILC_5D!=0) delete PDF_ILC_5D;
  if (PDF_ILC_3D!=0) delete PDF_ILC_3D;
}


