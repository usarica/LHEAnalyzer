#include "../include/ScalarPdfFactory_VH.h"

ScalarPdfFactory_VH::ScalarPdfFactory_VH(RooSpinZero::modelMeasurables measurables_, double sqrts_, int PDFType_, int VHmode_) :
ScalarPdfFactory(measurables_, false, VHmode_, VHmode_),
sqrts(sqrts_),
PDFType(PDFType_)
{
  makeParamsConst(true);
  initPDF();
}
ScalarPdfFactory_VH::ScalarPdfFactory_VH(RooSpinZero::modelMeasurables measurables_, double gRatio_[4][8], double sqrts_, int PDFType_, bool pmf_applied_, int VHmode_) :
ScalarPdfFactory(measurables_, gRatio_, pmf_applied_, false, VHmode_, VHmode_),
sqrts(sqrts_),
PDFType(PDFType_)
{
  makeParamsConst(true);
  initPDF();
}
ScalarPdfFactory_VH::~ScalarPdfFactory_VH(){
  destroyPDF();
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
}

void ScalarPdfFactory_VH::initPDF(){
  PDF_ILC_5D=0;
  PDF_LHC_3D=0;
  if (PDFType==2) PDF_ILC_5D = new RooSpinZero_5D_VH(
    "PDF", "PDF",
    measurables,
    parameters
    );
  else if (PDFType==1) PDF_LHC_3D = new RooSpinZero_3D_pp_VH(
    "PDF", "PDF",
    measurables,
    parameters,
    sqrts
    );
}

RooSpinZero* ScalarPdfFactory_VH::getPDF(){
  if (PDFType==2) return (RooSpinZero*)PDF_ILC_5D;
  else if (PDFType==1) return (RooSpinZero*)PDF_LHC_3D;
  else return 0;
}

void ScalarPdfFactory_VH::destroyPDF(){
  if (PDF_ILC_5D!=0) delete PDF_ILC_5D;
  if (PDF_LHC_3D!=0) delete PDF_LHC_3D;
}


