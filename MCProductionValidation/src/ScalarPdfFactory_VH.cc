#ifdef _def_melatools_
#include <ZZMatrixElement/MELA/interface/ScalarPdfFactory_VH.h>
#else
#include "../include/ScalarPdfFactory_VH.h"
#endif


ScalarPdfFactory_VH::ScalarPdfFactory_VH(RooSpinZero::modelMeasurables measurables_, double sqrts_, int VHmode1_, int VHmode2_) :
ScalarPdfFactory(measurables_, false, VHmode1_, VHmode2_),
sqrts(sqrts_)
{
  if (VHmode1_==-1 || VHmode1_==3 || VHmode1_==4 || VHmode1_==5) PDFType = 1;
  else PDFType = 2;
  if (PDFType==2) measurables.Y=0;

  makeParamsConst(true);
  initPDF();
}
ScalarPdfFactory_VH::ScalarPdfFactory_VH(RooSpinZero::modelMeasurables measurables_, double gRatio_[4][8], double gZGsRatio_[4][1], double gGsGsRatio_[3][1], double sqrts_, bool pmf_applied_, int VHmode1_, int VHmode2_) :
ScalarPdfFactory(measurables_, gRatio_, gZGsRatio_, gGsGsRatio_, pmf_applied_, false, VHmode1_, VHmode2_),
sqrts(sqrts_)
{
  if (VHmode1_==-1 || VHmode1_==3 || VHmode1_==4 || VHmode1_==5) PDFType = 1;
  else PDFType = 2;
  if (PDFType==2) measurables.Y=0;

  makeParamsConst(true);
  initPDF();
}
ScalarPdfFactory_VH::~ScalarPdfFactory_VH(){
  destroyPDF();
}

void ScalarPdfFactory_VH::makeParamsConst(bool yesNo){
  parameters.Lambda->setConstant(true);
  parameters.Lambda_zgs1->setConstant(true);
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
  if (PDFType==2){
    PDF_ILC_5D = new RooSpinZero_5D_VH(
      "PDF", "PDF",
      measurables,
      parameters,
      V1decay, V2decay
      );
    PDF_base = (RooSpinZero*)PDF_ILC_5D;
  }
  else if (PDFType==1){
    PDF_LHC_3D = new RooSpinZero_3D_pp_VH(
      "PDF", "PDF",
      measurables,
      parameters,
      sqrts,
      V1decay, V2decay
      );
    PDF_base = (RooSpinZero*)PDF_LHC_3D;
  }
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


