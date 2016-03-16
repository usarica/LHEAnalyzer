#ifndef SCALAR_PDF_FACTORY_VH
#define SCALAR_PDF_FACTORY_VH

#include "RooSpinZero_5D_VH.h"
#include "RooSpinZero_3D_pp_VH.h"
#include "ScalarPdfFactory.h"

class ScalarPdfFactory_VH : public ScalarPdfFactory {
public:

  ScalarPdfFactory_VH(RooSpinZero::modelMeasurables measurables_, double sqrts_, int VHmode1_=3, int VHmode2_=3);
  ScalarPdfFactory_VH(RooSpinZero::modelMeasurables measurables_, double gRatio_[4][8], double gZGsRatio_[4][1], double gGsGsRatio_[3][1], double sqrts_, bool pmf_applied_=false, int VHmode1_=3, int VHmode2_=3);
  ~ScalarPdfFactory_VH();

  void makeParamsConst(bool yesNo=true);
  RooSpinZero* getPDF();

protected:
  RooSpinZero_5D_VH* PDF_ILC_5D;
  RooSpinZero_3D_pp_VH* PDF_LHC_3D;
  double sqrts;
  int PDFType;

  void initPDF();
  void destroyPDF();
};

#endif



