#ifndef SCALAR_PDF_FACTORY_VH
#define SCALAR_PDF_FACTORY_VH

#include "RooSpinZero_5D_VH.h"
#include "RooSpinZero_3D_withAccep_VH.h"
#include "ScalarPdfFactory.h"

class ScalarPdfFactory_VH : public ScalarPdfFactory {
public:
  RooSpinZero_3D_withAccep_VH::accepParameters accepParams;

  ScalarPdfFactory_VH(RooSpinZero::modelMeasurables measurables_, int PDFType_=0, bool acceptance_=false, int V1decay_=1, int V2decay_=1);
  ScalarPdfFactory_VH(RooSpinZero::modelMeasurables measurables_, double gRatio_[4][8], int PDFType_=0, bool pmf_applied_=false, bool acceptance_=false, int V1decay_=1, int V2decay_=1);
  ~ScalarPdfFactory_VH();

  void makeParamsConst(bool yesNo=true);
  RooSpinZero* getPDF();

protected:
  RooSpinZero_5D_VH* PDF_ILC_5D;
  RooSpinZero_3D_withAccep_VH* PDF_ILC_3D;
  int PDFType;

  virtual void initAcceptanceParams();
  virtual void destroyAcceptanceParams();

  void initPDF();
  void destroyPDF();
};




#endif



