#ifndef SCALAR_PDF_FACTORY_GGH
#define SCALAR_PDF_FACTORY_GGH

#include "RooSpinZero_7DComplex_withAccep_ggH.h"
#include "ScalarPdfFactory.h"

class ScalarPdfFactory_ggH : public ScalarPdfFactory {
public:
  RooSpinZero_7DComplex_withAccep_ggH::accepParameters accepParams;

  ScalarPdfFactory_ggH(RooSpinZero::modelMeasurables measurables_, bool acceptance_=false, int V1decay_=1, int V2decay_=1);
  ScalarPdfFactory_ggH(RooSpinZero::modelMeasurables measurables_, double gRatio_[4][8], double gZGsRatio_[4][1], double gGsGsRatio_[3][1], bool pmf_applied_=false, bool acceptance_=false, int V1decay_=1, int V2decay_=1);
  ~ScalarPdfFactory_ggH();

  void makeParamsConst(bool yesNo=true);
  RooSpinZero* getPDF(){ return (RooSpinZero*)PDF; }

protected:
  RooSpinZero_7DComplex_withAccep_ggH* PDF;

  virtual void initAcceptanceParams();
  virtual void destroyAcceptanceParams();

  void initPDF();
  void destroyPDF(){ delete PDF; }
};


#endif



