#ifndef TENSOR_PDF_FACTORY_HZZ
#define TENSOR_PDF_FACTORY_HZZ

#include "RooSpinTwo_7DComplex_HVV.h"
#include "TensorPdfFactory.h"

class TensorPdfFactory_HVV : public TensorPdfFactory {
public:
  TensorPdfFactory_HVV(RooSpinTwo::modelMeasurables measurables_, int V1decay_=1, int V2decay_=1);
  ~TensorPdfFactory_HVV();

  void makeParamsConst(bool yesNo=true);
  RooSpinTwo* getPDF(){ return (RooSpinTwo*)PDF; }

protected:
  RooSpinTwo_7DComplex_HVV* PDF;

  void initPDF();
  void destroyPDF(){ delete PDF; }
};


#endif



