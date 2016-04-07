#ifndef TENSOR_PDF_FACTORY
#define TENSOR_PDF_FACTORY

#ifdef _def_melatools_
#include <ZZMatrixElement/MELA/interface/RooSpinTwo.h>
#else
#include "RooSpinTwo.h"
#endif
#include "TString.h"
#include "RooFormulaVar.h"


class TensorPdfFactory{
public:
  RooSpinTwo::modelMeasurables measurables;
  RooSpinTwo::modelParameters parameters;

  TensorPdfFactory(RooSpinTwo::modelMeasurables measurables_, int V1decay_=1, int V2decay_=1);
  virtual ~TensorPdfFactory();

  virtual void makeParamsConst(bool yesNo)=0;
  virtual void addHypothesis(int ig, double initval, double iphase=0);
  virtual void setTensorPolarization(int ig, double initval);
  virtual void resetHypotheses();
  virtual void resetVdecay(int V1decay_, int V2decay_);
  virtual RooSpinTwo* getPDF()=0;

protected:
  RooSpinTwo* PDF_base;

  int V1decay;
  int V2decay;

  virtual double getRValue(int Vdecay);

  virtual void initMeasurables(RooSpinTwo::modelMeasurables measurables_);
  virtual void initMassPole();
  virtual void initVdecayParams();
  virtual void initGVals();

  virtual void destroyMassPole();
  virtual void destroyVdecayParams();
  virtual void destroyGVals();

  virtual void initPDF()=0;
  virtual void destroyPDF()=0;
};




#endif



