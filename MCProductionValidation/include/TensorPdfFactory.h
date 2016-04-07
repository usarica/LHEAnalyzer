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
  RooSpin::modelMeasurables measurables;
  RooSpin::modelParameters parameters;
  RooSpinTwo::modelCouplings couplings;

  TensorPdfFactory(RooSpinTwo::modelMeasurables measurables_, RooSpin::VdecayType V1decay_=RooSpin::kVdecayType_Zll, RooSpin::VdecayType V2decay_=RooSpin::kVdecayType_Zll);
  virtual ~TensorPdfFactory();

  virtual void makeParamsConst(bool yesNo)=0;
  virtual void addHypothesis(int ig, double initval, double iphase=0);
  virtual void setTensorPolarization(int ig, double initval);
  virtual void resetHypotheses();
  virtual void resetVdecay(RooSpin::VdecayType V1decay_, RooSpin::VdecayType V2decay_);
  virtual RooSpinTwo* getPDF()=0;

protected:
  RooSpinTwo* PDF_base;

  RooSpin::VdecayType V1decay;
  RooSpin::VdecayType V2decay;

  virtual double getRValue(RooSpin::VdecayType Vdecay);

  virtual void initMeasurables(RooSpin::modelMeasurables measurables_);
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



