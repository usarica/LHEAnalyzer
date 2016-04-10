#ifndef SPIN_PDF_FACTORY
#define SPIN_PDF_FACTORY

#ifdef _def_melatools_
#include <ZZMatrixElement/MELA/interface/RooSpin.h>
#else
#include "RooSpin.h"
#endif
#include "TString.h"
#include "RooFormulaVar.h"


class SpinPdfFactory{
public:
  RooSpin::modelMeasurables measurables;
  RooSpin::modelParameters parameters;

  SpinPdfFactory(RooSpin::modelMeasurables measurables_, RooSpin::VdecayType V1decay_=RooSpin::kVdecayType_Zll, RooSpin::VdecayType V2decay_=RooSpin::kVdecayType_Zll);
  virtual ~SpinPdfFactory();

  virtual void getMVGamV(Double_t* mV=0, Double_t* gamV=0) const;

  virtual void makeParamsConst(bool yesNo)=0;
  virtual void resetHypotheses()=0;
  virtual void resetVdecay(RooSpin::VdecayType V1decay_, RooSpin::VdecayType V2decay_);

protected:
  RooSpin* PDF_base;

  RooSpin::VdecayType V1decay;
  RooSpin::VdecayType V2decay;

  virtual void initMeasurables(RooSpin::modelMeasurables measurables_);

  virtual void initMassPole();
  virtual void initVdecayParams();

  virtual void destroyMassPole();
  virtual void destroyVdecayParams();

  virtual void initGVals()=0;
  virtual void destroyGVals()=0;
  virtual void initPDF()=0;
  virtual void destroyPDF()=0;
};




#endif



