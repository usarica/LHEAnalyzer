#ifndef SCALAR_PDF_FACTORY
#define SCALAR_PDF_FACTORY

#include "RooSpinZero.h"
#include "TString.h"
#include "RooFormulaVar.h"

class ScalarPdfFactory{
public:
  RooSpinZero::modelMeasurables measurables;
  RooSpinZero::modelParameters parameters;

  RooRealVar* g1Frac[7]; // f_a1 = 1.-sum_{i, i>1}{fabs(f_ai)}
  RooRealVar* g2Frac[8];
  RooRealVar* g3Frac[8];
  RooRealVar* g4Frac[8];
  RooRealVar* g1Phase[7]; // phi_a1=0
  RooRealVar* g2Phase[8];
  RooRealVar* g3Phase[8];
  RooRealVar* g4Phase[8];
  RooFormulaVar* gFracSum; // sum_{i, i>1}{fabs(f_ai)}
  RooFormulaVar* g1FracInterp[8]; // f_a1^interp = (f_a1<0 ? 0 : f_ai)
  RooFormulaVar* g2FracInterp[8];
  RooFormulaVar* g3FracInterp[8];
  RooFormulaVar* g4FracInterp[8];

  RooRealVar* gRatioVal[4][8];

  ScalarPdfFactory(RooSpinZero::modelMeasurables measurables_, bool acceptance_=false, int V1decay_=1, int V2decay_=1);
  ScalarPdfFactory(RooSpinZero::modelMeasurables measurables_, double gRatio_[4][8], bool pmf_applied_=false, bool acceptance_=false, int V1decay_=1, int V2decay_=1);
  virtual ~ScalarPdfFactory();

  virtual void makeParamsConst(bool yesNo)=0;
  virtual RooSpinZero* getPDF()=0;

protected:
  int parameterization;
  bool pmf_applied;
  bool acceptance;

  int V1decay;
  int V2decay;

  double gRatio[4][8];

  virtual void initMeasurables(RooSpinZero::modelMeasurables measurables_);
  virtual void initMassPole();
  virtual void initVdecayParams();
  virtual void initFractionsPhases();
  virtual void initGVals();

  virtual void destroyMassPole();
  virtual void destroyVdecayParams();
  virtual void destroyFractionsPhases();
  virtual void destroyGVals();

  virtual void initPDF()=0;
  virtual void destroyPDF()=0;
};




#endif



