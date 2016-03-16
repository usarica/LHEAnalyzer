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

  RooRealVar* gzgs1Frac[1];
  RooRealVar* gzgs2Frac[1];
  RooRealVar* gzgs3Frac[1];
  RooRealVar* gzgs4Frac[1];
  RooRealVar* gzgs1Phase[1];
  RooRealVar* gzgs2Phase[1];
  RooRealVar* gzgs3Phase[1];
  RooRealVar* gzgs4Phase[1];

  RooRealVar* ggsgs2Frac[1];
  RooRealVar* ggsgs3Frac[1];
  RooRealVar* ggsgs4Frac[1];
  RooRealVar* ggsgs2Phase[1];
  RooRealVar* ggsgs3Phase[1];
  RooRealVar* ggsgs4Phase[1];

  RooFormulaVar* gFracSum; // sum_{i, i>1}{fabs(f_ai)}
  RooFormulaVar* g1FracInterp[8]; // f_a1^interp = (f_a1<0 ? 0 : f_ai)
  RooFormulaVar* g2FracInterp[8];
  RooFormulaVar* g3FracInterp[8];
  RooFormulaVar* g4FracInterp[8];

  RooFormulaVar* gzgs1FracInterp[1];
  RooFormulaVar* gzgs2FracInterp[1];
  RooFormulaVar* gzgs3FracInterp[1];
  RooFormulaVar* gzgs4FracInterp[1];
  RooFormulaVar* ggsgs2FracInterp[1];
  RooFormulaVar* ggsgs3FracInterp[1];
  RooFormulaVar* ggsgs4FracInterp[1];

  RooRealVar* gRatioVal[4][8];
  RooRealVar* gZGsRatioVal[4][1];
  RooRealVar* gGsGsRatioVal[3][1];

  ScalarPdfFactory(RooSpinZero::modelMeasurables measurables_, bool acceptance_=false, int V1decay_=1, int V2decay_=1);
  ScalarPdfFactory(RooSpinZero::modelMeasurables measurables_, double gRatio_[4][8], double gZGsRatio_[4][1], double gGsGsRatio_[3][1], bool pmf_applied_=false, bool acceptance_=false, int V1decay_=1, int V2decay_=1);
  virtual ~ScalarPdfFactory();

  virtual void makeParamsConst(bool yesNo)=0;
  virtual void addHypothesis(int ig, int ilam, double iphase=0, double altparam_fracval=0);
  virtual void resetHypotheses();
  virtual void resetVdecay(int V1decay_, int V2decay_);
  virtual RooSpinZero* getPDF()=0;

protected:
  RooSpinZero* PDF_base;

  int parameterization;
  bool pmf_applied;
  bool acceptance;

  int V1decay;
  int V2decay;

  double gRatio[4][8];
  double gZGsRatio[4][1];
  double gGsGsRatio[3][1];

  virtual double getRValue(int Vdecay);

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



