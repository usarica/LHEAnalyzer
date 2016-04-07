#ifndef ROOSPIN
#define ROOSPIN

#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooCategoryProxy.h"
#include "RooAbsReal.h"
#include "RooRealVar.h"
#include "RooFormulaVar.h"
#include "RooAbsCategory.h"
#include "Riostream.h" 
#include <cmath>
#include "TMath.h"

using namespace TMath;

class RooSpin : public RooAbsPdf {
public:

  enum VdecayType{
    kVdecayType_Wany=-1,
    kVdecayType_Zll=1,
    kVdecayType_Znn=2,
    kVdecayType_Zuu=3,
    kVdecayType_Zdd=4,
    kVdecayType_Zud=5
  };

  struct modelMeasurables{
    RooRealVar* h1;
    RooRealVar* h2;
    RooRealVar* hs;
    RooRealVar* Phi;
    RooRealVar* Phi1;
    RooRealVar* m1;
    RooRealVar* m2;
    RooRealVar* m12;
    RooRealVar* Y;
  };
  struct modelParameters{
    RooRealVar* mX;
    RooRealVar* gamX;
    RooRealVar* mV;
    RooRealVar* gamV;
    RooRealVar* R1Val;
    RooRealVar* R2Val;
  };

  RooSpin(){};
  RooSpin(
    const char* name, const char* title,
    modelMeasurables _measurables,
    modelParameters _parameters,
    RooSpin::VdecayType _Vdecay1=RooSpin::kVdecayType_Zll, RooSpin::VdecayType _Vdecay2=RooSpin::kVdecayType_Zll
    );
  RooSpin(const RooSpin& other, const char* name=0);
  inline virtual ~RooSpin(){}

  virtual TObject* clone(const char* newname) const = 0;

  virtual Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName=0) const = 0;
  virtual Double_t analyticalIntegral(Int_t code, const char* rangeName=0) const = 0;

  virtual void setProxy(RooRealProxy& proxy, RooAbsReal* objectPtr);
  virtual void setDecayModes(RooSpin::VdecayType Vdecay1_, RooSpin::VdecayType Vdecay2_){ Vdecay1=Vdecay1_; Vdecay2=Vdecay2_; }

protected:

  enum{
    prime_h1=2,
    prime_h2=3,
    prime_hs=5,
    prime_Phi=7,
    prime_Phi1=11,
    prime_m1=13,
    prime_m2=17,
    prime_m12=19,
    prime_Y=23
  };

  RooSpin::VdecayType Vdecay1;
  RooSpin::VdecayType Vdecay2;

  RooRealProxy h1;
  RooRealProxy h2;
  RooRealProxy Phi;
  RooRealProxy m1;
  RooRealProxy m2;
  RooRealProxy m12;
  RooRealProxy hs;
  RooRealProxy Phi1;
  RooRealProxy Y;

  RooRealProxy mX;
  RooRealProxy gamX;
  RooRealProxy mV;
  RooRealProxy gamV;
  RooRealProxy R1Val;
  RooRealProxy R2Val;

  virtual Double_t evaluate() const = 0;

  virtual void calculatePropagator(Double_t& propRe, Double_t& propIm, Double_t mass, bool useGamma=false) const;
  virtual void setProxies(modelMeasurables _measurables);
};

#endif
