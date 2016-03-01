#ifndef ROOSPINTWO_7DCOMPLEX_HVV
#define ROOSPINTWO_7DCOMPLEX_HVV

#include "RooSpinTwo.h"


class RooSpinTwo_7DComplex_HVV : public RooSpinTwo {

public:

  RooSpinTwo_7DComplex_HVV(){}
  RooSpinTwo_7DComplex_HVV(
    const char *name, const char *title,
    modelMeasurables _measurables,
    modelParameters _parameters
    );
  RooSpinTwo_7DComplex_HVV(const RooSpinTwo_7DComplex_HVV& other, const char* name=0);
  virtual TObject* clone(const char* newname) const { return new RooSpinTwo_7DComplex_HVV(*this, newname); }
  inline virtual ~RooSpinTwo_7DComplex_HVV(){}

  Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName=0) const;
  Double_t analyticalIntegral(Int_t code, const char* rangeName=0) const;

protected:

  Double_t evaluate() const;

};

#endif
