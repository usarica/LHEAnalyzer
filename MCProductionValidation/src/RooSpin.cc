#ifdef _def_melatools_
#include <ZZMatrixElement/MELA/interface/RooSpin.h>
#else
#include "../include/RooSpin.h"
#endif


RooSpin::RooSpin(
  const char* name, const char* title,
  modelMeasurables _measurables,
  modelParameters _parameters,
  RooSpin::VdecayType _Vdecay1, RooSpin::VdecayType _Vdecay2
  ) : RooAbsPdf(name, title),

  Vdecay1(_Vdecay1), Vdecay2(_Vdecay2),

  h1("h1", "h1", this),
  h2("h2", "h2", this),
  Phi("Phi", "Phi", this),
  m1("m1", "m1", this),
  m2("m2", "m2", this),
  m12("m12", "m12", this),
  hs("hs", "hs", this),
  Phi1("Phi1", "Phi1", this),
  Y("Y", "Y", this),

  mX("mX", "mX", this, (RooAbsReal&)*(_parameters.mX)),
  gamX("gamX", "gamX", this, (RooAbsReal&)*(_parameters.gamX)),
  mV("mV", "mV", this, (RooAbsReal&)*(_parameters.mV)),
  gamV("gamV", "gamV", this, (RooAbsReal&)*(_parameters.gamV)),
  R1Val("R1Val", "R1Val", this, (RooAbsReal&)*(_parameters.R1Val)),
  R2Val("R2Val", "R2Val", this, (RooAbsReal&)*(_parameters.R2Val))
{
  setProxies(_measurables);
}


RooSpin::RooSpin(const RooSpin& other, const char* name) :
RooAbsPdf(other, name),

Vdecay1(other.Vdecay1), Vdecay2(other.Vdecay2),

h1("h1", this, other.h1),
h2("h2", this, other.h2),
Phi("Phi", this, other.Phi),
m1("m1", this, other.m1),
m2("m2", this, other.m2),
m12("m12", this, other.m12),
hs("hs", this, other.hs),
Phi1("Phi1", this, other.Phi1),
Y("Y", this, other.Y),

mX("mX", this, other.mX),
gamX("gamX", this, other.gamX),
mV("mV", this, other.mV),
gamV("gamV", this, other.gamV),
R1Val("R1Val", this, other.R1Val),
R2Val("R2Val", this, other.R2Val)
{}

void RooSpin::calculatePropagator(Double_t& propRe, Double_t& propIm, Double_t mass, bool useGamma)const{
  // prop = -i / ((m**2-mV**2) + i*mV*GaV) = - ( mV*GaV + i*(m**2-mV**2) ) / ((m**2-mV**2)**2 + (mV*GaV)**2)
  if (useGamma){
    propRe = 0;
    propIm = -1./pow(mass, 2);
  }
  else{
    Double_t denominator = pow(mV*gamV, 2)+pow(pow(mass, 2)-pow(mV, 2), 2);
    propRe = -mV*gamV/denominator;
    propIm = -(pow(mass, 2)-pow(mV, 2))/denominator;
  }
}
void RooSpin::setProxies(modelMeasurables _measurables){
  setProxy(h1, (RooAbsReal*)_measurables.h1);
  setProxy(h2, (RooAbsReal*)_measurables.h2);
  setProxy(Phi, (RooAbsReal*)_measurables.Phi);
  setProxy(m1, (RooAbsReal*)_measurables.m1);
  setProxy(m2, (RooAbsReal*)_measurables.m2);
  setProxy(m12, (RooAbsReal*)_measurables.m12);
  setProxy(hs, (RooAbsReal*)_measurables.hs);
  setProxy(Phi1, (RooAbsReal*)_measurables.Phi1);
  setProxy(Y, (RooAbsReal*)_measurables.Y);
}
void RooSpin::setProxy(RooRealProxy& proxy, RooAbsReal* objectPtr){
  if (objectPtr!=0) proxy.setArg((RooAbsReal&)*objectPtr);
}
