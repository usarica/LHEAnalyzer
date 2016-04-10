#ifdef _def_melatools_
#include <ZZMatrixElement/MELA/interface/RooSpin.h>
#else
#include "../include/RooSpin.h"
#endif


using namespace std;

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
  mW("mW", "mW", this, (RooAbsReal&)*(_parameters.mW)),
  gamW("gamW", "gamW", this, (RooAbsReal&)*(_parameters.gamW)),
  mZ("mZ", "mZ", this, (RooAbsReal&)*(_parameters.mZ)),
  gamZ("gamZ", "gamZ", this, (RooAbsReal&)*(_parameters.gamZ)),
  Sin2ThetaW("Sin2ThetaW", "Sin2ThetaW", this, (RooAbsReal&)*(_parameters.Sin2ThetaW)),
  vev("vev", "vev", this, (RooAbsReal&)*(_parameters.vev))
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
mW("mW", this, other.mW),
gamW("gamW", this, other.gamW),
mZ("mZ", this, other.mZ),
gamZ("gamZ", this, other.gamZ),
Sin2ThetaW("Sin2ThetaW", this, other.Sin2ThetaW),
vev("vev", this, other.vev)
{}

void RooSpin::calculatePropagator(Double_t& propRe, Double_t& propIm, Double_t mass, Int_t propType)const{
  // prop = -i / ((m**2-mV**2) + i*mV*GaV) = - ( mV*GaV + i*(m**2-mV**2) ) / ((m**2-mV**2)**2 + (mV*GaV)**2)
  if (propType==0){
    propRe = 0.;
    propIm = (mass!=0. ? -1./pow(mass, 2) : 0.);
  }
  else if (propType==1){
    Double_t mV, gamV;
    getMVGamV(&mV, &gamV);
    if (gamV>0){
      Double_t denominator = pow(mV*gamV, 2)+pow(pow(mass, 2)-pow(mV, 2), 2);
      propRe = -mV*gamV/denominator;
      propIm = -(pow(mass, 2)-pow(mV, 2))/denominator;
    }
    else{
      propRe = (mass==mV ? 1. : 0.);
      propIm = 0.;
    }
  }
  else if (propType==2){ // Higgs prop = i / ((m**2-mX**2) + i*mX*GaX) = - ( mX*GaX + i*(m**2-mX**2) ) / ((m**2-mX**2)**2 + (mX*GaX)**2)
    if (gamX>0.){
      Double_t denominator = pow(mX*gamX, 2)+pow(pow(mass, 2)-pow(mX, 2), 2);
      propRe = mX*gamX/denominator;
      propIm = (pow(mass, 2)-pow(mX, 2))/denominator;
    }
    else{
      propRe = (mass==mX ? 1. : 0.);
      propIm = 0.;
    }
  }
  else{
    propRe = 1.;
    propIm = 0.;
  }
}
void RooSpin::calculateGVGA(Double_t& gV, Double_t& gA, RooSpin::VdecayType Vdecay, bool isGamma)const{
  const Double_t atomicT3 = 0.5;
  const Double_t atomicCharge = 1.;

  const Double_t gW = 2.*mW/vev;
  const Double_t overallFactorZ = gW*0.5/sqrt(1.-Sin2ThetaW);
  const Double_t overallFactorGamma = gW*sqrt(Sin2ThetaW);
  const Double_t overallFactorW = gW*0.5/sqrt(2.);

  const Double_t Q_up = 2.*atomicCharge/3.;
  const Double_t Q_dn = -atomicCharge/3.;
  const Double_t Q_l = -atomicCharge;
  const Double_t Q_nu = 0;

  // gV = T3 - 2*Qf*sintW**2
  const Double_t gV_up = atomicT3 - 2.*Q_up*Sin2ThetaW;
  const Double_t gV_dn = -atomicT3 - 2.*Q_dn*Sin2ThetaW;
  const Double_t gV_l = -atomicT3 - 2.*Q_l*Sin2ThetaW;
  const Double_t gV_nu = atomicT3 - 2.*Q_nu*Sin2ThetaW;

  // gA = T3
  const Double_t gA_up = atomicT3;
  const Double_t gA_dn = -atomicT3;
  const Double_t gA_l = -atomicT3;
  const Double_t gA_nu = atomicT3;

  if (Vdecay==RooSpin::kVdecayType_Zud){
    if (!isGamma){
      gV = overallFactorZ*(2.*gV_up + 3.*gV_dn)/5.;
      gA = overallFactorZ*(2.*gA_up + 3.*gA_dn)/5.;
    }
    else{
      gV = overallFactorGamma*(2.*Q_up + 3.*Q_dn)/5.;
      gA = 0;
    }
  }
  else if (Vdecay==RooSpin::kVdecayType_Zdd){
    if (!isGamma){
      gV = overallFactorZ*gV_dn;
      gA = overallFactorZ*gA_dn;
    }
    else{
      gV = overallFactorGamma*Q_dn;
      gA = 0;
    }
  }
  else if (Vdecay==RooSpin::kVdecayType_Zuu){
    if (!isGamma){
      gV = overallFactorZ*gV_up;
      gA = overallFactorZ*gA_up;
    }
    else{
      gV = overallFactorGamma*Q_up;
      gA = 0;
    }
  }
  else if (Vdecay==RooSpin::kVdecayType_Znn){
    if (!isGamma){
      gV = overallFactorZ*gV_nu;
      gA = overallFactorZ*gA_nu;
    }
    else{
      gV = overallFactorGamma*Q_nu;
      gA = 0;
    }
  }
  else if (Vdecay==RooSpin::kVdecayType_Zll){
    if (!isGamma){
      gV = overallFactorZ*gV_l;
      gA = overallFactorZ*gA_l;
    }
    else{
      gV = overallFactorGamma*Q_l;
      gA = 0;
    }
  }
  else if (Vdecay==RooSpin::kVdecayType_Wany){
    if (!isGamma){
      gV = overallFactorW;
      gA = overallFactorW;
    }
    else{
      gV = 0;
      gA = 0;
    }
  }
  else{
    gV = 1;
    gA = 0;
  }
}
void RooSpin::calculateR1R2(Double_t& R1Val, Double_t& R2Val, bool isGammaV1, bool isGammaV2)const{
  Double_t gV1, gV2, gA1, gA2;
  calculateGVGA(gV1, gA1, Vdecay1, isGammaV1);
  R1Val = 2.*gV1*gA1/(pow(gV1, 2) + pow(gA1, 2));
  calculateGVGA(gV2, gA2, Vdecay2, isGammaV2);
  R2Val = 2.*gV2*gA2/(pow(gV2, 2) + pow(gA2, 2));
}
void RooSpin::getMVGamV(Double_t* mV, Double_t* gamV)const{
  if (Vdecay1==RooSpin::kVdecayType_Wany){
    if (mV!=0) (*mV)=mW;
    if (gamV!=0) (*gamV)=gamW;
  }
  else if (!(Vdecay1==RooSpin::kVdecayType_GammaOnshell && Vdecay2==RooSpin::kVdecayType_GammaOnshell)){
    if (mV!=0) (*mV)=mZ;
    if (gamV!=0) (*gamV)=gamZ;
  }
  else{
    if (mV!=0) (*mV)=0;
    if (gamV!=0) (*gamV)=0;
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
Bool_t RooSpin::checkFundamentalType(const RooRealProxy& proxy)const{
  RooAbsArg* arg = proxy.absArg();
  return (dynamic_cast<RooRealVar*>(arg)!=0);
}
