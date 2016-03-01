#include "../include/RooSpinTwo_7DComplex_HVV.h"


RooSpinTwo_7DComplex_HVV::RooSpinTwo_7DComplex_HVV(
  const char *name, const char *title,
  modelMeasurables _measurables,
  modelParameters _parameters
  ) : RooSpinTwo(
  name, title,
  _measurables,
  _parameters
  )
{}


RooSpinTwo_7DComplex_HVV::RooSpinTwo_7DComplex_HVV(
  const RooSpinTwo_7DComplex_HVV& other, const char* name
  ) : RooSpinTwo(other, name)
{}



Double_t RooSpinTwo_7DComplex_HVV::evaluate() const{
  bool isZZ = true;
  if (mV < 90.) isZZ = false;
  if (isZZ) {
    if ((m1+m2) > m12 || fabs(m2-mV)<fabs(m1-mV) || m2 <= 0.0 || m1 <= 0.0) return 1e-15;
  }
  else if ((m1+m2) > m12 || m2 <= 0.0 || m1 <= 0.0) return 1e-15;

  Double_t betaValSq = (1.-(pow(m1-m2, 2)/pow(m12, 2)))*(1.-(pow(m1+m2, 2)/pow(m12, 2)));
  if (betaValSq<0) return 1e-15;
  Double_t betaVal = sqrt(betaValSq);

  Double_t term1Coeff = pow(m1, 3)/(pow(pow(m1, 2)-pow(mV, 2), 2)+pow(mV*gamV, 2));
  Double_t term2Coeff = pow(m2, 3)/(pow(pow(m2, 2)-pow(mV, 2), 2)+pow(mV*gamV, 2));
  Double_t f_spinz0 = 1. - f_spinz1 - f_spinz2;

  Double_t hsneg = -hs;
  Double_t value = 0;

  Double_t 
    A00Re, A00Im,
    AppRe, AppIm,
    A0pRe, A0pIm, Ap0Re, Ap0Im,
    AmmRe, AmmIm, A0mRe, A0mIm, Am0Re, Am0Im,
    ApmRe, ApmIm, AmpRe, AmpIm;

  calculateAmplitudes(
    A00Re, A00Im,
    AppRe, AppIm,
    A0pRe, A0pIm, Ap0Re, Ap0Im,
    AmmRe, AmmIm, A0mRe, A0mIm, Am0Re, Am0Im,
    ApmRe, ApmIm, AmpRe, AmpIm
    );

  Double_t A00 = A00Im*A00Im + A00Re*A00Re;
  Double_t App = AppIm*AppIm + AppRe*AppRe;
  Double_t Amm = AmmIm*AmmIm + AmmRe*AmmRe;
  Double_t Ap0 = Ap0Im*Ap0Im + Ap0Re*Ap0Re;
  Double_t A0p = A0pIm*A0pIm + A0pRe*A0pRe;
  Double_t Am0 = Am0Im*Am0Im + Am0Re*Am0Re;
  Double_t A0m = A0mIm*A0mIm + A0mRe*A0mRe;
  Double_t Apm=  ApmIm*ApmIm + ApmRe*ApmRe;
  Double_t Amp = AmpIm*AmpIm + AmpRe*AmpRe;

  Double_t phi00=atan2(A00Im, A00Re);
  Double_t phipp=atan2(AppIm, AppRe)-phi00;
  Double_t phimm=atan2(AmmIm, AmmRe)-phi00;
  Double_t phip0=atan2(Ap0Im, Ap0Re)-phi00;
  Double_t phi0p=atan2(A0pIm, A0pRe)-phi00;
  Double_t phim0=atan2(Am0Im, Am0Re)-phi00;
  Double_t phi0m=atan2(A0mIm, A0mRe)-phi00;
  Double_t phipm=atan2(ApmIm, ApmRe)-phi00;
  Double_t phimp=atan2(AmpIm, AmpRe)-phi00;

  value += (A00*(-1 + pow(h1, 2))*(-1 + pow(h2, 2))*(2*f_spinz0 + 3*f_spinz2 - 6*(2*f_spinz0 - 2*f_spinz1 + f_spinz2)*pow(hsneg, 2) +3*(6*f_spinz0 - 4*f_spinz1 + f_spinz2)*pow(hsneg, 4)))/32.;
  value += (Amm*(2*f_spinz0 + 3*f_spinz2 - 6*(2*f_spinz0 - 2*f_spinz1 + f_spinz2)*pow(hsneg, 2) +3*(6*f_spinz0 - 4*f_spinz1 + f_spinz2)*pow(hsneg, 4))*(-1 - pow(h1, 2) + 2*h1*R1Val)*(-1 - pow(h2, 2) + 2*h2*R2Val))/128.;
  value += (App*(2*f_spinz0 + 3*f_spinz2 - 6*(2*f_spinz0 - 2*f_spinz1 + f_spinz2)*pow(hsneg, 2) +3*(6*f_spinz0 - 4*f_spinz1 + f_spinz2)*pow(hsneg, 4))*(1 + pow(h1, 2) + 2*h1*R1Val)*(1 + pow(h2, 2) + 2*h2*R2Val))/128.;
  value += -(Ap0*(-1 + pow(h2, 2))*(f_spinz1 + f_spinz2 - 3*(-2*f_spinz0 + f_spinz1)*pow(hsneg, 2) -(6*f_spinz0 - 4*f_spinz1 + f_spinz2)*pow(hsneg, 4))*(1 + pow(h1, 2) + 2*h1*R1Val))/32.;
  value += (A0m*(-1 + pow(h1, 2))*(f_spinz1 + f_spinz2 - 3*(-2*f_spinz0 + f_spinz1)*pow(hsneg, 2) -(6*f_spinz0 - 4*f_spinz1 + f_spinz2)*pow(hsneg, 4))*(-1 - pow(h2, 2) + 2*h2*R2Val))/32.;
  value += -(A0p*(-1 + pow(h1, 2))*(f_spinz1 + f_spinz2 - 3*(-2*f_spinz0 + f_spinz1)*pow(hsneg, 2) -(6*f_spinz0 - 4*f_spinz1 + f_spinz2)*pow(hsneg, 4))*(1 + pow(h2, 2) + 2*h2*R2Val))/32.;
  value += (Am0*(-1 + pow(h2, 2))*(f_spinz1 + f_spinz2 - 3*(-2*f_spinz0 + f_spinz1)*pow(hsneg, 2) -(6*f_spinz0 - 4*f_spinz1 + f_spinz2)*pow(hsneg, 4))*(-1 - pow(h1, 2) + 2*h1*R1Val))/32.;
  value += -(Apm*(4*f_spinz1 + f_spinz2 + 6*f_spinz2*pow(hsneg, 2) + (-4*f_spinz1 + f_spinz2)*pow(hsneg, 4) +6*f_spinz0*pow(-1 + pow(hsneg, 2), 2))*(1 + pow(h1, 2) + 2*h1*R1Val)*(-1 - pow(h2, 2) + 2*h2*R2Val))/256.;
  value += -(Amp*(4*f_spinz1 + f_spinz2 + 6*f_spinz2*pow(hsneg, 2) + (-4*f_spinz1 + f_spinz2)*pow(hsneg, 4) +6*f_spinz0*pow(-1 + pow(hsneg, 2), 2))*(-1 - pow(h1, 2) + 2*h1*R1Val)*(1 + pow(h2, 2) + 2*h2*R2Val))/256.;
  value += (sqrt(A00)*sqrt(Amm)*sqrt(1 - pow(h1, 2))*sqrt(1 - pow(h2, 2))*(2*f_spinz0 + 3*f_spinz2 - 6*(2*f_spinz0 - 2*f_spinz1 + f_spinz2)*pow(hsneg, 2) +3*(6*f_spinz0 - 4*f_spinz1 + f_spinz2)*pow(hsneg, 4))*(h1 - R1Val)*(h2 - R2Val)*cos(Phi - phimm))/32.;
  value += (sqrt(A00)*sqrt(App)*sqrt(1 - pow(h1, 2))*sqrt(1 - pow(h2, 2))*(2*f_spinz0 + 3*f_spinz2 - 6*(2*f_spinz0 - 2*f_spinz1 + f_spinz2)*pow(hsneg, 2) +3*(6*f_spinz0 - 4*f_spinz1 + f_spinz2)*pow(hsneg, 4))*(h1 + R1Val)*(h2 + R2Val)*cos(Phi + phipp))/32.;
  value += -(sqrt(3)*sqrt(A00)*sqrt(Ap0)*sqrt(1 - pow(h1, 2))*(-1 + pow(h2, 2))*hsneg*sqrt(1 - pow(hsneg, 2))*(-2*f_spinz0 + 2*f_spinz1 - f_spinz2 +(6*f_spinz0 - 4*f_spinz1 + f_spinz2)*pow(hsneg, 2))*(h1 + R1Val)*cos(Phi1 - phip0))/16.;
  value += (sqrt(3)*sqrt(A00)*sqrt(A0m)*(-1 + pow(h1, 2))*sqrt(1 - pow(h2, 2))*hsneg*sqrt(1 - pow(hsneg, 2))*(-2*f_spinz0 + 2*f_spinz1 - f_spinz2 +(6*f_spinz0 - 4*f_spinz1 + f_spinz2)*pow(hsneg, 2))*(h2 - R2Val)*cos(Phi - phi0m + Phi1))/16.;
  value += (sqrt(3)*sqrt(A00)*sqrt(A0p)*(-1 + pow(h1, 2))*sqrt(1 - pow(h2, 2))*hsneg*sqrt(1 - pow(hsneg, 2))*(-2*f_spinz0 + 2*f_spinz1 - f_spinz2 +(6*f_spinz0 - 4*f_spinz1 + f_spinz2)*pow(hsneg, 2))*(h2 + R2Val)*cos(Phi + phi0p + Phi1))/16.;
  value += (sqrt(3)*sqrt(A00)*sqrt(Am0)*sqrt(1 - pow(h1, 2))*(-1 + pow(h2, 2))*hsneg*sqrt(1 - pow(hsneg, 2))*(-2*f_spinz0 + 2*f_spinz1 - f_spinz2 +(6*f_spinz0 - 4*f_spinz1 + f_spinz2)*pow(hsneg, 2))*(-h1 + R1Val)*cos(Phi1 + phim0))/16.;
  value += (sqrt(1.5)*sqrt(A00)*sqrt(Apm)*sqrt(1 - pow(h1, 2))*sqrt(1 - pow(h2, 2))*(-1 + hsneg)*(1 + hsneg)*(-2*f_spinz0 + f_spinz2 + (6*f_spinz0 - 4*f_spinz1 + f_spinz2)*pow(hsneg, 2))*(h1 + R1Val)*(h2 - R2Val)*cos(Phi + 2*Phi1 - phipm))/32.;
  value += (sqrt(1.5)*sqrt(A00)*sqrt(Amp)*sqrt(1 - pow(h1, 2))*sqrt(1 - pow(h2, 2))*(-1 + hsneg)*(1 + hsneg)*(-2*f_spinz0 + f_spinz2 + (6*f_spinz0 - 4*f_spinz1 + f_spinz2)*pow(hsneg, 2))*(h1 - R1Val)*(h2 + R2Val)*cos(Phi + 2*Phi1 + phimp))/32.;
  value += (sqrt(Amm)*sqrt(App)*(-1 + pow(h1, 2))*(-1 + pow(h2, 2))*(2*f_spinz0 + 3*f_spinz2 - 6*(2*f_spinz0 - 2*f_spinz1 + f_spinz2)*pow(hsneg, 2) +3*(6*f_spinz0 - 4*f_spinz1 + f_spinz2)*pow(hsneg, 4))*cos(2*Phi - phimm + phipp))/64.;
  value += (sqrt(3)*sqrt(Amm)*sqrt(Ap0)*(-1 + pow(h1, 2))*sqrt(1 - pow(h2, 2))*hsneg*sqrt(1 - pow(hsneg, 2))*(-2*f_spinz0 + 2*f_spinz1 - f_spinz2 +(6*f_spinz0 - 4*f_spinz1 + f_spinz2)*pow(hsneg, 2))*(h2 - R2Val)*cos(Phi - Phi1 - phimm + phip0))/32.;
  value += (sqrt(3)*sqrt(A0m)*sqrt(Amm)*sqrt(1 - pow(h1, 2))*hsneg*sqrt(1 - pow(hsneg, 2))*(-2*f_spinz0 + 2*f_spinz1 - f_spinz2 + (6*f_spinz0 - 4*f_spinz1 + f_spinz2)*pow(hsneg, 2))*(h1 - R1Val)*(-1 - pow(h2, 2) + 2*h2*R2Val)*cos(phi0m - Phi1 - phimm))/32.;
  value += (sqrt(3)*sqrt(A0p)*sqrt(Amm)*sqrt(1 - pow(h1, 2))*(-1 + pow(h2, 2))*hsneg*sqrt(1 - pow(hsneg, 2))*(-2*f_spinz0 + 2*f_spinz1 - f_spinz2 +(6*f_spinz0 - 4*f_spinz1 + f_spinz2)*pow(hsneg, 2))*(-h1 + R1Val)*cos(2*Phi + phi0p + Phi1 - phimm))/32.;
  value += (sqrt(3)*sqrt(Am0)*sqrt(Amm)*sqrt(1 - pow(h2, 2))*hsneg*sqrt(1 - pow(hsneg, 2))*(-2*f_spinz0 + 2*f_spinz1 - f_spinz2 + (6*f_spinz0 - 4*f_spinz1 + f_spinz2)*pow(hsneg, 2))*(-1 - pow(h1, 2) + 2*h1*R1Val)*(-h2 + R2Val)*cos(Phi + Phi1 + phim0 - phimm))/32.;
  value += -(sqrt(1.5)*sqrt(Amm)*sqrt(Apm)*(-1 + pow(h1, 2))*(-1 + hsneg)*(1 + hsneg)*(-2*f_spinz0 + f_spinz2 + (6*f_spinz0 - 4*f_spinz1 + f_spinz2)*pow(hsneg, 2))*(-1 - pow(h2, 2) + 2*h2*R2Val)*cos(2*Phi1 + phimm - phipm))/64.;
  value += -(sqrt(1.5)*sqrt(Amm)*sqrt(Amp)*(-1 + pow(h2, 2))*(-1 + hsneg)*(1 + hsneg)*(-2*f_spinz0 + f_spinz2 + (6*f_spinz0 - 4*f_spinz1 + f_spinz2)*pow(hsneg, 2))*(-1 - pow(h1, 2) + 2*h1*R1Val)*cos(2*Phi + 2*Phi1 - phimm + phimp))/64.;
  value += (sqrt(3)*sqrt(Ap0)*sqrt(App)*sqrt(1 - pow(h2, 2))*hsneg*sqrt(1 - pow(hsneg, 2))*(-2*f_spinz0 + 2*f_spinz1 - f_spinz2 + (6*f_spinz0 - 4*f_spinz1 + f_spinz2)*pow(hsneg, 2))*(1 + pow(h1, 2) + 2*h1*R1Val)*(h2 + R2Val)*cos(Phi + Phi1 - phip0 + phipp))/32.;
  value += -(sqrt(3)*sqrt(A0m)*sqrt(App)*sqrt(1 - pow(h1, 2))*(-1 + pow(h2, 2))*hsneg*sqrt(1 - pow(hsneg, 2))*(-2*f_spinz0 + 2*f_spinz1 - f_spinz2 +(6*f_spinz0 - 4*f_spinz1 + f_spinz2)*pow(hsneg, 2))*(h1 + R1Val)*cos(2*Phi - phi0m + Phi1 + phipp))/32.;
  value += -(sqrt(3)*sqrt(A0p)*sqrt(App)*sqrt(1 - pow(h1, 2))*hsneg*sqrt(1 - pow(hsneg, 2))*(-2*f_spinz0 + 2*f_spinz1 - f_spinz2 +(6*f_spinz0 - 4*f_spinz1 + f_spinz2)*pow(hsneg, 2))*(h1 + R1Val)*(1 + pow(h2, 2) + 2*h2*R2Val)*cos(phi0p + Phi1 - phipp))/32.;
  value += (sqrt(3)*sqrt(Am0)*sqrt(App)*(-1 + pow(h1, 2))*sqrt(1 - pow(h2, 2))*hsneg*sqrt(1 - pow(hsneg, 2))*(-2*f_spinz0 + 2*f_spinz1 - f_spinz2 +(6*f_spinz0 - 4*f_spinz1 + f_spinz2)*pow(hsneg, 2))*(h2 + R2Val)*cos(Phi - Phi1 - phim0 + phipp))/32.;
  value += (sqrt(1.5)*sqrt(Apm)*sqrt(App)*(-1 + pow(h2, 2))*(-1 + hsneg)*(1 + hsneg)*(-2*f_spinz0 + f_spinz2 + (6*f_spinz0 - 4*f_spinz1 + f_spinz2)*pow(hsneg, 2))*(1 + pow(h1, 2) + 2*h1*R1Val)*cos(2*Phi + 2*Phi1 - phipm + phipp))/64.;
  value += (sqrt(1.5)*sqrt(Amp)*sqrt(App)*(-1 + pow(h1, 2))*(-1 + hsneg)*(1 + hsneg)*(-2*f_spinz0 + f_spinz2 + (6*f_spinz0 - 4*f_spinz1 + f_spinz2)*pow(hsneg, 2))*(1 + pow(h2, 2) + 2*h2*R2Val)*cos(2*Phi1 + phimp - phipp))/64.;
  value += -(sqrt(A0m)*sqrt(Ap0)*sqrt(1 - pow(h1, 2))*sqrt(1 - pow(h2, 2))*(f_spinz1 + f_spinz2 - 3*(-2*f_spinz0 + f_spinz1)*pow(hsneg, 2) -(6*f_spinz0 - 4*f_spinz1 + f_spinz2)*pow(hsneg, 4))*(h1 + R1Val)*(h2 - R2Val)*cos(Phi - phi0m + phip0))/16.;
  value += -(sqrt(A0p)*sqrt(Ap0)*sqrt(1 - pow(h1, 2))*sqrt(1 - pow(h2, 2))*(-1 + pow(hsneg, 2))*(-f_spinz1 + f_spinz2 - (6*f_spinz0 - 4*f_spinz1 + f_spinz2)*pow(hsneg, 2))*(h1 + R1Val)*(h2 + R2Val)*cos(Phi + phi0p + 2*Phi1 - phip0))/16.;
  value += -(sqrt(Am0)*sqrt(Ap0)*(-1 + pow(h1, 2))*(-1 + pow(h2, 2))*(-1 + pow(hsneg, 2))*(-f_spinz1 + f_spinz2 - (6*f_spinz0 - 4*f_spinz1 + f_spinz2)*pow(hsneg, 2))*cos(2*Phi1 + phim0 - phip0))/16.;
  value += (sqrt(Ap0)*sqrt(Apm)*sqrt(1 - pow(h2, 2))*hsneg*sqrt(1 - pow(hsneg, 2))*(-6*f_spinz0 + 3*f_spinz2 + (6*f_spinz0 - 4*f_spinz1 + f_spinz2)*pow(hsneg, 2))*(1 + pow(h1, 2) + 2*h1*R1Val)*(h2 - R2Val)*cos(Phi + Phi1 + phip0 - phipm))/(32.*sqrt(2));
  value += -(sqrt(Amp)*sqrt(Ap0)*(6*f_spinz0 - 4*f_spinz1 + f_spinz2)*(-1 + pow(h1, 2))*sqrt(1 - pow(h2, 2))*hsneg*pow(1 - pow(hsneg, 2), 1.5)*(h2 + R2Val)*cos(Phi + 3*Phi1 + phimp - phip0))/(32.*sqrt(2));
  value += -(sqrt(A0m)*sqrt(A0p)*(-1 + pow(h1, 2))*(-1 + pow(h2, 2))*(-1 + pow(hsneg, 2))*(-f_spinz1 + f_spinz2 - (6*f_spinz0 - 4*f_spinz1 + f_spinz2)*pow(hsneg, 2))*cos(2*Phi - phi0m + phi0p + 2*Phi1))/16.;
  value += -(sqrt(A0m)*sqrt(Am0)*sqrt(1 - pow(h1, 2))*sqrt(1 - pow(h2, 2))*(-1 + pow(hsneg, 2))*(-f_spinz1 + f_spinz2 - (6*f_spinz0 - 4*f_spinz1 + f_spinz2)*pow(hsneg, 2))*(h1 - R1Val)*(h2 - R2Val)*cos(Phi - phi0m + 2*Phi1 + phim0))/16.;
  value += (sqrt(A0m)*sqrt(Apm)*sqrt(1 - pow(h1, 2))*hsneg*sqrt(1 - pow(hsneg, 2))*(-6*f_spinz0 + 3*f_spinz2 + (6*f_spinz0 - 4*f_spinz1 + f_spinz2)*pow(hsneg, 2))*(h1 + R1Val)*(-1 - pow(h2, 2) + 2*h2*R2Val)*cos(phi0m + Phi1 - phipm))/(32.*sqrt(2));
  value += (sqrt(A0m)*sqrt(Amp)*(6*f_spinz0 - 4*f_spinz1 + f_spinz2)*sqrt(1 - pow(h1, 2))*(-1 + pow(h2, 2))*hsneg*pow(1 - pow(hsneg, 2), 1.5)*(h1 - R1Val)*cos(2*Phi - phi0m + 3*Phi1 + phimp))/(32.*sqrt(2));
  value += -(sqrt(A0p)*sqrt(Am0)*sqrt(1 - pow(h1, 2))*sqrt(1 - pow(h2, 2))*(f_spinz1 + f_spinz2 - 3*(-2*f_spinz0 + f_spinz1)*pow(hsneg, 2) -(6*f_spinz0 - 4*f_spinz1 + f_spinz2)*pow(hsneg, 4))*(h1 - R1Val)*(h2 + R2Val)*cos(Phi + phi0p - phim0))/16.;
  value += (sqrt(A0p)*sqrt(Apm)*(6*f_spinz0 - 4*f_spinz1 + f_spinz2)*sqrt(1 - pow(h1, 2))*(-1 + pow(h2, 2))*hsneg*pow(1 - pow(hsneg, 2), 1.5)*(h1 + R1Val)*cos(2*Phi + phi0p + 3*Phi1 - phipm))/(32.*sqrt(2));
  value += (sqrt(A0p)*sqrt(Amp)*sqrt(1 - pow(h1, 2))*hsneg*sqrt(1 - pow(hsneg, 2))*(-6*f_spinz0 + 3*f_spinz2 + (6*f_spinz0 - 4*f_spinz1 + f_spinz2)*pow(hsneg, 2))*(-h1 + R1Val)*(1 + pow(h2, 2) + 2*h2*R2Val)*cos(phi0p - Phi1 - phimp))/(32.*sqrt(2));
  value += (sqrt(Am0)*sqrt(Apm)*(6*f_spinz0 - 4*f_spinz1 + f_spinz2)*(-1 + pow(h1, 2))*sqrt(1 - pow(h2, 2))*hsneg*pow(1 - pow(hsneg, 2), 1.5)*(-h2 + R2Val)*cos(Phi + 3*Phi1 + phim0 - phipm))/(32.*sqrt(2));
  value += -(sqrt(Am0)*sqrt(Amp)*sqrt(1 - pow(h2, 2))*hsneg*sqrt(1 - pow(hsneg, 2))*(-6*f_spinz0 + 3*f_spinz2 + (6*f_spinz0 - 4*f_spinz1 + f_spinz2)*pow(hsneg, 2))*(-1 - pow(h1, 2) + 2*h1*R1Val)*(h2 + R2Val)*cos(Phi + Phi1 - phim0 + phimp))/(32.*sqrt(2));
  value += (sqrt(Amp)*sqrt(Apm)*(6*f_spinz0 - 4*f_spinz1 + f_spinz2)*(-1 + pow(h1, 2))*(-1 + pow(h2, 2))*pow(-1 + pow(hsneg, 2), 2)*cos(2*Phi + 4*Phi1 + phimp - phipm))/128.;

  value *= term1Coeff*term2Coeff*betaVal;
  if (!(value==value)) cout << "Evaluate NaN=" << value << endl;
  return value;
}


Int_t RooSpinTwo_7DComplex_HVV::getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* /*rangeName*/) const{
  const int prime_h1=2;
  const int prime_h2=3;
  const int prime_hs=5;
  const int prime_Phi=7;
  const int prime_Phi1=11;

  if (matchArgs(allVars, analVars, RooArgSet(*hs.absArg(), *h1.absArg(), *h2.absArg(), *Phi.absArg(), *Phi1.absArg()))) return 6; // all integrated
  if (matchArgs(allVars, analVars, hs, h1, h2, Phi)) return 5; // No Phi1
  if (matchArgs(allVars, analVars, hs, h1, h2, Phi1)) return 4; // No Phi
  if (matchArgs(allVars, analVars, hs, h1, Phi, Phi1)) return 3; // No h2
  if (matchArgs(allVars, analVars, hs, h2, Phi, Phi1)) return 2; // No h1
  if (matchArgs(allVars, analVars, h1, h2, Phi, Phi1)) return 1; // No hs
  if (matchArgs(allVars, analVars, hs, Phi1)) return 7; // production angles integrated
  return 0;
}
Double_t RooSpinTwo_7DComplex_HVV::analyticalIntegral(Int_t code, const char* /*rangeName*/) const{
  bool isZZ = true;
  if (mV < 90.) isZZ = false;
  if (isZZ) {
    if ((m1+m2) > m12 || fabs(m2-mV)<fabs(m1-mV) || m2 <= 0.0 || m1 <= 0.0) return 1e-10;
  }
  else if ((m1+m2) > m12 || m2 <= 0.0 || m1 <= 0.0) return 1e-10;

  Double_t betaValSq = (1.-(pow(m1-m2, 2)/pow(m12, 2)))*(1.-(pow(m1+m2, 2)/pow(m12, 2)));
  if (betaValSq<0) return 1e-15;
  Double_t betaVal = sqrt(betaValSq);

  Double_t term1Coeff = pow(m1, 3)/(pow(pow(m1, 2)-pow(mV, 2), 2)+pow(mV*gamV, 2));
  Double_t term2Coeff = pow(m2, 3)/(pow(pow(m2, 2)-pow(mV, 2), 2)+pow(mV*gamV, 2));
  Double_t f_spinz0 = 1. - f_spinz1 - f_spinz2;

  Double_t Pi = TMath::Pi();

  Double_t hsneg = -hs;
  Double_t value = 0;

  Double_t
    A00Re, A00Im,
    AppRe, AppIm,
    A0pRe, A0pIm, Ap0Re, Ap0Im,
    AmmRe, AmmIm, A0mRe, A0mIm, Am0Re, Am0Im,
    ApmRe, ApmIm, AmpRe, AmpIm;

  calculateAmplitudes(
    A00Re, A00Im,
    AppRe, AppIm,
    A0pRe, A0pIm, Ap0Re, Ap0Im,
    AmmRe, AmmIm, A0mRe, A0mIm, Am0Re, Am0Im,
    ApmRe, ApmIm, AmpRe, AmpIm
    );

  Double_t A00 = A00Im*A00Im + A00Re*A00Re;
  Double_t App = AppIm*AppIm + AppRe*AppRe;
  Double_t Amm = AmmIm*AmmIm + AmmRe*AmmRe;
  Double_t Ap0 = Ap0Im*Ap0Im + Ap0Re*Ap0Re;
  Double_t A0p = A0pIm*A0pIm + A0pRe*A0pRe;
  Double_t Am0 = Am0Im*Am0Im + Am0Re*Am0Re;
  Double_t A0m = A0mIm*A0mIm + A0mRe*A0mRe;
  Double_t Apm=  ApmIm*ApmIm + ApmRe*ApmRe;
  Double_t Amp = AmpIm*AmpIm + AmpRe*AmpRe;

  Double_t phi00=atan2(A00Im, A00Re);
  Double_t phipp=atan2(AppIm, AppRe)-phi00;
  Double_t phimm=atan2(AmmIm, AmmRe)-phi00;
  Double_t phip0=atan2(Ap0Im, Ap0Re)-phi00;
  Double_t phi0p=atan2(A0pIm, A0pRe)-phi00;
  Double_t phim0=atan2(Am0Im, Am0Re)-phi00;
  Double_t phi0m=atan2(A0mIm, A0mRe)-phi00;
  Double_t phipm=atan2(ApmIm, ApmRe)-phi00;
  Double_t phimp=atan2(AmpIm, AmpRe)-phi00;

  switch (code)
  {
    // integrate all angles
  case 6:
  {
          value += (32*A00*(f_spinz0 + f_spinz1 + f_spinz2)*pow(Pi, 2))/45.;
          value += (32*Amm*(f_spinz0 + f_spinz1 + f_spinz2)*pow(Pi, 2))/45.;
          value += (32*App*(f_spinz0 + f_spinz1 + f_spinz2)*pow(Pi, 2))/45.;
          value += (32*Ap0*(f_spinz0 + f_spinz1 + f_spinz2)*pow(Pi, 2))/45.;
          value += (32*A0m*(f_spinz0 + f_spinz1 + f_spinz2)*pow(Pi, 2))/45.;
          value += (32*A0p*(f_spinz0 + f_spinz1 + f_spinz2)*pow(Pi, 2))/45.;
          value += (32*Am0*(f_spinz0 + f_spinz1 + f_spinz2)*pow(Pi, 2))/45.;
          value += (32*Apm*(f_spinz0 + f_spinz1 + f_spinz2)*pow(Pi, 2))/45.;
          value += (32*Amp*(f_spinz0 + f_spinz1 + f_spinz2)*pow(Pi, 2))/45.;
  }
    // projections onto Phi1, integrate all other angles
  case 5:
  {
          value += (16*A00*(f_spinz0 + f_spinz1 + f_spinz2)*Pi)/45.;
          value += (16*Amm*(f_spinz0 + f_spinz1 + f_spinz2)*Pi)/45.;
          value += (16*App*(f_spinz0 + f_spinz1 + f_spinz2)*Pi)/45.;
          value += (16*Ap0*(f_spinz0 + f_spinz1 + f_spinz2)*Pi)/45.;
          value += (16*A0m*(f_spinz0 + f_spinz1 + f_spinz2)*Pi)/45.;
          value += (16*A0p*(f_spinz0 + f_spinz1 + f_spinz2)*Pi)/45.;
          value += (16*Am0*(f_spinz0 + f_spinz1 + f_spinz2)*Pi)/45.;
          value += (16*Apm*(f_spinz0 + f_spinz1 + f_spinz2)*Pi)/45.;
          value += (16*Amp*(f_spinz0 + f_spinz1 + f_spinz2)*Pi)/45.;
          value += (-4*sqrt(2./3.)*sqrt(Amm)*sqrt(Apm)*(2*(f_spinz0 + f_spinz1) - 3*f_spinz2)*Pi*cos(2*Phi1 + phimm - phipm))/45.;
          value += (-4*sqrt(2./3.)*sqrt(Amp)*sqrt(App)*(2*(f_spinz0 + f_spinz1) - 3*f_spinz2)*Pi*cos(2*Phi1 + phimp - phipp))/45.;
          value += (-8*sqrt(Am0)*sqrt(Ap0)*(6*f_spinz0 + f_spinz1 - 4*f_spinz2)*Pi*cos(2*Phi1 + phim0 - phip0))/135.;
  }
    // projection to Phi, integrate all other angles
  case 4:
  {
          value += (16*A00*(f_spinz0 + f_spinz1 + f_spinz2)*Pi)/45.;
          value += (16*Amm*(f_spinz0 + f_spinz1 + f_spinz2)*Pi)/45.;
          value += (16*App*(f_spinz0 + f_spinz1 + f_spinz2)*Pi)/45.;
          value += (16*Ap0*(f_spinz0 + f_spinz1 + f_spinz2)*Pi)/45.;
          value += (16*A0m*(f_spinz0 + f_spinz1 + f_spinz2)*Pi)/45.;
          value += (16*A0p*(f_spinz0 + f_spinz1 + f_spinz2)*Pi)/45.;
          value += (16*Am0*(f_spinz0 + f_spinz1 + f_spinz2)*Pi)/45.;
          value += (16*Apm*(f_spinz0 + f_spinz1 + f_spinz2)*Pi)/45.;
          value += (16*Amp*(f_spinz0 + f_spinz1 + f_spinz2)*Pi)/45.;
          value += (sqrt(A00)*sqrt(Amm)*(f_spinz0 + f_spinz1 + f_spinz2)*pow(Pi, 3)*R1Val*R2Val*cos(Phi - phimm))/20.;
          value += (sqrt(A00)*sqrt(App)*(f_spinz0 + f_spinz1 + f_spinz2)*pow(Pi, 3)*R1Val*R2Val*cos(Phi + phipp))/20.;
          value += (8*sqrt(Amm)*sqrt(App)*(f_spinz0 + f_spinz1 + f_spinz2)*Pi*cos(2*Phi - phimm + phipp))/45.;
          value += (sqrt(A0m)*sqrt(Ap0)*(f_spinz0 + f_spinz1 + f_spinz2)*pow(Pi, 3)*R1Val*R2Val*cos(Phi - phi0m + phip0))/20.;
          value += (sqrt(A0p)*sqrt(Am0)*(f_spinz0 + f_spinz1 + f_spinz2)*pow(Pi, 3)*R1Val*R2Val*cos(Phi + phi0p - phim0))/20.;
  }
    // projections to h2, integrate over all others
  case 3:
  {
          value += (-8*A00*(f_spinz0 + f_spinz1 + f_spinz2)*(-1 + pow(h2, 2))*pow(Pi, 2))/15.;
          value += (-4*Amm*(f_spinz0 + f_spinz1 + f_spinz2)*pow(Pi, 2)*(-1 - pow(h2, 2) + 2*h2*R2Val))/15.;
          value += (4*App*(f_spinz0 + f_spinz1 + f_spinz2)*pow(Pi, 2)*(1 + pow(h2, 2) + 2*h2*R2Val))/15.;
          value += (-8*Ap0*(f_spinz0 + f_spinz1 + f_spinz2)*(-1 + pow(h2, 2))*pow(Pi, 2))/15.;
          value += (-4*A0m*(f_spinz0 + f_spinz1 + f_spinz2)*pow(Pi, 2)*(-1 - pow(h2, 2) + 2*h2*R2Val))/15.;
          value += (4*A0p*(f_spinz0 + f_spinz1 + f_spinz2)*pow(Pi, 2)*(1 + pow(h2, 2) + 2*h2*R2Val))/15.;
          value += (-8*Am0*(f_spinz0 + f_spinz1 + f_spinz2)*(-1 + pow(h2, 2))*pow(Pi, 2))/15.;
          value += (-4*Apm*(f_spinz0 + f_spinz1 + f_spinz2)*pow(Pi, 2)*(-1 - pow(h2, 2) + 2*h2*R2Val))/15.;
          value += (4*Amp*(f_spinz0 + f_spinz1 + f_spinz2)*pow(Pi, 2)*(1 + pow(h2, 2) + 2*h2*R2Val))/15.;
  }
    // projections to h1, integrate all others
  case 2:
  {
          value += (-8*A00*(f_spinz0 + f_spinz1 + f_spinz2)*(-1 + pow(h1, 2))*pow(Pi, 2))/15.;
          value += (-4*Amm*(f_spinz0 + f_spinz1 + f_spinz2)*pow(Pi, 2)*(-1 - pow(h1, 2) + 2*h1*R1Val))/15.;
          value += (4*App*(f_spinz0 + f_spinz1 + f_spinz2)*pow(Pi, 2)*(1 + pow(h1, 2) + 2*h1*R1Val))/15.;
          value += (4*Ap0*(f_spinz0 + f_spinz1 + f_spinz2)*pow(Pi, 2)*(1 + pow(h1, 2) + 2*h1*R1Val))/15.;
          value += (-8*A0m*(f_spinz0 + f_spinz1 + f_spinz2)*(-1 + pow(h1, 2))*pow(Pi, 2))/15.;
          value += (-8*A0p*(f_spinz0 + f_spinz1 + f_spinz2)*(-1 + pow(h1, 2))*pow(Pi, 2))/15.;
          value += (-4*Am0*(f_spinz0 + f_spinz1 + f_spinz2)*pow(Pi, 2)*(-1 - pow(h1, 2) + 2*h1*R1Val))/15.;
          value += (4*Apm*(f_spinz0 + f_spinz1 + f_spinz2)*pow(Pi, 2)*(1 + pow(h1, 2) + 2*h1*R1Val))/15.;
          value += (-4*Amp*(f_spinz0 + f_spinz1 + f_spinz2)*pow(Pi, 2)*(-1 - pow(h1, 2) + 2*h1*R1Val))/15.;
  }
    // projections to hs, integrate all others
  case 1:
  {
          value += (2*A00*(2*f_spinz0 + 3*f_spinz2 - 6*(2*f_spinz0 - 2*f_spinz1 + f_spinz2)*pow(hsneg, 2) +3*(6*f_spinz0 - 4*f_spinz1 + f_spinz2)*pow(hsneg, 4))*pow(Pi, 2))/9.;
          value += (2*Amm*(2*f_spinz0 + 3*f_spinz2 - 6*(2*f_spinz0 - 2*f_spinz1 + f_spinz2)*pow(hsneg, 2) +3*(6*f_spinz0 - 4*f_spinz1 + f_spinz2)*pow(hsneg, 4))*pow(Pi, 2))/9.;
          value += (2*App*(2*f_spinz0 + 3*f_spinz2 - 6*(2*f_spinz0 - 2*f_spinz1 + f_spinz2)*pow(hsneg, 2) +3*(6*f_spinz0 - 4*f_spinz1 + f_spinz2)*pow(hsneg, 4))*pow(Pi, 2))/9.;
          value += (4*Ap0*(f_spinz1 + f_spinz2 - 3*(-2*f_spinz0 + f_spinz1)*pow(hsneg, 2) -(6*f_spinz0 - 4*f_spinz1 + f_spinz2)*pow(hsneg, 4))*pow(Pi, 2))/9.;
          value += (4*A0m*(f_spinz1 + f_spinz2 - 3*(-2*f_spinz0 + f_spinz1)*pow(hsneg, 2) -(6*f_spinz0 - 4*f_spinz1 + f_spinz2)*pow(hsneg, 4))*pow(Pi, 2))/9.;
          value += (4*A0p*(f_spinz1 + f_spinz2 - 3*(-2*f_spinz0 + f_spinz1)*pow(hsneg, 2) -(6*f_spinz0 - 4*f_spinz1 + f_spinz2)*pow(hsneg, 4))*pow(Pi, 2))/9.;
          value += (4*Am0*(f_spinz1 + f_spinz2 - 3*(-2*f_spinz0 + f_spinz1)*pow(hsneg, 2) -(6*f_spinz0 - 4*f_spinz1 + f_spinz2)*pow(hsneg, 4))*pow(Pi, 2))/9.;
          value += (Apm*(4*f_spinz1 + f_spinz2 + 6*f_spinz2*pow(hsneg, 2) + (-4*f_spinz1 + f_spinz2)*pow(hsneg, 4) +6*f_spinz0*pow(-1 + pow(hsneg, 2), 2))*pow(Pi, 2))/9.;
          value += (Amp*(4*f_spinz1 + f_spinz2 + 6*f_spinz2*pow(hsneg, 2) + (-4*f_spinz1 + f_spinz2)*pow(hsneg, 4) +6*f_spinz0*pow(-1 + pow(hsneg, 2), 2))*pow(Pi, 2))/9.;
  }

    // production angles integrated out
  case 7:
  {
          value += (A00*(f_spinz0 + f_spinz1 + f_spinz2)*(-1 + pow(h1, 2))*(-1 + pow(h2, 2))*Pi)/5.;
          value += (Amm*(f_spinz0 + f_spinz1 + f_spinz2)*Pi*(-1 - pow(h1, 2) + 2*h1*R1Val)*
            (-1 - pow(h2, 2) + 2*h2*R2Val))/20.;
          value += (App*(f_spinz0 + f_spinz1 + f_spinz2)*Pi*(1 + pow(h1, 2) + 2*h1*R1Val)*
            (1 + pow(h2, 2) + 2*h2*R2Val))/20.;
          value += -(Ap0*(f_spinz0 + f_spinz1 + f_spinz2)*(-1 + pow(h2, 2))*Pi*
            (1 + pow(h1, 2) + 2*h1*R1Val))/10.;
          value += (A0m*(f_spinz0 + f_spinz1 + f_spinz2)*(-1 + pow(h1, 2))*Pi*
            (-1 - pow(h2, 2) + 2*h2*R2Val))/10.;
          value += -(A0p*(f_spinz0 + f_spinz1 + f_spinz2)*(-1 + pow(h1, 2))*Pi*
            (1 + pow(h2, 2) + 2*h2*R2Val))/10.;
          value += (Am0*(f_spinz0 + f_spinz1 + f_spinz2)*(-1 + pow(h2, 2))*Pi*
            (-1 - pow(h1, 2) + 2*h1*R1Val))/10.;
          value += -(Apm*(f_spinz0 + f_spinz1 + f_spinz2)*Pi*(1 + pow(h1, 2) + 2*h1*R1Val)*
            (-1 - pow(h2, 2) + 2*h2*R2Val))/20.;
          value += -(Amp*(f_spinz0 + f_spinz1 + f_spinz2)*Pi*(-1 - pow(h1, 2) + 2*h1*R1Val)*
            (1 + pow(h2, 2) + 2*h2*R2Val))/20.;
          value += (sqrt(A00)*sqrt(Amm)*(f_spinz0 + f_spinz1 + f_spinz2)*sqrt(1 - pow(h1, 2))*
            sqrt(1 - pow(h2, 2))*Pi*(h1 - R1Val)*(h2 - R2Val)*cos(Phi - phimm))/5.;
          value += (sqrt(A00)*sqrt(App)*(f_spinz0 + f_spinz1 + f_spinz2)*sqrt(1 - pow(h1, 2))*
            sqrt(1 - pow(h2, 2))*Pi*(h1 + R1Val)*(h2 + R2Val)*cos(Phi + phipp))/5.;

          value += (sqrt(Amm)*sqrt(App)*(f_spinz0 + f_spinz1 + f_spinz2)*(-1 + pow(h1, 2))*
            (-1 + pow(h2, 2))*Pi*cos(2*Phi - phimm + phipp))/10.;

          value += (sqrt(A0m)*sqrt(Ap0)*(f_spinz0 + f_spinz1 + f_spinz2)*sqrt(1 - pow(h1, 2))*
            sqrt(1 - pow(h2, 2))*Pi*(h1 + R1Val)*(-h2 + R2Val)*
            cos(Phi - phi0m + phip0))/5.;

          value += (sqrt(A0p)*sqrt(Am0)*(f_spinz0 + f_spinz1 + f_spinz2)*sqrt(1 - pow(h1, 2))*
            sqrt(1 - pow(h2, 2))*Pi*(-h1 + R1Val)*(h2 + R2Val)*
            cos(Phi + phi0p - phim0))/5.;
  }
  }

  value *= betaVal*term1Coeff*term2Coeff;

  if (!(value==value)){
    cout << "Integral NaN=" << value << " at "
      << "h1=" << h1 << '\t'
      << "h2=" << h2 << '\t'
      << "hs=" << hs << '\t'
      << "Phi1=" << Phi1 << '\t'
      << "Phi=" << Phi << '\t'
      << "m1=" << m1 << '\t'
      << "m2=" << m2 << '\t'
      << "m12=" << m12 << '\t'
      << endl;
    cout << "Possible sources:\n"
      << "betaVal=" << betaVal << '\t'
      << "term1Coeff=" << term1Coeff << '\t'
      << "term2Coeff=" << term2Coeff << '\t'
      << "A00=" << A00 << '\t'
      << "App=" << App << '\t'
      << "Amm=" << Amm << '\t'
      << "phi00=" << phi00 << '\t'
      << "phipp=" << phipp << '\t'
      << "phimm=" << phimm << '\t'
      << endl;
  }
  return value;
}

