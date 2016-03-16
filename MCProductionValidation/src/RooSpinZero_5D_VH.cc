#include "../include/RooSpinZero_5D_VH.h" 

RooSpinZero_5D_VH::RooSpinZero_5D_VH(
  const char *name, const char *title,
  modelMeasurables _measurables,
  modelParameters _parameters,
  int _Vdecay1, int _Vdecay2
  ) : RooSpinZero(
  name, title,
  _measurables,
  _parameters,
  _Vdecay1, _Vdecay2
  )
{}


RooSpinZero_5D_VH::RooSpinZero_5D_VH(
  const RooSpinZero_5D_VH& other, const char* name
  ) : RooSpinZero(other, name)
{}


Double_t RooSpinZero_5D_VH::evaluate() const{
  if (m1 <= 0.0 || (m2 <= 0.0 && Vdecay2!=0) || Vdecay1==0) return 1e-15;
  // No need to set m1_ or m2_

  Double_t A00Re, A00Im, AppRe, AppIm, AmmRe, AmmIm;
  calculateAmplitudes(A00Re, A00Im, AppRe, AppIm, AmmRe, AmmIm);

  Double_t A00 = A00Im*A00Im + A00Re*A00Re;
  Double_t App = AppIm*AppIm + AppRe*AppRe;
  Double_t Amm = AmmIm*AmmIm + AmmRe*AmmRe;

  Double_t phi00=atan2(A00Im, A00Re);
  Double_t phipp=atan2(AppIm, AppRe)-phi00;
  Double_t phimm=atan2(AmmIm, AmmRe)-phi00;

  Double_t value = 0;
  value += (A00*(-1 + pow(h1, 2))*(-1 + pow(h2, 2)))/4.;
  value += (App*(1 + pow(h1, 2) - 2*h1*R1Val)*(1 + pow(h2, 2) + 2*h2*R2Val))/16.;
  value += (Amm*(1 + pow(h1, 2) + 2*h1*R1Val)*(1 + pow(h2, 2) - 2*h2*R2Val))/16.;
  value += -(sqrt(A00)*sqrt(App)*sqrt(1 - pow(h1, 2))*sqrt(1 - pow(h2, 2))*(h1 - R1Val)*(h2 + R2Val)*cos(Phi + phipp))/4.;
  value += -(sqrt(A00)*sqrt(Amm)*sqrt(1 - pow(h1, 2))*sqrt(1 - pow(h2, 2))*(h1 + R1Val)*(h2 - R2Val)*cos(Phi - phimm))/4.;
  value += (sqrt(Amm)*sqrt(App)*(-1 + pow(h1, 2))*(-1 + pow(h2, 2))*cos(2*Phi - phimm + phipp))/8.;
  return value;
}

Int_t RooSpinZero_5D_VH::getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* /*rangeName*/) const{
  if (matchArgs(allVars, analVars, RooArgSet(*h1.absArg(), *h2.absArg(), *hs.absArg(), *Phi1.absArg(), *Phi.absArg()))) return 6;
  if (matchArgs(allVars, analVars, h1, h2, Phi1, Phi)) return 1;
  if (matchArgs(allVars, analVars, h1, h2, hs, Phi1)) return 5;
  if (matchArgs(allVars, analVars, h1, hs, Phi1, Phi)) return 3;
  if (matchArgs(allVars, analVars, h2, hs, Phi1, Phi)) return 4;
  if (matchArgs(allVars, analVars, h1, h2, hs, Phi)) return 2;
  return 0;
}

Double_t RooSpinZero_5D_VH::analyticalIntegral(Int_t code, const char* /*rangeName*/) const{
  if (m1 <= 0.0 || (m2 <= 0.0 && Vdecay2!=0) || Vdecay1==0) return 1e-10;
  // No need to set m1_ or m2_

  Double_t A00Re, A00Im, AppRe, AppIm, AmmRe, AmmIm;
  calculateAmplitudes(A00Re, A00Im, AppRe, AppIm, AmmRe, AmmIm);

  Double_t A00 = A00Im*A00Im + A00Re*A00Re;
  Double_t App = AppIm*AppIm + AppRe*AppRe;
  Double_t Amm = AmmIm*AmmIm + AmmRe*AmmRe;

  Double_t phi00=atan2(A00Im, A00Re);
  Double_t phipp=atan2(AppIm, AppRe)-phi00;
  Double_t phimm=atan2(AmmIm, AmmRe)-phi00;

  Double_t value = 0;
  switch (code)
  {
    // projections to hs
  case 1:
  {
          value += (16*A00*pow(TMath::Pi(), 2))/9.;
          value += (16*App*pow(TMath::Pi(), 2))/9.;
          value += (16*Amm*pow(TMath::Pi(), 2))/9.;
          return value;
  }
    // projections to Phi1
  case 2:
  {
          value += (16*A00*TMath::Pi())/9.;
          value += (16*App*TMath::Pi())/9.;
          value += (16*Amm*TMath::Pi())/9.;
          return value;
  }
    // projections to h2
  case 3:
  {
          value += (-8*A00*(-1 + pow(h2, 2))*pow(TMath::Pi(), 2))/3.;
          value += (4*App*pow(TMath::Pi(), 2)*(1 + pow(h2, 2) + 2*h2*R2Val))/3.;
          value += (4*Amm*pow(TMath::Pi(), 2)*(1 + pow(h2, 2) - 2*h2*R2Val))/3.;
          return value;
  }
    // projections to h1
  case 4:
  {
          value += (-8*A00*(-1 + pow(h1, 2))*pow(TMath::Pi(), 2))/3.;
          value += (4*App*pow(TMath::Pi(), 2)*(1 + pow(h1, 2) - 2*h1*R1Val))/3.;
          value += (4*Amm*pow(TMath::Pi(), 2)*(1 + pow(h1, 2) + 2*h1*R1Val))/3.;
  }
    // projections to Phi
  case 5:
  {
          value += (16*A00*TMath::Pi())/9.;
          value += (16*App*TMath::Pi())/9.;
          value += (16*Amm*TMath::Pi())/9.;
          value += (sqrt(A00)*sqrt(App)*pow(TMath::Pi(), 3)*R1Val*R2Val*cos(Phi + phipp))/4.;
          value += (sqrt(A00)*sqrt(Amm)*pow(TMath::Pi(), 3)*R1Val*R2Val*cos(Phi - phimm))/4.;
          value += (8*sqrt(Amm)*sqrt(App)*TMath::Pi()*cos(2*Phi - phimm + phipp))/9.;
          return value;
  }
    // projected everything
  case 6:
  {
          value += (32*A00*pow(TMath::Pi(), 2))/9.;
          value += (32*App*pow(TMath::Pi(), 2))/9.;
          value += (32*Amm*pow(TMath::Pi(), 2))/9.;
          return value;
  }
  }
  assert(0);
  return 0;
}




