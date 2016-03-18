#include "../include/RooSpinTwo.h"

RooSpinTwo::RooSpinTwo(
  const char* name, const char* title,
  modelMeasurables _measurables,
  modelParameters _parameters,
  int _Vdecay1, int _Vdecay2
  ) : RooAbsPdf(name, title),
  Vdecay1((Int_t)_Vdecay1), Vdecay2((Int_t)_Vdecay2),
  h1("h1", "h1", this),
  h2("h2", "h2", this),
  Phi("Phi", "Phi", this),
  m1("m1", "m1", this),
  m2("m2", "m2", this),
  m12("m12", "m12", this),
  hs("hs", "hs", this),
  Phi1("Phi1", "Phi1", this),
  //Y("Y", "Y", this),

  mX("mX", "mX", this, (RooAbsReal&)*(_parameters.mX)),
  gamX("gamX", "gamX", this, (RooAbsReal&)*(_parameters.gamX)),
  mV("mV", "mV", this, (RooAbsReal&)*(_parameters.mV)),
  gamV("gamV", "gamV", this, (RooAbsReal&)*(_parameters.gamV)),
  R1Val("R1Val", "R1Val", this, (RooAbsReal&)*(_parameters.R1Val)),
  R2Val("R2Val", "R2Val", this, (RooAbsReal&)*(_parameters.R2Val)),

  b1Val("b1Val", "b1Val", this, (RooAbsReal&)*(_parameters.bList[0][0])),
  b2Val("b2Val", "b2Val", this, (RooAbsReal&)*(_parameters.bList[1][0])),
  b3Val("b3Val", "b3Val", this, (RooAbsReal&)*(_parameters.bList[2][0])),
  b4Val("b4Val", "b4Val", this, (RooAbsReal&)*(_parameters.bList[3][0])),
  b5Val("b5Val", "b5Val", this, (RooAbsReal&)*(_parameters.bList[4][0])),
  b6Val("b6Val", "b6Val", this, (RooAbsReal&)*(_parameters.bList[5][0])),
  b7Val("b7Val", "b7Val", this, (RooAbsReal&)*(_parameters.bList[6][0])),
  b8Val("b8Val", "b8Val", this, (RooAbsReal&)*(_parameters.bList[7][0])),
  b9Val("b9Val", "b9Val", this, (RooAbsReal&)*(_parameters.bList[8][0])),
  b10Val("b10Val", "b10Val", this, (RooAbsReal&)*(_parameters.bList[9][0])),

  b1ValIm("b1ValIm", "b1ValIm", this, (RooAbsReal&)*(_parameters.bList[0][1])),
  b2ValIm("b2ValIm", "b2ValIm", this, (RooAbsReal&)*(_parameters.bList[1][1])),
  b3ValIm("b3ValIm", "b3ValIm", this, (RooAbsReal&)*(_parameters.bList[2][1])),
  b4ValIm("b4ValIm", "b4ValIm", this, (RooAbsReal&)*(_parameters.bList[3][1])),
  b5ValIm("b5ValIm", "b5ValIm", this, (RooAbsReal&)*(_parameters.bList[4][1])),
  b6ValIm("b6ValIm", "b6ValIm", this, (RooAbsReal&)*(_parameters.bList[5][1])),
  b7ValIm("b7ValIm", "b7ValIm", this, (RooAbsReal&)*(_parameters.bList[6][1])),
  b8ValIm("b8ValIm", "b8ValIm", this, (RooAbsReal&)*(_parameters.bList[7][1])),
  b9ValIm("b9ValIm", "b9ValIm", this, (RooAbsReal&)*(_parameters.bList[8][1])),
  b10ValIm("b10ValIm", "b10ValIm", this, (RooAbsReal&)*(_parameters.bList[9][1])),

  Lambda("Lambda", "Lambda", this, (RooAbsReal&)*(_parameters.Lambda)),

  f_spinz1("f_spinz1", "f_spinz1", this, (RooAbsReal&)*(_parameters.f_spinz1)),
  f_spinz2("f_spinz2", "f_spinz2", this, (RooAbsReal&)*(_parameters.f_spinz2))
{
  setProxies(_measurables);
}


RooSpinTwo::RooSpinTwo(const RooSpinTwo& other, const char* name) :
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
//Y("Y", this, other.Y),

mX("mX", this, other.mX),
gamX("gamX", this, other.gamX),
mV("mV", this, other.mV),
gamV("gamV", this, other.gamV),
R1Val("R1Val", this, other.R1Val),
R2Val("R2Val", this, other.R2Val),

b1Val("b1Val", this, other.b1Val),
b2Val("a2Val", this, other.b2Val),
b3Val("b3Val", this, other.b3Val),
b4Val("b4Val", this, other.b4Val),
b5Val("b5Val", this, other.b5Val),
b6Val("b6Val", this, other.b6Val),
b7Val("b7Val", this, other.b7Val),
b8Val("b8Val", this, other.b8Val),
b9Val("b9Val", this, other.b9Val),
b10Val("b10Val", this, other.b10Val),

b1ValIm("b1ValIm", this, other.b1ValIm),
b2ValIm("a2ValIm", this, other.b2ValIm),
b3ValIm("b3ValIm", this, other.b3ValIm),
b4ValIm("b4ValIm", this, other.b4ValIm),
b5ValIm("b5ValIm", this, other.b5ValIm),
b6ValIm("b6ValIm", this, other.b6ValIm),
b7ValIm("b7ValIm", this, other.b7ValIm),
b8ValIm("b8ValIm", this, other.b8ValIm),
b9ValIm("b9ValIm", this, other.b9ValIm),
b10ValIm("b10ValIm", this, other.b10ValIm),

Lambda("Lambda", this, other.Lambda),

f_spinz1("f_spinz1", this, other.f_spinz1),
f_spinz2("f_spinz2", this, other.f_spinz2)
{}

void RooSpinTwo::calculateCi(std::vector<Double_t>& ciRe, std::vector<Double_t>& ciIm, bool isGammaV1, bool isGammaV2) const{
  Double_t m1_=m1; if (Vdecay1==0) m1_=0;
  Double_t m2_=m2; if (Vdecay2==0) m2_=0;
  Double_t m1sq = pow(m1_, 2);
  Double_t m2sq = pow(m2_, 2);
  Double_t mVsq = pow(mV, 2);
  Double_t m12sq = pow(m12, 2);

  Double_t s = (m12sq-m1sq-m2sq)/2.;
  if (m1sq>m12sq || m2sq>m12sq) s = -s;
  Double_t kappa = s/pow(Lambda, 2);

  if (!isGammaV1 && !isGammaV2 && !(Vdecay1==0 || Vdecay2==0)){ // ZZ/WW
    Double_t c1Re = 2.*(b1Val + b2Val*kappa*(1.+m1sq/s)*(1.+m2sq/s) + b5Val*mVsq/s); ciRe.push_back(c1Re);
    Double_t c2Re = -0.5*b1Val + b3Val*kappa*(1.-(m1sq+m2sq)/(2*s)) + 2.*b4Val*kappa + b7Val*kappa*mVsq/s; ciRe.push_back(c2Re);
    Double_t c3Re = -(b2Val/2.+b3Val+2.*b4Val)*kappa*m12sq/s; ciRe.push_back(c3Re);
    Double_t c41Re = -b1Val - b2Val*kappa - (b2Val*m1sq+b3Val*m2sq+2.*b6Val*mVsq)*kappa/s; ciRe.push_back(c41Re);
    Double_t c42Re = -b1Val - b2Val*kappa - (b2Val*m2sq+b3Val*m1sq+2.*b6Val*mVsq)*kappa/s; ciRe.push_back(c42Re);
    Double_t c5Re = 2.*b8Val*kappa*(m12sq)/s; ciRe.push_back(c5Re);
    Double_t c6Re = b9Val*kappa*mVsq/s; ciRe.push_back(c6Re);
    Double_t c7Re = b10Val*m12sq*mVsq*pow(kappa/s, 2); ciRe.push_back(c7Re);

    Double_t c1Im = 2.*(b1ValIm + b2ValIm*kappa*(1.+m1sq/s)*(1.+m2sq/s) + b5ValIm*mVsq/s); ciIm.push_back(c1Im);
    Double_t c2Im = -0.5*b1ValIm + b3ValIm*kappa*(1.-(m1sq+m2sq)/(2*s)) + 2.*b4ValIm*kappa + b7ValIm*kappa*mVsq/s; ciIm.push_back(c2Im);
    Double_t c3Im = -(b2ValIm/2.+b3ValIm+2.*b4ValIm)*kappa*m12sq/s; ciIm.push_back(c3Im);
    Double_t c41Im = -b1ValIm - b2ValIm*kappa - (b2ValIm*m1sq+b3ValIm*m2sq+2.*b6ValIm*mVsq)*kappa/s; ciIm.push_back(c41Im);
    Double_t c42Im = -b1ValIm - b2ValIm*kappa - (b2ValIm*m2sq+b3ValIm*m1sq+2.*b6ValIm*mVsq)*kappa/s; ciIm.push_back(c42Im);
    Double_t c5Im = 2.*b8ValIm*kappa*(m12sq)/s; ciIm.push_back(c5Im);
    Double_t c6Im = b9ValIm*kappa*mVsq/s; ciIm.push_back(c6Im);
    Double_t c7Im = b10ValIm*m12sq*mVsq*pow(kappa/s, 2); ciIm.push_back(c7Im);
  }
  //else if ((!isGammaV1 || !isGammaV2) && !(Vdecay1==0 && Vdecay2==0)){ // ZGs/ZG
  //???
  //}
  //else{ // GG/GGs
  else{ // Z/G/Gs - G/Gs
    Double_t c1Re = 2.*(b1Val + b2Val*kappa*(1.+m1sq/s)*(1.+m2sq/s)); ciRe.push_back(c1Re);
    Double_t c2Re = -0.5*b1Val + b3Val*kappa*(1.-(m1sq+m2sq)/(2*s)) + 2.*b4Val*kappa; ciRe.push_back(c2Re);
    Double_t c3Re = -(b2Val/2.+b3Val+2.*b4Val)*kappa*m12sq/s; ciRe.push_back(c3Re);
    Double_t c41Re = -b1Val - b2Val*kappa - (b2Val*m1sq+b3Val*m2sq)*kappa/s; ciRe.push_back(c41Re);
    Double_t c42Re = -b1Val - b2Val*kappa - (b2Val*m2sq+b3Val*m1sq)*kappa/s; ciRe.push_back(c42Re);
    Double_t c5Re = 2.*b8Val*kappa*(m12sq)/s; ciRe.push_back(c5Re);
    Double_t c6Re = 0; ciRe.push_back(c6Re);
    Double_t c7Re = 0; ciRe.push_back(c7Re);

    Double_t c1Im = 2.*(b1ValIm + b2ValIm*kappa*(1.+m1sq/s)*(1.+m2sq/s)); ciIm.push_back(c1Im);
    Double_t c2Im = -0.5*b1ValIm + b3ValIm*kappa*(1.-(m1sq+m2sq)/(2*s)) + 2.*b4ValIm*kappa; ciIm.push_back(c2Im);
    Double_t c3Im = -(b2ValIm/2.+b3ValIm+2.*b4ValIm)*kappa*m12sq/s; ciIm.push_back(c3Im);
    Double_t c41Im = -b1ValIm - b2ValIm*kappa - (b2ValIm*m1sq+b3ValIm*m2sq)*kappa/s; ciIm.push_back(c41Im);
    Double_t c42Im = -b1ValIm - b2ValIm*kappa - (b2ValIm*m2sq+b3ValIm*m1sq)*kappa/s; ciIm.push_back(c42Im);
    Double_t c5Im = 2.*b8ValIm*kappa*(m12sq)/s; ciIm.push_back(c5Im);
    Double_t c6Im = 0; ciIm.push_back(c6Im);
    Double_t c7Im = 0; ciIm.push_back(c7Im);
  }
}
void RooSpinTwo::calculatePropagator(Double_t& propRe, Double_t& propIm, Double_t mass, bool useGamma)const{
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
void RooSpinTwo::calculateAmplitudeScale(bool isGammaV1, bool isGammaV2)const{

}
void RooSpinTwo::calculateAmplitudes(
  Double_t& A00Re, Double_t& A00Im,
  Double_t& AppRe, Double_t& AppIm, Double_t& A0pRe, Double_t& A0pIm, Double_t& Ap0Re, Double_t& Ap0Im,
  Double_t& AmmRe, Double_t& AmmIm, Double_t& A0mRe, Double_t& A0mIm, Double_t& Am0Re, Double_t& Am0Im,
  Double_t& ApmRe, Double_t& ApmIm, Double_t& AmpRe, Double_t& AmpIm,
  bool isGammaV1, bool isGammaV2
  )const{
  Double_t m1_=m1; if (Vdecay1==0) m1_=0;
  Double_t m2_=m2; if (Vdecay2==0) m2_=0;

  std::vector<Double_t> ciRe;
  std::vector<Double_t> ciIm;
  calculateCi(ciRe, ciIm, isGammaV1, isGammaV2);

  Double_t propV1Re=0, propV2Re=0;
  Double_t propV1Im=-1, propV2Im=-1;
  if (Vdecay1!=0) calculatePropagator(propV1Re, propV1Im, m1_, isGammaV1);
  if (Vdecay2!=0) calculatePropagator(propV2Re, propV2Im, m2_, isGammaV2);

  Double_t c1Re = ciRe.at(0);
  Double_t c2Re = ciRe.at(1);
  Double_t c3Re = ciRe.at(2);
  Double_t c41Re = ciRe.at(3);
  Double_t c42Re = ciRe.at(4);
  Double_t c5Re = ciRe.at(5);
  Double_t c6Re = ciRe.at(6);
  Double_t c7Re = ciRe.at(7);
  Double_t c1Im = ciIm.at(0);
  Double_t c2Im = ciIm.at(1);
  Double_t c3Im = ciIm.at(2);
  Double_t c41Im = ciIm.at(3);
  Double_t c42Im = ciIm.at(4);
  Double_t c5Im = ciIm.at(5);
  Double_t c6Im = ciIm.at(6);
  Double_t c7Im = ciIm.at(7);

  Double_t eta1 = m1_ / m12;
  Double_t eta2 = m2_ / m12;
  Double_t eta1p2 = eta1*eta2;

  Double_t eta1sq = eta1*eta1;
  Double_t eta2sq = eta2*eta2;
  Double_t eta1p2sq = pow(eta1p2, 2);

  Double_t etas = (1. - eta1sq - eta2sq)/2.;
  if (pow(eta1+eta2, 2)>1.) etas = -etas;
  Double_t x = etas;
  Double_t xsq = x*x;
  Double_t xxp = (pow(etas, 2)-eta1p2sq);
  if (xxp<0) xxp=0;

  Double_t A00Re_tmp=0, A00Im_tmp=0,
    AppRe_tmp=0, AppIm_tmp=0, A0pRe_tmp=0, A0pIm_tmp=0, Ap0Re_tmp=0, Ap0Im_tmp=0,
    AmmRe_tmp=0, AmmIm_tmp=0, A0mRe_tmp=0, A0mIm_tmp=0, Am0Re_tmp=0, Am0Im_tmp=0,
    ApmRe_tmp=0, ApmIm_tmp=0, AmpRe_tmp=0, AmpIm_tmp=0;

  A00Re_tmp =
    pow(m12, 4)*sqrt(2./3.)*
    (
    c1Re*(
    eta1p2sq * (xsq - eta1p2sq/4.)
    - (pow(eta1, 4)+pow(eta2, 4)) * xsq/2.
    + (pow(eta1, 8)+pow(eta2, 8))/8.
    + xsq/2.
    - (pow(eta1, 4)+pow(eta2, 4))/4.
    + 1.0/8.
    )
    + c2Re*2.*xxp*(
    (pow(eta1, 4) + pow(eta2, 4))
    - 2.*(eta1p2sq + 2.*xxp)
    - 1.
    )
    -c3Re*8.*pow(xxp, 2)
    + c41Re*2.*xxp*(1.+eta1sq-eta2sq)
    + c42Re*2.*xxp*(1.-eta1sq+eta2sq)
    );

  A00Im_tmp =
    pow(m12, 4)*sqrt(2./3.)*
    (
    c1Im*(
    eta1p2sq * (xsq - eta1p2sq/4.)
    - (pow(eta1, 4)+pow(eta2, 4)) * xsq/2.
    + (pow(eta1, 8)+pow(eta2, 8))/8.
    + xsq/2.
    - (pow(eta1, 4)+pow(eta2, 4))/4.
    + 1.0/8.
    )
    + c2Im*2.*xxp*(
    (pow(eta1, 4) + pow(eta2, 4))
    - 2.*(eta1p2sq + 2.*xxp)
    - 1.
    )
    -c3Im*8.*pow(xxp, 2)
    + c41Im*2.*xxp*(1.+eta1sq-eta2sq)
    + c42Im*2.*xxp*(1.-eta1sq+eta2sq)
    );

  //-----------------------------------------------------------------------
  // No m1_/m2_ singularities in A--

  AmmRe_tmp =
    pow(m12, 4)/sqrt(6.)*
    (
    c1Re/4.*(
    1. + 4.*xxp - pow(eta1sq-eta2sq, 2)
    )
    + c2Re*8.*xxp
    + c5Im*8.*pow(xxp, 1.5)
    - c6Im*4.*sqrt(xxp)
    );

  AmmIm_tmp =
    pow(m12, 4)/sqrt(6.)*
    (
    c1Im/4.*(
    1. + 4.*xxp - pow(eta1sq-eta2sq, 2)
    )
    + c2Im*8.*xxp
    - c5Re*8.*pow(xxp, 1.5)
    + c6Re*4.*sqrt(xxp)
    );

  if (m1_!=0) { AmmIm_tmp *= eta1; AmmRe_tmp *= eta1; }
  if (m2_!=0) { AmmIm_tmp *= eta2; AmmRe_tmp *= eta2; }

  //-----------------------------------------------------------------------
  // No m1_/m2_ singularities in A++

  AppRe_tmp =
    pow(m12, 4)/sqrt(6.)*
    (
    c1Re/4.*(
    1. + 4.*xxp - pow(eta1sq-eta2sq, 2)
    )
    + c2Re*8.*xxp
    - c5Im*8.*pow(xxp, 1.5)
    + c6Im*4.*sqrt(xxp)
    );

  AppIm_tmp =
    pow(m12, 4)/sqrt(6.)*
    (
    c1Im/4.*(
    1. + 4.*xxp - pow(eta1sq-eta2sq, 2)
    )
    + c2Im*8.*xxp
    + c5Re*8.*pow(xxp, 1.5)
    - c6Re*4.*sqrt(xxp)
    );

  if (m1_!=0) { AppIm_tmp *= eta1; AppRe_tmp *= eta1; }
  if (m2_!=0) { AppIm_tmp *= eta2; AppRe_tmp *= eta2; }

  //-----------------------------------------------------------------------

  Double_t A0m_0p_c1factor = (
    -(pow(eta1, 6)-pow(eta2, 6))/8.
    + (eta1sq-eta2sq)*(3.*eta1p2sq + 4.*xxp)/8.
    - pow(eta1sq-eta2sq, 2)/8.
    + xxp/2.
    + (1. + (eta1sq-eta2sq))/8.
    );
  Double_t Am0_p0_c1factor = (
    (pow(eta1, 6)-pow(eta2, 6))/8.
    - (eta1sq-eta2sq)*(3.*eta1p2sq + 4.*xxp)/8.
    - pow(eta1sq-eta2sq, 2)/8.
    + xxp/2.
    + (1. - (eta1sq-eta2sq))/8.
    );
  Double_t A0m_0p_m0_p0_c4factor = 2.*xxp;
  Double_t A0m_0p_c6factor = sqrt(xxp)*(1.+eta1sq-eta2sq); // x+-i
  Double_t Am0_p0_c6factor = sqrt(xxp)*(1.-eta1sq+eta2sq); // x+-i
  Double_t A0m_0p_m0_p0_c7factor = 4.*pow(xxp, 1.5); // x+-i

  A0mRe_tmp =
    pow(m12, 4)*eta2*
    (
    c1Re*A0m_0p_c1factor
    + c42Re*A0m_0p_m0_p0_c4factor
    - c6Im*A0m_0p_c6factor
    - c7Im*A0m_0p_m0_p0_c7factor
    );

  A0mIm_tmp =
    pow(m12, 4)*eta2*
    (
    c1Im*A0m_0p_c1factor
    + c42Im*A0m_0p_m0_p0_c4factor
    + c6Re*A0m_0p_c6factor
    + c7Re*A0m_0p_m0_p0_c7factor
    );

  //-----------------------------------------------------------------------

  Am0Re_tmp =
    pow(m12, 4)*eta1*
    (
    c1Re*Am0_p0_c1factor
    + c41Re*A0m_0p_m0_p0_c4factor
    - c6Im*Am0_p0_c6factor
    - c7Im*A0m_0p_m0_p0_c7factor
    );

  Am0Im_tmp =
    pow(m12, 4)*eta1*
    (
    c1Im*Am0_p0_c1factor
    + c41Im*A0m_0p_m0_p0_c4factor
    + c6Re*Am0_p0_c6factor
    + c7Re*A0m_0p_m0_p0_c7factor
    );

  //-----------------------------------------------------------------------

  A0pRe_tmp =
    pow(m12, 4)*eta2*
    (
    c1Re*A0m_0p_c1factor
    + c42Re*A0m_0p_m0_p0_c4factor
    + c6Im*A0m_0p_c6factor
    + c7Im*A0m_0p_m0_p0_c7factor
    );

  A0pIm_tmp =
    pow(m12, 4)*eta2*
    (
    c1Im*A0m_0p_c1factor
    + c42Im*A0m_0p_m0_p0_c4factor
    - c6Re*A0m_0p_c6factor
    - c7Re*A0m_0p_m0_p0_c7factor
    );

  //-----------------------------------------------------------------------

  Ap0Re_tmp =
    pow(m12, 4)*eta1*
    (
    c1Re*Am0_p0_c1factor
    + c41Re*A0m_0p_m0_p0_c4factor
    + c6Im*Am0_p0_c6factor
    + c7Im*A0m_0p_m0_p0_c7factor
    );

  Ap0Im_tmp =
    pow(m12, 4)*eta1*
    (
    c1Im*Am0_p0_c1factor
    + c41Im*A0m_0p_m0_p0_c4factor
    - c6Re*Am0_p0_c6factor
    - c7Re*A0m_0p_m0_p0_c7factor
    );

  //-----------------------------------------------------------------------
  // No m1_/m2_ singularities in A+- or A-+

  AmpRe_tmp = pow(m12, 4)*(c1Re/4.*(1.+4.*xxp-pow(eta1sq-eta2sq, 2)));
  AmpIm_tmp = pow(m12, 4)*(c1Im/4.*(1.+4.*xxp-pow(eta1sq-eta2sq, 2)));
  if (m1_!=0) { AmpIm_tmp *= eta1; AmpRe_tmp *= eta1; }
  if (m2_!=0) { AmpIm_tmp *= eta2; AmpRe_tmp *= eta2; }
  ApmRe_tmp = AmpRe_tmp;
  ApmIm_tmp = AmpIm_tmp;


  A00Re = ((A00Re_tmp*propV1Re - A00Im_tmp*propV1Im)*propV2Re - (A00Re_tmp*propV1Im + A00Im_tmp*propV1Re)*propV2Im);
  A00Im = ((A00Re_tmp*propV1Re - A00Im_tmp*propV1Im)*propV2Im + (A00Re_tmp*propV1Im + A00Im_tmp*propV1Re)*propV2Re);
  AmmRe = ((AmmRe_tmp*propV1Re - AmmIm_tmp*propV1Im)*propV2Re - (AmmRe_tmp*propV1Im + AmmIm_tmp*propV1Re)*propV2Im);
  AmmIm = ((AmmRe_tmp*propV1Re - AmmIm_tmp*propV1Im)*propV2Im + (AmmRe_tmp*propV1Im + AmmIm_tmp*propV1Re)*propV2Re);
  AppRe = ((AppRe_tmp*propV1Re - AppIm_tmp*propV1Im)*propV2Re - (AppRe_tmp*propV1Im + AppIm_tmp*propV1Re)*propV2Im);
  AppIm = ((AppRe_tmp*propV1Re - AppIm_tmp*propV1Im)*propV2Im + (AppRe_tmp*propV1Im + AppIm_tmp*propV1Re)*propV2Re);
  A0mRe = ((A0mRe_tmp*propV1Re - A0mIm_tmp*propV1Im)*propV2Re - (A0mRe_tmp*propV1Im + A0mIm_tmp*propV1Re)*propV2Im);
  A0mIm = ((A0mRe_tmp*propV1Re - A0mIm_tmp*propV1Im)*propV2Im + (A0mRe_tmp*propV1Im + A0mIm_tmp*propV1Re)*propV2Re);
  A0pRe = ((A0pRe_tmp*propV1Re - A0pIm_tmp*propV1Im)*propV2Re - (A0pRe_tmp*propV1Im + A0pIm_tmp*propV1Re)*propV2Im);
  A0pIm = ((A0pRe_tmp*propV1Re - A0pIm_tmp*propV1Im)*propV2Im + (A0pRe_tmp*propV1Im + A0pIm_tmp*propV1Re)*propV2Re);
  Am0Re = ((Am0Re_tmp*propV1Re - Am0Im_tmp*propV1Im)*propV2Re - (Am0Re_tmp*propV1Im + Am0Im_tmp*propV1Re)*propV2Im);
  Am0Im = ((Am0Re_tmp*propV1Re - Am0Im_tmp*propV1Im)*propV2Im + (Am0Re_tmp*propV1Im + Am0Im_tmp*propV1Re)*propV2Re);
  Ap0Re = ((Ap0Re_tmp*propV1Re - Ap0Im_tmp*propV1Im)*propV2Re - (Ap0Re_tmp*propV1Im + Ap0Im_tmp*propV1Re)*propV2Im);
  Ap0Im = ((Ap0Re_tmp*propV1Re - Ap0Im_tmp*propV1Im)*propV2Im + (Ap0Re_tmp*propV1Im + Ap0Im_tmp*propV1Re)*propV2Re);
  AmpRe = ((AmpRe_tmp*propV1Re - AmpIm_tmp*propV1Im)*propV2Re - (AmpRe_tmp*propV1Im + AmpIm_tmp*propV1Re)*propV2Im);
  AmpIm = ((AmpRe_tmp*propV1Re - AmpIm_tmp*propV1Im)*propV2Im + (AmpRe_tmp*propV1Im + AmpIm_tmp*propV1Re)*propV2Re);
  ApmRe = ((ApmRe_tmp*propV1Re - ApmIm_tmp*propV1Im)*propV2Re - (ApmRe_tmp*propV1Im + ApmIm_tmp*propV1Re)*propV2Im);
  ApmIm = ((ApmRe_tmp*propV1Re - ApmIm_tmp*propV1Im)*propV2Im + (ApmRe_tmp*propV1Im + ApmIm_tmp*propV1Re)*propV2Re);


  if (
    A00Re!=A00Re || A00Im!=A00Im ||
    AppRe!=AppRe || AppIm!=AppIm ||
    AmmRe!=AmmRe || AmmIm!=AmmIm ||
    A0pRe!=A0pRe || A0pIm!=A0pIm ||
    A0mRe!=A0mRe || A0mIm!=A0mIm ||
    Ap0Re!=Ap0Re || Ap0Im!=Ap0Im ||
    Am0Re!=Am0Re || Am0Im!=Am0Im ||
    ApmRe!=ApmRe || ApmIm!=ApmIm ||
    AmpRe!=AmpRe || AmpIm!=AmpIm ||
    (
    A00Re==0 && A00Im==0 &&
    AppRe==0 && AppIm==0 &&
    AmmRe==0 && AmmIm==0 &&
    A0pRe==0 && A0pIm==0 &&
    A0mRe==0 && A0mIm==0 &&
    Ap0Re==0 && Ap0Im==0 &&
    Am0Re==0 && Am0Im==0 &&
    ApmRe==0 && ApmIm==0 &&
    AmpRe==0 && AmpIm==0
    )
    ){
    std::cerr << "Some of the amplitudes are NaN or all are 0:" << endl;
    std::cerr << "A00Re=" << A00Re << ", A00Im=" << A00Im << endl;
    std::cerr << "AppRe=" << AppRe << ", AppIm=" << AppIm << endl;
    std::cerr << "AmmRe=" << AmmRe << ", AmmIm=" << AmmIm << endl;
    std::cerr << "A0pRe=" << A0pRe << ", A0pIm=" << A0pIm << endl;
    std::cerr << "A0mRe=" << A0mRe << ", A0mIm=" << A0mIm << endl;
    std::cerr << "Ap0Re=" << Ap0Re << ", Ap0Im=" << Ap0Im << endl;
    std::cerr << "Am0Re=" << Am0Re << ", Am0Im=" << Am0Im << endl;
    std::cerr << "ApmRe=" << ApmRe << ", ApmIm=" << ApmIm << endl;
    std::cerr << "AmpRe=" << AmpRe << ", AmpIm=" << AmpIm << endl;

    std::cerr << "Possible causes are" << endl;
    std::cerr << "m12=" << m12 << endl;
    std::cerr << "m1=" << m1_ << endl;
    std::cerr << "m2=" << m2_ << endl;
    std::cerr << "x=" << x << endl;
    std::cerr << "xxp=" << xxp << endl;

    std::cerr << "c1 = (" << c1Re << ", " << c1Im << ")" << endl;
    std::cerr << "c2 = (" << c2Re << ", " << c2Im << ")" << endl;
    std::cerr << "c3 = (" << c3Re << ", " << c3Im << ")" << endl;
    std::cerr << "c41 = (" << c41Re << ", " << c41Im << ")" << endl;
    std::cerr << "c42 = (" << c42Re << ", " << c42Im << ")" << endl;
    std::cerr << "c5 = (" << c5Re << ", " << c5Im << ")" << endl;
    std::cerr << "c6 = (" << c6Re << ", " << c6Im << ")" << endl;
    std::cerr << "c7 = (" << c7Re << ", " << c7Im << ")" << endl;
  }

  return;
}

void RooSpinTwo::setProxies(modelMeasurables _measurables){
  setProxy(h1, (RooAbsReal*)_measurables.h1);
  setProxy(h2, (RooAbsReal*)_measurables.h2);
  setProxy(Phi, (RooAbsReal*)_measurables.Phi);
  setProxy(m1, (RooAbsReal*)_measurables.m1);
  setProxy(m2, (RooAbsReal*)_measurables.m2);
  setProxy(m12, (RooAbsReal*)_measurables.m12);
  setProxy(hs, (RooAbsReal*)_measurables.hs);
  setProxy(Phi1, (RooAbsReal*)_measurables.Phi1);
  //setProxy(Y, (RooAbsReal*)_measurables.Y);
}
void RooSpinTwo::setProxy(RooRealProxy& proxy, RooAbsReal* objectPtr){
  if (objectPtr!=0) proxy.setArg((RooAbsReal&)*objectPtr);
}
