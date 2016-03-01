#include "../include/RooSpinTwo.h"

RooSpinTwo::RooSpinTwo(
  const char* name, const char* title,
  modelMeasurables _measurables,
  modelParameters _parameters
  ) : RooAbsPdf(name, title),
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

void RooSpinTwo::calculateCi(std::vector<Double_t>& ciRe, std::vector<Double_t>& ciIm) const{
  Double_t s = (pow(m12, 2) - pow(m1, 2) - pow(m2, 2))/2.;
  if (pow(m1, 2)>pow(m12, 2) || pow(m2, 2)>pow(m12, 2)) s = -s;
  Double_t kappa = s/pow(Lambda, 2);

  Double_t c1Re = 2.*(b1Val + b2Val*kappa*(1.+m1*m1/s)*(1.+m2*m2/s) + b5Val*(mV*mV)/s); ciRe.push_back(c1Re);
  Double_t c2Re = -0.5*b1Val + b3Val*kappa*(1.-(m1*m1+m2*m2)/(2*s)) + 2.*b4Val*kappa + b7Val*kappa*mV*mV/s; ciRe.push_back(c2Re);
  Double_t c3Re = -1.*(b2Val/2.+b3Val+2.*b4Val)*kappa*m12*m12/s; ciRe.push_back(c3Re);
  Double_t c41Re = -b1Val - b2Val*kappa - (b2Val*m1*m1+b3Val*m2*m2+2.*b6Val*mV*mV)*kappa/s; ciRe.push_back(c41Re);
  Double_t c42Re = -b1Val - b2Val*kappa - (b2Val*m2*m2+b3Val*m1*m1+2.*b6Val*mV*mV)*kappa/s; ciRe.push_back(c42Re);
  Double_t c5Re = 2.*b8Val*kappa*(m12*m12)/s; ciRe.push_back(c5Re);
  Double_t c6Re = b9Val*kappa*(mV*mV)/s; ciRe.push_back(c6Re);
  Double_t c7Re = b10Val*kappa*kappa*(m12*m12*mV*mV)/(s*s); ciRe.push_back(c7Re);

  Double_t c1Im = 2.*(b1ValIm + b2ValIm*kappa*(1.+m1*m1/s)*(1.+m2*m2/s) + b5ValIm*(mV*mV)/s); ciIm.push_back(c1Im);
  Double_t c2Im = -0.5*b1ValIm + b3ValIm*kappa*(1.-(m1*m1+m2*m2)/(2*s)) + 2.*b4ValIm*kappa + b7ValIm*kappa*mV*mV/s; ciIm.push_back(c2Im);
  Double_t c3Im = -1.*(b2ValIm/2.+b3ValIm+2.*b4ValIm)*kappa*m12*m12/s; ciIm.push_back(c3Im);
  Double_t c41Im = -b1ValIm - b2ValIm*kappa - (b2ValIm*m1*m1+b3ValIm*m2*m2+2.*b6ValIm*mV*mV)*kappa/s; ciIm.push_back(c41Im);
  Double_t c42Im = -b1ValIm - b2ValIm*kappa - (b2ValIm*m2*m2+b3ValIm*m1*m1+2.*b6ValIm*mV*mV)*kappa/s; ciIm.push_back(c42Im);
  Double_t c5Im = 2.*b8ValIm*kappa*(m12*m12)/s; ciIm.push_back(c5Im);
  Double_t c6Im = b9ValIm*kappa*(mV*mV)/s; ciIm.push_back(c6Im);
  Double_t c7Im = b10ValIm*kappa*kappa*(m12*m12*mV*mV)/(s*s); ciIm.push_back(c7Im);
}
void RooSpinTwo::calculateAmplitudes(
  Double_t& A00Re, Double_t& A00Im,
  Double_t& AppRe, Double_t& AppIm, Double_t& A0pRe, Double_t& A0pIm, Double_t& Ap0Re, Double_t& Ap0Im,
  Double_t& AmmRe, Double_t& AmmIm, Double_t& A0mRe, Double_t& A0mIm, Double_t& Am0Re, Double_t& Am0Im,
  Double_t& ApmRe, Double_t& ApmIm, Double_t& AmpRe, Double_t& AmpIm
  )const{
  std::vector<Double_t> ciRe;
  std::vector<Double_t> ciIm;
  calculateCi(ciRe, ciIm);

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

  Double_t x = (m12*m12-m1*m1-m2*m2)/(2.0*m1*m2);

  A00Re =
    pow(m12, -4)*pow(sqrt(6.), -1)*
    (
    c1Re*pow(m1, 3)*pow(m2, 3) * (3./4. + (x*x-1.))
    + c2Re*pow(m1, 3)*pow(m2, 3) * (-4.*(x*x-1.) - 8*pow(x*x-1., 2))
    + c3Re*pow(m1, 3)*pow(m2, 3) * (-8.*pow(x*x-1., 2))
    + c1Re*m1*pow(m2, 5) * (-1.0/2.0 - 1.0/2.0*(x*x-1.))
    + c2Re*m1*pow(m2, 5) * (2.*(x*x-1.))
    + c1Re*pow(m1, 5)*m2 * (-1.0/2.0 - 1.0/2.0*(x*x-1.))
    + c2Re*pow(m1, 5)*m2 * (2.*(x*x-1.))
    + c1Re*(1.0/m1)*pow(m2, 7) * (1.0/8.0)
    + c1Re*pow(m1, 7)*(1.0/m2) * (1.0/8.0)
    )
    + pow(sqrt(6.), -1)*
    (
    c1Re*m1*m2 * (1.0/2.0 + 1.0/2.0*(x*x-1.))
    + c2Re*m1*m2 * (-2.*(x*x-1.))
    + (c41Re+c42Re)*m1*m2 * (2.*(x*x-1.))
    + m1*m2*(m1*m1-m2*m2)*pow(m12, -2)*(c41Re-c42Re)*2.*(x*x-1.)
    + c1Re*(1.0/m1)*pow(m2, 3) * (-1.0/4.0)
    + c1Re*pow(m1, 3)*(1.0/m2) * (-1.0/4.0)
    + c1Re*pow(m12, 4)/(8.*m1*m2)
    );

  A00Im =
    pow(m12, -4)*pow(sqrt(6.), -1)*
    (
    c1Im*pow(m1, 3)*pow(m2, 3) * (3./4. + (x*x-1.))
    + c2Im*pow(m1, 3)*pow(m2, 3) * (-4.*(x*x-1.) - 8*pow(x*x-1., 2))
    + c3Im*pow(m1, 3)*pow(m2, 3) * (-8.*pow(x*x-1., 2))
    + c1Im*m1*pow(m2, 5) * (-1.0/2.0 - 1.0/2.0*(x*x-1.))
    + c2Im*m1*pow(m2, 5) * (2.*(x*x-1.))
    + c1Im*pow(m1, 5)*m2 * (-1.0/2.0 - 1.0/2.0*(x*x-1.))
    + c2Im*pow(m1, 5)*m2 * (2.*(x*x-1.))
    + c1Im*(1.0/m1)*pow(m2, 7) * (1.0/8.0)
    + c1Im*pow(m1, 7)*(1.0/m2) * (1.0/8.0)
    )
    + pow(sqrt(6.), -1)*
    (
    c1Im*m1*m2 * (1.0/2.0 + 1.0/2.0*(x*x-1.))
    + c2Im*m1*m2 * (-2.*(x*x-1.))
    + (c41Im+c42Im)*m1*m2 * (2.*(x*x-1.))
    + m1*m2*(m1*m1-m2*m2)*pow(m12, -2)*(c41Im-c42Im)*2.*(x*x-1.)
    + c1Im*(1.0/m1)*pow(m2, 3) * (-1.0/4.0)
    + c1Im*pow(m1, 3)*(1.0/m2) * (-1.0/4.0)
    + c1Im*pow(m12, 4)/(8.*m1*m2)
    );

  //-----------------------------------------------------------------------

  AppRe =
    pow(sqrt(6.), -1)*
    (
    pow(m12, 2)*c1Re * (1./4.)
    + pow(m12, -2)*
    (
    c1Re*pow(m2, 4) * (-1.0/4.0)
    + c1Re*pow(m1, 4) * (-1.0/4.0)
    + c1Re*pow(m1, 2)*pow(m2, 2) * (1.0/2.0 + (x*x-1.))
    + c2Re*pow(m1, 2)*pow(m2, 2) * (8.*(x*x-1.))
    )
    + pow(m12, -4)*c5Im*pow(m1, 3)*pow(m2, 3) * (8.*pow(sqrt(x*x-1.), 3))
    + c6Im*m1*m2 * (-4.*sqrt(x*x-1.))
    );

  AppIm =
    pow(sqrt(6.), -1)*
    (
    pow(m12, 2)*c1Im * (1./4.)
    + pow(m12, -2)*
    (
    c1Im*pow(m2, 4) * (-1.0/4.0)
    + c1Im*pow(m1, 4) * (-1.0/4.0)
    + c1Im*pow(m1, 2)*pow(m2, 2) * (1.0/2.0 + (x*x-1.))
    + c2Im*pow(m1, 2)*pow(m2, 2) * (8.*(x*x-1.))
    )
    + pow(m12, -4)*c5Re*pow(m1, 3)*pow(m2, 3) * (8.*pow(sqrt(x*x-1.), 3))
    + c6Re*m1*m2 * (-4.*sqrt(x*x-1.))
    );

  //-----------------------------------------------------------------------
  AmmRe =
    pow(sqrt(6.), -1)*
    (
    pow(m12, 2)*c1Re * (1.0/4.0)
    + pow(m12, -2)*
    (
    c1Re*pow(m1, 4) * (-1.0/4.0)
    + c1Re*pow(m2, 4) * (-1.0/4.0)
    + c1Re*pow(m1, 2)*pow(m2, 2) * (1.0/2.0 + (x*x-1.))
    + c2Re*pow(m1, 2)*pow(m2, 2) * (8.*(x*x-1.))
    )
    + pow(m12, -4)*c5Im*pow(m1, 3)*pow(m2, 3) * (-8.*pow(sqrt(x*x-1.), 3))
    + c6Im*m1*m2 * (4.*sqrt(x*x-1.))
    );

  AmmIm =
    pow(sqrt(6.), -1)*
    (
    pow(m12, 2)*c1Im * (1.0/4.0)
    + pow(m12, -2)*
    (
    c1Im*pow(m1, 4) * (-1.0/4.0)
    + c1Im*pow(m2, 4) * (-1.0/4.0)
    + c1Im*pow(m1, 2)*pow(m2, 2) * (1.0/2.0 + (x*x-1.))
    + c2Im*pow(m1, 2)*pow(m2, 2) * (8.*(x*x-1.))
    )
    + pow(m12, -4)*c5Re*pow(m1, 3)*pow(m2, 3) * (-8.*pow(sqrt(x*x-1.), 3))
    + c6Re*m1*m2 * (4.*sqrt(x*x-1.))
    );

  //-----------------------------------------------------------------------

  Ap0Re =
    (1.0/sqrt(2.))*
    (
    (1.0/pow(m12, 3))*
    (
    c1Re*pow(m2, 5) * (-1.0/8.0)
    + c1Re*pow(m1, 2)*pow(m2, 3) * (3.0/8.0 + 1.0/2.0*(x*x-1.))
    + c1Re*pow(m1, 4)*m2 * (-3.0/8.0 - 1.0/2.0*(x*x-1.))
    + c1Re*pow(m1, 6)*(1.0/m2) * (1.0/8.0)
    + c7Im*pow(m1, 3)*pow(m2, 2) * (-4.*pow(sqrt(x*x-1.), 3)
    )
    + (1.0/m12)*
    (
    c1Re*pow(m1, 2)*m2 * (1.0/4.0 + 1.0/2.0*(x*x-1.))
    + c41Re*pow(m1, 2)*m2 * (2.*(x*x-1.))
    + c1Re*pow(m2, 3) * (-1.0/8.0)
    + c1Re*pow(m1, 4)*(1.0/m2) * (-1.0/8.0)
    + c6Im*(pow(m1, 3)-m1*pow(m2, 2)) * sqrt(x*x-1.)
    )
    + m12/(8.*m2)*
    (
    c1Re* (pow(m2, 2) - pow(m1, 2) + pow(m12, 2))
    - c6Im* 8.*m1*m2 * sqrt(x*x-1.)
    )
    )
    );

  Ap0Im =
    (1.0/sqrt(2.))*
    (
    (1.0/pow(m12, 3))*
    (
    c1Im*pow(m2, 5) * (-1.0/8.0)
    + c1Im*pow(m1, 2)*pow(m2, 3) * (3.0/8.0 + 1.0/2.0*(x*x-1.))
    + c1Im*pow(m1, 4)*m2 * (-3.0/8.0 - 1.0/2.0*(x*x-1.))
    + c1Im*pow(m1, 6)*(1.0/m2) * (1.0/8.0)
    + c7Re*pow(m1, 3)*pow(m2, 2) * (-4.*pow(sqrt(x*x-1.), 3)
    )
    + (1.0/m12)*
    (
    c1Im*pow(m1, 2)*m2 * (1.0/4.0 + 1.0/2.0*(x*x-1.))
    + c41Im*pow(m1, 2)*m2 * (2.*(x*x-1.))
    + c1Im*pow(m2, 3) * (-1.0/8.0)
    + c1Im*pow(m1, 4)*(1.0/m2) * (-1.0/8.0)
    + c6Re*(pow(m1, 3)-m1*pow(m2, 2)) * sqrt(x*x-1.)
    )
    + m12/(8.*m2)*
    (
    c1Im* (pow(m2, 2) - pow(m1, 2) + pow(m12, 2))
    - c6Re* 8.*m1*m2 * sqrt(x*x-1.)
    )
    )
    );

  //-----------------------------------------------------------------------

  A0pRe =
    (1.0/sqrt(2.))*
    (
    (1.0/pow(m12, 3))*
    (
    c1Re*pow(m1, 5) * (-1.0/8.0)
    + c1Re*pow(m2, 2)*pow(m1, 3) * (3.0/8.0 + 1.0/2.0*(x*x-1.))
    + c1Re*pow(m2, 4)*m1 * (-3.0/8.0 - 1.0/2.0*(x*x-1.))
    + c1Re*pow(m2, 6)*(1.0/m1) * (1.0/8.0)
    + c7Im*pow(m2, 3)*pow(m1, 2) * (-4.*pow(sqrt(x*x-1.), 3)
    )
    + (1.0/m12)*
    (
    c1Re*pow(m2, 2)*m1 * (1.0/4.0 + 1.0/2.0*(x*x-1.))
    + c1Re*pow(m1, 3) * (-1.0/8.0)
    + c1Re*pow(m2, 4)*(1.0/m1) * (-1.0/8.0)
    + c41Re*pow(m2, 2)*m1 * (2.*(x*x-1.))
    + c6Im*(pow(m2, 3)-m2*pow(m1, 2)) * sqrt(x*x-1.)
    )
    + m12/(8.*m1)*
    (
    c1Re* (pow(m1, 2) - pow(m2, 2) + pow(m12, 2))
    - c6Im* 8.*m2*m1 * sqrt(x*x-1.)
    )
    )
    );

  A0pIm =
    (1.0/sqrt(2.))*
    (
    (1.0/pow(m12, 3))*
    (
    c1Im*pow(m1, 5) * (-1.0/8.0)
    + c1Im*pow(m2, 2)*pow(m1, 3) * (3.0/8.0 + 1.0/2.0*(x*x-1.))
    + c1Im*pow(m2, 4)*m1 * (-3.0/8.0 - 1.0/2.0*(x*x-1.))
    + c1Im*pow(m2, 6)*(1.0/m1) * (1.0/8.0)
    + c7Re*pow(m2, 3)*pow(m1, 2) * (-4.*pow(sqrt(x*x-1.), 3)
    )
    + (1.0/m12)*
    (
    c1Im*pow(m2, 2)*m1 * (1.0/4.0 + 1.0/2.0*(x*x-1.))
    + c1Im*pow(m1, 3) * (-1.0/8.0)
    + c1Im*pow(m2, 4)*(1.0/m1) * (-1.0/8.0)
    + c41Im*pow(m2, 2)*m1 * (2.*(x*x-1.))
    + c6Re*(pow(m2, 3)-m2*pow(m1, 2)) * sqrt(x*x-1.)
    )
    + m12/(8.*m1)*
    (
    c1Im* (pow(m1, 2) - pow(m2, 2) + pow(m12, 2))
    - c6Re* 8.*m2*m1 * sqrt(x*x-1.)
    )
    )
    );


  //-----------------------------------------------------------------------

  A0mRe =
    (1.0/sqrt(2.))*
    (
    (1.0/pow(m12, 3))*
    (
      c1Re*(1.0/m1)*pow(m2, 6) * (1.0/8.0)
    + c1Re*m1*pow(m2, 4) * (-3.0/8.0 - 1.0/2.0*(x*x-1.))
    + c1Re*pow(m1, 3)*pow(m2, 2) * (3.0/8.0 + 1.0/2.0*(x*x-1.))
    + c1Re*pow(m1, 5) * (-1.0/8.0)
    + c7Im*pow(m1, 2)*pow(m2, 3) * (4*pow(sqrt(x*x-1.), 3)
    )
    + (1.0/m12)*
    (
    c1Re*m1*pow(m2, 2) * (1.0/4.0 + 1.0/2.0*(x*x-1.))
    + c1Re*(1.0/m1)*pow(m2, 4) * (-1.0/8.0)
    + c1Re*pow(m1, 3) * (-1.0/8.0)
    + c42Re*m1*pow(m2, 2) * (2*(x*x-1.))
    + c6Im*(pow(m1, 2)*m2-pow(m2, 3))*sqrt(x*x-1.)
    )
    + m12/(8.*m1)*
    (
    c1Re* (pow(m1,2) - pow(m2, 2) + pow(m12,2) )
    + c6Im* 8.*m1*m2*sqrt(x*x-1.)
    )
    )
    );

  A0mIm =
    (1.0/sqrt(2.))*
    (
    (1.0/pow(m12, 3))*
    (
    c1Im*(1.0/m1)*pow(m2, 6) * (1.0/8.0)
    + c1Im*m1*pow(m2, 4) * (-3.0/8.0 - 1.0/2.0*(x*x-1.))
    + c1Im*pow(m1, 3)*pow(m2, 2) * (3.0/8.0 + 1.0/2.0*(x*x-1.))
    + c1Im*pow(m1, 5) * (-1.0/8.0)
    + c7Re*pow(m1, 2)*pow(m2, 3) * (4*pow(sqrt(x*x-1.), 3)
    )
    + (1.0/m12)*
    (
    c1Im*m1*pow(m2, 2) * (1.0/4.0 + 1.0/2.0*(x*x-1.))
    + c1Im*(1.0/m1)*pow(m2, 4) * (-1.0/8.0)
    + c1Im*pow(m1, 3) * (-1.0/8.0)
    + c42Im*m1*pow(m2, 2) * (2*(x*x-1.))
    + c6Re*(pow(m1, 2)*m2-pow(m2, 3))*sqrt(x*x-1.)
    )
    + m12/(8.*m1)*
    (
    c1Im* (pow(m1, 2) - pow(m2, 2) + pow(m12, 2))
    + c6Re* 8.*m1*m2*sqrt(x*x-1.)
    )
    )
    );

  //-----------------------------------------------------------------------

  Am0Re =
    (1.0/sqrt(2.))*
    (
    (1.0/pow(m12, 3))*
    (
    c1Re*(1.0/m2)*pow(m1, 6) * (1.0/8.0)
    + c1Re*m2*pow(m1, 4) * (-3.0/8.0 - 1.0/2.0*(x*x-1.))
    + c1Re*pow(m2, 3)*pow(m1, 2) * (3.0/8.0 + 1.0/2.0*(x*x-1.))
    + c1Re*pow(m2, 5) * (-1.0/8.0)
    + c7Im*pow(m2, 2)*pow(m1, 3) * (4*pow(sqrt(x*x-1.), 3)
    )
    + (1.0/m12)*
    (
    c1Re*m2*pow(m1, 2) * (1.0/4.0 + 1.0/2.0*(x*x-1.))
    + c1Re*(1.0/m2)*pow(m1, 4) * (-1.0/8.0)
    + c1Re*pow(m2, 3) * (-1.0/8.0)
    + c42Re*m2*pow(m1, 2) * (2*(x*x-1.))
    + c6Im*(pow(m2, 2)*m1-pow(m1, 3))*sqrt(x*x-1.)
    )
    + m12/(8.*m2)*
    (
    c1Re* (pow(m2, 2) - pow(m1, 2) + pow(m12, 2))
    + c6Im* 8.*m2*m1*sqrt(x*x-1.)
    )
    )
    );

  Am0Im =
    (1.0/sqrt(2.))*
    (
    (1.0/pow(m12, 3))*
    (
    c1Im*(1.0/m2)*pow(m1, 6) * (1.0/8.0)
    + c1Im*m2*pow(m1, 4) * (-3.0/8.0 - 1.0/2.0*(x*x-1.))
    + c1Im*pow(m2, 3)*pow(m1, 2) * (3.0/8.0 + 1.0/2.0*(x*x-1.))
    + c1Im*pow(m2, 5) * (-1.0/8.0)
    + c7Re*pow(m2, 2)*pow(m1, 3) * (4*pow(sqrt(x*x-1.), 3)
    )
    + (1.0/m12)*
    (
    c1Im*m2*pow(m1, 2) * (1.0/4.0 + 1.0/2.0*(x*x-1.))
    + c1Im*(1.0/m2)*pow(m1, 4) * (-1.0/8.0)
    + c1Im*pow(m2, 3) * (-1.0/8.0)
    + c42Im*m2*pow(m1, 2) * (2*(x*x-1.))
    + c6Re*(pow(m2, 2)*m1-pow(m1, 3))*sqrt(x*x-1.)
    )
    + m12/(8.*m2)*
    (
    c1Im* (pow(m2, 2) - pow(m1, 2) + pow(m12, 2))
    + c6Re* 8.*m2*m1*sqrt(x*x-1.)
    )
    )
    );


  //-----------------------------------------------------------------------

  ApmRe =
    pow(m12, -2)*c1Re/4.*
    (
    pow(m12, 4) - pow(m1, 4) - pow(m2, 4)
    + 4.*pow(m1, 2)*pow(m2, 2) * (x*x-0.5)
    );
  AmpRe = ApmRe;
  ApmIm =
    pow(m12, -2)*c1Im/4.*
    (
    pow(m12, 4) - pow(m1, 4) - pow(m2, 4)
    + 4.*pow(m1, 2)*pow(m2, 2) * (x*x-0.5)
    );
  AmpIm = ApmIm;

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
