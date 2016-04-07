#ifdef _def_melatools_
#include <ZZMatrixElement/MELA/interface/RooSpinZero.h>
#else
#include "../include/RooSpinZero.h"
#endif


RooSpinZero::RooSpinZero(
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
  Y("Y", "Y", this),

  mX("mX", "mX", this, (RooAbsReal&)*(_parameters.mX)),
  gamX("gamX", "gamX", this, (RooAbsReal&)*(_parameters.gamX)),
  mV("mV", "mV", this, (RooAbsReal&)*(_parameters.mV)),
  gamV("gamV", "gamV", this, (RooAbsReal&)*(_parameters.gamV)),
  R1Val("R1Val", "R1Val", this, (RooAbsReal&)*(_parameters.R1Val)),
  R2Val("R2Val", "R2Val", this, (RooAbsReal&)*(_parameters.R2Val)),

  g1Val("g1Val", "g1Val", this, (RooAbsReal&)*(_parameters.g1List[0][0])),
  g2Val("g2Val", "g2Val", this, (RooAbsReal&)*(_parameters.g2List[0][0])),
  g3Val("g3Val", "g3Val", this, (RooAbsReal&)*(_parameters.g3List[0][0])),
  g4Val("g4Val", "g4Val", this, (RooAbsReal&)*(_parameters.g4List[0][0])),

  g1_primeVal("g1_primeVal", "g1_primeVal", this, (RooAbsReal&)*(_parameters.g1List[1][0])),
  g2_primeVal("g2_primeVal", "g2_primeVal", this, (RooAbsReal&)*(_parameters.g2List[1][0])),
  g3_primeVal("g3_primeVal", "g3_primeVal", this, (RooAbsReal&)*(_parameters.g3List[1][0])),
  g4_primeVal("g4_primeVal", "g4_primeVal", this, (RooAbsReal&)*(_parameters.g4List[1][0])),

  g1_prime2Val("g1_prime2Val", "g1_prime2Val", this, (RooAbsReal&)*(_parameters.g1List[2][0])),
  g2_prime2Val("g2_prime2Val", "g2_prime2Val", this, (RooAbsReal&)*(_parameters.g2List[2][0])),
  g3_prime2Val("g3_prime2Val", "g3_prime2Val", this, (RooAbsReal&)*(_parameters.g3List[2][0])),
  g4_prime2Val("g4_prime2Val", "g4_prime2Val", this, (RooAbsReal&)*(_parameters.g4List[2][0])),

  g1_prime3Val("g1_prime3Val", "g1_prime3Val", this, (RooAbsReal&)*(_parameters.g1List[3][0])),
  g2_prime3Val("g2_prime3Val", "g2_prime3Val", this, (RooAbsReal&)*(_parameters.g2List[3][0])),
  g3_prime3Val("g3_prime3Val", "g3_prime3Val", this, (RooAbsReal&)*(_parameters.g3List[3][0])),
  g4_prime3Val("g4_prime3Val", "g4_prime3Val", this, (RooAbsReal&)*(_parameters.g4List[3][0])),

  g1_prime4Val("g1_prime4Val", "g1_prime4Val", this, (RooAbsReal&)*(_parameters.g1List[4][0])),
  g2_prime4Val("g2_prime4Val", "g2_prime4Val", this, (RooAbsReal&)*(_parameters.g2List[4][0])),
  g3_prime4Val("g3_prime4Val", "g3_prime4Val", this, (RooAbsReal&)*(_parameters.g3List[4][0])),
  g4_prime4Val("g4_prime4Val", "g4_prime4Val", this, (RooAbsReal&)*(_parameters.g4List[4][0])),

  g1_prime5Val("g1_prime5Val", "g1_prime5Val", this, (RooAbsReal&)*(_parameters.g1List[5][0])),
  g2_prime5Val("g2_prime5Val", "g2_prime5Val", this, (RooAbsReal&)*(_parameters.g2List[5][0])),
  g3_prime5Val("g3_prime5Val", "g3_prime5Val", this, (RooAbsReal&)*(_parameters.g3List[5][0])),
  g4_prime5Val("g4_prime5Val", "g4_prime5Val", this, (RooAbsReal&)*(_parameters.g4List[5][0])),

  g1_prime6Val("g1_prime6Val", "g1_prime6Val", this, (RooAbsReal&)*(_parameters.g1List[6][0])),
  g2_prime6Val("g2_prime6Val", "g2_prime6Val", this, (RooAbsReal&)*(_parameters.g2List[6][0])),
  g3_prime6Val("g3_prime6Val", "g3_prime6Val", this, (RooAbsReal&)*(_parameters.g3List[6][0])),
  g4_prime6Val("g4_prime6Val", "g4_prime6Val", this, (RooAbsReal&)*(_parameters.g4List[6][0])),

  g1_prime7Val("g1_prime7Val", "g1_prime7Val", this, (RooAbsReal&)*(_parameters.g1List[7][0])),
  g2_prime7Val("g2_prime7Val", "g2_prime7Val", this, (RooAbsReal&)*(_parameters.g2List[7][0])),
  g3_prime7Val("g3_prime7Val", "g3_prime7Val", this, (RooAbsReal&)*(_parameters.g3List[7][0])),
  g4_prime7Val("g4_prime7Val", "g4_prime7Val", this, (RooAbsReal&)*(_parameters.g4List[7][0])),

  gzgs1_prime2Val("gzgs1_prime2Val", "gzgs1_prime2Val", this, (RooAbsReal&)*(_parameters.gzgs1List[0][0])), // Special case!
  gzgs2Val("gzgs2Val", "gzgs2Val", this, (RooAbsReal&)*(_parameters.gzgs2List[0][0])),
  gzgs3Val("gzgs3Val", "gzgs3Val", this, (RooAbsReal&)*(_parameters.gzgs3List[0][0])),
  gzgs4Val("gzgs4Val", "gzgs4Val", this, (RooAbsReal&)*(_parameters.gzgs4List[0][0])),
  ggsgs2Val("ggsgs2Val", "ggsgs2Val", this, (RooAbsReal&)*(_parameters.ggsgs2List[0][0])),
  ggsgs3Val("ggsgs3Val", "ggsgs3Val", this, (RooAbsReal&)*(_parameters.ggsgs3List[0][0])),
  ggsgs4Val("ggsgs4Val", "ggsgs4Val", this, (RooAbsReal&)*(_parameters.ggsgs4List[0][0])),

  g1ValIm("g1ValIm", "g1ValIm", this, (RooAbsReal&)*(_parameters.g1List[0][1])),
  g2ValIm("g2ValIm", "g2ValIm", this, (RooAbsReal&)*(_parameters.g2List[0][1])),
  g3ValIm("g3ValIm", "g3ValIm", this, (RooAbsReal&)*(_parameters.g3List[0][1])),
  g4ValIm("g4ValIm", "g4ValIm", this, (RooAbsReal&)*(_parameters.g4List[0][1])),

  g1_primeValIm("g1_primeValIm", "g1_primeValIm", this, (RooAbsReal&)*(_parameters.g1List[1][1])),
  g2_primeValIm("g2_primeValIm", "g2_primeValIm", this, (RooAbsReal&)*(_parameters.g2List[1][1])),
  g3_primeValIm("g3_primeValIm", "g3_primeValIm", this, (RooAbsReal&)*(_parameters.g3List[1][1])),
  g4_primeValIm("g4_primeValIm", "g4_primeValIm", this, (RooAbsReal&)*(_parameters.g4List[1][1])),

  g1_prime2ValIm("g1_prime2ValIm", "g1_prime2ValIm", this, (RooAbsReal&)*(_parameters.g1List[2][1])),
  g2_prime2ValIm("g2_prime2ValIm", "g2_prime2ValIm", this, (RooAbsReal&)*(_parameters.g2List[2][1])),
  g3_prime2ValIm("g3_prime2ValIm", "g3_prime2ValIm", this, (RooAbsReal&)*(_parameters.g3List[2][1])),
  g4_prime2ValIm("g4_prime2ValIm", "g4_prime2ValIm", this, (RooAbsReal&)*(_parameters.g4List[2][1])),

  g1_prime3ValIm("g1_prime3ValIm", "g1_prime3ValIm", this, (RooAbsReal&)*(_parameters.g1List[3][1])),
  g2_prime3ValIm("g2_prime3ValIm", "g2_prime3ValIm", this, (RooAbsReal&)*(_parameters.g2List[3][1])),
  g3_prime3ValIm("g3_prime3ValIm", "g3_prime3ValIm", this, (RooAbsReal&)*(_parameters.g3List[3][1])),
  g4_prime3ValIm("g4_prime3ValIm", "g4_prime3ValIm", this, (RooAbsReal&)*(_parameters.g4List[3][1])),

  g1_prime4ValIm("g1_prime4ValIm", "g1_prime4ValIm", this, (RooAbsReal&)*(_parameters.g1List[4][1])),
  g2_prime4ValIm("g2_prime4ValIm", "g2_prime4ValIm", this, (RooAbsReal&)*(_parameters.g2List[4][1])),
  g3_prime4ValIm("g3_prime4ValIm", "g3_prime4ValIm", this, (RooAbsReal&)*(_parameters.g3List[4][1])),
  g4_prime4ValIm("g4_prime4ValIm", "g4_prime4ValIm", this, (RooAbsReal&)*(_parameters.g4List[4][1])),

  g1_prime5ValIm("g1_prime5ValIm", "g1_prime5ValIm", this, (RooAbsReal&)*(_parameters.g1List[5][1])),
  g2_prime5ValIm("g2_prime5ValIm", "g2_prime5ValIm", this, (RooAbsReal&)*(_parameters.g2List[5][1])),
  g3_prime5ValIm("g3_prime5ValIm", "g3_prime5ValIm", this, (RooAbsReal&)*(_parameters.g3List[5][1])),
  g4_prime5ValIm("g4_prime5ValIm", "g4_prime5ValIm", this, (RooAbsReal&)*(_parameters.g4List[5][1])),

  g1_prime6ValIm("g1_prime6ValIm", "g1_prime6ValIm", this, (RooAbsReal&)*(_parameters.g1List[6][1])),
  g2_prime6ValIm("g2_prime6ValIm", "g2_prime6ValIm", this, (RooAbsReal&)*(_parameters.g2List[6][1])),
  g3_prime6ValIm("g3_prime6ValIm", "g3_prime6ValIm", this, (RooAbsReal&)*(_parameters.g3List[6][1])),
  g4_prime6ValIm("g4_prime6ValIm", "g4_prime6ValIm", this, (RooAbsReal&)*(_parameters.g4List[6][1])),

  g1_prime7ValIm("g1_prime7ValIm", "g1_prime7ValIm", this, (RooAbsReal&)*(_parameters.g1List[7][1])),
  g2_prime7ValIm("g2_prime7ValIm", "g2_prime7ValIm", this, (RooAbsReal&)*(_parameters.g2List[7][1])),
  g3_prime7ValIm("g3_prime7ValIm", "g3_prime7ValIm", this, (RooAbsReal&)*(_parameters.g3List[7][1])),
  g4_prime7ValIm("g4_prime7ValIm", "g4_prime7ValIm", this, (RooAbsReal&)*(_parameters.g4List[7][1])),

  gzgs1_prime2ValIm("gzgs1_prime2ValIm", "gzgs1_prime2ValIm", this, (RooAbsReal&)*(_parameters.gzgs1List[0][1])), // Special case!
  gzgs2ValIm("gzgs2ValIm", "gzgs2ValIm", this, (RooAbsReal&)*(_parameters.gzgs2List[0][1])),
  gzgs3ValIm("gzgs3ValIm", "gzgs3ValIm", this, (RooAbsReal&)*(_parameters.gzgs3List[0][1])),
  gzgs4ValIm("gzgs4ValIm", "gzgs4ValIm", this, (RooAbsReal&)*(_parameters.gzgs4List[0][1])),
  ggsgs2ValIm("ggsgs2ValIm", "ggsgs2ValIm", this, (RooAbsReal&)*(_parameters.ggsgs2List[0][1])),
  ggsgs3ValIm("ggsgs3ValIm", "ggsgs3ValIm", this, (RooAbsReal&)*(_parameters.ggsgs3List[0][1])),
  ggsgs4ValIm("ggsgs4ValIm", "ggsgs4ValIm", this, (RooAbsReal&)*(_parameters.ggsgs4List[0][1])),

  Lambda("Lambda", "Lambda", this, (RooAbsReal&)*(_parameters.Lambda)),
  Lambda_zgs1("Lambda_zgs1", "Lambda_zgs1", this, (RooAbsReal&)*(_parameters.Lambda_zgs1)),
  Lambda_z1("Lambda_z1", "Lambda_z1", this, (RooAbsReal&)*(_parameters.Lambda_z1)),
  Lambda_z2("Lambda_z2", "Lambda_z2", this, (RooAbsReal&)*(_parameters.Lambda_z2)),
  Lambda_z3("Lambda_z3", "Lambda_z3", this, (RooAbsReal&)*(_parameters.Lambda_z3)),
  Lambda_z4("Lambda_z4", "Lambda_z4", this, (RooAbsReal&)*(_parameters.Lambda_z4)),
  Lambda_Q("Lambda_Q", "Lambda_Q", this, (RooAbsReal&)*(_parameters.Lambda_Q)),

  Lambda_z11("Lambda_z11", "Lambda_z11", this, (RooAbsReal&)*(_parameters.Lambda_z1qsq[0])),
  Lambda_z21("Lambda_z21", "Lambda_z21", this, (RooAbsReal&)*(_parameters.Lambda_z2qsq[0])),
  Lambda_z31("Lambda_z31", "Lambda_z31", this, (RooAbsReal&)*(_parameters.Lambda_z3qsq[0])),
  Lambda_z41("Lambda_z41", "Lambda_z41", this, (RooAbsReal&)*(_parameters.Lambda_z4qsq[0])),

  Lambda_z12("Lambda_z12", "Lambda_z12", this, (RooAbsReal&)*(_parameters.Lambda_z1qsq[1])),
  Lambda_z22("Lambda_z22", "Lambda_z22", this, (RooAbsReal&)*(_parameters.Lambda_z2qsq[1])),
  Lambda_z32("Lambda_z32", "Lambda_z32", this, (RooAbsReal&)*(_parameters.Lambda_z3qsq[1])),
  Lambda_z42("Lambda_z42", "Lambda_z42", this, (RooAbsReal&)*(_parameters.Lambda_z4qsq[1])),

  Lambda_z10("Lambda_z10", "Lambda_z10", this, (RooAbsReal&)*(_parameters.Lambda_z1qsq[2])),
  Lambda_z20("Lambda_z20", "Lambda_z20", this, (RooAbsReal&)*(_parameters.Lambda_z2qsq[2])),
  Lambda_z30("Lambda_z30", "Lambda_z30", this, (RooAbsReal&)*(_parameters.Lambda_z3qsq[2])),
  Lambda_z40("Lambda_z40", "Lambda_z40", this, (RooAbsReal&)*(_parameters.Lambda_z4qsq[2])),

  cz_q1sq("cz_q1sq", "cz_q1sq", this, (RooAbsReal&)*(_parameters.cLambda_qsq[0])),
  cz_q2sq("cz_q2sq", "cz_q2sq", this, (RooAbsReal&)*(_parameters.cLambda_qsq[1])),
  cz_q12sq("cz_q12sq", "cz_q12sq", this, (RooAbsReal&)*(_parameters.cLambda_qsq[2]))
{
  setProxies(_measurables);
}


RooSpinZero::RooSpinZero(const RooSpinZero& other, const char* name) :
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
R2Val("R2Val", this, other.R2Val),

g1Val("g1Val", this, other.g1Val),
g2Val("g2Val", this, other.g2Val),
g3Val("g3Val", this, other.g3Val),
g4Val("g4Val", this, other.g4Val),

g1_primeVal("g1_primeVal", this, other.g1_primeVal),
g2_primeVal("a2_primeVal", this, other.g2_primeVal),
g3_primeVal("g3_primeVal", this, other.g3_primeVal),
g4_primeVal("g4_primeVal", this, other.g4_primeVal),

g1_prime2Val("g1_prime2Val", this, other.g1_prime2Val),
g2_prime2Val("a2_prime2Val", this, other.g2_prime2Val),
g3_prime2Val("g3_prime2Val", this, other.g3_prime2Val),
g4_prime2Val("g4_prime2Val", this, other.g4_prime2Val),

g1_prime3Val("g1_prime3Val", this, other.g1_prime3Val),
g2_prime3Val("a2_prime3Val", this, other.g2_prime3Val),
g3_prime3Val("g3_prime3Val", this, other.g3_prime3Val),
g4_prime3Val("g4_prime3Val", this, other.g4_prime3Val),

g1_prime4Val("g1_prime4Val", this, other.g1_prime4Val),
g2_prime4Val("a2_prime4Val", this, other.g2_prime4Val),
g3_prime4Val("g3_prime4Val", this, other.g3_prime4Val),
g4_prime4Val("g4_prime4Val", this, other.g4_prime4Val),

g1_prime5Val("g1_prime5Val", this, other.g1_prime5Val),
g2_prime5Val("a2_prime5Val", this, other.g2_prime5Val),
g3_prime5Val("g3_prime5Val", this, other.g3_prime5Val),
g4_prime5Val("g4_prime5Val", this, other.g4_prime5Val),

g1_prime6Val("g1_prime6Val", this, other.g1_prime6Val),
g2_prime6Val("a2_prime6Val", this, other.g2_prime6Val),
g3_prime6Val("g3_prime6Val", this, other.g3_prime6Val),
g4_prime6Val("g4_prime6Val", this, other.g4_prime6Val),

g1_prime7Val("g1_prime7Val", this, other.g1_prime7Val),
g2_prime7Val("a2_prime7Val", this, other.g2_prime7Val),
g3_prime7Val("g3_prime7Val", this, other.g3_prime7Val),
g4_prime7Val("g4_prime7Val", this, other.g4_prime7Val),

gzgs1_prime2Val("gzgs1_prime2Val", this, other.gzgs1_prime2Val),
gzgs2Val("gzgs2Val", this, other.gzgs2Val),
gzgs3Val("gzgs3Val", this, other.gzgs3Val),
gzgs4Val("gzgs4Val", this, other.gzgs4Val),
ggsgs2Val("ggsgs2Val", this, other.ggsgs2Val),
ggsgs3Val("ggsgs3Val", this, other.ggsgs3Val),
ggsgs4Val("ggsgs4Val", this, other.ggsgs4Val),


g1ValIm("g1ValIm", this, other.g1ValIm),
g2ValIm("g2ValIm", this, other.g2ValIm),
g3ValIm("g3ValIm", this, other.g3ValIm),
g4ValIm("g4ValIm", this, other.g4ValIm),

g1_primeValIm("g1_primeValIm", this, other.g1_primeValIm),
g2_primeValIm("a2_primeValIm", this, other.g2_primeValIm),
g3_primeValIm("g3_primeValIm", this, other.g3_primeValIm),
g4_primeValIm("g4_primeValIm", this, other.g4_primeValIm),

g1_prime2ValIm("g1_prime2ValIm", this, other.g1_prime2ValIm),
g2_prime2ValIm("a2_prime2ValIm", this, other.g2_prime2ValIm),
g3_prime2ValIm("g3_prime2ValIm", this, other.g3_prime2ValIm),
g4_prime2ValIm("g4_prime2ValIm", this, other.g4_prime2ValIm),

g1_prime3ValIm("g1_prime3ValIm", this, other.g1_prime3ValIm),
g2_prime3ValIm("a2_prime3ValIm", this, other.g2_prime3ValIm),
g3_prime3ValIm("g3_prime3ValIm", this, other.g3_prime3ValIm),
g4_prime3ValIm("g4_prime3ValIm", this, other.g4_prime3ValIm),

g1_prime4ValIm("g1_prime4ValIm", this, other.g1_prime4ValIm),
g2_prime4ValIm("a2_prime4ValIm", this, other.g2_prime4ValIm),
g3_prime4ValIm("g3_prime4ValIm", this, other.g3_prime4ValIm),
g4_prime4ValIm("g4_prime4ValIm", this, other.g4_prime4ValIm),

g1_prime5ValIm("g1_prime5ValIm", this, other.g1_prime5ValIm),
g2_prime5ValIm("a2_prime5ValIm", this, other.g2_prime5ValIm),
g3_prime5ValIm("g3_prime5ValIm", this, other.g3_prime5ValIm),
g4_prime5ValIm("g4_prime5ValIm", this, other.g4_prime5ValIm),

g1_prime6ValIm("g1_prime6ValIm", this, other.g1_prime6ValIm),
g2_prime6ValIm("a2_prime6ValIm", this, other.g2_prime6ValIm),
g3_prime6ValIm("g3_prime6ValIm", this, other.g3_prime6ValIm),
g4_prime6ValIm("g4_prime6ValIm", this, other.g4_prime6ValIm),

g1_prime7ValIm("g1_prime7ValIm", this, other.g1_prime7ValIm),
g2_prime7ValIm("a2_prime7ValIm", this, other.g2_prime7ValIm),
g3_prime7ValIm("g3_prime7ValIm", this, other.g3_prime7ValIm),
g4_prime7ValIm("g4_prime7ValIm", this, other.g4_prime7ValIm),

gzgs1_prime2ValIm("gzgs1_prime2ValIm", this, other.gzgs1_prime2ValIm),
gzgs2ValIm("gzgs2ValIm", this, other.gzgs2ValIm),
gzgs3ValIm("gzgs3ValIm", this, other.gzgs3ValIm),
gzgs4ValIm("gzgs4ValIm", this, other.gzgs4ValIm),
ggsgs2ValIm("ggsgs2ValIm", this, other.ggsgs2ValIm),
ggsgs3ValIm("ggsgs3ValIm", this, other.ggsgs3ValIm),
ggsgs4ValIm("ggsgs4ValIm", this, other.ggsgs4ValIm),

Lambda("Lambda", this, other.Lambda),
Lambda_zgs1("Lambda_zgs1", this, other.Lambda_zgs1),
Lambda_z1("Lambda_z1", this, other.Lambda_z1),
Lambda_z2("Lambda_z2", this, other.Lambda_z2),
Lambda_z3("Lambda_z3", this, other.Lambda_z3),
Lambda_z4("Lambda_z4", this, other.Lambda_z4),
Lambda_Q("Lambda_Q", this, other.Lambda_Q),

Lambda_z11("Lambda_z11", this, other.Lambda_z11),
Lambda_z21("Lambda_z21", this, other.Lambda_z21),
Lambda_z31("Lambda_z31", this, other.Lambda_z31),
Lambda_z41("Lambda_z41", this, other.Lambda_z41),

Lambda_z12("Lambda_z12", this, other.Lambda_z12),
Lambda_z22("Lambda_z22", this, other.Lambda_z22),
Lambda_z32("Lambda_z32", this, other.Lambda_z32),
Lambda_z42("Lambda_z42", this, other.Lambda_z42),

Lambda_z10("Lambda_z10", this, other.Lambda_z10),
Lambda_z20("Lambda_z20", this, other.Lambda_z20),
Lambda_z30("Lambda_z30", this, other.Lambda_z30),
Lambda_z40("Lambda_z40", this, other.Lambda_z40),

cz_q1sq("cz_q1sq", this, other.cz_q1sq),
cz_q2sq("cz_q2sq", this, other.cz_q2sq),
cz_q12sq("cz_q12sq", this, other.cz_q12sq)
{}

void RooSpinZero::calculateAi(Double_t& a1Re, Double_t& a1Im, Double_t& a2Re, Double_t& a2Im, Double_t& a3Re, Double_t& a3Im, bool isGammaV1, bool isGammaV2)const{
  Double_t m1_=m1; if (Vdecay1==0) m1_=0;
  Double_t m2_=m2; if (Vdecay2==0) m2_=0;

  Double_t s = (pow(m12, 2) - pow(m1_, 2) - pow(m2_, 2))/2.;
  if (pow(m1_, 2)>pow(m12, 2) || pow(m2_, 2)>pow(m12, 2)) s = -s;

  Double_t g1_dyn=0;
  Double_t g2_dyn=0;
  Double_t g3_dyn=0;
  Double_t g4_dyn=0;
  Double_t g1_dynIm=0;
  Double_t g2_dynIm=0;
  Double_t g3_dynIm=0;
  Double_t g4_dynIm=0;

  if (!isGammaV1 && !isGammaV2 && !(Vdecay1==0 || Vdecay2==0)){
    if (g1_primeVal!=0) g1_dyn += g1_primeVal * pow(Lambda_z1, 4)/(pow(Lambda_z1, 2) + pow(m1_, 2))/(pow(Lambda_z1, 2) + pow(m2_, 2));
    if (g1_prime2Val!=0) g1_dyn += g1_prime2Val* (pow(m1_, 2) + pow(m2_, 2))/pow(Lambda_z1, 2);
    if (g1_prime3Val!=0) g1_dyn += g1_prime3Val* (pow(m1_, 2) - pow(m2_, 2))/pow(Lambda_z1, 2);
    if (g1_prime4Val!=0) g1_dyn += g1_prime4Val* (m12*m12)/pow(Lambda_Q, 2);
    if (g1_prime5Val!=0) g1_dyn += g1_prime5Val* (pow(m1_, 4) + pow(m2_, 4))/pow(Lambda_z1, 4);
    if (g1_prime6Val!=0) g1_dyn += g1_prime6Val* (pow(m1_, 4) - pow(m2_, 4))/pow(Lambda_z1, 4);
    if (g1_prime7Val!=0) g1_dyn += g1_prime7Val* (pow(m1_, 2) * pow(m2_, 2))/pow(Lambda_z1, 4);

    if (cz_q1sq!=0.) g1_dyn *= 1./(1.+ cz_q1sq*pow(m1_/Lambda_z11, 2));
    if (cz_q2sq!=0.) g1_dyn *= 1./(1.+ cz_q2sq*pow(m2_/Lambda_z12, 2));
    if (cz_q12sq!=0.) g1_dyn *= 1./(1.+ cz_q12sq*pow(m12/Lambda_z10, 2));
    g1_dyn += g1Val;

    g2_dyn = g2Val;
    if (g2_primeVal!=0) g2_dyn += g2_primeVal * pow(Lambda_z2, 4)/(pow(Lambda_z2, 2) + pow(m1_, 2))/(pow(Lambda_z2, 2) + pow(m2_, 2));
    if (g2_prime2Val!=0) g2_dyn += g2_prime2Val* (pow(m1_, 2) + pow(m2_, 2))/pow(Lambda_z2, 2);
    if (g2_prime3Val!=0) g2_dyn += g2_prime3Val* (pow(m1_, 2) - pow(m2_, 2))/pow(Lambda_z2, 2);
    if (g2_prime4Val!=0) g2_dyn += g2_prime4Val* (m12*m12)/pow(Lambda_Q, 2);
    if (g2_prime5Val!=0) g2_dyn += g2_prime5Val* (pow(m1_, 4) + pow(m2_, 4))/pow(Lambda_z2, 4);
    if (g2_prime6Val!=0) g2_dyn += g2_prime6Val* (pow(m1_, 4) - pow(m2_, 4))/pow(Lambda_z2, 4);
    if (g2_prime7Val!=0) g2_dyn += g2_prime7Val* (pow(m1_, 2) * pow(m2_, 2))/pow(Lambda_z2, 4);

    if (cz_q1sq!=0.) g2_dyn *= 1./(1.+ cz_q1sq*pow(m1_/Lambda_z21, 2));
    if (cz_q2sq!=0.) g2_dyn *= 1./(1.+ cz_q2sq*pow(m2_/Lambda_z22, 2));
    if (cz_q12sq!=0.) g2_dyn *= 1./(1.+ cz_q12sq*pow(m12/Lambda_z20, 2));

    g3_dyn = g3Val;
    if (g3_primeVal!=0) g3_dyn += g3_primeVal * pow(Lambda_z3, 4)/(pow(Lambda_z3, 2) + pow(m1_, 2))/(pow(Lambda_z3, 2) + pow(m2_, 2));
    if (g3_prime2Val!=0) g3_dyn += g3_prime2Val* (pow(m1_, 2) + pow(m2_, 2))/pow(Lambda_z3, 2);
    if (g3_prime3Val!=0) g3_dyn += g3_prime3Val* (pow(m1_, 2) - pow(m2_, 2))/pow(Lambda_z3, 2);
    if (g3_prime4Val!=0) g3_dyn += g3_prime4Val* (m12*m12)/pow(Lambda_Q, 2);
    if (g3_prime5Val!=0) g3_dyn += g3_prime5Val* (pow(m1_, 4) + pow(m2_, 4))/pow(Lambda_z3, 4);
    if (g3_prime6Val!=0) g3_dyn += g3_prime6Val* (pow(m1_, 4) - pow(m2_, 4))/pow(Lambda_z3, 4);
    if (g3_prime7Val!=0) g3_dyn += g3_prime7Val* (pow(m1_, 2) * pow(m2_, 2))/pow(Lambda_z3, 4);

    if (cz_q1sq!=0.) g3_dyn *= 1./(1.+ cz_q1sq*pow(m1_/Lambda_z31, 2));
    if (cz_q2sq!=0.) g3_dyn *= 1./(1.+ cz_q2sq*pow(m2_/Lambda_z32, 2));
    if (cz_q12sq!=0.) g3_dyn *= 1./(1.+ cz_q12sq*pow(m12/Lambda_z30, 2));


    g4_dyn = g4Val;
    if (g4_primeVal!=0) g4_dyn += g4_primeVal * pow(Lambda_z4, 4)/(pow(Lambda_z4, 2) + pow(m1_, 2))/(pow(Lambda_z4, 2) + pow(m2_, 2));
    if (g4_prime2Val!=0) g4_dyn += g4_prime2Val* (pow(m1_, 2) + pow(m2_, 2))/pow(Lambda_z4, 2);
    if (g4_prime3Val!=0) g4_dyn += g4_prime3Val* (pow(m1_, 2) - pow(m2_, 2))/pow(Lambda_z4, 2);
    if (g4_prime4Val!=0) g4_dyn += g4_prime4Val* (m12*m12)/pow(Lambda_Q, 2);
    if (g4_prime5Val!=0) g4_dyn += g4_prime5Val* (pow(m1_, 4) + pow(m2_, 4))/pow(Lambda_z4, 4);
    if (g4_prime6Val!=0) g4_dyn += g4_prime6Val* (pow(m1_, 4) - pow(m2_, 4))/pow(Lambda_z4, 4);
    if (g4_prime7Val!=0) g4_dyn += g4_prime7Val* (pow(m1_, 2) * pow(m2_, 2))/pow(Lambda_z4, 4);

    if (cz_q1sq!=0.) g4_dyn *= 1./(1.+ cz_q1sq*pow(m1_/Lambda_z41, 2));
    if (cz_q2sq!=0.) g4_dyn *= 1./(1.+ cz_q2sq*pow(m2_/Lambda_z42, 2));
    if (cz_q12sq!=0.) g4_dyn *= 1./(1.+ cz_q12sq*pow(m12/Lambda_z40, 2));

    g1_dynIm = 0;
    if (g1_primeValIm!=0) g1_dynIm += g1_primeValIm * pow(Lambda_z1, 4)/(pow(Lambda_z1, 2) + pow(m1_, 2))/(pow(Lambda_z1, 2) + pow(m2_, 2));
    if (g1_prime2ValIm!=0) g1_dynIm += g1_prime2ValIm* (pow(m1_, 2) + pow(m2_, 2))/pow(Lambda_z1, 2);
    if (g1_prime3ValIm!=0) g1_dynIm += g1_prime3ValIm* (pow(m1_, 2) - pow(m2_, 2))/pow(Lambda_z1, 2);
    if (g1_prime4ValIm!=0) g1_dynIm += g1_prime4ValIm* (m12*m12)/pow(Lambda_Q, 2);
    if (g1_prime5ValIm!=0) g1_dynIm += g1_prime5ValIm* (pow(m1_, 4) + pow(m2_, 4))/pow(Lambda_z1, 4);
    if (g1_prime6ValIm!=0) g1_dynIm += g1_prime6ValIm* (pow(m1_, 4) - pow(m2_, 4))/pow(Lambda_z1, 4);
    if (g1_prime7ValIm!=0) g1_dynIm += g1_prime7ValIm* (pow(m1_, 2) * pow(m2_, 2))/pow(Lambda_z1, 4);

    if (cz_q1sq!=0.) g1_dynIm *= 1./(1.+ cz_q1sq*pow(m1_/Lambda_z11, 2));
    if (cz_q2sq!=0.) g1_dynIm *= 1./(1.+ cz_q2sq*pow(m2_/Lambda_z12, 2));
    if (cz_q12sq!=0.) g1_dynIm *= 1./(1.+ cz_q12sq*pow(m12/Lambda_z10, 2));
    g1_dynIm += g1ValIm;

    g2_dynIm = g2ValIm;
    if (g2_primeValIm!=0) g2_dynIm += g2_primeValIm * pow(Lambda_z2, 4)/(pow(Lambda_z2, 2) + pow(m1_, 2))/(pow(Lambda_z2, 2) + pow(m2_, 2));
    if (g2_prime2ValIm!=0) g2_dynIm += g2_prime2ValIm* (pow(m1_, 2) + pow(m2_, 2))/pow(Lambda_z2, 2);
    if (g2_prime3ValIm!=0) g2_dynIm += g2_prime3ValIm* (pow(m1_, 2) - pow(m2_, 2))/pow(Lambda_z2, 2);
    if (g2_prime4ValIm!=0) g2_dynIm += g2_prime4ValIm* (m12*m12)/pow(Lambda_Q, 2);
    if (g2_prime5ValIm!=0) g2_dynIm += g2_prime5ValIm* (pow(m1_, 4) + pow(m2_, 4))/pow(Lambda_z2, 4);
    if (g2_prime6ValIm!=0) g2_dynIm += g2_prime6ValIm* (pow(m1_, 4) - pow(m2_, 4))/pow(Lambda_z2, 4);
    if (g2_prime7ValIm!=0) g2_dynIm += g2_prime7ValIm* (pow(m1_, 2) * pow(m2_, 2))/pow(Lambda_z2, 4);

    if (cz_q1sq!=0.) g2_dynIm *= 1./(1.+ cz_q1sq*pow(m1_/Lambda_z21, 2));
    if (cz_q2sq!=0.) g2_dynIm *= 1./(1.+ cz_q2sq*pow(m2_/Lambda_z22, 2));
    if (cz_q12sq!=0.) g2_dynIm *= 1./(1.+ cz_q12sq*pow(m12/Lambda_z20, 2));

    g3_dynIm = g3ValIm;
    if (g3_primeValIm!=0) g3_dynIm += g3_primeValIm * pow(Lambda_z3, 4)/(pow(Lambda_z3, 2) + pow(m1_, 2))/(pow(Lambda_z3, 2) + pow(m2_, 2));
    if (g3_prime2ValIm!=0) g3_dynIm += g3_prime2ValIm* (pow(m1_, 2) + pow(m2_, 2))/pow(Lambda_z3, 2);
    if (g3_prime3ValIm!=0) g3_dynIm += g3_prime3ValIm* (pow(m1_, 2) - pow(m2_, 2))/pow(Lambda_z3, 2);
    if (g3_prime4ValIm!=0) g3_dynIm += g3_prime4ValIm* (m12*m12)/pow(Lambda_Q, 2);
    if (g3_prime5ValIm!=0) g3_dynIm += g3_prime5ValIm* (pow(m1_, 4) + pow(m2_, 4))/pow(Lambda_z3, 4);
    if (g3_prime6ValIm!=0) g3_dynIm += g3_prime6ValIm* (pow(m1_, 4) - pow(m2_, 4))/pow(Lambda_z3, 4);
    if (g3_prime7ValIm!=0) g3_dynIm += g3_prime7ValIm* (pow(m1_, 2) * pow(m2_, 2))/pow(Lambda_z3, 4);

    if (cz_q1sq!=0.) g3_dynIm *= 1./(1.+ cz_q1sq*pow(m1_/Lambda_z31, 2));
    if (cz_q2sq!=0.) g3_dynIm *= 1./(1.+ cz_q2sq*pow(m2_/Lambda_z32, 2));
    if (cz_q12sq!=0.) g3_dynIm *= 1./(1.+ cz_q12sq*pow(m12/Lambda_z30, 2));


    g4_dynIm = g4ValIm;
    if (g4_primeValIm!=0) g4_dynIm += g4_primeValIm * pow(Lambda_z4, 4)/(pow(Lambda_z4, 2) + pow(m1_, 2))/(pow(Lambda_z4, 2) + pow(m2_, 2));
    if (g4_prime2ValIm!=0) g4_dynIm += g4_prime2ValIm* (pow(m1_, 2) + pow(m2_, 2))/pow(Lambda_z4, 2);
    if (g4_prime3ValIm!=0) g4_dynIm += g4_prime3ValIm* (pow(m1_, 2) - pow(m2_, 2))/pow(Lambda_z4, 2);
    if (g4_prime4ValIm!=0) g4_dynIm += g4_prime4ValIm* (m12*m12)/pow(Lambda_Q, 2);
    if (g4_prime5ValIm!=0) g4_dynIm += g4_prime5ValIm* (pow(m1_, 4) + pow(m2_, 4))/pow(Lambda_z4, 4);
    if (g4_prime6ValIm!=0) g4_dynIm += g4_prime6ValIm* (pow(m1_, 4) - pow(m2_, 4))/pow(Lambda_z4, 4);
    if (g4_prime7ValIm!=0) g4_dynIm += g4_prime7ValIm* (pow(m1_, 2) * pow(m2_, 2))/pow(Lambda_z4, 4);

    if (cz_q1sq!=0.) g4_dynIm *= 1./(1.+ cz_q1sq*pow(m1_/Lambda_z41, 2));
    if (cz_q2sq!=0.) g4_dynIm *= 1./(1.+ cz_q2sq*pow(m2_/Lambda_z42, 2));
    if (cz_q12sq!=0.) g4_dynIm *= 1./(1.+ cz_q12sq*pow(m12/Lambda_z40, 2));
  }
  else if ((!isGammaV1 || !isGammaV2) && !(Vdecay1==0 && Vdecay2==0)){
    if (gzgs1_prime2Val!=0 && !isGammaV1) g1_dyn += gzgs1_prime2Val* pow(m1_, 2)/pow(Lambda_zgs1, 2);
    if (gzgs1_prime2Val!=0 && !isGammaV2) g1_dyn += gzgs1_prime2Val* pow(m2_, 2)/pow(Lambda_zgs1, 2);
    g2_dyn = gzgs2Val;
    g3_dyn = gzgs3Val;
    g4_dyn = gzgs4Val;

    if (gzgs1_prime2ValIm!=0 && !isGammaV1) g1_dynIm += gzgs1_prime2ValIm* pow(m1_, 2)/pow(Lambda_zgs1, 2);
    if (gzgs1_prime2ValIm!=0 && !isGammaV2) g1_dynIm += gzgs1_prime2ValIm* pow(m2_, 2)/pow(Lambda_zgs1, 2);
    g2_dynIm = gzgs2ValIm;
    g3_dynIm = gzgs3ValIm;
    g4_dynIm = gzgs4ValIm;
  }
  else{
    g2_dyn = ggsgs2Val;
    g3_dyn = ggsgs3Val;
    g4_dyn = ggsgs4Val;

    g2_dynIm = ggsgs2ValIm;
    g3_dynIm = ggsgs3ValIm;
    g4_dynIm = ggsgs4ValIm;
  }
  /*
  cout << "g1: " << g1_dyn << " " << g1_dynIm << '\t';
  cout << "g2: " << g2_dyn << " " << g2_dynIm << '\t';
  cout << "g3: " << g3_dyn << " " << g3_dynIm << '\t';
  cout << "g4: " << g4_dyn << " " << g4_dynIm << '\t';
  cout << endl;
  */

  Double_t kappa = s/pow(Lambda,2);
  a3Re = -2.*g4_dyn*pow(m12, 2);
  a2Re = -(2.*g2_dyn + g3_dyn*kappa)*pow(m12, 2);
  a1Re = g1_dyn*pow(mV, 2) - a2Re*s/pow(m12, 2);
  a3Im = -2.*g4_dynIm*pow(m12, 2);
  a2Im = -(2.*g2_dynIm + g3_dynIm*kappa)*pow(m12, 2);
  a1Im = g1_dynIm*pow(mV, 2) - a2Im*s/pow(m12, 2);
}
void RooSpinZero::calculatePropagator(Double_t& propRe, Double_t& propIm, Double_t mass, bool useGamma)const{
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
void RooSpinZero::calculateAmplitudeScale(bool isGammaV1, bool isGammaV2)const{

}
void RooSpinZero::calculateAmplitudes(Double_t& A00Re, Double_t& A00Im, Double_t& AppRe, Double_t& AppIm, Double_t& AmmRe, Double_t& AmmIm, bool isGammaV1, bool isGammaV2)const{
  Double_t m1_=m1; if (Vdecay1==0) m1_=0;
  Double_t m2_=m2; if (Vdecay2==0) m2_=0;

  Double_t a1Re, a2Re, a3Re, a1Im, a2Im, a3Im;
  calculateAi(a1Re, a1Im, a2Re, a2Im, a3Re, a3Im, isGammaV1, isGammaV2);

  Double_t propV1Re=0, propV2Re=0;
  Double_t propV1Im=-1, propV2Im=-1;
  if (Vdecay1!=0) calculatePropagator(propV1Re, propV1Im, m1_, isGammaV1);
  if (Vdecay2!=0) calculatePropagator(propV2Re, propV2Im, m2_, isGammaV2);

  Double_t eta1 = m1_ / m12;
  Double_t eta2 = m2_ / m12;
  Double_t eta1p2 = 1.;
  if (Vdecay1!=0) eta1p2 *= eta1;
  if (Vdecay2!=0) eta1p2 *= eta2;

  Double_t etas = (1. - pow(eta1, 2) - pow(eta2, 2))/2.;
  if (pow(eta1+eta2, 2)>1.) etas = -etas;
  Double_t etasp = pow(etas, 2) - pow(eta1*eta2, 2); // Notice how eta1p2 is not used. The second set of multiplications below is the reason to it!
  if (etasp<0) etasp=0;

  Double_t A00Re_tmp, A00Im_tmp, AppRe_tmp, AppIm_tmp, AmmRe_tmp, AmmIm_tmp;
  A00Re_tmp = -(a1Re*etas + a2Re*etasp);
  A00Im_tmp = -(a1Im*etas + a2Im*etasp);
  AppRe_tmp = (a1Re - a3Im*sqrt(etasp));
  AppIm_tmp = (a1Im + a3Re*sqrt(etasp));
  AmmRe_tmp = (a1Re + a3Im*sqrt(etasp));
  AmmIm_tmp = (a1Im - a3Re*sqrt(etasp));

  //A00Re /= eta1p2;
  //A00Im /= eta1p2;
  AppRe_tmp *= eta1p2;
  AppIm_tmp *= eta1p2;
  AmmRe_tmp *= eta1p2;
  AmmIm_tmp *= eta1p2;
  
  // A_old = ARe_old+i*AIm_old => A_new = ARe_new + i*AIm_new = A_old*propV1*propV2
  A00Re = ((A00Re_tmp*propV1Re - A00Im_tmp*propV1Im)*propV2Re - (A00Re_tmp*propV1Im + A00Im_tmp*propV1Re)*propV2Im)/2.;
  A00Im = ((A00Re_tmp*propV1Re - A00Im_tmp*propV1Im)*propV2Im + (A00Re_tmp*propV1Im + A00Im_tmp*propV1Re)*propV2Re)/2.;
  AppRe = ((AppRe_tmp*propV1Re - AppIm_tmp*propV1Im)*propV2Re - (AppRe_tmp*propV1Im + AppIm_tmp*propV1Re)*propV2Im)/4.;
  AppIm = ((AppRe_tmp*propV1Re - AppIm_tmp*propV1Im)*propV2Im + (AppRe_tmp*propV1Im + AppIm_tmp*propV1Re)*propV2Re)/4.;
  AmmRe = ((AmmRe_tmp*propV1Re - AmmIm_tmp*propV1Im)*propV2Re - (AmmRe_tmp*propV1Im + AmmIm_tmp*propV1Re)*propV2Im)/4.;
  AmmIm = ((AmmRe_tmp*propV1Re - AmmIm_tmp*propV1Im)*propV2Im + (AmmRe_tmp*propV1Im + AmmIm_tmp*propV1Re)*propV2Re)/4.;

  /*
  cout << "A00 = " << A00Re << ", " << A00Im << endl;
  cout << "App = " << AppRe << ", " << AppIm << endl;
  cout << "Amm = " << AmmRe << ", " << AmmIm << endl;
  */
  if (!(A00Re==A00Re) || !(AppRe==AppRe) || !(AmmRe==AmmRe) || !(A00Im==A00Im) || !(AppIm==AppIm) || !(AmmIm==AmmIm)){
    cout << eta1 << '\t' << eta2 << '\t' << etas << '\t' << etasp << '\t' << eta1p2 << endl;
  }
}

void RooSpinZero::setProxies(modelMeasurables _measurables){
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
void RooSpinZero::setProxy(RooRealProxy& proxy, RooAbsReal* objectPtr){
  if (objectPtr!=0) proxy.setArg((RooAbsReal&)*objectPtr);
}
