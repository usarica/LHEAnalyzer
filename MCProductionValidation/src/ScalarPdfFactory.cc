#include "../include/ScalarPdfFactory.h"

ScalarPdfFactory::ScalarPdfFactory(RooSpinZero::modelMeasurables measurables_, bool acceptance_, int V1decay_, int V2decay_) :
parameterization(0),
pmf_applied(false),
acceptance(acceptance_),
V1decay(V1decay_),
V2decay(V2decay_)
{
  initMeasurables(measurables_);
  initMassPole();
  initVdecayParams();
  initGVals();
}
ScalarPdfFactory::ScalarPdfFactory(RooSpinZero::modelMeasurables measurables_, double gRatio_[4][8], bool pmf_applied_, bool acceptance_, int V1decay_, int V2decay_) :
parameterization(1),
pmf_applied(pmf_applied_),
acceptance(acceptance_),
V1decay(V1decay_),
V2decay(V2decay_)
{
  for (int v=0; v<4; v++){
    for (int k=0; k<8; k++) gRatio[v][k] = gRatio_[v][k];
  }
  initMeasurables(measurables_);
  initMassPole();
  initVdecayParams();
  initGVals();
}

ScalarPdfFactory::~ScalarPdfFactory(){
  destroyGVals();
  destroyVdecayParams();
  destroyMassPole();
}
void ScalarPdfFactory::initMeasurables(RooSpinZero::modelMeasurables measurables_){
  measurables.h1 = measurables_.h1;
  measurables.h2 = measurables_.h2;
  measurables.Phi = measurables_.Phi;
  measurables.m1 = measurables_.m1;
  measurables.m2 = measurables_.m2;
  measurables.m12 = measurables_.m12;
  measurables.hs = measurables_.hs;
  measurables.Phi1 = measurables_.Phi1;
  measurables.Y = measurables_.Y;
}
void ScalarPdfFactory::initMassPole(){
  parameters.mX = new RooRealVar("mX", "mX", (measurables.m12)->getVal());
  parameters.gamX = new RooRealVar("gamX", "gamX", 4.07e-3);
}
void ScalarPdfFactory::initVdecayParams(){
  if ((V1decay>0 && V2decay<0) || (V1decay<0 && V2decay>0)) cerr << "ScalarPdfFactory::initVdecayParams: V1 and V2 decays are inconsistent!" << endl;
  if (V1decay>0){
    parameters.mV = new RooRealVar("mV", "mV", 91.1876);
    parameters.gamV = new RooRealVar("gamV", "gamV", 2.4952);
  }
  else if (V1decay<0){
    parameters.mV = new RooRealVar("mV", "mV", 80.399);
    parameters.gamV = new RooRealVar("gamV", "gamV", 2.085);
  }
  else{
    parameters.mV = new RooRealVar("mV", "mV", 0);
    parameters.gamV = new RooRealVar("gamV", "gamV", 0);
  }

  /*
  For Zs,
  gV = T3 - 2 Q sin**2 thetaW
  gA = T3
  T3 = 1/2 for nu, u; -1/2 for l, d
  Q = 0, -1, 2/3, -1/3 for nu, l, u, d
  Values below are for sin**2 thetaW = 0.23119
  For Ws, gV=gA=1.
  R1, R2 = 2 gV gA / (gV**2 + gA**2) for Z1, Z2
  */

  double atomicT3 = 0.5;
  double atomicCharge = 1.;
  double sin2t = 0.23119;

  double gV_up = atomicT3 - 2.*(2.*atomicCharge/3.)*sin2t;
  double gV_dn = -atomicT3 - 2.*(-atomicCharge/3.)*sin2t;
  double gV_l = -atomicT3 - 2.*(-atomicCharge)*sin2t;
  double gV_nu = atomicT3;

  double gA_up = atomicT3;
  double gA_dn = -atomicT3;
  double gA_l = -atomicT3;
  double gA_nu = atomicT3;

  switch (V1decay){
  case 3:
    parameters.R1Val = new RooRealVar("R1Val", "R1Val", 0.5*((2.*gV_up*gA_up)/(pow(gV_up, 2)+pow(gA_up, 2))+(2.*gV_dn*gA_dn)/(pow(gV_dn, 2)+pow(gA_dn, 2)))); // Z->uu+dd avg.
    break;
  case 5:
    parameters.R1Val = new RooRealVar("R1Val", "R1Val", (2.*gV_dn*gA_dn)/(pow(gV_dn, 2)+pow(gA_dn, 2))); // Z->dd
    break;
  case 4:
    parameters.R1Val = new RooRealVar("R1Val", "R1Val", (2.*gV_up*gA_up)/(pow(gV_up, 2)+pow(gA_up, 2))); // Z->uu
    break;
  case 2:
    parameters.R1Val = new RooRealVar("R1Val", "R1Val", (2.*gV_nu*gA_nu)/(pow(gV_nu, 2)+pow(gA_nu, 2))); // Z->nunu
    break;
  case 1:
    parameters.R1Val = new RooRealVar("R1Val", "R1Val", (2.*gV_l*gA_l)/(pow(gV_l, 2)+pow(gA_l, 2))); // Z->ll
    break;
  case -1:
    parameters.R1Val = new RooRealVar("R1Val", "R1Val", 1); // W
    break;
  default:
    parameters.R1Val = new RooRealVar("R1Val", "R1Val", 0); // gamma
    break;
  }
  switch (V2decay){
  case 3:
    parameters.R2Val = new RooRealVar("R2Val", "R2Val", 0.5*((2.*gV_up*gA_up)/(pow(gV_up, 2)+pow(gA_up, 2))+(2.*gV_dn*gA_dn)/(pow(gV_dn, 2)+pow(gA_dn, 2)))); // Z->uu+dd avg.
    break;
  case 5:
    parameters.R2Val = new RooRealVar("R2Val", "R2Val", (2.*gV_dn*gA_dn)/(pow(gV_dn, 2)+pow(gA_dn, 2))); // Z->dd
    break;
  case 4:
    parameters.R2Val = new RooRealVar("R2Val", "R2Val", (2.*gV_up*gA_up)/(pow(gV_up, 2)+pow(gA_up, 2))); // Z->uu
    break;
  case 2:
    parameters.R2Val = new RooRealVar("R2Val", "R2Val", (2.*gV_nu*gA_nu)/(pow(gV_nu, 2)+pow(gA_nu, 2))); // Z->nunu
    break;
  case 1:
    parameters.R2Val = new RooRealVar("R2Val", "R2Val", (2.*gV_l*gA_l)/(pow(gV_l, 2)+pow(gA_l, 2))); // Z->ll
    break;
  case -1:
    parameters.R2Val = new RooRealVar("R2Val", "R2Val", 1); // W
    break;
  default:
    parameters.R2Val = new RooRealVar("R2Val", "R2Val", 0); // gamma
    break;
  }

  parameters.mV->removeRange();
  parameters.gamV->removeRange();
  parameters.R1Val->removeRange();
  parameters.R2Val->removeRange();
}
void ScalarPdfFactory::resetVdecay(int V1decay_, int V2decay_){
  if ((V1decay_>0 && V2decay_<0) || (V1decay_<0 && V2decay_>0)){ cerr << "ScalarPdfFactory::resetVdecayParams: V1 and V2 decays are inconsistent!" << endl; return; }

  V1decay=V1decay_;
  V2decay=V2decay_;

  bool is_mvconst = parameters.mV->isConstant();
  bool is_gamvconst = parameters.gamV->isConstant();
  bool is_r1const = parameters.R1Val->isConstant();
  bool is_r2const = parameters.R2Val->isConstant();
  parameters.mV->setConstant(false);
  parameters.gamV->setConstant(false);
  parameters.R1Val->setConstant(false);
  parameters.R2Val->setConstant(false);

  if (V1decay_>0){
    parameters.mV->setVal(91.1876);
    parameters.gamV->setVal(2.4952);
  }
  else if (V1decay_<0){
    parameters.mV->setVal(80.399);
    parameters.gamV->setVal(2.085);
  }
  else{
    parameters.mV->setVal(0);
    parameters.gamV->setVal(0);
  }

  double atomicT3 = 0.5;
  double atomicCharge = 1.;
  double sin2t = 0.23119;

  double gV_up = atomicT3 - 2.*(2.*atomicCharge/3.)*sin2t;
  double gV_dn = -atomicT3 - 2.*(-atomicCharge/3.)*sin2t;
  double gV_l = -atomicT3 - 2.*(-atomicCharge)*sin2t;
  double gV_nu = atomicT3;

  double gA_up = atomicT3;
  double gA_dn = -atomicT3;
  double gA_l = -atomicT3;
  double gA_nu = atomicT3;

  switch (V1decay_){
  case 3:
    parameters.R1Val->setVal(0.5*((2.*gV_up*gA_up)/(pow(gV_up, 2)+pow(gA_up, 2))+(2.*gV_dn*gA_dn)/(pow(gV_dn, 2)+pow(gA_dn, 2)))); // Z->uu+dd avg.
    break;
  case 5:
    parameters.R1Val->setVal((2.*gV_dn*gA_dn)/(pow(gV_dn, 2)+pow(gA_dn, 2))); // Z->dd
    break;
  case 4:
    parameters.R1Val->setVal((2.*gV_up*gA_up)/(pow(gV_up, 2)+pow(gA_up, 2))); // Z->uu
    break;
  case 2:
    parameters.R1Val->setVal((2.*gV_nu*gA_nu)/(pow(gV_nu, 2)+pow(gA_nu, 2))); // Z->nunu
    break;
  case 1:
    parameters.R1Val->setVal((2.*gV_l*gA_l)/(pow(gV_l, 2)+pow(gA_l, 2))); // Z->ll
    break;
  case -1:
    parameters.R1Val->setVal(1); // W
    break;
  default:
    parameters.R1Val->setVal(0); // gamma
    break;
  }
  switch (V2decay_){
  case 3:
    parameters.R2Val->setVal(0.5*((2.*gV_up*gA_up)/(pow(gV_up, 2)+pow(gA_up, 2))+(2.*gV_dn*gA_dn)/(pow(gV_dn, 2)+pow(gA_dn, 2)))); // Z->uu+dd avg.
    break;
  case 5:
    parameters.R2Val->setVal((2.*gV_dn*gA_dn)/(pow(gV_dn, 2)+pow(gA_dn, 2))); // Z->dd
    break;
  case 4:
    parameters.R2Val->setVal((2.*gV_up*gA_up)/(pow(gV_up, 2)+pow(gA_up, 2))); // Z->uu
    break;
  case 2:
    parameters.R2Val->setVal((2.*gV_nu*gA_nu)/(pow(gV_nu, 2)+pow(gA_nu, 2))); // Z->nunu
    break;
  case 1:
    parameters.R2Val->setVal((2.*gV_l*gA_l)/(pow(gV_l, 2)+pow(gA_l, 2))); // Z->ll
    break;
  case -1:
    parameters.R2Val->setVal(1); // W
    break;
  default:
    parameters.R2Val->setVal(0); // gamma
    break;
  }

  parameters.mV->setConstant(is_mvconst);
  parameters.gamV->setConstant(is_gamvconst);
  parameters.R1Val->setConstant(is_r1const);
  parameters.R2Val->setConstant(is_r2const);
}
void ScalarPdfFactory::initFractionsPhases(){
  for (int v=0; v<4; v++){
    for (int k=0; k<8; k++) gRatioVal[v][k] = new RooRealVar(Form("g%i_%iMixVal", v+1, k), Form("g%i_%iMixVal", v+1, k), gRatio[v][k]);
  }

  RooArgList gBareFracList;
  for (int v=0; v<8; v++){
    TString strcore;
    TString strapp;
    if (v>1) strapp.Prepend(Form("%i", v));
    if (v>0) strapp.Prepend("_prime");
    TString strappPhase = strapp;
    strapp.Append("Frac");
    strappPhase.Append("Phase");

    double firstFracVal = 0;
    double phaseBound = TMath::Pi();
    if (pmf_applied){
      firstFracVal = -1;
      phaseBound = 0;
    }

    if (v>0){
      strcore = "g1";
      strcore.Append(strapp);
      if (gRatio[0][v]!=0) g1Frac[v-1] = new RooRealVar(strcore, strcore, 0, firstFracVal, 1);
      else g1Frac[v-1] = new RooRealVar(strcore, strcore, 0);
      gBareFracList.add(*(g1Frac[v-1]));
    }

    strcore = "g2";
    strcore.Append(strapp);
    if (gRatio[1][v]!=0) g2Frac[v] = new RooRealVar(strcore, strcore, 0, firstFracVal, 1);
    else g2Frac[v] = new RooRealVar(strcore, strcore, 0);
    gBareFracList.add(*(g2Frac[v]));

    strcore = "g3";
    strcore.Append(strapp);
    if (gRatio[2][v]!=0) g3Frac[v] = new RooRealVar(strcore, strcore, 0, firstFracVal, 1);
    else g3Frac[v] = new RooRealVar(strcore, strcore, 0);
    gBareFracList.add(*(g3Frac[v]));

    strcore = "g4";
    strcore.Append(strapp);
    if (gRatio[3][v]!=0) g4Frac[v] = new RooRealVar(strcore, strcore, 0, firstFracVal, 1);
    else g4Frac[v] = new RooRealVar(strcore, strcore, 0);
    gBareFracList.add(*(g4Frac[v]));


    if (v>0){
      strcore = "g1";
      strcore.Append(strappPhase);
      if (gRatio[0][v]!=0) g1Phase[v-1] = new RooRealVar(strcore, strcore, 0, -phaseBound, phaseBound);
      else g1Phase[v-1] = new RooRealVar(strcore, strcore, 0);
    }

    strcore = "g2";
    strcore.Append(strappPhase);
    if (gRatio[1][v]!=0) g2Phase[v] = new RooRealVar(strcore, strcore, 0, -phaseBound, phaseBound);
    else g2Phase[v] = new RooRealVar(strcore, strcore, 0);

    strcore = "g3";
    strcore.Append(strappPhase);
    if (gRatio[2][v]!=0) g3Phase[v] = new RooRealVar(strcore, strcore, 0, -phaseBound, phaseBound);
    else g3Phase[v] = new RooRealVar(strcore, strcore, 0);

    strcore = "g4";
    strcore.Append(strappPhase);
    if (gRatio[3][v]!=0) g4Phase[v] = new RooRealVar(strcore, strcore, 0, -phaseBound, phaseBound);
    else g4Phase[v] = new RooRealVar(strcore, strcore, 0);
  }
  TString sumFormula = "(";
  for (int gg=0; gg<gBareFracList.getSize(); gg++){
    TString bareFormula = "abs(@";
    bareFormula.Append(Form("%i", gg));
    bareFormula.Append(")");
    if (gg!=(gBareFracList.getSize()-1)) bareFormula.Append("+");
    sumFormula.Append(bareFormula);
  }
  sumFormula.Append(")");
  TString gFracSumName = "sum_gjAbsFrac";
  gFracSum = new RooFormulaVar(gFracSumName, gFracSumName, sumFormula.Data(), gBareFracList);
  for (int v=0; v<8; v++){
    TString strcore;
    TString strapp;
    if (v>1) strapp.Prepend(Form("%i", v));
    if (v>0) strapp.Prepend("_prime");
    strapp.Append("InterpFrac");

    if (v>0){
      strcore = "g1";
      strcore.Append(strapp);
      RooArgList tmpArg1(*(g1Frac[v-1]), *gFracSum);
      g1FracInterp[v] = new RooFormulaVar(strcore, strcore, "(@1>1 ? 0. : @0)", tmpArg1);
    }
    else{
      strcore = "g1";
      strcore.Append(strapp);
      RooArgList tmpGArg(*gFracSum);
      g1FracInterp[v] = new RooFormulaVar(strcore, strcore, "(@0>1 ? 0. : abs(1.-@0))", tmpGArg);
    }

    strcore = "g2";
    strcore.Append(strapp);
    RooArgList tmpArg2(*(g2Frac[v]), *gFracSum);
    g2FracInterp[v] = new RooFormulaVar(strcore, strcore, "(@1>1 ? 0. : @0)", tmpArg2);

    strcore = "g3";
    strcore.Append(strapp);
    RooArgList tmpArg3(*(g3Frac[v]), *gFracSum);
    g3FracInterp[v] = new RooFormulaVar(strcore, strcore, "(@1>1 ? 0. : @0)", tmpArg3);

    strcore = "g4";
    strcore.Append(strapp);
    RooArgList tmpArg4(*(g4Frac[v]), *gFracSum);
    g4FracInterp[v] = new RooFormulaVar(strcore, strcore, "(@1>1 ? 0. : @0)", tmpArg4);
  }
  for (int v=0; v<8; v++){
    for (int im=0; im<2; im++){
      TString strcore;
      TString strapp = "Val";
      if (im==1) strapp.Append("Im");
      if (v>1) strapp.Prepend(Form("%i", v));
      if (v>0) strapp.Prepend("_prime");

      strcore = "g1";
      strcore.Append(strapp);
      RooArgList tmpArg1(*(gRatioVal[0][v]), *(g1FracInterp[v]));
      TString strFormulag1 = "@0*sqrt(@1)";
      if (v>0){
        tmpArg1.add(*(g1Phase[v-1]));
        if (im==0) strFormulag1.Append("*cos(@2)");
        else strFormulag1.Append("*sin(@2)");
      }
      RooFormulaVar* g1Val = new RooFormulaVar(strcore, strcore, strFormulag1, tmpArg1);
      parameters.g1List[v][im] = (RooAbsReal*)g1Val;

      strcore = "g2";
      strcore.Append(strapp);
      RooArgList tmpArg2(*(gRatioVal[1][v]), *(g2FracInterp[v]), *(g2Phase[v]));
      TString strFormulag2 = "@0*sqrt(@1)";
      if (im==0) strFormulag2.Append("*cos(@2)");
      else strFormulag2.Append("*sin(@2)");
      RooFormulaVar* g2Val = new RooFormulaVar(strcore, strcore, strFormulag2, tmpArg2);
      parameters.g2List[v][im] = (RooAbsReal*)g2Val;

      strcore = "g3";
      strcore.Append(strapp);
      RooArgList tmpArg3(*(gRatioVal[2][v]), *(g3FracInterp[v]), *(g3Phase[v]));
      TString strFormulag3 = "@0*sqrt(@1)";
      if (im==0) strFormulag3.Append("*cos(@2)");
      else strFormulag3.Append("*sin(@2)");
      RooFormulaVar* g3Val = new RooFormulaVar(strcore, strcore, strFormulag3, tmpArg3);
      parameters.g3List[v][im] = (RooAbsReal*)g3Val;

      strcore = "g4";
      strcore.Append(strapp);
      RooArgList tmpArg4(*(gRatioVal[3][v]), *(g4FracInterp[v]), *(g4Phase[v]));
      TString strFormulag4 = "@0*sqrt(@1)";
      if (im==0) strFormulag4.Append("*cos(@2)");
      else strFormulag4.Append("*sin(@2)");
      RooFormulaVar* g4Val = new RooFormulaVar(strcore, strcore, strFormulag4, tmpArg4);
      parameters.g4List[v][im] = (RooAbsReal*)g4Val;
    }
  }
}
void ScalarPdfFactory::initGVals(){
  parameters.Lambda = new RooRealVar("Lambda", "Lambda", 1000.);
  parameters.Lambda_z1 = new RooRealVar("Lambda_z1", "Lambda_z1", 10000.);
  parameters.Lambda_z2 = new RooRealVar("Lambda_z2", "Lambda_z2", 10000.);
  parameters.Lambda_z3 = new RooRealVar("Lambda_z3", "Lambda_z3", 10000.);
  parameters.Lambda_z4 = new RooRealVar("Lambda_z4", "Lambda_z4", 10000.);
  parameters.Lambda_Q = new RooRealVar("Lambda_Q", "Lambda_Q", 10000.);

  for (int v=0; v<3; v++){
    TString strcore;
    double initval = 100.;

    strcore = "Lambda_z1";
    strcore.Append(Form("%i", (v!=3 ? v+1 : 0)));
    parameters.Lambda_z1qsq[v] = new RooRealVar(strcore, strcore, initval, 0., 1e15);
    parameters.Lambda_z1qsq[v]->removeMax();
    strcore = "Lambda_z2";
    strcore.Append(Form("%i", (v!=3 ? v+1 : 0)));
    parameters.Lambda_z2qsq[v] = new RooRealVar(strcore, strcore, initval, 0., 1e15);
    parameters.Lambda_z2qsq[v]->removeMax();
    strcore = "Lambda_z3";
    strcore.Append(Form("%i", (v!=3 ? v+1 : 0)));
    parameters.Lambda_z3qsq[v] = new RooRealVar(strcore, strcore, initval, 0., 1e15);
    parameters.Lambda_z3qsq[v]->removeMax();
    strcore = "Lambda_z4";
    strcore.Append(Form("%i", (v!=3 ? v+1 : 0)));
    parameters.Lambda_z4qsq[v] = new RooRealVar(strcore, strcore, initval, 0., 1e15);
    parameters.Lambda_z4qsq[v]->removeMax();

    initval=0;
    strcore = "cz_q";
    strcore.Append(Form("%i%s", (v!=3 ? v+1 : 12),"sq"));
    parameters.cLambda_qsq[v] = new RooRealVar(strcore, strcore, initval, -1., 1.);
  }

  if (parameterization==0){
    for (int v=0; v<8; v++){
      for (int im=0; im<2; im++){
        TString strcore;
        double initval = 0;
        TString strapp = "Val";
        if (im==1) strapp.Append("Im");
        if (v>1) strapp.Prepend(Form("%i", v));
        if (v>0) strapp.Prepend("_prime");

        strcore = "g1";
        strcore.Append(strapp);
        if (v==0 && im==0) initval = 1;
        RooRealVar* g1Val = new RooRealVar(strcore, strcore, initval, -1e15, 1e15);
        g1Val->removeMin();
        g1Val->removeMax();
        parameters.g1List[v][im] = (RooAbsReal*)g1Val;

        strcore = "g2";
        strcore.Append(strapp);
        RooRealVar* g2Val = new RooRealVar(strcore, strcore, 0, -1e15, 1e15);
        g2Val->removeMin();
        g2Val->removeMax();
        parameters.g2List[v][im] = (RooAbsReal*)g2Val;

        strcore = "g3";
        strcore.Append(strapp);
        RooRealVar* g3Val = new RooRealVar(strcore, strcore, 0, -1e15, 1e15);
        g3Val->removeMin();
        g3Val->removeMax();
        parameters.g3List[v][im] = (RooAbsReal*)g3Val;

        strcore = "g4";
        strcore.Append(strapp);
        RooRealVar* g4Val = new RooRealVar(strcore, strcore, 0, -1e15, 1e15);
        g4Val->removeMin();
        g4Val->removeMax();
        parameters.g4List[v][im] = (RooAbsReal*)g4Val;
      }
    }
  }
  else initFractionsPhases();
}

void ScalarPdfFactory::destroyMassPole(){
  delete parameters.mX;
  delete parameters.gamX;
}
void ScalarPdfFactory::destroyVdecayParams(){
  delete parameters.R1Val;
  delete parameters.R2Val;
  delete parameters.mV;
  delete parameters.gamV;
}
void ScalarPdfFactory::destroyFractionsPhases(){
  for (int v=0; v<8; v++){
    delete g1FracInterp[v];
    delete g2FracInterp[v];
    delete g3FracInterp[v];
    delete g4FracInterp[v];
  }
  delete gFracSum;
  for (int v=0; v<8; v++){
    delete g2Frac[v];
    delete g3Frac[v];
    delete g4Frac[v];
    delete g2Phase[v];
    delete g3Phase[v];
    delete g4Phase[v];
    if (v>0){
      delete g1Frac[v-1];
      delete g1Phase[v-1];
    }
  }
  for (int gg=0; gg<4; gg++){
    for (int v=0; v<8; v++) delete gRatioVal[gg][v];
  }
}
void ScalarPdfFactory::destroyGVals(){
  for (int v=0; v<8; v++){
    for (int im=0; im<2; im++){
      delete parameters.g1List[v][im];
      delete parameters.g2List[v][im];
      delete parameters.g3List[v][im];
      delete parameters.g4List[v][im];
    }
  }
  if (parameterization!=0) destroyFractionsPhases();
  for (int v=0; v<3; v++){
    delete parameters.Lambda_z1qsq[v];
    delete parameters.Lambda_z2qsq[v];
    delete parameters.Lambda_z3qsq[v];
    delete parameters.Lambda_z4qsq[v];
    delete parameters.cLambda_qsq[v];
  }
  delete parameters.Lambda;
  delete parameters.Lambda_z1;
  delete parameters.Lambda_z2;
  delete parameters.Lambda_z3;
  delete parameters.Lambda_z4;
  delete parameters.Lambda_Q;
}

void ScalarPdfFactory::addHypothesis(int ig, int ilam, double iphase, double altparam_fracval){
  if (ig>=4 || ig<0) cerr << "Invalid g" << ig << endl;
  if (ilam>=8 || ilam<0) cerr << "Out-of-range g" << ig << "_prime" << ilam << endl;
  if (parameterization==0){
    // Good guesses of c-constants
    Double_t initval=1.;
    if (ilam>0){
      if (ig==0){ // g1_dyn
        if (ilam>=1 && ilam<=3) initval = pow(parameters.Lambda_z1->getVal()/parameters.mV->getVal(), 2);
        if (ilam>=5 && ilam<=7) initval = pow(parameters.Lambda_z1->getVal()/parameters.mV->getVal(), 4);
        if (ilam==4) initval = pow(parameters.Lambda_Q->getVal()/parameters.mX->getVal(), 2);
      }
      else if (ig==1){ // g2_dyn
        if (ilam>=1 && ilam<=3) initval = pow(parameters.Lambda_z2->getVal()/parameters.mV->getVal(), 2);
        if (ilam>=5 && ilam<=7) initval = pow(parameters.Lambda_z2->getVal()/parameters.mV->getVal(), 4);
        if (ilam==4) initval = pow(parameters.Lambda_Q->getVal()/parameters.mX->getVal(), 2);
      }
      else if (ig==2){ // g3_dyn
        if (ilam>=1 && ilam<=3) initval = pow(parameters.Lambda_z3->getVal()/parameters.mV->getVal(), 2);
        if (ilam>=5 && ilam<=7) initval = pow(parameters.Lambda_z3->getVal()/parameters.mV->getVal(), 4);
        if (ilam==4) initval = pow(parameters.Lambda_Q->getVal()/parameters.mX->getVal(), 2);
      }
      else if (ig==3){ // g4_dyn
        if (ilam>=1 && ilam<=3) initval = pow(parameters.Lambda_z4->getVal()/parameters.mV->getVal(), 2);
        if (ilam>=5 && ilam<=7) initval = pow(parameters.Lambda_z4->getVal()/parameters.mV->getVal(), 4);
        if (ilam==4) initval = pow(parameters.Lambda_Q->getVal()/parameters.mX->getVal(), 2);
      }
    }
    if (ig==2) initval *= fabs(pow(parameters.Lambda->getVal(), 2)/(pow(parameters.mX->getVal(), 2) - pow(parameters.mV->getVal(), 2)));

    if (ig==0){
      ((RooRealVar*)parameters.g1List[ilam][0])->setVal(initval*cos(iphase));
      ((RooRealVar*)parameters.g1List[ilam][1])->setVal(initval*sin(iphase));
    }
    else if (ig==1){
      ((RooRealVar*)parameters.g2List[ilam][0])->setVal(initval*cos(iphase));
      ((RooRealVar*)parameters.g2List[ilam][1])->setVal(initval*sin(iphase));
    }
    else if (ig==2){
      ((RooRealVar*)parameters.g3List[ilam][0])->setVal(initval*cos(iphase));
      ((RooRealVar*)parameters.g3List[ilam][1])->setVal(initval*sin(iphase));
    }
    else if (ig==3){
      ((RooRealVar*)parameters.g4List[ilam][0])->setVal(initval*cos(iphase));
      ((RooRealVar*)parameters.g4List[ilam][1])->setVal(initval*sin(iphase));
    }
  }
  else{
    if (ig==0 && ilam==0) cerr << "Cannot set fa1! Try to set everything else." << endl;
    else{
      if (ig==0){
        g1Frac[ilam-1]->setVal(altparam_fracval);
        g1Phase[ilam-1]->setVal(iphase);
      }
      else if (ig==1){
        g2Frac[ilam]->setVal(altparam_fracval);
        g2Phase[ilam]->setVal(iphase);
      }
      else if (ig==2){
        g3Frac[ilam]->setVal(altparam_fracval);
        g3Phase[ilam]->setVal(iphase);
      }
      else if (ig==3){
        g4Frac[ilam]->setVal(altparam_fracval);
        g4Phase[ilam]->setVal(iphase);
      }
    }
  }
}
void ScalarPdfFactory::resetHypotheses(){
  for (int ilam=0; ilam<8; ilam++){
    if (parameterization==0){
      ((RooRealVar*)parameters.g1List[ilam][0])->setVal(0.); // This is different from setting the fractions!
      ((RooRealVar*)parameters.g1List[ilam][1])->setVal(0.);
      ((RooRealVar*)parameters.g2List[ilam][0])->setVal(0.);
      ((RooRealVar*)parameters.g2List[ilam][1])->setVal(0.);
      ((RooRealVar*)parameters.g3List[ilam][0])->setVal(0.);
      ((RooRealVar*)parameters.g3List[ilam][1])->setVal(0.);
      ((RooRealVar*)parameters.g4List[ilam][0])->setVal(0.);
      ((RooRealVar*)parameters.g4List[ilam][1])->setVal(0.);
    }
    else{
      if (ilam>0){
        g1Frac[ilam-1]->setVal(0.);
        g1Phase[ilam-1]->setVal(0.);
      }
      else{
        g2Frac[ilam]->setVal(0.);
        g2Phase[ilam]->setVal(0.);
        g3Frac[ilam]->setVal(0.);
        g3Phase[ilam]->setVal(0.);
        g4Frac[ilam]->setVal(0.);
        g4Phase[ilam]->setVal(0.);
      }
    }
  }
}


