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
  switch (V1decay){
  case 2:
    parameters.R1Val = new RooRealVar("R1Val", "R1Val", 1);
    break;
  case 1:
    parameters.R1Val = new RooRealVar("R1Val", "R1Val", 0.15);
    break;
  case -1:
    parameters.R1Val = new RooRealVar("R1Val", "R1Val", 1);
    break;
  default:
    parameters.R1Val = new RooRealVar("R1Val", "R1Val", 0);
    break;
  }
  switch (V2decay){
  case 2:
    parameters.R2Val = new RooRealVar("R2Val", "R2Val", 1);
    break;
  case 1:
    parameters.R2Val = new RooRealVar("R2Val", "R2Val", 0.15);
    break;
  case -1:
    parameters.R2Val = new RooRealVar("R2Val", "R2Val", 1);
    break;
  default:
    parameters.R2Val = new RooRealVar("R2Val", "R2Val", 0);
    break;
  }
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
  delete parameters.Lambda;
  delete parameters.Lambda_z1;
  delete parameters.Lambda_z2;
  delete parameters.Lambda_z3;
  delete parameters.Lambda_z4;
  delete parameters.Lambda_Q;
}


