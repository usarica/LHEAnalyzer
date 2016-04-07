#ifdef _def_melatools_
#include <ZZMatrixElement/MELA/interface/TensorPdfFactory.h>
#else
#include "../include/TensorPdfFactory.h"
#endif


TensorPdfFactory::TensorPdfFactory(RooSpinTwo::modelMeasurables measurables_, RooSpin::VdecayType V1decay_, RooSpin::VdecayType V2decay_) :
V1decay(V1decay_),
V2decay(V2decay_)
{
  initMeasurables(measurables_);
  initMassPole();
  initVdecayParams();
  initGVals();
}

TensorPdfFactory::~TensorPdfFactory(){
  destroyGVals();
  destroyVdecayParams();
  destroyMassPole();
}
void TensorPdfFactory::initMeasurables(RooSpinTwo::modelMeasurables measurables_){
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
void TensorPdfFactory::initMassPole(){
  parameters.mX = new RooRealVar("mX", "mX", (measurables.m12)->getVal());
  parameters.gamX = new RooRealVar("gamX", "gamX", 4.07e-3);
}
void TensorPdfFactory::initVdecayParams(){
  if ((((int)V1decay)>0 && ((int)V2decay)<0) || (((int)V1decay)<0 && ((int)V2decay)>0)) cerr << "TensorPdfFactory::initVdecayParams: V1 and V2 decays are inconsistent!" << endl;
  if (V1decay==RooSpin::kVdecayType_GammaOnshell){
    parameters.mV = new RooRealVar("mV", "mV", 0);
    parameters.gamV = new RooRealVar("gamV", "gamV", 0);
  }
  else if (V1decay==RooSpin::kVdecayType_Wany){
    parameters.mV = new RooRealVar("mV", "mV", 80.399);
    parameters.gamV = new RooRealVar("gamV", "gamV", 2.085);
  }
  else{
    parameters.mV = new RooRealVar("mV", "mV", 91.1876);
    parameters.gamV = new RooRealVar("gamV", "gamV", 2.4952);
  }
  parameters.R1Val = new RooRealVar("R1Val", "R1Val", getRValue(V1decay));
  parameters.R2Val = new RooRealVar("R2Val", "R2Val", getRValue(V2decay));

  parameters.mV->removeRange();
  parameters.gamV->removeRange();
  parameters.R1Val->removeRange();
  parameters.R2Val->removeRange();
}
void TensorPdfFactory::resetVdecay(RooSpin::VdecayType V1decay_, RooSpin::VdecayType V2decay_){
  if ((((int)V1decay)>0 && ((int)V2decay)<0) || (((int)V1decay)<0 && ((int)V2decay)>0)) cerr << "TensorPdfFactory::resetVdecay: V1 and V2 decays are inconsistent!" << endl;

  V1decay=V1decay_;
  V2decay=V2decay_;
  PDF_base->setDecayModes(V1decay, V2decay);

  bool is_mvconst = parameters.mV->isConstant();
  bool is_gamvconst = parameters.gamV->isConstant();
  bool is_r1const = parameters.R1Val->isConstant();
  bool is_r2const = parameters.R2Val->isConstant();
  parameters.mV->setConstant(false);
  parameters.gamV->setConstant(false);
  parameters.R1Val->setConstant(false);
  parameters.R2Val->setConstant(false);

  if (V1decay_==RooSpin::kVdecayType_GammaOnshell){
    parameters.mV = new RooRealVar("mV", "mV", 0);
    parameters.gamV = new RooRealVar("gamV", "gamV", 0);
  }
  else if (V1decay_==RooSpin::kVdecayType_Wany){
    parameters.mV = new RooRealVar("mV", "mV", 80.399);
    parameters.gamV = new RooRealVar("gamV", "gamV", 2.085);
  }
  else{
    parameters.mV = new RooRealVar("mV", "mV", 91.1876);
    parameters.gamV = new RooRealVar("gamV", "gamV", 2.4952);
  }

  parameters.R1Val->setVal(getRValue(V1decay));
  parameters.R2Val->setVal(getRValue(V2decay));

  parameters.mV->setConstant(is_mvconst);
  parameters.gamV->setConstant(is_gamvconst);
  parameters.R1Val->setConstant(is_r1const);
  parameters.R2Val->setConstant(is_r2const);
}
double TensorPdfFactory::getRValue(RooSpin::VdecayType Vdecay){
  double atomicT3 = 0.5;
  double atomicCharge = 1.;
  double sin2t = 0.23119;

  // gV = T3 - 2*Qf*sintW**2
  double gV_up = atomicT3 - 2.*(2.*atomicCharge/3.)*sin2t;
  double gV_dn = -atomicT3 - 2.*(-atomicCharge/3.)*sin2t;
  double gV_l = -atomicT3 - 2.*(-atomicCharge)*sin2t;
  double gV_nu = atomicT3;

  // gA = T3
  double gA_up = atomicT3;
  double gA_dn = -atomicT3;
  double gA_l = -atomicT3;
  double gA_nu = atomicT3;

  double Rval;
  switch (Vdecay){
  case RooSpin::kVdecayType_Zud:
    Rval = (2.*2.*gV_up*gA_up+3.*2.*gV_dn*gA_dn)/(2.*(pow(gV_up, 2)+pow(gA_up, 2))+3.*(pow(gV_dn, 2)+pow(gA_dn, 2))); // Z->uu+dd avg.
    break;
  case RooSpin::kVdecayType_Zdd:
    Rval = (2.*gV_dn*gA_dn)/(pow(gV_dn, 2)+pow(gA_dn, 2)); // Z->dd
    break;
  case RooSpin::kVdecayType_Zuu:
    Rval = (2.*gV_up*gA_up)/(pow(gV_up, 2)+pow(gA_up, 2)); // Z->uu
    break;
  case RooSpin::kVdecayType_Znn:
    Rval = (2.*gV_nu*gA_nu)/(pow(gV_nu, 2)+pow(gA_nu, 2)); // Z->nunu
    break;
  case RooSpin::kVdecayType_Zll:
    Rval = (2.*gV_l*gA_l)/(pow(gV_l, 2)+pow(gA_l, 2)); // Z->ll
    break;
  case RooSpin::kVdecayType_Wany:
    Rval = 1; // W
    break;
  default:
    Rval = 0; // gamma
    break;
  }
  return Rval;
}
void TensorPdfFactory::initGVals(){
  couplings.Lambda = new RooRealVar("Lambda", "Lambda", 1000.);

  for (int v=0; v<10; v++){
    for (int im=0; im<2; im++){
      TString strcore;
      double initval = 0;
      TString strapp = "Val";
      if (im==1) strapp.Append("Im");
      strapp.Prepend(Form("%i", v+1));

      strcore = "g";
      strcore.Append(strapp);
      RooRealVar* gVal = new RooRealVar(strcore, strcore, initval, -1e15, 1e15);
      gVal->removeMin();
      gVal->removeMax();
      couplings.bList[v][im] = (RooAbsReal*)gVal;
    }
  }

  for (int f=0; f<2; f++){
    TString strcore = Form("f_spinz%i", f+1);
    RooRealVar* fVal = new RooRealVar(strcore, strcore, 0., 0., 1.);
    if (f==0) couplings.f_spinz1 = (RooRealVar*)fVal;
    else couplings.f_spinz2 = (RooRealVar*)fVal;
  }
}

void TensorPdfFactory::destroyMassPole(){
  delete parameters.mX;
  delete parameters.gamX;
}
void TensorPdfFactory::destroyVdecayParams(){
  delete parameters.R1Val;
  delete parameters.R2Val;
  delete parameters.mV;
  delete parameters.gamV;
}
void TensorPdfFactory::destroyGVals(){
  for (int v=0; v<10; v++){
    for (int im=0; im<2; im++){
      delete couplings.bList[v][im];
    }
  }
  delete couplings.Lambda;

  delete couplings.f_spinz1;
  delete couplings.f_spinz2;
}

void TensorPdfFactory::addHypothesis(int ig, double initval, double iphase){
  if (ig>=10 || ig<0) cerr << "Invalid g" << ig << endl;
  else{
    ((RooRealVar*)couplings.bList[ig][0])->setVal(initval*cos(iphase));
    ((RooRealVar*)couplings.bList[ig][1])->setVal(initval*sin(iphase));
  }
}
void TensorPdfFactory::setTensorPolarization(int ig, double initval){
  if (ig>2 || ig<=0) cerr << "Cannot set f_spinz" << ig << ". Please st f_spinz1 or f_spinz2 only." << endl;
  else{
    if (ig==1) ((RooRealVar*)couplings.f_spinz1)->setVal(initval);
    else ((RooRealVar*)couplings.f_spinz2)->setVal(initval);
  }
}
void TensorPdfFactory::resetHypotheses(){
  for (int ig=0; ig<10; ig++){
    for (int im=0; im<2; im++) ((RooRealVar*)couplings.bList[ig][im])->setVal(0.);
  }
  ((RooRealVar*)couplings.f_spinz1)->setVal(0.);
  ((RooRealVar*)couplings.f_spinz2)->setVal(0.);
}


