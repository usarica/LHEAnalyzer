#ifdef _def_melatools_
#include <ZZMatrixElement/MELA/interface/SpinPdfFactory.h>
#else
#include "../include/SpinPdfFactory.h"
#endif


SpinPdfFactory::SpinPdfFactory(RooSpin::modelMeasurables measurables_, RooSpin::VdecayType V1decay_, RooSpin::VdecayType V2decay_) :
V1decay(V1decay_),
V2decay(V2decay_)
{
  initMeasurables(measurables_);
  initMassPole();
  initVdecayParams();
}

SpinPdfFactory::~SpinPdfFactory(){
  destroyVdecayParams();
  destroyMassPole();
}
void SpinPdfFactory::initMeasurables(RooSpin::modelMeasurables measurables_){
  measurables.h1 = (RooAbsReal*)measurables_.h1;
  measurables.h2 = (RooAbsReal*)measurables_.h2;
  measurables.Phi = (RooAbsReal*)measurables_.Phi;
  measurables.m1 = (RooAbsReal*)measurables_.m1;
  measurables.m2 = (RooAbsReal*)measurables_.m2;
  measurables.m12 = (RooAbsReal*)measurables_.m12;
  measurables.hs = (RooAbsReal*)measurables_.hs;
  measurables.Phi1 = (RooAbsReal*)measurables_.Phi1;
  measurables.Y = (RooAbsReal*)measurables_.Y;
}
void SpinPdfFactory::initMassPole(){
  parameters.mX = new RooRealVar("mX", "mX", (measurables.m12)->getVal());
  parameters.gamX = new RooRealVar("gamX", "gamX", 0);
}
void SpinPdfFactory::initVdecayParams(){
  if ((((int)V1decay)>0 && ((int)V2decay)<0) || (((int)V1decay)<0 && ((int)V2decay)>0)) cerr << "SpinPdfFactory::initVdecayParams: V1 and V2 decays are inconsistent!" << endl;

  const Double_t GfVal = 1.16639e-5;
  const Double_t vevVal = 1./sqrt(GfVal*sqrt(2.));

  parameters.mW = new RooRealVar("mW", "mW", 80.399);
  parameters.gamW = new RooRealVar("gamW", "gamW", 2.085);
  parameters.mZ = new RooRealVar("mZ", "mZ", 91.1876);
  parameters.gamZ = new RooRealVar("gamZ", "gamZ", 2.4952);
  parameters.Sin2ThetaW = new RooRealVar("Sin2ThetaW", "Sin2ThetaW", 0.23119);
  parameters.vev = new RooRealVar("vev", "vev", vevVal);
}
void SpinPdfFactory::getMVGamV(Double_t* mV, Double_t* gamV)const{
  if (V1decay==RooSpin::kVdecayType_Wany){
    if (mV!=0) (*mV)=(parameters.mW)->getVal();
    if (gamV!=0) (*gamV)=(parameters.gamW)->getVal();
  }
  else if (!(V1decay==RooSpin::kVdecayType_GammaOnshell && V2decay==RooSpin::kVdecayType_GammaOnshell)){
    if (mV!=0) (*mV)=(parameters.mZ)->getVal();
    if (gamV!=0) (*gamV)=(parameters.gamZ)->getVal();
  }
  else{
    if (mV!=0) (*mV)=0;
    if (gamV!=0) (*gamV)=0;
  }
}
void SpinPdfFactory::resetVdecay(RooSpin::VdecayType V1decay_, RooSpin::VdecayType V2decay_){
  if ((((int)V1decay)>0 && ((int)V2decay)<0) || (((int)V1decay)<0 && ((int)V2decay)>0)) cerr << "SpinPdfFactory::resetVdecay: V1 and V2 decays are inconsistent!" << endl;

  V1decay=V1decay_;
  V2decay=V2decay_;
  PDF_base->setDecayModes(V1decay, V2decay);
}

void SpinPdfFactory::destroyMassPole(){
  delete parameters.gamX;
  delete parameters.mX;
}
void SpinPdfFactory::destroyVdecayParams(){
  delete parameters.vev;
  delete parameters.Sin2ThetaW;
  delete parameters.gamZ;
  delete parameters.mZ;
  delete parameters.gamW;
  delete parameters.mW;
}


