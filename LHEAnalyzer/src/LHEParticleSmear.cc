#include "LHEParticleSmear.h"

namespace LHEParticleSmear{
  TRandom randomForSmearing;
}

MELAParticle* LHEParticleSmear::smearParticle(MELAParticle* myParticle){
  TLorentzVector pOriginal, pSmeared;
  pOriginal.SetXYZT(myParticle->x(), myParticle->y(), myParticle->z(), myParticle->t());

  if (PDGHelpers::isALepton(myParticle->id)){
    pSmeared = LHEParticleSmear::smearLepton(pOriginal);
  }
  else if (PDGHelpers::isAQuark(myParticle->id) || PDGHelpers::isAGluon(myParticle->id)){
    pSmeared = LHEParticleSmear::smearJet(pOriginal);
  }
  else pSmeared = pOriginal;

  MELAParticle* smearedParticle = new MELAParticle(myParticle->id, pSmeared);
  return smearedParticle;
}

TLorentzVector LHEParticleSmear::smearLepton(TLorentzVector l_gen){
  float l_Perp, l_Theta, l_Phi;

  if (LHEParticleSmear::randomForSmearing.Uniform()<.9){
    l_Perp = l_gen.Pt()*pow(1.+0.012*1.15, LHEParticleSmear::randomForSmearing.Gaus(0., 1.));
    l_Theta = l_gen.Theta()+(LHEParticleSmear::randomForSmearing.Gaus(0, 0.001));
    l_Phi = l_gen.Phi()+(LHEParticleSmear::randomForSmearing.Gaus(0, 0.001));
  }
  else{
    l_Perp = l_gen.Pt()*(1.-0.04)*pow(1.+0.08, LHEParticleSmear::randomForSmearing.Gaus(0., 1.));
    l_Theta = l_gen.Theta()+(LHEParticleSmear::randomForSmearing.Gaus(0, 0.001));
    l_Phi = l_gen.Phi()+(LHEParticleSmear::randomForSmearing.Gaus(0, 0.001));
  }
  float l_Px = l_Perp*cos(l_Phi);
  float l_Py = l_Perp*sin(l_Phi);
  float l_Pz = l_Perp/tan(l_Theta);
  float l_E = sqrt(l_Px*l_Px+l_Py*l_Py+l_Pz*l_Pz+l_gen.M2());

  TLorentzVector final_l;
  final_l.SetXYZT(l_Px, l_Py, l_Pz, l_E);
  return final_l;
}
TLorentzVector LHEParticleSmear::smearJet(TLorentzVector l_gen){
  float l_Perp, l_Theta, l_Phi;

  // Assume the same smearing for pT, theta and phi
  if (LHEParticleSmear::randomForSmearing.Uniform()<.9){
    l_Perp = l_gen.Pt()*pow(1.+0.012*1.15, LHEParticleSmear::randomForSmearing.Gaus(0., 1.));
    l_Theta = l_gen.Theta()+(LHEParticleSmear::randomForSmearing.Gaus(0, 0.001));
    l_Phi = l_gen.Phi()+(LHEParticleSmear::randomForSmearing.Gaus(0, 0.001));
  }
  else{
    l_Perp = l_gen.Pt()*(1.-0.04)*pow(1.+0.08, LHEParticleSmear::randomForSmearing.Gaus(0., 1.));
    l_Theta = l_gen.Theta()+(LHEParticleSmear::randomForSmearing.Gaus(0, 0.001));
    l_Phi = l_gen.Phi()+(LHEParticleSmear::randomForSmearing.Gaus(0, 0.001));
  }

  float l_E = l_gen.T()*pow(1.15, LHEParticleSmear::randomForSmearing.Gaus(0., 1.));
  float l_Px = l_Perp*cos(l_Phi);
  float l_Py = l_Perp*sin(l_Phi);
  float l_Pz = l_Perp/tan(l_Theta);

  TLorentzVector final_l;
  final_l.SetXYZT(l_Px, l_Py, l_Pz, l_E);
  return final_l;
}


