#include <cmath>
#include "LHEParticleSmear.h"
#include "ParticleComparators.h"
#include "TRandom3.h"


namespace LHEParticleSmear{
  TRandom3 randomForSmearing;
}


MELAParticle* LHEParticleSmear::smearParticle(MELAParticle const* myParticle){
  if (!myParticle) return nullptr;

  randomForSmearing.SetSeed(std::abs(static_cast<int>(std::sin(myParticle->phi())*100000.)));
  TLorentzVector pOriginal(myParticle->p4), pSmeared(myParticle->p4);

  if (PDGHelpers::isALepton(myParticle->id)) pSmeared = LHEParticleSmear::smearLepton(pOriginal);
  else if (PDGHelpers::isAQuark(myParticle->id) || PDGHelpers::isAGluon(myParticle->id)) pSmeared = LHEParticleSmear::smearKnownJet(pOriginal);
  else if (PDGHelpers::isAnUnknownJet(myParticle->id)) pSmeared = LHEParticleSmear::smearUnknownJet(pOriginal);

  MELAParticle* smearedParticle = new MELAParticle(myParticle->id, pSmeared);
  return smearedParticle;
}

TLorentzVector LHEParticleSmear::smearLepton(TLorentzVector l_gen){
  float l_Perp, l_Theta, l_Phi;

  if (randomForSmearing.Uniform()<.9){
    l_Perp = l_gen.Pt()*pow(1.+0.012*1.15, randomForSmearing.Gaus(0., 1.));
    l_Theta = l_gen.Theta()+(randomForSmearing.Gaus(0, 0.001));
    l_Phi = l_gen.Phi()+(randomForSmearing.Gaus(0, 0.001));
  }
  else{
    l_Perp = l_gen.Pt()*(1.-0.04)*pow(1.+0.08, randomForSmearing.Gaus(0., 1.));
    l_Theta = l_gen.Theta()+(randomForSmearing.Gaus(0, 0.001));
    l_Phi = l_gen.Phi()+(randomForSmearing.Gaus(0, 0.001));
  }
  float l_Px = l_Perp*cos(l_Phi);
  float l_Py = l_Perp*sin(l_Phi);
  float l_Pz = l_Perp/tan(l_Theta);
  float l_E = sqrt(l_Px*l_Px+l_Py*l_Py+l_Pz*l_Pz+l_gen.M2());

  TLorentzVector final_l;
  final_l.SetXYZT(l_Px, l_Py, l_Pz, l_E);
  return final_l;
}
TLorentzVector LHEParticleSmear::smearKnownJet(TLorentzVector l_gen){
  double l_Pt, l_Eta, l_Phi, l_Mass;
  l_Pt = l_gen.Pt();
  l_Eta = l_gen.Eta();
  l_Phi = l_gen.Phi();
  l_Mass = l_gen.M();

  const double& dRjet = ParticleComparators::jetDeltaR;
  double etaSmear=0, phiSmear=0;
  double etaphiSmear = std::abs(randomForSmearing.Gaus(0, dRjet/2.));
  randomForSmearing.Circle(etaSmear, phiSmear, etaphiSmear);

  constexpr double massSmearConst = 6.; // GeV
  double massSmear = randomForSmearing.Exp(massSmearConst); // GeV

  constexpr double ptSmearConst = 0.05; // Ratio
  double ptSmear = randomForSmearing.Gaus(0, ptSmearConst);

  l_Pt *= std::max(0., 1.+ptSmear);
  l_Eta += etaSmear;
  l_Phi += phiSmear;
  l_Mass += massSmear;

  TLorentzVector final_l; final_l.SetPtEtaPhiM(l_Pt, l_Eta, l_Phi, l_Mass);
  final_l = smearUnknownJet(final_l);
  return final_l;
}
TLorentzVector LHEParticleSmear::smearUnknownJet(TLorentzVector l_gen){
  double l_Pt, l_Eta, l_Phi, l_Mass;
  l_Pt = l_gen.Pt();
  l_Eta = l_gen.Eta();
  l_Phi = l_gen.Phi();
  l_Mass = l_gen.M();

  /*
  const double& dRjet = ParticleComparators::jetDeltaR;
  double etaSmear=0, phiSmear=0;
  double etaphiSmear = std::abs(randomForSmearing.Gaus(0, dRjet));
  randomForSmearing.Circle(etaSmear, phiSmear, etaphiSmear);
  */

  constexpr double massSmearConst = 0.1; // Ratio
  double massSmear = randomForSmearing.Gaus(0, massSmearConst);

  constexpr double ptSmearConst = 0.1; // Ratio
  double ptSmear = randomForSmearing.Gaus(0, ptSmearConst);

  l_Pt *= std::max(0., 1.+ptSmear);
  //l_Eta += etaSmear;
  //l_Phi += phiSmear;
  l_Mass *= std::max(0., 1.+massSmear);

  TLorentzVector final_l; final_l.SetPtEtaPhiM(l_Pt, l_Eta, l_Phi, l_Mass);
  return final_l;
}


