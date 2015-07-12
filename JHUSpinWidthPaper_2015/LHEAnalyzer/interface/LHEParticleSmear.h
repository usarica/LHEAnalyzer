#ifndef LHEPARTICLESMEAR_H
#define LHEPARTICLESMEAR_H

#include "Particle.h"
#include "TRandom.h"

namespace LHEParticleSmear{
  TRandom randomForSmearing;


  Particle* smearParticle(Particle* myParticle);
  TLorentzVector smearLepton(TLorentzVector l_gen);
  TLorentzVector smearJet(TLorentzVector l_gen);

}

#endif
