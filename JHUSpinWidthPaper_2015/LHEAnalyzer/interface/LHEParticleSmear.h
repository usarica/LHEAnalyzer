#ifndef LHEPARTICLESMEAR_H
#define LHEPARTICLESMEAR_H

#include "MELAParticle.h"
#include "TRandom.h"

namespace LHEParticleSmear{
  extern TRandom randomForSmearing;


  MELAParticle* smearParticle(MELAParticle* myParticle);
  TLorentzVector smearLepton(TLorentzVector l_gen);
  TLorentzVector smearJet(TLorentzVector l_gen);

}

#endif
