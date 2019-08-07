#ifndef LHEPARTICLESMEAR_H
#define LHEPARTICLESMEAR_H

#include "MELAParticle.h"


namespace LHEParticleSmear{

  MELAParticle* smearParticle(MELAParticle const* myParticle);
  TLorentzVector smearLepton(TLorentzVector l_gen);
  TLorentzVector smearKnownJet(TLorentzVector l_gen);
  TLorentzVector smearUnknownJet(TLorentzVector l_gen);

}

#endif
