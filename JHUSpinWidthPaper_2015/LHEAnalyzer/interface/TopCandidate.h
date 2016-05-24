#ifndef TOPCANDIDATE_H
#define TOPCANDIDATE_H

#include "Particle.h"

class TopCandidate : public Particle{
public:
  TopCandidate(int id_, TLorentzVector p4_) : Particle(id_, p4_), lightQuark(0), Wferm(0), Wfermbar(0) {}
  TopCandidate(Particle* lightQuark_, Particle* Wferm_, Particle* Wfermbar_);
  ~TopCandidate(){}

  Particle* setLightQuark(Particle* myParticle){ lightQuark=myParticle; }
  Particle* setWFermion(Particle* myParticle){ Wferm=myParticle; }
  Particle* setWAntifermion(Particle* myParticle){ Wfermbar=myParticle; }

  Particle* getLightQuark(){ return lightQuark; }
  Particle* getWFermion(){ return Wferm; }
  Particle* getWAntifermion(){ return Wfermbar; }

protected:
  Particle* lightQuark;
  Particle* Wferm;
  Particle* Wfermbar;

};


#endif
