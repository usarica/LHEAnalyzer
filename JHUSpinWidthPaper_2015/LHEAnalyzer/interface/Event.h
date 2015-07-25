#ifndef EVENTBASE_H
#define EVENTBASE_H

#include <vector>
#include "TLorentzVector.h"
#include "ZZCandidate.h"
#include "ParticleComparators.h"

class Event{
public:

  // Constructors

  Event(){ weight=0; xsec=0; };
  ~Event(){ wipeAll(); };

  // Data

  // Member functions
  double getXSec() const{ return xsec; }
  double getWeight() const{ return weight; }
  std::vector<double>& getExtraWeights(){ return extraWeight; }
  void setXSec(double xsec_){ xsec=xsec_; }
  void setWeight(double weight_){ weight=weight_; }

  void constructVVCandidates(bool isZZ=true, int fstype=0);
  void applyParticleSelection();
  void addVVCandidateAppendages();


  int getNZZCandidates() const{ return ZZcandidates.size(); };
  int getNLeptons() const{ return leptons.size(); };
  int getNNeutrinos() const{ return neutrinos.size(); };
  int getNJets() const{ return jets.size(); };
  int getNParticles() const{ return particles.size(); };

  ZZCandidate* getZZCandidate(int index) const;
  Particle* getLepton(int index) const;
  Particle* getNeutrino(int index) const;
  Particle* getJet(int index) const;
  Particle* getParticle(int index) const;

  void addParticle(Particle* myParticle){ particles.push_back(myParticle); };
  void addLepton(Particle* myParticle, bool genuineParticle=true){ leptons.push_back(myParticle); if (genuineParticle) addParticle(myParticle); };
  void addNeutrino(Particle* myParticle, bool genuineParticle=true){ neutrinos.push_back(myParticle); if (genuineParticle) addParticle(myParticle); };
  void addJet(Particle* myParticle, bool genuineParticle=true){ jets.push_back(myParticle); if (genuineParticle) addParticle(myParticle); };
  void addZZCandidate(ZZCandidate* myParticle){ ZZcandidates.push_back(myParticle); };

  TLorentzVector missingP() const;

protected:
  std::vector<Particle*> particles;
  std::vector<Particle*> leptons;
  std::vector<Particle*> neutrinos;
  std::vector<Particle*> jets;
  std::vector<ZZCandidate*> ZZcandidates;

  template<typename ParticleType> void wipeArray(std::vector<ParticleType*>& particleArray, bool doDelete=true){ if (doDelete){ for (int i=0; i<particleArray.size(); i++){ ParticleType* delpar = particleArray.at(i); delete delpar; } } particleArray.clear(); };
  void wipeAll(){ leptons.clear(); neutrinos.clear(); jets.clear(); wipeArray(ZZcandidates, true); wipeArray(particles, false); };

  void applyLeptonSelection();
  void applyNeutrinoSelection();
  void applyJetSelection();
  void applyZZSelection();

  double xsec;
  double weight;
  std::vector<double> extraWeight;
};


#endif
