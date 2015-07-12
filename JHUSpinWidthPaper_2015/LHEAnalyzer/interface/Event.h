#ifndef EVENTBASE_H
#define EVENTBASE_H

#include <vector>
#include "TLorentzVector.h"
#include "ZZCandidate.h"


class Event{
public:

  // Constructors

  Event(){};
  ~Event(){ wipeAll(); };

  // Data

  // Member functions
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

  template<typename ParticleType> void wipeArray(std::vector<ParticleType*>& particleArray){ for (int i=0; i<particleArray.size(); i++) delete particleArray.at(i); particleArray.clear(); };
  void wipeAll(){ leptons.clear(); jets.clear(); wipeArray(ZZcandidates); wipeArray(particles); };

  void applyLeptonSelection();
  void applyNeutrinoSelection();
  void applyJetSelection();
  void applyZZSelection();
};


#endif
