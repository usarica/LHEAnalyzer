#ifndef ZZCANDIDATE_H
#define ZZCANDIDATE_H

#include "Particle.h"
#include "TopCandidate.h"

class ZZCandidate : public Particle{
public:
  ZZCandidate(int id_, TLorentzVector p4_, bool associatedByHighestPt_=true);
  ~ZZCandidate();
  ZZCandidate* shallowCopy();

  // Member functions

  Particle* getSortedDaughter(int index)const;
  Particle* getSortedV(int index)const;
  Particle* getAssociatedLepton(int index)const;
  Particle* getAssociatedNeutrino(int index)const;
  Particle* getAssociatedPhoton(int index)const;
  Particle* getAssociatedJet(int index)const;
  TopCandidate* getAssociatedTop(int index)const;
  TLorentzVector getAlternativeVMomentum(int index)const;

  virtual std::vector<int> getDaughterIds()const;
  std::vector<int> getAssociatedParticleIds()const;

  int getNAssociatedLeptons()const{ return associatedLeptons.size(); }
  int getNAssociatedNeutrinos()const{ return associatedNeutrinos.size(); }
  int getNAssociatedPhotons()const{ return associatedPhotons.size(); }
  int getNAssociatedJets()const{ return associatedJets.size(); }
  int getNAssociatedTops()const{ return associatedTops.size(); }
  int getNSortedVs()const{ return sortedVs.size(); }

  void addAssociatedLeptons(Particle* myParticle);
  void addAssociatedNeutrinos(Particle* myParticle);
  void addAssociatedPhotons(Particle* myParticle);
  void addAssociatedJets(Particle* myParticle);
  void addAssociatedTops(TopCandidate* myParticle);

  void addSortedV(Particle* myParticle){ sortedVs.push_back(myParticle); }
  void addAssociatedVs();

  void sortDaughters();
  void testPreSelectedDaughters();
  bool testShallowCopy();

  bool daughtersInterfere()const;
  void setAddAssociatedByHighestPt(bool associatedByHighestPt_);
  void setShallowCopy(bool flag);

protected:
  bool associatedByHighestPt;
  bool isShallowCopy;

  std::vector<Particle*> associatedLeptons;
  std::vector<Particle*> associatedNeutrinos;
  std::vector<Particle*> associatedPhotons;
  std::vector<Particle*> associatedJets;
  std::vector<TopCandidate*> associatedTops;

  std::vector<Particle*> sortedDaughters;
  std::vector<Particle*> sortedVs;

  void sortDaughtersInitial();
  void sortDaughtersByBestZ1();
  void createSortedVs();
  bool checkDaughtership(Particle* myParticle)const;
  void createAssociatedVs(std::vector<Particle*>& particleArray);
  void addByHighestPt(Particle* myParticle, std::vector<Particle*>& particleArray);
  void addByHighestPt(TopCandidate* myParticle, std::vector<TopCandidate*>& particleArray);

  bool checkTopCandidateExists(TopCandidate* myParticle, std::vector<TopCandidate*>& particleArray);
};



#endif
