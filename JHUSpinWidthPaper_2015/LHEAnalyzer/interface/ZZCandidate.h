#include "Particle.h"

class ZZCandidate : public Particle{
public:
  ~ZZCandidate(){ for (int i=0; i<sortedVs.size(); i++) delete sortedVs.at(i); sortedVs.clear(); sortedDaughters.clear(); associatedJets.clear(); associatedLeptons.clear(); };


  std::vector<Particle*> associatedLeptons;
  std::vector<Particle*> associatedJets;

  std::vector<Particle*> sortedDaughters;
  std::vector<Particle*> sortedVs;


  // Member functions

  Particle* getSortedDaughter(int index);
  Particle* getSortedV(int index);
  void addSortedV(Particle* myParticle){ sortedVs.push_back(myParticle); };


private:
  void sortDaughters();
  void sortDaughtersInitial();
  void sortDaughtersByBestZ1();
  void createSortedVs();
}
