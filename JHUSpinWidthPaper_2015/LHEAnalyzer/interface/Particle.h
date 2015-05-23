#include <vector>
#include "TLorentzVector.h"
//#include <string>
//#include "TString.h"

class Particle{
public:

  // Constructors

  Particle();
  Particle(int id_, TLorentzVector p4_);
  Particle(const Particle& particle_);
  Particle& operator= (const Particle& particle_);
  ~Particle(){ daughters.clear(); mothers.clear(); };

  // Data

  int id;
  TLorentzVector p4;
  bool passSelection;


  // Member functions
  void setSelected(bool isSelected=true){ passSelection = isSelected; }

  void addMother(Particle* myParticle){ mothers.push_back(myParticle); };
  void addDaughter(Particle* myParticle){ daughters.push_back(myParticle); };

  int getNMothers(){ return mothers.size() };
  int getNDaughters(){ return daughters.size() };

  Particle* getMother(int index);
  Particle* getDaughter(int index);

  double charge();
  double m(){ return p4.M(); }
  double x(){ return p4.X(); }
  double y(){ return p4.Y(); }
  double z(){ return p4.Z(); }
  double t(){ return p4.T(); }
  double pt(){ return p4.Pt(); }
  double eta(){ return p4.Eta(); }
  double phi(){ return p4.Phi(); }
  double rapidity(){ return p4.Rapidity(); }
  double deltaR(const TLorentzVector& v){ return p4.DeltaR(v); }


protected:
  std::vector<Particle*> mothers;
  std::vector<Particle*> daughters;
}

