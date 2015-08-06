#ifndef PARTICLE_H
#define PARTICLE_H

#include <vector>
#include "TVector3.h"
#include "TLorentzVector.h"
#include "PDGHelpers.h"

class Particle{
public:

  // Constructors

  Particle();
  Particle(int id_, TLorentzVector p4_);
  Particle(const Particle& particle_);
  Particle& operator=(const Particle& particle_);
  virtual ~Particle(){};

  // Data

  int id;
  TLorentzVector p4;
  bool passSelection;
  int genStatus;
  double lifetime;


  // Member functions
  void setSelected(bool isSelected=true){ passSelection = isSelected; }
  void setGenStatus(int status_){ genStatus=status_; }
  void setLifetime(int life_){ lifetime=life_; }

  void addMother(Particle* myParticle);
  void addDaughter(Particle* myParticle);

  int getNMothers() const{ return mothers.size(); };
  int getNDaughters() const{ return daughters.size(); };

  Particle* getMother(int index) const;
  Particle* getDaughter(int index) const;

  double charge()const;
  double m()const{ return p4.M(); }
  double x()const{ return p4.X(); }
  double y()const{ return p4.Y(); }
  double z()const{ return p4.Z(); }
  double t()const{ return p4.T(); }
  double pt()const{ return p4.Pt(); }
  double eta()const{ return p4.Eta(); }
  double phi()const{ return p4.Phi(); }
  double rapidity()const{ return p4.Rapidity(); }
  double deltaR(const TLorentzVector& v)const{ return p4.DeltaR(v); }


protected:
  std::vector<Particle*> mothers;
  std::vector<Particle*> daughters;

  bool checkParticleExists(Particle* myParticle, std::vector<Particle*>& particleArray);
};


#endif
