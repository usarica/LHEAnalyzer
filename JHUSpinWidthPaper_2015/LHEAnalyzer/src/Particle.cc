#include "Particle.h"

Particle::Particle():
id(0),
passSelection(false)
{
  p4.SetXYZT(0, 0, 0, 0);
}

Particle::Particle(int id_, TLorentzVector p4_) :
id(id_),
passSelection(false)
{
  p4.SetXYZT(p4_.X(), p4_.Y(), p4_.Z(), p4_.T());
}
Particle::Particle(const Particle& particle_)
id(particle_.id),
passSelection(particle_.passSelection)
{
  p4.SetXYZT(particle_.p4.X(), particle_.p4.Y(), particle_.p4.Z(), particle_.p4.T());
  for (int index=0; index<particle_.getNMothers(); index++) addMother(particle_.getMother(index));
  for (int index=0; index<particle_.getNDaughters(); index++) addDaughter(particle_.getDaughter(index));
}
Particle::Particle& operator= (const Particle& particle_){
  id=particle_.id;
  passSelection=particle_.passSelection;
  for (int index=0; index<particle_.getNMothers(); index++) addMother(particle_.getMother(index));
  for (int index=0; index<particle_.getNDaughters(); index++) addDaughter(particle_.getDaughter(index));
}



Particle* Particle::getMother(int index){
  if ((int)mothers.size()>index) return mothers.at(index);
  else return 0;
}
Particle* Particle::getDaughter(int index){
  if ((int)daughters.size()>index) return daughters.at(index);
  else return 0;
}

double Particle::charge(){
  double cpos=0;
  if (abs(id)==24 || abs(id)==37 || abs(id)==2212 || abs(id)==211 || abs(id)==321 || abs(id)==411 || abs(id)==521) cpos = 1.;
  else if (abs(id)==11 || abs(id)==13 || abs(id)==15) cpos = -1.;
  else if (abs(id)==2 || abs(id)==4 || abs(id)==6) cpos = 2./3.;
  else if (abs(id)==1 || abs(id)==3 || abs(id)==5) cpos = -1./3.;
  if (id<0) cpos *= -1.;

  return cpos;
}