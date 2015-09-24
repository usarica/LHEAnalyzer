#include "../interface/ZZCandidate.h"

void ZZCandidate::sortDaughters(){
  sortDaughtersInitial();
  sortDaughtersByBestZ1();
  createSortedVs();
}

Particle* ZZCandidate::getSortedDaughter(int index) const{
  if ((int)sortedDaughters.size()>index) return sortedDaughters.at(index);
  else return 0;
}
Particle* ZZCandidate::getSortedV(int index) const{
  if ((int)sortedVs.size()>index) return sortedVs.at(index);
  else return 0;
}
Particle* ZZCandidate::getAssociatedLepton(int index)const{
  if ((int)associatedLeptons.size()>index) return associatedLeptons.at(index);
  else return 0;
}
Particle* ZZCandidate::getAssociatedNeutrino(int index)const{
  if ((int)associatedNeutrinos.size()>index) return associatedNeutrinos.at(index);
  else return 0;
}
Particle* ZZCandidate::getAssociatedJet(int index)const{
  if ((int)associatedJets.size()>index) return associatedJets.at(index);
  else return 0;
}
void ZZCandidate::sortDaughtersInitial(){
  int tmpDindex[2]={ 0 };
  Particle* df[2] = { getDaughter(0), 0 };
  for (int j=1; j<getNDaughters(); j++){
    Particle* dtmp = getDaughter(j);
    if ((dtmp->charge()+df[0]->charge()==0 && PDGHelpers::HVVmass==PDGHelpers::Zmass) || (std::abs(dtmp->charge()+df[0]->charge())==1 && PDGHelpers::HVVmass==PDGHelpers::Wmass)){
      df[1] = dtmp;
      tmpDindex[1] = j;
      break;
    }
  }
  Particle* ds[2] ={ 0 };
  int sindex=0;
  for (int j=1; j<getNDaughters(); j++){
    if (j==tmpDindex[1]) continue;
    Particle* dtmp = getDaughter(j);
    ds[sindex] = dtmp;
    sindex++;
  }
  if (
    (df[0]!=0 && df[1]!=0)
    &&
    (
    (df[0]->id<df[1]->id && PDGHelpers::HVVmass==PDGHelpers::Zmass) || (std::abs(df[0]->id)<std::abs(df[1]->id) && PDGHelpers::HVVmass==PDGHelpers::Wmass)
    )
    ){
    Particle* dtmp = df[0];
    df[0] = df[1];
    df[1] = dtmp;
  }
  if (
    (ds[0]!=0 && ds[1]!=0)
    &&
    (
    (ds[0]->id<ds[1]->id && PDGHelpers::HVVmass==PDGHelpers::Zmass) || (std::abs(ds[0]->id)<std::abs(ds[1]->id) && PDGHelpers::HVVmass==PDGHelpers::Wmass)
    )
    ){
    Particle* dtmp = ds[0];
    ds[0] = ds[1];
    ds[1] = dtmp;
  }
  for (int i=0; i<2; i++){
    if (df[i]!=0) sortedDaughters.push_back(df[i]);
  }
  for (int i=0; i<2; i++){
    if (ds[i]!=0) sortedDaughters.push_back(ds[i]);
  }
}
void ZZCandidate::sortDaughtersByBestZ1(){
  Particle* orderedDs[2][2]={ { 0 } };

  TLorentzVector pZ1(0, 0, 0, 0);
  TLorentzVector pZ2(0, 0, 0, 0);
  for (int d=0; d<2; d++){
    if (sortedDaughters.at(d)!=0) pZ1 = pZ1 + sortedDaughters.at(d)->p4;
  }
  for (int d=2; d<4; d++){
    if (sortedDaughters.at(d)!=0) pZ2 = pZ2 + sortedDaughters.at(d)->p4;
  }
  if (std::abs(pZ1.M() - PDGHelpers::HVVmass)<std::abs(pZ2.M() - PDGHelpers::HVVmass)){
    orderedDs[0][0]=sortedDaughters.at(0);
    orderedDs[0][1]=sortedDaughters.at(1);
    orderedDs[1][0]=sortedDaughters.at(2);
    orderedDs[1][1]=sortedDaughters.at(3);
  }
  else{
    orderedDs[0][0]=sortedDaughters.at(2);
    orderedDs[0][1]=sortedDaughters.at(3);
    orderedDs[1][0]=sortedDaughters.at(0);
    orderedDs[1][1]=sortedDaughters.at(1);
    TLorentzVector ptmp = pZ1;
    pZ1 = pZ2;
    pZ2 = ptmp;
  }
  if (
    (
    (orderedDs[0][1]!=0 && orderedDs[1][1]!=0)
    &&
    (orderedDs[0][1]->id == orderedDs[1][1]->id)
    )
    ||
    (
    (orderedDs[1][0]!=0 && orderedDs[0][0]!=0)
    &&
    (orderedDs[1][0]->id == orderedDs[0][0]->id)
    )
    ){
    Particle* orderedDps[2][2]={ { 0 } };

    TLorentzVector pZ1p(0, 0, 0, 0);
    TLorentzVector pZ2p(0, 0, 0, 0);
    for (int d=0; d<2; d++){
      if (orderedDs[d][d]!=0) pZ1p = pZ1p + orderedDs[d][d]->p4;
    }
    for (int d=0; d<2; d++){
      if (orderedDs[1-d][d]!=0) pZ2p = pZ2p + orderedDs[1-d][d]->p4;
    }

    if (std::abs(pZ1p.M() - PDGHelpers::HVVmass)<std::abs(pZ2p.M() - PDGHelpers::HVVmass)){
      orderedDps[0][0]=orderedDs[0][0];
      orderedDps[0][1]=orderedDs[1][1];
      orderedDps[1][0]=orderedDs[1][0];
      orderedDps[1][1]=orderedDs[0][1];
    }
    else{
      orderedDps[0][0]=orderedDs[1][0];
      orderedDps[0][1]=orderedDs[0][1];
      orderedDps[1][0]=orderedDs[0][0];
      orderedDps[1][1]=orderedDs[1][1];
      TLorentzVector ptmp = pZ1p;
      pZ1p = pZ2p;
      pZ2p = ptmp;
    }
    if (std::abs(pZ1p.M() - PDGHelpers::HVVmass)<std::abs(pZ1.M() - PDGHelpers::HVVmass) || (std::abs(pZ1p.M() - PDGHelpers::HVVmass)==std::abs(pZ1.M() - PDGHelpers::HVVmass) && pZ2p.Pt()>pZ2.Pt()) ){
      for (int i=0; i<2; i++){
        for (int j=0; j<2; j++) orderedDs[i][j] = orderedDps[i][j];
      }
    }
  }
  sortedDaughters.clear();
  for (int i=0; i<2; i++){
    for (int j=0; j<2; j++){
      if (orderedDs[i][j]!=0) sortedDaughters.push_back(orderedDs[i][j]);
    }
  }
}
void ZZCandidate::createSortedVs(){
  int VID = 23;
  if (PDGHelpers::HVVmass==PDGHelpers::Wmass) VID = 24;

  TLorentzVector pZ1(0, 0, 0, 0);
  TLorentzVector pZ2(0, 0, 0, 0);
  for (int d=0; d<2; d++){
    if (sortedDaughters.at(d)!=0) pZ1 = pZ1 + sortedDaughters.at(d)->p4;
  }
  for (int d=2; d<4; d++){
    if (sortedDaughters.at(d)!=0) pZ2 = pZ2 + sortedDaughters.at(d)->p4;
  }

  // If the number of Zs is less than 2, should still create empty particles
  Particle* Z1 = new Particle(VID, pZ1);
  Z1->addMother(this);
  for (int d=0; d<2; d++){
    if (sortedDaughters.at(d)!=0) Z1->addDaughter(sortedDaughters.at(d));
  }
  addSortedV(Z1);

  Particle* Z2 = new Particle(VID, pZ2);
  Z2->addMother(this);
  for (int d=2; d<4; d++){
    if (sortedDaughters.at(d)!=0) Z2->addDaughter(sortedDaughters.at(d));
  }
  addSortedV(Z2);
}
TLorentzVector ZZCandidate::getAlternativeVMomentum(int index)const{
  TLorentzVector nullFourVector(0, 0, 0, 0);
  if (sortedDaughters.size()>3){
    TLorentzVector pZ1 = sortedDaughters.at(0)->p4;
    if (sortedDaughters.size()>=3) pZ1 = pZ1 + sortedDaughters.at(3)->p4;
    TLorentzVector pZ2 = sortedDaughters.at(2)->p4+sortedDaughters.at(1)->p4;
    if (std::abs(pZ1.M() - PDGHelpers::HVVmass)>std::abs(pZ2.M() - PDGHelpers::HVVmass)){
      TLorentzVector pZtmp = pZ1;
      pZ1 = pZ2;
      pZ2 = pZtmp;
    }
    return (index==0 ? pZ1 : pZ2);
  }
  else if (sortedDaughters.size()<3) return nullFourVector;
}

bool ZZCandidate::checkDaughtership(Particle* myParticle)const{
  for (int dd=0; dd<getNDaughters(); dd++){
    if (myParticle==getDaughter(dd)) return true;
  }
  return false;
}

void ZZCandidate::addAssociatedLeptons(Particle* myParticle){
  if (!checkDaughtership(myParticle)) addByHighestPt(myParticle, associatedLeptons);
}
void ZZCandidate::addAssociatedNeutrinos(Particle* myParticle){
  if (!checkDaughtership(myParticle)){
    addByHighestPt(myParticle, associatedLeptons); // Neutrinos are leptons at the ZZ candidate level
    addByHighestPt(myParticle, associatedNeutrinos);
  }
}
void ZZCandidate::addAssociatedJets(Particle* myParticle){
  if (!checkDaughtership(myParticle)) addByHighestPt(myParticle, associatedJets);
}
void ZZCandidate::addByHighestPt(Particle* myParticle, std::vector<Particle*>& particleArray){
  bool inserted = checkParticleExists(myParticle, particleArray); // Test if the particle is already in the vector
  if (!inserted){
    for (std::vector<Particle*>::iterator it = particleArray.begin(); it<particleArray.end(); it++){
      if ((*it)->pt()<myParticle->pt()){
        inserted=true;
        particleArray.insert(it, myParticle);
        break;
      }
    }
    if (!inserted) particleArray.push_back(myParticle);
  }
}
void ZZCandidate::addAssociatedVs(){
  createAssociatedVs(associatedJets);
  createAssociatedVs(associatedLeptons);
}
void ZZCandidate::createAssociatedVs(std::vector<Particle*>& particleArray){
  for (int i = 0; i<particleArray.size(); i++){
    double Qi = particleArray.at(i)->charge();
    int id_i = particleArray.at(i)->id;

    for (int j = 1; j<particleArray.size(); j++){
      if (j<=i) continue;
      double Qj = particleArray.at(j)->charge();
      int id_j = particleArray.at(j)->id;

      int bosonId=-1;
      if ((Qi+Qj)==0 && (id_i+id_j)==0) bosonId = (id_i==0 ? 0 : 23);
      else if (
          abs(Qi+Qj)==1 // W boson
          &&
          (
            ((PDGHelpers::isALepton(id_i) || PDGHelpers::isALepton(id_j)) && std::abs(id_i+id_j)==1) // Require SF lnu particle-antiparticle pairs
          ||
            (PDGHelpers::isUpTypeQuark(id_i) && PDGHelpers::isDownTypeQuark(id_j)) // Require ud- or du-type pairs, qqbar requirement is satisfied with charge.
          ||
            (PDGHelpers::isDownTypeQuark(id_i) && PDGHelpers::isUpTypeQuark(id_j))
          )
        ) bosonId=24*(Qi+Qj);

      if (bosonId!=-1){
        TLorentzVector pV = particleArray.at(i)->p4+particleArray.at(j)->p4;
        Particle* boson = new Particle(bosonId, pV);
        int firstdaughter = i, seconddaughter = j;
        if (
          (particleArray.at(firstdaughter)->id<particleArray.at(seconddaughter)->id && !PDGHelpers::isAWBoson(bosonId))
          ||
          (std::abs(particleArray.at(firstdaughter)->id)<std::abs(particleArray.at(seconddaughter)->id) && PDGHelpers::isAWBoson(bosonId))
          ){
          firstdaughter = j; seconddaughter = i;
        }
        boson->addDaughter(particleArray.at(firstdaughter));
        boson->addDaughter(particleArray.at(seconddaughter));
        addSortedV(boson);
      }
    }
  }
}

void ZZCandidate::testPreSelectedDaughters(){
  for (int i = 0; i<getNDaughters(); i++){
    if (!(daughters.at(i)->passSelection)){
      passSelection=false;
      break;
    }
  }
}


