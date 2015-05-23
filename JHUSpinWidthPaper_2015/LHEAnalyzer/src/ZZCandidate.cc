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
Particle* ZZCandidate::getAssociatedJet(int index)const{
  if ((int)associatedJets.size()>index) return associatedJets.at(index);
  else return 0;
}
void ZZCandidate::sortDaughtersInitial(){
  int tmpDindex[2]={ 0 };
  Particle* df[2] = { getDaughter(0), 0 };
  for (int j=1; j<getNDaughters(); j++){
    Particle* dtmp = getDaughter(j);
    if (dtmp->charge()==-df[0]->charge()/* && (df[0]->getMother(0)==dtmp->getMother(0) || df[0]->getMother(0)==dtmp->getMother(1))*/){
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
  if (df[0]->charge()<df[1]->charge()){
    Particle* dtmp = df[0];
    df[0] = df[1];
    df[1] = dtmp;
  }
  if (ds[0]->charge()<ds[1]->charge()){
    Particle* dtmp = ds[0];
    ds[0] = ds[1];
    ds[1] = dtmp;
  }
  for (int i=0; i<2; i++) sortedDaughters.push_back(df[i]);
  for (int i=0; i<2; i++) sortedDaughters.push_back(ds[i]);
}
void ZZCandidate::sortDaughtersByBestZ1(){
  Particle* orderedDs[2][2]={ { 0 } };

  TLorentzVector pZ1 = sortedDaughters.at(0)->p4+sortedDaughters.at(1)->p4;
  TLorentzVector pZ2 = sortedDaughters.at(2)->p4+sortedDaughters.at(3)->p4;
  TLorentzVector pZ1p = sortedDaughters.at(0)->p4+sortedDaughters.at(3)->p4;
  TLorentzVector pZ2p = sortedDaughters.at(2)->p4+sortedDaughters.at(1)->p4;
  if (std::abs(pZ1.M() - 91.1876)<std::abs(pZ2.M() - 91.1876)){
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
  if (orderedDs[0][1]->id == orderedDs[1][1]->id){
    Particle* orderedDps[2][2]={ { 0 } };
    TLorentzVector pZ1p = orderedDs[0][0]->p4+orderedDs[1][1]->p4;
    TLorentzVector pZ2p = orderedDs[1][0]->p4+orderedDs[0][1]->p4;
    if (std::abs(pZ1p.M() - 91.1876)<std::abs(pZ2p.M() - 91.1876)){
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
    if (std::abs(pZ1p.M() - 91.1876)<std::abs(pZ1.M() - 91.1876) || (std::abs(pZ1p.M() - 91.1876)==std::abs(pZ1.M() - 91.1876) && pZ2p.Pt()>pZ2.Pt()) ){
      for (int i=0; i<2; i++){
        for (int j=0; j<2; j++) orderedDs[i][j] = orderedDps[i][j];
      }
    }
  }
  sortedDaughters.clear();
  for (int i=0; i<2; i++){
    for (int j=0; j<2; j++) sortedDaughters.push_back(orderedDs[i][j]);
  }
}
void ZZCandidate::createSortedVs(){
  TLorentzVector pZ1 = sortedDaughters.at(0)->p4+sortedDaughters.at(1)->p4;
  Particle* Z1 = new Particle(23, pZ1);
  Z1->addDaughter(sortedDaughters.at(0));
  Z1->addDaughter(sortedDaughters.at(1));
  addSortedV(Z1);

  TLorentzVector pZ2 = sortedDaughters.at(2)->p4+sortedDaughters.at(3)->p4;
  Particle* Z2 = new Particle(23, pZ2);
  Z2->addDaughter(sortedDaughters.at(2));
  Z2->addDaughter(sortedDaughters.at(3));
  addSortedV(Z2);
}
void ZZCandidate::addAssociatedLeptons(Particle* myParticle){
  addByHighestPt(myParticle, associatedLeptons);
}
void ZZCandidate::addAssociatedJets(Particle* myParticle){
  addByHighestPt(myParticle, associatedJets);
}
void ZZCandidate::addByHighestPt(Particle* myParticle, std::vector<Particle*>& particleArray){
  if (particleArray.size()==0) particleArray.push_back(myParticle);
  else{
    for (std::vector<Particle*>::iterator it = particleArray.begin(); it<particleArray.begin(); it++){
      if ((*it)->pt()<myParticle->pt()){
        particleArray.insert(it, myParticle);
        break;
      }
    }
  }
}


