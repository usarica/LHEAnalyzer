#include <iostream>
#include "../interface/Event.h"

using namespace PDGHelpers;
using namespace ParticleComparators;

void Event::applyParticleSelection(){
  applyLeptonSelection();
  applyNeutrinoSelection();
  applyJetSelection();
  applyZZSelection(); // Order matters here
}
void Event::applyLeptonSelection(){
  for (std::vector<Particle*>::iterator it = leptons.begin(); it<leptons.end(); it++){
    // Trigger and acceptance
    bool passAcceptance = true;
    if (std::abs((*it)->id)==11 && ((*it)->pt()<=electronPTCut || std::abs((*it)->eta())>=electronEtaAcceptanceCut)) passAcceptance = false;
    else if (std::abs((*it)->id)==13 && ((*it)->pt()<=muonPTCut || std::abs((*it)->eta())>=muonEtaAcceptanceCut)) passAcceptance = false;
    else if (std::abs((*it)->id)==15) passAcceptance = false;
    for (std::vector<Particle*>::iterator it2 = leptons.begin(); it2<leptons.end(); it2++){
      if ((*it2)==(*it)) continue; // Every particle is their own ghost.
      else if ((*it)->deltaR((*it2)->p4)<=ghostDeltaRCut) passAcceptance = false; // Ghost removal
    }
    (*it)->setSelected(passAcceptance);
  }
}
void Event::applyNeutrinoSelection(){
  for (std::vector<Particle*>::iterator it = neutrinos.begin(); it<neutrinos.end(); it++) (*it)->setSelected(false);
}
void Event::applyJetSelection(){
  for (std::vector<Particle*>::iterator it = jets.begin(); it<jets.end(); it++){
    bool passAcceptance = true;
    if ((*it)->pt()<=jetPTCut || std::abs((*it)->eta())>=jetEtaAcceptanceCut) passAcceptance = false; // ZZ4l selection and acceptance
    for (std::vector<Particle*>::iterator it2 = leptons.begin(); it2<leptons.end(); it2++){ // Clean from selected leptons
      if ((*it2)->passSelection){ // If it is not selected at all, why would I care?
        if ((*it)->deltaR((*it2)->p4)<=jetDeltaR) passAcceptance = false;
      }
    }
    (*it)->setSelected(passAcceptance);
  }
}
void Event::applyZZSelection(){
  for (std::vector<ZZCandidate*>::iterator it = ZZcandidates.begin(); it<ZZcandidates.end(); it++){
    (*it)->testPreSelectedDaughters();
    if (!(*it)->passSelection) continue;

    bool passAcceptance = true;
    if ((*it)->getSortedV(0)->m()<=mV1LowCut || (*it)->getSortedV(0)->m()>=mV12HighCut){
      passAcceptance = false; (*it)->getSortedV(0)->setSelected(passAcceptance);
    } // Z1 selection
    if ((*it)->getSortedV(1)->m()<=mV2LowCut || (*it)->getSortedV(1)->m()>=mV12HighCut){
      passAcceptance = false; (*it)->getSortedV(1)->setSelected(passAcceptance);
    } // Z2 selection
    for (int iZ=2; iZ<(*it)->getNSortedVs(); iZ++){
      Particle* extraV = (*it)->getSortedV(iZ);
      if (!isAZBoson(extraV->id)) continue;
      else{
        if (extraV->m()<=mllLowCut || extraV->m()>=mV12HighCut || (extraV->getDaughter(0)!=0 && isANeutrino(extraV->getDaughter(0)->id))) extraV->setSelected(false); // Extra Z selection, no effect on ZZ candidate
      }
    }
    TLorentzVector pLOC[2];
    pLOC[0]=(*it)->getAlternativeVMomentum(0);
    pLOC[1]=(*it)->getAlternativeVMomentum(1);
    if (pLOC[0].M()<=mllLowCut || pLOC[1].M()<=mllLowCut) passAcceptance=false;

    (*it)->setSelected(passAcceptance);
  }
}

void Event::addZZCandidate(ZZCandidate* myParticle){
  bool isIdentical = (getNZZCandidates()>0);
  for (int cc=0; cc<getNZZCandidates(); cc++){
    ZZCandidate* testCand = ZZcandidates.at(cc);

    for (int i=0; i<2; i++){
      Particle* testV = testCand->getSortedV(i);
      Particle* partV = myParticle->getSortedV(i);
      if (partV==0 || testV==0){
        if (partV==0 && testV==0) continue; // Check the next intermediate V if there is any.
        else { isIdentical=false; break; } // They are definitely not the same.
      }
      if (testV->getNDaughters() != partV->getNDaughters()) { isIdentical=false;  break; }; // Again, they cannot be the same.
      for (int j=0; j<testV->getNDaughters(); j++){
        Particle* testD = testV->getDaughter(j);
        Particle* partD = partV->getDaughter(j);
        isIdentical = isIdentical && (testD==partD);
        if (!isIdentical) break;
      }
      if (!isIdentical) break;
    }
    if (!isIdentical) break;
  }
  if (!isIdentical) ZZcandidates.push_back(myParticle);
  else { delete myParticle; myParticle=0; }
}

void Event::constructVVCandidates(bool isZZ, int fstype){
  /*
  ZZ / WW
  fstype=0: 4l / lnulnu
  fstype=1: 4q / 4q
  fstype=2: 2l2q / lnu2q
  fstype=3: 2l2nu / -
  fstype=4: 2q2nu / -
  fstype=5: 4nu / -
  fstype=-1: Any / Any
  */

  if (!isZZ && fstype>2){
    std::cerr << "No " << (isZZ ? "ZZ" : "WW") << " candidate with final state " << fstype << " is possible!" << std::endl;
    return;
  }

  std::vector<Particle*> lepPlusMinus[3][2];
  std::vector<Particle*> lepNu[3][2];
  std::vector<Particle*> quarkPlusMinus[7][2];

  for (std::vector<Particle*>::iterator it = leptons.begin(); it<leptons.end(); it++){ // Leptons
    int iFirst=0;
    int iSecond=0;

    if (abs((*it)->id)==11) iFirst = 0;
    else if (abs((*it)->id)==13) iFirst = 1;
    else if (abs((*it)->id)==15) iFirst = 2;
    if ((*it)->id>0) iSecond=1;
    lepPlusMinus[iFirst][iSecond].push_back(*it);
  }
  for (std::vector<Particle*>::iterator it = neutrinos.begin(); it<neutrinos.end(); it++){ // Neutrinos
    int iFirst=0;
    int iSecond=0;

    if (abs((*it)->id)==12) iFirst = 0;
    else if (abs((*it)->id)==14) iFirst = 1;
    else if (abs((*it)->id)==16) iFirst = 2;
    if ((*it)->id>0) iSecond=1;
    lepNu[iFirst][iSecond].push_back(*it);
  }
  for (std::vector<Particle*>::iterator it = jets.begin(); it<jets.end(); it++){ // Jets
    int iFirst=abs((*it)->id); // Yes, 0-6, 0 being unknown
    if (PDGHelpers::isAGluon(iFirst)) continue;
    int iSecond=0;

    if ((*it)->id>0) iSecond=1;
    quarkPlusMinus[iFirst][iSecond].push_back(*it);
  }

  std::vector<Particle*> tmpVhandle;

  if (isZZ){ // ZZ

    if (fstype==-1 || fstype==0 || fstype==2 || fstype==3){ // Z->2l
      for (int c=0; c<3; c++){
        for (int i=0; i<lepPlusMinus[c][0].size(); i++){
          for (int j=0; j<lepPlusMinus[c][1].size(); j++){
            TLorentzVector pV = lepPlusMinus[c][0].at(i)->p4+lepPlusMinus[c][1].at(j)->p4;
            Particle* V = new Particle(23, pV);
            V->addDaughter(lepPlusMinus[c][0].at(i));
            V->addDaughter(lepPlusMinus[c][1].at(j));
            tmpVhandle.push_back(V);
          }
        }
      }
    }
    if (fstype==-1 || fstype==3 || fstype==4 || fstype==5){ // Z->2nu
      for (int c=0; c<3; c++){
        for (int i=0; i<lepNu[c][0].size(); i++){
          for (int j=0; j<lepNu[c][1].size(); j++){
            TLorentzVector pV = lepNu[c][0].at(i)->p4+lepNu[c][1].at(j)->p4;
            Particle* V = new Particle(23, pV);
            V->addDaughter(lepNu[c][0].at(i));
            V->addDaughter(lepNu[c][1].at(j));
            tmpVhandle.push_back(V);
          }
        }
      }
    }
    if (fstype==-1 || fstype==1 || fstype==2 || fstype==4){ // Z->2q
      for (int c=1; c<7; c++){
        for (int i=0; i<quarkPlusMinus[c][0].size(); i++){
          for (int j=0; j<quarkPlusMinus[c][1].size(); j++){
            TLorentzVector pV = quarkPlusMinus[c][0].at(i)->p4+quarkPlusMinus[c][1].at(j)->p4;
            Particle* V = new Particle(23, pV);
            V->addDaughter(quarkPlusMinus[c][0].at(i));
            V->addDaughter(quarkPlusMinus[c][1].at(j));
            tmpVhandle.push_back(V);
          }
        }
      }
    }

  }
  else{ // WW

    if (fstype==-1 || fstype==0 || fstype==2){ // W->lnu
      for (int c=0; c<3; c++){
        for (int i=0; i<lepPlusMinus[c][0].size(); i++){
          for (int j=0; j<lepNu[c][1].size(); j++){
            TLorentzVector pV = lepPlusMinus[c][0].at(i)->p4+lepNu[c][1].at(j)->p4;
            Particle* V = new Particle(24, pV);
            V->addDaughter(lepPlusMinus[c][0].at(i));
            V->addDaughter(lepNu[c][1].at(j));
            tmpVhandle.push_back(V);
          }
        }
      }
      for (int c=0; c<3; c++){
        for (int i=0; i<lepPlusMinus[c][0].size(); i++){
          for (int j=0; j<lepNu[c][1].size(); j++){
            TLorentzVector pV = lepPlusMinus[c][0].at(i)->p4+lepNu[c][1].at(j)->p4;
            Particle* V = new Particle(-24, pV);
            V->addDaughter(lepPlusMinus[c][0].at(i));
            V->addDaughter(lepNu[c][1].at(j));
            tmpVhandle.push_back(V);
          }
        }
      }
    }
    if (fstype==-1 || fstype==1 || fstype==2){ // W->2q
      for (int c=1; c<7; c++){
        for (int d=1; d<7; d++){
          if (d==c) continue;
          for (int i=0; i<quarkPlusMinus[c][0].size(); i++){
            for (int j=0; j<quarkPlusMinus[d][1].size(); j++){
              int totalcharge = quarkPlusMinus[c][0].at(i)->charge() + quarkPlusMinus[d][1].at(j)->charge();
              if (abs(totalcharge)!=1) continue;

              TLorentzVector pV = quarkPlusMinus[c][0].at(i)->p4+quarkPlusMinus[d][1].at(j)->p4;
              Particle* V = new Particle(24*totalcharge, pV);
              V->addDaughter(quarkPlusMinus[c][0].at(i));
              V->addDaughter(quarkPlusMinus[d][1].at(j));
              tmpVhandle.push_back(V);
            }
          }
        }
      }
    }

  }
  if (fstype==-1 || fstype==1 || fstype==2 || fstype==4){ // Z/W->2j reco.-level
    for (int i=0; i<quarkPlusMinus[0][0].size(); i++){
      if (quarkPlusMinus[0][0].at(i)->id!=0) continue;
      for (int j=i+1; j<quarkPlusMinus[0][0].size(); j++){
        if (quarkPlusMinus[0][0].at(j)->id!=0) continue;
        TLorentzVector pV = quarkPlusMinus[0][0].at(i)->p4+quarkPlusMinus[0][0].at(j)->p4;
        Particle* V = new Particle(0, pV);
        V->addDaughter(quarkPlusMinus[0][0].at(i));
        V->addDaughter(quarkPlusMinus[0][0].at(j));
        tmpVhandle.push_back(V);
      }
    }
  }



  for (int i=0; i<tmpVhandle.size(); i++){
    for (int j=i; j<tmpVhandle.size(); j++){
      if ((tmpVhandle.at(i)->charge()+tmpVhandle.at(j)->charge())!=0) continue;
      if (tmpVhandle.at(i)==tmpVhandle.at(j)) continue;
      Particle* Vi1 = tmpVhandle.at(i)->getDaughter(0);
      Particle* Vi2 = tmpVhandle.at(i)->getDaughter(1);
      Particle* Vj1 = tmpVhandle.at(j)->getDaughter(0);
      Particle* Vj2 = tmpVhandle.at(j)->getDaughter(1);
      if (Vi1==Vj1 || Vi2==Vj2) continue;

      TLorentzVector pH = Vi1->p4+Vi2->p4+Vj1->p4+Vj2->p4;
      ZZCandidate* cand = new ZZCandidate(25, pH);
      cand->addDaughter(Vi1);
      cand->addDaughter(Vi2);
      cand->addDaughter(Vj1);
      cand->addDaughter(Vj2);

      double defaultHVVmass = HVVmass;
      if (isZZ){
        setHVVmass(Zmass);
      }
      else{
        setHVVmass(Wmass);
      }
      cand->sortDaughters();
      setHVVmass(defaultHVVmass);

      addZZCandidate(cand);
    }
  }

  for (int i=0; i<tmpVhandle.size(); i++) delete tmpVhandle.at(i);
  tmpVhandle.clear();
}

ZZCandidate* Event::getZZCandidate(int index)const{
  if ((int)ZZcandidates.size()>index) return ZZcandidates.at(index);
  else return 0;
}
Particle* Event::getLepton(int index)const{
  if ((int)leptons.size()>index) return leptons.at(index);
  else return 0;
}
Particle* Event::getNeutrino(int index)const{
  if ((int)neutrinos.size()>index) return neutrinos.at(index);
  else return 0;
}
Particle* Event::getJet(int index)const{
  if ((int)jets.size()>index) return jets.at(index);
  else return 0;
}
Particle* Event::getParticle(int index)const{
  if ((int)particles.size()>index) return particles.at(index);
  else return 0;
}

TLorentzVector Event::missingP() const{
  TLorentzVector totalP(0, 0, 0, 0);
  for (int pp=0; pp<particles.size();pp++){
    Particle* part = getParticle(pp);
    if (part->passSelection) totalP = totalP + part->p4;
  }
  totalP.SetT(totalP.P());
  totalP.SetVect(-totalP.Vect());
  return totalP;
}

void Event::addVVCandidateMother(Particle* mother){
  for (std::vector<ZZCandidate*>::iterator it = ZZcandidates.begin(); it<ZZcandidates.end(); it++) (*it)->addMother(mother);
}
void Event::addVVCandidateAppendages(){
  for (std::vector<ZZCandidate*>::iterator it = ZZcandidates.begin(); it<ZZcandidates.end(); it++){
    for (std::vector<Particle*>::iterator iL = leptons.begin(); iL<leptons.end(); iL++){ if ((*iL)->passSelection) (*it)->addAssociatedLeptons(*iL); }
    for (std::vector<Particle*>::iterator iL = neutrinos.begin(); iL<neutrinos.end(); iL++){ if ((*iL)->passSelection) (*it)->addAssociatedNeutrinos(*iL); }
    for (std::vector<Particle*>::iterator iJ = jets.begin(); iJ<jets.end(); iJ++){ if ((*iJ)->passSelection) (*it)->addAssociatedJets(*iJ); }
    (*it)->addAssociatedVs();
  }
}





