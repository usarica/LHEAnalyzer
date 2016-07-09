#include "../interface/HiggsComparators.h"

ZZCandidate* HiggsComparators::matchAHiggsToParticle(Event& ev, Particle* genH){
  ZZCandidate* cand=0;
  for (int t=0; t<ev.getNZZCandidates(); t++){
    ZZCandidate* tmpCand = ev.getZZCandidate(t);
    double genhmassquant = genH->m()+genH->pt()+fabs(genH->z());
    double massdiff = fabs(genhmassquant-tmpCand->m()-tmpCand->pt()-fabs(tmpCand->z()));
    double massratio = 0;
    if (genhmassquant>0) massratio = massdiff / genhmassquant;
    if (massratio<0.001){
      if (cand==0) cand = tmpCand;
      else{
        TLorentzVector vGen = genH->p4;
        TLorentzVector vTmp = tmpCand->p4;
        TLorentzVector vCur = cand->p4;

        double dot_tmp = vTmp.Dot(vGen);
        double dot_curr = vCur.Dot(vGen);
        if (fabs(dot_tmp-vGen.M2())<fabs(dot_curr - vGen.M2())) cand = tmpCand;
      }
    }
  }
  return cand;
}

ZZCandidate* HiggsComparators::candidateSelector(Event& ev, HiggsComparators::CandidateSelection scheme, int isZZ){
  ZZCandidate* cand=0;
  for (int t=0; t<ev.getNZZCandidates(); t++){
    ZZCandidate* tmpCand = ev.getZZCandidate(t);
    if (!tmpCand->passSelection) continue;
    if (cand==0) cand=tmpCand;
    else cand = HiggsComparators::candComparator(cand, tmpCand, scheme, isZZ);
  }
  return cand;
}

ZZCandidate* HiggsComparators::candComparator(ZZCandidate* cand1, ZZCandidate* cand2, HiggsComparators::CandidateSelection scheme, int isZZ){
  ZZCandidate* theChosenOne=0;

  double defaultHVVmass = PDGHelpers::HVVmass;
  if (isZZ==0){
    PDGHelpers::setHVVmass(PDGHelpers::Wmass);
  }
  else if (isZZ==1){
    PDGHelpers::setHVVmass(PDGHelpers::Zmass);
  }

  if (isZZ==-1){
    if (
      (
      (cand2->m()>cand1->m())
      ||
      (cand2->m()==cand1->m() && cand2->pt()>cand1->pt())
      ) && PDGHelpers::isAHiggs(cand2->id)
      ) theChosenOne = cand2;
    else if (PDGHelpers::isAHiggs(cand1->id)) theChosenOne = cand1;
  }
  else if (scheme==HiggsComparators::BestZ1ThenZ2ScSumPt){
    double diffmass1 = fabs(cand1->getSortedV(0)->m()-PDGHelpers::HVVmass);
    double diffmass2 = fabs(cand2->getSortedV(0)->m()-PDGHelpers::HVVmass);
    double Z2scsumpt_cand1=0, Z2scsumpt_cand2=0;
    Particle* c11 = cand1->getSortedV(1)->getDaughter(0);
    Particle* c12 = cand1->getSortedV(1)->getDaughter(1);
    Particle* c21 = cand2->getSortedV(1)->getDaughter(0);
    Particle* c22 = cand2->getSortedV(1)->getDaughter(1);
    if (c11!=0) Z2scsumpt_cand1 += c11->pt();
    if (c12!=0) Z2scsumpt_cand1 += c12->pt();
    if (c21!=0) Z2scsumpt_cand2 += c21->pt();
    if (c22!=0) Z2scsumpt_cand2 += c22->pt();
    if (
      (diffmass1>diffmass2)
      ||
      (diffmass1==diffmass2 && Z2scsumpt_cand2>Z2scsumpt_cand1)
      ) theChosenOne = cand2;
    else theChosenOne = cand1;
  }

  PDGHelpers::setHVVmass(defaultHVVmass);
  return theChosenOne;
}


