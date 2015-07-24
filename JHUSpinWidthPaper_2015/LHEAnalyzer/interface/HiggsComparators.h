#ifndef HIGGSCOMPARATORS_H
#define HIGGSCOMPARATORS_H

#include <iostream>
#include "Event.h"

namespace HiggsComparators{
  enum CandidateSelection{
    BestZ1ThenZ2ScSumPt,
    nCandidateSelections
  };

  ZZCandidate* matchAHiggsToParticle(Event& ev, Particle* genH);
  ZZCandidate* candidateSelector(Event& ev, HiggsComparators::CandidateSelection scheme, bool isZZ);
  ZZCandidate* candComparator(ZZCandidate* cand1, ZZCandidate* cand2, HiggsComparators::CandidateSelection scheme, bool isZZ);
}

#endif
