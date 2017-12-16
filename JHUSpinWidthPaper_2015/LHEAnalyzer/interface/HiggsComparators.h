#ifndef HIGGSCOMPARATORS_H
#define HIGGSCOMPARATORS_H

#include <iostream>
#include "Event.h"

namespace HiggsComparators{
  enum CandidateSelection{
    BestZ1ThenZ2ScSumPt,
    nCandidateSelections
  };

  MELACandidate* matchAHiggsToParticle(Event& ev, MELAParticle* genH);
  MELACandidate* candidateSelector(Event& ev, HiggsComparators::CandidateSelection scheme, int isZZ);
  MELACandidate* candComparator(MELACandidate* cand1, MELACandidate* cand2, HiggsComparators::CandidateSelection scheme, int isZZ);
}

#endif
