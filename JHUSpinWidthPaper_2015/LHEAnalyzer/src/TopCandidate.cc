#include "../interface/TopCandidate.h"

using namespace std;

TopCandidate::TopCandidate(
  Particle* lightQuark_,
  Particle* Wferm_,
  Particle* Wfermbar_
  ) : Particle(),
  lightQuark(lightQuark_),
  Wferm(Wferm_),
  Wfermbar(Wfermbar_)
{
  p4 = lightQuark->p4 + Wferm->p4 + Wfermbar->p4;
  if (
    (PDGHelpers::isDownTypeQuark(lightQuark->id) && lightQuark->id>0)
    ||
    ((Wferm->charge()+Wfermbar->charge())>0)
    ) id = 6;
  else if (
    (PDGHelpers::isDownTypeQuark(lightQuark->id) && lightQuark->id<0)
    ||
    ((Wferm->charge()+Wfermbar->charge())<0)
    ) id = -6;
  else id=0;
}
