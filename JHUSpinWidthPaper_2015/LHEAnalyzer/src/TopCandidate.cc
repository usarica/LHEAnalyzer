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
  if (lightQuark!=0){
    p4 = p4 + lightQuark->p4;
    addDaughter(lightQuark);
  }
  double Wcharge=0;
  if (Wferm!=0){
    p4 = p4 + Wferm->p4;
    Wcharge += Wferm->charge();
    addDaughter(Wferm);
  }
  if (Wfermbar!=0){
    p4 = p4 + Wfermbar->p4;
    Wcharge += Wfermbar->charge();
    addDaughter(Wfermbar);
  }
  if (
    (lightQuark!=0 && PDGHelpers::isDownTypeQuark(lightQuark->id) && lightQuark->id>0)
    ||
    (Wcharge>0)
    ) id = 6;
  else if (
    (lightQuark!=0 && PDGHelpers::isDownTypeQuark(lightQuark->id) && lightQuark->id<0)
    ||
    (Wcharge<0)
    ) id = -6;
  else id=0;
}

Particle* TopCandidate::setLightQuark(Particle* myParticle){ lightQuark=myParticle; if (lightQuark!=0) addDaughter(lightQuark); }
Particle* TopCandidate::setWFermion(Particle* myParticle){ Wferm=myParticle; if (Wferm!=0) addDaughter(Wferm); }
Particle* TopCandidate::setWAntifermion(Particle* myParticle){ Wfermbar=myParticle; if (Wfermbar!=0) addDaughter(Wfermbar); }

