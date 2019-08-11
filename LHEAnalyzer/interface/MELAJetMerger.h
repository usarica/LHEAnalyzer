#ifndef MELAJETMERGER_H
#define MELAJETMERGER_H

#include <string>
#include "MELAParticle.h"
#include "TVar.hh"


namespace MELAJetMerger{

  void mergeParticlesToJets(std::string const& algo, std::vector<MELAParticle const*> const& inparts, std::vector<MELAParticle>& outparts);

  void mergeParticlesToJets(int ipow, std::vector<MELAParticle const*> const& inparts, std::vector<MELAParticle>& outparts);

  void mergeMomenta(int ipow, double dR, SimpleParticleCollection_t& inoutvecs, SimpleParticleCollection_t& removedvecs);

}

#endif
