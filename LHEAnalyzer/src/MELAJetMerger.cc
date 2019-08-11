#include <cassert>
#include <utility>
#include "MELAJetMerger.h"
#include "ParticleComparators.h"
#include "MELAStreamHelpers.hh"


using namespace std;
using namespace MELAStreamHelpers;


void MELAJetMerger::mergeParticlesToJets(std::string const& algo, std::vector<MELAParticle const*> const& inparts, std::vector<MELAParticle>& outparts){
  int ipow;
  if (algo.find("ak")!=std::string::npos) ipow=-2;
  else if (algo.find("kt")!=std::string::npos) ipow=+2;
  else if (algo.find("ca")!=std::string::npos) ipow=0;
  else{
    MELAerr << "MELAJetMerger::mergeParticlesToJets: Jet algorithm (currently " << algo << ") can only be ak, kt or ca." << endl;
    assert(0);
  }
  mergeParticlesToJets(ipow, inparts, outparts);
}
void MELAJetMerger::mergeParticlesToJets(int ipow, std::vector<MELAParticle const*> const& inparts, std::vector<MELAParticle>& outparts){
  if (inparts.empty()) return;

  double const& dRconst = ParticleComparators::jetDeltaR;

  SimpleParticleCollection_t inoutvecs; inoutvecs.reserve(inparts.size());
  SimpleParticleCollection_t removedvecs; removedvecs.reserve(inparts.size());
  for (MELAParticle const* part:inparts) inoutvecs.emplace_back(part->id, part->p4);
  mergeMomenta(ipow, dRconst, inoutvecs, removedvecs);
  if (!inoutvecs.empty() || !removedvecs.empty()) outparts.reserve(inoutvecs.size()+removedvecs.size()+outparts.size());
  for (SimpleParticle_t const& vec:inoutvecs) outparts.emplace_back(vec.first, vec.second);
  for (SimpleParticle_t const& vec:removedvecs) outparts.emplace_back(vec.first, vec.second);
}

void MELAJetMerger::mergeMomenta(int ipow, double dR, SimpleParticleCollection_t& inoutvecs, SimpleParticleCollection_t& removedvecs){
  if (inoutvecs.size()<2) return;

  //MELAout << "Begin MELAJetMerger::mergeMomenta." << endl;

  SimpleParticleCollection_t inoutvecs_old(inoutvecs);
  // Do not copy removedvecs!

  //MELAout << "\t- Particles to merge:" << endl;
  //for (auto const& part:inoutvecs_old) MELAout << part << " (address: " << &(part.second) << ")" << endl;

  // Find min. dij
  std::pair<TLorentzVector const*, TLorentzVector const*> vecpair_dmin(nullptr, nullptr);
  double dminval=-1;

  TLorentzVector const* minETvec=nullptr;
  double minETval=-1;

  for (SimpleParticleCollection_t::const_iterator it1=inoutvecs_old.cbegin(); it1!=inoutvecs_old.cend(); it1++){
    //int const& id1 = it1->first;
    TLorentzVector const& p1 = it1->second;
    double qt1 = std::pow(p1.Pt(), ipow);

    for (SimpleParticleCollection_t::const_iterator it2=it1+1; it2!=inoutvecs_old.cend(); it2++){
      //if (it1==it2) continue;

      //int const& id2 = it2->first;
      TLorentzVector const& p2 = it2->second;

      double deltaR = std::pow(p1.DeltaR(p2)/dR, 2);
      double qt2 = std::pow(p2.Pt(), ipow);

      double dval = deltaR*std::min(qt1, qt2);
      if (dminval<0. || dval<dminval){
        dminval = dval;
        vecpair_dmin.first = &p1;
        vecpair_dmin.second = &p2;
      }
    }

    if (minETval<0. || qt1<minETval){
      minETval = qt1;
      minETvec = &p1;
    }
  }

  if (dminval<=minETval){
    //MELAout << "\t- Particles identified for merging: " << vecpair_dmin << endl;
    // Merge vecpair_dmin.first and vecpair_dmin.second, remake the inoutvecs and call again
    inoutvecs.clear(); inoutvecs.reserve(inoutvecs_old.size()-1);
    for (SimpleParticle_t const& part:inoutvecs_old){
      //int const& id = part.first;
      TLorentzVector const& vec = part.second;

      if ((&vec)==vecpair_dmin.second) continue;
      else if ((&vec)==vecpair_dmin.first) inoutvecs.emplace_back(0, vec + *(vecpair_dmin.second));
      else inoutvecs.emplace_back(part);
    }

    mergeMomenta(ipow, dR, inoutvecs, removedvecs); // Recursively call until no unmerged particles remain
  }
  else if (minETvec){
    inoutvecs.clear(); inoutvecs.reserve(inoutvecs_old.size()-1);
    //MELAout << "\t- Particle to be removed as a single jet: " << minETvec << endl;
    for (SimpleParticle_t const& part:inoutvecs_old){
      int const& id = part.first;
      TLorentzVector const& vec = part.second;
      if ((&vec)==minETvec) removedvecs.emplace_back(id, vec);
      else inoutvecs.emplace_back(part);
    }
    mergeMomenta(ipow, dR, inoutvecs, removedvecs); // Recursively call until no unmerged particles remain
  }

  //MELAout << "End MELAJetMerger::mergeMomenta." << endl;
}
