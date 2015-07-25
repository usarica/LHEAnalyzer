#ifndef CONVERT_PYTHIA_H
#define CONVERT_PYTHIA_H

#include "converter.h"
#include <cstdlib>
#include <sstream>
#include <utility>
//#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
//#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
//#include "DataFormats/Common/interface/Wrapper.h"
//#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
//#include "DataFormats/VertexReco/interface/Vertex.h"
//#include "DataFormats/VertexReco/interface/VertexFwd.h"
//#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
//#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
//#include "DataFormats/JetReco/interface/GenJet.h"
//#include "DataFormats/JetReco/interface/GenJetCollection.h"
//#include "DataFormats/METReco/interface/GenMET.h"
//#include "DataFormats/METReco/interface/MET.h"
//#include "DataFormats/PatCandidates/interface/GenericParticle.h"
//#include "DataFormats/PatCandidates/interface/Muon.h"
//#include "DataFormats/PatCandidates/interface/Electron.h"
//#include "DataFormats/PatCandidates/interface/Jet.h"
//#include "DataFormats/PatCandidates/interface/Photon.h"
//#include "DataFormats/PatCandidates/interface/MET.h"
//#include "DataFormats/PatCandidates/interface/Tau.h"
//#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"

class convertPythia : public converter{
public:
  convertPythia(OptionParser* options_);
  ~convertPythia(){};
  void run();

protected:
  void configure(); // Set output file, tree
  void finalizeRun();
  void readEvent(TTree* tin, int ev, vector<Particle*>& genCollection, bool& genSuccess, vector<Particle*>& recoCollection, bool& smearedSuccess, double& eventWeight);

  TFile* getIntermediateFile(string cinput);
};
#endif
