#include <iostream>
#include <utility>
#include <vector>
#include "TFile.h"
#include "TLorentzVector.h"
#include "TString.h"
#include "TTree.h"
#include <DataFormats/HepMCCandidate/interface/GenParticle.h>
#include <DataFormats/JetReco/interface/GenJet.h>
#include <FWCore/Utilities/interface/TypeWithDict.h>
#include <SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h>
#include <DataFormats/Common/interface/Wrapper.h>


template<typename T> using CMSEDMWrapper = edm::Wrapper<T>;
using namespace std;


void trimPythia(TString cinput, TString outdir, int pythiaStep, TString jetAlgorithm){
  TString coutput = outdir;
  coutput.Append("pythiaTemp.root");

  TFile ftemp(coutput, "recreate");
  TTree* tmpTree = new TTree("tmpTree", "");

  vector<double> geneventinfoweights;

  vector<double> reco_GenJet_FV[4];
  vector<int> reco_GenJet_id;
  vector<int> reco_GenJet_status;

  vector<double> reco_GenParticle_FV[4];
  vector<int> reco_GenParticle_id;
  vector<int> reco_GenParticle_status;

  tmpTree->Branch("genWeights", &geneventinfoweights);

  tmpTree->Branch("reco_GenJet_X", reco_GenJet_FV);
  tmpTree->Branch("reco_GenJet_Y", reco_GenJet_FV+1);
  tmpTree->Branch("reco_GenJet_Z", reco_GenJet_FV+2);
  tmpTree->Branch("reco_GenJet_E", reco_GenJet_FV+3);
  tmpTree->Branch("reco_GenJet_id", &reco_GenJet_id);
  tmpTree->Branch("reco_GenJet_status", &reco_GenJet_status);

  tmpTree->Branch("reco_GenParticle_X", reco_GenParticle_FV);
  tmpTree->Branch("reco_GenParticle_Y", reco_GenParticle_FV+1);
  tmpTree->Branch("reco_GenParticle_Z", reco_GenParticle_FV+2);
  tmpTree->Branch("reco_GenParticle_E", reco_GenParticle_FV+3);
  tmpTree->Branch("reco_GenParticle_id", &reco_GenParticle_id);
  tmpTree->Branch("reco_GenParticle_status", &reco_GenParticle_status);

  TFile f(cinput, "read");
  if (f.IsOpen() && !f.IsZombie()){
    TTree* events = (TTree*) f.Get("Events");
    //events->SetBranchStatus("*", 0);

    CMSEDMWrapper< vector<reco::GenJet> >* jetWrapper = nullptr;
    CMSEDMWrapper< vector<reco::GenParticle> >* genparticleWrapper = nullptr;
    CMSEDMWrapper< GenEventInfoProduct >* geneventinfoWrapper = nullptr;

    TString suffix;
    if (pythiaStep == 1) suffix = "SIM";
    else if (pythiaStep == 0) suffix = "GEN";
    else{ cout << "trimPythia should not be called with pythiaStep=" << pythiaStep << endl; assert(0); }

    events->SetBranchStatus("recoGenJets_"+jetAlgorithm+"GenJets__"+suffix+"*", 1);
    events->SetBranchStatus("recoGenParticles_genParticles__"+suffix+"*", 1);
    events->SetBranchStatus("GenEventInfoProduct_generator__"+suffix+"*", 1);

    events->SetBranchAddress("recoGenJets_"+jetAlgorithm+"GenJets__"+suffix+".", &jetWrapper);
    events->SetBranchAddress("recoGenParticles_genParticles__"+suffix+".", &genparticleWrapper);
    events->SetBranchAddress("GenEventInfoProduct_generator__"+suffix+".", &geneventinfoWrapper);

    for (int ev=0; ev<events->GetEntries(); ev++){
      events->GetEntry(ev);

      geneventinfoweights.clear();

      for (int v=0; v<4; v++){
        reco_GenJet_FV[v].clear();
        reco_GenParticle_FV[v].clear();
      }
      reco_GenJet_id.clear();
      reco_GenParticle_id.clear();
      reco_GenJet_status.clear();
      reco_GenParticle_status.clear();

      GenEventInfoProduct const* geneventinfo = geneventinfoWrapper->product();
      if (geneventinfo){
        //gen::PdfInfo const* genpdf = geneventinfo->pdf();
        geneventinfoweights = geneventinfo->weights();
      }

      vector<reco::GenParticle> const* particles = genparticleWrapper->product();
      if (particles){
        for (reco::GenParticle const& part:(*particles)){
          if (
            part.status()==21 || part.status()==23 || part.status()==22 || part.status()==1 || part.status()==2
            ){
            reco_GenParticle_FV[0].push_back(part.px());
            reco_GenParticle_FV[1].push_back(part.py());
            reco_GenParticle_FV[2].push_back(part.pz());
            reco_GenParticle_FV[3].push_back(part.energy());
            reco_GenParticle_id.push_back(part.pdgId());
            reco_GenParticle_status.push_back(part.status());
            cout << "\t- Particle " << reco_GenParticle_id.back() << " with status " << reco_GenParticle_status.back() << endl;
          }
        }
      }

      vector<reco::GenJet> const* jets = jetWrapper->product();
      if (jets){
        for (reco::GenJet const& jet:(*jets)){
          reco_GenJet_FV[0].push_back(jet.px());
          reco_GenJet_FV[1].push_back(jet.py());
          reco_GenJet_FV[2].push_back(jet.pz());
          reco_GenJet_FV[3].push_back(jet.energy());
          reco_GenJet_id.push_back(jet.pdgId());
          reco_GenJet_status.push_back(jet.status());
        }
      }

      tmpTree->Fill();
    }
    f.Close();
  }
  else if (f.IsOpen()) f.Close();

  ftemp.WriteTObject(tmpTree); delete tmpTree;
  ftemp.Close();
}
