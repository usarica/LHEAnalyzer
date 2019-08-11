#include <cassert>
#include <iostream>
#include <utility>
#include <vector>
#include "TFile.h"
#include "TLorentzVector.h"
#include "TString.h"
#include "TTree.h"
#include <DataFormats/Common/interface/Wrapper.h>
#include <DataFormats/HepMCCandidate/interface/GenParticle.h>
#include <DataFormats/JetReco/interface/GenJet.h>
#include <SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h>


template<typename T> using CMSEDMWrapper = edm::Wrapper<T>;
using namespace std;


struct BasicParticle{
  int id;
  int status;
  TLorentzVector p4;
  std::vector<BasicParticle const*> mothers;
  std::vector<BasicParticle const*> daughters;

  BasicParticle(int const& id_, int const& status_, TLorentzVector const& p4_) : id(id_), status(status_), p4(p4_) {}
  BasicParticle(BasicParticle const& other) : id(other.id), status(other.status), p4(other.p4), mothers(other.mothers), daughters(other.daughters) {}
  ~BasicParticle(){}

  double px() const{ return p4.X(); }
  double py() const{ return p4.Y(); }
  double pz() const{ return p4.Z(); }
  double energy() const{ return p4.T(); }
  double mass() const{ return p4.M(); }

  void addMother(BasicParticle const* part);
  void addMother(BasicParticle const& part){ addMother(&part); }

  void addDaughter(BasicParticle const* part);
  void addDaughter(BasicParticle const& part){ addDaughter(&part); }

};

void BasicParticle::addMother(BasicParticle const* part){
  if (!part) return;
  bool found=false;
  for (BasicParticle const* pp:mothers){ if (pp==part){ found=true; break; } }
  if (!found) mothers.push_back(part);
}
void BasicParticle::addDaughter(BasicParticle const* part){
  if (!part) return;
  bool found=false;
  for (BasicParticle const* pp:daughters){ if (pp==part){ found=true; break; } }
  if (!found) daughters.push_back(part);
}


void trimPythia(TString cinput, TString outdir, int pythiaStep, TString jetAlgorithm){
  TString coutput = outdir;
  coutput.Append("pythiaTemp.root");

  TFile* ftemp = TFile::Open(coutput, "recreate");
  TTree* tmpTree = new TTree("TrimmedTree", "");

  vector<double> geneventinfoweights;
  tmpTree->Branch("genWeights", &geneventinfoweights);

  vector<float> GenParticles_FV[4];
  vector<int> GenParticles_id;
  vector<int> GenParticles_status;
  tmpTree->Branch("GenParticles_X", GenParticles_FV);
  tmpTree->Branch("GenParticles_Y", GenParticles_FV+1);
  tmpTree->Branch("GenParticles_Z", GenParticles_FV+2);
  tmpTree->Branch("GenParticles_E", GenParticles_FV+3);
  tmpTree->Branch("GenParticles_id", &GenParticles_id);
  tmpTree->Branch("GenParticles_status", &GenParticles_status);

  vector<float> FinalParticles_FV[4];
  vector<int> FinalParticles_id;
  vector<int> FinalParticles_status;
  tmpTree->Branch("FinalParticles_X", FinalParticles_FV);
  tmpTree->Branch("FinalParticles_Y", FinalParticles_FV+1);
  tmpTree->Branch("FinalParticles_Z", FinalParticles_FV+2);
  tmpTree->Branch("FinalParticles_E", FinalParticles_FV+3);
  tmpTree->Branch("FinalParticles_id", &FinalParticles_id);
  tmpTree->Branch("FinalParticles_status", &FinalParticles_status);

  vector<float> GenJets_FV[4];
  vector<int> GenJets_id;
  vector<int> GenJets_status;
  tmpTree->Branch("GenJets_X", GenJets_FV);
  tmpTree->Branch("GenJets_Y", GenJets_FV+1);
  tmpTree->Branch("GenJets_Z", GenJets_FV+2);
  tmpTree->Branch("GenJets_E", GenJets_FV+3);
  tmpTree->Branch("GenJets_id", &GenJets_id);
  tmpTree->Branch("GenJets_status", &GenJets_status);

  TFile* f=nullptr; f = TFile::Open(cinput, "read");
  if (f && f->IsOpen() && !f->IsZombie()){
    TTree* events = (TTree*) f->Get("Events");

    CMSEDMWrapper< vector<reco::GenJet> >* jetWrapper = nullptr;
    CMSEDMWrapper< vector<reco::GenParticle> >* genparticleWrapper = nullptr;
    CMSEDMWrapper< GenEventInfoProduct >* geneventinfoWrapper = nullptr;

    TString suffix;
    if (pythiaStep == 1) suffix = "SIM";
    else if (pythiaStep == 0) suffix = "GEN";
    else if (pythiaStep == 2) suffix = "PAT";
    else{ cout << "trimPythia should not be called with pythiaStep=" << pythiaStep << endl; assert(0); }

    events->SetBranchStatus("*", 0);
    if (pythiaStep<2){
      events->SetBranchStatus("recoGenJets_"+jetAlgorithm+"GenJets__"+suffix+"*", 1);
      events->SetBranchAddress("recoGenJets_"+jetAlgorithm+"GenJets__"+suffix+".", &jetWrapper);

      events->SetBranchStatus("recoGenParticles_genParticles__"+suffix+"*", 1);
      events->SetBranchAddress("recoGenParticles_genParticles__"+suffix+".", &genparticleWrapper);

      events->SetBranchStatus("GenEventInfoProduct_generator__"+suffix+"*", 1);
      events->SetBranchAddress("GenEventInfoProduct_generator__"+suffix+".", &geneventinfoWrapper);
    }
    else if (pythiaStep==2){
      if (jetAlgorithm.Contains("ak4")){
        events->SetBranchStatus("recoGenJets_slimmedGenJets__"+suffix+"*", 1);
        events->SetBranchAddress("recoGenJets_slimmedGenJets__"+suffix+".", &jetWrapper);
      }
      else if (jetAlgorithm.Contains("ak8")){
        events->SetBranchStatus("recoGenJets_slimmedGenJetsAK8__"+suffix+"*", 1);
        events->SetBranchAddress("recoGenJets_slimmedGenJetsAK8__"+suffix+".", &jetWrapper);
      }

      events->SetBranchStatus("recoGenParticles_prunedGenParticles__"+suffix+"*", 1);
      events->SetBranchAddress("recoGenParticles_prunedGenParticles__"+suffix+".", &genparticleWrapper);

      events->SetBranchStatus("GenEventInfoProduct_generator__SIM*", 1);
      events->SetBranchAddress("GenEventInfoProduct_generator__SIM.", &geneventinfoWrapper);
    }

    for (int ev=0; ev<events->GetEntries(); ev++){
      events->GetEntry(ev);
      if (ev%10000 == 0) cout << "Event " << ev << '/' << events->GetEntries() << "..." << endl;
      //if (ev>0) break;

      geneventinfoweights.clear();

      for (int v=0; v<4; v++){
        GenParticles_FV[v].clear();
        FinalParticles_FV[v].clear();
        GenJets_FV[v].clear();
      }
      GenParticles_id.clear();
      GenParticles_status.clear();
      FinalParticles_id.clear();
      FinalParticles_status.clear();
      GenJets_id.clear();
      GenJets_status.clear();

      GenEventInfoProduct const* geneventinfo = geneventinfoWrapper->product();
      if (geneventinfo){
        //gen::PdfInfo const* genpdf = geneventinfo->pdf();
        geneventinfoweights = geneventinfo->weights();
      }

      vector<BasicParticle> hardlist;
      vector<BasicParticle> finallist;
      vector<reco::GenParticle> const* particles = genparticleWrapper->product();
      if (particles){
        hardlist.reserve(particles->size());
        finallist.reserve(particles->size());
        for (reco::GenParticle const& part:(*particles)){
          int i_st = part.status();
          int i_id = part.pdgId();
          int abs_id = std::abs(i_id);
          bool isHardProcess = part.isHardProcess();
          bool isPromptFinal = part.isPromptFinalState() || part.isDirectPromptTauDecayProductFinalState();
          bool isFinal = (
            i_st==1
            &&
            ((abs_id>=11 && abs_id<=16) || (i_id==22 && isPromptFinal))
            );
          if (isHardProcess){
            if (i_st==2) i_st=1;

            hardlist.emplace_back(
              i_id, i_st,
              TLorentzVector(part.px(), part.py(), part.pz(), part.energy())
            );
          }
          if (isFinal) finallist.emplace_back(
            i_id, i_st,
            TLorentzVector(part.px(), part.py(), part.pz(), part.energy())
          );
        }
      }

      for (BasicParticle const& part:hardlist){
        GenParticles_FV[0].push_back(part.px());
        GenParticles_FV[1].push_back(part.py());
        GenParticles_FV[2].push_back(part.pz());
        GenParticles_FV[3].push_back(part.energy());
        GenParticles_id.push_back(part.id);
        GenParticles_status.push_back(part.status);
      }
      for (BasicParticle const& part:finallist){
        FinalParticles_FV[0].push_back(part.px());
        FinalParticles_FV[1].push_back(part.py());
        FinalParticles_FV[2].push_back(part.pz());
        FinalParticles_FV[3].push_back(part.energy());
        FinalParticles_id.push_back(part.id);
        FinalParticles_status.push_back(part.status);
      }

      vector<reco::GenJet> const* jets = jetWrapper->product();
      if (jets){
        for (reco::GenJet const& jet:(*jets)){
          GenJets_FV[0].push_back(jet.px());
          GenJets_FV[1].push_back(jet.py());
          GenJets_FV[2].push_back(jet.pz());
          GenJets_FV[3].push_back(jet.energy());
          GenJets_id.push_back(jet.pdgId());
          GenJets_status.push_back(jet.status());
        }
      }

      tmpTree->Fill();
    }
    f->Close();
  }
  else if (f && f->IsOpen()) f->Close();

  ftemp->WriteTObject(tmpTree); delete tmpTree;
  ftemp->Close();
}

int main(int argc, char** argv){
  TString cinput;
  TString outdir;
  int pythiaStep=-1;
  TString jetAlgorithm;

  if (argc!=5){
    cerr << "trimPythia: The argument count is " << argc << " != 5." << endl;
    assert(0);
  }
  for (int a=0; a<argc; a++){
    switch (a){
    case 0:
      break;
    case 1:
      cinput = argv[a];
      break;
    case 2:
      outdir = argv[a];
      break;
    case 3:
      pythiaStep = atoi(argv[a]);
      break;
    case 4:
      jetAlgorithm = argv[a];
      break;
    default:
      cerr << "trimPythia: Arguments " << a << " is not recognized." << endl;
      assert(0);
      break;
    }
  }

  trimPythia(cinput, outdir, pythiaStep, jetAlgorithm);
}
