#include <iostream>
#include <vector>
#include <utility>
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TLorentzVector.h"

using namespace std;

vector<pair<int, int>> findDuplicates(const vector<double>* fourvector, vector<int> id, vector<int> status);

void trimPythia(TString cinput, TString outdir="./"){
  TString coutput = outdir;
  coutput.Append("pythiaTemp.root");

  TFile* ftemp = new TFile(coutput, "recreate");
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
  if (!f.IsZombie()){
    TTree* events = (TTree*)f.Get("Events");

    edm::Wrapper< vector<reco::GenJet> >* jetWrapper;
    edm::Wrapper< vector<reco::GenParticle> >* genparticleWrapper;
    edm::Wrapper< GenEventInfoProduct >* geneventinfoWrapper;

    events->SetBranchAddress("recoGenJets_ak5GenJets__SIM.", &jetWrapper);
    events->SetBranchAddress("recoGenParticles_genParticles__SIM.", &genparticleWrapper);
    events->SetBranchAddress("GenEventInfoProduct_generator__SIM.", &geneventinfoWrapper);

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

      gen::PdfInfo* genpdf = geneventinfoWrapper->product()->pdf();
      const vector<double> genweights = geneventinfoWrapper->product()->weights();
      for (int g=0; g<genweights.size(); g++) geneventinfoweights.push_back(genweights.at(g));

      const vector<reco::GenJet>* reco_GenJets = jetWrapper->product();
      const vector<reco::GenParticle>* reco_GenParticles = genparticleWrapper->product();

//      int nlep=0;
      for (int p=0; p<reco_GenParticles->size(); p++){
        reco::GenParticle& part = reco_GenParticles->at(p);
        if (
          ((part.pdgId()==25 || part.pdgId()==32) && part.status()==22) // Generated Higgs
          ||
          ((std::abs(part.pdgId())>=11 && std::abs(part.pdgId())<=16) && (part.status()==23 || part.status()==1)) // Generated leptons
          ||
          ((std::abs(part.pdgId())<=6 || std::abs(part.pdgId())==21) && (part.status()==21 || part.status()==23 || part.status()==22)) // Generated partons
          ){
          reco_GenParticle_FV[0].push_back(part.px());
          reco_GenParticle_FV[1].push_back(part.py());
          reco_GenParticle_FV[2].push_back(part.pz());
          reco_GenParticle_FV[3].push_back(part.energy());
          reco_GenParticle_id.push_back(part.pdgId());
          reco_GenParticle_status.push_back(part.status());
          if ((std::abs(part.pdgId())>=11 && std::abs(part.pdgId())<=16) && std::abs(part.pdgId()) % 2 == 1){
//            nlep++;

            if (part.status()==1){
              reco_GenJet_FV[0].push_back(part.px());
              reco_GenJet_FV[1].push_back(part.py());
              reco_GenJet_FV[2].push_back(part.pz());
              reco_GenJet_FV[3].push_back(part.energy());
              reco_GenJet_id.push_back(part.pdgId());
              reco_GenJet_status.push_back(part.status());
            }

          }
        }
      }
//      cout << "Nleps: " << nlep << endl;
      vector<pair<int, int>> reco_GenParticle_duplicates = findDuplicates(reco_GenParticle_FV, reco_GenParticle_id, reco_GenParticle_status);
      vector<int> removalArray;
      for (int dd=0; dd<reco_GenParticle_duplicates.size(); dd++){
        pair<int, int> duplicate = reco_GenParticle_duplicates.at(dd);
        int iOriginal = duplicate.first;
        int iTransfer = duplicate.second;

        bool inserted=false;
        for (int it = 0; it<removalArray.size(); it++){
          int iIndex = removalArray.at(it);
          if (iTransfer > iIndex){
            removalArray.insert(removalArray.begin()+it, iTransfer);
            inserted=true;
            break;
          }
        }
        if (!inserted) removalArray.push_back(iTransfer);
      }
//      cout << "Nremoved: " << removalArray.size() << endl;
      for (int dd=0; dd<removalArray.size(); dd++){
        int iTransfer = removalArray.at(dd);
        for (int fv=0; fv<4; fv++){
          reco_GenParticle_FV[fv].erase(reco_GenParticle_FV[fv].begin()+iTransfer);
        }
        reco_GenParticle_id.erase(reco_GenParticle_id.begin()+iTransfer);
        reco_GenParticle_status.erase(reco_GenParticle_status.begin()+iTransfer);
      }

      for (int p=0; p<reco_GenJets->size(); p++){
        reco::GenJet& jet = reco_GenJets->at(p);
        reco_GenJet_FV[0].push_back(jet.px());
        reco_GenJet_FV[1].push_back(jet.py());
        reco_GenJet_FV[2].push_back(jet.pz());
        reco_GenJet_FV[3].push_back(jet.energy());
        reco_GenJet_id.push_back(jet.pdgId());
        reco_GenJet_status.push_back(jet.status());
      }
/*
      int nlepreco = 0;
      for (int r=0; r<reco_GenJet_id.size(); r++){
        int pid = reco_GenJet_id.at(r);
        if ((std::abs(pid)>=11 && std::abs(pid)<=16) && std::abs(pid) % 2 == 1) nlepreco++;
      }
      cout << "Nrecolep: " << nlepreco << endl;
*/
      tmpTree->Fill();
    }
    f.Close();
  }
  else if (f.IsOpen()) f.Close();
  ftemp->WriteTObject(tmpTree);
  delete tmpTree;
  ftemp->Close();
}

vector<pair<int, int>> findDuplicates(const vector<double>* fourvector, vector<int> id, vector<int> status){
  vector<pair<int, int>> duplicates;

  for (int xx=0; xx<id.size(); xx++){
    TLorentzVector p_tm(fourvector[0].at(xx), fourvector[1].at(xx), fourvector[2].at(xx), fourvector[3].at(xx));
    int id_tm = id.at(xx);
    int st_tm = status.at(xx);
    if (st_tm==23 || st_tm==1){
      for (int yy=xx+1; yy<id.size(); yy++){
        TLorentzVector p_tbm(fourvector[0].at(yy), fourvector[1].at(yy), fourvector[2].at(yy), fourvector[3].at(yy));
        int id_tbm = id.at(yy);
        int st_tbm = status.at(yy);

        if (id_tbm==id_tm && st_tbm!=st_tm && (st_tbm==1 || st_tbm==23)){
          double dot_tm = p_tbm.Dot(p_tm);
          double diff_sqmass = dot_tm - p_tm.M2();
          diff_sqmass = fabs((double)diff_sqmass);
          if (diff_sqmass<0.005*fabs(p_tm.M2())){
            int iFirst = (st_tm==23 ? xx : yy); // Should be order-independent
            int iSecond = (st_tbm==1 ? yy : xx);
            pair<int, int> dupinst(iFirst, iSecond);
            duplicates.push_back(dupinst);
          }
        }
      }
    }
  }

  return duplicates;
}

