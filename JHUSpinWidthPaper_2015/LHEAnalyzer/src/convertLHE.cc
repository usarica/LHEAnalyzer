#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <vector>
#include <iomanip>
#include <cstdlib>
#include "TSystem.h"
#include "TInterpreter.h"
#include "TFile.h"
#include "TList.h"
#include "TNtuple.h"
#include "TTree.h"
#include "TRandom.h"
#include "TLorentzVector.h"
#include "../interface/convertLHE.h"

using namespace PDGHelpers;
using namespace LHEParticleSmear;


int main(int argc, char ** argv){
  const int minArgsExpected=6;
  if (argc<minArgsExpected){
    cerr << "Minimum number of arguments should be " << (minArgsExpected-1) << ". Please review." << endl;
    return 0;
  }
  double mPOLE = (double)atof(argv[1]);
  int erg_tev = atoi(argv[2]);
  string cdir = argv[3];
  string coutput = argv[4];
  if (cdir.find(".lhe")!=string::npos || coutput.find(".root")==string::npos || mPOLE==0 || erg_tev==0){
    cerr << "Arguments should follow the order \"[mH] [sqrts] [main directory] [output file name] [file1] [optional: file2] ...\"." << endl;
    return 0;
  }
  coutput = cdir + "/" + coutput;
  char TREE_NAME[] = "SelectedTree";

  vector<string> filename;
  for (int a=5; a<argc; a++){
    string ftmp = argv[a];
    if (ftmp.find(".lhe")!=string::npos){
      ftmp = cdir + "/" + ftmp;
      filename.push_back(ftmp);
    }
    else cerr << "Argument " << a << " is not a LHE file! Skipping this file." << endl;
  }
  if (filename.size()==0){
    cerr << "No valid LHE files are passed." << endl;
    return 0;
  }

  cout << "Creating file " << coutput << endl;

  TFile* foutput = new TFile(coutput.c_str(), "recreate");
  foutput->cd();
  TTree* tree = new TTree(TREE_NAME, TREE_NAME);
  tree->SetAutoSave(5000000000);

  Float_t GenhelcosthetaZ1, GenhelcosthetaZ2, Genhelphi, Gencosthetastar, GenphistarZ1;
  Float_t GenHMass, GenZ1Mass, GenZ2Mass, GenZaMass, GenZbMass;
  Float_t GenHPt, GenZ1Pt, GenZ2Pt, GenZaPt, GenZbPt;
  Float_t GenHPhi, GenZ1Phi, GenZ2Phi, GenZaPhi, GenZbPhi;
  Float_t GenHPz, GenZ1Eta, GenZ2Eta, GenZaEta, GenZbEta;

  Float_t GenLep1Mass, GenLep2Mass, GenLep3Mass, GenLep4Mass;
  Float_t GenLep1Pt, GenLep2Pt, GenLep3Pt, GenLep4Pt;
  Float_t GenLep1Eta, GenLep2Eta, GenLep3Eta, GenLep4Eta;
  Float_t GenLep1Phi, GenLep2Phi, GenLep3Phi, GenLep4Phi;
  Int_t GenLep1Id, GenLep2Id, GenLep3Id, GenLep4Id;

  vector<double> GenMotherMass;
  vector<double> GenMotherPt;
  vector<double> GenMotherPz;
  vector<double> GenMotherPhi;
  vector<int> GenMotherId;

  int NGenAssociatedVs=0;
  vector<double> GenAssociatedParticleMass;
  vector<double> GenAssociatedParticlePt;
  vector<double> GenAssociatedParticleEta;
  vector<double> GenAssociatedParticlePhi;
  vector<int> GenAssociatedParticleId;
  vector<double> GenAssociatedVMass;
  vector<double> GenAssociatedVPt;
  vector<double> GenAssociatedVEta;
  vector<double> GenAssociatedVPhi;
  vector<int> GenAssociatedVId;
  vector<int> GenAssociatedV_Particle1Index;
  vector<int> GenAssociatedV_Particle2Index;

  Float_t helcosthetaZ1, helcosthetaZ2, helphi, costhetastar, phistarZ1;
  Float_t ZZMass, Z1Mass, Z2Mass, ZaMass, ZbMass;
  Float_t ZZPt, ZZPz, Z1Pt, Z2Pt, ZaPt, ZbPt;
  Float_t ZZPhi, Z1Phi, Z2Phi, ZaPhi, ZbPhi;
  Float_t ZZEta, Z1Eta, Z2Eta, ZaEta, ZbEta;

  Float_t Lep1Mass, Lep2Mass, Lep3Mass, Lep4Mass;
  Float_t Lep1Pt, Lep2Pt, Lep3Pt, Lep4Pt;
  Float_t Lep1Eta, Lep2Eta, Lep3Eta, Lep4Eta;
  Float_t Lep1Phi, Lep2Phi, Lep3Phi, Lep4Phi;
  Int_t Lep1Id, Lep2Id, Lep3Id, Lep4Id;

  int NAssociatedVs=0;
  vector<double> AssociatedParticleMass;
  vector<double> AssociatedParticlePt;
  vector<double> AssociatedParticleEta;
  vector<double> AssociatedParticlePhi;
  vector<int> AssociatedParticleId;
  vector<double> AssociatedVMass;
  vector<double> AssociatedVPt;
  vector<double> AssociatedVEta;
  vector<double> AssociatedVPhi;
  vector<int> AssociatedVId;
  vector<int> AssociatedV_Particle1Index;
  vector<int> AssociatedV_Particle2Index;


  double weight;
  Float_t MC_weight;
  Int_t genFinalState=-1;
  Int_t isSelected=0;

  tree->Branch("GenHMass", &GenHMass);
  tree->Branch("GenHPt", &GenHPt);
  tree->Branch("GenHPz", &GenHPz);
  tree->Branch("GenHPhi", &GenHPhi);

  tree->Branch("GenZ1Mass", &GenZ1Mass);
  tree->Branch("GenZ1Pt", &GenZ1Pt);
  tree->Branch("GenZ1Phi", &GenZ1Phi);
  tree->Branch("GenZ1Eta", &GenZ1Eta);

  tree->Branch("GenZ2Mass", &GenZ2Mass);
  tree->Branch("GenZ2Pt", &GenZ2Pt);
  tree->Branch("GenZ2Phi", &GenZ2Phi);
  tree->Branch("GenZ2Eta", &GenZ2Eta);

  tree->Branch("GenZaMass", &GenZaMass);
  tree->Branch("GenZaPt", &GenZaPt);
  tree->Branch("GenZaPhi", &GenZaPhi);
  tree->Branch("GenZaEta", &GenZaEta);

  tree->Branch("GenZbMass", &GenZbMass);
  tree->Branch("GenZbPt", &GenZbPt);
  tree->Branch("GenZbPhi", &GenZbPhi);
  tree->Branch("GenZbEta", &GenZbEta);

  tree->Branch("GenMotherMass", &GenMotherMass);
  tree->Branch("GenMotherPt", &GenMotherPt);
  tree->Branch("GenMotherPz", &GenMotherPz);
  tree->Branch("GenMotherPhi", &GenMotherPhi);
  tree->Branch("GenMotherId", &GenMotherId);

  tree->Branch("NGenAssociatedVs", &NGenAssociatedVs);
  tree->Branch("GenAssociatedParticleMass", &GenAssociatedParticleMass);
  tree->Branch("GenAssociatedParticlePt", &GenAssociatedParticlePt);
  tree->Branch("GenAssociatedParticleEta", &GenAssociatedParticleEta);
  tree->Branch("GenAssociatedParticlePhi", &GenAssociatedParticlePhi);
  tree->Branch("GenAssociatedParticleId", &GenAssociatedParticleId);
  tree->Branch("GenAssociatedVMass", &GenAssociatedVMass);
  tree->Branch("GenAssociatedVPt", &GenAssociatedVPt);
  tree->Branch("GenAssociatedVEta", &GenAssociatedVEta);
  tree->Branch("GenAssociatedVPhi", &GenAssociatedVPhi);
  tree->Branch("GenAssociatedVId", &GenAssociatedVId);
  tree->Branch("GenAssociatedV_Particle1Index", &GenAssociatedV_Particle1Index);
  tree->Branch("GenAssociatedV_Particle2Index", &GenAssociatedV_Particle2Index);

  tree->Branch("GenhelcosthetaZ1", &GenhelcosthetaZ1);
  tree->Branch("GenhelcosthetaZ2", &GenhelcosthetaZ2);
  tree->Branch("Genhelphi", &Genhelphi);
  tree->Branch("Gencosthetastar", &Gencosthetastar);
  tree->Branch("GenphistarZ1", &GenphistarZ1);

  tree->Branch("GenLep1Mass", &GenLep1Mass);
  tree->Branch("GenLep2Mass", &GenLep2Mass);
  tree->Branch("GenLep3Mass", &GenLep3Mass);
  tree->Branch("GenLep4Mass", &GenLep4Mass);
  tree->Branch("GenLep1Pt", &GenLep1Pt);
  tree->Branch("GenLep2Pt", &GenLep2Pt);
  tree->Branch("GenLep3Pt", &GenLep3Pt);
  tree->Branch("GenLep4Pt", &GenLep4Pt);
  tree->Branch("GenLep1Eta", &GenLep1Eta);
  tree->Branch("GenLep2Eta", &GenLep2Eta);
  tree->Branch("GenLep3Eta", &GenLep3Eta);
  tree->Branch("GenLep4Eta", &GenLep4Eta);
  tree->Branch("GenLep1Phi", &GenLep1Phi);
  tree->Branch("GenLep2Phi", &GenLep2Phi);
  tree->Branch("GenLep3Phi", &GenLep3Phi);
  tree->Branch("GenLep4Phi", &GenLep4Phi);
  tree->Branch("GenLep1Id", &GenLep1Id);
  tree->Branch("GenLep2Id", &GenLep2Id);
  tree->Branch("GenLep3Id", &GenLep3Id);
  tree->Branch("GenLep4Id", &GenLep4Id);

  tree->Branch("genFinalState", &genFinalState);
  tree->Branch("MC_weight", &MC_weight);


  tree->Branch("isSelected", &isSelected);

  tree->Branch("ZZMass", &ZZMass);
  tree->Branch("ZZPt", &ZZPt);
  tree->Branch("ZZPz", &ZZPz);
  tree->Branch("ZZPhi", &ZZPhi);

  tree->Branch("Z1Mass", &Z1Mass);
  tree->Branch("Z1Pt", &Z1Pt);
  tree->Branch("Z1Phi", &Z1Phi);
  tree->Branch("Z1Eta", &Z1Eta);

  tree->Branch("Z2Mass", &Z2Mass);
  tree->Branch("Z2Pt", &Z2Pt);
  tree->Branch("Z2Phi", &Z2Phi);
  tree->Branch("Z2Eta", &Z2Eta);

  tree->Branch("ZaMass", &ZaMass);
  tree->Branch("ZaPt", &ZaPt);
  tree->Branch("ZaPhi", &ZaPhi);
  tree->Branch("ZaEta", &ZaEta);

  tree->Branch("ZbMass", &ZbMass);
  tree->Branch("ZbPt", &ZbPt);
  tree->Branch("ZbPhi", &ZbPhi);
  tree->Branch("ZbEta", &ZbEta);

  tree->Branch("NAssociatedVs", &NAssociatedVs);
  tree->Branch("AssociatedParticleMass", &AssociatedParticleMass);
  tree->Branch("AssociatedParticlePt", &AssociatedParticlePt);
  tree->Branch("AssociatedParticleEta", &AssociatedParticleEta);
  tree->Branch("AssociatedParticlePhi", &AssociatedParticlePhi);
  tree->Branch("AssociatedParticleId", &AssociatedParticleId);
  tree->Branch("AssociatedVMass", &AssociatedVMass);
  tree->Branch("AssociatedVPt", &AssociatedVPt);
  tree->Branch("AssociatedVEta", &AssociatedVEta);
  tree->Branch("AssociatedVPhi", &AssociatedVPhi);
  tree->Branch("AssociatedVId", &AssociatedVId);
  tree->Branch("AssociatedV_Particle1Index", &AssociatedV_Particle1Index);
  tree->Branch("AssociatedV_Particle2Index", &AssociatedV_Particle2Index);

  tree->Branch("helcosthetaZ1", &helcosthetaZ1);
  tree->Branch("helcosthetaZ2", &helcosthetaZ2);
  tree->Branch("helphi", &helphi);
  tree->Branch("costhetastar", &costhetastar);
  tree->Branch("phistarZ1", &phistarZ1);

  tree->Branch("Lep1Mass", &Lep1Mass);
  tree->Branch("Lep2Mass", &Lep2Mass);
  tree->Branch("Lep3Mass", &Lep3Mass);
  tree->Branch("Lep4Mass", &Lep4Mass);
  tree->Branch("Lep1Pt", &Lep1Pt);
  tree->Branch("Lep2Pt", &Lep2Pt);
  tree->Branch("Lep3Pt", &Lep3Pt);
  tree->Branch("Lep4Pt", &Lep4Pt);
  tree->Branch("Lep1Eta", &Lep1Eta);
  tree->Branch("Lep2Eta", &Lep2Eta);
  tree->Branch("Lep3Eta", &Lep3Eta);
  tree->Branch("Lep4Eta", &Lep4Eta);
  tree->Branch("Lep1Phi", &Lep1Phi);
  tree->Branch("Lep2Phi", &Lep2Phi);
  tree->Branch("Lep3Phi", &Lep3Phi);
  tree->Branch("Lep4Phi", &Lep4Phi);
  tree->Branch("Lep1Id", &Lep1Id);
  tree->Branch("Lep2Id", &Lep2Id);
  tree->Branch("Lep3Id", &Lep3Id);
  tree->Branch("Lep4Id", &Lep4Id);


  for (int f=0; f<filename.size(); f++){
    string cinput = filename.at(f);
    ifstream fin;
    fin.open(cinput.c_str());
    if (fin.good()){
      while (!fin.eof()){
        GenMotherMass.clear();
        GenMotherPt.clear();
        GenMotherPz.clear();
        GenMotherPhi.clear();
        GenMotherId.clear();

        GenAssociatedParticleMass.clear();
        GenAssociatedParticlePt.clear();
        GenAssociatedParticleEta.clear();
        GenAssociatedParticlePhi.clear();
        GenAssociatedParticleId.clear();
        GenAssociatedVMass.clear();
        GenAssociatedVPt.clear();
        GenAssociatedVEta.clear();
        GenAssociatedVPhi.clear();
        GenAssociatedVId.clear();
        GenAssociatedV_Particle1Index.clear();
        GenAssociatedV_Particle2Index.clear();

        AssociatedParticleMass.clear();
        AssociatedParticlePt.clear();
        AssociatedParticleEta.clear();
        AssociatedParticlePhi.clear();
        AssociatedParticleId.clear();
        AssociatedVMass.clear();
        AssociatedVPt.clear();
        AssociatedVEta.clear();
        AssociatedVPhi.clear();
        AssociatedVId.clear();
        AssociatedV_Particle1Index.clear();
        AssociatedV_Particle2Index.clear();

        vector<Particle*> particleList = readLHEEvent(fin, weight);
        vector<Particle*> smearedParticleList; // Bookkeeping
        vector<ZZCandidate*> candList; // Bookkeeping
        vector<ZZCandidate*> smearedCandList; // Bookkeeping

        if (particleList.size()==0 && weight!=0) weight=0;
        if (weight!=0){
          Event genEvent;
          genEvent.setWeight(weight);
          bool hasGenHiggs=false;
          Event smearedEvent;
          smearedEvent.setWeight(weight);
          for (int p=0; p<particleList.size(); p++){
            Particle* genPart = particleList.at(p); // Has mother info from LHE reading
            if (isAHiggs(genPart->id)) hasGenHiggs=true;

            if (genPart->genStatus==1){
              if (isALepton(genPart->id)) genEvent.addLepton(genPart);
              else if (isANeutrino(genPart->id)) genEvent.addNeutrino(genPart);
              else if (isAGluon(genPart->id) || isAQuark(genPart->id)) genEvent.addJet(genPart);

              Particle* smearedPart = smearParticle(genPart); // Has no mother info
              smearedParticleList.push_back(smearedPart);
              if (isALepton(smearedPart->id)) smearedEvent.addLepton(smearedPart);
              else if (isANeutrino(smearedPart->id)) smearedEvent.addNeutrino(smearedPart);
              else if (isAGluon(smearedPart->id) || isAQuark(smearedPart->id)) smearedEvent.addJet(smearedPart);
              else smearedEvent.addParticle(smearedPart);
            }
            else if (genPart->genStatus==-1){
              GenMotherMass.push_back(genPart->m());
              GenMotherPt.push_back(genPart->pt());
              GenMotherPz.push_back(genPart->z());
              GenMotherPhi.push_back(genPart->phi());
              GenMotherId.push_back(genPart->id);
            }
          }

          //for (int p=0; p<genEvent.getNLeptons(); p++) cout << "Lepton " << p << " (x, y, z, t): " << genEvent.getLepton(p)->x() << '\t' << genEvent.getLepton(p)->y() << '\t' << genEvent.getLepton(p)->z() << '\t' << genEvent.getLepton(p)->t() << endl;

          genEvent.constructVVCandidates();
          genEvent.addVVCandidateAppendages();
          ZZCandidate* genCand=0;
          for (int t=0; t<genEvent.getNZZCandidates(); t++){
            ZZCandidate* tmpCand = genEvent.getZZCandidate(t);
            if (hasGenHiggs){
              if (!(isAHiggs(tmpCand->getSortedV(0)->getDaughter(0)->getMother(0)->id) || isAHiggs(tmpCand->getSortedV(0)->getDaughter(0)->getMother(0)->getMother(0)->id))) continue;
            }
            if (genCand==0) genCand=tmpCand;
            else if (fabs(genCand->getSortedV(0)->m()-PDGHelpers::HVVmass)>fabs(tmpCand->getSortedV(0)->m()-PDGHelpers::HVVmass)) genCand=tmpCand;
          }
          if (genCand!=0){
            Particle* gZ1=genCand->getSortedV(0);
            Particle* gZ2=genCand->getSortedV(1);

            GenHMass=genCand->m();
            GenHPt=genCand->pt();
            GenHPz=genCand->z();
            GenHPhi=genCand->phi();

            GenZ1Mass=gZ1->m();
            GenZ1Pt=gZ1->pt();
            GenZ1Eta=gZ1->eta();
            GenZ1Phi=gZ1->phi();

            GenZ2Mass=gZ2->m();
            GenZ2Pt=gZ2->pt();
            GenZ2Eta=gZ2->eta();
            GenZ2Phi=gZ2->phi();

            TLorentzVector pZ1alt = genCand->getAlternativeVMomentum(0);
            TLorentzVector pZ2alt = genCand->getAlternativeVMomentum(1);

            GenZaMass=pZ1alt.M();
            GenZaPt=pZ1alt.Pt();
            GenZaEta=pZ1alt.Eta();
            GenZaPhi=pZ1alt.Phi();

            GenZbMass=pZ2alt.M();
            GenZbPt=pZ2alt.Pt();
            GenZbEta=pZ2alt.Eta();
            GenZbPhi=pZ2alt.Phi();

            calculateAngles(
              genCand->p4,
              gZ1->getDaughter(0)->p4, gZ1->getDaughter(1)->p4,
              gZ2->getDaughter(0)->p4, gZ2->getDaughter(1)->p4,
              GenhelcosthetaZ1, GenhelcosthetaZ2, Genhelphi, Gencosthetastar, GenphistarZ1
              );

            GenLep1Id = gZ1->getDaughter(0)->id;
            GenLep1Mass = gZ1->getDaughter(0)->m();
            GenLep1Pt = gZ1->getDaughter(0)->pt();
            GenLep1Eta = gZ1->getDaughter(0)->eta();
            GenLep1Phi = gZ1->getDaughter(0)->phi();

            GenLep2Id = gZ1->getDaughter(1)->id;
            GenLep2Mass = gZ1->getDaughter(1)->m();
            GenLep2Pt = gZ1->getDaughter(1)->pt();
            GenLep2Eta = gZ1->getDaughter(1)->eta();
            GenLep2Phi = gZ1->getDaughter(1)->phi();

            GenLep3Id = gZ2->getDaughter(0)->id;
            GenLep3Mass = gZ2->getDaughter(0)->m();
            GenLep3Pt = gZ2->getDaughter(0)->pt();
            GenLep3Eta = gZ2->getDaughter(0)->eta();
            GenLep3Phi = gZ2->getDaughter(0)->phi();

            GenLep4Id = gZ2->getDaughter(1)->id;
            GenLep4Mass = gZ2->getDaughter(1)->m();
            GenLep4Pt = gZ2->getDaughter(1)->pt();
            GenLep4Eta = gZ2->getDaughter(1)->eta();
            GenLep4Phi = gZ2->getDaughter(1)->phi();

            for (int aa=0; aa<genCand->getNAssociatedJets(); aa++){
              Particle* apart = genCand->getAssociatedJet(aa);
              GenAssociatedParticleMass.push_back(apart->m());
              GenAssociatedParticlePt.push_back(apart->pt());
              GenAssociatedParticleEta.push_back(apart->eta());
              GenAssociatedParticlePhi.push_back(apart->phi());
              GenAssociatedParticleId.push_back(apart->id);
            }
            for (int aa=0; aa<genCand->getNAssociatedLeptons(); aa++){
              Particle* apart = genCand->getAssociatedLepton(aa);
              GenAssociatedParticleMass.push_back(apart->m());
              GenAssociatedParticlePt.push_back(apart->pt());
              GenAssociatedParticleEta.push_back(apart->eta());
              GenAssociatedParticlePhi.push_back(apart->phi());
              GenAssociatedParticleId.push_back(apart->id);
            }
            for (int aa=0; aa<genCand->getNAssociatedNeutrinos(); aa++){
              Particle* apart = genCand->getAssociatedNeutrino(aa);
              GenAssociatedParticleMass.push_back(apart->m());
              GenAssociatedParticlePt.push_back(apart->pt());
              GenAssociatedParticleEta.push_back(apart->eta());
              GenAssociatedParticlePhi.push_back(apart->phi());
              GenAssociatedParticleId.push_back(apart->id);
            }
            NGenAssociatedVs = genCand->getNSortedVs()-2;
            for (int av=2; av<genCand->getNSortedVs(); av++){
              Particle* associatedV = genCand->getSortedV(av);
              GenAssociatedVMass.push_back(associatedV->m());
              GenAssociatedVPt.push_back(associatedV->pt());
              GenAssociatedVEta.push_back(associatedV->eta());
              GenAssociatedVPhi.push_back(associatedV->phi());
              GenAssociatedVId.push_back(associatedV->id);

              Particle* avd1 = associatedV->getDaughter(0);
              Particle* avd2 = associatedV->getDaughter(1);
              for (int aa=0; aa<genCand->getNAssociatedJets(); aa++){
                if (avd1==genCand->getAssociatedJet(aa)) GenAssociatedV_Particle1Index.push_back(aa);
                else if (avd2==genCand->getAssociatedJet(aa)) GenAssociatedV_Particle2Index.push_back(aa);
              }
              for (int aa=0; aa<genCand->getNAssociatedLeptons(); aa++){
                if (avd1==genCand->getAssociatedLepton(aa)) GenAssociatedV_Particle1Index.push_back(aa+genCand->getNAssociatedJets());
                else if (avd2==genCand->getAssociatedLepton(aa)) GenAssociatedV_Particle2Index.push_back(aa+genCand->getNAssociatedJets());
              }
              for (int aa=0; aa<genCand->getNAssociatedNeutrinos(); aa++){
                if (avd1==genCand->getAssociatedNeutrino(aa)) GenAssociatedV_Particle1Index.push_back(aa+genCand->getNAssociatedJets()+genCand->getNAssociatedLeptons());
                else if (avd2==genCand->getAssociatedNeutrino(aa)) GenAssociatedV_Particle2Index.push_back(aa+genCand->getNAssociatedJets()+genCand->getNAssociatedLeptons());
              }
            }
          }
          else{
            cout << "No gen. level Higgs candidate was found!" << endl;

            GenHMass=0;
            GenHPt=0;
            GenHPz=0;
            GenHPhi=0;

            GenZ1Mass=0;
            GenZ1Pt=0;
            GenZ1Eta=0;
            GenZ1Phi=0;

            GenZ2Mass=0;
            GenZ2Pt=0;
            GenZ2Eta=0;
            GenZ2Phi=0;

            GenZaMass=0;
            GenZaPt=0;
            GenZaEta=0;
            GenZaPhi=0;

            GenZbMass=0;
            GenZbPt=0;
            GenZbEta=0;
            GenZbPhi=0;

            GenhelcosthetaZ1=0;
            GenhelcosthetaZ2=0;
            Genhelphi=0;
            Gencosthetastar=0;
            GenphistarZ1=0;

            GenLep1Id = 0;
            GenLep1Mass = 0;
            GenLep1Pt = 0;
            GenLep1Eta = 0;
            GenLep1Phi = 0;

            GenLep2Id = 0;
            GenLep2Mass = 0;
            GenLep2Pt = 0;
            GenLep2Eta = 0;
            GenLep2Phi = 0;

            GenLep3Id = 0;
            GenLep3Mass = 0;
            GenLep3Pt = 0;
            GenLep3Eta = 0;
            GenLep3Phi = 0;

            GenLep4Id = 0;
            GenLep4Mass = 0;
            GenLep4Pt = 0;
            GenLep4Eta = 0;
            GenLep4Phi = 0;

            NGenAssociatedVs = 0;
          }

          smearedEvent.constructVVCandidates();
          smearedEvent.applyParticleSelection();
          smearedEvent.addVVCandidateAppendages();
          ZZCandidate* rCand=0;
          for (int t=0; t<smearedEvent.getNZZCandidates(); t++){
            ZZCandidate* tmpCand = smearedEvent.getZZCandidate(t);
            if (!tmpCand->passSelection) continue;
            if (rCand==0) rCand=tmpCand;
            else if (fabs(rCand->getSortedV(0)->m()-PDGHelpers::HVVmass)>fabs(tmpCand->getSortedV(0)->m()-PDGHelpers::HVVmass)) rCand=tmpCand;
          }
          if (rCand!=0){
            isSelected=1;

            Particle* rZ1=rCand->getSortedV(0);
            Particle* rZ2=rCand->getSortedV(1);

            ZZMass=rCand->m();
            ZZPt=rCand->pt();
            ZZPz=rCand->z();
            ZZEta=rCand->eta();

            Z1Mass=rZ1->m();
            Z1Pt=rZ1->pt();
            Z1Eta=rZ1->eta();
            Z1Phi=rZ1->phi();

            Z2Mass=rZ2->m();
            Z2Pt=rZ2->pt();
            Z2Eta=rZ2->eta();
            Z2Phi=rZ2->phi();

            TLorentzVector pZ1alt = rCand->getAlternativeVMomentum(0);
            TLorentzVector pZ2alt = rCand->getAlternativeVMomentum(1);

            ZaMass=pZ1alt.M();
            ZaPt=pZ1alt.Pt();
            ZaEta=pZ1alt.Eta();
            ZaPhi=pZ1alt.Phi();

            ZbMass=pZ2alt.M();
            ZbPt=pZ2alt.Pt();
            ZbEta=pZ2alt.Eta();
            ZbPhi=pZ2alt.Phi();

            calculateAngles(
              rCand->p4,
              rZ1->getDaughter(0)->p4, rZ1->getDaughter(1)->p4,
              rZ2->getDaughter(0)->p4, rZ2->getDaughter(1)->p4,
              helcosthetaZ1, helcosthetaZ2, helphi, costhetastar, phistarZ1
              );

            Lep1Id = rZ1->getDaughter(0)->id;
            Lep1Mass = rZ1->getDaughter(0)->m();
            Lep1Pt = rZ1->getDaughter(0)->pt();
            Lep1Eta = rZ1->getDaughter(0)->eta();
            Lep1Phi = rZ1->getDaughter(0)->phi();

            Lep2Id = rZ1->getDaughter(1)->id;
            Lep2Mass = rZ1->getDaughter(1)->m();
            Lep2Pt = rZ1->getDaughter(1)->pt();
            Lep2Eta = rZ1->getDaughter(1)->eta();
            Lep2Phi = rZ1->getDaughter(1)->phi();

            Lep3Id = rZ2->getDaughter(0)->id;
            Lep3Mass = rZ2->getDaughter(0)->m();
            Lep3Pt = rZ2->getDaughter(0)->pt();
            Lep3Eta = rZ2->getDaughter(0)->eta();
            Lep3Phi = rZ2->getDaughter(0)->phi();

            Lep4Id = rZ2->getDaughter(1)->id;
            Lep4Mass = rZ2->getDaughter(1)->m();
            Lep4Pt = rZ2->getDaughter(1)->pt();
            Lep4Eta = rZ2->getDaughter(1)->eta();
            Lep4Phi = rZ2->getDaughter(1)->phi();

            for (int aa=0; aa<rCand->getNAssociatedJets(); aa++){
              Particle* apart = rCand->getAssociatedJet(aa);
              AssociatedParticleMass.push_back(apart->m());
              AssociatedParticlePt.push_back(apart->pt());
              AssociatedParticleEta.push_back(apart->eta());
              AssociatedParticlePhi.push_back(apart->phi());
              AssociatedParticleId.push_back(apart->id);
            }
            for (int aa=0; aa<rCand->getNAssociatedLeptons(); aa++){
              Particle* apart = rCand->getAssociatedLepton(aa);
              AssociatedParticleMass.push_back(apart->m());
              AssociatedParticlePt.push_back(apart->pt());
              AssociatedParticleEta.push_back(apart->eta());
              AssociatedParticlePhi.push_back(apart->phi());
              AssociatedParticleId.push_back(apart->id);
            }
            for (int aa=0; aa<rCand->getNAssociatedNeutrinos(); aa++){
              Particle* apart = rCand->getAssociatedNeutrino(aa);
              AssociatedParticleMass.push_back(apart->m());
              AssociatedParticlePt.push_back(apart->pt());
              AssociatedParticleEta.push_back(apart->eta());
              AssociatedParticlePhi.push_back(apart->phi());
              AssociatedParticleId.push_back(apart->id);
            }
            NAssociatedVs = rCand->getNSortedVs()-2;
            for (int av=2; av<rCand->getNSortedVs(); av++){
              Particle* associatedV = rCand->getSortedV(av);
              AssociatedVMass.push_back(associatedV->m());
              AssociatedVPt.push_back(associatedV->pt());
              AssociatedVEta.push_back(associatedV->eta());
              AssociatedVPhi.push_back(associatedV->phi());
              AssociatedVId.push_back(associatedV->id);

              Particle* avd1 = associatedV->getDaughter(0);
              Particle* avd2 = associatedV->getDaughter(1);
              for (int aa=0; aa<rCand->getNAssociatedJets(); aa++){
                if (avd1==rCand->getAssociatedJet(aa)) AssociatedV_Particle1Index.push_back(aa);
                else if (avd2==rCand->getAssociatedJet(aa)) AssociatedV_Particle2Index.push_back(aa);
              }
              for (int aa=0; aa<rCand->getNAssociatedLeptons(); aa++){
                if (avd1==rCand->getAssociatedLepton(aa)) AssociatedV_Particle1Index.push_back(aa+rCand->getNAssociatedJets());
                else if (avd2==rCand->getAssociatedLepton(aa)) AssociatedV_Particle2Index.push_back(aa+rCand->getNAssociatedJets());
              }
              for (int aa=0; aa<rCand->getNAssociatedNeutrinos(); aa++){
                if (avd1==rCand->getAssociatedNeutrino(aa)) AssociatedV_Particle1Index.push_back(aa+rCand->getNAssociatedJets()+rCand->getNAssociatedLeptons());
                else if (avd2==rCand->getAssociatedNeutrino(aa)) AssociatedV_Particle2Index.push_back(aa+rCand->getNAssociatedJets()+rCand->getNAssociatedLeptons());
              }
            }
          }
          else{
            isSelected=0;

            ZZMass=0;
            ZZPt=0;
            ZZPz=0;
            ZZPhi=0;

            Z1Mass=0;
            Z1Pt=0;
            Z1Eta=0;
            Z1Phi=0;

            Z2Mass=0;
            Z2Pt=0;
            Z2Eta=0;
            Z2Phi=0;

            ZaMass=0;
            ZaPt=0;
            ZaEta=0;
            ZaPhi=0;

            ZbMass=0;
            ZbPt=0;
            ZbEta=0;
            ZbPhi=0;

            helcosthetaZ1=0;
            helcosthetaZ2=0;
            helphi=0;
            costhetastar=0;
            phistarZ1=0;

            Lep1Id = 0;
            Lep1Mass = 0;
            Lep1Pt = 0;
            Lep1Eta = 0;
            Lep1Phi = 0;

            Lep2Id = 0;
            Lep2Mass = 0;
            Lep2Pt = 0;
            Lep2Eta = 0;
            Lep2Phi = 0;

            Lep3Id = 0;
            Lep3Mass = 0;
            Lep3Pt = 0;
            Lep3Eta = 0;
            Lep3Phi = 0;

            Lep4Id = 0;
            Lep4Mass = 0;
            Lep4Pt = 0;
            Lep4Eta = 0;
            Lep4Phi = 0;

            NAssociatedVs = 0;
          }

          MC_weight = (float)weight;
          tree->Fill();
        }

        for (int p=0; p<smearedCandList.size(); p++){ // Bookkeeping
          ZZCandidate* tmpCand = (ZZCandidate*)smearedCandList.at(p);
          if (tmpCand!=0) delete tmpCand;
        }
        for (int p=0; p<smearedParticleList.size(); p++){ // Bookkeeping
          Particle* tmpPart = (Particle*)smearedParticleList.at(p);
          if (tmpPart!=0) delete tmpPart;
        }

        for (int p=0; p<candList.size(); p++){ // Bookkeeping
          ZZCandidate* tmpCand = (ZZCandidate*)candList.at(p);
          if (tmpCand!=0) delete tmpCand;
        }
        for (int p=0; p<particleList.size(); p++){ // Bookkeeping
          Particle* tmpPart = (Particle*)particleList.at(p);
          if (tmpPart!=0) delete tmpPart;
        }

        // Bookkeeping
        smearedCandList.clear();
        smearedParticleList.clear();
        candList.clear();
        particleList.clear();
      }
      fin.close();
    }
  }

  cout << "Number of recorded events: " << tree->GetEntries() << endl;
  foutput->WriteTObject(tree);
  delete tree;
  foutput->Close();
  return 0;
}

vector<Particle*> readLHEEvent(ifstream& input_lhe, double& weight){
  string event_beginning = "<event>";
  string event_end = "</event>";
  string file_closing = "</LesHouchesEvents>";
  string str_in="";

  vector<Particle*> collection;
  vector<int> motherIDs_first;
  vector<int> motherIDs_second;

// Test whether the string read is the beginning of the event for a valid file
  while (str_in.find(event_beginning)==string::npos){
    if (input_lhe.eof()){
      weight=0;
      return collection;
    }
    getline(input_lhe, str_in);
    if (str_in.find(file_closing)!=string::npos){
      weight=0;
      return collection;
    }
  }

  int nparticle, para;
  double m_V, alpha_qed, alpha_s;

  input_lhe >> nparticle >> para >> weight >> m_V >> alpha_qed >> alpha_s;
  for (int a = 0; a < nparticle; a++){
    int idup, istup, mothup[2], icolup[2];
    double pup[5], vtimup, spinup;
    TLorentzVector partFourVec;

    input_lhe >> idup >> istup >> mothup[0] >> mothup[1] >> icolup[0] >> icolup[1];
    for (int i = 0; i < 5; i++){
      input_lhe >> pup[i];
    }
    input_lhe >> vtimup >> spinup;

    motherIDs_first.push_back(mothup[0]);
    motherIDs_second.push_back(mothup[1]);

    partFourVec.SetXYZT(pup[0], pup[1], pup[2], pup[3]);
    Particle* onePart = new Particle(idup, partFourVec);
    onePart->setGenStatus(istup);
    onePart->setLifetime(spinup);
    collection.push_back(onePart);
  }

// Test whether the end of event is reached indeed
  for(int t=0;t<2;t++) getline(input_lhe, str_in); // Read twice to get rid of the end-of-line
  while (str_in.find("#")!=string::npos) getline(input_lhe, str_in);
  if (str_in.find(event_end)==string::npos){
    cerr << "End of event not reached! string is " << str_in << endl;
    weight=0;
    for (int a = 0; a < collection.size(); a++){
      Particle* tmpPart = collection.at(a);
      delete tmpPart;
    }
    collection.clear();
    motherIDs_first.clear();
    motherIDs_second.clear();
    weight=0;
    return collection;
  }

// Assign the mothers
  for (int a = 0; a < nparticle; a++){
    if (motherIDs_first.at(a)>0) collection.at(a)->addMother(collection.at(motherIDs_first.at(a)-1));
    if (motherIDs_second.at(a)>0 && motherIDs_first.at(a) != motherIDs_second.at(a)) collection.at(a)->addMother(collection.at(motherIDs_second.at(a)-1));
  }

  return collection;
}


void calculateAngles(TLorentzVector thep4H, TLorentzVector thep4M11, TLorentzVector thep4M12, TLorentzVector thep4M21, TLorentzVector thep4M22, float& costheta1, float& costheta2, float& phi, float& costhetastar, float& phistar1){
	
  TLorentzVector thep4Z1 = thep4M11+thep4M12;
  TLorentzVector thep4Z2 = thep4M21+thep4M22;

  float norm;
  
  TVector3 boostX = -(thep4H.BoostVector());
  TLorentzVector thep4Z1inXFrame( thep4Z1 );
  TLorentzVector thep4Z2inXFrame( thep4Z2 );	
  thep4Z1inXFrame.Boost( boostX );
  thep4Z2inXFrame.Boost( boostX );
  TVector3 theZ1X_p3 = TVector3( thep4Z1inXFrame.X(), thep4Z1inXFrame.Y(), thep4Z1inXFrame.Z() );
  TVector3 theZ2X_p3 = TVector3( thep4Z2inXFrame.X(), thep4Z2inXFrame.Y(), thep4Z2inXFrame.Z() );

  
  ///////////////////////////////////////////////
  // check for z1/z2 convention, redefine all 4 vectors with convention
  ///////////////////////////////////////////////	
  TLorentzVector p4H, p4Z1, p4M11, p4M12, p4Z2, p4M21, p4M22;
  p4H = thep4H;
   
  p4Z1 = thep4Z1; p4M11 = thep4M11; p4M12 = thep4M12;
  p4Z2 = thep4Z2; p4M21 = thep4M21; p4M22 = thep4M22;
  costhetastar = theZ1X_p3.CosTheta();
	
  // now helicity angles................................
  // ...................................................
  TVector3 boostZ1 = -(p4Z1.BoostVector());
  TLorentzVector p4Z2Z1(p4Z2);
  p4Z2Z1.Boost(boostZ1);
  // find the decay axis
  TVector3 unitx_1( -p4Z2Z1.X(), -p4Z2Z1.Y(), -p4Z2Z1.Z() );
  norm = 1/(unitx_1.Mag());
  unitx_1*=norm;
  // boost daughters of z2
  TLorentzVector p4M21Z1(p4M21);
  TLorentzVector p4M22Z1(p4M22);
  p4M21Z1.Boost(boostZ1);
  p4M22Z1.Boost(boostZ1);
  // create z and y axes
  TVector3 p4M21Z1_p3( p4M21Z1.X(), p4M21Z1.Y(), p4M21Z1.Z() );
  TVector3 p4M22Z1_p3( p4M22Z1.X(), p4M22Z1.Y(), p4M22Z1.Z() );
  TVector3 unitz_1 = p4M21Z1_p3.Cross( p4M22Z1_p3 );
  norm = 1/(unitz_1.Mag());
  unitz_1 *= norm;
  TVector3 unity_1 = unitz_1.Cross(unitx_1);
  
  // calculate theta1
  TLorentzVector p4M11Z1(p4M11);
  p4M11Z1.Boost(boostZ1);
  TVector3 p3M11( p4M11Z1.X(), p4M11Z1.Y(), p4M11Z1.Z() );
  TVector3 unitM11 = p3M11.Unit();
  float x_m11 = unitM11.Dot(unitx_1); float y_m11 = unitM11.Dot(unity_1); float z_m11 = unitM11.Dot(unitz_1);
  TVector3 M11_Z1frame(y_m11, z_m11, x_m11);
  costheta1 = M11_Z1frame.CosTheta();

  //////-----------------------old way of calculating phi---------------/////////
  phi = M11_Z1frame.Phi();
  
  // set axes for other system
  TVector3 boostZ2 = -(p4Z2.BoostVector());
  TLorentzVector p4Z1Z2(p4Z1);
  p4Z1Z2.Boost(boostZ2);
  TVector3 unitx_2( -p4Z1Z2.X(), -p4Z1Z2.Y(), -p4Z1Z2.Z() );
  norm = 1/(unitx_2.Mag());
  unitx_2*=norm;
  // boost daughters of z2
  TLorentzVector p4M11Z2(p4M11);
  TLorentzVector p4M12Z2(p4M12);
  p4M11Z2.Boost(boostZ2);
  p4M12Z2.Boost(boostZ2);
  TVector3 p4M11Z2_p3( p4M11Z2.X(), p4M11Z2.Y(), p4M11Z2.Z() );
  TVector3 p4M12Z2_p3( p4M12Z2.X(), p4M12Z2.Y(), p4M12Z2.Z() );
  TVector3 unitz_2 = p4M11Z2_p3.Cross( p4M12Z2_p3 );
  norm = 1/(unitz_2.Mag());
  unitz_2*=norm;
  TVector3 unity_2 = unitz_2.Cross(unitx_2);
  // calcuate theta2
  TLorentzVector p4M21Z2(p4M21);
  p4M21Z2.Boost(boostZ2);
  TVector3 p3M21( p4M21Z2.X(), p4M21Z2.Y(), p4M21Z2.Z() );
  TVector3 unitM21 = p3M21.Unit();
  float x_m21 = unitM21.Dot(unitx_2); float y_m21 = unitM21.Dot(unity_2); float z_m21 = unitM21.Dot(unitz_2);
  TVector3 M21_Z2frame(y_m21, z_m21, x_m21);
  costheta2 = M21_Z2frame.CosTheta();
  
  // calculate phi
  // calculating phi_n
  TLorentzVector n_p4Z1inXFrame( p4Z1 );
  TLorentzVector n_p4M11inXFrame( p4M11 );
  n_p4Z1inXFrame.Boost( boostX );
  n_p4M11inXFrame.Boost( boostX );        
  TVector3 n_p4Z1inXFrame_unit = n_p4Z1inXFrame.Vect().Unit();
  TVector3 n_p4M11inXFrame_unit = n_p4M11inXFrame.Vect().Unit();  
  TVector3 n_unitz_1( n_p4Z1inXFrame_unit );
  //// y-axis is defined by neg lepton cross z-axis
  //// the subtle part is here...
  TVector3 n_unity_1 = n_unitz_1.Cross( n_p4M11inXFrame_unit );
  TVector3 n_unitx_1 = n_unity_1.Cross( n_unitz_1 );
  
  TLorentzVector n_p4M21inXFrame( p4M21 );
  n_p4M21inXFrame.Boost( boostX );
  TVector3 n_p4M21inXFrame_unit = n_p4M21inXFrame.Vect().Unit();
  //rotate into other plane
  TVector3 n_p4M21inXFrame_unitprime( n_p4M21inXFrame_unit.Dot(n_unitx_1), n_p4M21inXFrame_unit.Dot(n_unity_1), n_p4M21inXFrame_unit.Dot(n_unitz_1) );
  
  ///////-----------------new way of calculating phi-----------------///////
  // float phi_n =  n_p4M21inXFrame_unitprime.Phi();
  /// and then calculate phistar1
  TVector3 n_p4PartoninXFrame_unit( 0.0, 0.0, 1.0 );
  TVector3 n_p4PartoninXFrame_unitprime( n_p4PartoninXFrame_unit.Dot(n_unitx_1), n_p4PartoninXFrame_unit.Dot(n_unity_1), n_p4PartoninXFrame_unit.Dot(n_unitz_1) );
  // negative sign is for arrow convention in paper
  phistar1 = (n_p4PartoninXFrame_unitprime.Phi());
}





