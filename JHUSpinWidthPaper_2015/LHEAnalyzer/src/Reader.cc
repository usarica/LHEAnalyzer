#include "TSystem.h"
#include "TInterpreter.h"
#include "TList.h"
#include "TRandom.h"
#include "TLorentzVector.h"
#include "../interface/Reader.h"

using namespace PDGHelpers;

Reader::Reader(OptionParser* options_) : converter(options_){
  configure();
  run();
}

void Reader::configure(){}
void Reader::finalizeRun(){}

template<typename returnType> bool Reader::setVariable(const Event* ev, string& branchname, returnType(*evalVar)(const Event*, string&)){
  returnType result = evalVar(ev, branchname);
  tree->setVal(branchname, result);
  return true;
}

void Reader::bindInputBranches(HVVTree* tin){
  vector<string> inputBranches = tin->getBranchList();
  for (int br=0; br<inputBranches.size(); br++){
    int pos=-1;
    string branchname = inputBranches.at(br);
    bool outfound=false;
    BaseTree::BranchTypes oBT = tree->searchArray(branchname, pos);
    if (!(pos==-1 || oBT==BaseTree::nBranchTypes)) outfound = true;
    pos=-1;
    BaseTree::BranchTypes iBT = tin->searchArray(branchname, pos);
    if (iBT!=oBT) outfound = false;
    if (!outfound) continue;
    else{
      void* iBVRef = tin->getBranchHandleRef(branchname);
      void* oBVRef = tree->getBranchHandleRef(branchname);
      if (iBVRef==0 || oBVRef==0) continue;
      if (iBT == BaseTree::bInt){
        Int_t* iBRef = (Int_t*)iBVRef;
        Int_t* oBRef = (Int_t*)oBVRef;
        pair<Int_t*, Int_t*> refpair(iBRef, oBRef);
        intBranchMap.push_back(refpair);
      }
      else if (iBT == BaseTree::bFloat){
        Float_t* iBRef = (Float_t*)iBVRef;
        Float_t* oBRef = (Float_t*)oBVRef;
        pair<Float_t*, Float_t*> refpair(iBRef, oBRef);
        floatBranchMap.push_back(refpair);
      }
      else if (iBT == BaseTree::bVectorInt){
        vector<int>** iBRef = (vector<int>**)iBVRef;
        vector<int>* oBRef = (vector<int>*)oBVRef;
        pair<vector<int>**, vector<int>*> refpair(iBRef, oBRef);
        vectorIntBranchMap.push_back(refpair);
      }
      else if (iBT == BaseTree::bVectorDouble){
        vector<double>** iBRef = (vector<double>**)iBVRef;
        vector<double>* oBRef = (vector<double>*)oBVRef;
        pair<vector<double>**, vector<double>*> refpair(iBRef, oBRef);
        vectorDoubleBranchMap.push_back(refpair);
      }
    }
  }
}
void Reader::resetBranchBinding(){
  intBranchMap.clear();
  floatBranchMap.clear();
  vectorIntBranchMap.clear();
  vectorDoubleBranchMap.clear();
}


void Reader::synchMappedBranches(){
  for (int b=0; b<intBranchMap.size(); b++) *(intBranchMap.at(b).second) = *(intBranchMap.at(b).first);
  for (int b=0; b<floatBranchMap.size(); b++) *(floatBranchMap.at(b).second) = *(floatBranchMap.at(b).first);
  for (int b=0; b<vectorIntBranchMap.size(); b++){
    vector<int>* inhandle = *(vectorIntBranchMap.at(b).first);
    vector<int>* outhandle = vectorIntBranchMap.at(b).second;
    for (int el=0; el<inhandle->size(); el++){
      outhandle->push_back(inhandle->at(el));
    }
  }
  for (int b=0; b<vectorDoubleBranchMap.size(); b++){
    vector<double>* inhandle = *(vectorDoubleBranchMap.at(b).first);
    vector<double>* outhandle = vectorDoubleBranchMap.at(b).second;
    for (int el=0; el<inhandle->size(); el++){
      outhandle->push_back(inhandle->at(el));
    }
  }
}



void Reader::run(){/*
  Float_t MC_weight=0;
  Int_t isSelected=0;

  tree->bookAllBranches(false);

  for (int f=0; f<filename.size(); f++){
    string cinput = filename.at(f);
    cout << "Processing " << cinput << "..." << endl;
    TFile* fin = 0;
    fin = getIntermediateFile(cinput);
    if (fin!=0 && fin->IsZombie()){
      if (fin->IsOpen()) fin->Close();
      delete fin;
      fin=0;
    }
    else if (fin!=0 && !fin->IsOpen()){ // Separate if-statement on purpose
      delete fin;
      fin=0;
    }
    else if (fin!=0){
      TTree* tin = (TTree*)fin->Get("tmpTree");
      foutput->cd();

      int nInputEvents = tin->GetEntries();
      for (int ev=0; ev<nInputEvents; ev++){
        double weight;
        bool genSuccess=false, smearedSuccess=false;

        vector<Particle*> genParticleList;
        vector<Particle*> smearedParticleList;
        vector<ZZCandidate*> genCandList; // Bookkeeping
        vector<ZZCandidate*> smearedCandList; // Bookkeeping

        tree->initializeBranches();

        readEvent(tin, ev, genParticleList, genSuccess, smearedParticleList, smearedSuccess, weight);

        if (weight!=0){
          MC_weight = (float)weight;

          Event genEvent;
          if (genSuccess){
            genEvent.setWeight(weight);
            vector<int> hasGenHiggs;
            for (int p=0; p<genParticleList.size(); p++){
              Particle* genPart = genParticleList.at(p); // Has mother info from Pythia reading
              if (isAHiggs(genPart->id)) hasGenHiggs.push_back(p);

              if (genPart->genStatus==1){
                if (isALepton(genPart->id)) genEvent.addLepton(genPart);
                else if (isANeutrino(genPart->id)) genEvent.addNeutrino(genPart);
                else if (isAGluon(genPart->id) || isAQuark(genPart->id)) genEvent.addJet(genPart);
              }
              else if (genPart->genStatus==-1) tree->fillMotherInfo(genPart);
            }

            genEvent.constructVVCandidates(options->doGenHZZdecay(), options->genDecayProducts());
            genEvent.addVVCandidateAppendages();
            ZZCandidate* genCand=0;
            if (hasGenHiggs.size()>0){
              for (int gk=0; gk<hasGenHiggs.size(); gk++){
                ZZCandidate* tmpCand = HiggsComparators::matchAHiggsToParticle(genEvent, genParticleList.at(hasGenHiggs.at(gk)));
                if (tmpCand!=0){
                  if (genCand==0) genCand=tmpCand;
                  else genCand = HiggsComparators::candComparator(genCand, tmpCand, options->getHiggsCandidateSelectionScheme(true), options->doGenHZZdecay());
                }
              }
            }
            else genCand = HiggsComparators::candidateSelector(genEvent, options->getHiggsCandidateSelectionScheme(true), options->doGenHZZdecay());
            if (genCand!=0) tree->fillCandidate(genCand, true);
            else cout << cinput << " (" << ev << "): No gen. level Higgs candidate was found!" << endl;
          }

          Event smearedEvent;
          if (smearedSuccess){
            smearedEvent.setWeight(weight);
            for (int p=0; p<smearedParticleList.size(); p++){
              Particle* smearedPart = smearedParticleList.at(p);
              if (isALepton(smearedPart->id)) smearedEvent.addLepton(smearedPart);
              else if (isANeutrino(smearedPart->id)) smearedEvent.addNeutrino(smearedPart);
              else if (smearedPart->id==0) smearedEvent.addJet(smearedPart);
              else smearedEvent.addParticle(smearedPart);
            }
            smearedEvent.constructVVCandidates(options->doRecoHZZdecay(), options->recoDecayProducts());
            smearedEvent.applyParticleSelection();
            smearedEvent.addVVCandidateAppendages();
            ZZCandidate* rCand = HiggsComparators::candidateSelector(smearedEvent, options->getHiggsCandidateSelectionScheme(false), options->doRecoHZZdecay());
            if (rCand!=0){
              isSelected=1;
              tree->fillCandidate(rCand, false);
            }
            else isSelected=0;
          }
        }

        tree->fillEventVariables(MC_weight, isSelected);
        if ((genSuccess || smearedSuccess) && weight!=0) tree->record();

        for (int p=0; p<smearedCandList.size(); p++){ // Bookkeeping
          ZZCandidate* tmpCand = (ZZCandidate*)smearedCandList.at(p);
          if (tmpCand!=0) delete tmpCand;
        }
        for (int p=0; p<smearedParticleList.size(); p++){ // Bookkeeping
          Particle* tmpPart = (Particle*)smearedParticleList.at(p);
          if (tmpPart!=0) delete tmpPart;
        }

        for (int p=0; p<genCandList.size(); p++){ // Bookkeeping
          ZZCandidate* tmpCand = (ZZCandidate*)genCandList.at(p);
          if (tmpCand!=0) delete tmpCand;
        }
        for (int p=0; p<genParticleList.size(); p++){ // Bookkeeping
          Particle* tmpPart = (Particle*)genParticleList.at(p);
          if (tmpPart!=0) delete tmpPart;
        }

        // Bookkeeping
        smearedCandList.clear();
        smearedParticleList.clear();
        genCandList.clear();
        genParticleList.clear();
      }

      fin->Close();
    }
  }
  finalizeRun();
*/}



void Reader::readEvent(TTree* tin, int ev, bool isGen, Event& outEvent){/*
  int nEvents = tin->GetEntries();
  vector<double> weights;
  if (ev>=nEvents){
    double weight = 0;
    genSuccess=false;
    smearedSuccess=false;
    weights.push_back(0);
  }
  else{
    vector<double>* geneventinfoweights=0;
    TBranch* b_geneventinfoweights=0;

    vector<double>* reco_GenParticle_FV[4]={ 0 };
    vector<int>* reco_GenParticle_id=0;
    vector<int>* reco_GenParticle_status=0;
    TBranch* b_reco_GenParticle_FV[4]={ 0 };
    TBranch* b_reco_GenParticle_id=0;
    TBranch* b_reco_GenParticle_status=0;

    vector<double>* reco_GenJet_FV[4]={ 0 };
    vector<int>* reco_GenJet_id=0;
    vector<int>* reco_GenJet_status=0;
    TBranch* b_reco_GenJet_FV[4]={ 0 };
    TBranch* b_reco_GenJet_id=0;
    TBranch* b_reco_GenJet_status=0;

    tin->SetBranchAddress("genWeights", &geneventinfoweights, &b_geneventinfoweights);

    tin->SetBranchAddress("reco_GenParticle_X", reco_GenParticle_FV, b_reco_GenParticle_FV);
    tin->SetBranchAddress("reco_GenParticle_Y", reco_GenParticle_FV+1, b_reco_GenParticle_FV+1);
    tin->SetBranchAddress("reco_GenParticle_Z", reco_GenParticle_FV+2, b_reco_GenParticle_FV+2);
    tin->SetBranchAddress("reco_GenParticle_E", reco_GenParticle_FV+3, b_reco_GenParticle_FV+3);
    tin->SetBranchAddress("reco_GenParticle_id", &reco_GenParticle_id, &b_reco_GenParticle_id);
    tin->SetBranchAddress("reco_GenParticle_status", &reco_GenParticle_status, &b_reco_GenParticle_status);

    tin->SetBranchAddress("reco_GenJet_X", reco_GenJet_FV, b_reco_GenJet_FV);
    tin->SetBranchAddress("reco_GenJet_Y", reco_GenJet_FV+1, b_reco_GenJet_FV+1);
    tin->SetBranchAddress("reco_GenJet_Z", reco_GenJet_FV+2, b_reco_GenJet_FV+2);
    tin->SetBranchAddress("reco_GenJet_E", reco_GenJet_FV+3, b_reco_GenJet_FV+3);
    tin->SetBranchAddress("reco_GenJet_id", &reco_GenJet_id, &b_reco_GenJet_id);
    tin->SetBranchAddress("reco_GenJet_status", &reco_GenJet_status, &b_reco_GenJet_status);

    tin->GetEntry(ev);

    // Gen. particles
    int motherID[2];
    int mctr=0;
    for (int a = 0; a < reco_GenParticle_id->size(); a++){
      int istup = reco_GenParticle_status->at(a);
      if (istup==21 && mctr<2){
        motherID[mctr] = a;
        mctr++;
      }
      int idup = reco_GenParticle_id->at(a);
      TLorentzVector partFourVec(reco_GenParticle_FV[0]->at(a), reco_GenParticle_FV[1]->at(a), reco_GenParticle_FV[2]->at(a), reco_GenParticle_FV[3]->at(a));

      Particle* onePart = new Particle(idup, partFourVec);
      onePart->setGenStatus(PDGHelpers::ReaderStatus(istup));
      onePart->setLifetime(0);
      genCollection.push_back(onePart);
    }
    // Assign the mothers
    for (int a = 0; a < genCollection.size(); a++){
      for(int m=0;m<2;m++) genCollection.at(a)->addMother(genCollection.at(motherID[m]));
    }
    genSuccess=(genCollection.size()>0);
    // Reco. particles
    for (int a = 0; a < reco_GenJet_id->size(); a++){
      int istup = reco_GenJet_status->at(a);
      int idup = reco_GenJet_id->at(a);
      TLorentzVector partFourVec(reco_GenJet_FV[0]->at(a), reco_GenJet_FV[1]->at(a), reco_GenJet_FV[2]->at(a), reco_GenJet_FV[3]->at(a));

      Particle* onePart = new Particle(idup, partFourVec);
      onePart->setGenStatus(PDGHelpers::ReaderStatus(istup));
      onePart->setLifetime(0);
      recoCollection.push_back(onePart);
    }
    smearedSuccess=(recoCollection.size()>0);

    tin->ResetBranchAddresses();

    if (geneventinfoweights!=0 && geneventinfoweights->size()>0){
      for (int w=0; w<geneventinfoweights->size(); w++) weights.push_back(geneventinfoweights->at(w));
    }
    else if (!smearedSuccess && !genSuccess) weights.push_back(0);
    else weights.push_back(1.);
  }
  eventWeight = weights.at(0);
*/}
