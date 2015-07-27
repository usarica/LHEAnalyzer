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
  returnType* resultPtr = &result;
  if (dynamic_cast<vectorInt*>(resultPtr)!=0 || dynamic_cast<vectorDouble*>(resultPtr)!=0){
    for (int el=0; el<result.size(); el++) tree->setVal(branchname, result.at(el));
    return true;
  }
  else if (dynamic_cast<Int_t*>(resultPtr)!=0 || dynamic_cast<Float_t*>(resultPtr)!=0){
    tree->setVal(branchname, result);
    return true;
  }
  else return false;
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
        vectorInt** iBRef = (vectorInt**)iBVRef;
        vectorInt** oBRef = (vectorInt**)oBVRef;
        pair<vectorInt**, vectorInt**> refpair(iBRef, oBRef);
        vectorIntBranchMap.push_back(refpair);
      }
      else if (iBT == BaseTree::bVectorDouble){
        vectorDouble** iBRef = (vectorDouble**)iBVRef;
        vectorDouble** oBRef = (vectorDouble**)oBVRef;
        pair<vectorDouble**, vectorDouble**> refpair(iBRef, oBRef);
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
    vectorInt* inhandle = *(vectorIntBranchMap.at(b).first);
    vectorInt* outhandle = *(vectorIntBranchMap.at(b).second);
    for (int el=0; el<inhandle->size(); el++){
      outhandle->push_back(inhandle->at(el));
    }
  }
  for (int b=0; b<vectorDoubleBranchMap.size(); b++){
    vectorDouble* inhandle = *(vectorDoubleBranchMap.at(b).first);
    vectorDouble* outhandle = *(vectorDoubleBranchMap.at(b).second);
    for (int el=0; el<inhandle->size(); el++){
      outhandle->push_back(inhandle->at(el));
    }
  }
}



void Reader::run(){
  Float_t MC_weight=0;
  Int_t isSelected=0;

  tree->bookAllBranches(false);

  for (int f=0; f<filename.size(); f++){
    string cinput = filename.at(f);
    cout << "Processing " << cinput << "..." << endl;
    TFile* fin = new TFile(cinput.c_str(), "read");
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
      HVVTree* tin = new HVVTree("SelectedTree", fin);
      if (tin->getTree()!=0){
        tin->setOptions(options);
        tin->bookAllBranches(true);
        cout << "Input tree branches booked...";
        bindInputBranches(tin);
        cout << "bound..." << endl;
        foutput->cd();

        int nInputEvents = tin->getTree()->GetEntries();
        cout << "Number of input events to process: " << nInputEvents << endl;
        for (int ev=0; ev<nInputEvents; ev++){
          vector<Particle*> genParticleList;
          vector<Particle*> recoParticleList;
          vector<ZZCandidate*> genCandList; // Bookkeeping
          vector<ZZCandidate*> recoCandList; // Bookkeeping

          tree->initializeBranches();

          tin->getTree()->GetEntry(ev);
          synchMappedBranches();

          Event genEvent, recoEvent;

          readEvent(genEvent, true);
          readEvent(recoEvent, false);


          // Do magic


          tree->record();

          for (int p=0; p<recoCandList.size(); p++){ // Bookkeeping
            ZZCandidate* tmpCand = (ZZCandidate*)recoCandList.at(p);
            if (tmpCand!=0) delete tmpCand;
          }
          for (int p=0; p<recoParticleList.size(); p++){ // Bookkeeping
            Particle* tmpPart = (Particle*)recoParticleList.at(p);
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
          recoCandList.clear();
          recoParticleList.clear();
          genCandList.clear();
          genParticleList.clear();
        }
        resetBranchBinding();
      }
      delete tin;
      fin->Close();
    }
  }
  finalizeRun();
}



void Reader::readEvent(Event& outEvent, bool isGen){/*
  int nEvents = tin->GetEntries();
  vectorDouble weights;
  if (ev>=nEvents){
    double weight = 0;
    genSuccess=false;
    recoSuccess=false;
    weights.push_back(0);
  }
  else{
    vectorDouble* geneventinfoweights=0;
    TBranch* b_geneventinfoweights=0;

    vectorDouble* reco_GenParticle_FV[4]={ 0 };
    vectorInt* reco_GenParticle_id=0;
    vectorInt* reco_GenParticle_status=0;
    TBranch* b_reco_GenParticle_FV[4]={ 0 };
    TBranch* b_reco_GenParticle_id=0;
    TBranch* b_reco_GenParticle_status=0;

    vectorDouble* reco_GenJet_FV[4]={ 0 };
    vectorInt* reco_GenJet_id=0;
    vectorInt* reco_GenJet_status=0;
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
    recoSuccess=(recoCollection.size()>0);

    tin->ResetBranchAddresses();

    if (geneventinfoweights!=0 && geneventinfoweights->size()>0){
      for (int w=0; w<geneventinfoweights->size(); w++) weights.push_back(geneventinfoweights->at(w));
    }
    else if (!recoSuccess && !genSuccess) weights.push_back(0);
    else weights.push_back(1.);
  }
  eventWeight = weights.at(0);
*/}

