#include "TSystem.h"
#include "TInterpreter.h"
#include "TList.h"
#include "TRandom.h"
#include "TLorentzVector.h"
#include "Reader.h"

using namespace PDGHelpers;

Reader::Reader(OptionParser* options_) : converter(options_){
  configure();
  run();
}

void Reader::configure(){

}
void Reader::finalizeRun(){}

template<typename returnType> bool Reader::setVariable(const MELAEvent* ev, string& branchname, returnType(*evalVar)(const MELAEvent*, string&)){
  returnType result = evalVar(ev, branchname);
  returnType* resultPtr = &result;
  if (dynamic_cast<vectorInt*>(resultPtr) || dynamic_cast<vectorDouble*>(resultPtr) || dynamic_cast<vectorFloat*>(resultPtr)){
    for (unsigned int el=0; el<result.size(); el++) tree->setVal(branchname, result.at(el));
    return true;
  }
  else if (dynamic_cast<Int_t*>(resultPtr) || dynamic_cast<Float_t*>(resultPtr)){
    tree->setVal(branchname, result);
    return true;
  }
  else return false;
}

void Reader::bindInputBranches(HVVTree* tin){
  vector<string> inputBranches = tin->getBranchList();
  vector<pair<string, BaseTree::BranchTypes>> unreservedBranches;
  for (unsigned int br=0; br<inputBranches.size(); br++){
    int pos=-1;
    string branchname = inputBranches.at(br);

    BaseTree::BranchTypes iBT = tin->searchArray(branchname, pos);
    pos=-1;
    BaseTree::BranchTypes oBT = tree->searchArray(branchname, pos);
    bool outfound = !(pos==-1 || oBT==BaseTree::nBranchTypes || iBT!=oBT);
    if (!outfound && tree->getBranchList().empty()) unreservedBranches.push_back(pair<string, BaseTree::BranchTypes>(branchname, iBT));
  }
  if (!unreservedBranches.empty()){
    for (pair<string, BaseTree::BranchTypes> const& tmp_branch:unreservedBranches) tree->reserveBranch(tmp_branch.first, tmp_branch.second, false);
    tree->actuateBranches(false);
  }

  for (unsigned int br=0; br<inputBranches.size(); br++){
    int pos=-1;
    string branchname = inputBranches.at(br);

    BaseTree::BranchTypes iBT = tin->searchArray(branchname, pos);
    pos=-1;
    BaseTree::BranchTypes oBT = tree->searchArray(branchname, pos);
    bool outfound = !(pos==-1 || oBT==BaseTree::nBranchTypes || iBT!=oBT);
    if (!outfound) continue; // If still not found, forget it!
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
      else if (iBT == BaseTree::bVectorFloat){
        vectorFloat** iBRef = (vectorFloat**) iBVRef;
        vectorFloat** oBRef = (vectorFloat**) oBVRef;
        pair<vectorFloat**, vectorFloat**> refpair(iBRef, oBRef);
        vectorFloatBranchMap.push_back(refpair);
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
  vectorFloatBranchMap.clear();
  vectorDoubleBranchMap.clear();
}
void Reader::synchMappedBranches(){
  for (unsigned int b=0; b<intBranchMap.size(); b++) *(intBranchMap.at(b).second) = *(intBranchMap.at(b).first);
  for (unsigned int b=0; b<floatBranchMap.size(); b++) *(floatBranchMap.at(b).second) = *(floatBranchMap.at(b).first);
  for (unsigned int b=0; b<vectorIntBranchMap.size(); b++){
    vectorInt* inhandle = *(vectorIntBranchMap.at(b).first);
    vectorInt* outhandle = *(vectorIntBranchMap.at(b).second);
    for (unsigned int el=0; el<inhandle->size(); el++){
      outhandle->push_back(inhandle->at(el));
    }
  }
  for (unsigned int b=0; b<vectorFloatBranchMap.size(); b++){
    vectorFloat* inhandle = *(vectorFloatBranchMap.at(b).first);
    vectorFloat* outhandle = *(vectorFloatBranchMap.at(b).second);
    for (unsigned int el=0; el<inhandle->size(); el++){
      outhandle->push_back(inhandle->at(el));
    }
  }
  for (unsigned int b=0; b<vectorDoubleBranchMap.size(); b++){
    vectorDouble* inhandle = *(vectorDoubleBranchMap.at(b).first);
    vectorDouble* outhandle = *(vectorDoubleBranchMap.at(b).second);
    for (unsigned int el=0; el<inhandle->size(); el++){
      outhandle->push_back(inhandle->at(el));
    }
  }
}

void Reader::run(){
  int globalNEvents = 0;
  int maxProcEvents = options->maxEventsToProcess();
  vector < pair<Int_t, Int_t> > eventSkipList = options->getSkippedEvents();

  bool firstFile = true;
  for (string const& cinput:filename){
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
      int nProcessed = 0;

      HVVTree* tin = new HVVTree("SelectedTree", fin);
      if (tin->getTree()){
        tin->setOptions(options);
        tin->bookAllBranches(true);
        cout << "Input tree branches booked...";
        bindInputBranches(tin);
        cout << "bound..." << endl;

        foutput->cd();
        if (firstFile){
          tree->bookAllBranches(false);
          firstFile=false;
        }

        int nInputEvents = tin->getTree()->GetEntries();
        cout << "Number of input events to process: " << nInputEvents << endl;
        for (int ev=0; ev<nInputEvents; ev++){
          if (globalNEvents>=maxProcEvents && maxProcEvents>=0) break;
          bool doSkipEvent = false;
          for (unsigned int es=0; es<eventSkipList.size(); es++){
            if (
              (eventSkipList.at(es).first<=globalNEvents && eventSkipList.at(es).second>=globalNEvents)
              ||
              (eventSkipList.at(es).first<=globalNEvents && eventSkipList.at(es).second<0)
              )doSkipEvent=true;
          }
          if (doSkipEvent){ globalNEvents++; continue; }

          vector<MELAParticle*> genParticleList;
          vector<MELAParticle*> recoParticleList;
          vector<MELACandidate*> genCandList; // Bookkeeping
          vector<MELACandidate*> recoCandList; // Bookkeeping

          tree->initializeBranches();

          tin->getTree()->GetEntry(ev);
          synchMappedBranches();

          MELAEvent genEvent, recoEvent;
          MELACandidate* genCand=0;
          MELACandidate* recoCand=0;
          if (options->processGenInfo()){
            readEvent(genEvent, genParticleList, true);
            if (genEvent.getNCandidates()<=1) genCand = genEvent.getCandidate(0);
            else genCand = HiggsComparators::candidateSelector(genEvent, options->getHiggsCandidateSelectionScheme(true), options->doGenHZZdecay());
          }
          if (options->processRecoInfo()){
            readEvent(recoEvent, recoParticleList, false);
            if (recoEvent.getNCandidates()<=1) recoCand = recoEvent.getCandidate(0);
            else recoCand = HiggsComparators::candidateSelector(recoEvent, options->getHiggsCandidateSelectionScheme(false), options->doRecoHZZdecay());
          }

          if (melaHelpers::melaHandle && genCand){
            melaHelpers::melaHandle->setCurrentCandidate(genCand);

            if (options->doComputeDecayAngles()) tree->fillDecayAngles(true);
            if (options->doComputeVBFAngles()) tree->fillVBFProductionAngles(true);
            if (options->doComputeVHAngles()) tree->fillVHProductionAngles(true);
            tree->fillMELAProbabilities(true);

            melaHelpers::melaHandle->resetInputEvent();
          }
          if (recoCand){
            if (options->recoSelectionMode()!=0){
              tree->fillCandidate(recoCand, false);
              tree->fillEventVariables(*((Float_t*)tree->getBranchHandleRef("MC_weight")), 0 /*isSelected*/);
            }
            else if (melaHelpers::melaHandle){
              melaHelpers::melaHandle->setCurrentCandidate(recoCand);

              if (options->doComputeDecayAngles()) tree->fillDecayAngles(false);
              if (options->doComputeVBFAngles()) tree->fillVBFProductionAngles(false);
              if (options->doComputeVHAngles()) tree->fillVHProductionAngles(false);
              tree->fillMELAProbabilities(false); // Do it at the last step

              melaHelpers::melaHandle->resetInputEvent();
            }
          }
          else if (options->recoSelectionMode()!=0) tree->fillEventVariables(*((Float_t*)tree->getBranchHandleRef("MC_weight")), 0 /*isSelected*/);

          if ((recoCand && options->processRecoInfo()) || (genCand && options->processGenInfo())){
            tree->record();
            nProcessed++;
          }

          for (unsigned int p=0; p<recoCandList.size(); p++){ // Bookkeeping
            MELACandidate* tmpCand = (MELACandidate*)recoCandList.at(p);
            if (tmpCand!=0) delete tmpCand;
          }
          for (unsigned int p=0; p<recoParticleList.size(); p++){ // Bookkeeping
            MELAParticle* tmpPart = (MELAParticle*)recoParticleList.at(p);
            if (tmpPart!=0) delete tmpPart;
          }

          for (unsigned int p=0; p<genCandList.size(); p++){ // Bookkeeping
            MELACandidate* tmpCand = (MELACandidate*)genCandList.at(p);
            if (tmpCand!=0) delete tmpCand;
          }
          for (unsigned int p=0; p<genParticleList.size(); p++){ // Bookkeeping
            MELAParticle* tmpPart = (MELAParticle*)genParticleList.at(p);
            if (tmpPart!=0) delete tmpPart;
          }

          // Bookkeeping
          recoCandList.clear();
          recoParticleList.clear();
          genCandList.clear();
          genParticleList.clear();

          globalNEvents++;
        }
        resetBranchBinding();
        cout << "Processed number of events from the input file (recorded events / sample size observed / cumulative traversed): " << nProcessed << " / " << nInputEvents << " / " << globalNEvents << endl;
      }
      delete tin;
      fin->Close();
    }
  }
  finalizeRun();
}

void Reader::readEvent(MELAEvent& outEvent, vector<MELAParticle*>& particles, bool isGen){
  typedef vectorFloat VectorPrecision_t;

  string varname;

  string strId = "Id";
  vector<string> strMassPtEtaPhi; strMassPtEtaPhi.push_back("Pt"); strMassPtEtaPhi.push_back("Eta"); strMassPtEtaPhi.push_back("Phi"); strMassPtEtaPhi.push_back("Mass");
  vector<string> strMassPtPzPhi; strMassPtPzPhi.push_back("Pt"); strMassPtPzPhi.push_back("Pz"); strMassPtPzPhi.push_back("Phi"); strMassPtPzPhi.push_back("Mass");
  vector<string> strPzE; strPzE.push_back("Pz"); strPzE.push_back("E");

  string strLepCore = "CandDau";
  if (isGen) strLepCore.insert(0, "Gen");

  // Search the tree for the Higgs daughters
  vector<MELAParticle*> candFinalDaughters;
  VectorPrecision_t** ref_dauFV[4];
  bool dau_success=true;
  for (int fv=0; fv<4; fv++){
    varname = strLepCore + strMassPtEtaPhi.at(fv);
    ref_dauFV[fv] = (VectorPrecision_t**) tree->getBranchHandleRef(varname);
    if (!ref_dauFV[fv]) { dau_success=false; break; }
    else if (!(*ref_dauFV[fv])) { dau_success=false; break; }
  }
  varname = strLepCore + strId;
  vectorInt** ref_dauId = (vectorInt**) tree->getBranchHandleRef(varname);
  if (!ref_dauId) dau_success=false;
  if (!(*ref_dauId)) dau_success=false;
  if (dau_success){
    vectorInt* idhandle = *ref_dauId;
    VectorPrecision_t* fvhandle[4];

    for (int fv=0; fv<4; fv++){
      fvhandle[fv]=*(ref_dauFV[fv]);
      if (fvhandle[fv]->size()!=idhandle->size()){ dau_success=false; break; }
    }
    for (unsigned int el=0; el<idhandle->size(); el++){
      Int_t partId = idhandle->at(el);
      TLorentzVector partFV; partFV.SetPtEtaPhiM(fvhandle[0]->at(el), fvhandle[1]->at(el), fvhandle[2]->at(el), fvhandle[3]->at(el));
      if (partId==-9000 || (partFV.P()==0. && partFV.T()==0.)) continue; // Invalid particle
      MELAParticle* part = new MELAParticle((int) partId, partFV);
      if (isGen) part->setGenStatus(1);
      particles.push_back(part);
      candFinalDaughters.push_back(part);
    }
  }

  // Search the tree for the associated particles
  vector<MELAParticle*> associatedParticles;
  string strAPCore = "AssociatedParticle";
  if (isGen) strAPCore.insert(0, "Gen");
  VectorPrecision_t** ref_apartFV[4];
  bool apart_success=true;
  for (int fv=0; fv<4; fv++){
    varname = strAPCore + strMassPtEtaPhi.at(fv);
    ref_apartFV[fv] = (VectorPrecision_t**) tree->getBranchHandleRef(varname);
    if (!ref_apartFV[fv]) { apart_success=false; break; }
    else if (!(*ref_apartFV[fv])) { apart_success=false; break; }
  }
  varname = strAPCore + strId;
  vectorInt** ref_apartId = (vectorInt**) tree->getBranchHandleRef(varname);
  if (!ref_apartId) apart_success=false;
  if (!(*ref_apartId)) apart_success=false;
  if (apart_success){
    vectorInt* idhandle = *ref_apartId;
    VectorPrecision_t* fvhandle[4];

    for (int fv=0; fv<4; fv++){
      fvhandle[fv]=*(ref_apartFV[fv]);
      if (fvhandle[fv]->size()!=idhandle->size()) { apart_success=false; break; }
    }
    for (unsigned int el=0; el<idhandle->size(); el++){
      Int_t partId = idhandle->at(el);
      TLorentzVector partFV; partFV.SetPtEtaPhiM(fvhandle[0]->at(el), fvhandle[1]->at(el), fvhandle[2]->at(el), fvhandle[3]->at(el));
      if (partId==-9000 || (partFV.P()==0. && partFV.T()==0.)) continue; // Invalid particle
      MELAParticle* part = new MELAParticle((int) partId, partFV);
      if (isGen) part->setGenStatus(1);
      particles.push_back(part);
      associatedParticles.push_back(part);
    }
  }

  // Reconstruct the gen and reco candidates
  if (isGen){ // Do not find the best combination, leave it to the input tree

    // First search the tree for the mothers
    vector<MELAParticle*> motherParticles;
    string strMPCore = "Mother";
    strMPCore.insert(0, "Gen");
    VectorPrecision_t** ref_mpartFV[2];
    bool mpart_success=true;
    for (int fv=0; fv<2; fv++){
      varname = strMPCore + strPzE.at(fv);
      ref_mpartFV[fv] = (VectorPrecision_t**) tree->getBranchHandleRef(varname);
      if (!ref_mpartFV[fv]) { mpart_success=false; break; }
      else if (!(*ref_mpartFV[fv])) { mpart_success=false; break; }
    }
    varname = strMPCore + strId;
    vectorInt** ref_mpartId = (vectorInt**) tree->getBranchHandleRef(varname);
    if (!ref_mpartId) mpart_success=false;
    if (!(*ref_mpartId)) mpart_success=false;
    if (mpart_success){
      vectorInt* idhandle = *ref_mpartId;
      VectorPrecision_t* fvhandle[2];

      for (int fv=0; fv<2; fv++){
        fvhandle[fv]=*(ref_mpartFV[fv]);
        if (fvhandle[fv]->size()!=idhandle->size()) { mpart_success=false; break; }
      }
      for (unsigned int el=0; el<idhandle->size(); el++){
        Int_t partId = idhandle->at(el);
        TLorentzVector partFV; partFV.SetXYZT(0, 0, fvhandle[0]->at(el), fvhandle[1]->at(el));
        MELAParticle* part = new MELAParticle((int) partId, partFV);
        part->setGenStatus(-1);
        particles.push_back(part);
        motherParticles.push_back(part);
      }
    }

    for (unsigned int d=0; d<candFinalDaughters.size(); d++){
      MELAParticle* part = candFinalDaughters.at(d);
      if (isALepton(part->id)) outEvent.addLepton(part);
      else if (isANeutrino(part->id)) outEvent.addNeutrino(part);
      else if (isAPhoton(part->id)) outEvent.addPhoton(part);
      else if (isAGluon(part->id) || isAQuark(part->id)) outEvent.addJet(part);
      else if (part->id==0) outEvent.addJet(part);
      else outEvent.addParticle(part);
    }
    outEvent.constructVVCandidates(options->doGenHZZdecay(), options->genDecayProducts());
    for (unsigned int p=0; p<motherParticles.size(); p++){
      MELAParticle* part = motherParticles.at(p);
      outEvent.addMother(part);
    }
    for (unsigned int d=0; d<associatedParticles.size(); d++){
      MELAParticle* part = associatedParticles.at(d);
      if (isALepton(part->id)) outEvent.addLepton(part);
      else if (isANeutrino(part->id)) outEvent.addNeutrino(part);
      else if (isAPhoton(part->id)) outEvent.addPhoton(part);
      else if (isAGluon(part->id) || isAQuark(part->id)) outEvent.addJet(part);
      else if (part->id==0) outEvent.addJet(part);
      else outEvent.addParticle(part);
    }
  }
  else{
    for (unsigned int d=0; d<candFinalDaughters.size(); d++){
      MELAParticle* part = candFinalDaughters.at(d);
      if (isALepton(part->id)) outEvent.addLepton(part);
      else if (isANeutrino(part->id)) outEvent.addNeutrino(part);
      else if (isAPhoton(part->id)) outEvent.addPhoton(part);
      else if (isAGluon(part->id) || isAQuark(part->id)) outEvent.addJet(part);
      else if (part->id==0) outEvent.addJet(part);
      else outEvent.addParticle(part);
    }

    if (options->recoSelectionMode()==0){ // Do not find the best combination, leave it to the input tree
      outEvent.constructVVCandidates(options->doRecoHZZdecay(), options->recoDecayProducts());
      for (unsigned int d=0; d<associatedParticles.size(); d++){
        MELAParticle* part = associatedParticles.at(d);
        if (isALepton(part->id)) outEvent.addLepton(part);
        else if (isANeutrino(part->id)) outEvent.addNeutrino(part);
        else if (isAPhoton(part->id)) outEvent.addPhoton(part);
        else if (isAGluon(part->id) || isAQuark(part->id)) outEvent.addJet(part);
        else if (part->id==0) outEvent.addJet(part);
        else outEvent.addParticle(part);
      }
    }
    else{ // Find the best combination
      for (unsigned int d=0; d<associatedParticles.size(); d++){
        MELAParticle* part = associatedParticles.at(d);
        if (isALepton(part->id)) outEvent.addLepton(part);
        else if (isANeutrino(part->id)) outEvent.addNeutrino(part);
        else if (isAPhoton(part->id)) outEvent.addPhoton(part);
        else if (isAGluon(part->id) || isAQuark(part->id)) outEvent.addJet(part);
        else if (part->id==0) outEvent.addJet(part);
        else outEvent.addParticle(part);
      }
      outEvent.constructVVCandidates(options->doRecoHZZdecay(), options->recoDecayProducts());
      outEvent.applyParticleSelection();
    }
  }
  outEvent.addVVCandidateAppendages();

  /*
    if (isGen) cout << "NGenCandidates: " << outEvent.getNCandidates() << endl;
    else cout << "NRecoCandidates: " << outEvent.getNCandidates() << endl;
    for (int cc=0; cc<outEvent.getNCandidates(); cc++){
    cout << outEvent.getMELACandidate(cc)->m()
    << " ("
    << outEvent.getMELACandidate(cc)->getSortedV(0)->m()
    << " -> "
    << "pT = " << outEvent.getMELACandidate(cc)->getSortedV(0)->getDaughter(0)->pt() << " + " << outEvent.getMELACandidate(cc)->getSortedV(0)->getDaughter(1)->pt()
    << ", "
    << outEvent.getMELACandidate(cc)->getSortedV(1)->m()
    << " -> "
    << "pT = " << outEvent.getMELACandidate(cc)->getSortedV(1)->getDaughter(0)->pt() << " + " << outEvent.getMELACandidate(cc)->getSortedV(1)->getDaughter(1)->pt()
    << ")\n";
    }
    if (isGen) cout << "Recorded (m1, m2): (" << *((Float_t*)tree->getBranchHandleRef("GenZ1Mass")) << ", " << *((Float_t*)tree->getBranchHandleRef("GenZ2Mass")) << ")\n";
    else cout << "Recorded (m1, m2): (" << *((Float_t*)tree->getBranchHandleRef("Z1Mass")) << ", " << *((Float_t*)tree->getBranchHandleRef("Z2Mass")) << ")\n";
    cout << endl;
  */
}

