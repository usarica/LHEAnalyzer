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
#include "TRandom.h"
#include "TLorentzVector.h"
#include "LHEConverter.h"
#include "MELAJetMerger.h"
#include "MELAStreamHelpers.hh"


using namespace std;
using namespace PDGHelpers;
using namespace LHEParticleSmear;
using namespace MELAStreamHelpers;


LHEConverter::LHEConverter(OptionParser* options_) : converter(options_){
  configure();
  run();
}

void LHEConverter::configure(){}
void LHEConverter::finalizeRun(){}
void LHEConverter::run(){
  constexpr bool useNewJetMerging = true;

  double weight;
  Float_t MC_weight=0;
  Int_t isSelected=0;

  int globalNEvents = 0;
  int maxProcEvents = options->maxEventsToProcess();
  vector<pair<Int_t, Int_t>> eventSkipList = options->getSkippedEvents();

  tree->bookAllBranches(false);

  for (string const& cinput:filename){
    MELAout << "Processing " << cinput << "..." << endl;
    ifstream fin;
    fin.open(cinput.c_str());
    if (fin.good()){
      int nProcessed = 0;
      int ev = 0;
      int fileLine = 0; // Debugging
      while (!fin.eof()){
        // Bookkeeping vectors
        vector<MELAParticle*> particleList = readEvent(fin, fileLine, weight);
        vector<MELAParticle*> smearedParticleList;
        vector<MELACandidate*> candList;
        vector<MELACandidate*> smearedCandList;

        if (globalNEvents>=maxProcEvents && maxProcEvents>=0) break;
        if (particleList.empty()) continue;
        bool doSkipEvent = false;
        for (auto const& es:eventSkipList) doSkipEvent |= (
          (es.first<=globalNEvents && es.second>=globalNEvents)
          ||
          (es.first<=globalNEvents && es.second<0)
          );
        if (doSkipEvent){ globalNEvents++; ev++; continue; }

        if (particleList.empty() && weight!=0.) weight=0;
        if (weight!=0.){
          tree->initializeBranches();

          MELAEvent genEvent;
          genEvent.setWeight(weight);
          vector<MELAParticle*> writtenGenCands;
          vector<MELAParticle*> writtenGenTopCands;
          vector<MELAParticle const*> injetparts;

          MELAEvent smearedEvent;
          smearedEvent.setWeight(weight);

          for (MELAParticle* genPart:particleList){
            if (isATopQuark(genPart->id)){
              writtenGenTopCands.push_back(genPart);
              if (genPart->genStatus==1) genEvent.addIntermediate(genPart);
            }
            if (isAHiggs(genPart->id)){
              writtenGenCands.push_back(genPart);
              if (options->doGenHZZdecay()==MELAEvent::UndecayedMode && (genPart->genStatus==1 || genPart->genStatus==2)) genEvent.addIntermediate(genPart);
            }
            if (genPart->genStatus==1){
              bool addJetConstituent = false;
              bool skipSmearingGen = false;
              if (isALepton(genPart->id)){ genEvent.addLepton(genPart); addJetConstituent=useNewJetMerging; }
              else if (isANeutrino(genPart->id)) genEvent.addNeutrino(genPart);
              else if (isAPhoton(genPart->id)){ genEvent.addPhoton(genPart); addJetConstituent=useNewJetMerging; }
              else if (isAKnownJet(genPart->id) && !isATopQuark(genPart->id)){ genEvent.addJet(genPart); addJetConstituent=useNewJetMerging; skipSmearingGen=useNewJetMerging; }

              MELAParticle* smearedPart=nullptr; // Has no mother info
              if (!skipSmearingGen && options->recoSmearingMode()==0) smearedPart = smearParticle(genPart);
              else smearedPart = new MELAParticle(*genPart);
              smearedParticleList.push_back(smearedPart);
              if (isALepton(smearedPart->id)) smearedEvent.addLepton(smearedPart);
              else if (isANeutrino(smearedPart->id)) smearedEvent.addNeutrino(smearedPart);
              else if (isAPhoton(smearedPart->id)) smearedEvent.addPhoton(smearedPart);
              else if (isAKnownJet(smearedPart->id) && !addJetConstituent){
                if (!isATopQuark(smearedPart->id)) smearedEvent.addJet(smearedPart);
                else smearedEvent.addIntermediate(smearedPart);
                smearedPart->id=0; // Wipe id from reco. quark/gluon
              }
              else smearedEvent.addParticle(smearedPart);

              if (addJetConstituent) injetparts.push_back(smearedPart);
            }
            else if (genPart->genStatus==-1){
              genEvent.addMother(genPart);
              tree->fillMotherInfo(genPart);
            }
          }

          genEvent.constructTopCandidates();
          // Disable tops unmatched to a gen. top
          {
            vector<MELATopCandidate_t*> matchedTops;
            for (auto* writtenGenTopCand:writtenGenTopCands){
              MELATopCandidate_t* tmpCand = TopComparators::matchATopToParticle(genEvent, writtenGenTopCand);
              if (tmpCand) matchedTops.push_back(tmpCand);
            }
            for (MELATopCandidate_t* tmpCand:genEvent.getTopCandidates()){
              if (std::find(matchedTops.begin(), matchedTops.end(), tmpCand)==matchedTops.end()) tmpCand->setSelected(false);
            }
          }

          if (debugVars::debugFlag) MELAout << "Starting to construct gen. VV candidates." << endl;
          genEvent.constructVVCandidates(options->doGenHZZdecay(), options->genDecayProducts());
          if (debugVars::debugFlag) MELAout << "Successfully constructed gen. VV candidates." << endl;
          if (debugVars::debugFlag) MELAout << "Starting to add gen. VV candidate appendages." << endl;
          genEvent.addVVCandidateAppendages();
          MELACandidate* genCand=nullptr;
          if (debugVars::debugFlag) MELAout << "Number of gen. Higgs candidates directly from the LHE: " << writtenGenCands.size() << endl;
          if (!writtenGenCands.empty()){
            for (auto* writtenGenCand:writtenGenCands){
              MELACandidate* tmpCand = HiggsComparators::matchAHiggsToParticle(genEvent, writtenGenCand);
              if (tmpCand){
                if (!genCand) genCand = tmpCand;
                else genCand = HiggsComparators::candComparator(genCand, tmpCand, options->getHiggsCandidateSelectionScheme(true), options->doGenHZZdecay());
              }
            }
          }
          if (!genCand) genCand = HiggsComparators::candidateSelector(genEvent, options->getHiggsCandidateSelectionScheme(true), options->doGenHZZdecay());
          if (genCand) tree->fillCandidate(genCand, true);
          else MELAout << cinput << " (" << ev << "): No gen. level Higgs candidate was found!" << endl;

          // Get merged jets and add htem to the list of particles
          if (useNewJetMerging){ // Not really necessary to put the if-statement, but better to allow compiler optimization
            vector<MELAParticle> outjets;
            MELAJetMerger::mergeParticlesToJets(options->jetAlgorithm(), injetparts, outjets);
            for (MELAParticle const& part:outjets){
              MELAParticle* smearedPart=nullptr;
              if (options->recoSmearingMode()==0 && isAJet(part.id)) smearedPart = smearParticle(&part);
              else smearedPart = new MELAParticle(part);
              smearedPart->id=0;
              smearedEvent.addJet(smearedPart);
              smearedParticleList.push_back(smearedPart);
            }
          }
          smearedEvent.constructTopCandidates();
          smearedEvent.constructVVCandidates(options->doRecoHZZdecay(), options->recoDecayProducts());
          if (options->recoSelectionMode()==0) smearedEvent.applyParticleSelection();
          smearedEvent.addVVCandidateAppendages();

          MELACandidate* rCand = HiggsComparators::candidateSelector(smearedEvent, options->getHiggsCandidateSelectionScheme(false), options->doRecoHZZdecay());
          if (rCand){
            isSelected=1;
            tree->fillCandidate(rCand, false);
          }
          else isSelected=0;
          MC_weight = (float) weight;
          tree->fillEventVariables(MC_weight, isSelected);

          tree->fillXsec(options->get_xsec(), options->get_xsecerr());

          if ((rCand && options->processRecoInfo()) || (genCand && options->processGenInfo())){
            tree->record();
            nProcessed++;
          }
        }
        else MELAerr << "Weight=0 at event " << ev << endl;
        ev++;

        for (auto*& tmpPart:smearedCandList) delete tmpPart;
        for (auto*& tmpPart:smearedParticleList) delete tmpPart;
        for (auto*& tmpPart:candList) delete tmpPart;
        for (auto*& tmpPart:particleList) delete tmpPart;

        globalNEvents++;
        if (globalNEvents % 100000 == 0) MELAout << "Event " << globalNEvents << "..." << endl;
      }
      fin.close();
      MELAout << "Processed number of events from the input file (recorded events / sample size observed / cumulative traversed): " << nProcessed << " / " << ev << " / " << globalNEvents << endl;
    }
  }
  finalizeRun();
}

std::vector<MELAParticle*> LHEConverter::readEvent(std::ifstream& input_lhe, int& fline, double& weight){
  string event_beginning = "<event>";
  string event_end = "</event>";
  string file_closing = "</LesHouchesEvents>";
  string str_in="";

  vector<MELAParticle*> collection;
  vectorInt motherIDs_first;
  vectorInt motherIDs_second;

  // Test whether the string read is the beginning of the event for a valid file
  while (str_in.find(event_beginning)==string::npos){
    if (input_lhe.eof()){
      weight=0;
      return collection;
    }
    getline(input_lhe, str_in); fline++;
    if (str_in.find(file_closing)!=string::npos){
      weight=0;
      return collection;
    }
  }

  int nparticle, para;
  double m_V, alpha_qed, alpha_s;

  fline++;
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
    MELAParticle* onePart = new MELAParticle(idup, partFourVec);
    onePart->setGenStatus(istup);
    onePart->setLifetime(spinup);
    collection.push_back(onePart);
  }

  // Test whether the end of event is reached indeed
  str_in = "";
  while (str_in==""){ getline(input_lhe, str_in); } // Do not count empty lines or e.o.l. in the middle of events
  while (str_in.find("#")!=string::npos){ getline(input_lhe, str_in); fline++; }
  if (str_in.find(event_end)==string::npos){
    MELAerr << "End of event not reached! string is " << str_in << " on line " << fline << endl;
    weight=0;
    for (unsigned int a = 0; a < collection.size(); a++){
      MELAParticle* tmpPart = collection.at(a);
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

