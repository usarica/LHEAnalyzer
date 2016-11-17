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
#include "../interface/LHEConverter.h"

using namespace PDGHelpers;
using namespace LHEParticleSmear;

LHEConverter::LHEConverter(OptionParser* options_) : converter(options_){
  configure();
  run();
}

void LHEConverter::configure(){}
void LHEConverter::finalizeRun(){}
void LHEConverter::run(){
  double weight;
  Float_t MC_weight=0;
  Int_t isSelected=0;

  int globalNEvents = 0;
  int maxProcEvents = options->maxEventsToProcess();
  vector < pair<Int_t, Int_t> > eventSkipList = options->getSkippedEvents();

  tree->bookAllBranches(false);

  for (unsigned int f=0; f<filename.size(); f++){
    string cinput = filename.at(f);
    cout << "Processing " << cinput << "..." << endl;
    ifstream fin;
    fin.open(cinput.c_str());
    if (fin.good()){
      int nProcessed = 0;
      int ev = 0;
      int fileLine = 0; // Debugging
      while (!fin.eof()){
        vector<Particle*> particleList = readEvent(fin, fileLine, weight);
        vector<Particle*> smearedParticleList; // Bookkeeping
        vector<ZZCandidate*> candList; // Bookkeeping
        vector<ZZCandidate*> smearedCandList; // Bookkeeping

        if (globalNEvents>=maxProcEvents && maxProcEvents>=0) break;
        if (particleList.empty()) continue;
        bool doSkipEvent = false;
        for (unsigned int es=0; es<eventSkipList.size(); es++){
          if (
            (eventSkipList.at(es).first<=globalNEvents && eventSkipList.at(es).second>=globalNEvents)
            ||
            (eventSkipList.at(es).first<=globalNEvents && eventSkipList.at(es).second<0)
            )doSkipEvent=true;
        }
        if (doSkipEvent){ globalNEvents++; ev++; continue; }

        if (particleList.size()==0 && weight!=0) weight=0;
        if (weight!=0){
          tree->initializeBranches();

          Event genEvent;
          genEvent.setWeight(weight);
          vectorInt hasGenHiggs;
          Event smearedEvent;
          smearedEvent.setWeight(weight);
          for (unsigned int p=0; p<particleList.size(); p++){
            Particle* genPart = particleList.at(p); // Has mother info from LHE reading
            if (isAHiggs(genPart->id)){
              hasGenHiggs.push_back(p);
              if (options->doGenHZZdecay()==-1 && (genPart->genStatus==1 || genPart->genStatus==2)) genEvent.addIntermediate(genPart);
            }
            if (genPart->genStatus==1){
              if (isALepton(genPart->id)) genEvent.addLepton(genPart);
              else if (isANeutrino(genPart->id)) genEvent.addNeutrino(genPart);
              else if (isAPhoton(genPart->id)) genEvent.addPhoton(genPart);
              else if (isAGluon(genPart->id) || isAQuark(genPart->id)) genEvent.addJet(genPart);

              Particle* smearedPart; // Has no mother info
              if (options->recoSmearingMode()==0) smearedPart = smearParticle(genPart);
              else{
                int recoid = genPart->id;
                TLorentzVector recofv;
                recofv.SetXYZT(genPart->x(), genPart->y(), genPart->z(), genPart->t());
                smearedPart = new Particle(recoid, recofv);
              }
              smearedParticleList.push_back(smearedPart);
              if (isALepton(smearedPart->id)) smearedEvent.addLepton(smearedPart);
              else if (isANeutrino(smearedPart->id)) smearedEvent.addNeutrino(smearedPart);
              else if (isAPhoton(smearedPart->id)) smearedEvent.addPhoton(smearedPart);
              else if (isAGluon(smearedPart->id) || isAQuark(smearedPart->id)){
                smearedPart->id=0; // Wipe id from reco. quark/gluon
                smearedEvent.addJet(smearedPart);
              }
              else smearedEvent.addParticle(smearedPart);
            }
            else if (genPart->genStatus==-1) tree->fillMotherInfo(genPart);
          }

          if (debugVars::debugFlag) cout << "Starting to construct gen. VV candidates." << endl;
          genEvent.constructVVCandidates(options->doGenHZZdecay(), options->genDecayProducts());
          if (debugVars::debugFlag) cout << "Successfully constructed gen. VV candidates." << endl;
          for (unsigned int p=0; p<particleList.size(); p++){
            Particle* genPart = particleList.at(p);
            if (genPart->genStatus==-1) genEvent.addVVCandidateMother(genPart);
          }
          if (debugVars::debugFlag) cout << "Starting to add gen. VV candidate appendages." << endl;
          genEvent.addVVCandidateAppendages();
          ZZCandidate* genCand=0;
          if (debugVars::debugFlag) cout << "Number of gen. Higgs candidates directly from the LHE: " << hasGenHiggs.size() << endl;
          if (hasGenHiggs.size()>0){
            for (unsigned int gk=0; gk<hasGenHiggs.size(); gk++){
              ZZCandidate* tmpCand = HiggsComparators::matchAHiggsToParticle(genEvent, particleList.at(hasGenHiggs.at(gk)));
              if (tmpCand!=0){
                if (genCand==0) genCand=tmpCand;
                else genCand = HiggsComparators::candComparator(genCand, tmpCand, options->getHiggsCandidateSelectionScheme(true), options->doGenHZZdecay());
              }
            }
          }
          if (genCand==0) genCand = HiggsComparators::candidateSelector(genEvent, options->getHiggsCandidateSelectionScheme(true), options->doGenHZZdecay());
          if (genCand!=0) tree->fillCandidate(genCand, true);
          else cout << cinput << " (" << ev << "): No gen. level Higgs candidate was found!" << endl;

          smearedEvent.constructVVCandidates(options->doRecoHZZdecay(), options->recoDecayProducts());
          if (options->recoSelectionMode()==0) smearedEvent.applyParticleSelection();
          smearedEvent.addVVCandidateAppendages();
          ZZCandidate* rCand = HiggsComparators::candidateSelector(smearedEvent, options->getHiggsCandidateSelectionScheme(false), options->doRecoHZZdecay());

          if (rCand!=0){
            isSelected=1;
            tree->fillCandidate(rCand, false);
          }
          else isSelected=0;
          MC_weight = (float)weight;
          tree->fillEventVariables(MC_weight, isSelected);

          if ((rCand!=0 && options->processRecoInfo()) || (genCand!=0 && options->processGenInfo())){
            tree->record();
            nProcessed++;
          }
        }
        else cerr << "Weight=0 at event " << ev << endl;
        ev++;

        for (unsigned int p=0; p<smearedCandList.size(); p++){ // Bookkeeping
          ZZCandidate* tmpCand = (ZZCandidate*)smearedCandList.at(p);
          if (tmpCand!=0) delete tmpCand;
        }
        for (unsigned int p=0; p<smearedParticleList.size(); p++){ // Bookkeeping
          Particle* tmpPart = (Particle*)smearedParticleList.at(p);
          if (tmpPart!=0) delete tmpPart;
        }

        for (unsigned int p=0; p<candList.size(); p++){ // Bookkeeping
          ZZCandidate* tmpCand = (ZZCandidate*)candList.at(p);
          if (tmpCand!=0) delete tmpCand;
        }
        for (unsigned int p=0; p<particleList.size(); p++){ // Bookkeeping
          Particle* tmpPart = (Particle*)particleList.at(p);
          if (tmpPart!=0) delete tmpPart;
        }

        // Bookkeeping
        smearedCandList.clear();
        smearedParticleList.clear();
        candList.clear();
        particleList.clear();

        globalNEvents++;
        if (globalNEvents % 100000 == 0) cout << "Event " << globalNEvents << "..." << endl;
      }
      fin.close();
      cout << "Processed number of events from the input file (recorded events / sample size observed / cumulative traversed): " << nProcessed << " / " << ev << " / " << globalNEvents << endl;
    }
  }
  finalizeRun();
}

vector<Particle*> LHEConverter::readEvent(ifstream& input_lhe, int& fline, double& weight){
  string event_beginning = "<event>";
  string event_end = "</event>";
  string file_closing = "</LesHouchesEvents>";
  string str_in="";

  vector<Particle*> collection;
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
    Particle* onePart = new Particle(idup, partFourVec);
    onePart->setGenStatus(istup);
    onePart->setLifetime(spinup);
    collection.push_back(onePart);
  }

// Test whether the end of event is reached indeed
  str_in = "";
  while (str_in==""){ getline(input_lhe, str_in); } // Do not count empty lines or e.o.l. in the middle of events
  while (str_in.find("#")!=string::npos){ getline(input_lhe, str_in); fline++; }
  if (str_in.find(event_end)==string::npos){
    cerr << "End of event not reached! string is " << str_in << " on line " << fline << endl;
    weight=0;
    for (unsigned int a = 0; a < collection.size(); a++){
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

