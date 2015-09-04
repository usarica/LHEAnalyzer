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
#include "../interface/convertLHE.h"

using namespace PDGHelpers;
using namespace LHEParticleSmear;

convertLHE::convertLHE(OptionParser* options_) : converter(options_){
  configure();
  run();
}

void convertLHE::configure(){}
void convertLHE::finalizeRun(){}
void convertLHE::run(){
  double weight;
  Float_t MC_weight=0;
  Int_t isSelected=0;

  tree->bookAllBranches(false);

  for (int f=0; f<filename.size(); f++){
    string cinput = filename.at(f);
    cout << "Processing " << cinput << "..." << endl;
    ifstream fin;
    fin.open(cinput.c_str());
    if (fin.good()){
      int nProcessed = 0;
      int ev = 0;
      while (!fin.eof()){
        vector<Particle*> particleList = readEvent(fin, weight);
        vector<Particle*> smearedParticleList; // Bookkeeping
        vector<ZZCandidate*> candList; // Bookkeeping
        vector<ZZCandidate*> smearedCandList; // Bookkeeping

        if (particleList.size()==0 && weight!=0) weight=0;
        if (weight!=0){
          tree->initializeBranches();

          Event genEvent;
          genEvent.setWeight(weight);
          vectorInt hasGenHiggs;
          Event smearedEvent;
          smearedEvent.setWeight(weight);
          for (int p=0; p<particleList.size(); p++){
            Particle* genPart = particleList.at(p); // Has mother info from LHE reading
            if (isAHiggs(genPart->id)) hasGenHiggs.push_back(p);

            if (genPart->genStatus==1){
              if (isALepton(genPart->id)) genEvent.addLepton(genPart);
              else if (isANeutrino(genPart->id)) genEvent.addNeutrino(genPart);
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
              else if (isAGluon(smearedPart->id) || isAQuark(smearedPart->id)){
                smearedPart->id=0; // Wipe id from reco. quark/gluon
                smearedEvent.addJet(smearedPart);
              }
              else smearedEvent.addParticle(smearedPart);
            }
            else if (genPart->genStatus==-1) tree->fillMotherInfo(genPart);
          }

          genEvent.constructVVCandidates(options->doGenHZZdecay(), options->genDecayProducts());
          for (int p=0; p<particleList.size(); p++){
            Particle* genPart = particleList.at(p);
            if (genPart->genStatus==-1) genEvent.addVVCandidateMother(genPart);
          }
          genEvent.addVVCandidateAppendages();
          ZZCandidate* genCand=0;
          if (hasGenHiggs.size()>0){
            for (int gk=0; gk<hasGenHiggs.size(); gk++){
              ZZCandidate* tmpCand = HiggsComparators::matchAHiggsToParticle(genEvent, particleList.at(hasGenHiggs.at(gk)));
              if (tmpCand!=0){
                if (genCand==0) genCand=tmpCand;
                else genCand = HiggsComparators::candComparator(genCand, tmpCand, options->getHiggsCandidateSelectionScheme(true), options->doGenHZZdecay());
              }
            }
          }
          else genCand = HiggsComparators::candidateSelector(genEvent, options->getHiggsCandidateSelectionScheme(true), options->doGenHZZdecay());
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
      cout << "Processed number of events from the input file: " << nProcessed << " / " << ev << endl;
    }
  }
  finalizeRun();
}

vector<Particle*> convertLHE::readEvent(ifstream& input_lhe, double& weight){
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

