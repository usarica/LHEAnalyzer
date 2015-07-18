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
#include "../interface/HVVTree.h"

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
  HVVTree* tree = new HVVTree(TREE_NAME, TREE_NAME);

  double weight;
  Float_t MC_weight;
  Int_t isSelected;
  Int_t genFinalState=-1;

  tree->bookAllBranches();

  for (int f=0; f<filename.size(); f++){
    string cinput = filename.at(f);
    ifstream fin;
    fin.open(cinput.c_str());
    if (fin.good()){
      while (!fin.eof()){
        vector<Particle*> particleList = readLHEEvent(fin, weight);
        vector<Particle*> smearedParticleList; // Bookkeeping
        vector<ZZCandidate*> candList; // Bookkeeping
        vector<ZZCandidate*> smearedCandList; // Bookkeeping

        if (particleList.size()==0 && weight!=0) weight=0;
        if (weight!=0){
          tree->initializeBranches();

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
            else if (genPart->genStatus==-1) tree->fillMotherInfo(genPart);
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
          if (genCand!=0) tree->fillCandidate(genCand, true);
          else{
            cout << "No gen. level Higgs candidate was found!" << endl;
//            tree->fillCandidate(genCand, true);
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
            tree->fillCandidate(rCand, false);
          }
          else{
            isSelected=0;
//            tree->fillCandidate(rCand, false);
          }
          MC_weight = (float)weight;

          tree->fillEventVariables(MC_weight, isSelected, genFinalState);

          tree->record();
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

  cout << "Number of recorded events: " << tree->getTree()->GetEntries() << endl;
  tree->writeTree(foutput);
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

