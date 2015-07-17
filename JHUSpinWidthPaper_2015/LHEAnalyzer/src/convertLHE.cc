#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <vector>
#include <iomanip>
#include <cstdlib>
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
  Float_t GenHPt, GenHPz, GenZ1Pt, GenZ2Pt, GenZaPt, GenZbPt;
  Float_t GenHPhi, GenZ1Phi, GenZ2Phi, GenZaPhi, GenZbPhi;
  Float_t GenHEta, GenZ1Eta, GenZ2Eta, GenZaEta, GenZbEta;

  Float_t GenLep1Mass, GenLep2Mass, GenLep3Mass, GenLep4Mass;
  Float_t GenLep1Pt, GenLep2Pt, GenLep3Pt, GenLep4Pt;
  Float_t GenLep1Eta, GenLep2Eta, GenLep3Eta, GenLep4Eta;
  Float_t GenLep1Phi, GenLep2Phi, GenLep3Phi, GenLep4Phi;
  Int_t GenLep1Id, GenLep2Id, GenLep3Id, GenLep4Id;

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

  Float_t MC_weight;
  double weight;
  int interf=-1;
  int genFinalState=-1;

  tree->Branch("GenZZMass", &GenHMass);
  tree->Branch("GenZZPt", &GenHPt);
  tree->Branch("GenZZPz", &GenHPz);
  tree->Branch("GenZZPhi", &GenHPhi);

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
  tree->Branch("MC_weight", &MC_weight, "MC_weight/F");

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
        vector<Particle*> particleList = readLHEEvent(fin, weight);
        vector<Particle*> smearedParticleList; // Bookkeeping
        vector<ZZCandidate*> candList; // Bookkeeping
        vector<ZZCandidate*> smearedCandList; // Bookkeeping

        if (particleList.size()==0 && weight!=0) weight=0;
        if (weight!=0){
          Event genEvent;
          genEvent.setWeight(weight);
          Event smearedEvent;
          smearedEvent.setWeight(weight);
          for (int p=0; p<particleList.size(); p++){
            Particle* genPart = particleList.at(p); // Has mother info from LHE reading
            // No ZZCandidate at this moment
            if (isALepton(genPart->id)) genEvent.addLepton(genPart);
            else if (isANeutrino(genPart->id)) genEvent.addNeutrino(genPart);
            else if (isAGluon(genPart->id) || isAQuark(genPart->id)) genEvent.addJet(genPart);
            else genEvent.addParticle(genPart);

            if (genPart->genStatus==1){
              Particle* smearedPart = smearParticle(genPart); // Has no mother info
              smearedParticleList.push_back(smearedPart);
              if (isALepton(smearedPart->id)) smearedEvent.addLepton(smearedPart);
              else if (isANeutrino(smearedPart->id)) smearedEvent.addNeutrino(smearedPart);
              else if (isAGluon(smearedPart->id) || isAQuark(smearedPart->id)) smearedEvent.addJet(smearedPart);
              else smearedEvent.addParticle(smearedPart);
            }
          }

          //for (int p=0; p<genEvent.getNLeptons(); p++) cout << "Lepton " << p << " (x, y, z, t): " << genEvent.getLepton(p)->x() << '\t' << genEvent.getLepton(p)->y() << '\t' << genEvent.getLepton(p)->z() << '\t' << genEvent.getLepton(p)->t() << endl;

          genEvent.constructVVCandidates();
          ZZCandidate* genCand=0;
          for (int t=0; t<genEvent.getNZZCandidates(); t++){
            ZZCandidate* tmpCand = genEvent.getZZCandidate(t);
            if (genCand==0) genCand=tmpCand;
            else if (fabs(genCand->getSortedV(0)->m()-PDGHelpers::HVVmass)>fabs(tmpCand->getSortedV(0)->m()-PDGHelpers::HVVmass)) genCand=tmpCand;
          }
          Particle* gZ1=genCand->getSortedV(0);
          Particle* gZ2=genCand->getSortedV(1);

          GenHMass=genCand->m();
          GenZ1Mass=gZ1->m();
          GenZ2Mass=gZ2->m();

          smearedEvent.constructVVCandidates();
          ZZCandidate* rCand=0;
          for (int t=0; t<smearedEvent.getNZZCandidates(); t++){
            ZZCandidate* tmpCand = smearedEvent.getZZCandidate(t);
            if (rCand==0) rCand=tmpCand;
            else if (fabs(rCand->getSortedV(0)->m()-PDGHelpers::HVVmass)>fabs(tmpCand->getSortedV(0)->m()-PDGHelpers::HVVmass)) rCand=tmpCand;
          }
          Particle* rZ1=rCand->getSortedV(0);
          Particle* rZ2=rCand->getSortedV(1);

          ZZMass=rCand->m();
          Z1Mass=rZ1->m();
          Z2Mass=rZ2->m();

          MC_weight = (float)weight;
          tree->Fill();
        }

        for (int p=0; p<smearedCandList.size(); p++){ // Bookkeeping
          ZZCandidate* tmpCand = (ZZCandidate*)smearedCandList.at(p);
          if (tmpCand!=0) delete tmpCand;
        }
        for (int p=0; p<smearedParticleList.size(); p++){ // Bookkeeping
          Particle* tmpPart = (Particle*)smearedParticleList.at(p);
//          cout << "smeared "  << tmpPart->genStatus << '\t'  << tmpPart->id << '\t' << tmpPart->m() << endl;
          if (tmpPart!=0) delete tmpPart;
        }

        for (int p=0; p<candList.size(); p++){ // Bookkeeping
          ZZCandidate* tmpCand = (ZZCandidate*)candList.at(p);
          if (tmpCand!=0) delete tmpCand;
        }
        for (int p=0; p<particleList.size(); p++){ // Bookkeeping
          Particle* tmpPart = (Particle*)particleList.at(p);
//          cout << "gen "  << tmpPart->genStatus << '\t'  << tmpPart->id << '\t' << tmpPart->m() << endl;
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


void calculateAngles(TLorentzVector thep4H, TLorentzVector thep4Z1, TLorentzVector thep4M11, TLorentzVector thep4M12, TLorentzVector thep4Z2, TLorentzVector thep4M21, TLorentzVector thep4M22, double& costheta1, double& costheta2, double& phi, double& costhetastar, double& phistar1, double& phistar2, double& phistar12, double& phi1, double& phi2){
	
  float norm;
  
  TVector3 boostX = -(thep4H.BoostVector());
  TLorentzVector thep4Z1inXFrame( thep4Z1 );
  TLorentzVector thep4Z2inXFrame( thep4Z2 );	
  thep4Z1inXFrame.Boost( boostX );
  thep4Z2inXFrame.Boost( boostX );
  TVector3 theZ1X_p3 = TVector3( thep4Z1inXFrame.X(), thep4Z1inXFrame.Y(), thep4Z1inXFrame.Z() );
  TVector3 theZ2X_p3 = TVector3( thep4Z2inXFrame.X(), thep4Z2inXFrame.Y(), thep4Z2inXFrame.Z() );
  
  // calculate phi1, phi2, costhetastar
  phi1 = theZ1X_p3.Phi();
  phi2 = theZ2X_p3.Phi();
  
  ///////////////////////////////////////////////
  // check for z1/z2 convention, redefine all 4 vectors with convention
  ///////////////////////////////////////////////	
  TLorentzVector p4H, p4Z1, p4M11, p4M12, p4Z2, p4M21, p4M22;
  p4H = thep4H;
  
  /* ORDER OF Z1 AND Z2 ALREADY CHOSEN IN MAIN FUNCTION!!!!!! - - - - - - 
     if ((phi1 < 0)&&(phi1 >= -TMath::Pi())){   // old convention based on phi
     p4Z1 = thep4Z2; p4M11 = thep4M21; p4M12 = thep4M22;
     p4Z2 = thep4Z1; p4M21 = thep4M11; p4M22 = thep4M12;		
     costhetastar = theZ2X_p3.CosTheta();
     }
     else{
     p4Z1 = thep4Z1; p4M11 = thep4M11; p4M12 = thep4M12;
     p4Z2 = thep4Z2; p4M21 = thep4M21; p4M22 = thep4M22;
     costhetastar = theZ1X_p3.CosTheta();
     }
     - - - - - - - - - - - - - - - - - - - - - - - - -  - - - - - - - - -*/
  
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
  
  // and the calculate phistar2
  TLorentzVector n_p4Z2inXFrame( p4Z2 );
  n_p4Z2inXFrame.Boost( boostX );
  TVector3 n_p4Z2inXFrame_unit = n_p4Z2inXFrame.Vect().Unit();
  TVector3 n_unitz_2( n_p4Z2inXFrame_unit );
  //// y-axis is defined by neg lepton cross z-axis
  //// the subtle part is here...
  TVector3 n_unity_2 = n_unitz_2.Cross( n_p4M21inXFrame_unit );
  TVector3 n_unitx_2 = n_unity_2.Cross( n_unitz_2 );
  TVector3 n_p4PartoninZ2PlaneFrame_unitprime( n_p4PartoninXFrame_unit.Dot(n_unitx_2), n_p4PartoninXFrame_unit.Dot(n_unity_2), n_p4PartoninXFrame_unit.Dot(n_unitz_2) );
  phistar2 = (n_p4PartoninZ2PlaneFrame_unitprime.Phi());
  
  float phistar12_0 = phistar1 + phistar2;
  if (phistar12_0 > TMath::Pi()) phistar12 = phistar12_0 - 2*TMath::Pi();
  else if (phistar12_0 < (-1.)*TMath::Pi()) phistar12 = phistar12_0 + 2*TMath::Pi();
  else phistar12 = phistar12_0;
	
}





