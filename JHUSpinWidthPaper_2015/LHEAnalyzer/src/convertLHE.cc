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

  Float_t m_costheta1, m_costheta2, m_phi, m_costhetastar, m_phistar1;
  Float_t m_phistar2, m_phistar12, m_phi1, m_phi2;
  Float_t m_zzmass, m_z1mass, m_z2mass, m_zamass, m_zbmass;
  Float_t m_zzpt, m_zzpz, m_z1pt, m_z2pt, m_zapt, m_zbpt;
  Float_t m_zzphi, m_z1phi, m_z2phi, m_zaphi, m_zbphi;
  Float_t m_zzeta, m_z1eta, m_z2eta, m_zaeta, m_zbeta;

  Float_t m_l1minus_mass, m_l1plus_mass, m_l2minus_mass, m_l2plus_mass;
  Float_t m_l1minus_pT, m_l1plus_pT, m_l2minus_pT, m_l2plus_pT;
  Float_t m_l1minus_eta, m_l1plus_eta, m_l2minus_eta, m_l2plus_eta;
  Float_t m_l1minus_phi, m_l1plus_phi, m_l2minus_phi, m_l2plus_phi;
  Int_t m_l1minus_id, m_l1plus_id, m_l2minus_id, m_l2plus_id;

  Float_t m_wt;
  double weight;
  int interf=-1;
  int genFinalState=-1;

  tree->Branch("GenZZMass", &m_zzmass);
  tree->Branch("GenZZPt", &m_zzpt);
  tree->Branch("GenZZPz", &m_zzpz);
  tree->Branch("GenZZPhi", &m_zzphi);

  tree->Branch("GenZ1Mass", &m_z1mass);
  tree->Branch("GenZ1Pt", &m_z1pt);
  tree->Branch("GenZ1Phi", &m_z1phi);
  tree->Branch("GenZ1Eta", &m_z1eta);

  tree->Branch("GenZ2Mass", &m_z2mass);
  tree->Branch("GenZ2Pt", &m_z2pt);
  tree->Branch("GenZ2Phi", &m_z2phi);
  tree->Branch("GenZ2Eta", &m_z2eta);

  tree->Branch("GenZaMass", &m_zamass);
  tree->Branch("GenZaPt", &m_zapt);
  tree->Branch("GenZaPhi", &m_zaphi);
  tree->Branch("GenZaEta", &m_zaeta);

  tree->Branch("GenZbMass", &m_zbmass);
  tree->Branch("GenZbPt", &m_zbpt);
  tree->Branch("GenZbPhi", &m_zbphi);
  tree->Branch("GenZbEta", &m_zbeta);

  tree->Branch("GenhelcosthetaZ1", &m_costheta1);
  tree->Branch("GenhelcosthetaZ2", &m_costheta2);
  tree->Branch("Genhelphi", &m_phi);
  tree->Branch("Gencosthetastar", &m_costhetastar);
  tree->Branch("GenphistarZ1", &m_phistar1);

  tree->Branch("GenLep1Mass", &m_l1minus_mass);
  tree->Branch("GenLep2Mass", &m_l1plus_mass);
  tree->Branch("GenLep3Mass", &m_l2minus_mass);
  tree->Branch("GenLep4Mass", &m_l2plus_mass);
  tree->Branch("GenLep1Pt", &m_l1minus_pT);
  tree->Branch("GenLep2Pt", &m_l1plus_pT);
  tree->Branch("GenLep3Pt", &m_l2minus_pT);
  tree->Branch("GenLep4Pt", &m_l2plus_pT);
  tree->Branch("GenLep1Eta", &m_l1minus_eta);
  tree->Branch("GenLep2Eta", &m_l1plus_eta);
  tree->Branch("GenLep3Eta", &m_l2minus_eta);
  tree->Branch("GenLep4Eta", &m_l2plus_eta);
  tree->Branch("GenLep1Phi", &m_l1minus_phi);
  tree->Branch("GenLep2Phi", &m_l1plus_phi);
  tree->Branch("GenLep3Phi", &m_l2minus_phi);
  tree->Branch("GenLep4Phi", &m_l2plus_phi);
  tree->Branch("GenLep1Id", &m_l1minus_id);
  tree->Branch("GenLep2Id", &m_l1plus_id);
  tree->Branch("GenLep3Id", &m_l2minus_id);
  tree->Branch("GenLep4Id", &m_l2plus_id);
  tree->Branch("genFinalState", &genFinalState);
  tree->Branch("MC_weight", &m_wt, "MC_weight/F");

  for (int f=0; f<filename.size(); f++){
    string cinput = filename.at(f);
    ifstream fin;
    fin.open(cinput.c_str());
    if (fin.good()){
      while (!fin.eof()){
        vector<Particle*> particleList = readLHEEvent(fin, weight);
        if (particleList.size()==0 && weight!=0) weight=0;
        if (weight!=0){

          m_wt = (float)weight;
          tree->Fill();
        }
        for (int p=0; p<particleList.size(); p++){
          Particle* tmpPart = (Particle*) particleList.at(p);
          delete tmpPart;
        }
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





