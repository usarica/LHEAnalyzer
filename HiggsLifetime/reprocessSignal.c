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
#include "TRandom3.h"
#include "TLorentzVector.h"
#include <ZZMatrixElement/MELA/interface/Mela.h>
#include <ZZMatrixElement/MELA/src/computeAngles.h>
#include "./data/ZZ4l_125_6_Samples.h"

using namespace std;

float getMCFMMELAWeight(Mela& myMela, int lepId[4], float angularOrdered[8], double ggcoupl[2],int useConstant=0){
	float myprob=1.0;
	int myflavor=-1;
	if(abs(lepId[0])==abs(lepId[1]) &&
		abs(lepId[0])==abs(lepId[2]) &&
		abs(lepId[0])==abs(lepId[3])){
			if(abs(lepId[0])==11) myflavor=1;
			else myflavor=2;
	}
	else myflavor=3;

	if(myflavor>=0) myMela.computeP(angularOrdered[0],angularOrdered[1],angularOrdered[2],angularOrdered[3],
		angularOrdered[4],angularOrdered[5],angularOrdered[6],angularOrdered[7],
	    myflavor,
	    myprob,
		useConstant
		);
	return myprob;
}
float getJHUGenMELAWeight(Mela& myMela, int lepId[4], float angularOrdered[8], double selfDHvvcoupl[SIZE_HVV][2]){
	float myprob=1.0;
	int myflavor=-1;
	if(abs(lepId[0])==abs(lepId[1]) &&
		abs(lepId[0])==abs(lepId[2]) &&
		abs(lepId[0])==abs(lepId[3])){
			if(abs(lepId[0])==11) myflavor=1;
			else myflavor=2;
	}
	else myflavor=3;

	if(myflavor>=0) myMela.computeP(angularOrdered[0],angularOrdered[1],angularOrdered[2],angularOrdered[3],
		angularOrdered[4],angularOrdered[5],angularOrdered[6],angularOrdered[7],
	    myflavor,
	    selfDHvvcoupl,
	    myprob
		);
	return myprob;
}
float getSuperMELA(Mela& myMela, int lepId[4], float mZZ, TVar::SuperMelaSyst syst){
	float myprob=1.0;
	TVar::LeptonFlavor myflavor=TVar::Flavor_Dummy;
	if(abs(lepId[0])==abs(lepId[1]) &&
		abs(lepId[0])==abs(lepId[2]) &&
		abs(lepId[0])==abs(lepId[3])){
			if(abs(lepId[0])==11) myflavor=TVar::Flavor_4e;
			else myflavor=TVar::Flavor_4mu;
	}
	else myflavor=TVar::Flavor_2e2mu;

	if(myflavor!=TVar::Flavor_Dummy) myMela.computePM4l(mZZ,
	    myflavor,
	    syst,
	    myprob
		);
	return myprob;
}

void reprocessSignal(int erg_tev, int smp_min = 0, int smp_max = 3){
	int maxEvents = 20000000;
	bool debug = false;
	string input_main_folder;
	int EnergyIndex=1;
	if (erg_tev == 7){ input_main_folder = "LHC_7TeV"; EnergyIndex = 0; }
	else if (erg_tev == 8){ input_main_folder = "LHC_8TeV"; EnergyIndex = 1; }
	string user_main = user_dir_hep + input_main_folder;


	float mPOLE=125.6;
	float wPOLE=4.15e-3;
	Mela mela(erg_tev,mPOLE);


	Float_t m_costheta1, m_costheta2, m_phi, m_costhetastar, m_phistar1;
	Float_t m_zzmass, m_z1mass, m_z2mass;
	Float_t m_l1minus_pT, m_l1plus_pT, m_l2minus_pT, m_l2plus_pT;
	Float_t m_l1minus_eta, m_l1plus_eta, m_l2minus_eta, m_l2plus_eta;
	Float_t m_l1minus_phi, m_l1plus_phi, m_l2minus_phi, m_l2plus_phi;
	Int_t m_l1minus_id, m_l1plus_id, m_l2minus_id, m_l2plus_id;
	Float_t Y4l, eta4l, phi4l, pT4l, m_wt, interf;
	Float_t CandVtx_x,CandVtx_y,CandVtx_z;
	Int_t isSelected=0,passVtx=0;
	Int_t passLeptIso=0,passLepPT=0,passLepEta=0,passZLepSel=0,passZMassSel=0;
	Float_t m_gen_costheta1, m_gen_costheta2, m_gen_phi, m_gen_costhetastar, m_gen_phistar1;
	Float_t m_gen_zzmass, m_gen_z1mass, m_gen_z2mass;
	Float_t m_gen_l1minus_pT, m_gen_l1plus_pT, m_gen_l2minus_pT, m_gen_l2plus_pT;
	Float_t m_gen_l1minus_eta, m_gen_l1plus_eta, m_gen_l2minus_eta, m_gen_l2plus_eta;
	Float_t m_gen_l1minus_phi, m_gen_l1plus_phi, m_gen_l2minus_phi, m_gen_l2plus_phi;
	Float_t GenIntVtx_x,GenIntVtx_y,GenIntVtx_z;
	Int_t m_gen_l1minus_id, m_gen_l1plus_id, m_gen_l2minus_id, m_gen_l2plus_id;
	Float_t gen_eta4l, gen_phi4l, gen_pT4l, gen_ctau;

	float p0plus_VAJHU;
	float p0plus_VAMCFM;
	float bkg_VAMCFM;
	float p0plus_m4l, p0plus_m4l_ScaleUp, p0plus_m4l_ScaleDown, p0plus_m4l_ResUp, p0plus_m4l_ResDown;
	float bkg_m4l, bkg_m4l_ScaleUp, bkg_m4l_ScaleDown, bkg_m4l_ResUp, bkg_m4l_ResDown;

	double PeakYield[2][3] = {
 		{1.0902689,0.6087736,1.4452079},
 		{5.1963998,2.6944639,6.6562629}
	};

	for (int s = smp_min; s < smp_max; s++){
		string cinput_main = user_main + "/" + sample_suffix[s] + "/";
		char cinput[1000];
		sprintf(cinput, "%s%s%s%s", cinput_main.c_str(), "/HZZ4lTree_powheg15jhuGenV3-CTau", sample_suffix[s], "-0PMH125.6_GenLevel.root");
		char coutput[1000];
		sprintf(coutput, "%s%s%s%s", cinput_main.c_str(), "/HZZ4lTree_powheg15jhuGenV3-CTau", sample_suffix[s], "-0PMH125.6_GenLevel_Reprocessed.root");

		TFile* finput = new TFile(cinput,"read");
		TTree* tinput = (TTree*) finput->Get("SelectedTree");

		tinput->SetBranchAddress("ZZMass", &m_zzmass);
		tinput->SetBranchAddress("Z1Mass", &m_z1mass);
		tinput->SetBranchAddress("Z2Mass", &m_z2mass);
		tinput->SetBranchAddress("helcosthetaZ1", &m_costheta1);
		tinput->SetBranchAddress("helcosthetaZ2", &m_costheta2);
		tinput->SetBranchAddress("helphi", &m_phi);
		tinput->SetBranchAddress("costhetastar", &m_costhetastar);
		tinput->SetBranchAddress("phistarZ1", &m_phistar1);
		tinput->SetBranchAddress("Lep1Pt", &m_l1minus_pT);
		tinput->SetBranchAddress("Lep2Pt", &m_l1plus_pT);
		tinput->SetBranchAddress("Lep3Pt", &m_l2minus_pT);
		tinput->SetBranchAddress("Lep4Pt", &m_l2plus_pT);
		tinput->SetBranchAddress("Lep1Eta", &m_l1minus_eta);
		tinput->SetBranchAddress("Lep2Eta", &m_l1plus_eta);
		tinput->SetBranchAddress("Lep3Eta", &m_l2minus_eta);
		tinput->SetBranchAddress("Lep4Eta", &m_l2plus_eta);
		tinput->SetBranchAddress("Lep1Phi", &m_l1minus_phi);
		tinput->SetBranchAddress("Lep2Phi", &m_l1plus_phi);
		tinput->SetBranchAddress("Lep3Phi", &m_l2minus_phi);
		tinput->SetBranchAddress("Lep4Phi", &m_l2plus_phi);
		tinput->SetBranchAddress("Lep1ID", &m_l1minus_id);
		tinput->SetBranchAddress("Lep2ID", &m_l1plus_id);
		tinput->SetBranchAddress("Lep3ID", &m_l2minus_id);
		tinput->SetBranchAddress("Lep4ID", &m_l2plus_id);
		tinput->SetBranchAddress("ZZPt", &pT4l);
		tinput->SetBranchAddress("ZZEta", &eta4l);
		tinput->SetBranchAddress("ZZPhi", &phi4l);
		tinput->SetBranchAddress("CandVtx_x", &CandVtx_x);
		tinput->SetBranchAddress("CandVtx_y", &CandVtx_y);
		tinput->SetBranchAddress("CandVtx_z", &CandVtx_z);
		tinput->SetBranchAddress("isSelected", &isSelected);
		tinput->SetBranchAddress("passLeptIso", &passLeptIso);
		tinput->SetBranchAddress("passLepPT", &passLepPT);
		tinput->SetBranchAddress("passLepEta", &passLepEta);
		tinput->SetBranchAddress("passZLepSel", &passZLepSel);
		tinput->SetBranchAddress("passZMassSel", &passZMassSel);
		tinput->SetBranchAddress("passVtx", &passVtx);

		tinput->SetBranchAddress("GenHMass", &m_gen_zzmass);
		tinput->SetBranchAddress("GenZ1Mass", &m_gen_z1mass);
		tinput->SetBranchAddress("GenZ2Mass", &m_gen_z2mass);
		tinput->SetBranchAddress("GenhelcosthetaZ1", &m_gen_costheta1);
		tinput->SetBranchAddress("GenhelcosthetaZ2", &m_gen_costheta2);
		tinput->SetBranchAddress("Genhelphi", &m_gen_phi);
		tinput->SetBranchAddress("Gencosthetastar", &m_gen_costhetastar);
		tinput->SetBranchAddress("GenphistarZ1", &m_gen_phistar1);
		tinput->SetBranchAddress("GenLep1Pt", &m_gen_l1minus_pT);
		tinput->SetBranchAddress("GenLep2Pt", &m_gen_l1plus_pT);
		tinput->SetBranchAddress("GenLep3Pt", &m_gen_l2minus_pT);
		tinput->SetBranchAddress("GenLep4Pt", &m_gen_l2plus_pT);
		tinput->SetBranchAddress("GenLep1Eta", &m_gen_l1minus_eta);
		tinput->SetBranchAddress("GenLep2Eta", &m_gen_l1plus_eta);
		tinput->SetBranchAddress("GenLep3Eta", &m_gen_l2minus_eta);
		tinput->SetBranchAddress("GenLep4Eta", &m_gen_l2plus_eta);
		tinput->SetBranchAddress("GenLep1Phi", &m_gen_l1minus_phi);
		tinput->SetBranchAddress("GenLep2Phi", &m_gen_l1plus_phi);
		tinput->SetBranchAddress("GenLep3Phi", &m_gen_l2minus_phi);
		tinput->SetBranchAddress("GenLep4Phi", &m_gen_l2plus_phi);
		tinput->SetBranchAddress("GenLep1ID", &m_l1minus_id);
		tinput->SetBranchAddress("GenLep2ID", &m_l1plus_id);
		tinput->SetBranchAddress("GenLep3ID", &m_l2minus_id);
		tinput->SetBranchAddress("GenLep4ID", &m_l2plus_id);
		tinput->SetBranchAddress("GenHPt", &gen_pT4l);
		tinput->SetBranchAddress("GenHEta", &gen_eta4l);
		tinput->SetBranchAddress("GenHPhi", &gen_phi4l);
		tinput->SetBranchAddress("GenIntVtx_x", &GenIntVtx_x);
		tinput->SetBranchAddress("GenIntVtx_y", &GenIntVtx_y);
		tinput->SetBranchAddress("GenIntVtx_z", &GenIntVtx_z);

		tinput->SetBranchAddress("GenCTau", &gen_ctau); // Will be in microns

		TFile* foutput = new TFile(coutput,"recreate");
		TTree* tree = new TTree("SelectedTree", "SelectedTree");

		tree->Branch("ZZMass", &m_zzmass);
		tree->Branch("Z1Mass", &m_z1mass);
		tree->Branch("Z2Mass", &m_z2mass);
		tree->Branch("helcosthetaZ1", &m_costheta1);
		tree->Branch("helcosthetaZ2", &m_costheta2);
		tree->Branch("helphi", &m_phi);
		tree->Branch("costhetastar", &m_costhetastar);
		tree->Branch("phistarZ1", &m_phistar1);
		tree->Branch("Lep1Pt", &m_l1minus_pT);
		tree->Branch("Lep2Pt", &m_l1plus_pT);
		tree->Branch("Lep3Pt", &m_l2minus_pT);
		tree->Branch("Lep4Pt", &m_l2plus_pT);
		tree->Branch("Lep1Eta", &m_l1minus_eta);
		tree->Branch("Lep2Eta", &m_l1plus_eta);
		tree->Branch("Lep3Eta", &m_l2minus_eta);
		tree->Branch("Lep4Eta", &m_l2plus_eta);
		tree->Branch("Lep1Phi", &m_l1minus_phi);
		tree->Branch("Lep2Phi", &m_l1plus_phi);
		tree->Branch("Lep3Phi", &m_l2minus_phi);
		tree->Branch("Lep4Phi", &m_l2plus_phi);
		tree->Branch("Lep1ID", &m_l1minus_id);
		tree->Branch("Lep2ID", &m_l1plus_id);
		tree->Branch("Lep3ID", &m_l2minus_id);
		tree->Branch("Lep4ID", &m_l2plus_id);
		tree->Branch("ZZPt", &pT4l);
		tree->Branch("ZZEta", &eta4l);
		tree->Branch("ZZPhi", &phi4l);
		tree->Branch("CandVtx_x", &CandVtx_x);
		tree->Branch("CandVtx_y", &CandVtx_y);
		tree->Branch("CandVtx_z", &CandVtx_z);
		tree->Branch("MC_weight", &m_wt);

		tree->Branch("p0plus_VAJHU", &p0plus_VAJHU);
		tree->Branch("p0plus_VAMCFM", &p0plus_VAMCFM);
		tree->Branch("bkg_VAMCFM", &bkg_VAMCFM);
		tree->Branch("p0plus_m4l", &p0plus_m4l);
		tree->Branch("p0plus_m4l_ScaleUp", &p0plus_m4l_ScaleUp);
		tree->Branch("p0plus_m4l_ScaleDown", &p0plus_m4l_ScaleDown);
		tree->Branch("p0plus_m4l_ResUp", &p0plus_m4l_ResUp);
		tree->Branch("p0plus_m4l_ResDown", &p0plus_m4l_ResDown);
		tree->Branch("bkg_m4l", &bkg_m4l);
		tree->Branch("bkg_m4l_ScaleUp", &bkg_m4l_ScaleUp);
		tree->Branch("bkg_m4l_ScaleDown", &bkg_m4l_ScaleDown);
		tree->Branch("bkg_m4l_ResUp", &bkg_m4l_ResUp);
		tree->Branch("bkg_m4l_ResDown", &bkg_m4l_ResDown);

		tree->Branch("GenHMass", &m_gen_zzmass);
		tree->Branch("GenZ1Mass", &m_gen_z1mass);
		tree->Branch("GenZ2Mass", &m_gen_z2mass);
		tree->Branch("GenhelcosthetaZ1", &m_gen_costheta1);
		tree->Branch("GenhelcosthetaZ2", &m_gen_costheta2);
		tree->Branch("Genhelphi", &m_gen_phi);
		tree->Branch("Gencosthetastar", &m_gen_costhetastar);
		tree->Branch("GenphistarZ1", &m_gen_phistar1);
		tree->Branch("GenLep1Pt", &m_gen_l1minus_pT);
		tree->Branch("GenLep2Pt", &m_gen_l1plus_pT);
		tree->Branch("GenLep3Pt", &m_gen_l2minus_pT);
		tree->Branch("GenLep4Pt", &m_gen_l2plus_pT);
		tree->Branch("GenLep1Eta", &m_gen_l1minus_eta);
		tree->Branch("GenLep2Eta", &m_gen_l1plus_eta);
		tree->Branch("GenLep3Eta", &m_gen_l2minus_eta);
		tree->Branch("GenLep4Eta", &m_gen_l2plus_eta);
		tree->Branch("GenLep1Phi", &m_gen_l1minus_phi);
		tree->Branch("GenLep2Phi", &m_gen_l1plus_phi);
		tree->Branch("GenLep3Phi", &m_gen_l2minus_phi);
		tree->Branch("GenLep4Phi", &m_gen_l2plus_phi);
		tree->Branch("GenLep1ID", &m_l1minus_id);
		tree->Branch("GenLep2ID", &m_l1plus_id);
		tree->Branch("GenLep3ID", &m_l2minus_id);
		tree->Branch("GenLep4ID", &m_l2plus_id);
		tree->Branch("GenHPt", &gen_pT4l);
		tree->Branch("GenHEta", &gen_eta4l);
		tree->Branch("GenHPhi", &gen_phi4l);
		tree->Branch("GenIntVtx_x", &GenIntVtx_x);
		tree->Branch("GenIntVtx_y", &GenIntVtx_y);
		tree->Branch("GenIntVtx_z", &GenIntVtx_z);

		tree->Branch("GenCTau", &gen_ctau); // Will be in microns

		double selfDHvvcoupl[SIZE_HVV][2] = { { 0 } };
		double ggvvcoupl[2] = { 0 };

		int myflavor=0;
		double nSelected = tinput->GetEntries("isSelected>0 && ZZMass>=105.6 && ZZMass<140.6");
		for(int ev=0;ev<tinput->GetEntries();ev++){
			tinput->GetEntry(ev);
			if(isSelected==0) continue;
			if(abs(m_l1minus_id)==abs(m_l2minus_id) && abs(m_l1minus_id)==13) myflavor=0;
			else if(abs(m_l1minus_id)==abs(m_l2minus_id) && abs(m_l1minus_id)==11) myflavor=1;
			else if(abs(m_l1minus_id)!=abs(m_l2minus_id) && abs(m_l1minus_id)!=15 && abs(m_l2minus_id)!=15) myflavor=2;
			else myflavor=4;
			if(myflavor>3) continue;

			int lepIdOrdered[4]={ m_l1minus_id,m_l1plus_id,m_l2minus_id,m_l2plus_id };
			float angularOrdered[8] = { m_zzmass, m_z1mass, m_z2mass, m_costhetastar, m_costheta1, m_costheta2, m_phi, m_phistar1 };

			mela.setProcess(TVar::SelfDefine_spin0, TVar::JHUGen, TVar::ZZGG);

			m_wt = PeakYield[EnergyIndex][myflavor]/nSelected;

			mela.setProcess(TVar::SelfDefine_spin0, TVar::JHUGen, TVar::ZZGG);
			selfDHvvcoupl[0][0]=1;
			p0plus_VAJHU = getJHUGenMELAWeight(mela, lepIdOrdered, angularOrdered, selfDHvvcoupl);

			mela.setProcess(TVar::HSMHiggs, TVar::MCFM, TVar::ZZGG);
			p0plus_VAMCFM = getMCFMMELAWeight(mela, lepIdOrdered, angularOrdered, ggvvcoupl);
			mela.setProcess(TVar::bkgZZ, TVar::MCFM, TVar::ZZQQB);
			bkg_VAMCFM = getMCFMMELAWeight(mela, lepIdOrdered, angularOrdered, ggvvcoupl,1); // |qqZZ|**2

			mela.setProcess(TVar::HSMHiggs, TVar::JHUGen, TVar::ZZGG);
			p0plus_m4l = getSuperMELA(mela,lepIdOrdered,angularOrdered[0],TVar::SMSyst_None);
			mela.setProcess(TVar::HSMHiggs, TVar::JHUGen, TVar::ZZGG);
			p0plus_m4l_ScaleUp = getSuperMELA(mela,lepIdOrdered,angularOrdered[0],TVar::SMSyst_ScaleUp);
			mela.setProcess(TVar::HSMHiggs, TVar::JHUGen, TVar::ZZGG);
			p0plus_m4l_ScaleDown = getSuperMELA(mela,lepIdOrdered,angularOrdered[0],TVar::SMSyst_ScaleDown);
			mela.setProcess(TVar::HSMHiggs, TVar::JHUGen, TVar::ZZGG);
			p0plus_m4l_ResUp = getSuperMELA(mela,lepIdOrdered,angularOrdered[0],TVar::SMSyst_ResUp);
			mela.setProcess(TVar::HSMHiggs, TVar::JHUGen, TVar::ZZGG);
			p0plus_m4l_ResDown = getSuperMELA(mela,lepIdOrdered,angularOrdered[0],TVar::SMSyst_ResDown);

			mela.setProcess(TVar::bkgZZ, TVar::JHUGen, TVar::ZZGG);
			bkg_m4l = getSuperMELA(mela,lepIdOrdered,angularOrdered[0],TVar::SMSyst_None);
			mela.setProcess(TVar::bkgZZ, TVar::JHUGen, TVar::ZZGG);
			bkg_m4l_ScaleUp = getSuperMELA(mela,lepIdOrdered,angularOrdered[0],TVar::SMSyst_ScaleUp);
			mela.setProcess(TVar::bkgZZ, TVar::JHUGen, TVar::ZZGG);
			bkg_m4l_ScaleDown = getSuperMELA(mela,lepIdOrdered,angularOrdered[0],TVar::SMSyst_ScaleDown);
			mela.setProcess(TVar::bkgZZ, TVar::JHUGen, TVar::ZZGG);
			bkg_m4l_ResUp = getSuperMELA(mela,lepIdOrdered,angularOrdered[0],TVar::SMSyst_ResUp);
			mela.setProcess(TVar::bkgZZ, TVar::JHUGen, TVar::ZZGG);
			bkg_m4l_ResDown = getSuperMELA(mela,lepIdOrdered,angularOrdered[0],TVar::SMSyst_ResDown);

			tree->Fill();
		}

		foutput->WriteTObject(tree);
		foutput->Close();
	}
}

