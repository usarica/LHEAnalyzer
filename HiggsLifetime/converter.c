#include <iostream>
#include <fstream>
#include <cstdio>
#include <string>
#include <iomanip>
#include "./data/ZZ4l_125_6_Samples.h"

using namespace std;

void converter_single(int erg_tev, int smp);
void converter_main(ifstream& input_lhe, ofstream& output_txt);

void converter(){
	for (int erg_tev = 7; erg_tev < 9; erg_tev++){
		for (int smp = 0; smp < kNumSamples; smp++) converter_single(erg_tev, smp);
	}
}


void converter_single(int erg_tev, int smp){
	string input_main_folder = "LHC_";
	string cergtev;
	if(erg_tev==7) cergtev = "7TeV";
	else if(erg_tev==8) cergtev = "8TeV";
	input_main_folder = input_main_folder + cergtev;
	string user_main_dir = user_dir + input_main_folder;
	string user_main_dir_hep = user_dir_hep + input_main_folder;
	string user_main_dir_hhpc = user_dir_hhpc + input_main_folder + "/" + input_main_folder + "_0";

	string user_input_dir_used = user_main_dir_hep;
	if(smp==0) user_input_dir_used = user_main_dir_hhpc;
	string user_output_dir_used = user_main_dir_hep;

	string cinput_core = user_input_dir_used + "/";
	if(smp>0) cinput_core = cinput_core + sample_suffix[smp] + "/";
	string coutput_core = user_output_dir_used + "/";
	coutput_core = coutput_core + sample_suffix[smp] + "/";

	string filename = "HZZ4l-";
	filename = filename + cergtev;
	filename = filename + "_0+m_";
	string outfilename = filename;
	outfilename = outfilename + sample_suffix[smp];

	string coutput_main = coutput_core;
	coutput_main += outfilename;

	char coutput[1000];
	sprintf(coutput,"%s%s",coutput_main.c_str(),".txt");
	ofstream output_txt(coutput,ios::out);
	cout << coutput << endl;
	
	int nFiles=1;
	if(smp==0) nFiles=5;
	int fileCtr=0;
	while (fileCtr<nFiles){
		string infilename = filename;
		char indexFile[3];
		sprintf(indexFile,"%i", fileCtr);
		if (smp>0) infilename = infilename + sample_suffix[smp];
		else infilename = infilename + indexFile;

		string cinput_main = cinput_core;
		cinput_main += infilename;

		char cinput[1000];
		sprintf(cinput, "%s%s", cinput_main.c_str(), ".lhe");
		ifstream input_lhe(cinput, ios::in);
		cout << cinput << endl;
		if (!input_lhe.good()) assert(0);

		converter_main(input_lhe, output_txt);
		input_lhe.close();
		fileCtr++;
	}

	output_txt.close();
	cout << "End of program..." << endl;
}

void converter_main(ifstream& input_lhe, ofstream& output_txt){
	bool specifications=true;
	int event_ctr=0;

	string event_beginning = "<event>";
	string event_end = "</event>";
	string file_closing = "</LesHouchesEvents>";
	while(!input_lhe.eof()){
		string str_in;
		string str_out;
		getline(input_lhe,str_in);
		size_t event_begins=str_in.find(event_beginning);
		if(event_begins!=string::npos && specifications){
			specifications=false;
			getline(input_lhe,str_in);
			++event_ctr;
		}

		size_t event_ends=str_in.find(event_end);
		if(event_ends!=string::npos && !specifications) specifications=true;
		size_t file_ends=str_in.find(file_closing);
		if(file_ends!=string::npos && !specifications) specifications=true;
// 11, 13, 15: e, mu, tau; 25: Higgs
		bool entry_is_good = false;
		size_t higgs = str_in.find("25");
		size_t antihiggs = str_in.find("-25");
		size_t electron = str_in.find("11");
		size_t muon = str_in.find("13");
		size_t tau = str_in.find("15");
		size_t antielectron = str_in.find("-11");
		size_t antimuon = str_in.find("-13");
		size_t antitau = str_in.find("-15");

		if(higgs==1 || antihiggs==1) entry_is_good=true;
		if(electron==1 || muon==1 || tau==1) entry_is_good=true;
		if(antielectron==1 || antimuon==1 || antitau==1) entry_is_good=true;
		if(!specifications && entry_is_good){
			str_out = str_in;
			output_txt << str_out << endl;
		}
	}
}
