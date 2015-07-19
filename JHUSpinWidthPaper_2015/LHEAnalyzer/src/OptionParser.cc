#include "../interface/OptionParser.h"


OptionParser::OptionParser(int argc, char** argv):
mPOLE(125.), // mH
wPOLE(4.07), // GammaH
erg_tev(13), // C.o.M. energy in TeV
includeGenInfo(1), // Record gen. level quantities
includeRecoInfo(1), // Record reco. level quantities
removeDaughterMasses(1), // Force lepton masses to 0 in MELA
fileLevel(0), // LHE, Pythia
isHZZ(1), // H->ZZ or H->WW
decayMode(0), // 4l with HZZ, 2l2nu with HWW, see Event::constructVVCandidates(bool isZZ, int fstype).
indir("./"),
outdir("./"),
coutput("tmp.root")
{
  for (int a=0; a<argc; a++){
    string tmpArg(argv[a]);
    rawOptions.push_back(tmpArg);
  }
  analyze();
}
void OptionParser::analyze(){
  bool redefinedOutputFile=false;
  for (int opt=1; opt<rawOptions.size(); opt++){
    string wish, value;
    splitOption(rawOptions.at(opt), wish, value);
    interpretOption(wish, value);
    if (wish=="outfile") redefinedOutputFile=true;
  }

  if (filename.size()==0){ cerr << "You have to specify the input files." << endl; assert(0); }
  else{
    for (int f=0; f<filename.size(); f++){
      if ((filename.at(f).find(".lhe")!=string::npos && fileLevel==1) || (filename.at(f).find(".root")!=string::npos && fileLevel==0)){
        cerr << "Inconsistent fila name " << filename.at(f) << " and fileLevel option " << fileLevel << "!" << endl;
        assert(0);
      }
    }
  }
  if (includeGenInfo==0 && includeRecoInfo==0){ cerr << "Cannot omit both reco. and gen. level info." << endl; assert(0); }
  if (mPOLE==0 || wPOLE==0 || erg_tev==0){ cerr << "Cannot have mH, GammaH or sqrts == 0" << endl; assert(0); }
  if (!redefinedOutputFile) cout << "WARNING: No output file specified. Defaulting to " << coutput << "." << endl;
}
void OptionParser::splitOption(string rawoption, string& wish, string& value){
  size_t posEq = rawoption.find('=');
  if (posEq!=string::npos){
    wish=rawoption;
    value=rawoption.substr(posEq+1);
    wish.erase(wish.begin()+posEq,wish.end());
  }
  else{
    wish="";
    value=rawoption;
  }
}

void OptionParser::interpretOption(string wish, string value){
  if (wish.empty()){
    if (value.find(".lhe")!=string::npos || value.find(".root")!=string::npos) filename.push_back(value);
    else if (value.find("help")!= string::npos) printOptionsHelp();
    else cerr << "Unknown unspecified argument: " << value << endl;
  }
  else if (wish=="mH" || wish=="MH" || wish=="mPOLE") mPOLE = (double)atof(value.c_str());
  else if (wish=="GaH" || wish=="GammaH" || wish=="wPOLE") wPOLE = (double)atof(value.c_str());
  else if (wish=="sqrts") erg_tev = (int)atoi(value.c_str());
  else if (wish=="includeGenInfo") includeGenInfo = (int)atoi(value.c_str());
  else if (wish=="includeRecoInfo") includeRecoInfo = (int)atoi(value.c_str());
  else if (wish=="removeDaughterMasses") removeDaughterMasses = (int)atoi(value.c_str());
  else if (wish=="fileLevel") fileLevel = (int)atoi(value.c_str());
  else if (wish=="isHZZ") {
    isHZZ = (int)atoi(value.c_str());
    if (isHZZ==0) PDGHelpers::setHVVmass(PDGHelpers::Wmass);
    else PDGHelpers::setHVVmass(PDGHelpers::Zmass);
  }
  else if (wish=="decayMode") decayMode = (int)atoi(value.c_str());
  else if (wish=="indir") indir = value;
  else if (wish=="outdir") outdir = value;
  else if (wish=="outfile") coutput = value;
  else cerr << "Unknown specified argument: " << value << " with specifier " << wish << endl;
}

void OptionParser::printOptionsHelp(){

  assert(0);
}

void OptionParser::printOptionSummary(){

}

