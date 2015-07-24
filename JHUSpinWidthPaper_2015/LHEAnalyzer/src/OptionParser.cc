#include "../interface/OptionParser.h"


OptionParser::OptionParser(int argc, char** argv) :
mPOLE(125.), // mH
wPOLE(4.07e-3), // GammaH
erg_tev(13), // C.o.M. energy in TeV
includeGenInfo(1), // Record gen. level quantities
includeRecoInfo(1), // Record reco. level quantities
removeDaughterMasses(1), // Force lepton masses to 0 in MELA
fileLevel(0), // LHE, Pythia
isGenHZZ(1), // H->ZZ or H->WW
isRecoHZZ(1), // H->ZZ or H->WW
genDecayMode(0), // 4l with HZZ, 2l2nu with HWW, see Event::constructVVCandidates(bool isZZ, int fstype)
recoDecayMode(0), // 4l with HZZ, 2l2nu with HWW, see Event::constructVVCandidates(bool isZZ, int fstype)
genHiggsCandidateSelectionScheme(HiggsComparators::BestZ1ThenZ2ScSumPt),
recoHiggsCandidateSelectionScheme(HiggsComparators::BestZ1ThenZ2ScSumPt),
indir("./"),
outdir("./"),
tmpDir("./tmpStore/"),
coutput("tmp.root")
{
  for (int a=0; a<argc; a++){
    string tmpArg(argv[a]);
    rawOptions.push_back(tmpArg);
  }
  analyze();
}
void OptionParser::analyze(){
  bool hasInvalidOption=false;
  bool redefinedOutputFile=false;
  for (int opt=1; opt<rawOptions.size(); opt++){
    string wish, value;
    splitOption(rawOptions.at(opt), wish, value);
    interpretOption(wish, value);
    if (wish=="outfile") redefinedOutputFile=true;
  }

  if (filename.size()==0){ cerr << "You have to specify the input files." << endl; if(!hasInvalidOption) hasInvalidOption=true; }
  else{
    for (int f=0; f<filename.size(); f++){
      if ((filename.at(f).find(".lhe")!=string::npos && fileLevel==1) || (filename.at(f).find(".root")!=string::npos && fileLevel==0)){
        cerr << "Inconsistent fila name " << filename.at(f) << " and fileLevel option " << fileLevel << "!" << endl;
        if(!hasInvalidOption) hasInvalidOption=true;
      }
    }
  }
  if (includeGenInfo==0 && includeRecoInfo==0){ cerr << "Cannot omit both reco. and gen. level info." << endl; if(!hasInvalidOption) hasInvalidOption=true; }
  if (mPOLE==0 || wPOLE==0 || erg_tev==0){ cerr << "Cannot have mH, GammaH or sqrts == 0" << endl; if(!hasInvalidOption) hasInvalidOption=true; }
  if (genHiggsCandidateSelectionScheme>=HiggsComparators::nCandidateSelections){ cerr << "Gen. H selection scheme is invalid!" << endl; if(!hasInvalidOption) hasInvalidOption=true; }
  if (recoHiggsCandidateSelectionScheme>=HiggsComparators::nCandidateSelections){ cerr << "Reco. H selection scheme is invalid!" << endl; if(!hasInvalidOption) hasInvalidOption=true; }
  if (!redefinedOutputFile) cout << "WARNING: No output file specified. Defaulting to " << coutput << "." << endl;

  if (hasInvalidOption) printOptionsHelp();
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
  else if (wish=="GH" || wish=="GaH" || wish=="GammaH" || wish=="wPOLE") wPOLE = (double)atof(value.c_str());
  else if (wish=="sqrts") erg_tev = (int)atoi(value.c_str());
  else if (wish=="includeGenInfo") includeGenInfo = (int)atoi(value.c_str());
  else if (wish=="includeRecoInfo") includeRecoInfo = (int)atoi(value.c_str());
  else if (wish=="removeDaughterMasses") removeDaughterMasses = (int)atoi(value.c_str());
  else if (wish=="fileLevel") fileLevel = (int)atoi(value.c_str());
  else if (wish=="isGenHZZ") {
    isGenHZZ = (int)atoi(value.c_str());
    if (isGenHZZ==0) PDGHelpers::setHVVmass(PDGHelpers::Wmass);
    else PDGHelpers::setHVVmass(PDGHelpers::Zmass);
  }
  else if (wish=="isRecoHZZ") isRecoHZZ = (int)atoi(value.c_str());
  else if (wish=="genDecayMode") genDecayMode = (int)atoi(value.c_str());
  else if (wish=="recoDecayMode") recoDecayMode = (int)atoi(value.c_str());
  else if (wish=="genCandidateSelection" || wish=="genCandSel"){
    if (value=="BestZ1ThenZ2" || value=="BestZ1ThenZ2ScSumPt") genHiggsCandidateSelectionScheme = HiggsComparators::BestZ1ThenZ2ScSumPt;
  }
  else if (wish=="recoCandidateSelection" || wish=="recoCandSel"){
    if (value=="BestZ1ThenZ2" || value=="BestZ1ThenZ2ScSumPt") recoHiggsCandidateSelectionScheme = HiggsComparators::BestZ1ThenZ2ScSumPt;
  }
  else if (wish=="indir") indir = value;
  else if (wish=="outdir") outdir = value;
  else if (wish=="outfile") coutput = value;
  else if (wish=="tmpDir" || wish=="tempDir") tmpDir = value;
  else cerr << "Unknown specified argument: " << value << " with specifier " << wish << endl;
}

void OptionParser::printOptionsHelp(){
  cout << endl;
  cout << "The options implemented for the LHEAnalyzer (format: specifier=value):\n";

  cout << "No option specifier: Input files with extension .lhe or .root. Multiple input files can be passed as different arguments.\n";
  cout << "indir: Location of input files. Default=\"./\"\n";
  cout << "fileLevel: 0==LHE (no .root extensions), 1==Pythia8 (no .lhe extensions). Default=0\n";
  cout << "outfile: Output file name. Default=\"tmp.root\"\n";
  cout << "outdir: Location of the output file. Default=\"./\"\n";
  cout << "tmpDir: Location of temporary files. Default=\"./tmpStore/\"\n";

  cout << "mH / MH / mPOLE: Mass of the Higgs. Used in common for generator and reco. objects. Default=125 (GeV)\n";
  cout << "GH / GaH / GammaH / wPOLE: Width of the generated Higgs. Used in generator objects. Default=4.07 (MeV)\n";
  cout << "sqrts: Width of the generated Higgs. Used in generator objects. Default=13 (TeV)\n";
  cout << "includeGenInfo, includeRecoInfo: Flags to control the writing of gen. and reco. info., respctively. Cannot be both false (0). Default=(1, 1)\n";
  cout << "removeDaughterMasses: Flag to control the removal of lepton masses in the angle computation. Default=1\n";
  cout << "isGenHZZ, isRecoHZZ: Gen. or reco. H->VV decay. 0==H->ZZ decay, 1==H->WW decay. Defaults=(0, 0)\n";
  cout << "genDecayMode, recoDecayMode: Gen. or reco. H->VV->final states. Defaults=(0, 0)\n\tIf H->ZZ decay is specified, 0-5==4l, 4q, 2l2q, 2l2nu, 2q2nu, 4nu.\n\tIf H->WW decay is specified, 0-2==2l2nu, 4nu, lnu2q.\n";

  cout << "genCandidateSelection, recoCandidateSelection: Higgs candidate selection algorithm. Values accepted are\n\t->BestZ1ThenZ2 (=BestZ1ThenZ2ScSumPt).\n\tDefaults==(BestZ1ThenZ2, BestZ1ThenZ2)\n";


  cout << endl;
  assert(0);
}

void OptionParser::printOptionSummary(){

}

