#include <fstream>
#include "TUtilHelpers.hh"
#include "OptionParser.h"
#include "HostHelpersCore.h"
#include "HelperFunctions.h"
#include "SampleHelpersCore.h"


using namespace std;
using namespace HelperFunctions;


OptionParser::OptionParser(int argc, char** argv) :
  mPOLE(125.), // mH
  wPOLE(4.07e-3), // GammaH
  wPOLEStandard(4.07e-3), // GHSM

  sampleXsec(-1), // xsec value
  sampleXsecErr(-1), // xsec error

  erg_tev(13), // C.o.M. energy in TeV
  includeGenInfo(1), // Record gen. level quantities
  includeRecoInfo(1), // Record reco. level quantities
  removeDaughterMasses(1), // Force lepton masses to 0 in MELA
  computeDecayAngles(1), // Decay angles
  computeVBFAngles(0), // VBF production angles
  computeVHAngles(0), // VH production angles
  computeTTHAngles(0), // VH production angles
  fileLevel(0), // -1: ReadMode, 0: LHE, 1: Pythia,
  pythiaStep(1), //0: GEN, 1: GEN-SIM
  isGenHZZ(MELAEvent::ZZMode), // H->ZZ or H->WW
  isRecoHZZ(MELAEvent::ZZMode), // H->ZZ or H->WW
  genDecayMode(0), // 4l with HZZ, 2l2nu with HWW, see Event::constructVVCandidates(bool isZZ, int fstype)
  recoDecayMode(0), // 4l with HZZ, 2l2nu with HWW, see Event::constructVVCandidates(bool isZZ, int fstype)
  recoSelBehaviour(0),
  recoSmearBehaviour(0),

  genHiggsCandidateSelectionScheme(HiggsComparators::BestZ1ThenZ2ScSumPt),
  recoHiggsCandidateSelectionScheme(HiggsComparators::BestZ1ThenZ2ScSumPt),

  jetDeltaRIso(0.4),
  jetAlgo("ak"),

  indir("./"),
  outdir("./"),
  coutput("tmp.root"),
  tmpDir("./tmpStore/"),
  maxEvents(-1)
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
  bool hasDecayAngles=false;
  bool hasJetAlgo=false;
  char rawdelimiter = '=';
  for (unsigned int opt=1; opt<rawOptions.size(); opt++){
    string wish, value;
    splitOption(rawOptions.at(opt), wish, value, rawdelimiter);
    interpretOption(wish, value);
    if (wish=="outfile") redefinedOutputFile=true;
    if (wish=="computeDecayAngles") hasDecayAngles=true;
    if (wish=="JetAlgorithm" || wish=="jetAlgorithm" || wish=="jetalgorithm") hasJetAlgo=true;
  }

  if (filename.empty()){ cerr << "You have to specify the input files." << endl; if (!hasInvalidOption) hasInvalidOption=true; }
  else{
    for (unsigned int f=0; f<filename.size(); f++){
      if ((filename.at(f).find(".lhe")!=string::npos && fileLevel!=0) || (filename.at(f).find(".root")!=string::npos && fileLevel==0)){
        cerr << "Inconsistent file name " << filename.at(f) << " and fileLevel option " << fileLevel << "!" << endl;
        if (!hasInvalidOption) hasInvalidOption=true;
      }
    }
  }

  if (maxEvents>=0){
    for (unsigned int es=0; es<eventSkipRanges.size(); es++) maxEvents += (eventSkipRanges.at(es).second-eventSkipRanges.at(es).first+1);
  }

  // Check for any invalid options and print an error
  if (isGenHZZ==-1){
    cout << "Gen. Higgs decay is disabled. Disabling MEs, decay angles and anything reco. as well." << endl;
    includeRecoInfo=0;
    hasDecayAngles=false; computeDecayAngles=0;
  }
  if (isRecoHZZ==-1){ cerr << "Reco. Higgs decay cannot be disabled." << endl; if (!hasInvalidOption) hasInvalidOption=true; }
  if (includeGenInfo==0 && includeRecoInfo==0){ cerr << "Cannot omit both reco. and gen. level info." << endl; if (!hasInvalidOption) hasInvalidOption=true; }
  if (mPOLE==0 || wPOLE==0 || erg_tev==0){ cerr << "Cannot have mH, GammaH or sqrts == 0" << endl; if (!hasInvalidOption) hasInvalidOption=true; }
  if (genHiggsCandidateSelectionScheme>=HiggsComparators::nCandidateSelections){ cerr << "Gen. H selection scheme is invalid!" << endl; if (!hasInvalidOption) hasInvalidOption=true; }
  if (recoHiggsCandidateSelectionScheme>=HiggsComparators::nCandidateSelections){ cerr << "Reco. H selection scheme is invalid!" << endl; if (!hasInvalidOption) hasInvalidOption=true; }
  if (!(((jetDeltaRIso==0.4 || jetDeltaRIso==0.5 || jetDeltaRIso==0.8) && jetAlgo=="ak") || ((jetDeltaRIso==0.4 || jetDeltaRIso==0.6) && jetAlgo=="kt"))){ cerr << "Jet algorithm can only be used with object isolations 0.4, 0.5 or 0.8 for ak, and 0.4 or 0.6 for kt jets at this moment." << endl; if (!hasInvalidOption) hasInvalidOption=true; }
  else{ jetAlgo.append(std::to_string((int) (10*jetDeltaRIso))); if (hasJetAlgo) cout << "Jet algorithm string " << jetAlgo << " has the isolation appended." << endl; }

  // Warnings-only
  if (!redefinedOutputFile) cout << "WARNING: No output file specified. Defaulting to " << coutput << "." << endl;
  if (!hasDecayAngles && fileLevel<0 && recoSelBehaviour==0) { cout << "Disabling the re-calculation of decay angles in ReadMode by default since no relevant option is specified." << endl; computeDecayAngles=0; }
  if (fileLevel<0 && recoSelBehaviour!=0) { cout << "Enabling the (re-)calculation of all angles in ReadMode since the re-selection option is specified." << endl; computeVBFAngles=1; computeVHAngles=1; computeTTHAngles=1; computeDecayAngles=1; }

  // Print help if needed and abort at this point, nowhere later
  if (hasInvalidOption) printOptionsHelp();

  // Append extra "/" if they do not exist.
  unsigned int tlen=(unsigned int) indir.length();
  if (tlen>1 && indir[tlen-1]!='/') indir.append("/");
  tlen=(unsigned int) outdir.length();
  if (tlen>1 && outdir[tlen-1]!='/') outdir.append("/");
  tlen=(unsigned int) tmpDir.length();
  if (tlen>1 && tmpDir[tlen-1]!='/') tmpDir.append("/");

  // Set isolation
  ParticleComparators::setJetDeltaR(jetDeltaRIso);

  // Initialize the global Mela if needed
  extractMElines();
  extractXsec();
  configureMela();
}
Bool_t OptionParser::isAnExcludedBranch(std::string const& branchname)const{
  bool isExcluded=false;
  for (auto const& exbranch:excludedBranch){
    if (branchname.find(exbranch)!=string::npos && !(branchname.find("Gen")!=string::npos && exbranch.find("Gen")==string::npos)){
      isExcluded=true;
      break;
    }
  }
  return isExcluded;
}
void OptionParser::extractSkippedEvents(std::string const& rawoption){
  vector<string> skipPairs;
  splitOptionRecursive(rawoption, skipPairs, ',');
  for (string const& skipPair:skipPairs){
    string strlow, strhigh;
    splitOption(skipPair, strlow, strhigh, '.');

    bool firstInclusive = true;
    size_t posFirstInc=strlow.find("[");
    size_t posFirstExc=strlow.find("(");
    if (posFirstInc!=string::npos || posFirstExc!=string::npos){
      if (posFirstInc!=string::npos && posFirstExc!=string::npos){
        cerr << "Invalid skipEvents range. Ignoring..." << endl;
        continue;
      }
      if (posFirstExc!=string::npos){
        firstInclusive=false;
        strlow = strlow.substr(posFirstExc+1);
      }
      else strlow = strlow.substr(posFirstInc+1);
    }
    Int_t firstId=(Int_t) atoi(strlow.c_str());
    if (!firstInclusive) firstId++;

    bool lastInclusive = false;
    size_t posLastInc=strhigh.find("]");
    size_t posLastExc=strhigh.find(")");
    if (posLastInc!=string::npos || posLastExc!=string::npos){
      if (posLastInc!=string::npos && posLastExc!=string::npos){
        cerr << "Invalid skipEvents range. Ignoring..." << endl;
        continue;
      }
      if (posLastExc!=string::npos) strhigh.erase(strhigh.begin()+posLastExc, strhigh.end());
      else{
        lastInclusive=true;
        strhigh.erase(strhigh.begin()+posLastInc, strhigh.end());
      }
    }
    Int_t lastId=(Int_t) atoi(strhigh.c_str());
    if (!lastInclusive) lastId--;

    if ((lastId>=0 && lastId<firstId) || (firstId<0 && lastId<0)){
      cerr << "Invalid skipEvents range. Ignoring..." << endl;
      continue;
    }
    pair<Int_t, Int_t> tmpPair(firstId, lastId);
    cout << "OptionParser: Will skip events " << tmpPair.first << " - " << tmpPair.second << endl;
    eventSkipRanges.push_back(tmpPair);
  }
}
void OptionParser::extractGlobalRecordSet(std::string const& rawoption){
  vector<string> compositePair;
  splitOptionRecursive(rawoption, compositePair, ',');
  for (string const& cpair:compositePair){
    string strname, strvalue;
    splitOption(cpair, strname, strvalue, ':');

    size_t posFirstBrac=strname.find("[");
    size_t posLastBrac=strvalue.find("]");
    if (posFirstBrac==string::npos || posLastBrac==string::npos){
      cerr << "Invalid global record branch specification. Ignoring..." << endl;
      continue;
    }

    strname = strname.substr(posFirstBrac+1);
    strvalue.erase(strvalue.begin()+posLastBrac, strvalue.end());
    pair<string, string> tmpPair(strname, strvalue);
    cout << "OptionParser: Will create branch " << tmpPair.first << " = " << tmpPair.second << " in the globals tree." << endl;
    globalRecordSet.push_back(tmpPair);
  }
}

void OptionParser::configureMela()const{
  Bool_t needMela = doMEComputations() || doComputeDecayAngles() || doComputeVBFAngles() || doComputeVHAngles() || doComputeTTHAngles();
  if (needMela) melaHelpers::melaHandle = new Mela((int)erg_tev, (float)mPOLE);
  melaHelpers::setSamplePoleWidth(wPOLE);
  melaHelpers::setStandardPoleWidth(wPOLEStandard);
  TUtil::applyLeptonMassCorrection(doRemoveLepMasses()); // Remains fixed, so not a problem to set it here
}
void OptionParser::deconfigureMela()const{
  if (melaHelpers::melaHandle) delete melaHelpers::melaHandle;
}
void OptionParser::extractMElines(){
  OptionParser::extractMElines(lheMEfile, lheMElist);
  OptionParser::extractMElines(recoMEfile, recoMElist);
}
void OptionParser::extractMElines(std::string const& sfile, std::vector<std::string>& llist){
  using namespace HostHelpers;
  if (!sfile.empty()){
    cout << "OptionParser::extractMElines: Attempting to read the ME file " << sfile << '.' << endl;
    if (FileReadable(sfile.c_str())){
      ifstream fin(sfile.c_str());
      if (fin.good()){
        while (!fin.eof()){
          string strline;
          getline(fin, strline);
          string strlinestrip=strline;
          lstrip(strlinestrip, " \"");
          rstrip(strlinestrip, " \",");
          if (strlinestrip.find('#')==0) continue;
          else if (!strlinestrip.empty()){
            cout << "\t- Adding ME option " << strlinestrip << "..." << endl;
            llist.push_back(strlinestrip);
          }
        }
      }
      fin.close();
    }
    else cerr << "OptionParser::extractMElines: ME file " << sfile << " is not readable." << endl;
  }
}
void OptionParser::extractXsec(){
  if (sampleXsec!=-1.f|| sampleXsecErr!=-1.f) return;

  std::vector<std::pair<float, float>> xsec_val_err; xsec_val_err.reserve(filename.size());
  for (auto fname:filename){
    fname = this->inputDir() + fname;
    cout << "OptionParser::extractXsec: Checking file " << fname  << " for xsec and xsecerr..." << endl;

    if (fileLevel==-1){
      TFile* fin = TFile::Open(fname.c_str(), "read");
      if (fin && fin->IsZombie()){
        if (fin->IsOpen()) fin->Close();
        delete fin;
      }
      else if (fin && !fin->IsOpen()) delete fin;
      else if (fin){ // File is open and good
        TTree* tin = (TTree*) fin->Get("SelectedTree");
        if (tin){
          using namespace SampleHelpers;
          tin->SetBranchStatus("*", 0);
          if (!branchExists(tin, "xsec") || !branchExists(tin, "xsecerr")) continue;
          Float_t xsecval=0, xsecerrval=0;
          bookBranch(tin, "xsec", &xsecval);
          bookBranch(tin, "xsecerr", &xsecerrval);
          tin->GetEntry(0);
          if (xsecerrval>0.f && xsecval!=0.f && std::isfinite(xsecval) && std::isfinite(xsecerrval)){
            bool isUnique = true;
            for (std::pair<float, float> const& tmp_pair:xsec_val_err){ if (tmp_pair.first == xsecval && tmp_pair.second == xsecerrval){ isUnique = false; break; } }
            if (isUnique) xsec_val_err.emplace_back(xsecval, xsecerrval);
          }
        }
        fin->Close();
      }
    }
    else if (fileLevel==0){
      ifstream fin;
      fin.open(fname.c_str());
      bool xsec_line_found = false;
      size_t init_line = 0;
      if (fin.good()){
        while (!fin.eof()){
          string const xsec_MCFM = "Cross-section is";

          string strline;
          getline(fin, strline);

          // Replace lines
          if (
            strline.find(xsec_MCFM)==string::npos
            && (strline.find("<init>")==string::npos && init_line==0)
            ) continue;
          else if (!xsec_line_found && (strline.find("<init>")!=string::npos || init_line>0)){
            //cout << "init block line: '" << strline << "' (" << init_line << ")" << endl;
            if (init_line==2){
              string strlinestrip=strline;
              lstrip(strlinestrip);
              rstrip(strlinestrip);
              stringstream ss(strlinestrip);
              float xsecval=0, xsecerrval=0;
              ss >> xsecval >> xsecerrval;
              if (xsecerrval>0.f && xsecval!=0.f && std::isfinite(xsecval) && std::isfinite(xsecerrval)){
                bool isUnique = true;
                for (std::pair<float, float> const& tmp_pair:xsec_val_err){ if (tmp_pair.first == xsecval && tmp_pair.second == xsecerrval){ isUnique = false; break; } }
                if (isUnique) xsec_val_err.emplace_back(xsecval, xsecerrval);
              }
              break;
            }
            init_line++;
          }
          else{
            replaceString<std::string, const std::string>(strline, xsec_MCFM, "");

            string strlinestrip=strline;
            lstrip(strlinestrip, " \":,()#");
            rstrip(strlinestrip, " \":,()#");
            if (!strlinestrip.empty()){
              replaceString<std::string, const char*>(strlinestrip, " ", "");
              replaceString<std::string, const char*>(strlinestrip, "+-", "|");
              {
                string xsec, xsecerr;
                splitOption(strlinestrip, xsec, xsecerr, '|');
                float xsecval=0, xsecerrval=0;
                try{ xsecval = stoi(xsec.c_str()); }
                catch (std::invalid_argument& e){
                  cerr << "OptionParser::extractXsec: Could not interpret the cross section string '" << xsec << "'" << endl;
                  xsecval=0;
                }
                try{ xsecerrval = stoi(xsecerr.c_str()); }
                catch (std::invalid_argument& e){
                  cerr << "OptionParser::extractXsec: Could not interpret the cross section error string '" << xsecerr << "'" << endl;
                  xsecerrval=0;
                }
                if (xsecerrval>0.f && xsecval!=0.f && std::isfinite(xsecval) && std::isfinite(xsecerrval)){
                  bool isUnique = true;
                  for (std::pair<float, float> const& tmp_pair:xsec_val_err){ if (tmp_pair.first == xsecval && tmp_pair.second == xsecerrval){ isUnique = false; break; } }
                  if (isUnique) xsec_val_err.emplace_back(xsecval, xsecerrval);
                }
              }
            }
            xsec_line_found = true;
            break;
          }
        }
      }
      fin.close();
    }
  }
  if (!xsec_val_err.empty()){
    sampleXsec=sampleXsecErr=0;
    for (std::pair<float, float> const& pp:xsec_val_err){
      float const wgt = pow(pp.second, -2);
      sampleXsec += pp.first*wgt;
      sampleXsecErr += wgt;
    }
    sampleXsec /= sampleXsecErr;
    sampleXsecErr = 1.f/sqrt(sampleXsecErr);
    cout << "OptionParser::extractXsec: Final cross section estimate is " << sampleXsec << " +- " << sampleXsecErr << " pb." << endl;
  }
}
Bool_t OptionParser::checkListVariable(std::vector<std::string> const& list, std::string const& var)const{ return TUtilHelpers::checkElementExists(var, list); }

void OptionParser::interpretOption(std::string const& wish, std::string const& value){
  if (wish.empty()){
    if (value.find(".lhe")!=string::npos || value.find(".root")!=string::npos) filename.push_back(value);
    else if (value.find("help")!= string::npos) printOptionsHelp();
    else cerr << "Unknown unspecified argument: " << value << endl;
  }
  else if (wish=="includeGenInfo") includeGenInfo = (int)atoi(value.c_str());
  else if (wish=="includeRecoInfo") includeRecoInfo = (int)atoi(value.c_str());
  else if (wish=="fileLevel") fileLevel = (int)atoi(value.c_str());
  else if (wish=="pythiaStep" || wish=="pythialevel") pythiaStep = (int)atoi(value.c_str());
  else if (wish=="isGenHZZ"){
    isGenHZZ = MELAEvent::getCandidateVVModeFromString(value);

    if (isGenHZZ==MELAEvent::WWMode) PDGHelpers::setCandidateDecayMode(TVar::CandidateDecay_WW);
    else if (isGenHZZ==MELAEvent::ZZMode) PDGHelpers::setCandidateDecayMode(TVar::CandidateDecay_ZZ);
    else if (isGenHZZ==MELAEvent::ZGammaMode) PDGHelpers::setCandidateDecayMode(TVar::CandidateDecay_ZG);
    else if (isGenHZZ==MELAEvent::GammaGammaMode) PDGHelpers::setCandidateDecayMode(TVar::CandidateDecay_GG);
    else PDGHelpers::setCandidateDecayMode(TVar::CandidateDecay_ff);
  }
  else if (wish=="isRecoHZZ") isRecoHZZ = MELAEvent::getCandidateVVModeFromString(value);
  else if (wish=="genDecayMode") genDecayMode = (int)atoi(value.c_str());
  else if (wish=="recoDecayMode") recoDecayMode = (int)atoi(value.c_str());
  else if (wish=="recoSelBehaviour" || wish=="recoSelBehavior") recoSelBehaviour = (int)atoi(value.c_str());
  else if (wish=="recoSmearBehaviour" || wish=="recoSmearBehavior") recoSmearBehaviour = (int)atoi(value.c_str());
  else if (wish=="genCandidateSelection" || wish=="genCandSel"){
    if (value=="BestZ1ThenZ2" || value=="BestZ1ThenZ2ScSumPt") genHiggsCandidateSelectionScheme = HiggsComparators::BestZ1ThenZ2ScSumPt;
  }
  else if (wish=="recoCandidateSelection" || wish=="recoCandSel"){
    if (value=="BestZ1ThenZ2" || value=="BestZ1ThenZ2ScSumPt") recoHiggsCandidateSelectionScheme = HiggsComparators::BestZ1ThenZ2ScSumPt;
  }

  else if (wish=="indir") indir = value;
  else if (wish=="outdir") outdir = value;
  else if (wish=="tmpDir" || wish=="tempDir") tmpDir = value;

  else if (wish=="outfile") coutput = value;
  else if (wish=="lheMEfile") lheMEfile = value;
  else if (wish=="recoMEfile") recoMEfile = value;

  else if (wish=="excludeBranch") splitOptionRecursive(value, excludedBranch, ',');
  else if (wish=="maxevents" || wish=="maxEvents") maxEvents = (int)atoi(value.c_str());
  else if (wish=="skipevents" || wish=="skipEvents") extractSkippedEvents(value);

  else if (wish=="globalRecord") extractGlobalRecordSet(value);

  else if (wish=="mH" || wish=="MH" || wish=="mPOLE") mPOLE = (double)atof(value.c_str());
  else if (wish=="GH" || wish=="GaH" || wish=="GammaH" || wish=="wPOLE") wPOLE = (double)atof(value.c_str());
  else if (wish=="GHSM" || wish=="GaHSM" || wish=="GammaHSM" || wish=="wPOLEStandard") wPOLEStandard = (double) atof(value.c_str());
  else if (wish=="xsec") sampleXsec = (float) atof(value.c_str());
  else if (wish=="xsecerr") sampleXsecErr = (float) atof(value.c_str());
  else if (wish=="sqrts") erg_tev = (int)atoi(value.c_str());
  else if (wish=="JetAlgorithm" || wish=="jetAlgorithm" || wish=="jetalgorithm") jetAlgo = value;
  else if (wish=="jetDeltaR" || wish=="jetIso" || wish=="jetIsolation" || wish=="jetDeltaRIso" || wish=="jetDeltaRIsolation") jetDeltaRIso = (double)atof(value.c_str());
  else if (wish=="removeDaughterMasses") removeDaughterMasses = (int) atoi(value.c_str());
  else if (wish=="computeDecayAngles") computeDecayAngles = (int)atoi(value.c_str());
  else if (wish=="computeVBFProdAngles") computeVBFAngles = (int)atoi(value.c_str());
  else if (wish=="computeVHProdAngles") computeVHAngles = (int)atoi(value.c_str());
  else if (wish=="computeTTHProdAngles") computeTTHAngles = (int) atoi(value.c_str());

  else cerr << "Unknown specified argument: " << value << " with specifier " << wish << endl;
}

void OptionParser::printOptionsHelp()const{
  cout << endl;
  cout << "The options implemented for the LHEAnalyzer (format: specifier=value):\n\n";

  cout << "- No option specifier: Input files with extension .lhe or .root. Multiple input files can be passed as different arguments.\n\n";
  cout << "- indir: Location of input files. Default=\"./\"\n\n";
  cout << "- fileLevel: -1==ReadMode, 0==LHE, 1==Pythia8. \".lhe\" extension only allowed for 0, and \".root\" is the only format for the others. Default=0\n\n";
  cout << "- pythiaStep: Level of reconstruction to use for fileLevel==1: 0==GEN, 1==GEN-SIM, 2==MINIAODSIM. A value of -1 interprets the inputs as in other fileLevel options. Default=1\n\n";
  cout << "- outfile: Output file name. Default=\"tmp.root\"\n\n";
  cout << "- outdir: Location of the output file. Default=\"./\"\n\n";
  cout << "- tmpDir: Location of temporary files. Default=\"./tmpStore/\"\n\n";
  cout << "- maxevents / maxEvents: Maximum number of events to process. Default=-1 (all events)\n\n";
  cout << "- skipevents / skipEvents: Events to skip at the beginning. Default=none.\n\tAssignment is made in the form [ev1.ev2],[ev3.ev4] to allow multiple ranges. Use [ or ] for inclusive, ( or ) for exclusive ranges. Counting the events begins from 0.\n\n";

  cout << "- sqrts: pp collision c.o.m. energy. Default=13 (TeV)\n\n";
  cout << "- removeDaughterMasses: Switch to control the removal of lepton masses in the angle computation. Default=1\n\n";
  cout << "- computeDecayAngles: Switch to control the decay angles computation. Default=1 in LHE or Pythia modes, 0 in ReadMode.\n\n";
  cout << "- computeVBFProdAngles: Switch to control the VBF production angles computation. Default=0\n\n";
  cout << "- computeVHProdAngles: Switch to control the VH production angles computation. Default=0.\n\tPossible values are 0 (==disable), 1 (==compute from jets only), 2 (==compute from leptons only), 3 (==compute from jets and leptons, with separate variable suffixes \"*_VHhadronic\" and \"*_VHleptonic\".).\n\n";
  cout << "- computeTTHProdAngles: Switch to control the ttH production angles computation. Default=0.\n\tPossible values are 0 (==disable), 1 (==compute from jets and leptons, with same variable suffixes\".).\n\n";

  cout << "- mH / MH / mPOLE: Mass of the Higgs. Used in common for generator and reco. objects. Default=125 (GeV)\n\n";
  cout << "- GH / GaH / GammaH / wPOLE: Width of the generated Higgs. Used in generator objects. Default=4.07 (MeV)\n\n";
  cout << "- GHSM / GaHSM / GammaHSM / wPOLEStandard: Standard SM width. Used in scaling Mela probabilities properly. Default=4.07 (MeV).\n\n";
  cout << "- xsec: Assign a cross section to the sample. Will override the OptionParser determination. Default=-1 (pb)\n\n";
  cout << "- xsecerr: Assign a cross section error to the sample. Will override the OptionParser determination. Default=-1 (pb)\n\n";
  cout << "- includeGenInfo, includeRecoInfo: Flags to control the writing of gen. and reco. info., respectively. Cannot be both false (0). Default=(1, 1)\n\n";
  cout << "- isGenHZZ, isRecoHZZ: Gen. or reco. candidate decay hypothesis. Undecayed (Higgs, gen.-only), WW, ZZ, ff (or ffb), Zgamma (or Zgam), gammagamma (or gamgam), Z (->ffb). isGenHZZ also (re)sets the default V mass in H->VV decay. Defaults=(ZZ, ZZ)\n\n";
  cout << "- genDecayMode, recoDecayMode: Gen. or reco. H->VV->final states. Behavior changes based on isGenHZZ and isRecoHZZ. Defaults=(0, 0)\n";
  MELAEvent::printCandidateDecayModeDescriptions(); cout << "\n";

  cout << "- recoSelBehavior / recoSelBehaviour: Selection behaviour on all reco. final states. Default=0.\n\t0==Apply selection in LHE and Pythia modes, apply no re-selection in ReadMode.\n\t1==Opposite of 0. Also enables the computation of all angles in ReadMode, overriding the relevant command line options.\n\n";
  cout << "- recoSmearBehavior / recoSmearBehaviour: Smearing behaviour on all reco. final states. Does not apply to ReadMode. Default=0.\n\t0==Apply smearing in LHE mode, no smearing in Pythia mode\n\t1==Opposite of 0\n\n";
  cout << "- genCandidateSelection, recoCandidateSelection: Higgs candidate selection algorithm. Values accepted are\n\t->BestZ1ThenZ2 (=BestZ1ThenZ2ScSumPt).\n\tDefaults==(BestZ1ThenZ2, BestZ1ThenZ2)\n\n";

  cout << "- JetAlgorithm / jetAlgorithm / jetalgorithm: Jet algorithm to use if available in the input tree. Isolation needs to be set separately if different from the default value. Default=ak\n\n";
  cout << "- jetDeltaR / jetIso / jetIsolation / jetDeltaRIso / jetDeltaRIsolation: deltaR_jet isolation cut used in selecting jets. This value (x10) is appended to the jet algorithm string and can only be 0.4, 0.5 or 0.8 for ak, and 0.4 or 0.6 for kt jets at the moment. Default=0.5\n\n";

  cout << "- excludeBranch: Comma-separated list of excluded branches. Default is to include all branches called via HVVTree::bookAllBranches.\n\n";

  cout << "- lheMEfile: File listing truth-level MEs. Default=\"\"\n\n";
  cout << "- recoMEfile: File listing reco-level MEs. Default=\"\"\n\n";

  cout << "- globalRecord: Global values to set (e.g. cross section). Creates SelectedTree_Globals with a single event. Default=none.\n\tBranches are assigned in the form [name:type_value], and multiple forms can be specified with comma separation. Use C++ type names (e.g. [xsec:float_0.001].\n\n";

  cout << endl;
  assert(0);
}

void OptionParser::printOptionSummary()const{

}

