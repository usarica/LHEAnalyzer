#ifndef READER_H
#define READER_H

#include "converter.h"


class Reader : public converter{
public:
  Reader(OptionParser* options_);
  ~Reader(){ cleanFunctions(); resetBranchBinding(); }
  void run();

  // This function obtains the result of type "returnType" for the MELAEvent ev in the context of Branch "branchname"
  // from an operation "evalVar" that receives a pointer to an MELAEvent (most general scenario)
  // and the address of a string "branchname" (so that it can be either set or passed inside) as argument,
  // records the result to the corresponding branch,
  // and returns a bool to indicate the success or failure of the operation.
  // The function prohibits the modification of the MELAEvent object and checks the return type via dynamic_cast to a pointer to it.
  template<typename returnType> bool setVariable(const MELAEvent* ev, std::string& branchname, returnType(*evalVar)(const MELAEvent*, std::string&));

  // Add the branch name - function pointer pairs to the arrays for easier access.
  template<typename returnType> void addFunction(std::string& branchname, returnType(*evalVar)(const MELAEvent*, std::string&, int)){
    using namespace std;

    int varposition=-1;
    BaseTree::BranchTypes varbranchtype = tree->searchArray(branchname, varposition);
    if (varposition==-1 || varbranchtype==BaseTree::nBranchTypes){
      cerr << "Reader::addFunction -> Could not find the branch called " << branchname << "!" << endl;
    }
    else{
      if (varbranchtype==BaseTree::bInt) { pair<string, Int_t(*)(const MELAEvent*, string&, int)> varPair(branchname, evalVar); intFunctions.push_back(varPair); }
      else if (varbranchtype==BaseTree::bFloat) { pair<string, Float_t(*)(const MELAEvent*, string&, int)> varPair(branchname, evalVar); floatFunctions.push_back(varPair); }
      else if (varbranchtype==BaseTree::bVectorInt) { pair<string, vectorInt(*)(const MELAEvent*, string&, int)> varPair(branchname, evalVar); vectorIntFunctions.push_back(varPair); }
      else if (varbranchtype==BaseTree::bVectorFloat) { pair<string, vectorFloat(*)(const MELAEvent*, string&, int)> varPair(branchname, evalVar); vectorFloatFunctions.push_back(varPair); }
      else if (varbranchtype==BaseTree::bVectorDouble) { pair<string, vectorDouble(*)(const MELAEvent*, string&, int)> varPair(branchname, evalVar); vectorDoubleFunctions.push_back(varPair); }
    }
  }

  // Clear function arrays
  void cleanFunctions(){ intFunctions.clear(); floatFunctions.clear(); vectorIntFunctions.clear(); vectorFloatFunctions.clear(); vectorDoubleFunctions.clear(); }

protected:
  void configure();
  void finalizeRun();
  void readEvent(MELAEvent& outEvent, std::vector<MELAParticle*>& particles, bool isGen);

  void bindInputBranches(HVVTree* tin);
  void synchMappedBranches();
  void resetBranchBinding();

  // Vectors of branch name - evaluation function pointer pairs for streamlined evaluation of multiple functions and their recording to the output tree
  std::vector < std::pair<std::string, Int_t(*)(const MELAEvent*, std::string&, int)> > intFunctions;
  std::vector < std::pair<std::string, Float_t(*)(const MELAEvent*, std::string&, int)> > floatFunctions;
  std::vector < std::pair<std::string, vectorInt(*)(const MELAEvent*, std::string&, int)> > vectorIntFunctions;
  std::vector < std::pair<std::string, vectorFloat(*)(const MELAEvent*, std::string&, int)> > vectorFloatFunctions;
  std::vector < std::pair<std::string, vectorDouble(*)(const MELAEvent*, std::string&, int)> > vectorDoubleFunctions;

  // Vectors of branch pairs between input HVVTree (first) and output HVVTree (second)
  std::vector < std::pair<Int_t*, Int_t*> > intBranchMap;
  std::vector < std::pair<Float_t*, Float_t*> > floatBranchMap;
  std::vector < std::pair<vectorInt**, vectorInt**> > vectorIntBranchMap; // Map double pointer to pointer due to the object SetBranchAddress in the first HVVTree receives
  std::vector < std::pair<vectorFloat**, vectorFloat**> > vectorFloatBranchMap; // Map double pointer to pointer due to the object SetBranchAddress in the first HVVTree receives
  std::vector < std::pair<vectorDouble**, vectorDouble**> > vectorDoubleBranchMap; // Map double pointer to pointer due to the object SetBranchAddress in the first HVVTree receives
};
#endif
