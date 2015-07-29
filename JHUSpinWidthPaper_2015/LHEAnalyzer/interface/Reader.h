#ifndef READER_H
#define READER_H

#include "converter.h"


class Reader : public converter{
public:
  Reader(OptionParser* options_);
  ~Reader(){ cleanFunctions(); resetBranchBinding(); }
  void run();

  // This function obtains the result of type "returnType" for the Event ev in the context of Branch "branchname"
  // from an operation "evalVar" that receives a pointer to an Event (most general scenario)
  // and the address of a string "branchname" (so that it can be either set or passed inside) as argument,
  // records the result to the corresponding branch,
  // and returns a bool to indicate the success or failure of the operation.
  // The function prohibits the modification of the Event object and checks the return type via dynamic_cast to a pointer to it.
  template<typename returnType> bool setVariable(const Event* ev, string& branchname, returnType(*evalVar)(const Event*, string&));

  // Add the branch name - function pointer pairs to the arrays for easier access.
  template<typename returnType> void addFunction(string branchname, returnType(*evalVar)(const Event*, string&, int)){
    int varposition=-1;
    BaseTree::BranchTypes varbranchtype = tree->searchArray(branchname, varposition);
    if (varposition==-1 || varbranchtype==BaseTree::nBranchTypes){
      cerr << "Reader::addFunction -> Could not find the branch called " << branchname << "!" << endl;
    }
    else{
      if (varbranchtype==BaseTree::bInt) { pair<string, Int_t(*)(const Event*, string&, int)> varPair(branchname, evalVar); intFunctions.push_back(varPair); }
      else if (varbranchtype==BaseTree::bFloat) { pair<string, Float_t(*)(const Event*, string&, int)> varPair(branchname, evalVar); floatFunctions.push_back(varPair); }
      else if (varbranchtype==BaseTree::bVectorInt) { pair<string, vectorInt(*)(const Event*, string&, int)> varPair(branchname, evalVar); vectorIntFunctions.push_back(varPair); }
      else if (varbranchtype==BaseTree::bVectorDouble) { pair<string, vectorDouble(*)(const Event*, string&, int)> varPair(branchname, evalVar); vectorDoubleFunctions.push_back(varPair); }
    }
  }

  // Clear function arrays
  void cleanFunctions(){ intFunctions.clear(); floatFunctions.clear(); vectorIntFunctions.clear(); vectorDoubleFunctions.clear(); }

protected:
  void configure();
  void finalizeRun();
  void readEvent(Event& outEvent, vector<Particle*>& particles, bool isGen);

  void bindInputBranches(HVVTree* tin);
  void synchMappedBranches();
  void resetBranchBinding();

  // Vectors of branch name - evaluation function pointer pairs for streamlined evaluation of multiple functions and their recording to the output tree
  vector < pair<string, Int_t(*)(const Event*, string&, int)> > intFunctions;
  vector < pair<string, Float_t(*)(const Event*, string&, int)> > floatFunctions;
  vector < pair<string, vectorInt(*)(const Event*, string&, int)> > vectorIntFunctions;
  vector < pair<string, vectorDouble(*)(const Event*, string&, int)> > vectorDoubleFunctions;

  // Vectors of branch pairs between input HVVTree (first) and output HVVTree (second)
  vector < pair<Int_t*, Int_t*> > intBranchMap;
  vector < pair<Float_t*, Float_t*> > floatBranchMap;
  vector < pair<vectorInt**, vectorInt**> > vectorIntBranchMap; // Map double pointer to pointer due to the object SetBranchAddress in the first HVVTree receives
  vector < pair<vectorDouble**, vectorDouble**> > vectorDoubleBranchMap; // Map double pointer to pointer due to the object SetBranchAddress in the first HVVTree receives
};
#endif
