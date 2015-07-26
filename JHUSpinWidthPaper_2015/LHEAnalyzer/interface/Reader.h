#ifndef READER_H
#define READER_H

#include "converter.h"


class Reader : public converter{
public:
  Reader(OptionParser* options_);
  ~Reader(){};
  void run();

protected:
  void configure();
  void finalizeRun();
  void readEvent(TTree* tin, int ev, bool isGen, Event& outEvent);

  // This function obtains the result of type "returnType" for the Event ev in the context of Branch "branchname"
  // from an operation "evalVar" that receives a pointer to an Event (most general scenario)
  // and the address of a string "branchname" (so that it can be either set or passed inside) as argument,
  // records the result to the corresponding branch,
  // and returns a bool to indicate the success or failure of the operation.
  // The function prohibits the modification of the Event object.
  template<typename returnType> bool setVariable(const Event* ev, string& branchname, returnType(*evalVar)(const Event*, string&));

  void bindInputBranches(HVVTree* tin);
  void synchMappedBranches();
  void resetBranchBinding();

  // Vectors of branch name - evaluation function pointer pairs for streamlined evaluation of multiple functions and their recording to the output tree
  vector < pair<string, Int_t(*)(const Event*, string&)> > intFunctions;
  vector < pair<string, Float_t(*)(const Event*, string&)> > floatFunctions;
  vector < pair<string, vector<int>(*)(const Event*, string&)> > vectorIntFunctions;
  vector < pair<string, vector<double>(*)(const Event*, string&)> > vectorDoubleFunctions;

  // Vectors of branch pairs between input HVVTree (first) and output HVVTree (second)
  vector < pair<Int_t*, Int_t*> > intBranchMap;
  vector < pair<Float_t*, Float_t*> > floatBranchMap;
  vector < pair<vector<int>**, vector<int>*> > vectorIntBranchMap; // Map double pointer to pointer due to the object SetBranchAddress in the first HVVTree receives
  vector < pair<vector<double>**, vector<double>*> > vectorDoubleBranchMap; // Map double pointer to pointer due to the object SetBranchAddress in the first HVVTree receives
};
#endif
