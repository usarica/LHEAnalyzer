#ifndef BASETREE_H
#define BASETREE_H

#include <iostream>
#include <cstdio>
#include <string>
#include <vector>
#include <utility>
#include "TFile.h"
#include "TString.h"
#include "TTree.h"

using namespace std;

class BaseTree{
public:
  BaseTree(){ hvvtree=0; };
  BaseTree(string treename);
  BaseTree(string treename, string treetitle);
  ~BaseTree(){ if (hvvtree!=0) delete hvvtree; cleanBranches(); };

  // Innocuous functions
  void initTree(string treename, string treetitle){ hvvtree = new TTree(treename.c_str(), treetitle.c_str()); hvvtree->SetAutoSave(5000000000); };
  TTree* getTree(){ return hvvtree; };
  void record(){ hvvtree->Fill(); }
  void writeTree(TFile* foutput){ foutput->cd(); foutput->WriteTObject(hvvtree); }

  // This is where things get complicated pretty quickly
  enum BranchTypes{ bInt, bFloat, bVectorInt, bVectorDouble, nBranchTypes };
  bool bookBranch(string branchname, BaseTree::BranchTypes bType, bool doSetAddress);
  bool actuateBranches(bool doSetAddress);
  BranchTypes searchArray(string branchname, int& position);
  template<typename varType> void setVal(string branchname, varType value){
    int varposition=-1;
    BaseTree::BranchTypes varbranchtype = searchArray(branchname, varposition);
    if (varposition==-1 || varbranchtype==BaseTree::nBranchTypes){
      cerr << "Could not find the branch called " << branchname << "!" << endl;
    }
    else{
      if (varbranchtype==BaseTree::bInt) *(intBranches.at(varposition).second)=value;
      else if (varbranchtype==BaseTree::bFloat) *(floatBranches.at(varposition).second)=value;
      else if (varbranchtype==BaseTree::bVectorInt) vectorIntBranches.at(varposition).second->push_back(value);
      else if (varbranchtype==BaseTree::bVectorDouble) vectorDoubleBranches.at(varposition).second->push_back(value);
    }
  }
  void initializeBranches();
  void cleanBranches();



protected:
  TTree* hvvtree;

  vector < pair<string, Int_t*> > intBranches;
  vector < pair<string, Float_t*> > floatBranches;
  vector < pair<string, vector<int>*> > vectorIntBranches;
  vector < pair<string, vector<double>*> > vectorDoubleBranches;
};

#endif
