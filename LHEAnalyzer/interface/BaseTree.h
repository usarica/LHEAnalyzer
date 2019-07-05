#ifndef BASETREE_H
#define BASETREE_H

#include <iostream>
#include <cstdio>
#include <string>
#include <vector>
#include <utility>
#include <algorithm>
#include "TFile.h"
#include "TString.h"
#include "TTree.h"


typedef std::vector<int> vectorInt;
typedef std::vector<float> vectorFloat;
typedef std::vector<double> vectorDouble;


class BaseTree{
public:
  BaseTree(){ hvvtree=nullptr; }
  BaseTree(std::string treename);
  BaseTree(std::string treename, std::string treetitle);
  BaseTree(std::string treename, TFile* fin);
  virtual ~BaseTree(){ delete hvvtree; cleanBranches(); }

  // Innocuous functions
  void initTree(std::string treename, std::string treetitle){ hvvtree = new TTree(treename.c_str(), treetitle.c_str()); hvvtree->SetAutoSave(5000000000); }
  void getTreeFromFile(std::string treename, TFile* fin);
  TTree* getTree(){ return hvvtree; }
  void record(){ hvvtree->Fill(); }
  void writeTree(TFile* foutput){ foutput->cd(); foutput->WriteTObject(hvvtree); }

  // This is where things get complicated pretty quickly
  enum BranchTypes{
    bInt,
    bFloat,
    bVectorInt,
    bVectorFloat,
    bVectorDouble,
    nBranchTypes
  };
  bool bookBranch(std::string& branchname, const BaseTree::BranchTypes& bType, const bool& doSetAddress);
  bool actuateBranches(const bool& doSetAddress);
  std::vector<std::string> getBranchList()const;
  BranchTypes searchArray(const std::string& branchname, int& position);
  template<typename varType> void setVal(const std::string& branchname, const varType& value){
    int varposition=-1;
    BaseTree::BranchTypes varbranchtype = searchArray(branchname, varposition);
    if (varposition==-1 || varbranchtype==BaseTree::nBranchTypes){
      std::cerr << "BaseTree::setVal -> Could not find the branch called " << branchname << "!" << std::endl;
    }
    else{
      if (varbranchtype==BaseTree::bInt) *(intBranches.at(varposition).second)=value;
      else if (varbranchtype==BaseTree::bFloat) *(floatBranches.at(varposition).second)=value;
      else if (varbranchtype==BaseTree::bVectorInt) vectorIntBranches.at(varposition).second->push_back(value);
      else if (varbranchtype==BaseTree::bVectorFloat) vectorFloatBranches.at(varposition).second->push_back(value);
      else if (varbranchtype==BaseTree::bVectorDouble) vectorDoubleBranches.at(varposition).second->push_back(value);
    }
  }
  void* getBranchHandleRef(const std::string& branchname){
    int varposition=-1;
    BaseTree::BranchTypes varbranchtype = searchArray(branchname, varposition);
    if (varposition==-1 || varbranchtype==BaseTree::nBranchTypes){
      std::cerr << "Could not find the branch called " << branchname << "!" << std::endl;
      return 0;
    }
    else{
      if (varbranchtype==BaseTree::bInt) return intBranches.at(varposition).second;
      else if (varbranchtype==BaseTree::bFloat) return floatBranches.at(varposition).second;
      else if (varbranchtype==BaseTree::bVectorInt) return &(vectorIntBranches.at(varposition).second);
      else if (varbranchtype==BaseTree::bVectorFloat) return &(vectorFloatBranches.at(varposition).second);
      else if (varbranchtype==BaseTree::bVectorDouble) return &(vectorDoubleBranches.at(varposition).second);
      else return 0;
    }
  }
  void initializeBranches();
  void cleanBranches();
  void printEntry(int ev); // Mainly for debugging


protected:
  TTree* hvvtree;

  std::vector < std::pair<std::string, Int_t*> > intBranches;
  std::vector < std::pair<std::string, Float_t*> > floatBranches;
  std::vector < std::pair<std::string, vectorInt*> > vectorIntBranches;
  std::vector < std::pair<std::string, vectorFloat*> > vectorFloatBranches;
  std::vector < std::pair<std::string, vectorDouble*> > vectorDoubleBranches;
};

#endif
