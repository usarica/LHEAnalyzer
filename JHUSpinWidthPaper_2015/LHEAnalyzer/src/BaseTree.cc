#include "../interface/BaseTree.h"

BaseTree::BaseTree(string treename){ initTree(treename, ""); }
BaseTree::BaseTree(string treename, string treetitle){ initTree(treename, treetitle); }

bool BaseTree::bookBranch(string branchname, BaseTree::BranchTypes bType){
  bool success=1;
  if (hvvtree!=0){
    if (bType==BaseTree::bInt){
      Int_t* container = new Int_t;
      pair<string, Int_t*> varPair(branchname, container);
      intBranches.push_back(varPair);
      hvvtree->Branch(varPair.first.c_str(), varPair.second);
    }
    else if (bType==BaseTree::bFloat){
      Float_t* container = new Float_t;
      pair<string, Float_t*> varPair(branchname, container);
      floatBranches.push_back(varPair);
      hvvtree->Branch(varPair.first.c_str(), varPair.second);
    }
    else if (bType==BaseTree::bVectorInt){
      vector<int>* container = new vector<int>;
      pair<string, vector<int>*> varPair(branchname, container);
      vectorIntBranches.push_back(varPair);
      hvvtree->Branch(varPair.first.c_str(), varPair.second);
    }
    else if (bType==BaseTree::bVectorDouble){
      vector<double>* container = new vector<double>;
      pair<string, vector<double>*> varPair(branchname, container);
      vectorDoubleBranches.push_back(varPair);
      hvvtree->Branch(varPair.first.c_str(), varPair.second);
    }
    success=true;
  }
  else success=false;
  if (!success) cerr << "BaseTree::bookBranch: Failed to create branch named " << branchname << "!" << endl;
  return success;
}

BaseTree::BranchTypes BaseTree::searchArray(string branchname, int& position){
  for (int el=0; el<intBranches.size(); el++){
    if (branchname==intBranches.at(el).first){
      position = el;
      return BranchTypes::bInt;
    }
  }
  for (int el=0; el<floatBranches.size(); el++){
    if (branchname==floatBranches.at(el).first){
      position = el;
      return BranchTypes::bFloat;
    }
  }
  for (int el=0; el<vectorIntBranches.size(); el++){
    if (branchname==vectorIntBranches.at(el).first){
      position = el;
      return BranchTypes::bVectorInt;
    }
  }
  for (int el=0; el<vectorDoubleBranches.size(); el++){
    if (branchname==vectorDoubleBranches.at(el).first){
      position = el;
      return BranchTypes::bVectorDouble;
    }
  }
  return BaseTree::nBranchTypes;
}
void BaseTree::cleanBranches(){
  for (int el=0; el<intBranches.size(); el++){
    if (intBranches.at(el).second!=0) delete intBranches.at(el).second;
    intBranches.at(el).second=0;
  }
  intBranches.clear();
  for (int el=0; el<floatBranches.size(); el++){
    if (floatBranches.at(el).second!=0) delete floatBranches.at(el).second;
    floatBranches.at(el).second=0;
  }
  floatBranches.clear();
  for (int el=0; el<vectorIntBranches.size(); el++){
    vectorIntBranches.at(el).second->clear();
    if (vectorIntBranches.at(el).second!=0) delete vectorIntBranches.at(el).second;
    vectorIntBranches.at(el).second=0;
  }
  vectorIntBranches.clear();
  for (int el=0; el<vectorDoubleBranches.size(); el++){
    vectorDoubleBranches.at(el).second->clear();
    if (vectorDoubleBranches.at(el).second!=0) delete vectorDoubleBranches.at(el).second;
    vectorDoubleBranches.at(el).second=0;
  }
  vectorDoubleBranches.clear();
}
void BaseTree::initializeBranches(){
  for (int el=0; el<intBranches.size(); el++){
    *(intBranches.at(el).second)=0;
  }
  for (int el=0; el<floatBranches.size(); el++){
    *(floatBranches.at(el).second)=0;
  }
  for (int el=0; el<vectorIntBranches.size(); el++){
    vectorIntBranches.at(el).second->clear();
  }
  for (int el=0; el<vectorDoubleBranches.size(); el++){
    vectorDoubleBranches.at(el).second->clear();
  }
}




