#include "../interface/BaseTree.h"

BaseTree::BaseTree(string treename){ initTree(treename, ""); }
BaseTree::BaseTree(string treename, string treetitle){ initTree(treename, treetitle); }

bool BaseTree::bookBranch(string branchname, BaseTree::BranchTypes bType, bool doSetAddress){
  bool success=true;
  if (hvvtree!=0){
    if (bType==BaseTree::bInt){
      Int_t* container = 0;
      if (!doSetAddress || hvvtree->GetBranchStatus(branchname.c_str())){
        container = new Int_t; // Cannot be a null pointer even if setting branch address
        pair<string, Int_t*> varPair(branchname, container);
        intBranches.push_back(varPair);
      }
      else success=false;
    }
    else if (bType==BaseTree::bFloat){
      Float_t* container = 0;
      if (!doSetAddress || hvvtree->GetBranchStatus(branchname.c_str())){
        container = new Float_t; // Cannot be a null pointer even if setting branch address
        pair<string, Float_t*> varPair(branchname, container);
        floatBranches.push_back(varPair);
      }
      else success=false;
    }
    else if (bType==BaseTree::bVectorInt){
      vector<int>* container = 0;
      if (!doSetAddress) container = new vector<int>;
      if (!doSetAddress || hvvtree->GetBranchStatus(branchname.c_str())){
        pair<string, vector<int>*> varPair(branchname, container);
        vectorIntBranches.push_back(varPair);
      }
      else success=false;
    }
    else if (bType==BaseTree::bVectorDouble){
      vector<double>* container = 0;
      if (!doSetAddress) container = new vector<double>;
      if (!doSetAddress || hvvtree->GetBranchStatus(branchname.c_str())){
        pair<string, vector<double>*> varPair(branchname, container);
        vectorDoubleBranches.push_back(varPair);
      }
      else success=false;
    }
    success=true;
  }
  else success=false;
  if (!success) cerr << "BaseTree::bookBranch: Failed to book branch named " << branchname << "!" << endl;
  return success;
}
bool BaseTree::actuateBranches(bool doSetAddress){
  bool success=true;
  if (hvvtree!=0){
    if (!doSetAddress){
      for (int el=0; el<intBranches.size(); el++) hvvtree->Branch(intBranches.at(el).first.c_str(), intBranches.at(el).second);
      for (int el=0; el<floatBranches.size(); el++) hvvtree->Branch(floatBranches.at(el).first.c_str(), floatBranches.at(el).second);
      for (int el=0; el<vectorIntBranches.size(); el++) hvvtree->Branch(vectorIntBranches.at(el).first.c_str(), vectorIntBranches.at(el).second);
      for (int el=0; el<vectorDoubleBranches.size(); el++) hvvtree->Branch(vectorDoubleBranches.at(el).first.c_str(), vectorDoubleBranches.at(el).second);
    }
    else{
      for (int el=0; el<intBranches.size(); el++) hvvtree->SetBranchAddress(intBranches.at(el).first.c_str(), intBranches.at(el).second); // Already a pointer
      for (int el=0; el<floatBranches.size(); el++) hvvtree->SetBranchAddress(floatBranches.at(el).first.c_str(), floatBranches.at(el).second); // Already a pointer
      for (int el=0; el<vectorIntBranches.size(); el++) hvvtree->SetBranchAddress(vectorIntBranches.at(el).first.c_str(), &(vectorIntBranches.at(el).second)); // Need to pass the address of the pointer to std::vector
      for (int el=0; el<vectorDoubleBranches.size(); el++) hvvtree->SetBranchAddress(vectorDoubleBranches.at(el).first.c_str(), &(vectorDoubleBranches.at(el).second)); // Need to pass the address of the pointer to std::vector
    }
  }
  else success=false;
  if (!success) cerr << "BaseTree::actuateBranch: Failed to actuate the branches!" << endl;
  return success;
}
vector<string> BaseTree::getBranchList(){
  vector<string> branchlist;
  for (int el=0; el<intBranches.size(); el++) branchlist.push_back(intBranches.at(el).first);
  for (int el=0; el<floatBranches.size(); el++) branchlist.push_back(floatBranches.at(el).first);
  for (int el=0; el<vectorIntBranches.size(); el++) branchlist.push_back(vectorIntBranches.at(el).first);
  for (int el=0; el<vectorDoubleBranches.size(); el++) branchlist.push_back(vectorDoubleBranches.at(el).first);
  return branchlist;
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
    if (vectorIntBranches.at(el).second!=0){
      vectorIntBranches.at(el).second->clear();
      delete vectorIntBranches.at(el).second;
    }
    vectorIntBranches.at(el).second=0;
  }
  vectorIntBranches.clear();
  for (int el=0; el<vectorDoubleBranches.size(); el++){
    if (vectorDoubleBranches.at(el).second!=0){
      vectorDoubleBranches.at(el).second->clear();
      delete vectorDoubleBranches.at(el).second;
    }
    vectorDoubleBranches.at(el).second=0;
  }
  vectorDoubleBranches.clear();
}
void BaseTree::initializeBranches(){
  for (int el=0; el<intBranches.size(); el++){
    if (intBranches.at(el).second!=0) *(intBranches.at(el).second)=0;
  }
  for (int el=0; el<floatBranches.size(); el++){
    if (floatBranches.at(el).second!=0) *(floatBranches.at(el).second)=0;
  }
  for (int el=0; el<vectorIntBranches.size(); el++){
    if (vectorIntBranches.at(el).second!=0) vectorIntBranches.at(el).second->clear();
  }
  for (int el=0; el<vectorDoubleBranches.size(); el++){
    if (vectorDoubleBranches.at(el).second!=0) vectorDoubleBranches.at(el).second->clear();
  }
}




