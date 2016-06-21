#ifndef GLOBALRECORDTREE_H
#define GLOBALRECORDTREE_H

#include "BaseTree.h"
#include "OptionParser.h"

class GlobalRecordTree : public BaseTree{
public:
  GlobalRecordTree():BaseTree(), options(0){}
  GlobalRecordTree(string treename) : BaseTree(treename), options(0){}
  GlobalRecordTree(string treename, string treetitle) : BaseTree(treename, treetitle), options(0){}
  GlobalRecordTree(string treename, TFile* fin) : BaseTree(treename, fin), options(0){}

  void setOptions(OptionParser* options_){ options=options_; }

  bool reserveBranch(string branchname, BaseTree::BranchTypes branchtype, bool doSetAddress);
  void bookAllBranches(bool doSetAddress);

protected:
  void getGlobalRecordSet();
  BaseTree::BranchTypes interpretGlobalRecordType(string strRcdVal);


  OptionParser* options;
  vector < pair<string,string> > globalRecordSet;
};

#endif
