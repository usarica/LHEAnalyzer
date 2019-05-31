#include "GlobalRecordTree.h"


bool GlobalRecordTree::reserveBranch(string branchname, BaseTree::BranchTypes branchtype, bool doSetAddress){
  bool isAvailable = true;
  if (options->isAnExcludedBranch(branchname)){
    isAvailable=false;
    if (doSetAddress && hvvtree->GetBranchStatus(branchname.c_str())) hvvtree->SetBranchStatus(branchname.c_str(), 0);
  }
  else if (
    (doSetAddress && !hvvtree->GetBranchStatus(branchname.c_str()))
    ||
    (!doSetAddress && hvvtree->GetBranchStatus(branchname.c_str()))
    ) isAvailable=false; // If setAddress to a non-existing branch or branch to an existing address
  if (isAvailable) bookBranch(branchname, branchtype, doSetAddress);
  return isAvailable;
}
void GlobalRecordTree::bookAllBranches(bool doSetAddress){
  if (!options){
    cerr << "GlobalRecordTree::bookAllBranches -> No options are set for the GlobalRecordTree!" << endl;
    return;
  }

  if (options->hasGlobalRecord()){
    globalRecordSet = options->getGlobalRecordSet();
    for (int r=0; r<globalRecordSet.size(); r++){
      BaseTree::BranchTypes interpretGlobalRecordType(globalRecordSet.at(r).second);
      // LEFT HERE

    }

    reserveBranch("genFinalState", BranchTypes::bInt, doSetAddress);

    bookPtEtaPhiMassIdBranches("Mother", BranchTypes::bVectorDouble, doSetAddress, true, true, true);
    bookPtEtaPhiMassIdBranches("H", BranchTypes::bFloat, doSetAddress, false, true, true);
    bookPtEtaPhiMassIdBranches("Z1", BranchTypes::bFloat, doSetAddress, false, false, true);
    bookPtEtaPhiMassIdBranches("Z2", BranchTypes::bFloat, doSetAddress, false, false, true);
    bookPtEtaPhiMassIdBranches("Za", BranchTypes::bFloat, doSetAddress, false, false, true);
    bookPtEtaPhiMassIdBranches("Zb", BranchTypes::bFloat, doSetAddress, false, false, true);

    bookPtEtaPhiMassIdBranches("AssociatedParticle", BranchTypes::bVectorDouble, doSetAddress, true, false, true);
    reserveBranch("GenDijetMass", BranchTypes::bFloat, doSetAddress);
    reserveBranch("GenDileptonMass", BranchTypes::bFloat, doSetAddress);
    reserveBranch("GenDijetVVMass", BranchTypes::bFloat, doSetAddress);
    reserveBranch("GenDileptonVVMass", BranchTypes::bFloat, doSetAddress);
    reserveBranch("GenDRjet", BranchTypes::bFloat, doSetAddress);
    reserveBranch("GenDRlepton", BranchTypes::bFloat, doSetAddress);

    reserveBranch("NGenAssociatedVs", BranchTypes::bInt, doSetAddress);
    bookPtEtaPhiMassIdBranches("AssociatedV", BranchTypes::bVectorDouble, doSetAddress, true, false, true);
    reserveBranch("GenAssociatedV_Particle1Index", BranchTypes::bVectorInt, doSetAddress);
    reserveBranch("GenAssociatedV_Particle2Index", BranchTypes::bVectorInt, doSetAddress);

    bookPtEtaPhiMassIdBranches("Lep1", BranchTypes::bFloat, doSetAddress, true, false, true);
    bookPtEtaPhiMassIdBranches("Lep2", BranchTypes::bFloat, doSetAddress, true, false, true);
    bookPtEtaPhiMassIdBranches("Lep3", BranchTypes::bFloat, doSetAddress, true, false, true);
    bookPtEtaPhiMassIdBranches("Lep4", BranchTypes::bFloat, doSetAddress, true, false, true);
  }
  if (options->processRecoInfo()){
    reserveBranch("isSelected", BranchTypes::bInt, doSetAddress);

    bookPtEtaPhiMassIdBranches("ZZ", BranchTypes::bFloat, doSetAddress, false, true, false);
    bookPtEtaPhiMassIdBranches("Z1", BranchTypes::bFloat, doSetAddress, false, false, false);
    bookPtEtaPhiMassIdBranches("Z2", BranchTypes::bFloat, doSetAddress, false, false, false);
    bookPtEtaPhiMassIdBranches("Za", BranchTypes::bFloat, doSetAddress, false, false, false);
    bookPtEtaPhiMassIdBranches("Zb", BranchTypes::bFloat, doSetAddress, false, false, false);

    bookPtEtaPhiMassIdBranches("AssociatedParticle", BranchTypes::bVectorDouble, doSetAddress, true, false, false);
    reserveBranch("DijetMass", BranchTypes::bFloat, doSetAddress);
    reserveBranch("DileptonMass", BranchTypes::bFloat, doSetAddress);
    reserveBranch("DijetVVMass", BranchTypes::bFloat, doSetAddress);
    reserveBranch("DileptonVVMass", BranchTypes::bFloat, doSetAddress);
    reserveBranch("DRjet", BranchTypes::bFloat, doSetAddress);
    reserveBranch("DRlepton", BranchTypes::bFloat, doSetAddress);

    reserveBranch("NAssociatedVs", BranchTypes::bInt, doSetAddress);
    bookPtEtaPhiMassIdBranches("AssociatedV", BranchTypes::bVectorDouble, doSetAddress, true, false, false);
    reserveBranch("AssociatedV_Particle1Index", BranchTypes::bVectorInt, doSetAddress);
    reserveBranch("AssociatedV_Particle2Index", BranchTypes::bVectorInt, doSetAddress);

    bookPtEtaPhiMassIdBranches("Lep1", BranchTypes::bFloat, doSetAddress, true, false, false);
    bookPtEtaPhiMassIdBranches("Lep2", BranchTypes::bFloat, doSetAddress, true, false, false);
    bookPtEtaPhiMassIdBranches("Lep3", BranchTypes::bFloat, doSetAddress, true, false, false);
    bookPtEtaPhiMassIdBranches("Lep4", BranchTypes::bFloat, doSetAddress, true, false, false);
  }
  bookAngularBranches(doSetAddress);
  if (options->initializeMELABranches() || doSetAddress) bookMELABranches(doSetAddress);
  actuateBranches(doSetAddress);
}


