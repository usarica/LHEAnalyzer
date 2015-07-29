#include "../interface/HVVTree.h"


bool HVVTree::reserveBranch(string branchname, BaseTree::BranchTypes branchtype, bool doSetAddress){
  bool isAvailable = true;
  if (options->isAnExcludedBranch(branchname)){
    isAvailable=false;
    if (doSetAddress && hvvtree->GetBranchStatus(branchname.c_str())) hvvtree->SetBranchStatus(branchname.c_str(), 0);
  }
  else if (doSetAddress && !hvvtree->GetBranchStatus(branchname.c_str())) isAvailable=false;
  if (isAvailable) bookBranch(branchname, branchtype, doSetAddress);
  return isAvailable;
}
void HVVTree::bookAllBranches(bool doSetAddress){
  if (!options){
    cerr << "HVVTree::bookAllBranches -> No options are set for the HVVTree!" << endl;
    return;
  }

  reserveBranch("MC_weight", BranchTypes::bFloat, doSetAddress);
  if (options->processGenInfo()){
    reserveBranch("genFinalState", BranchTypes::bInt, doSetAddress);

    reserveBranch("GenMotherMass", BranchTypes::bVectorDouble, doSetAddress);
    reserveBranch("GenMotherPt", BranchTypes::bVectorDouble, doSetAddress);
    reserveBranch("GenMotherPz", BranchTypes::bVectorDouble, doSetAddress);
    reserveBranch("GenMotherPhi", BranchTypes::bVectorDouble, doSetAddress);
    reserveBranch("GenMotherId", BranchTypes::bVectorInt, doSetAddress);

    reserveBranch("GenHMass", BranchTypes::bFloat, doSetAddress);
    reserveBranch("GenHPt", BranchTypes::bFloat, doSetAddress);
    reserveBranch("GenHPz", BranchTypes::bFloat, doSetAddress);
    reserveBranch("GenHPhi", BranchTypes::bFloat, doSetAddress);

    reserveBranch("GenZ1Mass", BranchTypes::bFloat, doSetAddress);
    reserveBranch("GenZ1Pt", BranchTypes::bFloat, doSetAddress);
    reserveBranch("GenZ1Phi", BranchTypes::bFloat, doSetAddress);
    reserveBranch("GenZ1Eta", BranchTypes::bFloat, doSetAddress);

    reserveBranch("GenZ2Mass", BranchTypes::bFloat, doSetAddress);
    reserveBranch("GenZ2Pt", BranchTypes::bFloat, doSetAddress);
    reserveBranch("GenZ2Phi", BranchTypes::bFloat, doSetAddress);
    reserveBranch("GenZ2Eta", BranchTypes::bFloat, doSetAddress);

    reserveBranch("GenZaMass", BranchTypes::bFloat, doSetAddress);
    reserveBranch("GenZaPt", BranchTypes::bFloat, doSetAddress);
    reserveBranch("GenZaPhi", BranchTypes::bFloat, doSetAddress);
    reserveBranch("GenZaEta", BranchTypes::bFloat, doSetAddress);

    reserveBranch("GenZbMass", BranchTypes::bFloat, doSetAddress);
    reserveBranch("GenZbPt", BranchTypes::bFloat, doSetAddress);
    reserveBranch("GenZbPhi", BranchTypes::bFloat, doSetAddress);
    reserveBranch("GenZbEta", BranchTypes::bFloat, doSetAddress);

    reserveBranch("GenAssociatedParticleMass", BranchTypes::bVectorDouble, doSetAddress);
    reserveBranch("GenAssociatedParticlePt", BranchTypes::bVectorDouble, doSetAddress);
    reserveBranch("GenAssociatedParticleEta", BranchTypes::bVectorDouble, doSetAddress);
    reserveBranch("GenAssociatedParticlePhi", BranchTypes::bVectorDouble, doSetAddress);
    reserveBranch("GenAssociatedParticleId", BranchTypes::bVectorInt, doSetAddress);
    reserveBranch("GenDijetMass", BranchTypes::bFloat, doSetAddress);

    reserveBranch("GenNAssociatedVs", BranchTypes::bInt, doSetAddress);
    reserveBranch("GenAssociatedVMass", BranchTypes::bVectorDouble, doSetAddress);
    reserveBranch("GenAssociatedVPt", BranchTypes::bVectorDouble, doSetAddress);
    reserveBranch("GenAssociatedVEta", BranchTypes::bVectorDouble, doSetAddress);
    reserveBranch("GenAssociatedVPhi", BranchTypes::bVectorDouble, doSetAddress);
    reserveBranch("GenAssociatedVId", BranchTypes::bVectorInt, doSetAddress);
    reserveBranch("GenAssociatedV_Particle1Index", BranchTypes::bVectorInt, doSetAddress);
    reserveBranch("GenAssociatedV_Particle2Index", BranchTypes::bVectorInt, doSetAddress);

    reserveBranch("GenhelcosthetaZ1", BranchTypes::bFloat, doSetAddress);
    reserveBranch("GenhelcosthetaZ2", BranchTypes::bFloat, doSetAddress);
    reserveBranch("Genhelphi", BranchTypes::bFloat, doSetAddress);
    reserveBranch("Gencosthetastar", BranchTypes::bFloat, doSetAddress);
    reserveBranch("GenphistarZ1", BranchTypes::bFloat, doSetAddress);

    reserveBranch("GenLep1Mass", BranchTypes::bFloat, doSetAddress);
    reserveBranch("GenLep2Mass", BranchTypes::bFloat, doSetAddress);
    reserveBranch("GenLep3Mass", BranchTypes::bFloat, doSetAddress);
    reserveBranch("GenLep4Mass", BranchTypes::bFloat, doSetAddress);
    reserveBranch("GenLep1Pt", BranchTypes::bFloat, doSetAddress);
    reserveBranch("GenLep2Pt", BranchTypes::bFloat, doSetAddress);
    reserveBranch("GenLep3Pt", BranchTypes::bFloat, doSetAddress);
    reserveBranch("GenLep4Pt", BranchTypes::bFloat, doSetAddress);
    reserveBranch("GenLep1Eta", BranchTypes::bFloat, doSetAddress);
    reserveBranch("GenLep2Eta", BranchTypes::bFloat, doSetAddress);
    reserveBranch("GenLep3Eta", BranchTypes::bFloat, doSetAddress);
    reserveBranch("GenLep4Eta", BranchTypes::bFloat, doSetAddress);
    reserveBranch("GenLep1Phi", BranchTypes::bFloat, doSetAddress);
    reserveBranch("GenLep2Phi", BranchTypes::bFloat, doSetAddress);
    reserveBranch("GenLep3Phi", BranchTypes::bFloat, doSetAddress);
    reserveBranch("GenLep4Phi", BranchTypes::bFloat, doSetAddress);
    reserveBranch("GenLep1Id", BranchTypes::bInt, doSetAddress);
    reserveBranch("GenLep2Id", BranchTypes::bInt, doSetAddress);
    reserveBranch("GenLep3Id", BranchTypes::bInt, doSetAddress);
    reserveBranch("GenLep4Id", BranchTypes::bInt, doSetAddress);
  }
  if (options->processRecoInfo()){
    reserveBranch("isSelected", BranchTypes::bInt, doSetAddress);

    reserveBranch("ZZMass", BranchTypes::bFloat, doSetAddress);
    reserveBranch("ZZPt", BranchTypes::bFloat, doSetAddress);
    reserveBranch("ZZPz", BranchTypes::bFloat, doSetAddress);
    reserveBranch("ZZPhi", BranchTypes::bFloat, doSetAddress);

    reserveBranch("Z1Mass", BranchTypes::bFloat, doSetAddress);
    reserveBranch("Z1Pt", BranchTypes::bFloat, doSetAddress);
    reserveBranch("Z1Phi", BranchTypes::bFloat, doSetAddress);
    reserveBranch("Z1Eta", BranchTypes::bFloat, doSetAddress);

    reserveBranch("Z2Mass", BranchTypes::bFloat, doSetAddress);
    reserveBranch("Z2Pt", BranchTypes::bFloat, doSetAddress);
    reserveBranch("Z2Phi", BranchTypes::bFloat, doSetAddress);
    reserveBranch("Z2Eta", BranchTypes::bFloat, doSetAddress);

    reserveBranch("ZaMass", BranchTypes::bFloat, doSetAddress);
    reserveBranch("ZaPt", BranchTypes::bFloat, doSetAddress);
    reserveBranch("ZaPhi", BranchTypes::bFloat, doSetAddress);
    reserveBranch("ZaEta", BranchTypes::bFloat, doSetAddress);

    reserveBranch("ZbMass", BranchTypes::bFloat, doSetAddress);
    reserveBranch("ZbPt", BranchTypes::bFloat, doSetAddress);
    reserveBranch("ZbPhi", BranchTypes::bFloat, doSetAddress);
    reserveBranch("ZbEta", BranchTypes::bFloat, doSetAddress);

    reserveBranch("AssociatedParticleMass", BranchTypes::bVectorDouble, doSetAddress);
    reserveBranch("AssociatedParticlePt", BranchTypes::bVectorDouble, doSetAddress);
    reserveBranch("AssociatedParticleEta", BranchTypes::bVectorDouble, doSetAddress);
    reserveBranch("AssociatedParticlePhi", BranchTypes::bVectorDouble, doSetAddress);
    reserveBranch("AssociatedParticleId", BranchTypes::bVectorInt, doSetAddress);
    reserveBranch("DijetMass", BranchTypes::bFloat, doSetAddress);

    reserveBranch("NAssociatedVs", BranchTypes::bInt, doSetAddress);
    reserveBranch("AssociatedVMass", BranchTypes::bVectorDouble, doSetAddress);
    reserveBranch("AssociatedVPt", BranchTypes::bVectorDouble, doSetAddress);
    reserveBranch("AssociatedVEta", BranchTypes::bVectorDouble, doSetAddress);
    reserveBranch("AssociatedVPhi", BranchTypes::bVectorDouble, doSetAddress);
    reserveBranch("AssociatedVId", BranchTypes::bVectorInt, doSetAddress);
    reserveBranch("AssociatedV_Particle1Index", BranchTypes::bVectorInt, doSetAddress);
    reserveBranch("AssociatedV_Particle2Index", BranchTypes::bVectorInt, doSetAddress);

    reserveBranch("helcosthetaZ1", BranchTypes::bFloat, doSetAddress);
    reserveBranch("helcosthetaZ2", BranchTypes::bFloat, doSetAddress);
    reserveBranch("helphi", BranchTypes::bFloat, doSetAddress);
    reserveBranch("costhetastar", BranchTypes::bFloat, doSetAddress);
    reserveBranch("phistarZ1", BranchTypes::bFloat, doSetAddress);

    reserveBranch("Lep1Mass", BranchTypes::bFloat, doSetAddress);
    reserveBranch("Lep2Mass", BranchTypes::bFloat, doSetAddress);
    reserveBranch("Lep3Mass", BranchTypes::bFloat, doSetAddress);
    reserveBranch("Lep4Mass", BranchTypes::bFloat, doSetAddress);
    reserveBranch("Lep1Pt", BranchTypes::bFloat, doSetAddress);
    reserveBranch("Lep2Pt", BranchTypes::bFloat, doSetAddress);
    reserveBranch("Lep3Pt", BranchTypes::bFloat, doSetAddress);
    reserveBranch("Lep4Pt", BranchTypes::bFloat, doSetAddress);
    reserveBranch("Lep1Eta", BranchTypes::bFloat, doSetAddress);
    reserveBranch("Lep2Eta", BranchTypes::bFloat, doSetAddress);
    reserveBranch("Lep3Eta", BranchTypes::bFloat, doSetAddress);
    reserveBranch("Lep4Eta", BranchTypes::bFloat, doSetAddress);
    reserveBranch("Lep1Phi", BranchTypes::bFloat, doSetAddress);
    reserveBranch("Lep2Phi", BranchTypes::bFloat, doSetAddress);
    reserveBranch("Lep3Phi", BranchTypes::bFloat, doSetAddress);
    reserveBranch("Lep4Phi", BranchTypes::bFloat, doSetAddress);
    reserveBranch("Lep1Id", BranchTypes::bInt, doSetAddress);
    reserveBranch("Lep2Id", BranchTypes::bInt, doSetAddress);
    reserveBranch("Lep3Id", BranchTypes::bInt, doSetAddress);
    reserveBranch("Lep4Id", BranchTypes::bInt, doSetAddress);
  }
  if (options->initializeMELA()) setMELABranches(doSetAddress);
  actuateBranches(doSetAddress);
}

void HVVTree::constructMELABranchList(){
  vector<string> blist;

  // Add gen. prod. MEs
  TVar::Production prod = options->getSampleProductionId().first;
  TVar::MatrixElement me = options->getSampleProductionId().second;
  setupMELASignalMECases(blist, prod, me, true, true);

  cout << "constructMELABranchList: Gen prod MES done" << endl;

  // Reco. prod. MEs
  me = TVar::JHUGen;
  prod = TVar::ZZGG; setupMELASignalMECases(blist, prod, me, false, true);
  prod = TVar::JJGG; setupMELASignalMECases(blist, prod, me, false, true);
  prod = TVar::JJVBF; setupMELASignalMECases(blist, prod, me, false, true);
  prod = TVar::ZH; setupMELASignalMECases(blist, prod, me, false, true);
  prod = TVar::WH; setupMELASignalMECases(blist, prod, me, false, true);
  prod = TVar::JH; setupMELASignalMECases(blist, prod, me, false, true);
  //prod = TVar::ttH; setupMELASignalMECases(blist, prod, me, false, true);
  //prod = TVar::bbH; setupMELASignalMECases(blist, prod, me, false, true);

  cout << "constructMELABranchList: Reco prod MES done" << endl;

  // Gen. decay MEs
  me = TVar::JHUGen;
  prod = TVar::ZZGG; setupMELASignalMECases(blist, prod, me, true, false);
  me = TVar::MCFM;
  prod = TVar::ZZGG; setupMELASignalMECases(blist, prod, me, true, false);

  cout << "constructMELABranchList: Gen dec MES done" << endl;

  // Reco. decay MEs
  me = TVar::JHUGen;
  prod = TVar::ZZGG; setupMELASignalMECases(blist, prod, me, false, false);
  me = TVar::MCFM;
  prod = TVar::ZZGG; setupMELASignalMECases(blist, prod, me, false, false);

  cout << "constructMELABranchList: Reco dec MES done" << endl;

  // ggzz_VAMCFM
  blist.push_back("ggzz_VAMCFM");
  // VVzz_VAMCFM
  //blist.push_back("VVzz_VAMCFM");
  // bkg_VAMCFM
  blist.push_back("bkg_VAMCFM");
  // ggzz_VAMCFM
  blist.push_back("Gen_ggzz_VAMCFM");
  // VVzz_VAMCFM
  //blist.push_back("Gen_VVzz_VAMCFM");
  // bkg_VAMCFM
  blist.push_back("Gen_bkg_VAMCFM");

  cout << "constructMELABranchList: ggZZ MES done" << endl;

  // P(m4l)
  vector<string> ME_m4l;
  ME_m4l.push_back("p0plus_m4l");
  ME_m4l.push_back("p0plus_m4l_ScaleUp");
  ME_m4l.push_back("p0plus_m4l_ScaleDown");
  ME_m4l.push_back("p0plus_m4l_ResUp");
  ME_m4l.push_back("p0plus_m4l_ResDown");
  ME_m4l.push_back("bkg_m4l");
  ME_m4l.push_back("bkg_m4l_ScaleUp");
  ME_m4l.push_back("bkg_m4l_ScaleDown");
  ME_m4l.push_back("bkg_m4l_ResUp");
  ME_m4l.push_back("bkg_m4l_ResDown");
  string chvar = "m4l";
  if (options->hasRecoDecayME(chvar) || options->hasRecoDecayME("All") || options->hasRecoDecayME("all")){
    for (int b=0; b<ME_m4l.size(); b++) blist.push_back(ME_m4l.at(b));
  }
  cout << "constructMELABranchList: m4l MES done" << endl;

  for (int b=0; b<blist.size(); b++) melaProbBranches.push_back(blist.at(b));
}
void HVVTree::setupMELASignalMECases(vector<string>& accumulatedlist, TVar::Production prod, TVar::MatrixElement me, bool isGen, bool isProdME){
  vector<string> gList;
  gList.push_back("g1");
  gList.push_back("g2");
  gList.push_back("g4");
  if (!isProdME || prod==TVar::ZH || prod==TVar::WH) gList.push_back("g1_prime2");
  int sgList = gList.size();
  int** gCount = new int*[sgList];
  for (int gg=0; gg<sgList; gg++){
    gCount[gg] = new int[2];
    for (int im=0; im<2; im++) gCount[gg][im] = 0;
  }
  vector<int> v_gCount[2];

  // Check for any invalid productions
  bool invalidProduction = false;
  if (me==TVar::MCFM && !(prod==TVar::ZZGG/* || prod==TVar::JJVBF || prod==TVar::WH || prod==TVar::ZH*/)) invalidProduction=true;
  if (isProdME && me==TVar::MCFM && prod==TVar::ZZGG) invalidProduction=true;
  bool overridden = false;
  if (
    ((options->hasGenProdME("None") || options->hasGenProdME("none")) && isGen && isProdME)
    ||
    ((options->hasRecoProdME("None") || options->hasRecoProdME("none")) && !isGen && isProdME)
    ||
    ((options->hasGenDecayME("None") || options->hasGenDecayME("none")) && isGen && !isProdME)
    ||
    ((options->hasRecoDecayME("None") || options->hasRecoDecayME("none")) && !isGen && !isProdME)
    ) overridden=true;

  // Count number of gi occurences
  bool noInstance=true;
  for (int gg=0; gg<sgList; gg++){
    for (int im=0; im<2; im++){
      if (invalidProduction || overridden) break; // Overriding options and the invalid MEs
      if ((prod==TVar::JJGG || prod==TVar::JH || prod==TVar::ttH || prod==TVar::bbH) && gg==0) continue; // Unphysical ME, g1 does not exist for these
      if (prod==TVar::JH && gg!=0 && im!=0) continue; // Unimplemented ME
      string chvar = gList.at(gg);
      if (im==1) chvar.append("_pi2");
      if (
        ((options->hasGenProdME(chvar) || options->hasGenProdME("All") || options->hasGenProdME("all")) && isGen && isProdME)
        ||
        ((options->hasRecoProdME(chvar) || options->hasRecoProdME("All") || options->hasRecoProdME("all")) && !isGen && isProdME)
        ||
        ((options->hasGenDecayME(chvar) || options->hasGenDecayME("All") || options->hasGenDecayME("all")) && isGen && !isProdME)
        ||
        ((options->hasRecoDecayME(chvar) || options->hasRecoDecayME("All") || options->hasRecoDecayME("all")) && !isGen && !isProdME)
        ){
        gCount[gg][im]++;
        if (noInstance) noInstance=false;
      }
    }
  }
  for (int gg=0; gg<sgList; gg++){
    for (int im=0; im<2; im++) v_gCount[im].push_back(gCount[gg][im]);
  }
  if (!noInstance){
    vector<string> blist = getMELASignalMEBranches(prod, me, gList, v_gCount[0], v_gCount[1], isGen, isProdME);
    for (int b=0; b<blist.size(); b++) accumulatedlist.push_back(blist.at(b));
  }

  for (int gg=0; gg<sgList; gg++) delete[] gCount[gg];
  delete[] gCount;
}
vector<string> HVVTree::getMELASignalMEBranches(TVar::Production prod, TVar::MatrixElement me, vector<string> gList, vector<int> gCountRe, vector<int> gCountIm, bool isGen, bool isProdME){
  vector<string> blist;
  string strcore = "";

  if (prod==TVar::JJGG) strcore = "hjj_";
  else if (prod==TVar::JH) strcore = "hj_";
  else if (prod==TVar::ttH) strcore = "tth_";
  else if (prod==TVar::bbH) strcore = "bbh_";
  else if (prod==TVar::JJVBF) strcore = "vbf_";
  else if (prod==TVar::ZH) strcore = "zh_";
  else if (prod==TVar::WH) strcore = "wh_";
  if (isProdME && prod==TVar::ZZGG) strcore.insert(0, "prod_");
  if (isGen) strcore.insert(0, "Gen_");

  int sgList = gList.size();
  int** gCount = new int*[sgList];
  for (int gg=0; gg<sgList; gg++){
    gCount[gg] = new int[2];
    for (int im=0; im<2; im++) gCount[gg][im] = 0;
  }

  bool invalidProduction = false;
  if (me==TVar::MCFM && !(prod==TVar::ZZGG/* || prod==TVar::JJVBF || prod==TVar::WH || prod==TVar::ZH*/)) invalidProduction=true;
  if (isProdME && me==TVar::MCFM && prod==TVar::ZZGG) invalidProduction=true;

  // Count number of gi occurences
  for (int gg=0; gg<sgList; gg++){
    for (int im=0; im<2; im++){
      if (invalidProduction) break;
      if ((prod==TVar::JJGG || prod==TVar::JH || prod==TVar::ttH || prod==TVar::bbH) && gList.at(gg)=="g1") continue; // Unphysical ME, g1 does not exist for these
      if (prod==TVar::JH && !(gList.at(gg)=="g2" && im==0)) continue; // Unimplemented ME

      if (im==0) gCount[gg][im]=gCountRe.at(gg);
      else gCount[gg][im]=gCountIm.at(gg);
    }
  }
  for (int ai1=0; ai1<sgList; ai1++){
    if (gCount[ai1][0]==0 && gCount[ai1][1]==0) continue;
    for (int ai2=ai1; ai2<sgList; ai2++){
      if (gCount[ai2][0]==0 && gCount[ai2][1]==0) continue;
      if (ai1==ai2){ // ai1 MEs
        string tmpVarType;
        if (
          gList.at(ai1)=="g1"
          ||
          ((prod==TVar::JJGG || prod==TVar::JH || prod==TVar::ttH || prod==TVar::bbH) && gList.at(ai1)=="g2")
          ) tmpVarType = "0plus";
        else if (gList.at(ai1)=="g2") tmpVarType = "0hplus";
        else if (gList.at(ai1)=="g4") tmpVarType = "0minus";
        else if (gList.at(ai1)=="g1_prime2") tmpVarType = "0_p1prime2";
        tmpVarType.insert(0, "p");

        string strlist;

        strlist = tmpVarType;
        if (me==TVar::MCFM) strlist.append("_VAMCFM");
        else strlist.append("_VAJHU");
        strlist.insert(0, strcore);
        blist.push_back(strlist); // Sig

        if (me==TVar::MCFM){
          if (gCount[ai1][0]>0){
            strlist = tmpVarType;
            if (prod==TVar::ZZGG) strlist.insert(0, "ggzz_");
            else strlist.insert(0, "VVzz_");
            strlist.insert(0, strcore);
            strlist.append("_VAMCFM");
            blist.push_back(strlist); // BSI, re Sig component
          }
          if (gCount[ai1][1]>0){
            strlist = tmpVarType;
            if (prod==TVar::ZZGG) strlist.insert(0, "ggzz_");
            else strlist.insert(0, "VVzz_");
            strlist.insert(0, strcore);
            strlist.append("_pi2");
            strlist.append("_VAMCFM");
            blist.push_back(strlist); // BSI, im Sig component
          }
        }
      }
      else{ // ai1 - ai2 MEs
        string tmpVar[2] ={ gList.at(ai1), gList.at(ai2) };
        string strlist;
        for (int im1=0; im1<2; im1++){
          for (int im2=0; im2<2; im2++){
            if (gCount[ai1][im1]==0 || gCount[ai2][im2]==0) continue;
            string tmpVarType[2] ={ tmpVar[0], tmpVar[1] };

            if (im1==1) tmpVarType[0].append("_pi2_");
            if (im2==1) tmpVarType[1].append("_pi2");

            strlist = tmpVarType[0]+tmpVarType[1];
            if (me==TVar::MCFM) strlist.append("_VAMCFM");
            else strlist.append("_VAJHU");
            strlist.insert(0, "p");
            strlist.insert(0, strcore);
            if (
              !(im1==im2 && im1==1) // g2_pi2_g4_pi2 == g2g4 when no bkg is involved.
              &&
              !(im1>im2 && gCount[ai1][1-im1]==1 && gCount[ai2][1-im2]==1) // Take the form g2g4_pi2 etc., not g2_pi2_g4
              ) blist.push_back(strlist); // Sig, only

            if (me==TVar::MCFM){
              if (prod==TVar::ZZGG) strlist.insert(0, "ggzz_");
              else strlist.insert(0, "VVzz_");
              blist.push_back(strlist); // BSI, take all forms since bkg. sets the absolite phase.
            }
          }
        }
      }
    }
  }

  for (int gg=0; gg<sgList; gg++) delete[] gCount[gg];
  delete[] gCount;
  return blist;
}
void HVVTree::setMELABranches(bool doSetAddress){
  cout << "Starting MELA branches" << endl;
  constructMELABranchList();
  cout << "SConstructed..." << endl;
  for (int b=0; b<melaProbBranches.size(); b++){
    cout << melaProbBranches.at(b) << endl;
    reserveBranch(melaProbBranches.at(b), BranchTypes::bFloat, doSetAddress);
  }
}


void HVVTree::fillMotherInfo(Particle* mother){
  if (options!=0 && options->processGenInfo() && mother!=0){
    setVal("GenMotherMass", mother->m());
    setVal("GenMotherPt", mother->pt());
    setVal("GenMotherPz", mother->z());
    setVal("GenMotherPhi", mother->phi());
    setVal("GenMotherId", mother->id);
  }
}


void HVVTree::fillCandidate(ZZCandidate* pH, bool isGen){
  if (!options) return;
  if ((!options->processGenInfo() && isGen) || (!options->processRecoInfo() && !isGen)) return;

  string varname;
  string strcore = "ZZ";
  if (isGen) strcore = "GenH";

  varname = strcore + "Mass"; setVal(varname, (pH!=0 ? pH->m() : 0.));
  varname = strcore + "Pt"; setVal(varname, (pH!=0 ? pH->pt() : 0));
  varname = strcore + "Pz"; setVal(varname, (pH!=0 ? pH->z() : 0));
  varname = strcore + "Phi"; setVal(varname, (pH!=0 ? pH->phi() : 0));

  fillCandidateDaughters(pH, isGen);
  fillDaughterProducts(pH, isGen);
  fillAssociatedInfo(pH, isGen);
  fillDecayAngles(pH, isGen);
//  fillProductionAngles(pH, isGen);
}
void HVVTree::fillCandidateDaughters(ZZCandidate* pH, bool isGen){
  string varname;
  string strcore;

  Particle* pV1=(pH!=0 ? pH->getSortedV(0) : 0);
  Particle* pV2=(pH!=0 ? pH->getSortedV(1) : 0);

  strcore = "Z1";
  if (isGen) strcore.insert(0, "Gen");
  varname = strcore + "Mass"; setVal(varname, (pV1!=0 ? pV1->m() : 0));
  varname = strcore + "Pt"; setVal(varname, (pV1!=0 ? pV1->pt() : 0));
  varname = strcore + "Eta"; setVal(varname, (pV1!=0 ? pV1->eta() : 0));
  varname = strcore + "Phi"; setVal(varname, (pV1!=0 ? pV1->phi() : 0));
  strcore = "Z2";
  if (isGen) strcore.insert(0, "Gen");
  varname = strcore + "Mass"; setVal(varname, (pV2!=0 ? pV2->m() : 0));
  varname = strcore + "Pt"; setVal(varname, (pV2!=0 ? pV2->pt() : 0));
  varname = strcore + "Eta"; setVal(varname, (pV2!=0 ? pV2->eta() : 0));
  varname = strcore + "Phi"; setVal(varname, (pV2!=0 ? pV2->phi() : 0));

  TLorentzVector nullVector(0, 0, 0, 0);
  TLorentzVector pZ1alt = (pH!=0 ? pH->getAlternativeVMomentum(0) : nullVector);
  TLorentzVector pZ2alt = (pH!=0 ? pH->getAlternativeVMomentum(1) : nullVector);

  strcore = "Za";
  if (isGen) strcore.insert(0, "Gen");
  varname = strcore + "Mass"; setVal(varname, (pH!=0 ? pZ1alt.M() : 0));
  varname = strcore + "Pt"; setVal(varname, (pH!=0 ? pZ1alt.Pt() : 0));
  varname = strcore + "Eta"; setVal(varname, (pH!=0 ? pZ1alt.Eta() : 0));
  varname = strcore + "Phi"; setVal(varname, (pH!=0 ? pZ1alt.Phi() : 0));
  strcore = "Zb";
  if (isGen) strcore.insert(0, "Gen");
  varname = strcore + "Mass"; setVal(varname, (pH!=0 ? pZ2alt.M() : 0));
  varname = strcore + "Pt"; setVal(varname, (pH!=0 ? pZ2alt.Pt() : 0));
  varname = strcore + "Eta"; setVal(varname, (pH!=0 ? pZ2alt.Eta() : 0));
  varname = strcore + "Phi"; setVal(varname, (pH!=0 ? pZ2alt.Phi() : 0));
}
void HVVTree::fillDaughterProducts(ZZCandidate* pH, bool isGen){
  string varname;
  string strcore = "Lep";
  if (isGen) strcore.insert(0, "Gen");

  for (int v=0; v<2; v++){
    for (int d=0; d<2; d++){
      int iLep = 2*v+d+1;
      char cILep[2];
      sprintf(cILep, "%i", iLep);
      string strILep = string(cILep);

      Particle* lepton = (pH!=0 ? pH->getSortedV(v)->getDaughter(d) : 0);

      varname = strcore + strILep + "Mass"; setVal(varname, (lepton!=0 ? lepton->m() : 0));
      varname = strcore + strILep + "Pt"; setVal(varname, (lepton!=0 ? lepton->pt() : 0));
      varname = strcore + strILep + "Eta"; setVal(varname, (lepton!=0 ? lepton->eta() : 0));
      varname = strcore + strILep + "Phi"; setVal(varname, (lepton!=0 ? lepton->phi() : 0));
      varname = strcore + strILep + "Id"; setVal(varname, (lepton!=0 ? lepton->id : 0));
    }
  }

  if (isGen){
    Int_t genFinalState=-1;
    if (pH!=0){
      if (PDGHelpers::isAZBoson(pH->getSortedV(0)->id) && PDGHelpers::isAZBoson(pH->getSortedV(1)->id)){
        // 4l
        if (std::abs(pH->getSortedV(0)->getDaughter(0)->id)==13 && std::abs(pH->getSortedV(1)->getDaughter(0)->id)==13) genFinalState=0; // 4mu
        else if (std::abs(pH->getSortedV(0)->getDaughter(0)->id)==11 && std::abs(pH->getSortedV(1)->getDaughter(0)->id)==11) genFinalState=1; // 4e
        else if ((std::abs(pH->getSortedV(0)->getDaughter(0)->id)==11 && std::abs(pH->getSortedV(1)->getDaughter(0)->id)==13) || (std::abs(pH->getSortedV(0)->getDaughter(0)->id)==13 && std::abs(pH->getSortedV(1)->getDaughter(0)->id)==11)) genFinalState=2; // 2e2mu
        else if (std::abs(pH->getSortedV(0)->getDaughter(0)->id)==15 && std::abs(pH->getSortedV(1)->getDaughter(0)->id)==15) genFinalState=3; // 4tau
        else if ((std::abs(pH->getSortedV(0)->getDaughter(0)->id)==15 && std::abs(pH->getSortedV(1)->getDaughter(0)->id)==13) || (std::abs(pH->getSortedV(0)->getDaughter(0)->id)==13 && std::abs(pH->getSortedV(1)->getDaughter(0)->id)==15)) genFinalState=4; // 2mu2tau
        else if ((std::abs(pH->getSortedV(0)->getDaughter(0)->id)==11 && std::abs(pH->getSortedV(1)->getDaughter(0)->id)==15) || (std::abs(pH->getSortedV(0)->getDaughter(0)->id)==15 && std::abs(pH->getSortedV(1)->getDaughter(0)->id)==11)) genFinalState=5; // 2e2tau
        // 4nu, 4q
        else if (std::abs(pH->getSortedV(0)->getDaughter(0)->id)==std::abs(pH->getSortedV(1)->getDaughter(0)->id)) genFinalState=0;
        else genFinalState=2;
      }
      else genFinalState=2; // WW has no interference
    }
    setVal("genFinalState", genFinalState);
  }
}
void HVVTree::fillAssociatedInfo(ZZCandidate* pH, bool isGen){
  vector<Particle*> AssociatedParticle;
  vector<Particle*> tmpAssociatedParticle;
  Float_t DijetMass=-1;

  Int_t NAssociatedVs=0;
  vector<Particle*> AssociatedV;
  vectorInt AssociatedV_Particle1Index;
  vectorInt AssociatedV_Particle2Index;

  if (pH!=0){
    for (int aa=0; aa<pH->getNAssociatedJets(); aa++){
      Particle* apart = pH->getAssociatedJet(aa);
      tmpAssociatedParticle.push_back(apart);
    }
    for (int aa=0; aa<pH->getNAssociatedLeptons(); aa++){
      Particle* apart = pH->getAssociatedLepton(aa);
      tmpAssociatedParticle.push_back(apart);
    }
    for (int aa=0; aa<pH->getNAssociatedNeutrinos(); aa++){
      Particle* apart = pH->getAssociatedNeutrino(aa);
      tmpAssociatedParticle.push_back(apart);
    }
  }

  while (tmpAssociatedParticle.size()>0){ // Re-sort all associated particles by leading pT (categories are individually soreted, but mixing categories loses this sorting)
    Particle* tmpPart=0;
    int pos=0;
    for (int el=0; el<tmpAssociatedParticle.size(); el++){
      if (tmpPart==0){
        tmpPart = tmpAssociatedParticle.at(el); pos=el;
      }
      else if (tmpPart->pt()<tmpAssociatedParticle.at(el)->pt()){
        tmpPart = tmpAssociatedParticle.at(el); pos=el;
      }// Safer to do in two steps
    }
    AssociatedParticle.push_back(tmpPart);
    tmpAssociatedParticle.erase(tmpAssociatedParticle.begin()+pos);
  }

  NAssociatedVs = (pH!=0 ? pH->getNSortedVs()-2 : 0);
  for (int av=2; av<NAssociatedVs+2; av++){
    Particle* pAV = pH->getSortedV(av);
    AssociatedV.push_back(pAV);
    Particle* avd1 = pAV->getDaughter(0);
    Particle* avd2 = pAV->getDaughter(1);

    for (int aa=0; aa<AssociatedParticle.size(); aa++){
      if (avd1==AssociatedParticle.at(aa)) AssociatedV_Particle1Index.push_back(aa);
      else if (avd2==AssociatedParticle.at(aa)) AssociatedV_Particle2Index.push_back(aa);
    }
  }

  string varname;
  string strcore;

  if (pH->getNAssociatedJets()>1){
    DijetMass = (pH->getAssociatedJet(0)->p4+pH->getAssociatedJet(1)->p4).M();
    varname = "DijetMass";
    if (isGen) varname.insert(0, "Gen");
    setVal(varname, DijetMass);
  }

  varname = "NAssociatedVs";
  if (isGen) varname.insert(0, "Gen");
  setVal(varname, NAssociatedVs);

  strcore = "AssociatedParticle";
  if (isGen) strcore.insert(0, "Gen");
  for (int aa=0; aa<AssociatedParticle.size(); aa++){
    varname = strcore + "Mass"; setVal(varname, AssociatedParticle.at(aa)->m());
    varname = strcore + "Pt"; setVal(varname, AssociatedParticle.at(aa)->pt());
    varname = strcore + "Eta"; setVal(varname, AssociatedParticle.at(aa)->eta());
    varname = strcore + "Phi"; setVal(varname, AssociatedParticle.at(aa)->phi());
    varname = strcore + "Id"; setVal(varname, AssociatedParticle.at(aa)->id);
  }

  strcore = "AssociatedV";
  if (isGen) strcore.insert(0, "Gen");
  for (int aa=0; aa<NAssociatedVs; aa++){
    varname = strcore + "Mass"; setVal(varname, AssociatedV.at(aa)->m());
    varname = strcore + "Pt"; setVal(varname, AssociatedV.at(aa)->pt());
    varname = strcore + "Eta"; setVal(varname, AssociatedV.at(aa)->eta());
    varname = strcore + "Phi"; setVal(varname, AssociatedV.at(aa)->phi());
    varname = strcore + "Id"; setVal(varname, AssociatedV.at(aa)->id);
    varname = strcore + "_Particle1Index"; setVal(varname, AssociatedV_Particle1Index.at(aa));
    varname = strcore + "_Particle2Index"; setVal(varname, AssociatedV_Particle2Index.at(aa));
  }
}

void HVVTree::fillDecayAngles(ZZCandidate* pH, bool isGen){
  Float_t helcosthetaZ1=0, helcosthetaZ2=0, helphi=0, costhetastar=0, phistarZ1=0;

  if (pH!=0) mela::computeAngles(
    pH->getSortedV(0)->getDaughter(0)->p4, pH->getSortedV(0)->getDaughter(0)->id,
    pH->getSortedV(0)->getDaughter(1)->p4, pH->getSortedV(0)->getDaughter(1)->id,
    pH->getSortedV(1)->getDaughter(0)->p4, pH->getSortedV(1)->getDaughter(0)->id,
    pH->getSortedV(1)->getDaughter(1)->p4, pH->getSortedV(1)->getDaughter(1)->id,
    costhetastar, helcosthetaZ1, helcosthetaZ2, helphi, phistarZ1
    );
  // Protect against NaN
  if (!(costhetastar==costhetastar)) costhetastar=0;
  if (!(helcosthetaZ1==helcosthetaZ1)) helcosthetaZ1=0;
  if (!(helcosthetaZ2==helcosthetaZ2)) helcosthetaZ2=0;
  if (!(helphi==helphi)) helphi=0;
  if (!(phistarZ1==phistarZ1)) phistarZ1=0;

  string varname;
  varname = "helcosthetaZ1"; if (isGen) varname.insert(0, "Gen"); setVal(varname, helcosthetaZ1);
  varname = "helcosthetaZ2"; if (isGen) varname.insert(0, "Gen"); setVal(varname, helcosthetaZ2);
  varname = "helphi"; if (isGen) varname.insert(0, "Gen"); setVal(varname, helphi);
  varname = "costhetastar"; if (isGen) varname.insert(0, "Gen"); setVal(varname, costhetastar);
  varname = "phistarZ1"; if (isGen) varname.insert(0, "Gen"); setVal(varname, phistarZ1);
}
//  void HVVTree::fillProductionAngles(Particle* pH, bool isGen=false);

void HVVTree::fillEventVariables(Float_t weight, Int_t passSelection){
  setVal("MC_weight", weight);
  if (options!=0 && options->processRecoInfo()) setVal("isSelected", passSelection);
}

