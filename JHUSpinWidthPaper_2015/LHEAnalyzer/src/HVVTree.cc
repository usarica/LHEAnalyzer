#include "../interface/HVVTree.h"


bool HVVTree::reserveBranch(string branchname, BaseTree::BranchTypes branchtype, bool doSetAddress){
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
void HVVTree::bookAllBranches(bool doSetAddress){
  if (!options){
    cerr << "HVVTree::bookAllBranches -> No options are set for the HVVTree!" << endl;
    return;
  }

  reserveBranch("MC_weight", BranchTypes::bFloat, doSetAddress);
  if (options->processGenInfo()){
    reserveBranch("genFinalState", BranchTypes::bInt, doSetAddress);

    bookPtEtaPhiMassIdBranches("Mother", BranchTypes::bVectorDouble, doSetAddress, true, true, true);
    bookPtEtaPhiMassIdBranches("H", BranchTypes::bFloat, doSetAddress, false, true, true);
    bookPtEtaPhiMassIdBranches("Z1", BranchTypes::bFloat, doSetAddress, false, false, true);
    bookPtEtaPhiMassIdBranches("Z2", BranchTypes::bFloat, doSetAddress, false, false, true);
    bookPtEtaPhiMassIdBranches("Za", BranchTypes::bFloat, doSetAddress, false, false, true);
    bookPtEtaPhiMassIdBranches("Zb", BranchTypes::bFloat, doSetAddress, false, false, true);

    bookPtEtaPhiMassIdBranches("AssociatedParticle", BranchTypes::bVectorDouble, doSetAddress, true, false, true);
    reserveBranch("GenDijetMass", BranchTypes::bFloat, doSetAddress);

    reserveBranch("GenNAssociatedVs", BranchTypes::bInt, doSetAddress);
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
  if (options->initializeMELA() || doSetAddress) bookMELABranches(doSetAddress);
  actuateBranches(doSetAddress);
}


void HVVTree::bookPtEtaPhiMassIdBranches(string owner, BaseTree::BranchTypes btype, bool doSetAddress, bool addId, bool usePz, bool isGen){
  vector<string> tmpBranchList;
  getPtEtaPhiMIdBranches(tmpBranchList, owner, addId, usePz, isGen);
  for (int b=0; b<tmpBranchList.size(); b++){
    string branchname = tmpBranchList.at(b);
    if(!addId || branchname.find("Id")==string::npos) reserveBranch(tmpBranchList.at(b), btype, doSetAddress);
    else{
      BaseTree::BranchTypes bInttype = BranchTypes::bInt;
      if (btype==BranchTypes::bVectorDouble) bInttype = BranchTypes::bVectorInt;
      reserveBranch(tmpBranchList.at(b), bInttype, doSetAddress);
    }
  }
}
void HVVTree::getPtEtaPhiMIdBranches(vector<string>& blist, string owner, bool addId, bool usePz, bool isGen){
  string strGen = "Gen";
  vector<string> strtmp;

  strtmp.push_back("Pt");
  if (usePz) strtmp.push_back("Pz");
  else strtmp.push_back("Eta");
  strtmp.push_back("Phi");
  strtmp.push_back("Mass");
  if(addId) strtmp.push_back("Id");

  for (int b=0; b<strtmp.size(); b++){
    string varname = strtmp.at(b);
    varname.insert(0, owner);
    if (isGen) varname.insert(0, strGen);
    blist.push_back(varname);
  }
}
void HVVTree::bookAngularBranches(bool doSetAddress){
  vector<string> tmpBranchList;
  if (options->processGenInfo()){
    if (options->doComputeDecayAngles() || doSetAddress) getAngularBranches(tmpBranchList, 0, true);
    if (options->doComputeVBFAngles() || doSetAddress) getAngularBranches(tmpBranchList, 1, true);
    if (options->doComputeVHAngles() || doSetAddress) getAngularBranches(tmpBranchList, 2, true);
  }
  if (options->processRecoInfo()){
    if (options->doComputeDecayAngles() || doSetAddress) getAngularBranches(tmpBranchList, 0, false);
    if (options->doComputeVBFAngles() || doSetAddress) getAngularBranches(tmpBranchList, 1, false);
    if (options->doComputeVHAngles() || doSetAddress) getAngularBranches(tmpBranchList, 2, false);
  }
  for (int b=0; b<tmpBranchList.size(); b++){
    reserveBranch(tmpBranchList.at(b), BranchTypes::bFloat, doSetAddress);
  }
}
void HVVTree::getAngularBranches(vector<string>& blist, Int_t prodFlag /* 0: Decay, 1: VBF, 2: VH */, bool isGen){
  string strGen = "Gen";
  vector<string> strtmp;
  if (prodFlag==0){
    strtmp.push_back("costhetastar");
    strtmp.push_back("helcosthetaZ1");
    strtmp.push_back("helcosthetaZ2");
    strtmp.push_back("helphi");
    strtmp.push_back("phistarZ1");
  }
  else{
    strtmp.push_back("costhetastar");
    strtmp.push_back("helcosthetaV1");
    strtmp.push_back("helcosthetaV2");
    strtmp.push_back("helphi");
    strtmp.push_back("phistarV1");
  }
  for (int b=0; b<strtmp.size(); b++){
    string varname = strtmp.at(b);
    if (isGen) varname.insert(0, strGen);
    if (prodFlag==1) varname.append("_VBF");
    else if (prodFlag==2) varname.append("_VH");
    blist.push_back(varname);
  }
}


void HVVTree::bookMELABranches(bool doSetAddress){
  vector<string> tmpBranchList = constructMELABranchList(doSetAddress);
  for (int b=0; b<tmpBranchList.size(); b++){
    bool isReserved = reserveBranch(tmpBranchList.at(b), BranchTypes::bFloat, doSetAddress);
    if (isReserved){
      melaProbBranches.push_back(tmpBranchList.at(b));
    }
  }
}
vector<string> HVVTree::constructMELABranchList(bool doSetAddress){
  vector<string> blist;

  // Add gen. prod. MEs
  TVar::Production prod;
  TVar::MatrixElement me;
  if (!doSetAddress){
    prod = options->getSampleProductionId().first;
    me = options->getSampleProductionId().second;
    setupMELASignalMECases(blist, prod, me, true, true, doSetAddress);
  }
  else{
    vector<TVar::Production> prods;
    vector<TVar::MatrixElement> mes;
    mes.push_back(TVar::JHUGen);
    mes.push_back(TVar::MCFM);
    prods.push_back(TVar::JJGG);
    prods.push_back(TVar::JJVBF);
    prods.push_back(TVar::JH);
    prods.push_back(TVar::ZH);
    prods.push_back(TVar::WH);
    prods.push_back(TVar::ttH);
    prods.push_back(TVar::bbH);
    prods.push_back(TVar::ZZGG);
    for (int pp=0; pp<prods.size(); pp++){
      for (int mm=0; mm<mes.size(); mm++) setupMELASignalMECases(blist, prods.at(pp), mes.at(mm), true, true, doSetAddress);
    }
  }

  // Reco. prod. MEs
  me = TVar::JHUGen;
  prod = TVar::ZZGG; setupMELASignalMECases(blist, prod, me, false, true, doSetAddress);
  prod = TVar::JJGG; setupMELASignalMECases(blist, prod, me, false, true, doSetAddress);
  prod = TVar::JJVBF; setupMELASignalMECases(blist, prod, me, false, true, doSetAddress);
  prod = TVar::ZH; setupMELASignalMECases(blist, prod, me, false, true, doSetAddress);
  prod = TVar::WH; setupMELASignalMECases(blist, prod, me, false, true, doSetAddress);
  prod = TVar::JH; setupMELASignalMECases(blist, prod, me, false, true, doSetAddress);
  //prod = TVar::ttH; setupMELASignalMECases(blist, prod, me, false, true);
  //prod = TVar::bbH; setupMELASignalMECases(blist, prod, me, false, true);

  // Gen. decay MEs
  me = TVar::JHUGen;
  prod = TVar::ZZGG; setupMELASignalMECases(blist, prod, me, true, false, doSetAddress);
  me = TVar::MCFM;
  prod = TVar::ZZGG; setupMELASignalMECases(blist, prod, me, true, false, doSetAddress);

  // Reco. decay MEs
  me = TVar::JHUGen;
  prod = TVar::ZZGG; setupMELASignalMECases(blist, prod, me, false, false, doSetAddress);
  me = TVar::MCFM;
  prod = TVar::ZZGG; setupMELASignalMECases(blist, prod, me, false, false, doSetAddress);

  // Bkg MEs
  vector<string> ME_bkg;
  // ggzz_VAMCFM
  ME_bkg.push_back("ggzz_VAMCFM");
  // VVzz_VAMCFM
  //ME_bkg.push_back("VVzz_VAMCFM");
  // bkg_VAMCFM
  ME_bkg.push_back("bkg_VAMCFM");
  ME_bkg.push_back("bkg_VAMCFM_wconst");
  // ggzz_VAMCFM
  ME_bkg.push_back("Gen_ggzz_VAMCFM");
  // VVzz_VAMCFM
  //ME_bkg.push_back("Gen_VVzz_VAMCFM");
  // bkg_VAMCFM
  ME_bkg.push_back("Gen_bkg_VAMCFM");
  if (options->hasRecoDecayME("*") || doSetAddress){
    for (int b=0; b<ME_bkg.size(); b++) blist.push_back(ME_bkg.at(b));
  }

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
  if (options->hasRecoDecayME(chvar) || options->hasRecoDecayME("All") || options->hasRecoDecayME("all") || doSetAddress){
    for (int b=0; b<ME_m4l.size(); b++) blist.push_back(ME_m4l.at(b));
  }
  return blist;
}
void HVVTree::setupMELASignalMECases(vector<string>& accumulatedlist, TVar::Production prod, TVar::MatrixElement me, bool isGen, bool isProdME, bool doSetAddress){
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
    (
    ((options->hasGenProdME("None") || options->hasGenProdME("none")) && isGen && isProdME)
    ||
    ((options->hasRecoProdME("None") || options->hasRecoProdME("none")) && !isGen && isProdME)
    ||
    ((options->hasGenDecayME("None") || options->hasGenDecayME("none")) && isGen && !isProdME)
    ||
    ((options->hasRecoDecayME("None") || options->hasRecoDecayME("none")) && !isGen && !isProdME)
    ) && !doSetAddress // Override "None" options if reading a tree
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
        (
        ((options->hasGenProdME(chvar) || options->hasGenProdME("All") || options->hasGenProdME("all")) && isGen && isProdME)
        ||
        ((options->hasRecoProdME(chvar) || options->hasRecoProdME("All") || options->hasRecoProdME("all")) && !isGen && isProdME)
        ||
        ((options->hasGenDecayME(chvar) || options->hasGenDecayME("All") || options->hasGenDecayME("all")) && isGen && !isProdME)
        ||
        ((options->hasRecoDecayME(chvar) || options->hasRecoDecayME("All") || options->hasRecoDecayME("all")) && !isGen && !isProdME)
        ) || doSetAddress
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
    vector<string> blist = getMELASignalMEBranches(prod, me, gList, v_gCount[0], v_gCount[1], isGen, isProdME, doSetAddress);
    for (int b=0; b<blist.size(); b++) accumulatedlist.push_back(blist.at(b));
  }

  for (int gg=0; gg<sgList; gg++) delete[] gCount[gg];
  delete[] gCount;
}
vector<string> HVVTree::getMELASignalMEBranches(TVar::Production prod, TVar::MatrixElement me, vector<string> gList, vector<int> gCountRe, vector<int> gCountIm, bool isGen, bool isProdME, bool doSetAddress){
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
        else if (gList.at(ai1)=="g1_prime2") tmpVarType = "0_g1prime2";
        tmpVarType.insert(0, "p");

        string strlist;

        strlist = tmpVarType;
        if (me==TVar::MCFM) strlist.append("_VAMCFM");
        else strlist.append("_VAJHU");
        if ((prod==TVar::WH || prod==TVar::ZH) && me!=TVar::MCFM){
          string strlist1 = strlist;
          string strlist2 = strlist;
          strlist1.insert(0, "hadronic_");
          strlist1.insert(0, strcore);
          blist.push_back(strlist1); // Sig
          strlist2.insert(0, "leptonic_");
          strlist2.insert(0, strcore);
          blist.push_back(strlist2); // Sig
        }
        else{
          strlist.insert(0, strcore);
          blist.push_back(strlist); // Sig
        }

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

            if ((prod==TVar::WH || prod==TVar::ZH) && me!=TVar::MCFM){
              string strlist1 = strlist;
              string strlist2 = strlist;
              strlist1.insert(0, "hadronic_");
              strlist1.insert(0, strcore);
              strlist2.insert(0, "leptonic_");
              strlist2.insert(0, strcore);
              if (
                !(im1==im2 && im1==1) // g2_pi2_g4_pi2 == g2g4 when no bkg is involved.
                &&
                (!(im1>im2 && gCount[ai1][1-im1]==1 && gCount[ai2][1-im2]==1) || doSetAddress) // Take the form g2g4_pi2 etc., not g2_pi2_g4
                ){
                blist.push_back(strlist1); // Sig, only
                blist.push_back(strlist2); // Sig, only
              }
            }
            else{
              strlist.insert(0, strcore);
              if (
                !(im1==im2 && im1==1) // g2_pi2_g4_pi2 == g2g4 when no bkg is involved.
                &&
                (!(im1>im2 && gCount[ai1][1-im1]==1 && gCount[ai2][1-im2]==1) || doSetAddress) // Take the form g2g4_pi2 etc., not g2_pi2_g4
                ) blist.push_back(strlist); // Sig, only
            }

            if (me==TVar::MCFM){
              int insertPos = strlist.find("Gen_");
              if (insertPos==string::npos) insertPos=0;
              else insertPos+=4;
              if (prod==TVar::ZZGG) strlist.insert(insertPos, "ggzz_");
              else strlist.insert(insertPos, "VVzz_");
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
  if (pH==0) return;

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
  if (options->doComputeDecayAngles()) fillDecayAngles(pH, isGen);
//  fillProductionAngles(pH, isGen);

  if (melaProbBranches.size()>0) fillMELAProbabilities(pH, isGen); // Do it at the last step
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

  vector<string> varlist;
  getAngularBranches(varlist, 0, isGen);
  for (int b=0; b<varlist.size(); b++){
    string varname = varlist.at(b);
    if (varname.find("costhetastar")!=string::npos) setVal(varname, costhetastar);
    else if (varname.find("helcosthetaZ1")!=string::npos) setVal(varname, helcosthetaZ1);
    else if (varname.find("helcosthetaZ2")!=string::npos) setVal(varname, helcosthetaZ2);
    else if (varname.find("helphi")!=string::npos) setVal(varname, helphi);
    else if (varname.find("phistarZ1")!=string::npos) setVal(varname, phistarZ1);
    else cerr << "HVVTree::fillDecayAngles -> ERROR: Branch " << varname << " is invalid!" << endl;
  }
}
//  void HVVTree::fillProductionAngles(Particle* pH, bool isGen=false);


void HVVTree::fillMELAProbabilities(ZZCandidate* pH, bool isGen){
  if (pH==0) return;

  for (int b=0; b<melaProbBranches.size(); b++){
    string branchname = melaProbBranches.at(b);
    if ((isGen && branchname.find("Gen")==string::npos) || (!isGen && branchname.find("Gen")!=string::npos)) continue;
    Float_t prob = melaHelpers::melaBranchMEInterpreter(pH, branchname);
    if (prob!=0) setVal(branchname, prob);
  }
}


void HVVTree::fillEventVariables(Float_t weight, Int_t passSelection){
  setVal("MC_weight", weight);
  if (options!=0 && options->processRecoInfo()) setVal("isSelected", passSelection);
}

