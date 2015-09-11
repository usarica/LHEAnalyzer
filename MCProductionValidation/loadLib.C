{
  gSystem->AddIncludePath("-I$ROOFITSYS/include/");
  gROOT->LoadMacro("src/RooSpinZero.cc+");
  gROOT->LoadMacro("src/RooSpinZero_7DComplex_withAccep_ggH.cc+");
  gROOT->LoadMacro("src/RooSpinZero_5D_VH.cc+");
//  gROOT->LoadMacro("src/RooSpinZero_3D_withAccep_VH.cc+");
  gROOT->LoadMacro("src/ScalarPdfFactory.cc+");
  gROOT->LoadMacro("src/ScalarPdfFactory_ggH.cc+");
//  gROOT->LoadMacro("src/ScalarPdfFactory_VH.cc+");
}
