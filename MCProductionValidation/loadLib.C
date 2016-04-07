{
  gSystem->AddIncludePath("-I$ROOFITSYS/include/");

  // Higgs JCP mother class
  gROOT->LoadMacro("src/RooSpin.cc+");
  // Spin-0
  gROOT->LoadMacro("src/RooSpinZero.cc+");
  gROOT->LoadMacro("src/RooSpinZero_7DComplex_withAccep_ggH.cc+");
  gROOT->LoadMacro("src/RooSpinZero_5D_VH.cc+");
  gROOT->LoadMacro("src/RooSpinZero_3D_pp_VH.cc+");
  //gROOT->LoadMacro("src/RooSpinZero_3D_withAccep_VH.cc+");
  gROOT->LoadMacro("src/ScalarPdfFactory.cc+");
  gROOT->LoadMacro("src/ScalarPdfFactory_ggH.cc+");
  gROOT->LoadMacro("src/ScalarPdfFactory_VH.cc+");

  // Spin-2
  gROOT->LoadMacro("src/RooSpinTwo.cc+");
  gROOT->LoadMacro("src/RooSpinTwo_7DComplex_HVV.cc+");
  gROOT->LoadMacro("src/TensorPdfFactory.cc+");
  gROOT->LoadMacro("src/TensorPdfFactory_HVV.cc+");
}
