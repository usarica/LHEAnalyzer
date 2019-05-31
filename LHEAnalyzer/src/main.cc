#include "Reader.h"
#include "LHEConverter.h"
#include "PythiaConverter.h"

int main(int argc, char** argv){
  OptionParser options(argc, argv);
  options.printOptionSummary();
  if (options.analysisLevel()==-1) Reader analyzer(&options);
  else if (options.analysisLevel()==0) LHEConverter analyzer(&options);
  else if (options.analysisLevel()==1) PythiaConverter analyzer(&options);
}

