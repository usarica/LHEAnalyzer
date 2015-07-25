#include "../interface/convertLHE.h"
#include "../interface/convertPythia.h"

int main(int argc, char** argv){
  OptionParser options(argc, argv);
  options.printOptionSummary();
  if (options.analysisLevel()==0) convertLHE converter(&options);
  else if (options.analysisLevel()==1) convertPythia converter(&options);
}

