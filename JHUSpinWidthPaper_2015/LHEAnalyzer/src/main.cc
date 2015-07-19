#include "../interface/convertLHE.h"

int main(int argc, char** argv){
  OptionParser options(argc, argv);
  options.printOptionSummary();
  convertLHE converter(&options);
}

