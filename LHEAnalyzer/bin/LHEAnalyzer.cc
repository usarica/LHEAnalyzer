#include <exception>
#include "Reader.h"
#include "LHEConverter.h"
#include "PythiaConverter.h"

int main(int argc, char** argv){
  try{
    OptionParser options(argc, argv);
    if (!OptionParser::globalHelpFlag){
      //options.printOptionSummary();
      if (options.analysisLevel()==-1) Reader analyzer(&options);
      else if (options.analysisLevel()==0) LHEConverter analyzer(&options);
      else if (options.analysisLevel()==1) PythiaConverter analyzer(&options);
    }
    return 0;
  }
  catch (const std::exception&){
    return 1;
  }
}

