#ifndef READER_H
#define READER_H

#include "converter.h"


class Reader : public converter{
public:
  Reader(OptionParser* options_);
  ~Reader(){};
  void run();

protected:
  void configure();
  void finalizeRun();
  void readEvent(TTree* tin, int ev, bool isGen, Event& outEvent);

  // This function obtains the result of type "returnType" for the Event ev in the context of Branch "branchname"
  // from an operation "evalVar" that receives a pointer to an Event (most general scenario)
  // and the address of a string "branchname" (so that it can be either set or passed inside) as argument,
  // records the result to the corresponding branch,
  // and returns a bool to indicate the success or failure of the operation.
  // The function prohibits the modification of the Event object.
  template<typename returnType> bool setVariable(const Event* ev, string& branchname, returnType(*evalVar)(const Event*, string&));


};
#endif
