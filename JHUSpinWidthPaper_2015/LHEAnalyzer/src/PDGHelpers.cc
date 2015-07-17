#include "../interface/PDGHelpers.h"

namespace PDGHelpers{
   double HVVmass = Zmass;
}

bool PDGHelpers::isAQuark(int id){
  if (std::abs(id)<=6 && std::abs(id)>0) return true;
  else return false;
}
bool PDGHelpers::isUpTypeQuark(int id){
  if (std::abs(id)==2 || std::abs(id)==4 || std::abs(id)==6) return true;
  else return false;
}
bool PDGHelpers::isDownTypeQuark(int id){
  if (std::abs(id)==1 || std::abs(id)==3 || std::abs(id)==5) return true;
  else return false;
}
bool PDGHelpers::isALepton(int id){
  if (std::abs(id)==11 || std::abs(id)==13 || std::abs(id)==15) return true;
  else return false;
}
bool PDGHelpers::isANeutrino(int id){
  if (std::abs(id)==12 || std::abs(id)==14 || std::abs(id)==16) return true;
  else return false;
}
bool PDGHelpers::isAGluon(int id){
  if (std::abs(id)==21) return true;
  else return false;
}
bool PDGHelpers::isAZBoson(int id){
  if (std::abs(id)==23) return true;
  else return false;
}
bool PDGHelpers::isAWBoson(int id){
  if (std::abs(id)==24) return true;
  else return false;
}
bool PDGHelpers::isAHiggs(int id){
  if (std::abs(id)==25) return true;
  else return false;
}
void PDGHelpers::setHVVmass(double mymass){
  PDGHelpers::HVVmass=mymass;
}


