#ifndef PDGHELPERS_H
#define PDGHELPERS_H

#include <iostream>
#include <cmath>

namespace PDGHelpers{
  const double Wmass = 80.399;
  const double Zmass = 91.1876;
  const double Topmass = 173.2;
  const double Bottommass = 4.75;
  const double Zeromass = 0;
  const double m_el = 0.00051100;
  const double m_mu = 0.10566;
  const double m_tau = 1.7768;
  const double Wwidth = 2.085;
  const double Zwidth = 2.4952;
  const double Topwidth = 2.;

  extern double HVVmass;

  bool isALepton(int id);
  bool isANeutrino(int id);
  bool isAJet(int id);
  bool isAQuark(int id);
  bool isUpTypeQuark(int id);
  bool isDownTypeQuark(int id);
  bool isAGluon(int id);
  bool isAPhoton(int id);
  bool isAZBoson(int id);
  bool isAWBoson(int id);
  bool isAHiggs(int id);
  void setHVVmass(double mymass);

  int convertPythiaStatus(int pSt);
}

#endif
