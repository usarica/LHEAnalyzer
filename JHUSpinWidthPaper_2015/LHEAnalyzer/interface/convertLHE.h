#ifndef CONVERT_LHE_H
#define CONVERT_LHE_H

#include <iostream>
#include <fstream>
#include <cstdio>
#include <string>
#include <iomanip>
#include <vector>
#include "Event.h"
#include "LHEParticleSmear.h"

using namespace std;

vector<Particle*> readLHEEvent(ifstream& input_lhe, double& weight);

void calculateAngles(TLorentzVector p4H, TLorentzVector p4M11, TLorentzVector p4M12, TLorentzVector p4M21, TLorentzVector p4M22, float& costheta1, float& costheta2, float& phi, float& costhetastar, float& phistar1);

#endif
