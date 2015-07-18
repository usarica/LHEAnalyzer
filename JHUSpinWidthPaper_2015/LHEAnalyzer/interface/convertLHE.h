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

#endif
